#!/usr/bin/env python3
"""
Extract rRNA chains from a cryo-EM ribosome structure and build per-nucleotide
spatial distance matrices, plus the exact modelled sequence (so contact-map
coordinates and structure coordinates are 1:1).

    pip install gemmi numpy scipy matplotlib
    wget https://files.rcsb.org/download/6EK0.cif.gz && gunzip 6EK0.cif.gz
    python rrna_dist.py 6EK0.cif

Notes
-----
* An 80S entry contains 28S, 5.8S, 5S and 18S in the same file, so one
  download gives you every rRNA at once.
* Residues missing from the model (disordered expansion segments) are kept as
  NaN rows/cols so you can mask them out rather than silently shifting index.
  For mmCIF inputs, the `_pdbx_poly_seq_scheme` table is used when available.
* By default this writes unbinned matrices plus 8, 16, 64 and 128 nt binned
  matrices and plots.
* Coarse-bin distance matrices are block-reduced from nucleotide-level distances
  with a minimum, so a bin pair is close if any nucleotide pair is close.
"""

import argparse
import logging
import os
import shlex
import warnings
from pathlib import Path

# Rough residue-number spans of mature human rRNAs, for auto-labelling chains.
# Span is more stable than the number of modelled residues because structures
# often omit disordered expansion segments.
KNOWN = [(3400, 5600, "28S"), (1500, 1950, "18S"), (140, 200, "5.8S"), (100, 135, "5S")]

np = None
gemmi = None
cdist = None


def guess_name(n):
    for lo, hi, label in KNOWN:
        if lo <= n <= hi:
            return label
    return "other"


def parse_bins(value):
    bins = []
    for raw in value.split(","):
        raw = raw.strip()
        if not raw:
            continue
        try:
            bin_size = int(raw)
        except ValueError as exc:
            raise argparse.ArgumentTypeError(f"invalid bin size: {raw!r}") from exc
        if bin_size < 1:
            raise argparse.ArgumentTypeError("bin sizes must be >= 1")
        bins.append(bin_size)
    if not bins:
        raise argparse.ArgumentTypeError("at least one bin size is required")
    return sorted(set(bins))


def nucleotides(chain):
    """Residues with a C1' atom — catches modified bases (PSU, OMG, 1MA, ...)."""
    return [r for r in chain if r.find_atom("C1'", "*") is not None]


def one_letter(resname):
    info = gemmi.find_tabulated_residue(resname)
    if info is not None and info.one_letter_code:
        return info.one_letter_code.upper()
    return "N"


def scheme_one_letter(resname):
    return one_letter(resname) if resname not in {".", "?"} else "N"


def cif_value(value):
    if value in {".", "?"}:
        return ""
    return value.strip("'\"")


def insertion_code(seqid):
    icode = getattr(seqid, "icode", "")
    if callable(icode):
        icode = icode()
    return str(icode).strip()


def residue_span(residues):
    nums = [r.seqid.num for r in residues]
    return max(nums) - min(nums) + 1


def parse_poly_seq_scheme(cif_path):
    """Return authoritative full-polymer sequence rows keyed by mmCIF chain IDs."""
    tags = None
    rows = []
    in_loop = False
    collecting_rows = False
    wanted_prefix = "_pdbx_poly_seq_scheme."

    with open(cif_path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line == "loop_":
                tags = []
                rows = []
                in_loop = True
                collecting_rows = False
                continue
            if not in_loop:
                continue
            if line.startswith("_"):
                if collecting_rows:
                    break
                tags.append(line)
                continue
            if line.startswith("#"):
                if tags and all(tag.startswith(wanted_prefix) for tag in tags):
                    break
                in_loop = False
                tags = None
                rows = []
                collecting_rows = False
                continue
            if tags and all(tag.startswith(wanted_prefix) for tag in tags):
                collecting_rows = True
                values = shlex.split(line, posix=False)
                if len(values) == len(tags):
                    rows.append(dict(zip(tags, values)))

    if not rows:
        return {}

    by_pdb_strand = {}
    by_asym = {}
    for row in rows:
        seq_id = int(cif_value(row["_pdbx_poly_seq_scheme.seq_id"]))
        entry = {
            "seq_id": seq_id,
            "mon_id": cif_value(row["_pdbx_poly_seq_scheme.mon_id"]),
            "pdb_seq_num": cif_value(row["_pdbx_poly_seq_scheme.pdb_seq_num"]),
            "auth_seq_num": cif_value(row["_pdbx_poly_seq_scheme.auth_seq_num"]),
            "pdb_ins_code": cif_value(row["_pdbx_poly_seq_scheme.pdb_ins_code"]),
            "pdb_mon_id": cif_value(row["_pdbx_poly_seq_scheme.pdb_mon_id"]),
            "auth_mon_id": cif_value(row["_pdbx_poly_seq_scheme.auth_mon_id"]),
        }
        pdb_key = cif_value(row["_pdbx_poly_seq_scheme.pdb_strand_id"])
        asym_key = cif_value(row["_pdbx_poly_seq_scheme.asym_id"])
        if pdb_key:
            by_pdb_strand.setdefault(pdb_key, []).append(entry)
        if asym_key:
            by_asym.setdefault(asym_key, []).append(entry)

    return {
        "pdb_strand_id": normalize_scheme_lookup(by_pdb_strand),
        "asym_id": normalize_scheme_lookup(by_asym),
    }


def normalize_scheme_lookup(lookup):
    for key, chain_rows in lookup.items():
        unique = {}
        for row in chain_rows:
            unique[row["seq_id"]] = row
        lookup[key] = [unique[i] for i in sorted(unique)]
    return lookup


def coord_for_residue(residue, atom):
    if atom == "centroid":
        return np.mean([[a.pos.x, a.pos.y, a.pos.z] for a in residue], axis=0)

    a = residue.find_atom(atom, "*")
    if a is None:
        return None
    return [a.pos.x, a.pos.y, a.pos.z]


def residue_by_pdb_id(residues):
    out = {}
    for residue in residues:
        out[(residue.seqid.num, insertion_code(residue.seqid))] = residue
    return out


def chain_arrays_from_scheme(chain, residues, scheme_rows, atom="C1'"):
    """Return (first_seq_id, seq, coords) using mmCIF polymer sequence mapping."""
    res_by_id = residue_by_pdb_id(residues)
    coords = np.full((len(scheme_rows), 3), np.nan)
    seq = []

    for idx, row in enumerate(scheme_rows):
        seq.append(scheme_one_letter(row["mon_id"]))
        pdb_seq_num = row["pdb_seq_num"]
        if not pdb_seq_num:
            continue
        try:
            residue_num = int(pdb_seq_num)
        except ValueError:
            logging.debug("chain %s: non-integer pdb_seq_num %s", chain.name, pdb_seq_num)
            continue
        residue = res_by_id.get((residue_num, row["pdb_ins_code"]))
        if residue is None:
            continue
        pos = coord_for_residue(residue, atom)
        if pos is not None:
            coords[idx] = pos

    return scheme_rows[0]["seq_id"], "".join(seq), coords


def scheme_residue_labels(scheme_rows):
    labels = []
    for row in scheme_rows:
        label = row["pdb_seq_num"] or row["auth_seq_num"]
        if row["pdb_ins_code"]:
            label = f"{label}{row['pdb_ins_code']}"
        labels.append(label)
    return labels


def chain_arrays(residues, atom="C1'"):
    """Return (first_resid, seq, coords) padded to contiguous numeric residue IDs."""
    seqids = np.array([r.seqid.num for r in residues])
    insertion_codes = [insertion_code(r.seqid) for r in residues]
    if any(insertion_codes):
        raise ValueError("residue insertion codes are not supported by numeric matrix indexing")
    if len(np.unique(seqids)) != len(seqids):
        raise ValueError("duplicate numeric residue IDs are not supported by matrix indexing")
    lo, hi = seqids.min(), seqids.max()
    n = hi - lo + 1
    coords = np.full((n, 3), np.nan)
    seq = ["-"] * n
    for r in residues:
        i = r.seqid.num - lo
        seq[i] = one_letter(r.name)
        pos = coord_for_residue(r, atom)
        if pos is None:
            continue
        coords[i] = pos
    return lo, "".join(seq), coords


def block_min_distance(D, bin_size):
    """Reduce an nt-level distance matrix to bins by min distance per block."""
    pad = (-len(D)) % bin_size
    if pad:
        D = np.pad(D, ((0, pad), (0, pad)), constant_values=np.nan)

    nb = len(D) // bin_size
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmin(D.reshape(nb, bin_size, nb, bin_size), axis=(1, 3))


def binned_arrays(seq, coords, offset, bin_size, row_ids=None, resid_labels=None):
    """Return binned coordinates and row metadata for a given bin size."""
    n = len(coords)
    if row_ids is None:
        row_ids = [offset + i for i in range(n)]
    if resid_labels is None:
        resid_labels = [str(row_id) for row_id in row_ids]
    pad = (-n) % bin_size
    if pad:
        padded = np.vstack([coords, np.full((pad, 3), np.nan)])
    else:
        padded = coords

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        binned = np.nanmean(padded.reshape(-1, bin_size, 3), axis=1)

    rows = []
    for row_idx in range(len(binned)):
        start_idx = row_idx * bin_size
        end_idx = min(start_idx + bin_size, n)
        chunk = seq[start_idx:end_idx]
        window = coords[start_idx:end_idx]
        window_labels = [label for label in resid_labels[start_idx:end_idx] if label]
        rows.append(
            {
                "row": row_idx,
                "seq_id_start": row_ids[start_idx],
                "seq_id_end": row_ids[end_idx - 1],
                "resid_start": window_labels[0] if window_labels else "",
                "resid_end": window_labels[-1] if window_labels else "",
                "modelled": int(np.isfinite(window).all(axis=1).sum()),
                "seq": chunk,
            }
        )
    return binned, rows


def write_fasta(path, tag, offset, seq, modelled):
    with path.open("w") as fh:
        fh.write(f">{tag} first_seq_id={offset} modelled={modelled}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")


def write_bins_tsv(path, rows):
    with path.open("w") as fh:
        fh.write("row\tseq_id_start\tseq_id_end\tresid_start\tresid_end\tmodelled\tseq\n")
        for row in rows:
            fh.write(
                f"{row['row']}\t{row['seq_id_start']}\t{row['seq_id_end']}\t"
                f"{row['resid_start']}\t{row['resid_end']}\t{row['modelled']}\t{row['seq']}\n"
            )


def plot_distance_map(path, D, tag, bin_size, plot_max_size):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_d = D
    stride = 1
    if len(D) > plot_max_size:
        stride = int(np.ceil(len(D) / plot_max_size))
        plot_d = D[::stride, ::stride]

    finite = np.isfinite(plot_d)
    vmax = np.nanpercentile(plot_d[finite], 95) if finite.any() else 1.0

    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    im = ax.imshow(plot_d, cmap="viridis_r", origin="lower", vmin=0, vmax=vmax)
    title = f"{tag} physical distance map"
    if stride > 1:
        title += f" (shown every {stride} bins)"
    ax.set_title(title)
    ax.set_xlabel(f"matrix row ({bin_size} nt/bin)")
    ax.set_ylabel(f"matrix row ({bin_size} nt/bin)")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("distance (Angstrom)")
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_structure_trace(path, xyz, tag, bin_size):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    keep = np.isfinite(xyz).all(axis=1)
    coords = xyz[keep]
    if len(coords) == 0:
        return False

    colors = np.flatnonzero(keep)
    stride = max(1, int(np.ceil(len(coords) / 6000)))
    coords = coords[::stride]
    colors = colors[::stride]

    fig = plt.figure(figsize=(8, 7), constrained_layout=True)
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(xyz[:, 0], xyz[:, 1], xyz[:, 2], color="0.70", linewidth=0.6)
    sc = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        coords[:, 2],
        c=colors,
        cmap="plasma",
        s=4,
        depthshade=False,
    )
    ax.set_title(f"{tag} 3D coordinate trace")
    ax.set_xlabel("x (Angstrom)")
    ax.set_ylabel("y (Angstrom)")
    ax.set_zlabel("z (Angstrom)")
    cbar = fig.colorbar(sc, ax=ax, shrink=0.7)
    cbar.set_label(f"matrix row ({bin_size} nt/bin)")
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return True


def setup_logging(outdir, verbose=False):
    log_path = outdir / "rrna_dist.log"
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(message)s",
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler(log_path, mode="w"),
        ],
    )
    logging.info("writing log to %s", log_path)


def load_runtime_deps():
    global np, gemmi, cdist

    try:
        import numpy as _np
        import gemmi as _gemmi
        from scipy.spatial.distance import cdist as _cdist
    except ModuleNotFoundError as exc:
        raise SystemExit(
            f"missing dependency {exc.name!r}; create the environment with "
            "`conda env create -f nc_rna_benchmarking/environment.yml`"
        ) from exc

    np = _np
    gemmi = _gemmi
    cdist = _cdist


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cif")
    ap.add_argument("--outdir", default="rrna_dist_outputs")
    ap.add_argument("--atom", default="C1'", help="C1' (default), P, or 'centroid'")
    ap.add_argument("--min-len", type=int, default=100)
    ap.add_argument(
        "--bins",
        type=parse_bins,
        default=parse_bins("1,8,16,64,128"),
        help="comma-separated bin sizes in nt (default: 1,8,16,64,128)",
    )
    ap.add_argument(
        "--bin",
        type=int,
        action="append",
        help="legacy single-bin option; if supplied, overrides --bins",
    )
    ap.add_argument("--no-plots", action="store_true", help="skip PNG plot generation")
    ap.add_argument("--plot-max-size", type=int, default=2200, help="max plotted rows/cols")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    bins = sorted(set(args.bin)) if args.bin else args.bins
    if any(bin_size < 1 for bin_size in bins):
        ap.error("bin sizes must be >= 1")
    load_runtime_deps()
    outdir = Path(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir, args.verbose)

    logging.info("reading structure %s", args.cif)
    logging.info("atom=%s min_len=%d bins=%s plots=%s", args.atom, args.min_len, bins, not args.no_plots)
    schemes = parse_poly_seq_scheme(args.cif)
    if schemes:
        logging.info(
            "using _pdbx_poly_seq_scheme (%d pdb_strand_id, %d asym_id entries)",
            len(schemes["pdb_strand_id"]),
            len(schemes["asym_id"]),
        )
    else:
        logging.info("_pdbx_poly_seq_scheme not found; falling back to numeric residue spans")
    st = gemmi.read_structure(args.cif)
    st.setup_entities()
    model = st[0]

    seen = {}
    for chain in model:
        res = nucleotides(chain)
        if len(res) < args.min_len:
            logging.debug("skipping chain %s: %d nucleotide residues", chain.name, len(res))
            continue
        scheme_rows = schemes.get("pdb_strand_id", {}).get(chain.name)
        scheme_key = "pdb_strand_id"
        if scheme_rows is None:
            scheme_rows = schemes.get("asym_id", {}).get(chain.name)
            scheme_key = "asym_id"
        span = len(scheme_rows) if scheme_rows else residue_span(res)
        label = guess_name(span)
        if label in seen:
            logging.info("found another %s-like chain; chain name will keep outputs unique", label)
        seen[label] = True

        if scheme_rows:
            offset, seq, xyz = chain_arrays_from_scheme(chain, res, scheme_rows, args.atom)
            row_ids = [row["seq_id"] for row in scheme_rows]
            resid_labels = scheme_residue_labels(scheme_rows)
            axis_source = f"_pdbx_poly_seq_scheme.{scheme_key}"
        else:
            try:
                offset, seq, xyz = chain_arrays(res, args.atom)
            except ValueError as exc:
                raise ValueError(f"chain {chain.name}: {exc}") from exc
            row_ids = [offset + i for i in range(len(seq))]
            resid_labels = [str(row_id) for row_id in row_ids]
            axis_source = "numeric residue IDs"
        chain_tag = f"{label}_chain{chain.name}"
        D_full = cdist(xyz, xyz)               # NaN propagates for unmodelled residues
        write_fasta(outdir / f"{chain_tag}.fa", chain_tag, offset, seq, len(res))
        logging.info(
            "chain=%s label=%s span=%d modelled=%d axis=%s first_seq_id=%d gaps=%d",
            chain.name,
            label,
            span,
            len(res),
            axis_source,
            offset,
            seq.count("-"),
        )

        for bin_size in bins:
            bin_dir = outdir / f"bin_{bin_size:03d}nt"
            bin_dir.mkdir(parents=True, exist_ok=True)
            tag = f"{chain_tag}_bin{bin_size}nt"

            out_xyz, rows = binned_arrays(seq, xyz, offset, bin_size, row_ids, resid_labels)
            if bin_size == 1:
                D = D_full
            else:
                D = block_min_distance(D_full, bin_size)

            np.save(bin_dir / f"{tag}.dist.npy", D.astype(np.float32))
            np.save(bin_dir / f"{tag}.xyz.npy", out_xyz.astype(np.float32))
            write_bins_tsv(bin_dir / f"{tag}.bins.tsv", rows)

            if not args.no_plots:
                plot_distance_map(
                    bin_dir / f"{tag}.distance_map.png",
                    D,
                    tag,
                    bin_size,
                    args.plot_max_size,
                )
                made_trace = plot_structure_trace(
                    bin_dir / f"{tag}.structure_trace.png",
                    out_xyz,
                    tag,
                    bin_size,
                )
                if not made_trace:
                    logging.warning("no finite coordinates for %s; skipped structure trace", tag)

            logging.info("wrote %s matrix=%s", tag, D.shape)

    logging.info("done; wrote outputs to %s", outdir)
    print(f"\nwrote to {outdir}/")
    plot_note = " and PNG plots" if not args.no_plots else ""
    print(f"Root contains per-chain .fa files; each bin folder contains .dist.npy, .xyz.npy, .bins.tsv{plot_note}.")


if __name__ == "__main__":
    main()
