#!/usr/bin/env python3
"""
Extract rRNA chains from a cryo-EM ribosome structure and build per-nucleotide
and windowed spatial distance matrices, plus the exact modelled sequence, so
matrix coordinates and structure coordinates are 1:1.

    pip install gemmi numpy scipy matplotlib
    wget https://files.rcsb.org/download/4V6X.cif.gz && gunzip 4V6X.cif.gz
    python rrna_dist.py 4V6X.cif --karr-seq

KARR-seq convention
-------------------
Wu et al., Nat. Biotechnol. 2024 built their cryo-EM ground truth by taking the
centroid of each 5-nt window and computing Euclidean distances between window
centroids, then looking up values using the midpoint of each arm of a chimeric
read.  `--karr-seq` is a preset for exactly that:

    --bins 1,5 --atom centroid --centroid-weight atom --reduce centroid

They used 4V6X for 18S/28S (also 7QXB for TERC, 6QX9 for U1 snRNA).  Use 4V6X if you
want numbers directly comparable to the published benchmark; 6EK0 is newer and
better resolved if you only want the best available human structure.

Reduction modes (`--reduce`)
----------------------------
centroid  Distance between window centroids.  Matches KARR-seq.  Sensible while
          windows are small enough to be compact (roughly <= 10 nt).
min       Closest nucleotide pair in each block.  Answers "do these two regions
          touch anywhere".  Needed once windows get large, but values are NOT
          comparable across bin sizes (more pairs = more chances at a small
          minimum) and partially modelled windows are biased outward.
contact   Fraction of modelled nucleotide pairs below --contact-cutoff.  Closest
          analogue to what a ligation assay actually measures, and normalising
          by modelled pairs removes the partial-modelling bias.

Notes
-----
* An 80S entry contains 28S, 5.8S, 5S and 18S in one file, so a single download
  gives every rRNA at once.
* Residues missing from the model (disordered expansion segments) are kept as
  NaN rows/cols so you can mask them rather than silently shifting index.
  For mmCIF inputs, `_pdbx_poly_seq_scheme` provides the authoritative mapping.
* Each bin directory contains a `.lookup.tsv` giving residue -> matrix row, so
  chimeric-read midpoints can be mapped onto the matrix without recomputing
  offsets by hand:

      import numpy as np, csv
      row = {r["resid"]: int(r["row"]) for r in
             csv.DictReader(open("...lookup.tsv"), delimiter="\t")}
      D = np.load("...dist.npy")
      D[row["1055"], row["1462"]]
"""

import argparse
import logging
import os
import re
import shlex
import warnings
from pathlib import Path

try:
    from download_pdb_cifs import (
        default_outdir as default_structures_dir,
        download_record,
        read_accessions,
        safe_name,
        write_manifest,
    )
except ModuleNotFoundError:
    from nc_rna_benchmarking.download_pdb_cifs import (
        default_outdir as default_structures_dir,
        download_record,
        read_accessions,
        safe_name,
        write_manifest,
    )

# Rough residue-number spans of human RNAs. These are a fallback when mmCIF
# entity descriptions do not identify the RNA species.
# Caveat: bacterial 16S (~1540 nt) and mitochondrial 12S/16S fall inside or near
# these windows and will be mislabelled.  Use --label-from-chain for those.
KNOWN = [
    (3400, 5600, "28S"),
    (1500, 1950, "18S"),
    (400, 500, "TERC"),
    (140, 200, "5.8S"),
    (100, 135, "5S"),
]

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


def is_int(value):
    try:
        int(value)
    except ValueError:
        return False
    return True


def parse_entity_metadata(cif_path):
    """Read enough mmCIF metadata to classify polymer chains."""
    metadata = {"entities": {}, "poly_types": {}, "asym_to_entity": {}}
    tags = None
    in_loop = False
    category = None

    with open(cif_path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line == "loop_":
                tags = []
                in_loop = True
                category = None
                continue
            if not in_loop:
                continue
            if line.startswith("_"):
                tags.append(line)
                if line.startswith("_entity."):
                    category = "_entity"
                elif line.startswith("_entity_poly."):
                    category = "_entity_poly"
                elif line.startswith("_struct_asym."):
                    category = "_struct_asym"
                continue
            if line.startswith("#"):
                in_loop = False
                tags = None
                category = None
                continue
            if line.startswith(";") or category is None:
                continue

            try:
                values = shlex.split(line)
            except ValueError:
                continue
            if not values:
                continue

            if category == "_entity" and len(values) >= 4 and is_int(values[0]):
                metadata["entities"][values[0]] = {
                    "type": cif_value(values[1]),
                    "description": cif_value(values[3]),
                }
            elif category == "_entity_poly" and len(values) >= 2 and is_int(values[0]):
                metadata["poly_types"][values[0]] = cif_value(values[1])
            elif category == "_struct_asym" and len(values) >= 4:
                metadata["asym_to_entity"][cif_value(values[0])] = cif_value(values[3])

    return metadata


def description_label(description):
    desc = description.lower()
    if re.search(r"\bu1\s+(snrna|small nuclear rna)\b", desc):
        return "U1"
    if re.search(r"\bu2\s+(snrna|small nuclear rna)\b", desc):
        return "U2"
    if re.search(r"\bu4\s+(snrna|small nuclear rna)\b", desc):
        return "U4"
    if re.search(r"\bu5\s+(snrna|small nuclear rna)\b", desc):
        return "U5"
    if re.search(r"\bu6\s+(snrna|small nuclear rna)\b", desc):
        return "U6"
    if "terc" in desc or "telomerase rna" in desc or "telomerase rna component" in desc:
        return "TERC"
    if "5.8s" in desc or "5.8 s" in desc:
        return "5.8S"
    if "28s" in desc or "28 s" in desc:
        return "28S"
    if "18s" in desc or "18 s" in desc:
        return "18S"
    if re.search(r"\b5s\s+(rrna|ribosomal rna)\b", desc):
        return "5S"
    return None


def is_rna_entity(entity_id, metadata):
    if not entity_id:
        return None
    poly_type = metadata["poly_types"].get(entity_id, "").lower()
    description = metadata["entities"].get(entity_id, {}).get("description", "").lower()
    if "ribonucleotide" in poly_type or "rna" in description:
        return True
    if poly_type or description:
        return False
    return None


def entity_description(entity_id, metadata):
    if not entity_id:
        return ""
    return metadata["entities"].get(entity_id, {}).get("description", "")


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
            "entity_id": cif_value(row["_pdbx_poly_seq_scheme.entity_id"]),
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
    """Return (xyz, n_atoms). n_atoms is the weight used for atom-weighted centroids."""
    positions = [[a.pos.x, a.pos.y, a.pos.z] for a in residue]
    if atom == "centroid":
        if not positions:
            return None, 0
        return np.mean(positions, axis=0), len(positions)

    a = residue.find_atom(atom, "*")
    if a is None:
        return None, 0
    return [a.pos.x, a.pos.y, a.pos.z], len(positions)


def residue_by_pdb_id(residues):
    out = {}
    for residue in residues:
        out[(residue.seqid.num, insertion_code(residue.seqid))] = residue
    return out


def chain_arrays_from_scheme(chain, residues, scheme_rows, atom="C1'"):
    """Return (first_seq_id, seq, coords, weights) using mmCIF polymer sequence mapping."""
    res_by_id = residue_by_pdb_id(residues)
    coords = np.full((len(scheme_rows), 3), np.nan)
    weights = np.zeros(len(scheme_rows))
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
        pos, natoms = coord_for_residue(residue, atom)
        if pos is not None:
            coords[idx] = pos
            weights[idx] = natoms

    return scheme_rows[0]["seq_id"], "".join(seq), coords, weights


def scheme_residue_labels(scheme_rows):
    labels = []
    for row in scheme_rows:
        label = row["pdb_seq_num"] or row["auth_seq_num"]
        if row["pdb_ins_code"]:
            label = f"{label}{row['pdb_ins_code']}"
        labels.append(label)
    return labels


def chain_arrays(residues, atom="C1'"):
    """Return (first_resid, seq, coords, weights) padded to contiguous residue IDs."""
    seqids = np.array([r.seqid.num for r in residues])
    insertion_codes = [insertion_code(r.seqid) for r in residues]
    if any(insertion_codes):
        raise ValueError("residue insertion codes are not supported by numeric matrix indexing")
    if len(np.unique(seqids)) != len(seqids):
        raise ValueError("duplicate numeric residue IDs are not supported by matrix indexing")
    lo, hi = seqids.min(), seqids.max()
    n = hi - lo + 1
    coords = np.full((n, 3), np.nan)
    weights = np.zeros(n)
    seq = ["-"] * n
    for r in residues:
        i = r.seqid.num - lo
        seq[i] = one_letter(r.name)
        pos, natoms = coord_for_residue(r, atom)
        if pos is None:
            continue
        coords[i] = pos
        weights[i] = natoms
    return lo, "".join(seq), coords, weights


def window_centroids(coords, weights, bin_size, centroid_weight="atom"):
    """Centroid of each window.

    KARR-seq takes the centroid of each 5-nt window over its atoms, so purines
    (more atoms) pull the centroid slightly more than pyrimidines. Pass
    centroid_weight='residue' to weight every modelled residue equally instead.
    Windows with no modelled residue come back as NaN.
    """
    n = len(coords)
    pad = (-n) % bin_size
    c = coords
    w = np.ones(n) if centroid_weight == "residue" else np.asarray(weights, dtype=float)
    if pad:
        c = np.vstack([c, np.full((pad, 3), np.nan)])
        w = np.concatenate([w, np.zeros(pad)])

    finite = np.isfinite(c).all(axis=1)
    w = np.where(finite, w, 0.0)
    c = np.where(finite[:, None], c, 0.0)

    nb = len(c) // bin_size
    weighted = (c * w[:, None]).reshape(nb, bin_size, 3).sum(axis=1)
    totals = w.reshape(nb, bin_size).sum(axis=1)

    out = np.full((nb, 3), np.nan)
    ok = totals > 0
    out[ok] = weighted[ok] / totals[ok, None]
    return out


def mask_centroids_by_modelled_fraction(centroids, rows, min_fraction):
    if min_fraction <= 0:
        return centroids, 0

    out = centroids.copy()
    masked = 0
    for row in rows:
        if row["modelled_fraction"] < min_fraction and np.isfinite(out[row["row"]]).all():
            out[row["row"]] = np.nan
            masked += 1
    return out, masked


def block_min_distance(D, bin_size):
    """Reduce an nt-level distance matrix to bins by min distance per block."""
    pad = (-len(D)) % bin_size
    if pad:
        D = np.pad(D, ((0, pad), (0, pad)), constant_values=np.nan)

    nb = len(D) // bin_size
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        return np.nanmin(D.reshape(nb, bin_size, nb, bin_size), axis=(1, 3))


def block_contact_fraction(D, bin_size, cutoff):
    """Fraction of MODELLED nucleotide pairs in each block below `cutoff`.

    Normalising by modelled pairs rather than total pairs means a window with
    half its residues unmodelled is not automatically reported as low-contact.
    Blocks with no modelled pair at all come back as NaN.
    """
    pad = (-len(D)) % bin_size
    if pad:
        D = np.pad(D, ((0, pad), (0, pad)), constant_values=np.nan)

    valid = np.isfinite(D)
    hits = valid & (D < cutoff)

    nb = len(D) // bin_size
    shape = (nb, bin_size, nb, bin_size)
    n_valid = valid.reshape(shape).sum(axis=(1, 3))
    n_hits = hits.reshape(shape).sum(axis=(1, 3))

    out = np.full((nb, nb), np.nan)
    ok = n_valid > 0
    out[ok] = n_hits[ok] / n_valid[ok]
    return out


def mask_diagonal(D, min_separation):
    """NaN out |i-j| < min_separation.

    Regions close in sequence are trivially close in space, so the near-diagonal
    band carries no tertiary-structure information and will dominate any ROC or
    correlation against experimental data if left in.
    """
    if min_separation <= 0:
        return D
    n = len(D)
    i, j = np.indices((n, n))
    out = D.copy()
    out[np.abs(i - j) < min_separation] = np.nan
    return out


def mask_diagonal_band(D, half_width):
    out = D.copy()
    n = len(out)
    i, j = np.indices((n, n))
    out[np.abs(i - j) <= half_width] = np.nan
    return out


def bin_row_metadata(seq, coords, bin_size, row_ids, resid_labels):
    """Per-row provenance: which residues went in, and how many were modelled."""
    n = len(coords)
    nb = int(np.ceil(n / bin_size))
    rows = []
    for row_idx in range(nb):
        start_idx = row_idx * bin_size
        end_idx = min(start_idx + bin_size, n)
        window = coords[start_idx:end_idx]
        window_labels = [label for label in resid_labels[start_idx:end_idx] if label]
        n_modelled = int(np.isfinite(window).all(axis=1).sum())
        width = end_idx - start_idx
        rows.append(
            {
                "row": row_idx,
                "seq_id_start": row_ids[start_idx],
                "seq_id_end": row_ids[end_idx - 1],
                "resid_start": window_labels[0] if window_labels else "",
                "resid_end": window_labels[-1] if window_labels else "",
                "modelled": n_modelled,
                "window": width,
                "modelled_fraction": n_modelled / width,
                "seq": seq[start_idx:end_idx],
            }
        )
    return rows


def write_fasta(path, tag, offset, seq, modelled):
    with path.open("w") as fh:
        fh.write(f">{tag} first_seq_id={offset} modelled={modelled}\n")
        for i in range(0, len(seq), 60):
            fh.write(seq[i:i + 60] + "\n")


def write_bins_tsv(path, rows):
    with path.open("w") as fh:
        fh.write(
            "row\tseq_id_start\tseq_id_end\tresid_start\tresid_end\t"
            "modelled\twindow\tmodelled_fraction\tseq\n"
        )
        for row in rows:
            fh.write(
                f"{row['row']}\t{row['seq_id_start']}\t{row['seq_id_end']}\t"
                f"{row['resid_start']}\t{row['resid_end']}\t{row['modelled']}\t"
                f"{row['window']}\t{row['modelled_fraction']:.6g}\t{row['seq']}\n"
            )


def write_lookup_tsv(path, row_ids, resid_labels, coords, bin_size):
    """residue -> matrix row, for mapping chimeric-read midpoints onto the matrix."""
    modelled = np.isfinite(coords).all(axis=1)
    with path.open("w") as fh:
        fh.write("seq_id\tresid\trow\tmodelled\n")
        for idx, (seq_id, label) in enumerate(zip(row_ids, resid_labels)):
            fh.write(f"{seq_id}\t{label}\t{idx // bin_size}\t{int(modelled[idx])}\n")


def plot_distance_map(path, D, tag, bin_size, plot_max_size, reduce_mode, cutoff):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plot_d = D
    stride = 1
    if len(D) > plot_max_size:
        stride = int(np.ceil(len(D) / plot_max_size))
        plot_d = D[::stride, ::stride]

    finite = np.isfinite(plot_d)
    contact_mode = reduce_mode == "contact"

    if contact_mode:
        cmap_name, vmin = "magma", 0.0
        vmax = np.nanpercentile(plot_d[finite], 99) if finite.any() else 1.0
        vmax = max(vmax, 1e-6)
        cbar_label = f"fraction of modelled pairs < {cutoff:g} A"
        title = f"{tag} contact density"
    else:
        cmap_name, vmin = "viridis_r", 0.0
        vmax = np.nanpercentile(plot_d[finite], 95) if finite.any() else 1.0
        cbar_label = "distance (Angstrom)"
        title = f"{tag} physical distance map"

    if stride > 1:
        title += f" (shown every {stride} bins)"

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad("#f2f2f2")

    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    im = ax.imshow(
        np.ma.masked_invalid(plot_d),
        cmap=cmap,
        origin="lower",
        vmin=vmin,
        vmax=vmax,
        interpolation="nearest",
    )
    ax.set_title(title)
    ax.set_xlabel(f"matrix row ({bin_size} nt/bin)")
    ax.set_ylabel(f"matrix row ({bin_size} nt/bin)")
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(cbar_label)
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
        force=True,
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


def bins_slug(bins):
    return "bins-" + "-".join(str(bin_size) for bin_size in bins)


def mode_slug(args):
    if args.karr_seq:
        return "karrseq"
    return f"{args.reduce}_{args.atom.replace(chr(39), '')}"


def output_dir_for_run(base_outdir, cif_path, bins, args, metadata=""):
    stem = Path(cif_path).stem
    label = safe_name(metadata) if metadata else stem
    if label.upper() == stem.upper():
        prefix = stem
    else:
        prefix = f"{stem}_{label}"
    return Path(base_outdir) / f"{prefix}_{mode_slug(args)}_{bins_slug(bins)}"


def process_cif(cif_path, outdir, args, bins, logging_deferred=None):
    os.makedirs(outdir, exist_ok=True)
    setup_logging(outdir, args.verbose)

    logging.info("reading structure %s", cif_path)
    logging.info(
        "atom=%s reduce=%s centroid_weight=%s min_len=%d bins=%s min_sep=%d min_bin_modelled_fraction=%.3g plots=%s",
        args.atom, args.reduce, args.centroid_weight, args.min_len,
        bins, args.min_separation, args.min_bin_modelled_fraction, not args.no_plots,
    )
    if args.karr_seq:
        logging.info("KARR-seq preset active (Wu et al. 2024: 5-nt window centroids)")
    if logging_deferred:
        logging.warning(logging_deferred)
    if args.reduce == "contact":
        logging.info("contact cutoff = %.1f A (calibrate against known base pairs)", args.contact_cutoff)
    if args.reduce == "min" and len(bins) > 1:
        logging.warning(
            "min-reduced values are not comparable across bin sizes; "
            "do not put different bin sizes on a shared colour scale"
        )

    entity_metadata = parse_entity_metadata(cif_path)
    logging.info(
        "read entity metadata: %d entities, %d polymer types, %d asym mappings",
        len(entity_metadata["entities"]),
        len(entity_metadata["poly_types"]),
        len(entity_metadata["asym_to_entity"]),
    )
    schemes = parse_poly_seq_scheme(cif_path)
    if schemes:
        logging.info(
            "using _pdbx_poly_seq_scheme (%d pdb_strand_id, %d asym_id entries)",
            len(schemes["pdb_strand_id"]),
            len(schemes["asym_id"]),
        )
    else:
        logging.info("_pdbx_poly_seq_scheme not found; falling back to numeric residue spans")

    st = gemmi.read_structure(str(cif_path))
    st.setup_entities()
    model = st[0]

    seen = {}
    for chain in model:
        scheme_rows = schemes.get("pdb_strand_id", {}).get(chain.name)
        scheme_key = "pdb_strand_id"
        if scheme_rows is None:
            scheme_rows = schemes.get("asym_id", {}).get(chain.name)
            scheme_key = "asym_id"

        entity_id = scheme_rows[0]["entity_id"] if scheme_rows else entity_metadata["asym_to_entity"].get(chain.name)
        metadata_rna = is_rna_entity(entity_id, entity_metadata)
        description = entity_description(entity_id, entity_metadata)
        if metadata_rna is False:
            logging.debug("skipping chain %s: entity %s is %s", chain.name, entity_id, description)
            continue

        res = nucleotides(chain)
        if len(res) < args.min_len:
            logging.debug("skipping chain %s: %d nucleotide residues", chain.name, len(res))
            continue

        span = len(scheme_rows) if scheme_rows else residue_span(res)

        if args.label_from_chain:
            label = "chain"
            label_source = "chain ID"
        else:
            label = description_label(description)
            label_source = "entity description"
            if label is None:
                label = guess_name(span)
                label_source = "span fallback"
            if label == "other":
                logging.warning(
                    "chain %s (span %d) matches no known human rRNA; labelled 'other'",
                    chain.name, span,
                )
        seen[label] = seen.get(label, 0) + 1
        if seen[label] > 1:
            logging.info(
                "chain %s is %s-like #%d; chain ID keeps output filenames unique",
                chain.name, label, seen[label],
            )

        if scheme_rows:
            offset, seq, xyz, weights = chain_arrays_from_scheme(chain, res, scheme_rows, args.atom)
            row_ids = [row["seq_id"] for row in scheme_rows]
            resid_labels = scheme_residue_labels(scheme_rows)
            axis_source = f"_pdbx_poly_seq_scheme.{scheme_key}"
        else:
            try:
                offset, seq, xyz, weights = chain_arrays(res, args.atom)
            except ValueError as exc:
                raise ValueError(f"chain {chain.name}: {exc}") from exc
            row_ids = [offset + i for i in range(len(seq))]
            resid_labels = [str(row_id) for row_id in row_ids]
            axis_source = "numeric residue IDs"

        chain_tag = f"{label}_chain{chain.name}"
        n_placed = int(np.isfinite(xyz).all(axis=1).sum())
        write_fasta(outdir / f"{chain_tag}.fa", chain_tag, offset, seq, n_placed)
        logging.info(
            "chain=%s label=%s label_source=%s entity=%s span=%d rows=%d placed=%d unmodelled=%d axis=%s first_seq_id=%s",
            chain.name, label, label_source, entity_id or "", span, len(seq), n_placed, len(seq) - n_placed,
            axis_source, offset,
        )

        D_full = None
        if args.reduce in {"min", "contact"}:
            D_full = cdist(xyz, xyz)   # NaN propagates for unmodelled residues

        for bin_size in bins:
            bin_dir = outdir / f"bin_{bin_size:03d}nt"
            bin_dir.mkdir(parents=True, exist_ok=True)
            tag = f"{chain_tag}_bin{bin_size}nt"

            out_xyz = window_centroids(xyz, weights, bin_size, args.centroid_weight)
            rows = bin_row_metadata(seq, xyz, bin_size, row_ids, resid_labels)

            if args.reduce == "centroid":
                out_xyz, masked_bins = mask_centroids_by_modelled_fraction(
                    out_xyz, rows, args.min_bin_modelled_fraction
                )
                if masked_bins:
                    logging.info(
                        "%s masked %d centroid bins below modelled fraction %.3g",
                        tag, masked_bins, args.min_bin_modelled_fraction,
                    )
                D = cdist(out_xyz, out_xyz)
            elif args.reduce == "min":
                D = D_full if bin_size == 1 else block_min_distance(D_full, bin_size)
            else:
                D = block_contact_fraction(D_full, bin_size, args.contact_cutoff)

            D = mask_diagonal(D, args.min_separation)

            np.save(bin_dir / f"{tag}.dist.npy", D.astype(np.float32))
            np.save(bin_dir / f"{tag}.xyz.npy", out_xyz.astype(np.float32))
            write_bins_tsv(bin_dir / f"{tag}.bins.tsv", rows)
            write_lookup_tsv(bin_dir / f"{tag}.lookup.tsv", row_ids, resid_labels, xyz, bin_size)

            if not args.no_plots:
                plot_distance_map(
                    bin_dir / f"{tag}.distance_map.png",
                    D, tag, bin_size, args.plot_max_size,
                    args.reduce, args.contact_cutoff,
                )
                plot_distance_map(
                    bin_dir / f"{tag}.distance_map.masked_diagonal.png",
                    mask_diagonal_band(D, args.masked_diagonal_width),
                    f"{tag} diagonal masked (+/- {args.masked_diagonal_width})", bin_size,
                    args.plot_max_size, args.reduce, args.contact_cutoff,
                )
                made_trace = plot_structure_trace(
                    bin_dir / f"{tag}.structure_trace.png", out_xyz, tag, bin_size,
                )
                if not made_trace:
                    logging.warning("no finite coordinates for %s; skipped structure trace", tag)

            finite = int(np.isfinite(D).sum())
            logging.info(
                "wrote %s matrix=%s reduce=%s finite_cells=%d",
                tag, D.shape, args.reduce, finite,
            )

    logging.info("done; wrote outputs to %s", outdir)
    return outdir


def build_parser():
    ap = argparse.ArgumentParser(
        description="Windowed physical distance maps from cryo-EM rRNA structures.",
    )
    ap.add_argument("cif", nargs="?", help="single input CIF/mmCIF file")
    ap.add_argument(
        "--accessions-csv",
        "--csv",
        dest="accessions_csv",
        help="CSV with PDB accession ID in column 1 and optional metadata/name in column 2",
    )
    ap.add_argument("--outdir", default="rrna_dist_outputs", help="base output directory")
    ap.add_argument(
        "--structures-dir",
        type=Path,
        default=default_structures_dir(),
        help="where CSV-mode downloaded CIF files are stored",
    )
    ap.add_argument("--overwrite-downloads", action="store_true", help="re-download existing CSV-mode CIFs")
    ap.add_argument("--atom", default=None, help="C1' (default), P, or 'centroid'")
    ap.add_argument("--min-len", type=int, default=100)
    ap.add_argument(
        "--bins",
        type=parse_bins,
        default=None,
        help="comma-separated window sizes in nt (default: 1,5)",
    )
    ap.add_argument(
        "--bin",
        type=int,
        action="append",
        help="legacy single-bin option; if supplied, overrides --bins",
    )
    ap.add_argument(
        "--reduce",
        choices=["centroid", "min", "contact"],
        default=None,
        help="how nucleotide-level geometry is collapsed into a window pair value "
             "(default: centroid)",
    )
    ap.add_argument(
        "--centroid-weight",
        choices=["atom", "residue"],
        default=None,
        help="weight residues by atom count (KARR-seq) or equally (default: atom)",
    )
    ap.add_argument(
        "--contact-cutoff",
        type=float,
        default=16.0,
        help="Angstrom threshold for --reduce contact; CALIBRATE THIS, see --help notes",
    )
    ap.add_argument(
        "--min-bin-modelled-fraction",
        type=float,
        default=0.0,
        help="for --reduce centroid, NaN out bins with less than this resolved/modelled fraction",
    )
    ap.add_argument(
        "--min-separation",
        type=int,
        default=0,
        help="NaN out window pairs closer than this many windows to the diagonal",
    )
    ap.add_argument(
        "--karr-seq",
        action="store_true",
        help="preset: --bins 1,5 --atom centroid --centroid-weight atom --reduce centroid",
    )
    ap.add_argument(
        "--label-from-chain",
        action="store_true",
        help="skip human-rRNA span guessing; label outputs by chain ID only",
    )
    ap.add_argument("--no-plots", action="store_true", help="skip PNG plot generation")
    ap.add_argument("--plot-max-size", type=int, default=2200, help="max plotted rows/cols")
    ap.add_argument(
        "--masked-diagonal-width",
        type=int,
        default=1,
        help="half-width, in matrix rows, for plot-only diagonal masking (default: 1)",
    )
    ap.add_argument("--verbose", action="store_true")
    return ap


def main():
    ap = build_parser()
    args = ap.parse_args()

    # Presets fill in only what the user left unset, so an explicit flag always wins.
    preset = {"bins": parse_bins("1,5"), "atom": "centroid",
              "centroid_weight": "atom", "reduce": "centroid"}
    plain = {"bins": parse_bins("1,5"), "atom": "C1'",
             "centroid_weight": "atom", "reduce": "centroid"}
    defaults = preset if args.karr_seq else plain
    overridden = []
    for key, value in defaults.items():
        if getattr(args, key) is None:
            setattr(args, key, value)
        elif args.karr_seq and getattr(args, key) != value:
            overridden.append(key.replace("_", "-"))
    if overridden:
        logging_deferred = (
            "explicit flags override the KARR-seq preset: " + ", ".join(sorted(overridden))
        )
    else:
        logging_deferred = None

    bins = sorted(set(args.bin)) if args.bin else args.bins
    if any(bin_size < 1 for bin_size in bins):
        ap.error("bin sizes must be >= 1")
    if args.reduce == "contact" and args.contact_cutoff <= 0:
        ap.error("--contact-cutoff must be > 0")
    if not 0 <= args.min_bin_modelled_fraction <= 1:
        ap.error("--min-bin-modelled-fraction must be between 0 and 1")
    if args.masked_diagonal_width < 0:
        ap.error("--masked-diagonal-width must be >= 0")

    if bool(args.cif) == bool(args.accessions_csv):
        ap.error("provide exactly one of a single CIF positional argument or --accessions-csv")

    load_runtime_deps()
    base_outdir = Path(args.outdir)

    if args.accessions_csv:
        logging.basicConfig(
            level=logging.DEBUG if args.verbose else logging.INFO,
            format="%(asctime)s %(levelname)s %(message)s",
        )
        records = read_accessions(args.accessions_csv)
        args.structures_dir.mkdir(parents=True, exist_ok=True)
        base_outdir.mkdir(parents=True, exist_ok=True)

        download_rows = [
            download_record(record, args.structures_dir, args.overwrite_downloads)
            for record in records
        ]
        structures_manifest = args.structures_dir / "manifest.tsv"
        write_manifest(structures_manifest, download_rows)

        run_rows = []
        for record, download_row in zip(records, download_rows):
            cif_path = Path(download_row["cif"])
            run_outdir = output_dir_for_run(
                base_outdir,
                cif_path,
                bins,
                args,
                metadata=record["metadata"] or download_row["metadata"],
            )
            process_cif(cif_path, run_outdir, args, bins, logging_deferred)
            run_rows.append(
                {
                    "pdb_id": record["pdb_id"],
                    "metadata": record["metadata"],
                    "cif": str(cif_path),
                    "outdir": str(run_outdir),
                    "bins": ",".join(str(bin_size) for bin_size in bins),
                    "reduce": args.reduce,
                    "atom": args.atom,
                }
            )

        run_manifest = base_outdir / f"karrseq_runs_{bins_slug(bins)}.tsv"
        with run_manifest.open("w") as fh:
            fh.write("pdb_id\tmetadata\tcif\toutdir\tbins\treduce\tatom\n")
            for row in run_rows:
                fh.write(
                    f"{row['pdb_id']}\t{row['metadata']}\t{row['cif']}\t"
                    f"{row['outdir']}\t{row['bins']}\t{row['reduce']}\t{row['atom']}\n"
                )

        print(f"\nprocessed {len(run_rows)} structures")
        print(f"download manifest: {structures_manifest}")
        print(f"run manifest: {run_manifest}")
    else:
        cif_path = Path(args.cif)
        run_outdir = output_dir_for_run(base_outdir, cif_path, bins, args)
        process_cif(cif_path, run_outdir, args, bins, logging_deferred)
        print(f"\nwrote to {run_outdir}/")

    plot_note = " and PNG plots" if not args.no_plots else ""
    print(
        "Each run folder contains per-chain .fa files; each bin folder contains "
        f".dist.npy, .xyz.npy, .bins.tsv, .lookup.tsv{plot_note}."
    )


if __name__ == "__main__":
    main()
