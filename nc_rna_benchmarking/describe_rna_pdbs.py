#!/usr/bin/env python3
"""
Descriptive summary of an entity/representatives TSV from find_rna_pdbs.py,
meant to be run before download_pdb_cifs.py to sanity-check what a query
actually pulled in: which organisms, which kinds of RNA molecule, and the
length/resolution/release-date spread of the cryo-EM (or mixed-method) hits.

    python describe_rna_pdbs.py rna_pdbs_min800_cryoem_representatives.representatives.tsv
"""

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

try:
    from find_rna_pdbs import classify_molecule_type, first_resolution_value
except ImportError:
    from nc_rna_benchmarking.find_rna_pdbs import classify_molecule_type, first_resolution_value


def read_rows(path):
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def organism_counts(rows):
    counts = Counter()
    for row in rows:
        organisms = [o.strip() for o in (row.get("source_organisms") or "").split(";") if o.strip()]
        if not organisms:
            organisms = ["unknown/synthetic"]
        counts.update(organisms)
    return counts


def sequence_length(row):
    for key in ("entity_sequence_length", "canonical_sequence_length"):
        value = row.get(key)
        if value not in (None, ""):
            try:
                return float(value)
            except ValueError:
                continue
    return None


def release_year(row):
    value = row.get("entry_initial_release_date")
    if not value or len(value) < 4:
        return None
    try:
        return int(value[:4])
    except ValueError:
        return None


def bar_panel(ax, counts, title, top_n=None, other_label="other"):
    items = counts.most_common(top_n) if top_n else counts.most_common()
    if top_n and len(counts) > top_n:
        rest = sum(n for _, n in counts.most_common()[top_n:])
        items.append((other_label, rest))
    items = items[::-1]
    labels = [label for label, _ in items]
    values = [n for _, n in items]
    ax.barh(labels, values, color="#4C72B0")
    ax.set_title(title)
    ax.tick_params(axis="y", labelsize=8)
    for y, value in enumerate(values):
        ax.text(value, y, " %d" % value, va="center", fontsize=7)


def build_figure(rows, out_path, title, top_organisms):
    import matplotlib

    matplotlib.use("Agg")
    import numpy as np
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=(16, 9), constrained_layout=True)

    bar_panel(axes[0, 0], organism_counts(rows), "top organisms", top_n=top_organisms)

    type_counts = Counter(
        classify_molecule_type(
            row.get("description"), row.get("entry_title"),
            organism=row.get("source_organisms"), sequence_length=sequence_length(row),
        )
        for row in rows
    )
    bar_panel(axes[0, 1], type_counts, "molecule type")

    lengths = [v for v in (sequence_length(row) for row in rows) if v is not None]
    if lengths:
        axes[0, 2].hist(lengths, bins=40, color="#55A868")
        axes[0, 2].axvline(float(np.median(lengths)), color="black", linestyle="--", linewidth=1)
    axes[0, 2].set_title("sequence length (nt)")
    axes[0, 2].set_xlabel("nt")

    resolutions = [
        r for r in (first_resolution_value(row.get("entry_resolution")) for row in rows) if r != float("inf")
    ]
    if resolutions:
        axes[1, 0].hist(resolutions, bins=30, color="#C44E52")
        axes[1, 0].axvline(float(np.median(resolutions)), color="black", linestyle="--", linewidth=1)
    axes[1, 0].set_title("resolution (%d entries)" % len(resolutions))
    axes[1, 0].set_xlabel("Angstrom")

    years = Counter(y for y in (release_year(row) for row in rows) if y is not None)
    if years:
        sorted_years = sorted(years)
        axes[1, 1].bar(sorted_years, [years[y] for y in sorted_years], color="#8172B2")
        axes[1, 1].tick_params(axis="x", rotation=45)
    axes[1, 1].set_title("initial release year")

    method_counts = Counter()
    for row in rows:
        methods = [m.strip() for m in (row.get("entry_experimental_methods") or "").split(";") if m.strip()]
        method_counts.update(methods or ["unknown"])
    bar_panel(axes[1, 2], method_counts, "experimental method")

    fig.suptitle(title)
    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def build_selected_figure(rows, out_path, title):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    by_type = {}
    for row in rows:
        by_type.setdefault(row.get("molecule_type", "unknown"), []).append(row)
    ordered_types = sorted(by_type)

    fig_height = max(4.0, 0.7 * len(ordered_types) + 1.5)
    fig, ax = plt.subplots(figsize=(13, fig_height), constrained_layout=True)
    colors = {"priority": "#C44E52", "other": "#4C72B0"}
    row_span = 0.34

    for y, mtype in enumerate(ordered_types):
        entries = sorted(by_type[mtype], key=lambda row: int(row.get("selection_rank") or 0))
        n = len(entries)
        for i, row in enumerate(entries):
            resolution = first_resolution_value(row.get("entry_resolution"))
            if resolution == float("inf"):
                continue
            reason = row.get("selection_reason", "other")
            # Spread same-row points/labels vertically so close resolutions
            # don't stack their annotations on top of each other.
            yy = y + (0.0 if n == 1 else (i / (n - 1) - 0.5) * row_span)
            ax.scatter(resolution, yy, color=colors.get(reason, "#55A868"), s=60, zorder=3)
            label = "%s / %s / %.2gA" % (
                row.get("pdb_id"), row.get("source_organisms") or "?", resolution,
            )
            ax.annotate(
                label, (resolution, yy), textcoords="offset points", xytext=(8, 0),
                va="center", fontsize=7,
            )

    ax.set_yticks(range(len(ordered_types)))
    ax.set_yticklabels(ordered_types, fontsize=8)
    ax.set_ylim(-1, len(ordered_types))
    for y in range(len(ordered_types)):
        ax.axhline(y, color="0.9", linewidth=0.6, zorder=0)
    ax.set_xlabel("resolution (Angstrom, lower = better)")
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    ax.set_title(title)

    handles = [
        plt.Line2D([0], [0], marker="o", linestyle="", color=colors["priority"], label="priority species"),
        plt.Line2D([0], [0], marker="o", linestyle="", color=colors["other"], label="other species (best resolution)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8)

    fig.savefig(out_path, dpi=200)
    plt.close(fig)


def log_selected_summary(rows):
    by_type = {}
    for row in rows:
        by_type.setdefault(row.get("molecule_type", "unknown"), []).append(row)
    logging.info("%d curated rows across %d molecule types", len(rows), len(by_type))
    for mtype in sorted(by_type):
        parts = [
            "%s(%s,%.3gA,%s)" % (
                row.get("pdb_id"),
                row.get("source_organisms") or "?",
                first_resolution_value(row.get("entry_resolution")),
                row.get("selection_reason"),
            )
            for row in sorted(by_type[mtype], key=lambda row: int(row.get("selection_rank") or 0))
        ]
        logging.info("%s: %s", mtype, ", ".join(parts))


def log_summary(rows):
    n = len(rows)
    pdb_ids = {row.get("pdb_id") for row in rows}
    organisms = organism_counts(rows)
    type_counts = Counter(
        classify_molecule_type(
            row.get("description"), row.get("entry_title"),
            organism=row.get("source_organisms"), sequence_length=sequence_length(row),
        )
        for row in rows
    )
    lengths = [v for v in (sequence_length(row) for row in rows) if v is not None]
    resolutions = [
        r for r in (first_resolution_value(row.get("entry_resolution")) for row in rows) if r != float("inf")
    ]
    years = [y for y in (release_year(row) for row in rows) if y is not None]
    cryo_em = sum(1 for row in rows if row.get("cryo_em_method_present") == "1")

    logging.info("%d rows, %d unique PDB entries, %d unique organisms", n, len(pdb_ids), len(organisms))
    logging.info("molecule types: %s", dict(type_counts.most_common()))
    logging.info("top organisms: %s", organisms.most_common(5))
    if lengths:
        logging.info(
            "length (nt): median=%.0f min=%.0f max=%.0f", sorted(lengths)[len(lengths) // 2], min(lengths), max(lengths)
        )
    if resolutions:
        logging.info(
            "resolution (A): median=%.2f min=%.2f max=%.2f",
            sorted(resolutions)[len(resolutions) // 2], min(resolutions), max(resolutions),
        )
    if years:
        logging.info("release years: %d-%d", min(years), max(years))
    logging.info("cryo-EM method present: %d/%d", cryo_em, n)


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("manifest", type=Path, help="entities.tsv or representatives.tsv from find_rna_pdbs.py")
    ap.add_argument("--out", type=Path, help="output PNG path (default: <manifest-stem>.summary.png)")
    ap.add_argument("--top-organisms", type=int, default=15)
    ap.add_argument("--title", help="figure title (default: manifest filename)")
    ap.add_argument("--verbose", action="store_true")
    return ap


def main():
    ap = build_parser()
    args = ap.parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    rows = read_rows(args.manifest)
    if not rows:
        logging.error("no rows in %s", args.manifest)
        return 1

    out_path = args.out or args.manifest.with_name(args.manifest.stem + ".summary.png")
    title = args.title or args.manifest.name

    if "molecule_type" in rows[0]:
        log_selected_summary(rows)
        build_selected_figure(rows, out_path, title)
    else:
        log_summary(rows)
        build_figure(rows, out_path, title, args.top_organisms)
    print("wrote summary figure: %s" % out_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
