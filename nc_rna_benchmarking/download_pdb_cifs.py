#!/usr/bin/env python3
"""
Download PDB/mmCIF files listed in a CSV.

CSV format
----------
Header is optional. The first column is the PDB accession ID and the optional
second column is a human-readable label/metadata field.

Examples:

    pdb_id,name
    8QOI,human 80S ribosome
    9Q16,TERC / human telomerase RNA
    6QX9,U1 snRNA

or headerless:

    8QOI,human 80S ribosome
    9Q16,TERC / human telomerase RNA
    6QX9,U1 snRNA
"""

import argparse
import csv
import json
import logging
import re
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import urlretrieve


PDB_ID = re.compile(r"^[A-Za-z0-9]{4}$")
HEADER_NAMES = {"accession", "accession_id", "accesion", "pdb", "pdb_id", "id"}


def default_outdir():
    return Path(__file__).resolve().parent / "pdb_structures"


def safe_name(value):
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip())
    return cleaned.strip("._-") or "metadata"


def read_accessions(csv_path):
    records = []
    seen = set()

    with open(csv_path, newline="") as fh:
        reader = csv.reader(fh)
        for line_no, row in enumerate(reader, start=1):
            row = [cell.strip() for cell in row]
            if not row or not row[0] or row[0].startswith("#"):
                continue

            first = row[0].lower()
            if line_no == 1 and first in HEADER_NAMES:
                continue

            pdb_id = row[0].upper()
            if not PDB_ID.match(pdb_id):
                raise ValueError(f"{csv_path}:{line_no}: invalid PDB accession {row[0]!r}")
            if pdb_id in seen:
                logging.info("skipping duplicate accession %s from line %d", pdb_id, line_no)
                continue

            seen.add(pdb_id)
            records.append(
                {
                    "pdb_id": pdb_id,
                    "metadata": row[1] if len(row) > 1 else "",
                    "extra": row[2:] if len(row) > 2 else [],
                    "source_line": line_no,
                }
            )

    if not records:
        raise ValueError(f"no PDB accessions found in {csv_path}")
    return records


def download_record(record, outdir, overwrite=False):
    pdb_id = record["pdb_id"]
    label = safe_name(record["metadata"]) if record["metadata"] else pdb_id
    entry_dir = outdir / f"{pdb_id}_{label}"
    entry_dir.mkdir(parents=True, exist_ok=True)

    cif_path = entry_dir / f"{pdb_id}.cif"
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    if cif_path.exists() and not overwrite:
        logging.info("kept existing %s", cif_path)
        status = "exists"
    else:
        logging.info("downloading %s -> %s", url, cif_path)
        try:
            urlretrieve(url, cif_path)
        except (HTTPError, URLError) as exc:
            raise RuntimeError(f"failed to download {pdb_id} from {url}: {exc}") from exc
        status = "downloaded"

    metadata_path = entry_dir / "metadata.json"
    metadata_path.write_text(json.dumps(record | {"url": url, "cif": str(cif_path)}, indent=2) + "\n")

    return {
        "pdb_id": pdb_id,
        "metadata": record["metadata"],
        "status": status,
        "directory": str(entry_dir),
        "cif": str(cif_path),
        "url": url,
    }


def write_manifest(path, rows):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["pdb_id", "metadata", "status", "directory", "cif", "url"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def main():
    ap = argparse.ArgumentParser(description="Download PDB mmCIF files from a CSV of accessions.")
    ap.add_argument("csv", help="CSV with accession ID in column 1 and optional metadata/name in column 2")
    ap.add_argument("--outdir", type=Path, default=default_outdir(), help="output directory")
    ap.add_argument("--overwrite", action="store_true", help="re-download existing CIF files")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    try:
        records = read_accessions(args.csv)
        args.outdir.mkdir(parents=True, exist_ok=True)
        manifest_rows = [download_record(record, args.outdir, args.overwrite) for record in records]
        manifest_path = args.outdir / "manifest.tsv"
        write_manifest(manifest_path, manifest_rows)
    except Exception as exc:
        logging.error("%s", exc)
        return 1

    print(f"downloaded/checked {len(manifest_rows)} structures")
    print(f"wrote manifest: {manifest_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
