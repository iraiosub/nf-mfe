#!/usr/bin/env python3
"""
Find PDB entries containing RNA polymer entities of at least a given length.

The main output is a downstream-compatible accession CSV:

    pdb_id,name
    8QOI,RNA >=200 nt; entities 1 2

That CSV can be passed directly to download_pdb_cifs.py or to
rrna_dist_karrseq.py --accessions-csv.
"""

import argparse
import csv
import json
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
POLYMER_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
DEFAULT_ROWS = 1000
EXPERIMENTAL_METHODS = {
    "any": None,
    "cryo-em": "EM",
}


def default_prefix(min_length, experimental_method):
    suffix = "" if experimental_method == "any" else "_cryoem"
    return Path(__file__).resolve().parent / ("rna_pdbs_min%d%s" % (min_length, suffix))


def post_json(url, payload, timeout, retries):
    data = json.dumps(payload).encode("utf-8")
    req = Request(
        url,
        data=data,
        headers={"Content-Type": "application/json", "Accept": "application/json"},
    )
    return request_json(req, timeout, retries)


def get_json(url, timeout, retries):
    req = Request(url, headers={"Accept": "application/json"})
    return request_json(req, timeout, retries)


def request_json(req, timeout, retries):
    last_exc = None
    for attempt in range(retries + 1):
        try:
            with urlopen(req, timeout=timeout) as resp:
                if resp.status == 204:
                    return None
                return json.loads(resp.read().decode("utf-8"))
        except (HTTPError, URLError, TimeoutError) as exc:
            last_exc = exc
            if attempt >= retries:
                break
            time.sleep(1.5 * (attempt + 1))
    raise RuntimeError("RCSB request failed after %d attempts: %s" % (retries + 1, last_exc))


def build_search_query(min_length, experimental_method, start, rows):
    nodes = [
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "entity_poly.rcsb_entity_polymer_type",
                "operator": "exact_match",
                "value": "RNA",
            },
        },
        {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "entity_poly.rcsb_sample_sequence_length",
                "operator": "greater_or_equal",
                "value": min_length,
            },
        },
    ]

    method_value = EXPERIMENTAL_METHODS[experimental_method]
    if method_value:
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.experimental_method",
                    "operator": "exact_match",
                    "value": method_value,
                },
            }
        )

    return {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": nodes,
        },
        "return_type": "polymer_entity",
        "request_options": {
            "paginate": {"start": start, "rows": rows},
            "results_verbosity": "minimal",
        },
    }


def parse_polymer_entity_id(identifier):
    if "_" not in identifier:
        raise ValueError("unexpected polymer_entity identifier from RCSB: %r" % identifier)
    pdb_id, entity_id = identifier.split("_", 1)
    return pdb_id.upper(), entity_id


def find_rna_entities(min_length, experimental_method, rows, limit, timeout, retries):
    hits = []
    start = 0
    total = None

    while True:
        payload = build_search_query(min_length, experimental_method, start, rows)
        response = post_json(SEARCH_URL, payload, timeout, retries)
        if response is None:
            logging.info("RCSB returned no hits")
            return []

        if total is None:
            total = int(response.get("total_count", 0))
            method_label = "any method" if experimental_method == "any" else experimental_method
            logging.info(
                "RCSB reported %d RNA polymer entities >= %d nt (%s)",
                total,
                min_length,
                method_label,
            )

        result_set = response.get("result_set", [])
        if not result_set:
            break

        for item in result_set:
            identifier = item["identifier"]
            pdb_id, entity_id = parse_polymer_entity_id(identifier)
            hits.append(
                {
                    "pdb_id": pdb_id,
                    "entity_id": entity_id,
                    "rcsb_id": identifier,
                    "score": item.get("score", ""),
                    "query_min_length": min_length,
                    "query_experimental_method": experimental_method,
                }
            )
            if limit and len(hits) >= limit:
                logging.info("stopped at --limit %d entities", limit)
                return hits

        start += len(result_set)
        logging.info("retrieved %d/%d entities", min(start, total), total)
        if start >= total:
            break

    return hits


def collapse_accessions(entity_rows, min_length, experimental_method):
    grouped = {}
    ordered_pdb_ids = []
    for row in entity_rows:
        pdb_id = row["pdb_id"]
        if pdb_id not in grouped:
            grouped[pdb_id] = []
            ordered_pdb_ids.append(pdb_id)
        grouped[pdb_id].append(row["entity_id"])

    records = []
    method_prefix = "" if experimental_method == "any" else "cryo-EM "
    for pdb_id in ordered_pdb_ids:
        entity_ids = sorted(grouped[pdb_id], key=entity_sort_key)
        records.append(
            {
                "pdb_id": pdb_id,
                "metadata": "%sRNA >=%d nt; entities %s"
                % (method_prefix, min_length, " ".join(entity_ids)),
            }
        )
    return records


def entity_sort_key(value):
    try:
        return (0, int(value))
    except ValueError:
        return (1, value)


def fetch_entity_metadata(row, timeout, retries):
    url = POLYMER_ENTITY_URL.format(pdb_id=row["pdb_id"], entity_id=row["entity_id"])
    data = get_json(url, timeout, retries)
    entity_poly = data.get("entity_poly", {}) if data else {}
    entity = data.get("rcsb_polymer_entity", {}) if data else {}
    identifiers = data.get("rcsb_polymer_entity_container_identifiers", {}) if data else {}
    organisms = data.get("rcsb_entity_source_organism", []) if data else []

    if isinstance(organisms, dict):
        organisms = [organisms]

    enriched = dict(row)
    enriched.update(
        {
            "entity_polymer_type": entity_poly.get("rcsb_entity_polymer_type", ""),
            "entity_sequence_length": entity_poly.get("rcsb_sample_sequence_length", ""),
            "description": entity.get("pdbx_description", ""),
            "auth_asym_ids": ",".join(identifiers.get("auth_asym_ids", []) or []),
            "asym_ids": ",".join(identifiers.get("asym_ids", []) or []),
            "source_organisms": "; ".join(
                sorted(
                    {
                        organism.get("ncbi_scientific_name", "")
                        for organism in organisms
                        if organism.get("ncbi_scientific_name")
                    }
                )
            ),
            "metadata_status": "ok",
        }
    )
    return enriched


def enrich_metadata(entity_rows, workers, timeout, retries):
    enriched = [None] * len(entity_rows)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(fetch_entity_metadata, row, timeout, retries): idx
            for idx, row in enumerate(entity_rows)
        }
        for done, future in enumerate(as_completed(futures), start=1):
            idx = futures[future]
            try:
                enriched[idx] = future.result()
            except Exception as exc:
                row = dict(entity_rows[idx])
                row["metadata_status"] = "failed: %s" % exc
                enriched[idx] = row
                logging.warning("metadata failed for %s: %s", row["rcsb_id"], exc)
            if done % 100 == 0 or done == len(entity_rows):
                logging.info("metadata %d/%d entities", done, len(entity_rows))
    return enriched


def write_accessions_csv(path, records):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["pdb_id", "name"])
        writer.writeheader()
        for record in records:
            writer.writerow({"pdb_id": record["pdb_id"], "name": record["metadata"]})


def write_entity_manifest(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "pdb_id",
        "entity_id",
        "rcsb_id",
        "query_min_length",
        "query_experimental_method",
        "score",
        "entity_polymer_type",
        "entity_sequence_length",
        "description",
        "auth_asym_ids",
        "asym_ids",
        "source_organisms",
        "metadata_status",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def download_accessions(records, structures_dir, overwrite):
    try:
        from download_pdb_cifs import download_record, write_manifest
    except ImportError:
        from nc_rna_benchmarking.download_pdb_cifs import download_record, write_manifest

    structures_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for idx, record in enumerate(records, start=1):
        logging.info("download %d/%d %s", idx, len(records), record["pdb_id"])
        downloader_record = {
            "pdb_id": record["pdb_id"],
            "metadata": record["metadata"],
            "extra": [],
            "source_line": idx + 1,
        }
        rows.append(download_record(downloader_record, structures_dir, overwrite))

    manifest = structures_dir / "manifest.tsv"
    write_manifest(manifest, rows)
    return manifest


def build_parser():
    ap = argparse.ArgumentParser(
        description="Find PDB entries containing RNA polymer entities at least N nt long.",
    )
    ap.add_argument("--min-rna-length", type=int, default=200, help="minimum RNA polymer length in nt")
    ap.add_argument(
        "--experimental-method",
        choices=sorted(EXPERIMENTAL_METHODS),
        default="any",
        help="restrict by experimental method; currently supports 'cryo-em' or 'any'",
    )
    ap.add_argument(
        "--cryo-em-only",
        action="store_true",
        help="shortcut for --experimental-method cryo-em",
    )
    ap.add_argument("--accessions-out", type=Path, help="downstream-compatible PDB accession CSV")
    ap.add_argument("--manifest-out", type=Path, help="entity-level TSV manifest")
    ap.add_argument(
        "--out-prefix",
        type=Path,
        help="prefix for default outputs: PREFIX.accessions.csv and PREFIX.entities.tsv",
    )
    ap.add_argument("--rows", type=int, default=DEFAULT_ROWS, help="RCSB page size")
    ap.add_argument("--limit", type=int, help="stop after this many RNA entities, for smoke tests")
    ap.add_argument(
        "--fetch-metadata",
        action="store_true",
        help="enrich the TSV with exact entity lengths, descriptions, chains, and organisms",
    )
    ap.add_argument("--metadata-workers", type=int, default=4, help="parallel metadata requests")
    ap.add_argument("--timeout", type=float, default=30.0, help="network timeout in seconds")
    ap.add_argument("--retries", type=int, default=2, help="retry failed RCSB requests")
    ap.add_argument(
        "--download",
        action="store_true",
        help="also download the unique PDB/mmCIF files after writing the accession CSV",
    )
    ap.add_argument(
        "--structures-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "pdb_structures",
        help="download destination used with --download",
    )
    ap.add_argument("--overwrite", action="store_true", help="re-download existing CIFs with --download")
    ap.add_argument("--verbose", action="store_true")
    return ap


def main():
    ap = build_parser()
    args = ap.parse_args()

    if args.min_rna_length < 1:
        ap.error("--min-rna-length must be >= 1")
    if args.cryo_em_only:
        args.experimental_method = "cryo-em"
    if args.rows < 1:
        ap.error("--rows must be >= 1")
    if args.limit is not None and args.limit < 1:
        ap.error("--limit must be >= 1")
    if args.metadata_workers < 1:
        ap.error("--metadata-workers must be >= 1")

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    prefix = args.out_prefix or default_prefix(args.min_rna_length, args.experimental_method)
    accessions_out = args.accessions_out or prefix.with_suffix(".accessions.csv")
    manifest_out = args.manifest_out or prefix.with_suffix(".entities.tsv")

    try:
        entity_rows = find_rna_entities(
            args.min_rna_length,
            args.experimental_method,
            args.rows,
            args.limit,
            args.timeout,
            args.retries,
        )
        if args.fetch_metadata and entity_rows:
            entity_rows = enrich_metadata(
                entity_rows,
                args.metadata_workers,
                args.timeout,
                args.retries,
            )

        records = collapse_accessions(entity_rows, args.min_rna_length, args.experimental_method)
        write_accessions_csv(accessions_out, records)
        write_entity_manifest(manifest_out, entity_rows)

        download_manifest = None
        if args.download:
            download_manifest = download_accessions(records, args.structures_dir, args.overwrite)
    except Exception as exc:
        logging.error("%s", exc)
        return 1

    print("found %d RNA polymer entities" % len(entity_rows))
    print("collapsed to %d unique PDB accessions" % len(records))
    print("wrote downstream CSV: %s" % accessions_out)
    print("wrote entity manifest: %s" % manifest_out)
    if download_manifest:
        print("download manifest: %s" % download_manifest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
