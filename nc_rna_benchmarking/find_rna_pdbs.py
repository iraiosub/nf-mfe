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
import hashlib
import json
import logging
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import lru_cache
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
POLYMER_ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}"
ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
DEFAULT_ROWS = 1000
DEFAULT_METADATA_BATCH_SIZE = 100
EXPERIMENTAL_METHODS = {
    "any": None,
    "cryo-em": "ELECTRON MICROSCOPY",
}
GRAPHQL_METADATA_QUERY = """
query($ids: [String!]!) {
  polymer_entities(entity_ids: $ids) {
    rcsb_id
    entity_poly {
      rcsb_entity_polymer_type
      rcsb_sample_sequence_length
      pdbx_seq_one_letter_code_can
      pdbx_seq_one_letter_code
    }
    rcsb_polymer_entity { pdbx_description }
    rcsb_polymer_entity_container_identifiers {
      entry_id
      entity_id
      auth_asym_ids
      asym_ids
    }
    rcsb_entity_source_organism { ncbi_scientific_name }
    entry {
      struct { title }
      rcsb_accession_info { initial_release_date }
      rcsb_entry_info { experimental_method resolution_combined }
      exptl { method }
    }
  }
}
"""


def default_prefix(min_length, experimental_method, max_resolution=None, representatives=False):
    suffix = "" if experimental_method == "any" else "_cryoem"
    if max_resolution is not None:
        resolution = ("%g" % max_resolution).replace(".", "p")
        suffix += "_resle%sA" % resolution
    if representatives:
        suffix += "_representatives"
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


def build_search_query(min_length, experimental_method, max_resolution, start, rows):
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
                    "attribute": "exptl.method",
                    "operator": "exact_match",
                    "value": method_value,
                },
            }
        )
    if max_resolution is not None:
        nodes.append(
            {
                "type": "terminal",
                "service": "text",
                "parameters": {
                    "attribute": "rcsb_entry_info.resolution_combined",
                    "operator": "less_or_equal",
                    "value": max_resolution,
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


def find_rna_entities(min_length, experimental_method, max_resolution, rows, limit, timeout, retries):
    hits = []
    start = 0
    total = None

    while True:
        payload = build_search_query(min_length, experimental_method, max_resolution, start, rows)
        response = post_json(SEARCH_URL, payload, timeout, retries)
        if response is None:
            logging.info("RCSB returned no hits")
            return []

        if total is None:
            total = int(response.get("total_count", 0))
            method_label = "any method" if experimental_method == "any" else experimental_method
            resolution_label = (
                "any resolution" if max_resolution is None else "resolution <= %.3g A" % max_resolution
            )
            logging.info(
                "RCSB reported %d RNA polymer entities >= %d nt (%s, %s)",
                total, min_length, method_label, resolution_label,
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
                    "query_max_resolution": "" if max_resolution is None else max_resolution,
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


def collapse_accessions(entity_rows, min_length, experimental_method, representatives=False):
    grouped = {}
    ordered_pdb_ids = []
    for row in entity_rows:
        pdb_id = row["pdb_id"]
        if pdb_id not in grouped:
            grouped[pdb_id] = []
            ordered_pdb_ids.append(pdb_id)
        grouped[pdb_id].append(row)

    records = []
    method_prefix = "" if experimental_method == "any" else "cryo-EM "
    for pdb_id in ordered_pdb_ids:
        rows = grouped[pdb_id]
        entity_ids = sorted([row["entity_id"] for row in rows], key=entity_sort_key)
        if representatives:
            cluster_count = len({row.get("sequence_cluster_id", row["rcsb_id"]) for row in rows})
            metadata = (
                "%sRNA exact-sequence representatives; %d sequence groups; entities %s"
                % (method_prefix, cluster_count, " ".join(entity_ids))
            )
        else:
            metadata = "%sRNA >=%d nt; entities %s" % (
                method_prefix, min_length, " ".join(entity_ids)
            )
        records.append(
            {
                "pdb_id": pdb_id,
                "metadata": metadata,
            }
        )
    return records


def entity_sort_key(value):
    try:
        return (0, int(value))
    except ValueError:
        return (1, value)


def clean_sequence(value):
    if not value:
        return ""
    return "".join(str(value).split()).upper().replace("T", "U")


def sequence_digest(sequence):
    if not sequence:
        return ""
    return hashlib.sha1(sequence.encode("ascii", "ignore")).hexdigest()


def first_resolution_value(value):
    if value in (None, ""):
        return float("inf")
    if isinstance(value, (int, float)):
        return float(value)
    values = []
    for part in str(value).split(","):
        try:
            values.append(float(part))
        except ValueError:
            continue
    return min(values) if values else float("inf")


def release_timestamp(value):
    if not value:
        return 0.0
    from datetime import datetime

    try:
        return datetime.fromisoformat(str(value).replace("Z", "+00:00")).timestamp()
    except ValueError:
        return 0.0


@lru_cache(maxsize=None)
def fetch_entry_metadata(pdb_id, timeout, retries):
    data = get_json(ENTRY_URL.format(pdb_id=pdb_id), timeout, retries) or {}
    methods = [
        item.get("method", "")
        for item in data.get("exptl", []) or []
        if item.get("method")
    ]
    entry_info = data.get("rcsb_entry_info") or {}
    accession_info = data.get("rcsb_accession_info") or {}
    struct_info = data.get("struct") or {}
    resolutions = entry_info.get("resolution_combined", []) or []
    return {
        "entry_title": struct_info.get("title", ""),
        "entry_initial_release_date": accession_info.get("initial_release_date", ""),
        "entry_experimental_method_summary": entry_info.get("experimental_method", ""),
        "entry_experimental_methods": "; ".join(methods),
        "entry_resolution": ",".join("%g" % value for value in resolutions),
        "cryo_em_method_present": int("ELECTRON MICROSCOPY" in methods),
    }


def entity_metadata_from_graphql(row, entity):
    entity_poly = entity.get("entity_poly", {}) if entity else {}
    entity_poly = entity_poly or {}
    entity_info = entity.get("rcsb_polymer_entity", {}) if entity else {}
    entity_info = entity_info or {}
    identifiers = entity.get("rcsb_polymer_entity_container_identifiers", {}) if entity else {}
    identifiers = identifiers or {}
    organisms = entity.get("rcsb_entity_source_organism", []) if entity else []
    organisms = organisms or []
    entry = entity.get("entry", {}) if entity else {}
    entry = entry or {}

    if isinstance(organisms, dict):
        organisms = [organisms]

    methods = [
        item.get("method", "")
        for item in entry.get("exptl", []) or []
        if item.get("method")
    ]
    entry_info = entry.get("rcsb_entry_info") or {}
    accession_info = entry.get("rcsb_accession_info") or {}
    struct_info = entry.get("struct") or {}
    resolutions = entry_info.get("resolution_combined", []) or []
    canonical_sequence = clean_sequence(
        entity_poly.get("pdbx_seq_one_letter_code_can")
        or entity_poly.get("pdbx_seq_one_letter_code")
        or ""
    )

    enriched = dict(row)
    enriched.update(
        {
            "entry_title": struct_info.get("title", ""),
            "entry_initial_release_date": accession_info.get("initial_release_date", ""),
            "entry_experimental_method_summary": entry_info.get("experimental_method", ""),
            "entry_experimental_methods": "; ".join(methods),
            "entry_resolution": ",".join("%g" % value for value in resolutions),
            "cryo_em_method_present": int("ELECTRON MICROSCOPY" in methods),
            "entity_polymer_type": entity_poly.get("rcsb_entity_polymer_type", ""),
            "entity_sequence_length": entity_poly.get("rcsb_sample_sequence_length", ""),
            "canonical_sequence_length": len(canonical_sequence),
            "sequence_hash": sequence_digest(canonical_sequence),
            "description": entity_info.get("pdbx_description", ""),
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


def chunked(rows, size):
    for start in range(0, len(rows), size):
        yield start, rows[start:start + size]


def fetch_metadata_batch(batch, timeout, retries):
    payload = {
        "query": GRAPHQL_METADATA_QUERY,
        "variables": {"ids": [row["rcsb_id"] for row in batch]},
    }
    response = post_json(GRAPHQL_URL, payload, timeout, retries)
    if response.get("errors"):
        raise RuntimeError("RCSB GraphQL metadata query failed: %s" % response["errors"])

    entities = response.get("data", {}).get("polymer_entities", []) or []
    by_id = {entity["rcsb_id"]: entity for entity in entities if entity}
    out = []
    for row in batch:
        entity = by_id.get(row["rcsb_id"])
        if not entity:
            failed = dict(row)
            failed["metadata_status"] = "missing from GraphQL response"
            out.append(failed)
            continue
        out.append(entity_metadata_from_graphql(row, entity))
    return out


def fetch_entity_metadata(row, timeout, retries):
    url = POLYMER_ENTITY_URL.format(pdb_id=row["pdb_id"], entity_id=row["entity_id"])
    data = get_json(url, timeout, retries)
    entity_poly = data.get("entity_poly", {}) if data else {}
    entity_poly = entity_poly or {}
    entity = data.get("rcsb_polymer_entity", {}) if data else {}
    entity = entity or {}
    identifiers = data.get("rcsb_polymer_entity_container_identifiers", {}) if data else {}
    identifiers = identifiers or {}
    organisms = data.get("rcsb_entity_source_organism", []) if data else []
    organisms = organisms or []

    if isinstance(organisms, dict):
        organisms = [organisms]

    canonical_sequence = clean_sequence(
        entity_poly.get("pdbx_seq_one_letter_code_can")
        or entity_poly.get("pdbx_seq_one_letter_code")
        or ""
    )
    enriched = dict(row)
    enriched.update(fetch_entry_metadata(row["pdb_id"], timeout, retries))
    enriched.update(
        {
            "entity_polymer_type": entity_poly.get("rcsb_entity_polymer_type", ""),
            "entity_sequence_length": entity_poly.get("rcsb_sample_sequence_length", ""),
            "canonical_sequence_length": len(canonical_sequence),
            "sequence_hash": sequence_digest(canonical_sequence),
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


def enrich_metadata(entity_rows, workers, timeout, retries, batch_size=DEFAULT_METADATA_BATCH_SIZE):
    if not entity_rows:
        return []
    if batch_size > 1:
        return enrich_metadata_batched(entity_rows, workers, timeout, retries, batch_size)

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


def enrich_metadata_batched(entity_rows, workers, timeout, retries, batch_size):
    enriched = [None] * len(entity_rows)
    batches = list(chunked(entity_rows, batch_size))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {
            pool.submit(fetch_metadata_batch, batch, timeout, retries): (start, batch)
            for start, batch in batches
        }
        for done, future in enumerate(as_completed(futures), start=1):
            start, batch = futures[future]
            try:
                batch_rows = future.result()
            except Exception as exc:
                logging.warning(
                    "metadata batch failed for rows %d-%d: %s",
                    start + 1,
                    start + len(batch),
                    exc,
                )
                batch_rows = []
                for row in batch:
                    failed = dict(row)
                    failed["metadata_status"] = "failed: %s" % exc
                    batch_rows.append(failed)
            for offset, row in enumerate(batch_rows):
                enriched[start + offset] = row
            logging.info("metadata batches %d/%d", done, len(batches))
    return enriched


def representative_sort_key(row, ranking):
    timestamp = release_timestamp(row.get("entry_initial_release_date", ""))
    resolution = first_resolution_value(row.get("entry_resolution", ""))
    entity_key = entity_sort_key(str(row.get("entity_id", "")))
    if ranking == "best-resolution":
        return (resolution, -timestamp, row.get("pdb_id", ""), entity_key)
    return (-timestamp, resolution, row.get("pdb_id", ""), entity_key)


def sequence_cluster_key(row):
    organism = row.get("source_organisms", "") or "unknown organism"
    digest = row.get("sequence_hash", "")
    if digest:
        return "%s|sha1:%s" % (organism, digest)
    return "%s|entity:%s" % (organism, row["rcsb_id"])


def select_sequence_representatives(entity_rows, ranking):
    groups = {}
    ordered_keys = []
    for row in entity_rows:
        key = sequence_cluster_key(row)
        if key not in groups:
            groups[key] = []
            ordered_keys.append(key)
        groups[key].append(row)

    selected = []
    annotated = []
    for cluster_index, key in enumerate(ordered_keys, start=1):
        members = sorted(groups[key], key=lambda row: representative_sort_key(row, ranking))
        cluster_id = "seq%06d" % cluster_index
        cluster_size = len(members)
        winner_id = members[0]["rcsb_id"]
        for rank, row in enumerate(members, start=1):
            enriched = dict(row)
            enriched.update(
                {
                    "sequence_cluster_id": cluster_id,
                    "sequence_cluster_key": key,
                    "sequence_cluster_size": cluster_size,
                    "representative_rank": rank,
                    "selection_status": (
                        "selected_representative"
                        if row["rcsb_id"] == winner_id
                        else "skipped_duplicate_sequence"
                    ),
                    "representative_ranking": ranking,
                }
            )
            annotated.append(enriched)
            if row["rcsb_id"] == winner_id:
                selected.append(enriched)

    logging.info(
        "selected %d exact-sequence representatives from %d RNA entities",
        len(selected), len(entity_rows),
    )
    return annotated, selected


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
        "query_max_resolution",
        "score",
        "entry_title",
        "entry_initial_release_date",
        "entry_experimental_method_summary",
        "entry_experimental_methods",
        "entry_resolution",
        "cryo_em_method_present",
        "entity_polymer_type",
        "entity_sequence_length",
        "canonical_sequence_length",
        "sequence_hash",
        "sequence_cluster_id",
        "sequence_cluster_key",
        "sequence_cluster_size",
        "representative_rank",
        "representative_ranking",
        "selection_status",
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
    ap.add_argument(
        "--max-resolution",
        type=float,
        help="optional maximum entry resolution in Angstroms, using RCSB resolution_combined",
    )
    ap.add_argument(
        "--representatives",
        action="store_true",
        help="deduplicate exact deposited RNA sequences and output one representative entity per sequence",
    )
    ap.add_argument(
        "--representative-rank",
        choices=["newest", "best-resolution"],
        default="newest",
        help="ranking used within exact-sequence groups (default: newest)",
    )
    ap.add_argument("--accessions-out", type=Path, help="downstream-compatible PDB accession CSV")
    ap.add_argument("--manifest-out", type=Path, help="entity-level TSV manifest")
    ap.add_argument("--representatives-out", type=Path, help="selected representative entity TSV")
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
        help="enrich the TSV with entity details plus entry methods/resolution",
    )
    ap.add_argument("--metadata-workers", type=int, default=4, help="parallel metadata requests")
    ap.add_argument("--metadata-batch-size", type=int, default=DEFAULT_METADATA_BATCH_SIZE)
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
    if args.max_resolution is not None and args.max_resolution <= 0:
        ap.error("--max-resolution must be > 0")
    if args.rows < 1:
        ap.error("--rows must be >= 1")
    if args.limit is not None and args.limit < 1:
        ap.error("--limit must be >= 1")
    if args.metadata_workers < 1:
        ap.error("--metadata-workers must be >= 1")
    if args.metadata_batch_size < 1:
        ap.error("--metadata-batch-size must be >= 1")

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    prefix = args.out_prefix or default_prefix(
        args.min_rna_length,
        args.experimental_method,
        args.max_resolution,
        representatives=args.representatives,
    )
    accessions_out = args.accessions_out or prefix.with_suffix(".accessions.csv")
    manifest_out = args.manifest_out or prefix.with_suffix(".entities.tsv")
    representatives_out = args.representatives_out or prefix.with_suffix(".representatives.tsv")

    try:
        entity_rows = find_rna_entities(
            args.min_rna_length,
            args.experimental_method,
            args.max_resolution,
            args.rows,
            args.limit,
            args.timeout,
            args.retries,
        )
        if args.representatives and not args.fetch_metadata:
            logging.info("--representatives requires sequence metadata; fetching metadata")
        if (args.fetch_metadata or args.representatives) and entity_rows:
            entity_rows = enrich_metadata(
                entity_rows,
                args.metadata_workers,
                args.timeout,
                args.retries,
                args.metadata_batch_size,
            )

        representative_rows = []
        if args.representatives:
            entity_rows, representative_rows = select_sequence_representatives(
                entity_rows, args.representative_rank
            )
            records = collapse_accessions(
                representative_rows,
                args.min_rna_length,
                args.experimental_method,
                representatives=True,
            )
            write_entity_manifest(representatives_out, representative_rows)
        else:
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
    if args.representatives:
        print("selected %d exact-sequence representatives" % len(representative_rows))
        print("wrote representatives manifest: %s" % representatives_out)
    if download_manifest:
        print("download manifest: %s" % download_manifest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
