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
import re
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
                "%sRNA exact-sequence-per-organism representatives; %d sequence groups; entities %s"
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


RIBOSOMAL_KEYWORDS = ("ribosomal", "rrna", "ribosome", "ssu", "lsu")
RIBOSOMAL_SUBUNIT_RE = re.compile(
    r"(^|[^0-9])(5\.8s|16s|18s|23s|25s|28s|5s|30s|40s|50s|60s|70s|80s)([^a-z]|$)", re.I
)
# Deliberately no trailing \b: real-world descriptions run the S-value straight
# into the molecule name ("16SrRNA", "23SRNA"), so only the leading edge of the
# number is anchored.
RRNA_SVALUE_RE = re.compile(r"\b(\d+(?:\.\d+)?)S", re.I)
# These name the assembled ribosome particle (30S/50S -> bacterial SSU/LSU,
# 40S/60S -> eukaryotic SSU/LSU, 70S/80S -> whole ribosome), not a single RNA
# molecule, so a match here is never trusted as the molecule's own identity.
PARTICLE_SVALUES = {"30", "40", "50", "60", "70", "80"}
# "mitochondrial" alone misses real cases: depositors routinely title an
# entry "39S/55S mitoribosome" without ever writing "mitochondrial" in the
# individual chain's description (e.g. RCSB 9CN3, 8K2A).
MITO_COMPARTMENT_RE = re.compile(r"mitochondrial|mitoribosom|\bmt-", re.I)
# Pre-rRNA processing intermediates (internal/external transcribed spacers)
# are not any mature rRNA molecule -- their length can coincidentally land
# near a reference value (e.g. human ITS2 at 1167nt vs mt-12S's 954nt) and
# must never be force-matched into that bucket by the length fallback below.
SPACER_RE = re.compile(r"\bits[12]\b|external transcribed spacer|internal transcribed spacer|[35]['’]ets\b", re.I)
# Reference lengths (nt) let a known organism's own entity length overrule a
# depositor's free-text label rather than just trusting it -- catches both
# missing compartment words (label right, "(mt)" missing) and outright
# mislabels (e.g. RCSB 9RU9's human LSU rRNA entity is described "23S rRNA",
# copied from a bacterial deposition template; humans have no 23S rRNA at
# all -- its 5069nt length is human cytoplasmic 28S almost exactly).
# Scoped to Homo sapiens for now, since it is always a candidate priority
# species; extend this table if other organisms show the same failure.
ORGANISM_RRNA_REFERENCE_NT = {
    "Homo sapiens": {
        "cytoplasmic": {"18S": 1869.0, "5.8S": 156.0, "28S": 5070.0, "5S": 120.0},
        "mitochondrial": {"12S": 954.0, "16S": 1559.0},
    },
}
NON_RIBOSOMAL_TYPE_RULES = (
    ("spliceosomal RNA", lambda text: "snrna" in text or "spliceosom" in text),
    ("group I/II intron", lambda text: "group i intron" in text or "group ii intron" in text),
    ("telomerase RNA", lambda text: "telomerase" in text or "tlc1" in text),
    ("riboswitch", lambda text: "riboswitch" in text),
    ("ribozyme", lambda text: "ribozyme" in text),
    ("tRNA", lambda text: "trna" in text),
    ("mRNA", lambda text: "mrna" in text),
    ("viral/phage RNA", lambda text: any(kw in text for kw in ("phage", "viral", "virion", "virus"))),
)


def _split_organisms(organism):
    if not organism:
        return []
    return [name.strip() for name in organism.split(";") if name.strip()]


def _extract_rrna_svalue(text):
    values = RRNA_SVALUE_RE.findall(text)
    if not values:
        return None
    molecule_level = [v for v in values if v not in PARTICLE_SVALUES]
    if molecule_level:
        return "%sS" % molecule_level[0].upper(), "molecule"
    return values[0], "particle"


# How far an entity's own length may sit from a reference rRNA length and
# still be treated as that molecule (partial cryo-EM models regularly drop
# unresolved expansion segments, so some slack is needed) -- whichever of
# the two is more permissive: a flat floor for the shortest rRNAs (5S/5.8S,
# where a percentage alone would be unreasonably tight) or a percentage for
# the longer ones. Anything outside this is left unspecified rather than
# guessed: a bare "ribosomal RNA" at 300nt or 700nt is not obviously 5.8S or
# mt-12S just because those happen to be the nearest reference values.
MIN_ABSOLUTE_LENGTH_TOLERANCE_NT = 30.0
MAX_RELATIVE_LENGTH_DEVIATION = 0.20


def _nearest_reference_label(length, reference):
    if length is None or not reference:
        return None
    label, ref_length = min(reference.items(), key=lambda item: abs(item[1] - length))
    tolerance = max(MIN_ABSOLUTE_LENGTH_TOLERANCE_NT, MAX_RELATIVE_LENGTH_DEVIATION * ref_length)
    if abs(ref_length - length) > tolerance:
        return None
    return label


def _unspecified_ribosomal_label(source_text, particle_value, mito_suffix):
    if particle_value in ("30", "40"):
        return "SSU rRNA (unspecified)%s" % mito_suffix
    if particle_value in ("50", "60"):
        return "LSU rRNA (unspecified)%s" % mito_suffix
    if particle_value is not None:
        return "ribosomal RNA (unspecified)%s" % mito_suffix
    lowered = source_text.lower()
    if "small subunit" in lowered or "ssu" in lowered:
        return "SSU rRNA (unspecified)%s" % mito_suffix
    if "large subunit" in lowered or "lsu" in lowered:
        return "LSU rRNA (unspecified)%s" % mito_suffix
    return "ribosomal RNA (unspecified)%s" % mito_suffix


def _label_ribosomal(source_text, combined_lower, trust_svalue, organism, length):
    extraction = _extract_rrna_svalue(source_text) if trust_svalue else None
    svalue = extraction[0] if extraction and extraction[1] == "molecule" else None
    particle_value = extraction[0] if extraction and extraction[1] == "particle" else None

    if svalue is None and SPACER_RE.search(combined_lower):
        mito_suffix = " (mt)" if MITO_COMPARTMENT_RE.search(combined_lower) else ""
        return "pre-rRNA spacer (unspecified)%s" % mito_suffix

    # source_organisms can be a "; "-joined multi-organism string (chimeric
    # or co-purified constructs) -- match each named organism individually
    # rather than the whole joined string, which would never hit the table.
    organism_reference = None
    for name in _split_organisms(organism):
        if name in ORGANISM_RRNA_REFERENCE_NT:
            organism_reference = ORGANISM_RRNA_REFERENCE_NT[name]
            break

    if organism_reference is not None:
        # For an organism we have a reference table for, the S-value itself
        # is compartment-diagnostic (human 16S/12S are ALWAYS mitochondrial,
        # 18S/28S/5.8S/5S ALWAYS cytoplasmic) -- trust that directly rather
        # than depending on the title spelling "mitochondrial"/"mt-" in a
        # form our keyword regex happens to catch (real depositor variants
        # we've seen: "mitoribosome", "mtLSU", plain "16S rRNA" with the
        # compartment only implied by context).
        for compartment, reference in organism_reference.items():
            if svalue in reference:
                return "%s rRNA%s" % (svalue, " (mt)" if compartment == "mitochondrial" else "")

        # svalue is missing, or names something that isn't a real rRNA for
        # this organism at all (e.g. a bacterial-style "23S" on a human
        # entity, likely copy-pasted from another deposition's template).
        # This organism's rRNAs have well-separated lengths across BOTH
        # compartments together, so use the entity's own length instead of
        # trusting text at all -- this recovers the correct compartment too,
        # without needing a compartment keyword anywhere in the text.
        combined_reference = {}
        compartment_of = {}
        for compartment, reference in organism_reference.items():
            for label, ref_length in reference.items():
                combined_reference[label] = ref_length
                compartment_of[label] = compartment
        nearest = _nearest_reference_label(length, combined_reference)
        if nearest is not None:
            if svalue is not None:
                logging.warning(
                    "molecule-type override: %r is labelled %r but its length "
                    "(%.0fnt) matches %s %s rRNA -- using %s instead",
                    source_text, svalue, length, organism, nearest, nearest,
                )
            mito_suffix = " (mt)" if compartment_of[nearest] == "mitochondrial" else ""
            return "%s rRNA%s" % (nearest, mito_suffix)

        # Neither the text nor the length lands close enough to any known
        # rRNA for this organism. This organism's rRNA repertoire is fully
        # known, so an unmatched svalue (e.g. "23S" on a human entity) is
        # not trusted as a fallback label either -- it's already established
        # as invalid for this organism, and guessing would just trade one
        # wrong specific label for another.
        if svalue is not None:
            logging.warning(
                "molecule-type unresolved: %r (organism %s, %s) doesn't match any "
                "known rRNA closely enough -- leaving unspecified rather than guessing",
                source_text, organism, ("%.0fnt" % length) if length is not None else "no length",
            )
        mito_suffix = " (mt)" if MITO_COMPARTMENT_RE.search(combined_lower) else ""
        return _unspecified_ribosomal_label(source_text, particle_value, mito_suffix)

    # No reference table for this organism -- keep the older, purely
    # textual heuristic (organism-specific length validation isn't available).
    mito_suffix = " (mt)" if MITO_COMPARTMENT_RE.search(combined_lower) else ""
    if svalue is not None:
        return "%s rRNA%s" % (svalue, mito_suffix)
    return _unspecified_ribosomal_label(source_text, particle_value, mito_suffix)


def _classify_text(source_text, combined_lower, trust_svalue, organism, length):
    lowered = source_text.lower()
    is_ribosomal = (
        any(kw in lowered for kw in RIBOSOMAL_KEYWORDS)
        or bool(RIBOSOMAL_SUBUNIT_RE.search(lowered))
        or (trust_svalue and bool(RRNA_SVALUE_RE.search(source_text)))
    )
    if is_ribosomal:
        return _label_ribosomal(source_text, combined_lower, trust_svalue, organism, length)
    for label, matches in NON_RIBOSOMAL_TYPE_RULES:
        if matches(lowered):
            return label
    return None


def classify_molecule_type(description, entry_title, organism=None, sequence_length=None):
    """Classify an RNA polymer entity by molecule identity.

    The entity's own `description` is trusted first (it names that specific
    chain, and any S-value in it is taken as the molecule's own identity).
    `entry_title` describes the whole PDB entry/particle ("90S pre-ribosome",
    "80S ribosome bound to mRNA") and is only consulted when the description
    gives no signal at all -- and even then, a title S-value is never trusted
    as a specific molecule name (it almost always names the particle, not a
    single RNA chain within it), so title-derived hits land in an
    "(unspecified)" bucket rather than guessing a number.

    Mitochondrial vs cytoplasmic compartment is read from title+description
    keywords, independently of which one supplied the molecule name. For an
    organism in ORGANISM_RRNA_REFERENCE_NT, `sequence_length` is used to
    validate the text-derived label (or fill it in / correct it) against
    that organism's known rRNA lengths, since depositor free text is not
    always trustworthy (wrong S-value, or a missing "mitochondrial").
    """
    description = description or ""
    entry_title = entry_title or ""
    combined_lower = ("%s %s" % (description, entry_title)).lower()
    return (
        _classify_text(description, combined_lower, True, organism, sequence_length)
        or _classify_text(entry_title, combined_lower, False, organism, sequence_length)
        or "other/unclassified"
    )


def _row_sequence_length(row):
    for key in ("entity_sequence_length", "canonical_sequence_length"):
        value = row.get(key)
        if value not in (None, ""):
            try:
                return float(value)
            except (TypeError, ValueError):
                continue
    return None


def select_molecule_representatives(representative_rows, priority_species, other_species_count):
    by_type = {}
    ordered_types = []
    for row in representative_rows:
        mtype = classify_molecule_type(
            row.get("description"), row.get("entry_title"),
            organism=row.get("source_organisms"), sequence_length=_row_sequence_length(row),
        )
        if mtype not in by_type:
            by_type[mtype] = []
            ordered_types.append(mtype)
        by_type[mtype].append(row)

    def best_of(rows):
        return min(
            rows,
            key=lambda row: (
                first_resolution_value(row.get("entry_resolution")),
                -release_timestamp(row.get("entry_initial_release_date")),
            ),
        )

    selected = []
    for mtype in ordered_types:
        by_organism = {}
        for row in by_type[mtype]:
            organism = row.get("source_organisms") or "unknown/synthetic"
            by_organism.setdefault(organism, []).append(row)

        picks = []
        # Match against each individually-named organism, not the whole
        # (possibly "; "-joined, for chimeric/co-purified constructs) key,
        # so a joined string containing the priority species still counts.
        priority_keys = [
            key for key in by_organism if priority_species in _split_organisms(key)
        ]
        priority_rows = [row for key in priority_keys for row in by_organism.pop(key)]
        if priority_rows:
            picks.append(("priority", best_of(priority_rows)))

        other_best = sorted(
            (best_of(rows) for rows in by_organism.values()),
            key=lambda row: first_resolution_value(row.get("entry_resolution")),
        )
        for row in other_best[:other_species_count]:
            picks.append(("other", row))

        for rank, (reason, row) in enumerate(picks, start=1):
            enriched = dict(row)
            enriched.update(
                {
                    "molecule_type": mtype,
                    "selection_reason": reason,
                    "selection_rank": rank,
                }
            )
            selected.append(enriched)

    logging.info(
        "selected %d molecule-type representatives across %d molecule types",
        len(selected), len(ordered_types),
    )
    return selected


def collapse_selected_accessions(selected_rows):
    grouped = {}
    ordered_pdb_ids = []
    for row in selected_rows:
        pdb_id = row["pdb_id"]
        if pdb_id not in grouped:
            grouped[pdb_id] = []
            ordered_pdb_ids.append(pdb_id)
        grouped[pdb_id].append(row)

    records = []
    for pdb_id in ordered_pdb_ids:
        rows = grouped[pdb_id]
        entity_ids = sorted({row["entity_id"] for row in rows}, key=entity_sort_key)
        labels = sorted(
            {
                "%s (%s)" % (row["molecule_type"], row.get("source_organisms") or "unknown")
                for row in rows
            }
        )
        metadata = "curated representative; entities %s; %s" % (" ".join(entity_ids), "; ".join(labels))
        records.append({"pdb_id": pdb_id, "metadata": metadata})
    return records


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
        "selected %d exact-sequence-per-organism representatives from %d RNA entities",
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


def write_entity_manifest(path, rows, extra_fields=None):
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
    ] + list(extra_fields or [])
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
        help="deduplicate exact deposited RNA sequences per organism and output one "
        "representative entity per (organism, sequence) group -- the same exact "
        "sequence in two different organisms is kept as two representatives",
    )
    ap.add_argument(
        "--representative-rank",
        choices=["newest", "best-resolution"],
        default="newest",
        help="ranking used within exact-sequence-per-organism groups (default: newest)",
    )
    ap.add_argument(
        "--priority-species",
        default="Homo sapiens",
        help="organism guaranteed a slot (if present) in the --representatives-selected output "
        "(default: Homo sapiens); matched against RCSB's ncbi_scientific_name",
    )
    ap.add_argument(
        "--other-species-count",
        type=int,
        default=2,
        help="number of additional non-priority species kept per molecule type in the "
        "--representatives-selected output, chosen by best resolution (default: 2)",
    )
    ap.add_argument("--accessions-out", type=Path, help="downstream-compatible PDB accession CSV")
    ap.add_argument("--manifest-out", type=Path, help="entity-level TSV manifest")
    ap.add_argument("--representatives-out", type=Path, help="selected representative entity TSV")
    ap.add_argument(
        "--selected-out",
        type=Path,
        help="molecule-type-curated representative TSV (only written with --representatives)",
    )
    ap.add_argument(
        "--selected-accessions-out",
        type=Path,
        help="accession CSV for the curated set (only written with --representatives)",
    )
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
    if args.other_species_count < 0:
        ap.error("--other-species-count must be >= 0")

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
    selected_out = args.selected_out or prefix.with_suffix(".selected.tsv")
    selected_accessions_out = args.selected_accessions_out or prefix.with_suffix(".selected.accessions.csv")

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
        selected_rows = []
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

            selected_rows = select_molecule_representatives(
                representative_rows, args.priority_species, args.other_species_count
            )
            selected_records = collapse_selected_accessions(selected_rows)
            write_entity_manifest(
                selected_out, selected_rows,
                extra_fields=["molecule_type", "selection_reason", "selection_rank"],
            )
            write_accessions_csv(selected_accessions_out, selected_records)
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
        print("selected %d exact-sequence-per-organism representatives" % len(representative_rows))
        print("wrote representatives manifest: %s" % representatives_out)
        print("curated %d molecule-type representatives" % len(selected_rows))
        print("wrote selected manifest: %s" % selected_out)
        print("wrote selected accessions: %s" % selected_accessions_out)
    if download_manifest:
        print("download manifest: %s" % download_manifest)
    return 0


if __name__ == "__main__":
    sys.exit(main())
