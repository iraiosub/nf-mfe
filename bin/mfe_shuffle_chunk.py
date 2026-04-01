#!/usr/bin/env python3

import argparse
import csv
import subprocess
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
import RNA


def normalize_sequence(sequence):
    if pd.isna(sequence):
        return None
    sequence = str(sequence).strip().upper().replace("T", "U")
    return sequence if sequence else None


def compute_mfe(pair):
    """
    Compute RNA-RNA duplex MFE for one sequence pair.
    Expects a tuple: (lseq, rseq)
    Returns: (mfe, structure)
    """
    lseq, rseq = pair

    lseq = normalize_sequence(lseq)
    rseq = normalize_sequence(rseq)

    if not lseq or not rseq:
        return (None, None)

    try:
        duplex = RNA.duplexfold(lseq, rseq)
        return (duplex.energy, duplex.structure)
    except Exception:
        return (None, None)


def reverse_sequence(sequence):
    """
    Reverse sequence string only. Does not complement.
    """
    sequence = normalize_sequence(sequence)
    if not sequence:
        return None
    return sequence[::-1]


def shuffle_sequence(sequence, number=100, klet=2, seed=42):
    """
    Shuffle sequence using ushuffle.
    Returns a list of shuffled sequences.
    """
    if pd.isna(sequence):
        return []

    sequence = str(sequence).strip().upper()
    if not sequence:
        return []

    cmd = [
        "ushuffle",
        "-seed", str(seed),
        "-k", str(klet),
        "-n", str(number),
        "-s", sequence
    ]

    try:
        result = subprocess.run(
            cmd,
            check=True,
            capture_output=True,
            text=True
        )
        shuffles = [line.strip() for line in result.stdout.splitlines() if line.strip()]
        return shuffles
    except Exception:
        return []


def summarize_null(obs_mfe, null_mfes):
    """
    Compute summary statistics for null distribution.
    Lower/more negative MFE = stronger binding.
    empirical_p_lower = proportion of null MFEs <= observed MFE.
    """
    valid = [x for x in null_mfes if x is not None]

    if len(valid) == 0:
        return {
            "mean_shuffled_mfe": None,
            "sd_shuffled_mfe": None,
            "delta_mfe": None,
            "zscore_mfe": None,
            "empirical_p_lower": None,
            "n_shuffles_ok": 0
        }

    arr = np.array(valid, dtype=float)
    mean_shuffled = float(arr.mean())
    sd_shuffled = float(arr.std(ddof=1)) if len(arr) > 1 else 0.0

    delta_mfe = None if obs_mfe is None else float(obs_mfe - mean_shuffled)

    if obs_mfe is None or sd_shuffled == 0:
        zscore = None
    else:
        zscore = float((obs_mfe - mean_shuffled) / sd_shuffled)

    if obs_mfe is None:
        empirical_p_lower = None
    else:
        empirical_p_lower = float((1 + np.sum(arr <= obs_mfe)) / (1 + len(arr)))

    return {
        "mean_shuffled_mfe": mean_shuffled,
        "sd_shuffled_mfe": sd_shuffled,
        "delta_mfe": delta_mfe,
        "zscore_mfe": zscore,
        "empirical_p_lower": empirical_p_lower,
        "n_shuffles_ok": int(len(arr))
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input table path")
    parser.add_argument("--output-summary", required=True, help="Output summary table path")
    parser.add_argument("--output-detail", required=True, help="Output detailed shuffle table path")
    parser.add_argument(
        "--processes",
        type=int,
        required=True,
        help="Number of worker processes"
    )
    parser.add_argument(
        "--n-shuffles",
        type=int,
        default=100,
        help="Number of shuffles per row (default: 100)"
    )
    parser.add_argument(
        "--klet",
        type=int,
        default=2,
        help="Preserve k-let nucleotide frequencies in ushuffle (default: 2)"
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help="Input/output delimiter (default: tab)"
    )
    parser.add_argument(
        "--flipped-arm-mfe",
        action="store_true",
        help=(
            "Compute additional MFEs with one arm reversed at a time: "
            "reverse lseq only, and reverse rseq only. "
            "Results are saved to summary output columns only."
        )
    )
    args = parser.parse_args()

    print(f"Loading: {args.input}")
    df = pd.read_csv(args.input, sep=args.sep)

    required = {"lseq", "rseq"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if "row_id" not in df.columns:
        df["row_id"] = np.arange(len(df))

    if "name" not in df.columns:
        df["name"] = df["row_id"].astype(str)

    if "mfe" not in df.columns or "dot_bracket" not in df.columns:
        print("Observed mfe/dot_bracket columns not found; computing observed values first...")
        pairs = list(zip(df["lseq"], df["rseq"]))
        with ProcessPoolExecutor(max_workers=args.processes) as ex:
            observed_results = list(ex.map(compute_mfe, pairs, chunksize=200))
        df["mfe"] = [r[0] for r in observed_results]
        df["dot_bracket"] = [r[1] for r in observed_results]

    if args.flipped_arm_mfe:
        print("Computing flipped-arm MFEs (reverse lseq only; reverse rseq only)...")

        left_reversed_pairs = [
            (reverse_sequence(lseq), normalize_sequence(rseq))
            for lseq, rseq in zip(df["lseq"], df["rseq"])
        ]
        right_reversed_pairs = [
            (normalize_sequence(lseq), reverse_sequence(rseq))
            for lseq, rseq in zip(df["lseq"], df["rseq"])
        ]

        with ProcessPoolExecutor(max_workers=args.processes) as ex:
            left_reversed_results = list(ex.map(compute_mfe, left_reversed_pairs, chunksize=200))
            right_reversed_results = list(ex.map(compute_mfe, right_reversed_pairs, chunksize=200))

        df["mfe_lseq_reversed"] = [r[0] for r in left_reversed_results]
        df["flipped_lseq_dot_bracket"] = [r[1] for r in left_reversed_results]
        df["flipped_lseq_pair"] = [
            None if (l is None or r is None) else f"{l}&{r}"
            for l, r in left_reversed_pairs
        ]

        df["mfe_rseq_reversed"] = [r[0] for r in right_reversed_results]
        df["flipped_rseq_dot_bracket"] = [r[1] for r in right_reversed_results]
        df["flipped_rseq_pair"] = [
            None if (l is None or r is None) else f"{l}&{r}"
            for l, r in right_reversed_pairs
        ]

    print(f"Computing shuffled MFEs for {len(df):,} rows with {args.n_shuffles} shuffles per row...")

    summary_rows = []

    with open(args.output_detail, "w", newline="") as detail_handle:
        writer = csv.writer(detail_handle, delimiter=args.sep)
        writer.writerow([
            "row_id",
            "name",
            "shuffle_idx",
            "lseq_shuffle",
            "rseq_shuffle",
            "mfe_shuffle",
            "dot_bracket_shuffle"
        ])

        with ProcessPoolExecutor(max_workers=args.processes) as ex:
            for idx, row in df.iterrows():
                row_id = row["row_id"]
                name = row["name"]
                lseq = row["lseq"]
                rseq = row["rseq"]
                obs_mfe = row["mfe"]
                obs_db = row["dot_bracket"]

                if (idx + 1) % 100 == 0:
                    print(f"  {idx + 1}/{len(df)} rows completed")

                left_shuffles = shuffle_sequence(
                    lseq,
                    number=args.n_shuffles,
                    klet=args.klet,
                    seed=42 + int(row_id)
                )
                right_shuffles = shuffle_sequence(
                    rseq,
                    number=args.n_shuffles,
                    klet=args.klet,
                    seed=420000 + int(row_id)
                )

                n_pairs = min(len(left_shuffles), len(right_shuffles))
                if n_pairs == 0:
                    stats = summarize_null(obs_mfe, [])
                    summary_row = {
                        "row_id": row_id,
                        "name": name,
                        "mfe": obs_mfe,
                        "dot_bracket": obs_db,
                        **stats
                    }
                    if args.flipped_arm_mfe:
                        summary_row["mfe_lseq_reversed"] = row.get("mfe_lseq_reversed")
                        summary_row["flipped_lseq_dot_bracket"] = row.get("flipped_lseq_dot_bracket")
                        summary_row["flipped_lseq_pair"] = row.get("flipped_lseq_pair")
                        summary_row["mfe_rseq_reversed"] = row.get("mfe_rseq_reversed")
                        summary_row["flipped_rseq_dot_bracket"] = row.get("flipped_rseq_dot_bracket")
                        summary_row["flipped_rseq_pair"] = row.get("flipped_rseq_pair")
                    summary_rows.append(summary_row)
                    continue

                pairs = list(zip(left_shuffles[:n_pairs], right_shuffles[:n_pairs]))
                shuffle_results = list(ex.map(compute_mfe, pairs, chunksize=50))

                null_mfes = []
                for shuffle_idx, ((lshuf, rshuf), (mfe_shuffle, db_shuffle)) in enumerate(zip(pairs, shuffle_results)):
                    writer.writerow([
                        row_id,
                        name,
                        shuffle_idx,
                        lshuf,
                        rshuf,
                        mfe_shuffle,
                        db_shuffle
                    ])
                    null_mfes.append(mfe_shuffle)

                stats = summarize_null(obs_mfe, null_mfes)
                summary_row = {
                    "row_id": row_id,
                    "name": name,
                    "mfe": obs_mfe,
                    "dot_bracket": obs_db,
                    **stats
                }
                if args.flipped_arm_mfe:
                    summary_row["mfe_lseq_reversed"] = row.get("mfe_lseq_reversed")
                    summary_row["flipped_lseq_dot_bracket"] = row.get("flipped_lseq_dot_bracket")
                    summary_row["flipped_lseq_pair"] = row.get("flipped_lseq_pair")
                    summary_row["mfe_rseq_reversed"] = row.get("mfe_rseq_reversed")
                    summary_row["flipped_rseq_dot_bracket"] = row.get("flipped_rseq_dot_bracket")
                    summary_row["flipped_rseq_pair"] = row.get("flipped_rseq_pair")
                summary_rows.append(summary_row)

    summary_df = pd.DataFrame(summary_rows)

    df = df.merge(
        summary_df,
        on=["row_id", "name", "mfe", "dot_bracket"],
        how="left"
    )

    print(f"Writing summary: {args.output_summary}")
    df.to_csv(args.output_summary, sep=args.sep, index=False)

    n_ok = df["n_shuffles_ok"].fillna(0).astype(int)
    print(f"Done. Summary written for {len(df):,} rows.")
    print(f"Rows with >=1 successful shuffle MFE: {(n_ok > 0).sum():,}/{len(df):,}")


if __name__ == "__main__":
    main()