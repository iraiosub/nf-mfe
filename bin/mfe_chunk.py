#!/usr/bin/env python3

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import RNA


def compute_mfe(pair):
    """
    Compute RNA-RNA duplex MFE for one sequence pair.
    Expects a tuple: (lseq, rseq)
    """
    lseq, rseq = pair

    if pd.isna(lseq) or pd.isna(rseq):
        return None

    lseq = str(lseq).strip().upper().replace("T", "U")
    rseq = str(rseq).strip().upper().replace("T", "U")

    if not lseq or not rseq:
        return None

    try:
        duplex = RNA.duplexfold(lseq, rseq)
        return duplex.energy
    except Exception:
        return None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input table path")
    parser.add_argument("--output", required=True, help="Output table path")
    parser.add_argument(
        "--processes",
        type=int,
        required=True,
        help="Number of worker processes"
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help="Input/output delimiter (default: tab)"
    )
    args = parser.parse_args()

    print(f"Loading: {args.input}")
    df = pd.read_csv(args.input, sep=args.sep)

    required = {"lseq", "rseq"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    pairs = list(zip(df["lseq"], df["rseq"]))

    print(f"Computing MFEs for {len(pairs):,} rows using {args.processes} processes...")
    with ProcessPoolExecutor(max_workers=args.processes) as ex:
        mfes = list(ex.map(compute_mfe, pairs, chunksize=200))

    df["mfe"] = mfes

    print(f"Writing: {args.output}")
    df.to_csv(args.output, sep=args.sep, index=False)

    n_ok = df["mfe"].notna().sum()
    print(f"Done. Computed MFEs for {n_ok:,}/{len(df):,} rows.")


if __name__ == "__main__":
    main()