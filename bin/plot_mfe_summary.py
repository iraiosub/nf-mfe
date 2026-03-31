#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    from scipy.stats import gaussian_kde, wilcoxon
    SCIPY_AVAILABLE = True
except Exception:
    SCIPY_AVAILABLE = False


def ecdf(x):
    x = np.asarray(x)
    x = x[np.isfinite(x)]
    x = np.sort(x)
    y = np.arange(1, len(x) + 1) / len(x)
    return x, y


def bh_fdr(pvals):
    pvals = np.asarray(pvals, dtype=float)
    out = np.full_like(pvals, np.nan, dtype=float)

    valid = np.isfinite(pvals)
    if valid.sum() == 0:
        return out

    pv = pvals[valid]
    n = len(pv)
    order = np.argsort(pv)
    ranked = pv[order]

    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)

    tmp = np.empty_like(q)
    tmp[order] = q
    out[valid] = tmp
    return out


def density_line(ax, values, label, linewidth=2):
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return

    if SCIPY_AVAILABLE and len(np.unique(values)) > 1:
        kde = gaussian_kde(values)
        x = np.linspace(values.min(), values.max(), 500)
        y = kde(x)
        ax.plot(x, y, label=label, linewidth=linewidth)
    else:
        ax.hist(values, bins=50, density=True, histtype="step", label=label, linewidth=linewidth)


def safe_wilcoxon(x, y):
    if not SCIPY_AVAILABLE:
        return np.nan, np.nan
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) == 0:
        return np.nan, np.nan
    if np.allclose(x, y):
        return np.nan, np.nan
    try:
        stat, p = wilcoxon(x, y, zero_method="wilcox", alternative="two-sided")
        return stat, p
    except Exception:
        return np.nan, np.nan


def main():
    parser = argparse.ArgumentParser(description="Plot observed vs shuffled MFE summary distributions.")
    parser.add_argument("--input", required=True, help="Input summary TSV from shuffled MFE script")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument("--sample-name", default=None, help="Optional sample name for plot titles")
    parser.add_argument("--sep", default="\t", help="Delimiter for input/output tables (default: tab)")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.input, sep=args.sep)

    required = {
        "mfe",
        "mean_shuffled_mfe",
        "delta_mfe",
        "zscore_mfe",
        "empirical_p_lower"
    }
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    sample_name = args.sample_name if args.sample_name else Path(args.input).stem

    n_before = len(df)
    df = df[df["mfe"].isna() | (pd.to_numeric(df["mfe"], errors="coerce") <= 0)].copy()
    n_after = len(df)
    n_removed = n_before - n_after
    print(f"Filtered out {n_removed:,} rows with mfe > 0; retained {n_after:,}/{n_before:,} rows.")

    plot_df = df.copy()
    plot_df["fdr_bh"] = bh_fdr(plot_df["empirical_p_lower"].values)

    obs = plot_df["mfe"].astype(float).values
    null_mean = plot_df["mean_shuffled_mfe"].astype(float).values
    delta = plot_df["delta_mfe"].astype(float).values
    z = plot_df["zscore_mfe"].astype(float).values
    p = plot_df["empirical_p_lower"].astype(float).values
    q = plot_df["fdr_bh"].astype(float).values

    paired_mask = np.isfinite(obs) & np.isfinite(null_mean)
    obs_paired = obs[paired_mask]
    null_paired = null_mean[paired_mask]

    n_total = len(plot_df)
    n_paired = paired_mask.sum()
    frac_stronger = np.mean(delta[np.isfinite(delta)] < 0) if np.isfinite(delta).sum() else np.nan

    wilcox_stat, wilcox_p = safe_wilcoxon(obs, null_mean)

    stats = {
        "sample": sample_name,
        "n_rows_total": n_total,
        "n_rows_paired_valid": int(n_paired),
        "mean_observed_mfe": np.nanmean(obs),
        "median_observed_mfe": np.nanmedian(obs),
        "sd_observed_mfe": np.nanstd(obs, ddof=1),
        "mean_shuffled_mean_mfe": np.nanmean(null_mean),
        "median_shuffled_mean_mfe": np.nanmedian(null_mean),
        "sd_shuffled_mean_mfe": np.nanstd(null_mean, ddof=1),
        "mean_delta_mfe": np.nanmean(delta),
        "median_delta_mfe": np.nanmedian(delta),
        "sd_delta_mfe": np.nanstd(delta, ddof=1),
        "fraction_delta_mfe_lt_0": frac_stronger,
        "mean_zscore_mfe": np.nanmean(z),
        "median_zscore_mfe": np.nanmedian(z),
        "fraction_empirical_p_lt_0_05": np.mean(p[np.isfinite(p)] < 0.05) if np.isfinite(p).sum() else np.nan,
        "fraction_empirical_p_lt_0_01": np.mean(p[np.isfinite(p)] < 0.01) if np.isfinite(p).sum() else np.nan,
        "fraction_fdr_lt_0_05": np.mean(q[np.isfinite(q)] < 0.05) if np.isfinite(q).sum() else np.nan,
        "wilcoxon_statistic_observed_vs_shuffled_mean": wilcox_stat,
        "wilcoxon_pvalue_observed_vs_shuffled_mean": wilcox_p,
    }

    stats_df = pd.DataFrame([stats])
    stats_df.to_csv(outdir / f"{sample_name}.mfe_summary_stats.tsv", sep=args.sep, index=False)

    plt.rcParams.update({
        "figure.figsize": (16, 12),
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.2,
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "legend.frameon": False,
    })

    fig, axes = plt.subplots(2, 3, figsize=(16, 12))
    fig.suptitle(f"{sample_name}: observed vs shuffled MFE summary", fontsize=18, y=0.98)

    ax = axes[0, 0]
    density_line(ax, obs, "Observed MFE")
    density_line(ax, null_mean, "Mean shuffled MFE")
    ax.axvline(np.nanmedian(obs), linestyle="--", linewidth=1.5)
    ax.axvline(np.nanmedian(null_mean), linestyle="--", linewidth=1.5)
    ax.set_title("Density")
    ax.set_xlabel("MFE (kcal/mol)")
    ax.set_ylabel("Density")
    ax.legend()

    ax = axes[0, 1]
    x1, y1 = ecdf(obs)
    x2, y2 = ecdf(null_mean)
    if len(x1):
        ax.plot(x1, y1, linewidth=2, label="Observed MFE")
    if len(x2):
        ax.plot(x2, y2, linewidth=2, label="Mean shuffled MFE")
    ax.set_title("ECDF")
    ax.set_xlabel("MFE (kcal/mol)")
    ax.set_ylabel("Cumulative fraction")
    ax.legend()

    ax = axes[0, 2]
    delta_valid = delta[np.isfinite(delta)]
    if len(delta_valid):
        density_line(ax, delta_valid, "Delta MFE")
        ax.axvline(0, linestyle="--", linewidth=1.5)
        ax.axvline(np.nanmedian(delta_valid), linestyle=":", linewidth=1.5)
    ax.set_title("Delta MFE distribution")
    ax.set_xlabel("delta_mfe = observed - mean_shuffled")
    ax.set_ylabel("Density")
    ax.legend()

    ax = axes[1, 0]
    p_valid = p[np.isfinite(p)]
    if len(p_valid):
        bins = np.linspace(0, 1, 41)
        ax.hist(p_valid, bins=bins, edgecolor="black", linewidth=0.6)
        ax.axvline(0.05, linestyle="--", linewidth=1.5)
        ax.axvline(0.01, linestyle=":", linewidth=1.5)
    ax.set_title("Empirical p-value distribution")
    ax.set_xlabel("empirical_p_lower")
    ax.set_ylabel("Count")

    ax = axes[1, 1]
    z_valid = z[np.isfinite(z)]
    if len(z_valid):
        density_line(ax, z_valid, "Z-score")
        ax.axvline(0, linestyle="--", linewidth=1.5)
        ax.axvline(np.nanmedian(z_valid), linestyle=":", linewidth=1.5)
    ax.set_title("Z-score distribution")
    ax.set_xlabel("zscore_mfe")
    ax.set_ylabel("Density")
    ax.legend()

    ax = axes[1, 2]
    if len(obs_paired):
        ax.scatter(null_paired, obs_paired, s=12, alpha=0.35)
        lo = np.nanmin([obs_paired.min(), null_paired.min()])
        hi = np.nanmax([obs_paired.max(), null_paired.max()])
        ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.5)
    ax.set_title("Observed vs shuffled mean")
    ax.set_xlabel("Mean shuffled MFE")
    ax.set_ylabel("Observed MFE")

    text = (
        f"n = {n_paired:,}\n"
        f"mean obs = {np.nanmean(obs):.2f}\n"
        f"mean shuffled = {np.nanmean(null_mean):.2f}\n"
        f"mean delta = {np.nanmean(delta):.2f}\n"
        f"frac(delta < 0) = {frac_stronger:.3f}\n"
        f"frac(p < 0.05) = {stats['fraction_empirical_p_lt_0_05']:.3f}\n"
        f"frac(FDR < 0.05) = {stats['fraction_fdr_lt_0_05']:.3f}\n"
        f"Wilcoxon p = {wilcox_p:.3e}" if np.isfinite(wilcox_p) else
        f"n = {n_paired:,}\n"
        f"mean obs = {np.nanmean(obs):.2f}\n"
        f"mean shuffled = {np.nanmean(null_mean):.2f}\n"
        f"mean delta = {np.nanmean(delta):.2f}\n"
        f"frac(delta < 0) = {frac_stronger:.3f}\n"
        f"frac(p < 0.05) = {stats['fraction_empirical_p_lt_0_05']:.3f}\n"
        f"frac(FDR < 0.05) = {stats['fraction_fdr_lt_0_05']:.3f}"
    )
    fig.text(0.68, 0.08, text, fontsize=11, va="bottom")

    fig.tight_layout(rect=(0, 0.04, 1, 0.96))

    png_path = outdir / f"{sample_name}.mfe_summary_plots.png"
    pdf_path = outdir / f"{sample_name}.mfe_summary_plots.pdf"

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote plots: {png_path}")
    print(f"Wrote plots: {pdf_path}")
    print(f"Wrote stats: {outdir / f'{sample_name}.mfe_summary_stats.tsv'}")


if __name__ == "__main__":
    main()