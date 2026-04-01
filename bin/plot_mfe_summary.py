#!/usr/bin/env python3

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D

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


def density_line(ax, values, label, linewidth=2, color=None):
    values = np.asarray(values)
    values = values[np.isfinite(values)]
    if len(values) < 2:
        return

    if SCIPY_AVAILABLE and len(np.unique(values)) > 1:
        kde = gaussian_kde(values)
        x = np.linspace(values.min(), values.max(), 500)
        y = kde(x)
        ax.plot(x, y, label=label, linewidth=linewidth, color=color)
    else:
        ax.hist(values, bins=50, density=True, histtype="step", label=label, linewidth=linewidth, color=color)


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


def add_median_and_landmark_legend(ax, series_handles=None, median_label="Median", landmark_label="Landmark"):
    handles = []
    labels = []

    if series_handles:
        handles.extend(series_handles)
        labels.extend([h.get_label() for h in series_handles])

    handles.append(Line2D([0], [0], color="black", linestyle="--", linewidth=1.5, label=median_label))
    labels.append(median_label)

    handles.append(Line2D([0], [0], color="black", linestyle=":", linewidth=1.5, label=landmark_label))
    labels.append(landmark_label)

    ax.legend(handles=handles, labels=labels)


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
    obs_median = np.nanmedian(obs)
    null_median = np.nanmedian(null_mean)
    ax.axvline(obs_median, linestyle="--", linewidth=1.5, color="C0")
    ax.axvline(null_median, linestyle="--", linewidth=1.5, color="C1")
    ax.set_title(f"Density (n={n_paired:,})")
    ax.set_xlabel("MFE (kcal/mol)")
    ax.set_ylabel("Density")
    ax.legend()

    ax = axes[0, 1]
    x1, y1 = ecdf(obs)
    x2, y2 = ecdf(null_mean)
    if len(x1):
        ax.plot(x1, y1, linewidth=2, label="Observed MFE")
        ax.axvline(np.nanmedian(obs), linestyle="--", linewidth=1.5, color="C0")
    if len(x2):
        ax.plot(x2, y2, linewidth=2, label="Mean shuffled MFE")
        ax.axvline(np.nanmedian(null_mean), linestyle="--", linewidth=1.5, color="C1")
    ax.set_title(f"ECDF (n={n_paired:,})")
    ax.set_xlabel("MFE (kcal/mol)")
    ax.set_ylabel("Cumulative fraction")
    ax.legend()

    ax = axes[0, 2]
    delta_valid = delta[np.isfinite(delta)]
    if len(delta_valid):
        density_line(ax, delta_valid, "Delta MFE", color="gray")
        ax.axvline(np.nanmedian(delta_valid), linestyle="--", linewidth=1.5, color="gray")
        ax.axvline(0, linestyle=":", linewidth=1.5, color="gray")
    ax.set_title(f"Delta MFE distribution (n={len(delta_valid):,})")
    ax.set_xlabel("delta_mfe = observed - mean_shuffled")
    ax.set_ylabel("Density")
    ax.legend()

    ax = axes[1, 0]
    p_valid = p[np.isfinite(p)]
    if len(p_valid):
        bins = np.linspace(0, 1, 41)
        ax.hist(p_valid, bins=bins, edgecolor="black", linewidth=0.6, color="gray")
        ax.axvline(0.05, linestyle=":", linewidth=1.5, color="gray")
        ax.axvline(0.01, linestyle=":", linewidth=1.5, color="black")
    ax.set_title(f"Empirical p-value distribution (n={len(p_valid):,})")
    ax.set_xlabel("empirical_p_lower")
    ax.set_ylabel("Count")

    ax = axes[1, 1]
    z_valid = z[np.isfinite(z)]
    if len(z_valid):
        density_line(ax, z_valid, "Z-score", color="gray")
        ax.axvline(np.nanmedian(z_valid), linestyle="--", linewidth=1.5, color="gray")
        ax.axvline(0, linestyle=":", linewidth=1.5, color="gray")
    ax.set_title(f"Z-score distribution (n={len(z_valid):,})")
    ax.set_xlabel("zscore_mfe")
    ax.set_ylabel("Density")
    ax.legend()

    ax = axes[1, 2]
    if len(obs_paired):
        ax.scatter(null_paired, obs_paired, s=12, alpha=0.2, color="gray")
        lo = np.nanmin([obs_paired.min(), null_paired.min()])
        hi = np.nanmax([obs_paired.max(), null_paired.max()])
        ax.plot([lo, hi], [lo, hi], linestyle=":", linewidth=1.5, color="gray")
    ax.set_title(f"Observed vs shuffled mean (n={n_paired:,})")
    ax.set_xlabel("Mean shuffled MFE")
    ax.set_ylabel("Observed MFE")

    fig.tight_layout(rect=(0, 0.04, 1, 0.96))

    png_path = outdir / f"{sample_name}.shuffled_mfe_summary_plots.png"
    pdf_path = outdir / f"{sample_name}.shuffled_mfe_summary_plots.pdf"

    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    with PdfPages(pdf_path) as pdf:
        pdf.savefig(fig, bbox_inches="tight")

        flipped_cols_present = {"mfe_lseq_reversed", "mfe_rseq_reversed"}.issubset(plot_df.columns)
        if flipped_cols_present:
            flip_l = plot_df["mfe_lseq_reversed"].astype(float).values
            flip_r = plot_df["mfe_rseq_reversed"].astype(float).values

            delta_flip_l = flip_l - null_mean
            delta_flip_r = flip_r - null_mean

            flip_l_valid = flip_l[np.isfinite(flip_l)]
            flip_r_valid = flip_r[np.isfinite(flip_r)]
            delta_flip_l_valid = delta_flip_l[np.isfinite(delta_flip_l)]
            delta_flip_r_valid = delta_flip_r[np.isfinite(delta_flip_r)]

            fig2, axes2 = plt.subplots(1, 3, figsize=(18, 5))
            fig2.suptitle(f"{sample_name}: flipped-arm MFE summary", fontsize=18, y=0.98)

            ax = axes2[0]
            density_line(ax, flip_l_valid, "Reverse lseq MFE", color="C0")
            density_line(ax, flip_r_valid, "Reverse rseq MFE", color="C1")
            if len(flip_l_valid):
                ax.axvline(np.nanmedian(flip_l_valid), linestyle="--", linewidth=1.5, color="C0")
            if len(flip_r_valid):
                ax.axvline(np.nanmedian(flip_r_valid), linestyle="--", linewidth=1.5, color="C1")
            ax.axvline(0, linestyle=":", linewidth=1.5, color="gray")
            ax.set_title(f"Flipped arm density (nL={len(flip_l_valid):,}, nR={len(flip_r_valid):,})")
            ax.set_xlabel("MFE (kcal/mol)")
            ax.set_ylabel("Density")
            add_median_and_landmark_legend(
                ax,
                series_handles=[
                    Line2D([0], [0], color="C0", linewidth=2, label="Reverse lseq MFE"),
                    Line2D([0], [0], color="C1", linewidth=2, label="Reverse rseq MFE"),
                ],
                median_label="Median",
                landmark_label="0 landmark",
            )

            ax = axes2[1]
            if len(flip_l_valid):
                x, y = ecdf(flip_l_valid)
                ax.plot(x, y, linewidth=2, color="C0", label="Reverse lseq MFE")
                ax.axvline(np.nanmedian(flip_l_valid), linestyle="--", linewidth=1.5, color="C0")
            if len(flip_r_valid):
                x, y = ecdf(flip_r_valid)
                ax.plot(x, y, linewidth=2, color="C1", label="Reverse rseq MFE")
                ax.axvline(np.nanmedian(flip_r_valid), linestyle="--", linewidth=1.5, color="C1")
            ax.axvline(0, linestyle=":", linewidth=1.5, color="gray")
            ax.set_title("Flipped ECDF")
            ax.set_xlabel("MFE (kcal/mol)")
            ax.set_ylabel("Cumulative fraction")
            add_median_and_landmark_legend(
                ax,
                series_handles=[
                    Line2D([0], [0], color="C0", linewidth=2, label="Reverse lseq MFE"),
                    Line2D([0], [0], color="C1", linewidth=2, label="Reverse rseq MFE"),
                ],
                median_label="Median",
                landmark_label="0 landmark",
            )

            ax = axes2[2]
            density_line(ax, delta_flip_l_valid, "Delta reverse lseq", color="C0")
            density_line(ax, delta_flip_r_valid, "Delta reverse rseq", color="C1")
            if len(delta_flip_l_valid):
                ax.axvline(np.nanmedian(delta_flip_l_valid), linestyle="--", linewidth=1.5, color="C0")
            if len(delta_flip_r_valid):
                ax.axvline(np.nanmedian(delta_flip_r_valid), linestyle="--", linewidth=1.5, color="C1")
            ax.axvline(0, linestyle=":", linewidth=1.5, color="gray")
            ax.set_title(
                f"Flipped delta MFE (nL={len(delta_flip_l_valid):,}, nR={len(delta_flip_r_valid):,})"
            )
            ax.set_xlabel("delta_mfe = flipped - mean_shuffled")
            ax.set_ylabel("Density")
            add_median_and_landmark_legend(
                ax,
                series_handles=[
                    Line2D([0], [0], color="C0", linewidth=2, label="Delta reverse lseq"),
                    Line2D([0], [0], color="C1", linewidth=2, label="Delta reverse rseq"),
                ],
                median_label="Median",
                landmark_label="0 landmark",
            )

            fig2.tight_layout(rect=(0, 0.02, 1, 0.93))

            flipped_png_path = outdir / f"{sample_name}.flipped_mfe_summary_plots.png"
            flipped_pdf_path = outdir / f"{sample_name}.flipped_mfe_summary_plots.pdf"

            fig2.savefig(flipped_png_path, dpi=300, bbox_inches="tight")
            pdf.savefig(fig2, bbox_inches="tight")
            plt.close(fig2)

            print(f"Wrote plots: {flipped_png_path}")
            print(f"Wrote plots: {flipped_pdf_path}")

    plt.close(fig)

    print(f"Wrote plots: {png_path}")
    print(f"Wrote plots: {pdf_path}")
    print(f"Wrote stats: {outdir / f'{sample_name}.mfe_summary_stats.tsv'}")


if __name__ == "__main__":
    main()