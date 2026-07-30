#!/usr/bin/env python3
"""Export local-only SpaGE versus Tangram comparison tables and figures."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


COLORS = {"spage": "#3B6FB6", "tangram": "#E5863B"}


def parse_args() -> argparse.Namespace:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--metrics-dir",
        type=Path,
        default=root / "git_eval/metrics",
        help="Directory containing week1_meld_fold*_BACKEND.tsv files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=root / "import_result/week1_model_comparison",
        help="Local-only output directory (ignored by git).",
    )
    parser.add_argument("--dpi", type=int, default=300)
    return parser.parse_args()


def _read_backend(metrics_dir: Path, backend: str) -> pd.DataFrame:
    paths = sorted(metrics_dir.glob(f"week1_meld_fold*_{backend}.tsv"))
    if not paths:
        raise FileNotFoundError(
            f"no per-gene {backend} metrics under {metrics_dir}"
        )
    tables = []
    required = {
        "gene_id",
        "backend",
        "fold",
        "spearman",
        "z_nrmse",
        "truth_detected_spots",
        "truth_total",
    }
    for path in paths:
        table = pd.read_csv(path, sep="\t")
        missing = required - set(table.columns)
        if missing:
            raise ValueError(f"{path} is missing columns: {sorted(missing)}")
        if not table["backend"].eq(backend).all():
            raise ValueError(f"{path} contains an unexpected backend label")
        tables.append(table)
    result = pd.concat(tables, ignore_index=True)
    if result.duplicated(["fold", "gene_id"]).any():
        raise ValueError(f"duplicate {backend} fold/gene rows")
    return result


def build_comparison(metrics_dir: Path) -> pd.DataFrame:
    spage = _read_backend(metrics_dir, "spage")
    tangram = _read_backend(metrics_dir, "tangram")
    keys = ["fold", "gene_id"]
    truth = ["truth_detected_spots", "truth_total"]
    comparison = spage.merge(
        tangram,
        on=keys,
        how="inner",
        validate="one_to_one",
        suffixes=("_spage", "_tangram"),
    )
    if len(comparison) != len(spage) or len(comparison) != len(tangram):
        raise ValueError("SpaGE and Tangram do not contain identical fold/gene rows")
    for column in truth:
        left = comparison[f"{column}_spage"].to_numpy(dtype=float)
        right = comparison[f"{column}_tangram"].to_numpy(dtype=float)
        if not np.allclose(left, right, equal_nan=True):
            raise ValueError(f"truth mismatch between backends: {column}")
        comparison[column] = left

    comparison = comparison.rename(
        columns={
            "spearman_spage": "spage_spearman",
            "spearman_tangram": "tangram_spearman",
            "z_nrmse_spage": "spage_z_nrmse",
            "z_nrmse_tangram": "tangram_z_nrmse",
        }
    )
    comparison["delta_spearman_tangram_minus_spage"] = (
        comparison["tangram_spearman"] - comparison["spage_spearman"]
    )
    comparison["delta_z_nrmse_spage_minus_tangram"] = (
        comparison["spage_z_nrmse"] - comparison["tangram_z_nrmse"]
    )
    comparison["spearman_winner"] = np.where(
        comparison["tangram_spearman"] > comparison["spage_spearman"],
        "Tangram",
        np.where(
            comparison["tangram_spearman"] < comparison["spage_spearman"],
            "SpaGE",
            "Tie",
        ),
    )
    comparison["z_nrmse_winner"] = np.where(
        comparison["tangram_z_nrmse"] < comparison["spage_z_nrmse"],
        "Tangram",
        np.where(
            comparison["tangram_z_nrmse"] > comparison["spage_z_nrmse"],
            "SpaGE",
            "Tie",
        ),
    )
    comparison["detection_bin"] = pd.cut(
        comparison["truth_detected_spots"],
        bins=[0, 5, 10, np.inf],
        labels=["3-5 spots", "6-10 spots", ">10 spots"],
        include_lowest=True,
    )
    columns = [
        "fold",
        "gene_id",
        "truth_detected_spots",
        "truth_total",
        "detection_bin",
        "spage_spearman",
        "tangram_spearman",
        "delta_spearman_tangram_minus_spage",
        "spearman_winner",
        "spage_z_nrmse",
        "tangram_z_nrmse",
        "delta_z_nrmse_spage_minus_tangram",
        "z_nrmse_winner",
    ]
    return comparison[columns].sort_values(["fold", "gene_id"])


def summarize(comparison: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for backend in ("spage", "tangram"):
        rho = comparison[f"{backend}_spearman"]
        error = comparison[f"{backend}_z_nrmse"]
        rows.append(
            {
                "backend": backend,
                "targets": len(comparison),
                "median_spearman": rho.median(),
                "mean_spearman": rho.mean(),
                "positive_spearman": int(rho.gt(0).sum()),
                "spearman_ge_0.1": int(rho.ge(0.1).sum()),
                "median_z_nrmse": error.median(),
                "mean_z_nrmse": error.mean(),
            }
        )
    summary = pd.DataFrame(rows)

    strata = []
    for detection_bin, group in comparison.groupby(
        "detection_bin", observed=True
    ):
        for backend in ("spage", "tangram"):
            strata.append(
                {
                    "detection_bin": str(detection_bin),
                    "backend": backend,
                    "targets": len(group),
                    "median_detected_spots": group[
                        "truth_detected_spots"
                    ].median(),
                    "median_spearman": group[f"{backend}_spearman"].median(),
                    "median_z_nrmse": group[f"{backend}_z_nrmse"].median(),
                }
            )
    return summary, pd.DataFrame(strata)


def write_decision(
    path: Path,
    comparison: pd.DataFrame,
    summary: pd.DataFrame,
) -> None:
    s = summary.set_index("backend")
    tangram_rho_wins = int(comparison["spearman_winner"].eq("Tangram").sum())
    tangram_error_wins = int(comparison["z_nrmse_winner"].eq("Tangram").sum())
    provisional = (
        s.loc["tangram", "median_spearman"]
        > s.loc["spage", "median_spearman"]
        and s.loc["tangram", "median_z_nrmse"]
        < s.loc["spage", "median_z_nrmse"]
    )
    recommendation = (
        "Tangram is the provisional backbone; keep SpaGE as the comparator."
        if provisional
        else "No provisional winner; retain both backends for additional folds."
    )
    text = f"""# Week 1 model comparison

## Scope

- Targets compared: {len(comparison)}
- Folds present: {", ".join(map(str, sorted(comparison["fold"].unique())))}
- Median detected spots: {comparison["truth_detected_spots"].median():.1f}

## Results

| Model | Median Spearman | Median z-NRMSE | Positive Spearman |
|---|---:|---:|---:|
| SpaGE | {s.loc["spage", "median_spearman"]:.4f} | {s.loc["spage", "median_z_nrmse"]:.4f} | {int(s.loc["spage", "positive_spearman"])}/{len(comparison)} |
| Tangram | {s.loc["tangram", "median_spearman"]:.4f} | {s.loc["tangram", "median_z_nrmse"]:.4f} | {int(s.loc["tangram", "positive_spearman"])}/{len(comparison)} |

Tangram has higher per-gene Spearman for {tangram_rho_wins}/{len(comparison)}
targets and lower z-NRMSE for {tangram_error_wins}/{len(comparison)} targets.

## Decision

{recommendation}

This is an engineering smoke benchmark, not evidence that production lncRNA
mapping accuracy has passed. Complete all five folds and null-calibrated
evaluation before freezing the model choice.
"""
    path.write_text(text)


def plot_comparison(
    comparison: pd.DataFrame,
    output_png: Path,
    output_pdf: Path,
    dpi: int,
) -> None:
    import os

    matplotlib_config = output_png.parent / ".matplotlib"
    matplotlib_config.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(matplotlib_config))
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 8.2), constrained_layout=True)

    ax = axes[0, 0]
    x = comparison["spage_spearman"]
    y = comparison["tangram_spearman"]
    limits = [min(x.min(), y.min()) - 0.02, max(x.max(), y.max()) + 0.02]
    ax.scatter(x, y, c="#555555", s=32, alpha=0.8, edgecolors="white", linewidth=0.4)
    ax.plot(limits, limits, "--", color="#999999", linewidth=1)
    ax.axhline(0, color="#DDDDDD", linewidth=0.8)
    ax.axvline(0, color="#DDDDDD", linewidth=0.8)
    ax.set(xlim=limits, ylim=limits, xlabel="SpaGE Spearman", ylabel="Tangram Spearman")
    ax.set_title("A  Per-target rank correlation")

    ax = axes[0, 1]
    x = comparison["spage_z_nrmse"]
    y = comparison["tangram_z_nrmse"]
    limits = [min(x.min(), y.min()) - 0.01, max(x.max(), y.max()) + 0.01]
    ax.scatter(x, y, c="#555555", s=32, alpha=0.8, edgecolors="white", linewidth=0.4)
    ax.plot(limits, limits, "--", color="#999999", linewidth=1)
    ax.set(xlim=limits, ylim=limits, xlabel="SpaGE z-NRMSE", ylabel="Tangram z-NRMSE")
    ax.set_title("B  Per-target normalized error")

    ax = axes[1, 0]
    values = [
        comparison["spage_spearman"].dropna().to_numpy(),
        comparison["tangram_spearman"].dropna().to_numpy(),
    ]
    boxes = ax.boxplot(values, tick_labels=["SpaGE", "Tangram"], patch_artist=True)
    for patch, backend in zip(boxes["boxes"], ("spage", "tangram")):
        patch.set_facecolor(COLORS[backend])
        patch.set_alpha(0.55)
    rng = np.random.default_rng(0)
    for position, (backend, data) in enumerate(
        zip(("spage", "tangram"), values), start=1
    ):
        jitter = rng.normal(0, 0.035, len(data))
        ax.scatter(
            np.full(len(data), position) + jitter,
            data,
            s=16,
            color=COLORS[backend],
            alpha=0.65,
            edgecolors="none",
        )
    ax.axhline(0, color="#BBBBBB", linewidth=0.8)
    ax.set_ylabel("Spearman correlation")
    ax.set_title("C  Distribution across targets")

    ax = axes[1, 1]
    detected = comparison["truth_detected_spots"]
    for backend in ("spage", "tangram"):
        ax.scatter(
            detected,
            comparison[f"{backend}_spearman"],
            s=28,
            alpha=0.7,
            label=backend.capitalize(),
            color=COLORS[backend],
            edgecolors="white",
            linewidth=0.35,
        )
    ax.set_xscale("log")
    ax.axhline(0, color="#BBBBBB", linewidth=0.8)
    ax.set(
        xlabel="Truth-detected spots (log scale)",
        ylabel="Spearman correlation",
    )
    ax.legend(frameon=False)
    ax.set_title("D  Performance versus target support")

    fig.suptitle(
        "LncSpaceMap Week 1: SpaGE versus Tangram",
        fontsize=14,
        fontweight="bold",
    )
    fig.savefig(output_png, dpi=dpi, bbox_inches="tight", facecolor="white")
    fig.savefig(output_pdf, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    comparison = build_comparison(args.metrics_dir)
    summary, strata = summarize(comparison)
    comparison.to_csv(
        args.output_dir / "model_comparison_per_gene.tsv",
        sep="\t",
        index=False,
    )
    summary.to_csv(
        args.output_dir / "model_comparison_summary.tsv",
        sep="\t",
        index=False,
    )
    strata.to_csv(
        args.output_dir / "model_comparison_by_detection.tsv",
        sep="\t",
        index=False,
    )
    write_decision(
        args.output_dir / "model_comparison_decision.md",
        comparison,
        summary,
    )
    plot_comparison(
        comparison,
        args.output_dir / "model_comparison_4panel.png",
        args.output_dir / "model_comparison_4panel.pdf",
        args.dpi,
    )
    print(f"PASS_LOCAL_MODEL_COMPARISON output={args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
