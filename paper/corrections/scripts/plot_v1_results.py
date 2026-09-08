"""Redraw the paper's accuracy figures with the v1.0.0 estimates.

The style, the tests and the layout follow `paper/notebooks/results.ipynb` so the new figures can be
read against the published ones. What differs is the LRGE column: it is v1.0.0 rather than the
version the paper benchmarked, and the all-vs-all strategy is absent because it was not re-run.

The other three methods are the paper's own numbers. They were measured on the paper's read sets
rather than the sets rebuilt for these runs, which is sound because the two agree: 0.3.0 on the
rebuilt reads lands within 0.1 percentage points of the paper's own two-set figures on both
platforms, and that row is drawn in the version figure so the agreement is visible.
"""

from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd

DPI = 350
FIGSIZE = (6.4, 4.8)
ROOT = Path(__file__).resolve().parents[2]
FIGURES = ROOT / "corrections/figures"
CORRECTIONS = ROOT / "corrections"
PAPER_ESTIMATES = ROOT / "results/estimates/estimates.tsv"
V1_ESTIMATES = CORRECTIONS / "v1_estimates.tsv"

CUD_PALETTE = [
    "#000000", "#e69f00", "#56b4e9", "#009e73",
    "#f0e442", "#0072b2", "#d55e00", "#cc79a7",
]
LABELS = {"OXFORD_NANOPORE": "ONT", "PACBIO_SMRT": "PacBio"}
METHOD_ORDER = ["lrge 1.0.0", "genomescope", "mash", "raven"]
VERSION_ORDER = ["lrge 0.3.0", "lrge 1.0.0"]


def load():
    v1 = pd.read_csv(V1_ESTIMATES, sep="\t")
    versions = pd.concat(
        [
            v1[["run", "platform", column]]
            .rename(columns={column: "relative_error"})
            .assign(method=name)
            for name, column in (
                ("lrge 0.3.0", "v030_relative_error"),
                ("lrge 1.0.0", "v1_relative_error"),
            )
        ]
    )
    paper = pd.read_csv(PAPER_ESTIMATES, sep="\t")[
        ["run", "platform", "method", "relative_error"]
    ]
    methods = pd.concat(
        [
            v1[["run", "platform", "v1_relative_error"]]
            .rename(columns={"v1_relative_error": "relative_error"})
            .assign(method="lrge 1.0.0"),
            paper[paper["method"].isin(METHOD_ORDER)],
        ]
    )
    for frame in (versions, methods):
        frame["abs_relative_error"] = frame["relative_error"].abs()
    return versions, methods


def legend(ax, location):
    handles, labels = ax.get_legend_handles_labels()
    for handle in handles:
        handle.set_markersize(5)
        handle.set_alpha(1)
    ax.legend(
        handles,
        [LABELS.get(label, label) for label in labels],
        loc=location,
        ncol=2,
        title="",
        alignment="center",
    )


def significant_pairs(data, order):
    """The paper's test: one-way ANOVA per platform, then Tukey HSD on the pairs."""
    found = []
    for platform in data["platform"].unique():
        groups = [
            data[(data["platform"] == platform) & (data["method"] == method)][
                "abs_relative_error"
            ]
            for method in order
        ]
        if f_oneway(*groups).pvalue >= 0.05:
            continue
        subset = data[data["platform"] == platform]
        tukey = pairwise_tukeyhsd(subset["abs_relative_error"], subset["method"])
        for index, pair in enumerate(combinations(tukey.groupsunique, 2)):
            if tukey.reject[index]:
                found.append((pair[0], pair[1], platform, float(tukey.pvalues[index])))
    return found


def asterisks(pvalue):
    for threshold, mark in ((1e-4, "****"), (1e-3, "***"), (1e-2, "**"), (5e-2, "*")):
        if pvalue <= threshold:
            return mark
    return "ns"


def annotate(ax, pairs, order, hue_order):
    if not pairs:
        return
    brackets = []
    for first, second, platform, pvalue in pairs:
        offset = -0.2 if hue_order.index(platform) == 0 else 0.2
        low, high = sorted([order.index(first) + offset, order.index(second) + offset])
        brackets.append((low, high, CUD_PALETTE[hue_order.index(platform)], asterisks(pvalue)))
    heights = np.logspace(4, 6, num=len(brackets), base=10)
    for height, (low, high, colour, mark) in zip(heights, sorted(brackets)):
        cap = 0.1 * height
        ax.plot(
            [low, low, high, high],
            [height, height + cap, height + cap, height],
            lw=1,
            c=colour,
        )
        ax.text(
            (low + high) * 0.5, height * 0.85, mark,
            ha="center", va="bottom", color=colour, fontsize=8,
        )


def absolute_error_figure(data, order, name, annotated):
    fig, ax = plt.subplots(dpi=DPI, figsize=FIGSIZE)
    hue_order = sorted(data["platform"].unique())
    ax.set_yscale("symlog", linthresh=1)
    shared = dict(
        x="method", y="abs_relative_error", hue="platform", data=data, ax=ax,
        order=order, hue_order=hue_order,
    )
    sns.stripplot(**shared, dodge=True, alpha=0.1, jitter=0.2, size=2)
    sns.violinplot(
        **shared, fill=False, linewidth=1.5, density_norm="width", legend=False,
        inner="quart", cut=0, split=True,
    )
    ticks = [0, 1, 10, 100, 1000, 10000]
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticks)
    if annotated:
        annotate(ax, significant_pairs(data, order), order, hue_order)
    legend(ax, "upper left")
    ax.set_ylabel("Absolute relative error (%)")
    ax.set_xlabel("")
    ax.yaxis.grid(True, alpha=0.5)
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"{name}.{suffix}", bbox_inches="tight")
    plt.close(fig)


def signed_error_figure(data, order, name):
    fig, ax = plt.subplots(dpi=DPI, figsize=FIGSIZE)
    hue_order = sorted(data["platform"].unique())
    ax.set_yscale("symlog", linthresh=1)
    shared = dict(
        x="method", y="relative_error", hue="platform", data=data, ax=ax,
        order=order, hue_order=hue_order,
    )
    sns.violinplot(
        **shared, fill=False, linewidth=1.5, density_norm="width", legend=False,
        split=True, cut=0, inner="quart",
    )
    sns.stripplot(**shared, dodge=True, alpha=0.1, jitter=0.15, size=2)
    ticks = [-100, -10, -1, 0, 1, 10, 100, 1000, 10000, 100000]
    ax.set_yticks(ticks)
    ax.set_yticklabels(ticks)
    ax.set_xlabel("")
    ax.set_ylabel("Relative error (%)")
    legend(ax, "best")
    for suffix in ("png", "pdf"):
        fig.savefig(FIGURES / f"{name}.{suffix}", bbox_inches="tight")
    plt.close(fig)


def summarise(data, path):
    table = data.groupby(["method", "platform"])["abs_relative_error"].agg(
        n="size",
        median="median",
        mean="mean",
        within_5_pct=lambda column: (column <= 5).mean() * 100,
        within_10_pct=lambda column: (column <= 10).mean() * 100,
    )
    table.round(4).to_csv(path, sep="\t")
    return table


def main():
    sns.set_style("whitegrid")
    sns.set_palette(CUD_PALETTE)
    FIGURES.mkdir(parents=True, exist_ok=True)

    versions, methods = load()

    absolute_error_figure(methods, METHOD_ORDER, "method_absolute_relative_error", True)
    signed_error_figure(methods, METHOD_ORDER, "platform_relative_error")
    absolute_error_figure(versions, VERSION_ORDER, "version_absolute_relative_error", False)

    print(summarise(methods, CORRECTIONS / "v1_method_accuracy.tsv").round(2).to_string())
    print()
    print(summarise(versions, CORRECTIONS / "v1_version_accuracy.tsv").round(2).to_string())


if __name__ == "__main__":
    main()
