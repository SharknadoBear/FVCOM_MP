"""Analyze the Memo 11 FVCOM-MP case-B speed-test results.

The observed timings are stored here instead of in the memo body so speedup,
efficiency, and plots are generated reproducibly.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


@dataclass(frozen=True)
class TimingRun:
    run_id: str
    role: str
    nodes: int
    tasks_per_node: int
    ranks: int
    sec_per_it: float
    note: str


@dataclass(frozen=True)
class TimingMetric:
    run: TimingRun
    speedup: float
    efficiency: float
    iterations_per_second: float


ROOT_DIR = Path(__file__).resolve().parents[1]
MEMORY_DIR = ROOT_DIR / "Memory"
FIGURE_DIR = MEMORY_DIR / "memo_11_figures"
SUMMARY_TEX = MEMORY_DIR / "memo_11_speed_analysis.tex"


RUNS = [
    TimingRun("B01", "single-core baseline", 1, 1, 1, 11.61, "baseline timing"),
    TimingRun("B02", "small within-node", 1, 2, 2, 5.72, "corrected timing from user"),
    TimingRun("B03", "small within-node", 1, 4, 4, 2.99, "observed timing"),
    TimingRun("B04", "small within-node", 1, 8, 8, 1.45, "observed timing"),
    TimingRun("B05", "medium within-node", 1, 16, 16, 0.77, "observed timing"),
    TimingRun("B06", "medium within-node", 1, 32, 32, 0.41, "observed timing"),
    TimingRun("B07", "half node", 1, 52, 52, 0.25, "observed timing"),
    TimingRun("B08", "dense within-node", 1, 64, 64, 0.22, "observed timing"),
    TimingRun("B09", "dense within-node", 1, 78, 78, 0.18, "observed timing"),
    TimingRun("B10", "full node", 1, 104, 104, 0.15, "observed timing"),
    TimingRun("B11", "spread-vs-pack check", 2, 52, 104, 0.12, "observed timing"),
    TimingRun("B12", "two full nodes", 2, 104, 208, 0.074, "observed timing"),
    TimingRun("B13", "spread-vs-pack check", 4, 52, 208, 0.060, "observed timing"),
    TimingRun("B14", "three full nodes", 3, 104, 312, 0.052, "observed timing"),
    TimingRun("B15", "current full setup", 4, 104, 416, 0.0416, "first observed timing"),
    TimingRun("B16", "eight-node spread check", 8, 52, 416, 0.035, "observed timing"),
    TimingRun("B17", "dense eight-node test", 8, 78, 624, 0.025, "observed timing"),
    TimingRun("B18", "full eight-node test", 8, 104, 832, 0.032, "observed timing"),
    TimingRun("B19", "full five-node test", 5, 104, 520, 0.037, "observed timing"),
    TimingRun("B20", "full six-node test", 6, 104, 624, 0.039, "observed timing"),
    TimingRun("B21", "full seven-node test", 7, 104, 728, 0.034, "observed timing"),
    TimingRun("B22", "six-node spread check", 6, 70, 420, 0.043, "observed timing"),
]


NODE_COLORS = {
    1: "#4e79a7",
    2: "#f28e2b",
    3: "#e15759",
    4: "#76b7b2",
    5: "#59a14f",
    6: "#edc948",
    7: "#b07aa1",
    8: "#ff9da7",
}


def compute_metrics(runs: list[TimingRun]) -> list[TimingMetric]:
    baseline = next(run for run in runs if run.run_id == "B01").sec_per_it
    return [
        TimingMetric(
            run=run,
            speedup=baseline / run.sec_per_it,
            efficiency=(baseline / run.sec_per_it) / run.ranks,
            iterations_per_second=1.0 / run.sec_per_it,
        )
        for run in runs
    ]


def fmt_sec(value: float) -> str:
    if value < 0.05:
        return f"{value:.4f}".rstrip("0").rstrip(".")
    if value < 0.1:
        return f"{value:.3f}".rstrip("0").rstrip(".")
    return f"{value:.2f}".rstrip("0").rstrip(".")


def fmt(value: float, digits: int = 2) -> str:
    return f"{value:.{digits}f}"


def metric_by_id(metrics: list[TimingMetric], run_id: str) -> TimingMetric:
    return next(metric for metric in metrics if metric.run.run_id == run_id)


def percent_faster(reference_sec: float, candidate_sec: float) -> float:
    return 100.0 * (reference_sec - candidate_sec) / reference_sec


def annotate_selected(
    ax: plt.Axes,
    metrics: list[TimingMetric],
    y_name: str,
    offsets: dict[str, tuple[int, int]],
) -> None:
    selected = set(offsets)
    for metric in metrics:
        if metric.run.run_id not in selected:
            continue
        y_value = metric.run.sec_per_it if y_name == "sec_per_it" else getattr(metric, y_name)
        ax.annotate(
            metric.run.run_id,
            xy=(metric.run.ranks, y_value),
            xytext=offsets[metric.run.run_id],
            textcoords="offset points",
            fontsize=8,
            bbox={"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.72},
        )


def plot_sec_per_iteration(metrics: list[TimingMetric]) -> Path:
    output = FIGURE_DIR / "memo_11_sec_per_iteration.png"
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    sorted_metrics = sorted(metrics, key=lambda metric: metric.run.ranks)
    ax.plot(
        [metric.run.ranks for metric in sorted_metrics],
        [metric.run.sec_per_it for metric in sorted_metrics],
        color="#2f3542",
        linewidth=1.2,
        alpha=0.55,
        zorder=1,
    )
    for nodes, color in NODE_COLORS.items():
        group = [metric for metric in metrics if metric.run.nodes == nodes]
        if not group:
            continue
        ax.scatter(
            [metric.run.ranks for metric in group],
            [metric.run.sec_per_it for metric in group],
            s=[72 if metric.run.run_id == "B17" else 42 for metric in group],
            color=color,
            edgecolor="#1f2933",
            linewidth=0.7,
            label=f"{nodes} node{'s' if nodes != 1 else ''}",
            zorder=2,
        )
    annotate_selected(
        ax,
        metrics,
        "sec_per_it",
        {
            "B01": (5, 6),
            "B10": (7, 8),
            "B11": (7, -17),
            "B15": (-30, -18),
            "B16": (12, 10),
            "B17": (-8, 10),
            "B20": (8, -6),
        },
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Total MPI ranks")
    ax.set_ylabel("SEC/IT")
    ax.set_title("Case B observed seconds per model iteration")
    ax.grid(True, which="both", linewidth=0.45, alpha=0.35)
    ax.legend(ncol=4, fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)
    return output


def plot_speedup(metrics: list[TimingMetric]) -> Path:
    output = FIGURE_DIR / "memo_11_speedup.png"
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    sorted_metrics = sorted(metrics, key=lambda metric: metric.run.ranks)
    ranks = [metric.run.ranks for metric in sorted_metrics]
    speedups = [metric.speedup for metric in sorted_metrics]
    ax.plot(ranks, speedups, color="#2f3542", linewidth=1.2, alpha=0.55, zorder=1)
    ax.plot([1, max(ranks)], [1, max(ranks)], "--", color="#9aa5b1", linewidth=1.0, label="Ideal linear")
    for nodes, color in NODE_COLORS.items():
        group = [metric for metric in metrics if metric.run.nodes == nodes]
        if not group:
            continue
        ax.scatter(
            [metric.run.ranks for metric in group],
            [metric.speedup for metric in group],
            s=[72 if metric.run.run_id == "B17" else 42 for metric in group],
            color=color,
            edgecolor="#1f2933",
            linewidth=0.7,
            label=f"{nodes} node{'s' if nodes != 1 else ''}",
            zorder=2,
        )
    annotate_selected(
        ax,
        metrics,
        "speedup",
        {
            "B01": (7, 8),
            "B10": (7, 6),
            "B11": (7, -13),
            "B15": (8, -22),
            "B16": (8, 10),
            "B17": (8, 8),
            "B18": (-35, -18),
            "B20": (-35, -19),
            "B22": (8, 8),
        },
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Total MPI ranks")
    ax.set_ylabel("Speedup vs B01")
    ax.set_title("Case B observed speedup relative to one MPI rank")
    ax.grid(True, which="both", linewidth=0.45, alpha=0.35)
    ax.legend(ncol=3, fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(output, dpi=300)
    plt.close(fig)
    return output


def latex_metric_row(metric: TimingMetric) -> str:
    run = metric.run
    layout = f"{run.nodes} x {run.tasks_per_node}"
    return (
        f"{run.run_id} & {layout} & {run.ranks} & {fmt_sec(run.sec_per_it)} "
        f"& {fmt(metric.speedup, 1)} & {fmt(100.0 * metric.efficiency, 1)}\\% \\\\"
    )


def write_summary_tex(metrics: list[TimingMetric], sec_plot: Path, speedup_plot: Path) -> None:
    fastest = min(metrics, key=lambda metric: metric.run.sec_per_it)
    current = metric_by_id(metrics, "B15")
    best_416 = min(
        [metric for metric in metrics if 400 <= metric.run.ranks <= 430],
        key=lambda metric: metric.run.sec_per_it,
    )
    b10 = metric_by_id(metrics, "B10")
    b11 = metric_by_id(metrics, "B11")
    b12 = metric_by_id(metrics, "B12")
    b13 = metric_by_id(metrics, "B13")
    b15 = metric_by_id(metrics, "B15")
    b16 = metric_by_id(metrics, "B16")
    b17 = metric_by_id(metrics, "B17")
    b18 = metric_by_id(metrics, "B18")
    b20 = metric_by_id(metrics, "B20")

    focus_ids = ["B15", "B16", "B17", "B18", "B20"]
    focus_rows = "\n".join(latex_metric_row(metric_by_id(metrics, run_id)) for run_id in focus_ids)

    spread_rows = "\n".join(
        [
            (
                f"104 ranks & B10, 1 x 104 & B11, 2 x 52 "
                f"& {fmt(percent_faster(b10.run.sec_per_it, b11.run.sec_per_it), 1)}\\% faster \\\\"
            ),
            (
                f"208 ranks & B12, 2 x 104 & B13, 4 x 52 "
                f"& {fmt(percent_faster(b12.run.sec_per_it, b13.run.sec_per_it), 1)}\\% faster \\\\"
            ),
            (
                f"416 ranks & B15, 4 x 104 & B16, 8 x 52 "
                f"& {fmt(percent_faster(b15.run.sec_per_it, b16.run.sec_per_it), 1)}\\% faster \\\\"
            ),
            (
                f"624 ranks & B20, 6 x 104 & B17, 8 x 78 "
                f"& {fmt(percent_faster(b20.run.sec_per_it, b17.run.sec_per_it), 1)}\\% faster \\\\"
            ),
        ]
    )

    fastest_vs_current = percent_faster(current.run.sec_per_it, fastest.run.sec_per_it)
    b18_penalty = 100.0 * (b18.run.sec_per_it - b17.run.sec_per_it) / b17.run.sec_per_it

    content = rf"""
% Auto-generated by Analysis/analyze_memo_11_speed.py.
% Do not edit this file directly; update the Python timing data and rerun.

\subsection*{{Observed Timing Results}}

The observed SEC/IT values are stored in \texttt{{Analysis/analyze\_memo\_11\_speed.py}}.
That script computes the speedup relative to B01, parallel efficiency, same-rank
spread/pack comparisons, and the figures below.

\begin{{figure}}[htbp]
\centering
\includegraphics[width=0.92\linewidth]{{memo_11_sec_per_iteration.png}}
\caption{{Observed case-B SEC/IT versus total MPI ranks. Lower values are faster; B17 is the fastest observed point.}}
\end{{figure}}

\begin{{figure}}[htbp]
\centering
\includegraphics[width=0.92\linewidth]{{memo_11_speedup.png}}
\caption{{Observed speedup relative to the B01 one-rank baseline. The dashed line marks ideal linear speedup.}}
\end{{figure}}

\begin{{table}}[htbp]
\centering
\caption{{Key observed timing configurations generated from the analysis script.}}
\begin{{tabular}}{{lrrrrr}}
\toprule
Run ID & Layout & Ranks & SEC/IT & Speedup & Efficiency \\
\midrule
{focus_rows}
\bottomrule
\end{{tabular}}
\end{{table}}

\begin{{table}}[htbp]
\centering
\caption{{Same-rank or near-diagnostic node-spread comparisons.}}
\begin{{tabular}}{{llll}}
\toprule
Rank case & Packed/reference & Spread/candidate & SEC/IT change \\
\midrule
{spread_rows}
\bottomrule
\end{{tabular}}
\end{{table}}

\subsection*{{Analysis Interpretation}}

The fastest observed configuration is \textbf{{{fastest.run.run_id}}}, with
{fastest.run.nodes} nodes, {fastest.run.tasks_per_node} ranks per node,
{fastest.run.ranks} total ranks, and {fmt_sec(fastest.run.sec_per_it)} SEC/IT.
This is a {fmt(fastest.speedup, 1)}x speedup over B01 and is
{fmt(fastest_vs_current, 1)}\% faster than the current four-node B15 layout.

For the near-416-rank operating point, \textbf{{{best_416.run.run_id}}} is the
best observed choice: {fmt_sec(best_416.run.sec_per_it)} SEC/IT compared with
{fmt_sec(current.run.sec_per_it)} SEC/IT for B15.  The matched-rank checks show
that spreading ranks across more nodes improves timing for 104, 208, 416, and
624 total ranks in these observations.

The eight-node full-occupancy run B18 is slower than B17 by
{fmt(b18_penalty, 1)}\%, even though B18 uses more ranks.  This suggests that
the present sweet spot is not maximum rank count, but moderate per-node
occupancy on eight nodes.  The practical recommendation is to use B17 when the
queue allows eight nodes, use B16 for a 416-rank comparison run, and repeat
B16--B18 if final production timing needs tighter uncertainty bounds.
""".lstrip()
    SUMMARY_TEX.write_text(content, encoding="ascii")


def main() -> None:
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    metrics = compute_metrics(RUNS)
    sec_plot = plot_sec_per_iteration(metrics)
    speedup_plot = plot_speedup(metrics)
    write_summary_tex(metrics, sec_plot, speedup_plot)
    fastest = min(metrics, key=lambda metric: metric.run.sec_per_it)
    print(f"Wrote {SUMMARY_TEX.relative_to(ROOT_DIR)}")
    print(f"Wrote {sec_plot.relative_to(ROOT_DIR)}")
    print(f"Wrote {speedup_plot.relative_to(ROOT_DIR)}")
    print(
        f"Fastest: {fastest.run.run_id} "
        f"({fmt_sec(fastest.run.sec_per_it)} SEC/IT, {fmt(fastest.speedup, 1)}x speedup)"
    )


if __name__ == "__main__":
    main()
