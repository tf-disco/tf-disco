"""Module for creating the main score plots."""

if __name__ == "__main__":
    print("This module is not intended to be run directly.")
    print("Please import this as a module on the main Streamlit page (app.py).")
    exit()

from dataclasses import dataclass
from typing import NamedTuple, Any, Literal

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ============================================================================ #

ScoreName = Literal[
    "HI", "Aiupred-Disorder", "Aiupred-Binding",
    "Aiupred-Linker", "Fldpnn-Disorder", "Fldpnn-ProteinBinding",
    "Fldpnn-DnaBinding", "Fldpnn-RnaBinding", "Fldpnn-Linker",
    "Metapredict-Disorder"
]

@dataclass
class ScoreTemplate():
    name: ScoreName
    display_name: str
    threshold: float
    is_downward: bool
    yaxis_range: tuple[float, float]
    color: str
    fillcolor: str

@dataclass
class Score(ScoreTemplate):
    values: np.typing.NDArray[Any]

SCORE_PROPERTIES: dict[ScoreName, ScoreTemplate] = {
    "HI": ScoreTemplate("HI", "Hydropathy Index",
        0.0, True, (-3.0, 3.0), "red", "rgba(255, 0, 0, 0.3)"),
    "Aiupred-Disorder": ScoreTemplate("Aiupred-Disorder", "AIUPred - Disorder",
        0.5, False, (0.0, 1.0), "deepskyblue", "rgba(0, 191, 255, 0.3)"),
    "Aiupred-Binding": ScoreTemplate("Aiupred-Binding", "AIUPred - Binding",
        0.5, False, (0.0, 1.0), "deepskyblue", "rgba(0, 191, 255, 0.3)"),
    "Aiupred-Linker": ScoreTemplate("Aiupred-Linker", "AIUPred - Linker",
        0.5, False, (0.0, 1.0), "deepskyblue", "rgba(0, 191, 255, 0.3)"),
    "Fldpnn-Disorder": ScoreTemplate("Fldpnn-Disorder", "flDPnn - Disorder",
        0.3, False, (0.0, 1.0), "green", "rgba(0, 128, 0, 0.3)"),
    "Fldpnn-ProteinBinding": ScoreTemplate("Fldpnn-ProteinBinding", "flDPnn - Protein binding",
        0.5, False, (0.0, 1.0), "green", "rgba(0, 128, 0, 0.3)"),
    "Fldpnn-DnaBinding": ScoreTemplate("Fldpnn-DnaBinding", "flDPnn - DNA binding",
        0.5, False, (0.0, 1.0), "green", "rgba(0, 128, 0, 0.3)"),
    "Fldpnn-RnaBinding": ScoreTemplate("Fldpnn-RnaBinding", "flDPnn - RNA binding",
        0.5, False, (0.0, 1.0), "green", "rgba(0, 128, 0, 0.3)"),
    "Fldpnn-Linker": ScoreTemplate("Fldpnn-Linker", "flDPnn - Linker",
        0.5, False, (0.0, 1.0), "green", "rgba(0, 128, 0, 0.3)"),
    "Metapredict-Disorder": ScoreTemplate("Metapredict-Disorder", "Metapredict v3 - Disorder",
        0.5, False, (0.0, 1.0), "gold", "rgba(255, 215, 0, 0.3)"),
}

def make_score_renderable(score_name: ScoreName, df: pd.DataFrame) -> Score | None:
    """Generate a plot-ready object for plotting the given score name, extracted
    from the given TF df."""

    if score_name not in SCORE_PROPERTIES: return None
    if score_name not in df: return None

    result_template = SCORE_PROPERTIES[score_name]
    return Score(
        name=result_template.name,
        display_name=result_template.display_name,
        values=df[score_name].to_numpy(),
        threshold=result_template.threshold,
        is_downward=result_template.is_downward,
        yaxis_range=result_template.yaxis_range,
        color=result_template.color,
        fillcolor=result_template.fillcolor,
    )

def create_scores_plotly(sequence: str, scores_list: list[Score | None], dbd_ranges: list[tuple[int, int]], disprot_regions: pd.DataFrame):
    """
    Create a Plotly Figure for plotting the given list of scores.

    :param str sequence: the sequence to plot (for x-axis labels)
    :param list[Score|None] scores_list: list of scores you wanna plot. `None`s are skipped. Refer to `make_score_renderable()`.
    :param list[tuple[int, int]] dbd_ranges: list of (start, end) positions (1-based) for the DBD range.
    :param pd.DataFrame disprot_regions: df with columns: `Region_Id`, `Start`, `End`
    """

    length = len(sequence)
    x = np.arange(1, length + 1)

    scores = [s for s in scores_list if s is not None]
    is_disprot_available = len(disprot_regions) > 0

    # MARK: Initialize plot
    fig = make_subplots(
        rows=len(scores) + (1 if is_disprot_available else 0) + (1 if dbd_ranges else 0) + 1,
        cols=1,
        shared_xaxes=True,
        # vertical_spacing=(0.3 / len(scores)), # less spacing if no score rows
        row_heights=[
            *([10] * len(scores)),  # weight `10` for each big graph
            *([1] if is_disprot_available else []),  # small weight for DisProt row if it exists
            *([1] if dbd_ranges else []),  # small weight for DBD Ranges row if it exists
            2,
        ],
        subplot_titles=[
            *[s.display_name for s in scores],
            *(["DisProt"] if is_disprot_available else []),
            *(["DBD Ranges"] if dbd_ranges else []),
            # "Sequence",
        ],
    )

    # MARK: Each Disorder
    for row, score in enumerate(scores, start=1):
        # Plot main line
        fig.add_trace(
            go.Scatter(
                x=x,
                y=score.values,
                name=score.display_name, mode="lines", line=dict(color=score.color, width=2),
            ),
            row=row,
            col=1,
        )

        # Plot area under curve
        fig.add_trace(
            go.Scatter(
                x=x,
                y=np.repeat(score.threshold, length),
                mode="lines", line=dict(width=0), showlegend=False, hoverinfo="skip",
            ),
            row=row,
            col=1,
        )  # reference line for area-under-curve

        fig.add_trace(
            go.Scatter(
                x=x,
                y=(np.clip(score.values, a_min=score.threshold, a_max=None)
                   if not score.is_downward
                   else np.clip(score.values, a_min=None, a_max=score.threshold)),
                name=score.display_name,
                mode="lines", line=dict(color=score.color, width=0), fill="tonexty", fillcolor=score.fillcolor, showlegend=False, hoverinfo="skip",
            ),
            row=row,
            col=1,
        )  # area under curve

        # Configure axes
        fig.update_yaxes(
            range=[score.yaxis_range[0] - 0.05, score.yaxis_range[1] + 0.05],
            tickvals=[score.yaxis_range[0], score.threshold, score.yaxis_range[1]],  # [0, 0.2, 0.4, 0.6, 0.8, 1.0],
            showline=True,
            linewidth=1,
            linecolor="lightgrey",  # Draw left axis line
            row=row,
            col=1,
        )

    # MARK: DBD Ranges bars
    if dbd_ranges:
        for i, (start, end) in enumerate(dbd_ranges):
            start = max(1, int(start))
            end = min(length, int(end))
            if end < start:
                continue

            xs = list(range(start, end + 1))
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=[0] * len(xs),
                    mode="lines",
                    name="DBD Range",
                    line=dict(color="blue", width=15),
                    hovertemplate=f"DBD Range {i+1}<br>Start: {start}<br>End: {end}<br>Position: %{{x}}<extra></extra>",
                    showlegend=(i == 0),
                ),
                row=len(scores) + (1 if is_disprot_available else 0) + 1,
                col=1,
            )

        fig.update_yaxes(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            title_text="DBD",
            title_font=dict(size=12, color="grey"),
            title_standoff=0,
            showline=True,
            linewidth=1,
            range=[-0.1, 0.1],
            fixedrange=True,
            row=len(scores) + (1 if is_disprot_available else 0) + 1,
            col=1,
        )

    # MARK: DisProt bars
    if is_disprot_available: # DisProt track (bottom row) + build mask for overlaps
        disprot_mask = np.zeros(length, dtype=bool)
        for i, (region_id, start, end) in enumerate(disprot_regions[["Region_Id", "Start", "End"]].itertuples(index=False, name=None)):
            start = max(0, int(start))
            end = min(length - 1, int(end))
            if end < start:
                continue

            disprot_mask[start : end + 1] = True

            xs = list(range(start, end + 1))
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=[0] * len(xs),
                    mode="lines",
                    name="DisProt",
                    line=dict(color="purple", width=15),
                    # hoverinfo="skip",
                    hovertemplate=f"{region_id}<br>Start: {start}<br>End: {end}<br>Position: %{{x}}<extra></extra>",
                    showlegend=(i == 0),
                ),
                row=len(scores) + 1,
                col=1,
            )

        fig.update_yaxes(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            title_text="DisProt",
            title_font=dict(size=12, color="grey"),
            title_standoff=0,
            showline=True,
            linewidth=1,
            range=[-0.1, 0.1],
            fixedrange=True,
            row=len(scores) + 1,
            col=1,
        )

        # For each score row, add "in DisProt + above threshold" highlight segments
        for row, score in enumerate(scores, start=1):
            vals = np.asarray(score.values)
            n = min(len(vals), len(disprot_mask))

            hit_mask = ((vals[:n] >= score.threshold) if not score.is_downward else (vals[:n] <= score.threshold)) & disprot_mask[:n]
            if not np.any(hit_mask):
                continue

            idx = np.where(hit_mask)[0]
            segments = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)

            highlight_y = (
                score.yaxis_range[0] - 0.03
            )  # thin strip near the bottom of each score panel
            for seg in segments:
                s, e = int(seg[0]), int(seg[-1])
                fig.add_trace(
                    go.Scatter(
                        x=[s, e],
                        y=[highlight_y, highlight_y],
                        mode="lines",
                        line=dict(color=score.color, width=7),
                        name="",
                        showlegend=False,
                        hovertemplate=(
                            f"{score.display_name}<br>"
                            f"Segment: {s}-{e}<br>"
                            f"Condition: score >= {score.threshold} and in DisProt"
                        ),
                    ),
                    row=row,
                    col=1,
                )

    # MARK: AAs + position axis
    fig.add_trace(
        go.Scatter(
            x=np.arange(1, length + 1),
            y=[0] * length,
            mode="text",
            line=dict(color="red"),
            text=list(sequence),
            hoverinfo="skip",
            hovertemplate=f"Residue: %{{text}}<extra></extra>",
        ),
        row=len(scores) + (1 if is_disprot_available else 0) + 1 + (1 if dbd_ranges else 0),
        col=1,
    )

    fig.update_yaxes(
        showticklabels=False,
        showgrid=False,
        zeroline=False,
        title_text="",
        title_font=dict(size=12, color="grey"),
        title_standoff=0,
        showline=True,
        linewidth=1,
        range=[-0.5, 0.5],
        fixedrange=True,
        row=len(scores) + (1 if is_disprot_available else 0) + 1 + (1 if dbd_ranges else 0),
        col=1,
    )

    fig.update_xaxes(
        showgrid=False,
        showline=True,
        linewidth=1,
        row=len(scores) + (1 if is_disprot_available else 0) + 1 + (1 if dbd_ranges else 0),
        col=1,
    )

    # MARK: Overall config
    fig.update_layout(
        height=sum([
            125 * len(scores),
            50 if is_disprot_available else 0,
            50 if dbd_ranges else 0,
            50,
        ]),
        margin=dict(l=40, r=20, t=20, b=20),
        showlegend=False,
        hovermode="x unified",
        legend_traceorder="normal",
    )

    fig.update_xaxes(
        showspikes=True, spikemode="across", spikesnap="cursor", spikethickness=1
    )
    return fig
