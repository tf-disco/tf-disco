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
        0.0, True, (-3.0, 3.0), "slateblue", "rgba(106, 90, 238, 0.3)"),
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
        0.5, False, (0.0, 1.0), "orange", "rgba(255, 165, 0, 0.3)"),
}
COLOR_DBD = "magenta"
COLOR_DISPROT = "red"

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

# ============================================================================ #

def create_scores_plotly(sequence: str, scores_list: list[Score | None], dbd_ranges: list[tuple[int, int]], disprot_regions: pd.DataFrame, compare_with_disprot: bool=True):
    """
    Create a Plotly Figure for plotting the given list of scores.

    :param str sequence: the sequence to plot (for x-axis labels)
    :param list[Score|None] scores_list: list of scores you wanna plot. `None`s are skipped. Refer to `make_score_renderable()`.
    :param list[tuple[int, int]] dbd_ranges: list of (start, end) positions (1-based) for the DBD range.
    :param pd.DataFrame disprot_regions: df with columns: `Region_Id`, `Start`, `End`
    :param bool compare_with_disprot: whether to render DisProt compare bars below each graph
    """

    CAP_LENGTH = 0.25

    length = len(sequence)
    x = np.arange(1, length + 1)

    scores = [s for s in scores_list if s is not None]
    is_disprot_available = len(disprot_regions) > 0

    #region Initialize plot
    fig = make_subplots(
        rows=len(scores) + (1 if is_disprot_available else 0) + (1 if dbd_ranges else 0) + 1,
        cols=1,
        shared_xaxes=True,
        row_heights=[
            *([10] * len(scores)), # weight `10` for each big graph
            *([1] if is_disprot_available else []), # small weight for DisProt row if it exists
            *([1] if dbd_ranges else []), # small weight for DBD Ranges row if it exists
            2,
        ],
        subplot_titles=[
            *[s.display_name for s in scores],
            *(["DBD Range" + ("s" if len(dbd_ranges) > 1 else "")] if dbd_ranges else []),
            *(["DisProt"] if is_disprot_available else []),
            # "Sequence",
        ],
        font=dict(size=18),
    )

    fig.update_layout(
        height=sum([
            125 * len(scores),
            50 if is_disprot_available else 0,
            50 if dbd_ranges else 0,
            50,
        ]),
        margin=dict(l=20, r=20, t=20, b=20),
        showlegend=False,
        hovermode="x unified",
        legend_traceorder="normal",
    )

    fig.update_xaxes(
        showspikes=True, spikemode="across", spikesnap="data", spikethickness=1
    )

    #endregion

    #region Each Score
    for row, score in enumerate(scores, start=1):
        # Plot main line

        # hover_color = score.color
        hover_color = f"""color-mix(in srgb, currentcolor, {score.color} calc(100% * (
            sign(sign({ f'{score.threshold} - %{{y: }}' if score.is_downward else f'%{{y: }} - {score.threshold}' }) + 0.5) / 2 + 0.5
        )))""" # basically "score.color above threshold, default color below threshold"

        fig.add_trace(
            go.Scatter(
                x=x,
                y=score.values, name=score.display_name,
                mode="lines", line=dict(color=score.color, width=2),
                hovertemplate=f"{score.display_name}: <span style='color: {hover_color}; font-family: monospace; font-weight: bold;'>%{{y:+}}</span><extra></extra>",
            ),
            row=row,
            col=1,
        )

        shading_x = []
        shading_y = []
        idx = np.where(~np.isnan(score.values) & (score.values >= score.threshold if not score.is_downward else score.values <= score.threshold))[0] # indices where the score is not NaN
        segments = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
        for seg in segments:
            if seg.size == 0: continue

            # the problem is that, around each segment, the shading cuts off
            # abruptly. instead, it would be nice if we can have the shading
            # extend till the score's graph intercepts the threshold line.

            # to find out how far it extends, we will simply linear interpolate
            # (since the score's graph is ... well, linear), and add that point
            # to the list of points. you get the idea.
            # if the score's graph doesn't intercept (maybe if we don't have a
            # value before the start, likewise if we don't have a value after
            # the end), then we'll just extend to `CAP_LENGTH` amount.

            # naming convention for the positions:
            # assume scores = [0.2,  0.6,  0.8, 0.9,  0.7,  0.3]
            #                        ~~~~~~~~~~~~~~~~~~~~        = (segment to shade)
            #                  ~~~                               = beforestart (may be NaN)
            #                      ~                             = capstart!! to calculate
            #                        ~~~                         = start
            #                                         ~~~        = end
            #                                             ~      = capend!! to calculate
            #                                               ~~~  = afterend (may be NaN)
            # more naming conventions:
            #   cap = the extra shading (to be) added before the start &
            #         after the end of the segment.
            #   i   = 0-based index in sequence = (position - 1) = (xcoord - 1)

            start_i           = int(seg[0])
            start_value       = float(score.values[start_i])
            beforestart_i     = start_i-1
            beforestart_value = float(score.values[beforestart_i]) if beforestart_i >= 0 else np.nan
            capstart_i = (
                (start_i - CAP_LENGTH)
                if ( np.isnan(beforestart_value) or np.sign(start_value - score.threshold) == np.sign(beforestart_value - score.threshold) )
                # else (score.threshold - start_value) * (beforestart_i - start_i) / (beforestart_value - start_value) + start_i # simple intercept formula of a line
                else (start_value - score.threshold) / (beforestart_value - start_value) + start_i # same as above, but simplified
            )
            capstart_value = np.nan_to_num(np.interp(capstart_i, np.arange(length), score.values, left=np.nan), nan=score.threshold)

            end_i          = int(seg[-1])
            end_value      = float(score.values[end_i])
            afterend_i     = end_i + 1
            afterend_value = float(score.values[afterend_i]) if afterend_i < length else np.nan
            capend_i = (
                (end_i + CAP_LENGTH)
                if ( np.isnan(afterend_value) or np.sign(end_value - score.threshold) == np.sign(afterend_value - score.threshold) )
                # else (score.end_value - threshold) * (afterend_i - end_i) / (afterend_value - end_value) + end_i # simple intercept formula of a line
                else (score.threshold - end_value) / (afterend_value - end_value) + end_i # same as above, but simplified
            )
            capend_value = np.nan_to_num(np.interp(capend_i, np.arange(length), score.values, right=np.nan), nan=score.threshold)


            # capstart point
            shading_x.extend([capstart_i+1])
            shading_y.extend([capstart_value])

            # score points
            shading_x.extend(seg+1)
            shading_y.extend(score.values[seg])

            # capend point
            shading_x.extend([capend_i+1])
            shading_y.extend([capend_value])

        # Plot area under curve
        fig.add_trace(
            go.Scatter(
                x=shading_x,
                y=np.repeat(score.threshold, len(shading_x)),
                mode="lines", line=dict(width=0), showlegend=False, hoverinfo="skip",
            ),
            row=row,
            col=1,
        )  # reference line for area-under-curve

        fig.add_trace(
            go.Scatter(
                x=shading_x,
                y=shading_y,
                name=score.display_name,
                mode="lines", line=dict(color=score.color, width=0), fill="tonexty", fillcolor=score.fillcolor, showlegend=False, hoverinfo="skip",
            ),
            row=row,
            col=1,
        )  # area under curve

        # Configure y-axis
        fig.update_yaxes(
            range=[score.yaxis_range[0] - 0.05, score.yaxis_range[1] + 0.05],
            tickvals=[score.yaxis_range[0], score.threshold, score.yaxis_range[1]],
            showline=True,
            linewidth=1,
            linecolor="lightgrey",  # Draw left axis line
            row=row,
            col=1,
        )

        fig.layout.annotations[row - 1].font.color = score.color # type: ignore

    #endregion

    #region DBD Ranges bars
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
                    name="DBD Range",
                    mode="lines", line=dict(color=COLOR_DBD, width=15),
                    hovertemplate=f"DBD: {start}-{end}<extra></extra>",
                    showlegend=(i == 0),
                ),
                row=len(scores) + 1,
                col=1,
            )

        fig.update_yaxes(
            showticklabels=False, showgrid=False, zeroline=False,
            title_text="DBD", title_font=dict(size=12, color="grey"), title_standoff=0,
            showline=True, linewidth=1,
            range=[-0.1, 0.1],
            fixedrange=True,
            row=len(scores) + 1,
            col=1,
        )

        fig.layout.annotations[len(scores) + 1 - 1].font.color = COLOR_DBD # type: ignore

    #endregion

    #region DisProt bars
    if is_disprot_available: # DisProt track (bottom row) + build mask for overlaps
        disprot_mask = np.zeros(length, dtype=bool)
        for i, (region_id, start_pos, end_pos) in enumerate(disprot_regions[["Region_Id", "Start", "End"]].itertuples(index=False, name=None)):
            start_pos = max(0, int(start_pos))
            end_pos = min(length - 1, int(end_pos))
            if end_pos < start_pos:
                continue

            disprot_mask[start_pos-1 : end_pos + 1] = True

            xs = np.arange(start_pos, end_pos + 1)
            fig.add_trace(
                go.Scatter(
                    x=xs,
                    y=[0] * len(xs),
                    mode="lines",
                    name="DisProt",
                    line=dict(color=COLOR_DISPROT, width=15),
                    # hoverinfo="skip",
                    hovertemplate=f"DisProt <span style='color: {COLOR_DISPROT}'>{region_id}</span>: {start_pos}-{end_pos}<extra></extra>",
                    showlegend=(i == 0),
                ),
                row=len(scores) + 1 + 1,
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
            row=len(scores) + 1 + 1,
            col=1,
        )

        fig.layout.annotations[len(scores) + 1 + 1 - 1].font.color = COLOR_DISPROT # type: ignore

        # For each score row, add "in DisProt + above threshold" highlight segments
        for row, score in enumerate(scores if compare_with_disprot else [], start=1):
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
                s, e = int(seg[0]+1), int(seg[-1]+1) # 1-based position
                xs = np.arange(s, e + 1)
                fig.add_trace(
                    go.Scatter(
                        x=xs,
                        y=[highlight_y] * len(xs),
                        mode="lines",
                        line=dict(color=score.color, width=7),
                        name="",
                        showlegend=False,
                        hovertemplate=(
                            # f"{score.display_name}<br>"
                            f"In range ({s}-{e}): "
                            f"<span style='color: {score.color}'>metric</span> {'≥' if not score.is_downward else '≤'} {score.threshold} and in <span style='color: {COLOR_DISPROT}'>DisProt</span>"
                            "<extra></extra>"
                        ),
                    ),
                    row=row,
                    col=1,
                )

                # add "caps" to the left and right of the highlight segment
                fig.add_trace(go.Scatter(x=[s-CAP_LENGTH, s], y=[highlight_y] * 2, mode="lines", line=dict(color=score.color, width=7), showlegend=False, hoverinfo="skip",), row=row, col=1)
                fig.add_trace(go.Scatter(x=[e, e+CAP_LENGTH], y=[highlight_y] * 2, mode="lines", line=dict(color=score.color, width=7), showlegend=False, hoverinfo="skip",), row=row, col=1)

    #endregion

    #region AAs + position axis
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

    #endregion

    return fig
