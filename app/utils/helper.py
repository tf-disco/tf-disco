"""Module for helper functions."""

if __name__ == "__main__":
    print("This module is not intended to be run directly.")
    print("Please import this as a module on the main Streamlit page (app.py).")
    exit()

from random import randint
from typing import cast, Literal, Final
from collections import OrderedDict
import re

import streamlit as st
import streamlit.elements.arrow
import numpy as np
from numpy.typing import NDArray
import pandas as pd

from . import data_loading
from . import graph

# ============================================================================ #

#region Non-render
def convert_ranges_to_array(ranges: list[tuple[int, int]], length: int):
    """Convert a list of (start, end) tuples into a binary array of the given
    length. Provide 1-based start and end positions. The returned array will be
    0-based."""

    arr = np.zeros(length, dtype=int)
    for start, end in ranges:
        arr[start-1:end] = 1
    return arr

def calculate_disprot_perc(ranges: pd.DataFrame, seq_len: int) -> float:
    """For the given list of DisProt regions, calculate the disordered
    percentage.

    :param ranges: A DataFrame with `Start` and `End` columns.
    :param seq_len: The length of the full sequence.
    """

    if seq_len <= 0: return 0.0

    return convert_ranges_to_array(
        ranges=list(ranges[["Start", "End"]].astype(int).itertuples(index=False, name=None)),
        length=seq_len
    ).sum() / seq_len

def get_df_selected_rows(state: streamlit.elements.arrow.DataframeState) -> list[int]:
    """Shorthand for `state["selection"]["rows"]`, but with safety."""

    if "selection" not in state: return []
    if "rows" not in state["selection"]: return []
    return state["selection"]["rows"]

def parse_range_str(range_str: str) -> list[tuple[int, int]]:
    """Parse any range string into a list of (start, end) tuples.

    :param range_str: The range string, e.g. "10-20,30-40"
    :return: A list of (start, end) tuples, e.g. [(10, 20), (30, 40)]
    """

    ranges = []
    for part in range_str.split(","):
        part = part.strip()
        if not part: continue
        try:
            start, end = map(int, part.split("-"))
            ranges.append((start, end))
        except ValueError:
            continue
    return ranges

def genus_sort_key(genus_num: str, handle_zero: bool=True):
    """Split the given genus-number string into its parts. Useful for sorting.

    :param handle_zero: It is observed in the TFClass website that `0` is placed
        *after* the other numbers, instead of at the start. Set to `True` to
        respect this sorting logic.
    """

    key = list(map(int,
        genus_num.splitlines()[0].split(".")
    ))
    if handle_zero:
        key = [
            float("inf") if part == 0 else part
            for part in key
        ]
    return tuple(key)


#endregion

# ============================================================================ #

#region Render
FilterBy = Literal["Superclass", "Class", "Family", "Subfamily", "Name"]
# dict mapping `option` to its display name
FILTER_OPTIONS: Final[OrderedDict[FilterBy, str]] = OrderedDict([
    ("Superclass", "Superclass"),
    ("Class", "Class"),
    ("Family", "Family"),
    ("Subfamily", "Subfamily"),
    ("Name", "Genus number (or) UniProt Accession (or) Genus name"),
])

def render_search_dropdown(classes_df: pd.DataFrame, filter_by: FilterBy, st_key: str) -> list[str]:
    return st.multiselect(
        label=f"Select a {filter_by}:",
        label_visibility="collapsed",
        options=classes_df[filter_by].dropna().unique(),
        placeholder=f"Search for {FILTER_OPTIONS.get(filter_by, filter_by)}",
        key=st_key,
    )

def render_filter_controls(tfclasses_df: pd.DataFrame):
    """Filter the `tfclasses_df` based on user-input `filter_by` and `option`.
    Returns the following:
    1. `selected_option`
    2. `filter_by`
    3. `disprot_filter`
    """

    tab1, tab2 = st.tabs(["Basic search", "Advanced search"], on_change="rerun")

    disprot_filter = st.toggle(label="View only entries with DisProt data available", value=False, key="disprot_filter", help="This will filter the catalog to only show entries that have DisProt data available, which is required for pattern generation. (Currently, only a subset of entries have DisProt data available.)")

    with tab1:
        filter_by1 = "Name"
        selected_option1 = render_search_dropdown(tfclasses_df, filter_by1, "option1")

    with tab2:
        filter_by2: FilterBy = st.segmented_control(
            label="Filter by:",
            options=[k for k in FILTER_OPTIONS if k!="Name"],
            format_func=lambda k: FILTER_OPTIONS.get(k, f"{k}"),
            default="Subfamily",
            required=True,
            selection_mode="single",
            help="List all TFs which belong to a specific superclass/class/family/subfamily (refer to the [TFClass Resource](http://tfclass.bioinf.med.uni-goettingen.de/index.jsf) for details on TF classification)",
        ) # type: ignore because narrowing from str to `FilterBy`
        selected_option2 = render_search_dropdown(tfclasses_df, filter_by2, "option2")

    if tab1.open:
        return (selected_option1, filter_by1, disprot_filter)
    else: # tab2.open:
        return (selected_option2, filter_by2, disprot_filter)
    # else:
    #     return ([], "", disprot_filter)

# ============================================================================ #

# colors for match highlights
COLOR_MATCH_1_BG: Final = "#5f9ea0"
COLOR_MATCH_1_FG: Final = "white"
COLOR_MATCH_2_BG: Final = "#5f9ea07f"
COLOR_MATCH_2_FG: Final = "unset" # default (black/white as per theme)

def render_sequence(
    sequence: str,
    pattern: str|None=None,
    pattern_tooltip: bool=False,
    overline_residues: pd.Series|NDArray|None=None,
    overline_color: str="black",
    overline_tooltip: str|None=None,
    underline_residues: pd.Series|NDArray|None=None,
    underline_color: str="#d20f39",
    underline_tooltip: str|None=None,
    ) -> str:
    """
    Get the chunked sequence block in HTML.

    There are three types of annotations available:
    1. Highlights any matches of the provided pattern.
    2. Overline on residues (per-residue basis)
    3. Underline on residues (per-residue basis)

    Tooltips may be added on residues that have at least one annotation on them.
    To configure these tooltips, provide static strings for `pattern_tooltip`,
    `overline_tooltip`, and `underline_tooltip`.

    :param sequence: Amino acid sequence.
    :param pattern: The regex pattern to highlight in the sequence. If None, no highlighting is done.
    :param pattern_tooltip: Whether to include a tooltip for the residues within matches of the pattern.
    :param overline_residues: An array of `1`s and `0`s indicating positions to overline in the sequence.
    :param overline_color: Overline color. As given by `render_annotations_popover()`.
    :param overline_tooltip: Static tooltip for overlined residues.
    :param underline_residues: An array of `1`s and `0`s indicating positions to underline in the sequence.
    :param underline_color: Underline color. As given by `render_annotations_popover()`.
    :param underline_tooltip: Static tooltip for underlined residues.
    """

    GAP: Final = "10px"

    # this array will hold the index of the match for each residue, or 0 if not matched
    # use it to determine alternating colors
    residues_matched = np.zeros(len(sequence))
    overline_residues = overline_residues.astype(bool) if overline_residues is not None else np.zeros(len(sequence), dtype=bool)
    underline_residues = underline_residues.astype(bool) if underline_residues is not None else np.zeros(len(sequence), dtype=bool)

    if pattern:
        try:
            for i, match in enumerate(re.finditer(pattern, sequence), start=1):
                residues_matched[match.start():match.end()] = i
        except:
            pass

    html_chunks: list[str] = []

    for chunk_start in range(0, len(sequence), 10):
        chunk_chars = sequence[chunk_start:chunk_start+10]

        # Build HTML for this chunk with highlighting
        html_chunk = ""
        for index_within_chunk, char in enumerate(chunk_chars):
            index = chunk_start + index_within_chunk
            classes = []
            tooltips = []

            if residues_matched[index] > 0:
                classes.append("c1" if residues_matched[index] % 2 == 1 else "c2")
                if (index_within_chunk == len(chunk_chars) - 1) and (residues_matched[index+1:index+2].sum() == residues_matched[index]): classes.append("l")
                if pattern_tooltip: tooltips.append(f"▶️ Matches pattern: {pattern}")

            if overline_residues[index]:
                classes.append("o")
                if overline_tooltip: tooltips.append(f"▶️ {overline_tooltip}")

            if underline_residues[index]:
                classes.append("u")
                if underline_tooltip: tooltips.append(f"▶️ {underline_tooltip}")

            if tooltips: tooltips = [f"[{index+1}] {char}", ""] + tooltips
            tooltips_attr = f' title="{"&#013;".join(tooltips)}"' if tooltips else ""
            html_chunk += f"""<span class="{' '.join(classes)}"{tooltips_attr}>{char}</span>""" if classes else char

        pos = f"{(chunk_start + len(chunk_chars))}"
        # render chunk position only if the number fits the width of the chunk
        pos = f"{pos}" if len(chunk_chars) >= len(pos) else ""

        html_chunks.append(
            '    <span class="chunk">'
                f'<span class="pos">{pos}<br></span>'
                f'<span class="seq">{html_chunk}</span>'
            '</span>'
        )

    container_id = f"rendered-sequence-container-{randint(1, 100000)}"

    style = f"""
<style>
    #{container_id} {{
        --gap: {GAP}; --color-match-1-bg: {COLOR_MATCH_1_BG}; --color-match-1-fg: {COLOR_MATCH_1_FG}; --color-match-2-bg: {COLOR_MATCH_2_BG}; --color-match-2-fg: {COLOR_MATCH_2_FG}; --color-underline: {underline_color}; --color-overline: {overline_color};
        display: flex; flex-wrap: wrap; gap: var(--gap); align-items: flex-start; line-height: 1.6; padding: 0 0 1.5em 1.5em;
    }}
    #{container_id} .chunk {{display: flex; flex-direction: column; align-items: flex-end; font-family: monospace;}}

    #{container_id} .chunk .pos {{font-size: 0.8rem; opacity: 0.6; user-select: none;}}

    #{container_id} .chunk .seq > span {{padding: 1pt 0pt; background-clip: content-box;}}
    #{container_id} .chunk .seq .c1 {{background-color: var(--color-match-1-bg); color: var(--color-match-1-fg); font-weight: bold; text-decoration-color: var(--color-match-1-fg); position: relative;}}
    #{container_id} .chunk .seq .c1.l:after {{content: ''; display: block; box-sizing: border-box; position: absolute; top: 0; right: calc(-1 * var(--gap)); width: var(--gap); height: 100%; background-color: var(--color-match-1-bg); padding: 1pt 0pt; background-clip: content-box;}}
    #{container_id} .chunk .seq .c2 {{background-color: var(--color-match-2-bg); color: var(--color-match-2-fg); font-weight: bold; text-decoration-color: var(--color-match-2-fg); position: relative; /*box-shadow: blue 0pt 0pt 5pt 3pt; clip-path: inset(-8pt 0 -8pt 0);*/}}
    #{container_id} .chunk .seq .c2.l:after {{content: ''; display: block; box-sizing: border-box; position: absolute; top: 0; right: calc(-1 * var(--gap)); width: var(--gap); height: 100%; background-color: var(--color-match-2-bg); padding: 1pt 0pt; background-clip: content-box;}}
    #{container_id} .chunk .seq .u {{border-bottom: 2pt solid var(--color-underline);}}
    #{container_id} .chunk .seq .o {{border-top: 2pt solid var(--color-overline);}}
</style>
"""

    rendered = f"""{style}<div id="{container_id}">\n{"\n".join(html_chunks)}\n</div>"""
    return rendered

@st.cache_data()
def render_sequences_multiple_html(
    genus_nums: list[str], pattern: str|None=None,
    selected_pattern_df: pd.DataFrame|None=None,
    tooltips: bool=True,
    overline_type: Literal["", "DBD", "DPR"]="",
    overline_color: str="black",
    underline_type: graph.ScoreName|Literal[""]="",
    underline_color: str="#d20f39",
    render_pattern_heading: bool=True,
    make_relative_link: bool=False,
):
    """`render_sequence()` with same pattern and multiple sequences. Returns a
    full HTML document.

    :param genus_nums: List of genus numbers to render sequences for.
    :param pattern: The regex pattern to highlight in the sequences.
    :param selected_pattern_df: DataFrame containing the observed matches for the pattern, with columns `"Genus_Num"`, `"Observed"`.
    :param tooltips: Whether to include a tooltip for the residues within matches of the pattern, overlined residues, and/or underlined residues.
    :param overline_type: Overline selected annotation. One of `""`, `"DBD"`, `"DPR"`. As given by `render_annotations_popover()`.
    :param overline_color: Overline color. As given by `render_annotations_popover()`.
    :param underline_type: Underline selected annotation. Either `""` or the name of a disorder predictor (e.g. `"Aiupred-Disorder"`). As given by `render_annotations_popover()`.
    :param underline_color: Underline color. As given by `render_annotations_popover()`.
    :param render_pattern_heading: Whether to render a heading for the pattern at the top of the page.
    :param make_relative_link: Whether to make the genus number link relative (for use in the pattern explorer page) or absolute (for use in the help page).
    """

    html_list: list[str] = [f"<h1>Matches for pattern {pattern}</h1>"] if render_pattern_heading else []
    data = data_loading.init()

    for genus_num in genus_nums:
        tf_info = data.sequence_dict.get(genus_num, data_loading.TFInfo(Uniprot_Acc="", Genus_Name="", Sequence=""))
        tfclasses_row = data.tfclasses_df[data.tfclasses_df["Genus_Num"] == genus_num].iloc[0]
        observed_matches_count = selected_pattern_df[selected_pattern_df["Genus_Num"] == genus_num]["Observed"].sum() if selected_pattern_df is not None else None

        html = render_sequence(
            sequence=tf_info.Sequence, pattern=pattern,
            pattern_tooltip=tooltips,
            overline_residues=convert_ranges_to_array(parse_range_str(tfclasses_row["DBD"] if overline_type == "DBD" else (tfclasses_row["DPR"] if overline_type == "DPR" else "")), len(tf_info.Sequence)),
            overline_color=overline_color,
            overline_tooltip=(f"In {overline_type} region" if overline_type else None) if tooltips else None,
            underline_residues=data_loading.load_disorder_scores(genus_num)[underline_type] >= graph.SCORE_PROPERTIES[underline_type].threshold if underline_type else None,
            underline_color=underline_color,
            underline_tooltip=(graph.SCORE_PROPERTIES[underline_type].display_name if underline_type else None) if tooltips else None,
        )
        html_list.append("".join([
            f'<h3>Matches for <a href="https://uniprot.org/entry/{tf_info.Uniprot_Acc}" target="_blank">{tf_info.Uniprot_Acc}</a> | genus <a href="{'' if make_relative_link else 'https://tf-disco.com'}/tf_view?genus_num={genus_num}" target="_blank">{genus_num}</a> | {tf_info.Genus_Name}</h3>',
            f'<p>Observed matches: {observed_matches_count}</p>{html}' if observed_matches_count is not None else '',
        ]))

    return f"""<html><head><meta charset="UTF-8"><title>Matches for pattern {pattern}</title></head><body>\n\n{"\n\n<hr/>\n\n".join(html_list)}\n\n</body></html>"""

def render_annotations_popover():
    """Render a popover for selecting the annotations (underline and overline)
    to render on a sequence. Use the return value of this function to determine
    the inputs for `render_sequence()`.

    :return:
    A tuple containing 4 values:
        1. Overline selected annotation. One of `""`, `"DBD"`, `"DPR"`
        2. Overline color
        3. Underline selected annotation. Either `""` or the name of a disorder predictor (e.g. `"Aiupred-Disorder"`)
        4. Underline color
    """

    V_HIDE: Final = ":material/visibility_off: Hide"

    with st.popover(label=":material/format_ink_highlighter: Annotations"):
        overline_type = cast(Literal["", "DBD", "DPR"], st.pills(
            label="Overline annotations:",
            selection_mode="single",
            required=True,
            default="",
            options=["", "DBD", "DPR"],
            format_func=lambda value: {"": V_HIDE, "DBD": f':color[:material/format_overline:]{{foreground="{graph.COLOR_DBD}"}} DBD', "DPR": f':color[:material/format_overline:]{{foreground="{graph.COLOR_DPR}"}} DPR'}[value]
        ))
        if overline_type == "DBD":
            overline_color = graph.COLOR_DBD
            overline_name = "DNA Binding Domains"
        elif overline_type == "DPR":
            overline_color = graph.COLOR_DPR
            overline_name = "Droplet-Promoting Regions"
        else:
            overline_color = ""
            overline_name = "DNA Binding Domains or Droplet-Promoting Regions"

        st.write(f'Displays an **:color[overline]{{foreground="{overline_color}"}}** above the residues which are part of {overline_name}.\n\n---')

        underline_type = cast(graph.ScoreName|Literal[""], st.pills(
            label="Underline annotations - Disorder predictors:",
            selection_mode="single",
            required=True,
            default="",
            options=["", "Aiupred-Disorder", "Fldpnn-Disorder", "Metapredict-Disorder"],
            format_func=lambda score_name: f':color[:material/format_underlined:]{{foreground="{graph.SCORE_PROPERTIES[score_name].color}"}} {graph.SCORE_PROPERTIES[score_name].display_name}' if score_name else V_HIDE,
        ))

        if underline_type:
            score_properties = graph.SCORE_PROPERTIES[underline_type]
            underline_color = score_properties.color
            underline_name = score_properties.display_name.partition('-')[0].strip()
        else:
            underline_color = ""
            underline_name = "the predictor selected above"

        st.write(f'Displays an **:color[underline]{{foreground="{underline_color}"}}** below the residues predicted to be disordered by {underline_name}.')

    return (
        overline_type,
        cast(str, overline_color),
        underline_type,
        cast(str, underline_color),
    )

# ============================================================================ #

def render_vagueness_penalty_slider(key: str|None=None):
    """Render a vagueness penalty slider, and return the selected value.

    This is its own function because it's used in multiple places, and so that
    the huge ass help text is consistent."""

    return st.slider(
        "Vagueness Penalty", min_value=1.0, max_value=10.0, value=2.25, step=0.25,
        help="""
Vagueness is a measure of how specific a pattern is, with higher values
indicating "vague" patterns, and lower values indicating more "specific"
patterns (this value is internally calculated, and not shown in the table).

Using this slider, a :primary[**penalty**] may be applied to :primary[**highly
vague**] patterns. Higher penalty allows the user to focus on patterns that are
more likely to be biologically meaningful.

Set to 1.0 to not apply any penalty.

Note that this is only used to rank the patterns for convenience purposes, and
it doesn't affect the actual pattern matching, and it doesn't discard any
patterns.

For a detailed explanation, please visit the [Help page](/help#vagueness).
""",
        key=key,
    )

# ============================================================================ #

def container_indent(*args, **kwargs):
    """Render a st.container with an indent."""

    key = f"{randint(1,100000)}"
    st.html(f"<style>.st-key-{key} {{margin-left: 2em;}}</style>")
    return st.container(*args, key=key, **kwargs)

#endregion
