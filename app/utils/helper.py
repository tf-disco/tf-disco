"""Module for helper functions."""

if __name__ == "__main__":
    print("This module is not intended to be run directly.")
    print("Please import this as a module on the main Streamlit page (app.py).")
    exit()

from random import randint
from typing import Literal
from collections import OrderedDict
import re

import streamlit as st
import streamlit.elements.arrow
import numpy as np
from numpy.typing import NDArray
import pandas as pd

from . import data_loading

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

def parse_dbd_ranges(dbd_range_str: str) -> list[tuple[int, int]]:
    """Parse the DBD range string into a list of (start, end) tuples.

    :param dbd_range_str: The DBD range string, e.g. "10-20,30-40"
    :return: A list of (start, end) tuples, e.g. [(10, 20), (30, 40)]
    """

    ranges = []
    for part in dbd_range_str.split(","):
        part = part.strip()
        if not part: continue
        try:
            start, end = map(int, part.split("-"))
            ranges.append((start, end))
        except ValueError:
            continue
    return ranges

#endregion

# ============================================================================ #

#region Render
FilterBy = Literal["Superclass", "Class", "Family", "Subfamily", "Name"]
# dict mapping `option` to its display name
FILTER_OPTIONS: OrderedDict[FilterBy, str] = OrderedDict([
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

def __st_get_option(key: str, default=None):
    """Get a Streamlit option, with a default value if the option is not set."""
    try: return st.get_option(key)
    except KeyError: return default

def render_sequence(
    sequence: str, pattern: str|None=None,
    disorder_scores: pd.Series|NDArray|None=None,
    dbd_ranges: pd.Series|NDArray|None=None,
    theme: Literal["dark", "light"]|None=None,
    ) -> str:
    """
    Get the chunked sequence block in HTML. Also highlights any matches of the
    provided pattern.

    :param sequence: Amino acid sequence.
    :param pattern: The regex pattern to highlight in the sequence. If None, no highlighting is done.
    :param disorder_scores: A list of (start, end) tuples indicating regions to overline in the sequence.
    :param dbd_ranges: A list of (start, end) tuples indicating regions to underline in the sequence.
    :param theme: The theme to use for rendering. If `None`, the current Streamlit theme is used.
    """

    GAP = "10px"

    THEME = theme or st.context.theme.type or "light" # ! TODO: explore CSS media query instead

    COLOR_BACKGROUND = __st_get_option(f"theme.{THEME}.backgroundColor", "#eff1f5")

    COLOR_MATCH_1_BG = __st_get_option(f"theme.{THEME}.greenColor", "#40a02b")
    COLOR_MATCH_1_FG = COLOR_BACKGROUND

    COLOR_MATCH_2_BG = f"color-mix(in srgb, {COLOR_MATCH_1_BG} 50%, transparent)"
    COLOR_MATCH_2_FG = "unset" # default (black/white as per theme)

    COLOR_UNDERLINE = __st_get_option(f"theme.{THEME}.redColor", "#d20f39")

    # this array will hold the index of the match for each residue, or 0 if not matched
    # use it to determine alternating colors
    residues_matched = np.zeros(len(sequence))
    residues_overline = dbd_ranges.astype(bool) if dbd_ranges is not None else np.zeros(len(sequence), dtype=bool)
    residues_underline = disorder_scores.astype(bool) if disorder_scores is not None else np.zeros(len(sequence), dtype=bool)

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
            html_chunk_char = ""

            if residues_matched[index] > 0:
                class_name = "c1" if residues_matched[index] % 2 == 1 else "c2"
                last_child = " l" if (index_within_chunk == len(chunk_chars) - 1) and (residues_matched[index+1:index+2].sum() == residues_matched[index]) else ""
                html_chunk_char = f"""<span class="{class_name}{last_child}">{char}</span>"""
            else:
                html_chunk_char = char

            if residues_overline[index]:
                html_chunk_char = f"""<span class="o">{html_chunk_char}</span>"""

            if residues_underline[index]:
                html_chunk_char = f"""<span class="u">{html_chunk_char}</span>"""

            html_chunk += html_chunk_char

        pos = f"{(chunk_start + len(chunk_chars))}"
        # render chunk position only if the number fits the width of the chunk
        pos = f"{pos}" if len(chunk_chars) >= len(pos) else ""

        html_chunks.append(f"""
<span class="chunk">
    <span class="pos">{pos}<br/></span>
    <span class="seq">{html_chunk}</span>
</span>
""")

    style = f"""
<style>
    .rendered-sequence-container {{display: flex; flex-wrap: wrap; gap: {GAP}; align-items: flex-start; padding: 0 0 1.5em 1.5em;}}
    .rendered-sequence-container .chunk {{display: flex; flex-direction: column; align-items: flex-end; font-family: monospace;}}

    .rendered-sequence-container .chunk .pos {{font-size: 0.8rem; opacity: 0.6; user-select: none;}}

    .rendered-sequence-container .chunk .seq .c1 {{background-color: {COLOR_MATCH_1_BG}; color: {COLOR_MATCH_1_FG}; font-weight: bold; text-decoration-color: {COLOR_MATCH_1_BG}; position: relative;}}
    .rendered-sequence-container .chunk .seq .c1.l:after {{content: ''; display: block; position: absolute; top: 0; right: -{GAP}; width: {GAP}; height: 100%; background-color: {COLOR_MATCH_1_BG};}}
    .rendered-sequence-container .chunk .seq .c2 {{background-color: {COLOR_MATCH_2_BG}; color: {COLOR_MATCH_2_FG}; font-weight: bold; text-decoration-color: {COLOR_MATCH_2_BG}; position: relative;}}
    .rendered-sequence-container .chunk .seq .c2.l:after {{content: ''; display: block; position: absolute; top: 0; right: -{GAP}; width: {GAP}; height: 100%; background-color: {COLOR_MATCH_2_BG};}}
    .rendered-sequence-container .chunk .seq .u {{text-decoration: underline 2pt; text-underline-position: under; text-decoration-color: {COLOR_UNDERLINE}}}
    .rendered-sequence-container .chunk .seq .o {{text-decoration: overline 2pt;}}
</style>
"""

    rendered = f"""{style}<div class="rendered-sequence-container">{"\n".join(html_chunks)}</div>"""
    return rendered

@st.cache_data(persist="disk", show_spinner=True)
def render_sequences_multiple_html(
    genus_nums: list[str], pattern: str,
    selected_pattern_df: pd.DataFrame,
    theme: Literal["dark", "light"]|None=None,
):
    """`render_sequence()` with same pattern and multiple sequences. Returns a
    full HTML document."""

    html_list: list[str] = [
        f"<h1>Matches for pattern {pattern}</h1>",
    ]
    sequence_dict = data_loading.init().sequence_dict

    for genus_num in genus_nums:
        tf_info = sequence_dict.get(genus_num, data_loading.TFInfo(Uniprot_Acc="", Genus_Name="", Sequence=""))
        observed_matches_count = selected_pattern_df[selected_pattern_df["Genus_Num"] == genus_num]["Observed"].sum()

        html_list.append(f"""
<h3>Matches for <a href="https://uniprot.org/entry/{tf_info.Uniprot_Acc}" target="_blank">{tf_info.Uniprot_Acc}</a> | genus <a href="https://tf-disco.com/tf_view?genus_num={genus_num}" target="_blank">{genus_num}</a> | {tf_info.Genus_Name}</h3>
<p>Observed Matches: {observed_matches_count}</p>
{render_sequence(sequence=tf_info.Sequence, pattern=pattern, theme=theme)}
""")

    return f"""<html><head><meta charset="UTF-8"><title>Matches for pattern {pattern}</title></head><body>{"<hr/>".join(html_list)}</body></html>"""

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
