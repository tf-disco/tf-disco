"""Module for helper functions."""

if __name__ == "__main__":
    print("This module is not intended to be run directly.")
    print("Please import this as a module on the main Streamlit page (app.py).")
    exit()

from typing import Literal
from collections import OrderedDict
import streamlit as st
import streamlit.elements.arrow
import pandas as pd
import re

# ============================================================================ #

FilterBy = Literal["superclass", "class", "family", "subfamily", "genus"]
# dict mapping `option` to its display name
FILTER_OPTIONS: OrderedDict[FilterBy, str] = OrderedDict([
    ("superclass", "Superclass"),
    ("class", "Class"),
    ("family", "Family"),
    ("subfamily", "Subfamily"),
    ("genus", "Genus number (or) UniProt Accession (or) Genus name"),
])

def render_search_dropdown(classes_df: pd.DataFrame, filter_by: FilterBy, key: str) -> list[str]:
    return st.multiselect(
        label=f"Select a {filter_by}:",
        label_visibility="collapsed",
        options=classes_df[filter_by].dropna().unique(),
        placeholder=f"Search for {FILTER_OPTIONS.get(filter_by, filter_by)}",
        key=key,
    )

def render_filter_controls(tfclasses_df: pd.DataFrame):
    """Filter the `tfclasses_df` based on user-input `filter_by` and `option`.
    Returns the following:
    1. `selected_option`
    2. `filter_by`
    3. `disprot_filter`
    """

    tab1, tab2 = st.tabs(["Basic search", "Advanced search"], on_change="rerun")

    disprot_filter = st.toggle(label="View only entries with DisProt data available", value=False, key="disprot_filter")

    with tab1:
        filter_by1 = "genus"
        selected_option1 = render_search_dropdown(tfclasses_df, filter_by1, "option1")

    with tab2:
        filter_by2: FilterBy = st.segmented_control(
            label="Filter by:",
            options=[k for k in FILTER_OPTIONS if k!="genus"],
            format_func=lambda k: FILTER_OPTIONS.get(k, f"{k}"),
            default="subfamily",
            required=True,
            selection_mode="single",
        ) # type: ignore because narrowing from str to `FilterBy`
        selected_option2 = render_search_dropdown(tfclasses_df, filter_by2, "option2")

    if tab1.open:
        return (selected_option1, filter_by1, disprot_filter)
    else: # tab2.open:
        return (selected_option2, filter_by2, disprot_filter)
    # else:
    #     return ([], "", disprot_filter)

# ============================================================================ #

def render_sequence(sequence: str, pattern: str | None) -> str:
    """
    Render the chunked sequence block in HTML.
    Also highlights any matches of the provided pattern

    :param str sequence: the full sequence to render
    :param str pattern: the regex pattern to highlight in the sequence (optional)
    :return str: the rendered HTML string
    """

    merged_chunks = ""
    # Find all match positions if pattern is provided
    highlighted_positions = set()
    if pattern:
        try:
            regex = re.compile(pattern)
            for match in regex.finditer(sequence):
                for pos in range(match.start(), match.end()):
                    highlighted_positions.add(pos)
        except:
            pass

    for i in range(0, len(sequence), 10):
        chunk = sequence[i:i+10]
        index = i + len(chunk)

        # Build HTML for this chunk with highlighting
        chunk_html = ""
        for j, char in enumerate(chunk):
            pos = i + j
            if pos in highlighted_positions:
                chunk_html += f"<span style='background-color: orchid; font-weight: bold;'>{char}</span>"
            else:
                chunk_html += char

        html_chunk = f"""
<div style='display: flex; flex-direction: column; align-items: flex-end;'>
    <span style='font-family: monospace; font-size: 0.8rem; opacity: 0.6; user-select: none;'>{index}</span>
    <span style='font-family: monospace;'>{chunk_html}</span>
</div>
        """
        merged_chunks += html_chunk

    rendered = f"""
<div style='display: flex; flex-wrap: wrap; gap: 10px; align-items: flex-start;'>
    {merged_chunks}
</div>
        """
    return rendered

# ============================================================================ #

def get_df_selected_rows(state: streamlit.elements.arrow.DataframeState) -> list[int]:
    """Shorthand for `state["selection"]["rows"]`, but with safety."""

    if "selection" not in state: return []
    if "rows" not in state["selection"]: return []
    return state["selection"]["rows"]
