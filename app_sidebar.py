"""Module for rendering the elements in the sidebar of the Streamlit app."""

if __name__ == "__main__":
    print("This module is not intended to be run directly.")
    print("Please import this as a module on the main Streamlit page (app.py).")
    exit()

from typing import Literal
from collections import OrderedDict
import streamlit as st
import pandas as pd

import app_helper

# ============================================================================ #

def render_tf_summary(tfclass_row: pd.Series, length: int, disprot_regions: pd.DataFrame):
    selected_genus_num: str = tfclass_row["genus_num"]
    selected_uniprot: str = tfclass_row["uniprot_acc"]
    selected_genus_name: str = tfclass_row["genus_name"]
    disprot_id = disprot_regions.iloc[0]["region_id"].split("r")[0] if not disprot_regions.empty else "N/A"
    disprot_regions_unique = len(disprot_regions[["start", "end"]].drop_duplicates())

    st.header(f":material/select_check_box: Selected TF details")
    st.table({
        ":material/link: UniProt accession": f"[{selected_uniprot}](https://www.uniprot.org/uniprotkb/{selected_uniprot}/entry)",
        ":material/numbers: Genus number": selected_genus_num,
        ":material/immunology: Genus name": selected_genus_name,
        ":material/straighten: Length": f"{length} residues",
        ":material/error_med: DisProt": f"[{disprot_id}](https://disprot.org/{disprot_id}) ({disprot_regions_unique} distinct regions)" if not disprot_regions.empty else "N/A",
    }, border="horizontal", width="content", hide_header=True)

# ============================================================================ #

def render_pattern_selector(selected_patterns_df: pd.DataFrame, selected_genus_num: str) -> str | None:
    """Render pattern selector, and return the selected pattern's regex."""

    st.header(":material/view_timeline: Patterns to Display")

    patterns = selected_patterns_df.groupby(["Instance_accession", "Regex"])
    # print(*patterns.count().columns)

    pattern_options = (
        selected_patterns_df[["Instance_accession", "Regex"]].drop_duplicates()
        if not selected_patterns_df.empty
        else selected_patterns_df
    )

    if "selected_pattern_acc" not in st.session_state: st.session_state.selected_pattern_acc = ""

    pattern_options__state = st.dataframe(
        pattern_options,
        hide_index=True,
        key=f"pattern_options-{selected_genus_num}",
        column_config={
            "Instance_accession": "Instance acc.",
            "Regex": "Regex",
        },
        selection_mode="single-row",
        on_select="rerun",
    )

    # TODO: remove patterns with only 1 occurance
    pattern_options__sel_row_indices = app_helper.get_df_selected_rows(pattern_options__state)
    if len(pattern_options__sel_row_indices) == 1:
        selected_pattern__sel_row_index = pattern_options__sel_row_indices[0]
        selected_pattern__str: str = pattern_options.iloc[selected_pattern__sel_row_index]["Regex"]

        st.write(f"Selected pattern:\n\n`{selected_pattern__str}`")

        return selected_pattern__str

    return None
