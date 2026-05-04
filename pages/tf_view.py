import streamlit as st
import pandas as pd

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import app_data_loading
import app_sidebar
import app_helper
import app_graph

#region Config
st.set_page_config(
    page_title="Viewer | Janus",
    page_icon="🧬",
    initial_sidebar_state="expanded",
)
#endregion

#region Data loading
data: app_data_loading.InitResult = st.session_state["data"]
disprot_df = data.disprot_df
tfclasses_df = data.tfclasses_df
genus_num_name_map = data.genus_num_name_map
load_disorder_scores = app_data_loading.load_disorder_scores
load_pattern = app_data_loading.load_pattern

#endregion

@st.dialog(
    title="No TF selected",
    icon=":material/error:",
    dismissible=False
)
def blank_view():
    st.error("No TF selected. Please select a TF from the browser.", icon=":material/error:")
    st.page_link("pages/tf_browser.py", label="Go to TF Browser", icon=":material/arrow_back:", icon_position="left")

if "genus_num" not in st.query_params:
    blank_view()
    st.stop()

# region Selection details

selected_row = None

if "genus_num" in st.query_params:
    genus_num_from_url = st.query_params["genus_num"]
    matches = tfclasses_df[tfclasses_df["genus_num"] == genus_num_from_url]
    if not matches.empty:
        selected_row = matches.iloc[0]
elif "selected_tf" in st.session_state and "filtered_tfclasses_df" in st.session_state:
    selected_tf = st.session_state["selected_tf"]
    filtered_tfclasses_df = st.session_state["filtered_tfclasses_df"]
    if selected_tf["selection"]["rows"]:  # type: ignore
        selected_row_index: int = selected_tf["selection"]["rows"][0]  # type: ignore
        selected_row = filtered_tfclasses_df.iloc[selected_row_index]

if selected_row is not None:
    selected_genus_num: str = selected_row["genus_num"]
    selected_uniprot: str = selected_row["uniprot_acc"]
    selected_genus_name: str = selected_row["genus_name"]

    st.query_params["genus_num"] = selected_genus_num

    selected_disorder_df = load_disorder_scores(selected_genus_num)
    selected_patterns_df = load_pattern(selected_genus_num)
    length = len(selected_disorder_df)

    disprot_regions = disprot_df[disprot_df["acc"] == selected_uniprot]
    disprot_regions = disprot_regions.sort_values(by="end", ascending=False).sort_values(by="start")
    disprot_id = disprot_regions.iloc[0]["region_id"].split("r")[0] if not disprot_regions.empty else "N/A"

    # Sidebar details
    with st.sidebar:
        st.divider()
        app_sidebar.render_tf_summary(selected_row, length, disprot_regions)
        st.divider()
        selected_pattern = app_sidebar.render_pattern_selector(selected_patterns_df, selected_genus_num)

    # TF metrics cards
    st.header(f"{selected_genus_name}", anchor="selected_tf") #: {selected_genus_num} | {selected_genus_name} | [{selected_uniprot}](https://www.uniprot.org/uniprotkb/{selected_uniprot}/entry)")
    col1, col2, col3 = st.columns(3, border=True)
    col1.metric(":material/link: *UniProt Accession*", f"**[{selected_uniprot}](https://www.uniprot.org/uniprotkb/{selected_uniprot}/entry)**")
    col2.metric(":material/error_med: *DisProt ID*", f"**[{disprot_id}](https://disprot.org/{disprot_id})**" if not disprot_regions.empty else "N/A")
    col3.metric(":material/straighten: *Sequence length*", f"**{length} residues**")

    # Sequence section
    with st.expander("Show Sequence", icon=":material/genetics:", expanded=True):
        sequence = selected_disorder_df['Residue'].str.cat(sep='')
        sequence_html = app_helper.render_sequence(sequence=sequence, pattern=selected_pattern)
        sequence_fasta = f">{selected_genus_num}_{selected_uniprot}_{selected_genus_name}\n{sequence}\n"
        st.write(f"{sequence_html}", unsafe_allow_html=True)
        st.divider()
        with st.container(horizontal=True, vertical_alignment="center", horizontal_alignment="distribute"):
            st.text("FASTA:")
            st.download_button(
                label=":material/download: Download FASTA",
                data=sequence_fasta,
                file_name=f"{selected_genus_num}_{selected_uniprot}_{selected_genus_name}.fasta",
                mime="text/plain"
            )
        st.code(sequence_fasta, language="plaintext", line_numbers=False, wrap_lines=False)

        # TODO: weblogo
        # sigh

    st.divider()
    st.space(2)

    # Scores section
    score_cols: list[app_graph.ScoreName] = [
        "Aiupred",
        "AiupredLatest-Disorder",
        "AiupredLatest-Binding",
        "AiupredLatest-Linker",
        "HI",
        "Fldpnn-Disorder",
        "Fldpnn-ProteinBinding",
        "Fldpnn-DnaBinding",
        "Fldpnn-RnaBinding",
        "Fldpnn-Linker",
        "Metapredict-Disorder",
    ]

    default_cols: list[app_graph.ScoreName] = [
        "AiupredLatest-Disorder",
        "HI",
        "Fldpnn-Disorder",
        "Metapredict-Disorder",
    ]
    # persist selection across navigations
    st.session_state["scoresToDisplay"] = st.session_state.get("scoresToDisplay", default_cols)

    with st.expander("Score Plots", icon=":material/area_chart:", expanded=True):
        st.subheader("Score Plots:")
        display_scores = st.pills(
            label="Select scores to display:",
            options=score_cols,
            default=st.session_state["scoresToDisplay"],
            selection_mode="multi",
            key="scoresToDisplay",
        )
        # TODO: also add sequence in x axis
        fig = app_graph.create_scores_plotly(
            length=length,
            scores_list=[
                app_graph.make_score_renderable(score_name, selected_disorder_df)
                for score_name in score_cols
                if score_name in display_scores
            ],
            disprot_regions=disprot_regions,
        )
        st.plotly_chart(fig)

    # Disprot regions section
    if not disprot_regions.empty:
        with st.expander("DisProt Regions for selected TF", icon=":material/error_med:"):
            st.subheader(f"DisProt Regions for [{disprot_id}](https://disprot.org/{disprot_id}):")
            st.dataframe(
                disprot_regions[["region_id", "start", "end"]],
                # key=f"disprot_regions-{selected_uniprot}", # uncomment if you make the rows selectable later on
                hide_index=True,
            )

    # Patterns section
    patterns_reqd = load_pattern(selected_genus_num)
    # TODO: check if instances count works?
    # patterns_reqd = patterns_reqd[patterns_reqd["#Instances"] > 1]
    patterns_reqd = patterns_reqd[patterns_reqd.groupby("Regex")["Instance_accession"].transform("count") > 1]
    patterns_reqd = (
        patterns_reqd
            .groupby(['Instance_accession', 'Regex', "motif_start_in_query", "motif_end_in_query"], as_index=False)
            .first()
    )

    if not patterns_reqd.empty:
        with st.expander("ELM Patterns for selected TF", icon=":material/view_timeline:"):
            st.dataframe(
                patterns_reqd[["Instance_accession", "Regex", "Instances (Matched Sequence)" ,"motif_start_in_query", "motif_end_in_query"]],
                # patterns_reqd[["Instance_accession", "Regex", "Matched_sequence" ,"motif_start_in_query", "motif_end_in_query"]],
                # key=f"patterns-{selected_genus_num}", # uncomment if you make the rows selectable later on
                hide_index=False,
            )

    st.set_page_config(page_title=f"{selected_genus_name} | Janus")

# endregion
