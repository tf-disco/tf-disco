import streamlit as st

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import app_data_loading

#region Config
st.set_page_config(
    page_title="Janus | Home",
    page_icon="🧬",
    initial_sidebar_state="collapsed",
)

st.title("Welcome to 🧬 Janus!", text_alignment="center", anchor=False)
st.caption("A Consolidated Database of Human Transcription Factors Information", text_alignment="center")
st.divider()

#endregion

#region Data loading
data: app_data_loading.InitResult = st.session_state["data"]
disprot_df = data.disprot_df
tfclasses_df = data.tfclasses_df
genus_num_name_map = data.genus_num_name_map
load_disorder_scores = app_data_loading.load_disorder_scores
load_pattern = app_data_loading.load_pattern

#endregion

# region Search
with st.container(width="stretch", horizontal_alignment="center"):
    st.header(":material/search: Search for a TF", anchor=False)
    st.write("Search by UniProt accession, Genus number, or Genus name:")

    # TODO: try st.container ig
    col1, col2 = st.columns([2, 1], vertical_alignment="center", width="stretch")

    # with col1:
    searched_tf: str|None = (
        col1.selectbox(
            label="Search by UniProt accession, Genus number, or Genus name:",
            label_visibility="collapsed",
            options=tfclasses_df["genus_num"].dropna().unique(),
            format_func=lambda genus_num: genus_num_name_map.get(genus_num, genus_num), # type:ignore
            placeholder="Start typing...",
            index=None,
            # key="genus_num",
            # bind="query-params",
        )
    )

    col2.page_link(
        "pages/tf_view.py",
        icon=":material/search:", label="Go to Viewer",
        disabled=not searched_tf,
        query_params={
            "genus_num": searched_tf, # type:ignore
        }
    )
with st.container(width="stretch", horizontal=True):
    st.write("Or, try out an example:")
    examples = ["1.1.1.2.3", "1.1.7.1.1", "2.1.2.4.3"]
    for genus_num in examples:
        st.page_link(
            "pages/tf_view.py",
            icon=":material/search:", label=genus_num_name_map[genus_num],
            query_params={"genus_num": genus_num,}
        )

# endregion

st.divider()

# region Details
st.header(":material/info: What is Janus?", anchor=False)
st.markdown("""
    Janus is a comprehensive database and visualization tool for human
    transcription factors (TFs). It integrates data from multiple sources,
    including [UniProt](https://www.uniprot.org),
    [DisProt](https://disprot.org), and the [ELM Resource](https://elm.eu.org),
    to provide detailed information on selected TFs. Disorder and some relevent
    scores are precalculated from AIUPred, flDPnn, and Metapredict v3 and known
    patterns are retreived and stored from the above mentioned sources.

    Use the <a href="browser" target="_self" rel="noreferrer">Browser</a> page
    to explore the available TFs and their characteristics, or jump straight
    into the viewer for an in-depth look at a specific TF.

    More details can be found in the <a href="about" target="_self"
    rel="noreferrer">About</a> page.
""", unsafe_allow_html=True)
# endregion

st.divider()

# region Contact
st.header(":material/contact_mail: Contact", anchor=False)
st.markdown("""
    For inquiries, feedback, or contributions, please reach out to us at
    [f20220019@dubai.bits-pilani.ac.in](mailto:f20220019@dubai.bits-pilani.ac.in). We welcome your input and
    look forward to hearing from you!
""")

# endregion
