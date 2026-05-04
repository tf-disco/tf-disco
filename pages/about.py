import streamlit as st


#region Config
st.set_page_config(
    page_title="Janus | About",
    page_icon="🧬",
    initial_sidebar_state="collapsed",
)

#endregion

st.title("About 🧬 Janus!", text_alignment="center", anchor=False)
st.caption("A Consolidated Database of Human Transcription Factors Information", text_alignment="center")

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
