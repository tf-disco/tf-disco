import streamlit as st

import app_data_loading

#region Config
st.set_page_config(
    page_title="Janus",
    # page_title="...",
    page_icon="🧬",
    initial_sidebar_state="expanded",
    menu_items={
        # TODO: add about info (and repo link)
        "Get Help": "https://placecats.com/1920/108?fit=fill",
        "Report a bug": "https://placecats.com/1920/108?fit=fill",
        "About": "# This is a header. This is an *extremely* cool app!"
    },
    layout="wide",
)

st.session_state["data"] = app_data_loading.init()

#endregion

st.logo("🧬")

pages = [
    st.Page("pages/home.py", icon=":material/home:", title="Home"),
    st.Page("pages/tf_browser.py", icon=":material/view_list:", title="TF Browser"),
    st.Page("pages/tf_view.py", icon=":material/visibility:", visibility="hidden", title="TF Viewer"),
    st.Page("pages/tf_compare.py", icon=":material/visibility:", visibility="hidden", title="TF Compare"),
    st.Page("pages/pattern_explorer.py", icon=":material/regular_expression:", title="Pattern Explorer"),
    st.Page("pages/about.py", icon=":material/info:", title="About"),
]
pages = st.navigation(pages, position="top")
pages.run()

# TODO: Comparison page for comparing multiple TFs
# TODO: Filtering Algo for Patterns
