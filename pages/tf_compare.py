import streamlit as st

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
    title="No TFs selected",
    icon=":material/error:",
    dismissible=False
)
def blank_view():
    st.error("No TFs selected. Please select more than one TFs from the browser.", icon=":material/error:")
    st.page_link("pages/tf_browser.py", label="Go to TF Browser", icon=":material/arrow_back:", icon_position="left")

if "genus_nums" not in st.query_params:
    blank_view()
    st.stop()

# st.write(st.query_params)

# st.divider()

# =========================================================

selected_rows = None
genus_nums_from_url = st.query_params["genus_nums"].split(",") if "genus_nums" in st.query_params else []

matches = tfclasses_df[tfclasses_df["genus_num"].isin(genus_nums_from_url)]
selected_rows = matches.to_dict(orient="records")
