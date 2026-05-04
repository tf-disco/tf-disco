import streamlit as st
import pandas as pd
import time

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import app_data_loading
import app_helper

#region Config
st.set_page_config(
    page_title="Pattern Explorer | Janus",
    page_icon="🧬",
    initial_sidebar_state="expanded",
)

st.title("🧬 Janus Pattern Explorer", text_alignment="center", anchor=False)
st.caption("Explore and analyze transcription factor patterns", text_alignment="center")

with st.sidebar:
    st.empty()

#endregion

#region Data loading
data: app_data_loading.InitResult = st.session_state["data"]
disprot_df = data.disprot_df
tfclasses_df = data.tfclasses_df
genus_num_name_map = data.genus_num_name_map
load_disorder_scores = app_data_loading.load_disorder_scores
load_pattern = app_data_loading.load_pattern

#endregion

with st.expander(f"Select Patterns", expanded=True):

    #region Filter controls
    # TODO: send catalog_filt to filter instead of tfclasses_df? so that users can further filter down the already filtered results.
    catalog_filt = tfclasses_df
    option, filter_by, filter_disprot = app_helper.render_filter_controls(tfclasses_df)

    if filter_disprot:
        catalog_filt = catalog_filt[catalog_filt["disprot_available"] == "✅"] # 😭😭😭

    if option and filter_by:
        st.write(f"Filtering by {filter_by}: `{'`, `'.join(option)}`")

        catalog_filt = catalog_filt[
            catalog_filt[filter_by.lower()].isin(option)
        ]

    # Since we're passing tfclasses_df_filtered to st.dataframe, and st.dataframe
    # only accepts and returns *ROW NUMBERS*, we reset the index to count from 0.
    catalog_filt = catalog_filt.reset_index(drop=True)

    catalog_filt__before: pd.DataFrame = st.session_state.get("catalog_filt__before", catalog_filt)
    catalog_filt__is_refiltered = len(catalog_filt__before) != len(catalog_filt) or (catalog_filt__before["genus_num"] != catalog_filt["genus_num"]).any()
    st.session_state["catalog_filt__before"] = catalog_filt

    #endregion

    st.session_state["cart"] = st.session_state.get("cart", set())
    cart: set[str] = st.session_state["cart"]

    #region Pattern Catalog render
    # Add all variables to this list, such that when they change, it should trigger a re-render of the df:
    catalog_filt__dependencies: str = "-".join(map(str, [
        filter_by,
        option,
        filter_disprot,
        st.session_state.get("sidebar_cart__trigger", 0),
    ]))
    catalog_filt__state = st.dataframe(
        data=catalog_filt[["genus_num", "uniprot_acc", "genus_name", "disprot_available"]],
        selection_default={"selection": {"rows":
            catalog_filt.index[catalog_filt["genus_num"].isin(cart)].tolist()
        }},
        column_config={
            "genus_num": "Genus Number",
            "uniprot_acc": "UniProt Accession",
            "genus_name": "Genus Name",
            "disprot_available": "DisProt data available?",
        },
        height="stretch",
        hide_index=True,
        key=catalog_filt__dependencies,
        selection_mode="multi-row",
        on_select="rerun",
    )

    # if user re-filtered, then don't modify cart.
    # if only selection changed, that means user actively clicked items, hence update the cart.
    if not catalog_filt__is_refiltered:
        catalog_filt__sel_row_indices = app_helper.get_df_selected_rows(catalog_filt__state)
        catalog_filt__genus_num = catalog_filt["genus_num"]

        for i, genus_num in enumerate(catalog_filt__genus_num.tolist()):
            is_selected = i in catalog_filt__sel_row_indices

            if (is_selected) and (genus_num not in cart): cart.add(genus_num)
            elif (genus_num in cart) and (not is_selected): cart.remove(genus_num)


    #endregion

#TODO: time cart sidebar
#region Cart render
with st.sidebar:
    st.header(":material/shopping_cart: Cart contents")

    if not cart:
        st.info(":material/info: Cart is empty.")

    # TODO: update to multiselect in cart df, let users select tiems to remove.
    # add a "remove selected" button
    else:
        st.success(f":material/remove_shopping_cart: Items in cart: {len(cart)}")
        # st.info(":material/remove_shopping_cart: Select an item to remove it.")
        sidebar_cart = tfclasses_df[["genus_num", "uniprot_acc", "genus_name"]][tfclasses_df["genus_num"].isin(cart)]
        sidebar_cart__state = st.dataframe(
            data=sidebar_cart,
            selection_default={"selection": {"rows": list(range(len(sidebar_cart)))}},
            column_config={
                "genus_num": "Genus #",
                "uniprot_acc": "UniProt",
                "genus_name": "Genus Name",
            },
            hide_index=True,
            key=f"sidebar-cart-{cart}",
            on_select="rerun",
            # selection_mode="single-row",
            selection_mode="multi-row",
        )

        sidebar_cart__sel_row_indices = app_helper.get_df_selected_rows(sidebar_cart__state)
        sidebar_cart__sel_items = sidebar_cart["genus_num"].iloc[sidebar_cart__sel_row_indices] # type:ignore

        # if sidebar_cart__sel_row_indices:
        #     for sidebar_cart__row_to_remove in sidebar_cart__sel_row_indices:
        #         genus_num: str = sidebar_cart["genus_num"].iloc[sidebar_cart__row_to_remove]
        #         cart.remove(genus_num)
        #         # st.toast(f"Removed {genus_num}.")
        #     st.session_state["sidebar_cart__trigger"] = st.session_state.get("sidebar_cart__trigger", 0) + 1
        #     st.rerun() # to make the main catalog re-render


        with st.container(width="stretch", horizontal=False, horizontal_alignment="right"):
            if len(sidebar_cart__sel_items) == 1:
                if st.button("View selected item", icon=":material/visibility:"):
                    st.switch_page("pages/tf_view.py", query_params={"genus_num": list(cart)[0]})
            else:
                if st.button("Compare selected items", icon=":material/visibility:", disabled=len(sidebar_cart__sel_items)==0):
                    st.switch_page("pages/tf_compare.py", query_params={"genus_nums": ",".join(cart)})

            if st.button(f"Remove selected item" + ("" if len(sidebar_cart__sel_items)==1 else "s") + " from cart", icon=":material/delete:", type="primary", disabled=len(sidebar_cart__sel_items)==0):
                for sidebar_cart__row_to_remove in sidebar_cart__sel_row_indices:
                    genus_num: str = sidebar_cart["genus_num"].iloc[sidebar_cart__row_to_remove]
                    cart.remove(genus_num)
                    # st.toast(f"Removed {genus_num}.")
                st.session_state.update({"sidebar_cart__trigger": st.session_state.get("sidebar_cart__trigger", 0) + 1})
                st.rerun()

#endregion

st.divider()

st.header(":material/regular_expression: Pattern", anchor=False)
