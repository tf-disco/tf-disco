"""Module for loading data from disk.

Call the `init()` function to get the required variables."""

if __name__ == "__main__":
    print("This module is not intended to be run directly.")
    print("Please import this as a module on the main Streamlit page (app.py).")
    exit()

import streamlit as st
from pathlib import Path
import pandas as pd
from dataclasses import dataclass
from collections.abc import Callable
from time import perf_counter

#region Configuration
PATH_ROOT = Path(__file__).parents[1]
PATH_DISPROT = PATH_ROOT / "Datasets/disprot_data.tsv"
PATH_TF_CLASSES = PATH_ROOT / "Datasets/tf_classes_with_seq.tsv"
DIR_TFS_DISORDER = PATH_ROOT / "5. More metrics/"
DIR_TFS_PATTERNS = PATH_ROOT / "7. ELM Patterns/"
# DIR_TFS_PATTERNS = PATH_ROOT / "8. Pattern Statistics/"
COLS_DONT_REPLACE_X = {
    "Position", "Residue", "DbdFlag",
    "Disprot", "DisprotRegionIds",
    "AiupredAgreesWithDisprot", "HICorrelatesWithDisprot",
}

#endregion

# =========================================================================== #

@st.cache_data(show_spinner="Loading DisProt and TFClasses data...")
def load_disprot_tfclasses_dfs() -> tuple[pd.DataFrame, pd.DataFrame, dict[str, str]]:
    """Loads and returns:
    1. DisProt DataFrame
    2. TFClasses DataFrame
    3. A mapping from Genus number to Genus name
    """

    tfclasses_df = pd.read_csv(PATH_TF_CLASSES, sep="\t")
    disprot_df = pd.read_csv(PATH_DISPROT, sep="\t")

    # for each row (tf), `genus_num_parts` is a list of ints.

    # 1. sort on genus number
    tfclasses_df["genus_num_parts"] = tfclasses_df["genus_num"].str.split(".").map(lambda parts: [int(part) for part in parts]) # type:ignore
    tfclasses_df = tfclasses_df.sort_values(by="genus_num_parts")

    # 2. create columns superclass, class, family, subfamily, genus:
    tfclasses_df["superclass"] = tfclasses_df["genus_num_parts"].map(lambda parts: ".".join(str(part) for part in parts[:1])) # type: ignore
    tfclasses_df["class"] = tfclasses_df["genus_num_parts"].map(lambda parts: ".".join(str(part) for part in parts[:2])) # type: ignore
    tfclasses_df["family"] = tfclasses_df["genus_num_parts"].map(lambda parts: ".".join(str(part) for part in parts[:3])) # type: ignore
    tfclasses_df["subfamily"] = tfclasses_df["genus_num_parts"].map(lambda parts: ".".join(str(part) for part in parts[:4])) # type: ignore
    tfclasses_df = tfclasses_df.drop("genus_num_parts", axis=1)

    # 3. create genus column
    tfclasses_df["genus"] = tfclasses_df[["genus_num", "genus_name", "uniprot_acc"]].apply(
        lambda row: f"{row['genus_num']} | {row['uniprot_acc']} | {row['genus_name']}",
        axis=1
    )

    # 4. create new column "disprot_available"
    tfclasses_df["disprot_available"] = tfclasses_df["uniprot_acc"].isin(disprot_df["acc"]).map({True: "✅", False: "❌"})

    genus_num_name_map = tfclasses_df[["genus_num", "genus"]].drop_duplicates().set_index("genus_num").to_dict()["genus"]

    return disprot_df, tfclasses_df, genus_num_name_map

# =========================================================================== #

# @st.cache_data
# def load_files() -> tuple[dict[str, pd.DataFrame], dict[str, pd.DataFrame]]:
#     """:return disorder_map, pattern_map:
#     1. a `dict` mapping filenames to the disorder scores TSV (DataFrame)
#     2. a `dict` mapping filenames to the ELM patterns TSV (DataFrame)
#     """

#     disorder_map: dict[str, pd.DataFrame] = {}
#     pattern_map: dict[str, pd.DataFrame] = {}

#     print(f"{perf_counter():0.03f} - load_files start")

#     # Get Disorder scores from TF_DIR
#     for path in DIR_TFS_DISORDER.glob("*.tsv"):
#         genus_num = path.stem.split("_", maxsplit=1)[0]
#         df = pd.read_csv(path, sep="\t")

#         metric_cols = list( set(df.columns) - COLS_DONT_REPLACE_X )
#         df[metric_cols] = df[metric_cols].apply(pd.to_numeric, errors="coerce")

#         disorder_map[genus_num] = df

#     print(f"{perf_counter():0.03f} - load_files Disorder scores done")

#     # Get Patterns from PATTERNS_DIR
#     for path in DIR_TFS_PATTERNS.glob("*.tsv"):
#         genus_num = path.stem.split("_", maxsplit=1)[0]
#         df = pd.read_csv(path, sep="\t")
#         # pos_cols = ["motif_start_in_query", "motif_end_in_query"]
#         # available_pos_cols = [c for c in pos_cols if c in df.columns]
#         # if available_pos_cols:
#             # df[available_pos_cols] = df[available_pos_cols].apply(pd.to_numeric, errors="coerce")
#         pattern_map[genus_num] = df

#     print(f"{perf_counter():0.03f} - load_files Patterns done")

#     return disorder_map, pattern_map

@st.cache_data(show_spinner=False)
def load_paths_disorder() -> list[Path]:
    """Return list of available Paths for disorder scores."""
    # TODO: move ts to session_state
    return list(DIR_TFS_DISORDER.glob("*.tsv"))

# TODO: replace with the patterns in the `8.Pattern Statistics` folder, which are generated from the our list.
@st.cache_data(show_spinner=False)
def load_paths_patterns() -> list[Path]:
    """Return list of available Paths for ELM patterns."""
    # TODO: move ts to session_state
    return list(DIR_TFS_PATTERNS.glob("*.tsv"))

def load_df(paths: list[Path], genus_num: str) -> pd.DataFrame:
    """Load the DataFrame for the given genus number, from the given list of
    Paths."""

    paths_filtered = [
        path
        for path in paths
        if path.name.startswith(f"{genus_num}_")
    ]

    if not paths_filtered:
        return pd.DataFrame()
    if len(paths_filtered) > 1:
        raise ValueError(f"Ambiguous genus number: {genus_num}")
    path = paths_filtered[0]
    df = pd.read_csv(path, sep="\t")
    return df

@st.cache_data(show_spinner="Loading Disorder scores...")
def load_disorder_scores(genus_num: str) -> pd.DataFrame:
    """Load the disorder scores for the given genus number.

    Columns (except the ones listed in `COLS_DONT_REPLACE_X`) are coerced to
    numeric datatype. Non-numeric values (e.g. `X`) are replaced with `nan`."""

    df = load_df(load_paths_disorder(), genus_num)
    metric_cols = list( set(df.columns) - COLS_DONT_REPLACE_X )
    df[metric_cols] = df[metric_cols].apply(pd.to_numeric, errors="coerce")

    return df

@st.cache_data(show_spinner="Loading ELM patterns...")
def load_pattern(genus_num: str) -> pd.DataFrame:
    """Load the ELM patterns for the given genus number."""

    df = load_df(load_paths_patterns(), genus_num)
    return df

@dataclass
class InitResult:
    # disorder_map: dict[str, pd.DataFrame]
    # patterns_map: dict[str, pd.DataFrame]
    # lload_disorder_scores = load_disorder_scores#: Callable[[str], pd.DataFrame]
    # lload_pattern = load_pattern#: Callable[[str], pd.DataFrame]
    disprot_df: pd.DataFrame
    tfclasses_df: pd.DataFrame
    genus_num_name_map: dict[str, str]

def init() -> InitResult:
    """Initialize variables. Note that only the rows corresponding to the TFs
    available on disk are returned."""

    genus_numbers = {
        path.stem.partition("_")[0]
        for path in load_paths_disorder() + load_paths_patterns()
    }

    disprot_df, tfclasses_df, genus_num_name_map = load_disprot_tfclasses_dfs()
    # disorder_map, patterns_map = load_files()

    # Keep only those rows which have corresponding files on disk.
    tfclasses_df = tfclasses_df[tfclasses_df["genus_num"].isin(genus_numbers)].reset_index(drop=True)
    disprot_df = disprot_df[disprot_df["acc"].isin(tfclasses_df["uniprot_acc"])].reset_index(drop=True)
    genus_num_name_map = {k:v for k,v in genus_num_name_map.items() if k in genus_numbers}

    return InitResult(
        # disorder_map=disorder_map,
        # patterns_map=patterns_map,
        # load_disorder_scores=load_disorder_scores,
        # load_pattern=load_pattern,
        disprot_df=disprot_df,
        tfclasses_df=tfclasses_df,
        genus_num_name_map=genus_num_name_map,
    )
