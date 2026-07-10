"""Module for storing constants, such as page and data paths."""

if __name__ == "__main__":
    print("This module is not intended to be run directly.")
    print("Please import this as a module on the main Streamlit page (app.py).")
    exit()

from pathlib import Path
from os import getenv
from base64 import b64encode

from dotenv import load_dotenv

# ============================================================================ #

# Don't include "Browser" suffix here
APP_NAME = "TF-DISCO"
"""The name of the application, used for page titles and other places."""

PATH_ROOT = Path(__file__).parents[2]
"""Path to the root of project."""

def GET_APP_COMMIT_HASH() -> str:
    """Get the commit hash of the current version of the app."""
    try:
        __PATH_HEAD = PATH_ROOT / ".git/HEAD"
        __HEAD = __PATH_HEAD.read_text().strip() if __PATH_HEAD.exists() else ""
        return (PATH_ROOT / ".git/" / __HEAD[4:].strip()).read_text().strip() if __HEAD.startswith("ref:") else __HEAD
    except:
        return ""


#region Environment variables
if getenv("IS_DOCKER", "0") == "0":
    load_dotenv(PATH_ROOT / ".env")

__DATASET_PATH_OVERRIDE = getenv("DATASET_PATH_OVERRIDE")
ENV_DATASET_PATH_OVERRIDE = (PATH_ROOT / __DATASET_PATH_OVERRIDE).resolve() if __DATASET_PATH_OVERRIDE else None

__MUSCLE_PATH = getenv("MUSCLE_PATH")
ENV_MUSCLE_PATH = (PATH_ROOT / __MUSCLE_PATH).resolve() if __MUSCLE_PATH else None

#endregion


#region Page paths
PATH_ASSET_LOGO = PATH_ROOT / "assets/logo.png"
"""Path to the logo, 534x534."""

ASSET_LOGO_DATAURI = f"data:image/png;base64,{b64encode(PATH_ASSET_LOGO.read_bytes()).decode('utf-8')}"
"""Data URI for the logo, used for embedding the logo in markdown."""

PATH_PAGE_HOME = PATH_ROOT / "app/home.py"
"""Path to the **Home** page."""

PATH_PAGE_TF_BROWSER = PATH_ROOT / "app/tf_browser.py"
"""Path to the **TF Browser** page."""

PATH_PAGE_TF_VIEW = PATH_ROOT / "app/tf_view.py"
"""Path to the **TF Viewer** page."""

PATH_PAGE_PATTERN_EXPLORER = PATH_ROOT / "app/pattern_explorer.py"
"""Path to the **Pattern Explorer** page."""

PATH_PAGE_ABOUT = PATH_ROOT / "app/about.py"
"""Path to the **About** page."""

PATH_PAGE_HELP = PATH_ROOT / "app/about_help.py"
"""Path to the **Help** page."""

#endregion

#region Dataset
KAGGLE_HANDLE = "joejojoestar/tf-disco-datasets/versions/3"
"""The Kaggle handle for the datasets."""

class PathData:
    DISPROT: Path
    """Path to the DisProt dataset."""

    TFCLASSES: Path
    """Path to the TFClasses dataset."""

    TF_DBD_RANGES: Path
    """Path to the TSV file containing DBD regions for each TF."""

    TF_FASTA: Path
    """Path to the FASTA file containing all TF sequences."""

    TFS_DISORDER: Path
    """Path to the directory containing Disorder scores."""

    TF_PATTERNS: Path
    """Path to the TSV file containing a list of regex patterns."""

    TF_MATCHES: Path
    """Path to the TSV file with all matches of the patterns among all TFs."""

    def __init__(self, base_dir: Path):
        """Initialize paths to Data files, with the provided base directory."""
        self.set_base_dir(base_dir)

    def set_base_dir(self, base_dir: Path):
        """Re-initialize paths to Data files, with the provided base directory."""
        self.DISPROT = base_dir / "disprot_data.tsv"
        self.TFCLASSES = base_dir / "tfclasses.tsv"
        self.TF_DBD_RANGES = base_dir / "dbd_ranges.tsv"
        self.TF_FASTA = base_dir / "tfs_fasta/tfs_all.fasta"
        self.TFS_DISORDER = base_dir / "tfs_metrics/"
        self.TF_PATTERNS = base_dir / "patterns.tsv"
        self.TF_MATCHES = base_dir / "matches.tsv"

PATH_DATA = PathData(ENV_DATASET_PATH_OVERRIDE or Path("./data/")) # initial value

#endregion

#region Content
CONTENT_SUMMARY = f"""
{APP_NAME} (:primary[T]ranscription :primary[F]actor identity,
:primary[DIS]order, :primary[CO]ntext/motifs) is an integrated database and
visualization platform for transcription factors (TFs). It combines curated
information from [UniProt](https://www.uniprot.org),
[DisProt](https://disprot.org), and the
[Eukaryotic Linear Motif (ELM) Resource](http://elm.eu.org) with precomputed
disorder predictions from [AIUPred](https://aiupred.elte.hu),
[flDPnn](https://biomine.cs.vcu.edu/server-handler/?type=servers&target=flDPnn),
and [Metapredict v3](https://metapredict.net). {APP_NAME} lets users explore
known functional motifs in sequence context, examine their distribution within
DNA-binding and activation domains, identify shared motifs across TF sets, and
investigate relationships between motif occurrence and intrinsic disorder.
"""



CONTENT_CONTACT = f"""
For inquiries, feedback, or contributions, please reach out to us at
[debostuti@dubai.bits-pilani.ac.in](mailto:debostuti@dubai.bits-pilani.ac.in)
or make an [issue on GitHub](https://github.com/tf-disco/tf-disco/issues).
We welcome your input and look forward to hearing from you!
"""



CONTENT_HELP_TF_BROWSER = f"""
The user can use the :primary[:material/view_list: TF Browser] page to discover
and search amongst all the TFs available, which are :primary[listed in the
table].

To narrow down a search, :primary[filters] can be applied. There are *two* types
of search available:
1. :primary[**Basic search**]: Search for one specific TF, either by its
   **UniProt Accession**, **Genus number** or **Genus name**
2. :primary[**Advanced search**]: List all TFs which belong to a specific
   superclass/class/family/subfamily (refer to the [TFClass
   Resource](http://tfclass.bioinf.med.uni-goettingen.de/index.jsf) for details
   on TF classification)

To add a TF to the Cart, click on the :primary[:material/check_box:] icon.

Once TFs are selected, the Cart will be visible on the sidebar. One of the
following can be done:
- Select :primary[**1**] item in the Cart to :primary[**view its details**].
- Click on the :primary[:material/regular_expression: **Explore patterns**]
  button to explore the motifs occurring in all the TFs in the Cart.
"""



CONTENT_HELP_TF_VIEWER = f"""
This page shows the information of the selected TF. Basic info such as
UniProt accession number, DisProt ID (if exists), length of sequence, and the
ranges of positions marked as DNA Binding Domains (DBD) is shown.

The sequence is rendered 10-residue-long blocks. The sequence can be downloaded
in the FASTA format, if desired. There are two annotations which can be enabled
on the sequence:
1. DBD regions, denoted by a **Black/White** overline above the residues which
   are part of the DNA Binding Domains
2. Disorder predictors, denoted by a :red[Red] underline on residues which are
   predicted to be disordered. More on this below.

The Score Plots section displays a graph/plot of the metrics listed below. The
graph is interactive, and the user can make a selection with the mouse to zoom
in. Optionally, the scores can also be correlated with DisProt regions, by
showing bars under each graph to denote whether it agrees with DisProt or not.
The scores are available for download in tabular format (CSV/TSV), which
contains the per-residue scores for each of the metrics listed below.
- :primary[Hydropathy Index]: Kyte-Doolittle scale with window size 11
- :primary[AIUPred]: Disorder, Binding, Linker
- :primary[flDPnn]: Disorder, Protein binding, DNA binding, RNA binding, Linker
- :primary[Metapredict v3]: Disorder
- DisProt regions
- DBD ranges

The evidences of regions marked by DisProt can be viewed by the user in the
table below the graphs. The data from this table is obtained fully from the
DisProt database.

On the left sidebar, a table displays the ELM patterns that have matches
occurring in the currently selected TF. Click on the
:primary[:material/check_box:] checkbox next to a pattern to highlight its
matches in the sequence.
"""



def CONTENT_HELP_PATTERN_EXPLORER(is_on_help_page: bool=False):
    return f"""
The user can use the :primary[:material/regular_expression: Pattern Explorer]
page to explore the patterns/motifs occurring in a selection of TFs, or in the
entire dataset. The selection of TFs i.e. the cart can be made in the
:primary[:material/view_list: TF Browser] page first.

The table displays the ELM patterns which are available in the selected TFs. To
select a pattern for further analysis, click on the
:primary[:material/check_box:] checkbox next to the pattern.\\
These patterns are ranked with the help of the vagueness metric as described
{"[above](#vagueness)" if is_on_help_page else "in the [Help page](/help#vagueness)"}.

The sequence logo of the selected pattern is displayed below the table. Multiple
sequence alignment is performed using the [MUSCLE](https://drive5.com/muscle5)
tool. The sequence logo is generated using these
aligned sequences, which denotes the conservation of the pattern across the TFs.

The user can also visualize the sites in the TF's sequences where the pattern
occurs, where these sites are highlighted with a color. The user can choose to
download this data in the following formats:
1. HTML format, as shown on the page
2. Tabular format (CSV/TSV) where each row describes a matched site, along with
   whether the match is disordered or not, and whether it is part of the DNA Binding Domains (DBD).

Finally, the user can also view the frequency of the pattern's occurrences in
each TF, which is also available for download in tabular (CSV/TSV) format.
"""



CONTENT_DEVS = """
Website developed by:
- Joseph Cijo ([GitHub](https://github.com/joejo-joestar), [LinkedIn](https://linkedin.com/in/joseph-cijo), [:material/mail:](mailto:joecn2704+tfdisco@gmail.com))
- Sreenikethan Iyer ([GitHub](http://github.com/SreenikethanI), [LinkedIn](https://linkedin.com/in/sreenikethan-i), [:material/mail:](mailto:sreeni.s.iyer+tfdisco@gmail.com))
"""

#endregion
