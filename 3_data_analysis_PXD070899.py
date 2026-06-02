# ---
# jupyter:
#   jupytext:
#     cell_metadata_filter: tags,-all
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: base
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Exercise: Data Analysis PXD070899
#
# Plan
# - read data and log2 transform intensity values
# - aggregate peptide intensities to protein intensities
# - format data from long to wide format
# - remove contaminant proteins
# - check for missing values
# - Clustermap of sample and proteins
# - differential analysis (Volcano Plots)
# - Enrichment Analysis
# - check for maltose update pathway (Fig. 3 in paper)

# %% tags=["hide-output"]
# %pip install acore vuecore "pingouin<0.6.0"

# %%
from pathlib import Path

import acore.differential_regulation
import acore.enrichment_analysis
import acore.normalization
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import plotly.express as px
import scipy.stats
import seaborn as sns
import vuecore
from acore.io.uniprot import fetch_annotations, process_annotations
from vuecore.viz import get_enrichment_plots

# %% [markdown]
# # Paramters
# - `file_in`: input file with the quantified peptide data in MSstats format
#    as provided by quantms
# - `out_dir`: output directory for the results of the data analysis, 
#    which will be used later for the report generation with VueGen.
#
# The file will be loaded from the online repository if it is not present.

# %% tags=["parameters"]
file_in: str = Path(
    "data/PXD070899/processed/report.pg_matrix.tsv"
) # input file with the quantified peptide data in MSstats format as provided by quantms
out_dir = "data/PXD070899/report/" # output directory for the results of the data analysis, which will be used later for the report generation with VueGen.
min_obs_per_group:int = 3 # minimum number of observations per group for a protein to be included in the differential regulation analysis

# %% [markdown]
# Create output directory if it does not exist

# %%
out_dir = Path(out_dir)
out_dir.mkdir(parents=True, exist_ok=True)
print(f"Output directory: {out_dir}")

# %% [markdown]
# We have the following columns in the data:

# %%
if not file_in.exists():
    file_in = (
        "https://raw.githubusercontent.com/biosustain/dsp_course_proteomics_intro/HEAD"
        "/data/PXD070899/processed/report.pg_matrix.tsv"
    )
df = pd.read_csv(file_in, sep="\t", header=0, index_col=0)  # .set_index([])
df.head()

# %% [markdown]
# Potentiall clean filepath names in DIANN output

# %%
df.columns = df.columns.str.split(r"/|\\").str[-1]
df.head()

# %% [markdown]
# The first 6 columns contain the meta information about the peptides, 
# while the remaining columns contain the intensities.

# %%
proteins_meta = df.iloc[:,:5]
proteins_meta 

# %%
proteins = df.iloc[:,5:].T
proteins.index.name = "SampleID"
proteins.columns.name = "ProteinName"
proteins.head()

# %% [markdown]
# Log2 transform the intensity values and remove contaminant proteins

# %%
to_drop = proteins.filter(regex="cRAP-|CON_", axis=1).columns
to_drop

# %%
proteins = np.log2(proteins).drop(to_drop, axis=1)
proteins

# %% [markdown]
# Add label encoding


# %%
# label_encoding = {"WT": 0, "AYa2022": 1, "AYa18": 2, "AY": 3}
label_suf = pd.Series(
    proteins.index.str.split("_").str[-2],
    index=proteins.index,
    name="condition",
)
label_suf
# %% [markdown]
# # Homework
# Repeat the analysis based on the tutorial from the course.
#
# > Small adjustments are needed.
#
# - copy bit by bit
