# %% [markdown]
# # Protein Quantification for PXD040621
# - Using DirectLFQ: https://github.com/MannLabs/directlfq
# Python 3.9 environment with directlfq
# - `pip install 'directlfq[stable]'`

# %%
from pathlib import Path

import directlfq.lfq_manager as lfq_manager
import numpy as np
import pandas as pd

fname_peptides = Path("../processed/PXD040621.sdrf_openms_design_msstats_in.csv")


folder_out = Path(".")
folder_out.mkdir(exist_ok=True)

fname_out_peptides_directlfq = (
    folder_out / "PXD040621_peptides_directlfq_in.aq_reformat.tsv"
)


df_peptides = pd.read_csv(fname_peptides)
df_peptides.head()

# %%

df_peptides = df_peptides.rename(columns={"ProteinName": "protein"}).assign(
    ion=df_peptides["PeptideSequence"] + df_peptides["PrecursorCharge"].astype(str),
    run=df_peptides["Condition"] + "_run_" + df_peptides["Run"].astype("str"),
)
df_peptides[["protein", "ion", "Intensity", "run", "Reference", "Condition"]]

# %%
df_peptides_wide = df_peptides.pivot(
    index=["protein", "ion"],
    columns="run",
    values="Intensity",
)
df_peptides_wide.to_csv(fname_out_peptides_directlfq, sep="\t")
df_peptides_wide

# %%
lfq_manager.run_lfq(
    input_file=fname_out_peptides_directlfq,
    min_nonan=2,
    maximum_number_of_quadratic_ions_to_use_per_protein=10,  # default 10 if None
    number_of_quadratic_samples=50,  # default 50 if None
    num_cores=2,  # default all cores -1 if None
)

# %%
fname_in_proteins = fname_out_peptides_directlfq.with_suffix(
    ".tsv.protein_intensities.tsv"
)
df_proteins = pd.read_csv(fname_in_proteins, sep="\t", index_col=["protein"])
df_proteins.head()

# %%
df_proteins = df_proteins.loc[~df_proteins.index.str.contains("CON_")]
df_proteins = np.log2(df_proteins.replace(0, np.nan)).T
df_proteins

# %%
df_proteins.notna().sum(axis=1).describe()

# %%
df_proteins.insert(0, "Condition", df_proteins.index.str.split("_run_").str[0])

# %%
df_proteins.to_csv(
    folder_out / "PXD040621_proteins_directlfq_log2_intensities.tsv", sep="\t"
)

# %%
