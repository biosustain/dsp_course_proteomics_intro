# QuantMS Hands-On for PXD04621

- [slides](slides/quantms_and_data_analysis.pdf)

## Open GitHub codespace

Use the following link to open a GitHub codespace with most of the required software installed:

> ⚠️
> If you do it manually, make sure to select the bigger machine with 4 cores and 16GB RAM

[![Open in Codespace deeplink](https://github.com/codespaces/badge.svg)](https://github.com/codespaces/new?hide_repo_select=true&ref=main&repo=949944579&skip_quickstart=true&machine=standardLinux32gb&devcontainer_path=.devcontainer%2Fdevcontainer.json&geo=EuropeWest)

## Download the data

### Option 1: Open in Codespace

> Will only work during the course. You do not to do anything.

We are using a setup script in `.devcontainer/setup.sh` to download the data from
Azure Blob Storage using [azcopy](https://learn.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10).
Access to the data is restricted during the time of the course only.

### Option 2: Git LFS

> Could be added to git LFS storage via GitHub, but the billing is uncertain.

### Option 3: OneDrive Options

> If you want to run it outside of the setting of this course

Currently the data is stored on a
[Google Drive folder](https://drive.google.com/drive/folders/1gxUh9nMx9icFLrI0vn3zAB9dDjZf-1Nh).
It should be downloaded when you run the tutorial in a GitHub codespace automatically.

In case you do not see the `mzML` files in the `data/PXD040621/mzML/` folder,
you can manually use [gdown](https://github.com/wkentaro/gdown) to download
these files from Google Drive:

```bash
conda activate base
pip install gdown
python 0_download_PXD040621_data.py
```

## Run the analysis

```bash
# export NXF_VER=25.10.4
nextflow run bigbio/quantms \
         -revision 1.7.0 \
         -params-file PXD040621_w_contaminants-params.yaml \
         -profile docker \
         -resume
```

If you run **locally on a Mac with Apple Silicion (M-ships)**, you need to addtionally the `arm` profile:

```bash
# export NXF_VER=25.10.4
nextflow run bigbio/quantms \
         -revision 1.7.0 \
         -params-file PXD040621_w_contaminants-params.yaml \
         -profile docker,arm \
         -resume
```

After 25 to 35mins you should see the following steps to be executed in the terminal:


<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<!-- This file was created with the aha Ansi HTML Adapter. https://github.com/theZiz/aha -->
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta http-equiv="Content-Type" content="application/xml+xhtml; charset=UTF-8"/>
<title>stdin</title>
</head>
<body>
<pre>

 N E X T F L O W   ~  version 25.10.4

Launching `https://github.com/bigbio/quantms` [sleepy_celsius] DSL2 - revision: e719f43a9f [1.7.0]

[97/ee0b14] BIGBIO_QUANTMS:QUANTMS:INPUT_CHECK:SAMPLESHEET_CHECK (PXD040621.sdrf.tsv)                                  [100%] 1 of 1 ✔
[fb/4dfc8c] BIGBIO_QUANTMS:QUANTMS:CREATE_INPUT_CHANNEL:SDRF_PARSING (PXD040621.sdrf.tsv)                              [100%] 1 of 1 ✔
[-        ] BIGBIO_QUANTMS:QUANTMS:FILE_PREPARATION:DECOMPRESS                                                         -
[15/199cc7] BIGBIO_QUANTMS:QUANTMS:FILE_PREPARATION:MZML_INDEXING (20220830_JL-4884_Forster_Ecoli_Suf_rep3_EG-7)       [100%] 8 of 8 ✔
[-        ] BIGBIO_QUANTMS:QUANTMS:FILE_PREPARATION:THERMORAWFILEPARSER                                                -
[2f/7c857b] BIGBIO_QUANTMS:QUANTMS:FILE_PREPARATION:MZML_STATISTICS (20220830_JL-4884_Forster_Ecoli_Suf_rep3_EG-7)     [100%] 8 of 8 ✔
[22/cdd7b8] BIGBIO_QUANTMS:QUANTMS:GENERATE_DECOY_DATABASE (1)                                                         [100%] 1 of 1 ✔
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:ID:PEPTIDE_DATABASE_SEARCH:COMET                                                -
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:ID:PSM_RESCORING:PERCOLATOR                                                     -
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:ID:PSM_FDR_CONTROL:ID_FILTER                                                    -
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:FEATURE_MAPPER:ISOBARIC_ANALYZER                                                -
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:FEATURE_MAPPER:ID_MAPPER                                                        -
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:FILE_MERGE                                                                      -
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:PROTEIN_INFERENCE:PROTEIN_INFERENCE_GENERIC                                     -
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:PROTEIN_INFERENCE:ID_FILTER                                                     -
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:PROTEIN_QUANT:ID_CONFLICT_RESOLVER                                              -
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:PROTEIN_QUANT:PROTEIN_QUANTIFIER                                                -
[-        ] BIGBIO_QUANTMS:QUANTMS:TMT:PROTEIN_QUANT:MSSTATS_CONVERTER                                                 -
[ee/ed8ddd] BIGBIO_QUANTMS:QUANTMS:LFQ:ID:PEPTIDE_DATABASE_SEARCH:COMET (20220830_JL-4884_Forster_Ecoli_Suf_rep2_EG-6) [100%] 8 of 8 ✔
[a8/33e8de] BIGBIO_QUANTMS:QUANTMS:LFQ:ID:PSM_RESCORING:PERCOLATOR (20220830_JL-4884_Forster_Ecoli_Suf_rep4_EG-8)      [100%] 8 of 8 ✔
[12/197e9e] BIGBIO_QUANTMS:QUANTMS:LFQ:ID:PSM_FDR_CONTROL:ID_FILTER (20220830_JL-4884_Forster_Ecoli_Suf_rep2_EG-6)     [100%] 8 of 8 ✔
[71/513331] BIGBIO_QUANTMS:QUANTMS:LFQ:PROTEOMICSLFQ (PXD040621.sdrf_openms_design)                                    [100%] 1 of 1 ✔
[3b/44e4b4] BIGBIO_QUANTMS:QUANTMS:SUMMARY_PIPELINE                                                                    [100%] 1 of 1 ✔
Plus 7 more processes waiting for tasks…
-[bigbio/quantms] Pipeline completed successfully-


Completed at: 01-Jun-2026 12:47:20
Duration    : 1m 56s
CPU hours   : 1.9 (93.9% cached)
Succeeded   : 1
Cached      : 44


</pre>
</body>
</html>


## Quality control of the analysis

You can inspect the generated QC report (using pMulitQC which is an extension of MultiQC)
in the `results/PXD040621/pmultiqc/` folder.

Download a pre-created report and open it in your browser:

- [pMultiQC report for quantms analysis of PXD04621 (v1.7.0)](./pmultiqc_report/multiqc_report_pmultiqc.html)

## Copy files for further analysis

```bash
cp -aL results/PXD040621/proteomicslfq/. data/PXD040621/processed/
# And maybe save the parameters to reproduce the analysis (for 1.4.0 and above):
cp -aL results/PXD040621/pipeline_info/. data/PXD040621/processed/pipeline_info/
```

Now it is safe to delete the `results/PXD040621` and `work` folder.

```bash
rm -r results/PXD040621 work
```

## Clean up unused docker images

As we are running in a GitHub codespace, we have limited storage. Therefore let's
[clean up our docker images](https://docs.docker.com/engine/manage-resources/pruning/)
store after the analysis is done.

```bash
docker images # see all images
docker image prune -a
```

### Free up some more space

> Should not be necessary

```bash
# some cache files
rm -r  /.codespaces/bin/cache/bin/linux-x64/
ls /vscode/extensionsCache/
ls /vscode/serverCache/
```

### Check used disk space

Check the local storage usage (you have maximum of 32GB in a GitHub codespace)
in the root folder with the following command:

```bash
du -hd 1 /
```

## Run a different analysis

- feel free to run the project `PXD041301` in quantms if you want an exercise

> ⚠️
> not tested in codespace and probably too big for the codespace storage
