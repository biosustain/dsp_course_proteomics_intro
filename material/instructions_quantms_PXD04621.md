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
