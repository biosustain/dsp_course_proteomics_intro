# Massspectrometry-based Proteomcis introduction

The course will introduce MS-based proteomics data, it's basic processing and
the downstream analysis. The course will be a mix of lectures and hands-on session
where we will perform the basic processing steps using
[quantms](https://docs.quantms.org/en/latest/), a nextflow workflow, and
[acore, short for analytical core](https://analytics-core.readthedocs.io/latest/)
for the downstream analysis.

## Agenda

| Time          | Topic                                                 | lecturer        |
| ------------- | ----------------------------------------------------- | --------------- |
| 8.30 - 10.00  | Introduction with overview of all the components      | Marco Reverenna |
| 10.30 - 12.00 | Steps in data processing and running quantms hands-on | Henry Webel     |
| 12.00 - 13.00 | Lunch (sandwiches are provided)                       | -               |
| 13.00 - 14.30 | Steps in statistical analysis (lecture )              | Alberto Santos  |
| 15.00 - 16.30 | Steps in statistical analysis (Hands-On)              | Henry Webel     |

### Some details to the agenda points

#### 10.30-12.00 - details - data processing and hands-on quantms

Steps in data processing (using [quantms](https://docs.quantms.org/en/latest/))

- FASTA file to define search space
- Spectrum files from Mass-spectrometer
- Running quantms to process spectra to identified and quantified peptide sequences

#### 15.00 - 16.30 - details - hands-on statistical analysis

Basic statistical analysis of a two-group experiment with one timepoint (option 1) or four timepoints (option 2)

- Peptide to protein (group) aggregation
- Downstream data analysis of proteins (using [analytical core library](https://analytics-core.readthedocs.io/latest/)
  developed at biosustain and other Python libraries)
- Building a report with [vuegen](https://vuegen.readthedocs.io/en/latest/) reports (developed at biosustain)

## Run nextflow

```bash
nextflow run bigbio/quantms \
         -revision 1.3.0 \
         -params-file PXD040621_w_contaminants-params.yaml \
         -profile docker \
         -resume
```

## help commands

Check the local storage usage (you have maximum of 32GB in a GitHub codespace)

```bash
du -hd 1
```

See the downloaded docker images

```bash
docker images
```

### Free up some more space

```bash
# some cache files
rm -r  /.codespaces/bin/cache/bin/linux-x64/
ls /vscode/extensionsCache/
ls /vscode/serverCache/
```

## Tools

- [gdown](https://github.com/wkentaro/gdown) to download files from Google Drive

```
pip install gdown
```
