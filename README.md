# Mass spectrometry-based Proteomics introduction

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

Find the instruction [here](docs/instructions_quantms_PXD04621.md)

Steps in data processing (using [quantms](https://docs.quantms.org/en/latest/))

- FASTA file to define search space
- Spectrum files from Mass-spectrometer
- Running quantms to process spectra to identified and quantified peptide sequences

#### 15.00 - 16.30 - details - hands-on statistical analysis

Find the instruction [here](docs/instructions_statistics_PXD04621.md)

Basic statistical analysis of a two-group experiment with one timepoint (option 1) or four timepoints (option 2)

- Peptide to protein (group) aggregation
- Downstream data analysis of proteins (using [analytical core library](https://analytics-core.readthedocs.io/latest/)
  developed at biosustain and other Python libraries)
- Building a report with [vuegen](https://vuegen.readthedocs.io/en/latest/) reports (developed at biosustain)

## QuantMS help

- see the documentation for an overview: [docs.quantms.org](https://docs.quantms.org)
- ask question on the nf-core slack channel `quantms`: [https://nf-co.re/join/slack](https://nf-co.re/join/slack)
- submit an issue on the [GitHub repository](https://github.com/bigbio/quantms/issues)






