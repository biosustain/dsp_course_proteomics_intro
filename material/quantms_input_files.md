# QuantMS Input Files

See the documentation for an overview: [docs.quantms.org](https://docs.quantms.org/en/latest/formats.html)

## Sample Data Relationship Format (SDRF) file

- describes the relationship between the biological samples and the (raw) data files from
  the mass spectrometer

Example from PXD04621: [PXD040621.sdrf.tsv](../data/PXD040621/PXD040621.sdrf.tsv)

## Fasta File

- defines the search space: only these protein sequences will be considered for the
  spectrum database search

Example from PXD04621: [merged_ecoli_with_contaminants.fasta](../data/fasta/merged_ecoli_with_contaminants.fasta)

## Spectrum Files

- contains the spectra which are compared to the protein sequences in the fasta file to
  identify the proteins (actually peptides) in the sample

Example from PXD04621, see download instructions on the page
['QuantMS Hands-On for PXD04621'](instructions_quantms_PXD04621.md)

## QuantMS Parameters file

- defines the parameters for the analysis, e.g. the search engine to use, the quantification method,
  the database to use, etc.

Example from PXD04621: [PXD040621_w_contaminants-params.yaml](../PXD040621_w_contaminants-params.yaml)
