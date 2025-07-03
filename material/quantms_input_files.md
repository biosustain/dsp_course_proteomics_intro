# QuantMS Input Files

See the documentation for an overview: [docs.quantms.org](https://docs.quantms.org/en/latest/formats.html)

## Sample Data Relationship Format (SDRF) file

- describes the relationship between the biological samples and the (raw) data files from 
  the mass spectrometer

## Fasta File

- defines the search space: only these protein sequences will be considered for the spectrum database search

## Spectrum Files

- contains the spectra which are compared to the protein sequences in the fasta file to 
  identify the proteins (actually peptides) in the sample

## QuantMS Parameters file

- defines the parameters for the analysis, e.g. the search engine to use, the quantification method, 
  the database to use, etc.
