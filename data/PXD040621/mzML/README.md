# mzML files

mzML files are storing the mass spectra data in a text based format simialar to xml. 
The files are generated using the ThermorawFileParser, which take thermo raw files as input 
and convert them to mzML format. If you provide raw files to quantms, the raw files will 
be processed and stored as mzML files in the `thermorawFileReader` folder in the results
directory. The mzML files are then used for identification and quantification of the 
peptides recorded in the spectra.

Files expected in this folder are:

```bash
LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_01.mzML
LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_03.mzML
LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_02.mzML
LFQ_Orbitrap_DDA_Condition_A_Sample_Alpha_02.mzML
LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_01.mzML
LFQ_Orbitrap_DDA_Condition_B_Sample_Alpha_03.mzML
```

