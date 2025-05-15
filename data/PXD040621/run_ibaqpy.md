# IbaqPy

```bash
pip install git+https://github.com/bigbio/quantms.io
pip install ibaqpy
```

Then try to follow the steps as outlined in the 
[IbaqPy documentation](https://github.com/bigbio/ibaqpy/tree/master?tab=readme-aov-file#from-quantms-to-ibaq-values)

## Convert to quantms.io format

```bash
quantmsioc convert-feature --sdrf_file PXD040621.sdrf.tsv --msstats_file PXD040621.sdrf_openms_design_msstats_in.csv --mztab_file PXD040621.sdrf_openms_design_openms.mzTab --file_num 30 --output_folder res --duckdb_max_memory 8GB --output_prefix_file PXD040621
quantmsioc convert-ibaq --feature_file res/PXD040621-6c224f5a-7c1f-46f9-9dae-1541baeef8fe.feature.parquet --sdrf_file PXD040621.sdrf.tsv --output_folder ibaq --output_prefix_file PXD040621
```
