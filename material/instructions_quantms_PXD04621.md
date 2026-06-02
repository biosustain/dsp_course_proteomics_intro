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


### Setup parameters yourself

you can use a tool called nf-core to interactively select parameters

```bash
nf-core pipelines launch bigbio/quantms -r 1.7.0
```

or list all parameter with

```bash
# this is enabled in nextflow.config
nextflow run bigbio/quantms -r 1.7.0 --help
```

<details>
<summary>Help message for bigbio/quantms v1.7.0</summary>

```bash
>>> nextflow run bigbio/quantms -r 1.7.0 --help
N E X T F L O W   ~  version 25.10.4

Launching `https://github.com/bigbio/quantms` [chaotic_agnesi] DSL2 - revision: e719f43a9f [1.7.0]

--help                                    [boolean, string] Show the help message for all top level parameters. When a parameter is given to `--help`, the full help message of that parameter will be printed.
--helpFull                                [boolean]         Show the help message for all non-hidden parameters.
--showHidden                              [boolean]         Show all hidden parameters in the help message. This needs to be used in combination with `--help` or `--helpFull`.

Input/output options
  --input                                 [string]  URI/path to an [SDRF](https://github.com/bigbio/proteomics-metadata-standard/tree/master/annotated-projects) file (.sdrf.tsv) **OR** [OpenMS-style experimental
design](https://abibuilder.cs.uni-tuebingen.de/archive/openms/Documentation/release/latest/html/classOpenMS_1_1ExperimentalDesign.html#details) with paths to spectra files (.tsv)
  --outdir                                [string]  The output directory where the results will be saved. [default: ./results]
  --email                                 [string]  Email address for completion summary.
  --multiqc_title                         [string]  MultiQC report title. Printed as page header, used for filename if not otherwise specified.
  --root_folder                           [string]  Root folder in which the spectrum files specified in the SDRF/design are searched
  --local_input_type                      [string]  Overwrite the file type/extension of the filename as specified in the SDRF/design [default: mzML]
  --acquisition_method                    [string]  Proteomics data acquisition method  (accepted: dda, dia)
  --id_only                               [boolean] Only perform identification subworkflow.
  --export_decoy_psm                      [boolean] Whether export PSM from decoy in final identification results [default: true]

SDRF validation
  --validate_ontologies                   [boolean] Check that ontology terms in an input SDRF file exist. [default: true]
  --skip_ms_validation                    [boolean] Skip validation of mass spectrometry files.
  --skip_factor_validation                [boolean] Skip validation of factor columns. [default: true]
  --skip_experimental_design_validation   [boolean] Skip validation of experimental design.
  --use_ols_cache_only                    [boolean] Use cached version of the Ontology Lookup Service (OLS). [default: true]

Protein database
  --database                              [string]  The `fasta` protein database used during database search. *Note:* For DIA data, it must not contain decoys.
  --add_decoys                            [boolean] Generate and append decoys to the given protein database
  --decoy_string                          [string]  Pre- or suffix of decoy proteins in their accession [default: DECOY_]
  --decoy_string_position                 [string]  Location of the decoy marker string in the `fasta` accession. Before (prefix) or after (suffix) [default: prefix]
  --decoy_method                          [string]  Choose the method to produce decoys from input target database.  (accepted: reverse, shuffle) [default: reverse]
  --shuffle_max_attempts                  [integer] Maximum nr. of attempts to lower the amino acid sequence identity between target and decoy for the shuffle algorithm [default: 30]
  --shuffle_sequence_identity_threshold   [number]  Target-decoy amino acid sequence identity threshold for the shuffle algorithm. if the sequence identity is above this threshold, shuffling is repeated. In case of repeated failure, individual amino acids are 'mutated' to produce a difference amino acid
sequence. [default: 0.5]

Spectrum preprocessing
  --openms_peakpicking                    [boolean] Activate OpenMS-internal peak picking
  --peakpicking_inmemory                  [boolean] Perform peakpicking in memory
  --peakpicking_ms_levels                 [string]  Which MS levels to pick as comma separated list. Leave empty for auto-detection.
  --convert_dotd                          [boolean] Convert bruker .d files to mzML
  --reindex_mzml                          [boolean] Force initial re-indexing of input mzML files. Also fixes some common mistakes in slightly incomplete/outdated mzMLs. (Default: true for safety) [default: true]
  --mzml_features                         [boolean] Compute with mzmlstatistics step the features at MS1 level and output to a RAW file, only available for mzML files

Database search
  --search_engines                        [string]  A comma separated list of search engines to use (and combine). Valid: comet, msgf, sage [default: comet]
  --sage_processes                        [integer] Number of sage processes to be spawned. [default: 1]
  --enzyme                                [string]  The enzyme to be used for in-silico digestion, in 'OpenMS format' [default: Trypsin]
  --met_excision                          [boolean] Database searches accounted for N-terminal methionine excision, a common co-translational modification where the initial methionine is enzymatically removed from proteins. [default: true]
  --num_enzyme_termini                    [string]  Specify the amount of termini matching the enzyme cutting rules for a peptide to be considered. Valid values are `fully` (default), `semi`, or `none`  (accepted: fully, semi, none) [default: fully]
  --allowed_missed_cleavages              [integer] Specify the maximum number of allowed missed enzyme cleavages in a peptide. The parameter is not applied if `unspecific cleavage` is specified as enzyme. [default: 2]
  --precursor_mass_tolerance              [integer] Precursor mass tolerance used for database search. For High-Resolution instruments a precursor mass tolerance value of 5 ppm is recommended (i.e. 5). See also [`--precursor_mass_tolerance_unit`](#precursor_mass_tolerance_unit). [default: 5]
  --precursor_mass_tolerance_unit         [string]  Precursor mass tolerance unit used for database search. Possible values are 'ppm' (default) and 'Da'.  (accepted: Da, ppm) [default: ppm]
  --fragment_mass_tolerance               [number]  Fragment mass tolerance used for database search. The default of 0.03 Da is for high-resolution instruments. [default: 0.03]
  --fragment_mass_tolerance_unit          [string]  Fragment mass tolerance unit used for database search. Possible values are 'ppm' (default) and 'Da'.  (accepted: Da, ppm) [default: Da]
  --fixed_mods                            [string]  A comma-separated list of fixed modifications with their Unimod name to be searched during database search [default: Carbamidomethyl (C)]
  --variable_mods                         [string]  A comma-separated list of variable modifications with their Unimod name to be searched during database search [default: Oxidation (M)]
  --isotope_error_range                   [string]  Comma-separated range of integers with allowed isotope peak errors for precursor tolerance (like MS-GF+ parameter '-ti'). E.g. -1,3 [default: 0,1]
  --instrument                            [string]  Type of instrument that generated the data. 'low_res' or 'high_res' (default; refers to LCQ and LTQ instruments)
  --protocol                              [string]  MSGF only: Labeling or enrichment protocol used, if any (options: 'automatic', 'phospho', 'iTRAQ', 'iTRAQ_phospho', 'TMT', 'none') Default: automatic [default: automatic]
  --min_precursor_charge                  [integer] Minimum precursor ion charge. Omit the '+' [default: 2]
  --max_precursor_charge                  [integer] Maximum precursor ion charge. Omit the '+' [default: 4]
  --min_peptide_length                    [integer] Minimum peptide length to consider (works with MSGF and in newer Comet versions) [default: 6]
  --max_peptide_length                    [integer] Maximum peptide length to consider (works with MSGF and in newer Comet versions) [default: 40]
  --num_hits                              [integer] Specify the maximum number of top peptide candidates per spectrum to be reported by the search engine. Default: 1 [default: 1]
  --max_mods                              [integer] Maximum number of modifications per peptide. If this value is large, the search may take very long. [default: 3]
  --min_peaks                             [integer] Minimum number of peaks in the spectrum to be considered for the search engine. Default: 10 [default: 10]
  --min_pr_mz                             [number]  The minimum precursor m/z for the in silico library generation or library-free search [default: 400.0]
  --max_pr_mz                             [number]  The maximum precursor m/z for the in silico library generation or library-free search [default: 2400.0]
  --min_fr_mz                             [number]  The minimum fragment m/z for the in silico library generation or library-free search [default: 100.0]
  --max_fr_mz                             [number]  The maximum fragment m/z for the in silico library generation or library-free search [default: 1800.0]

Modification localization
  --enable_mod_localization               [boolean] Turn the mechanism on.
  --mod_localization                      [string]  Which variable modifications to use for scoring their localization. [default: Phospho (S),Phospho (T),Phospho (Y)]

Peptide re-indexing
  --unmatched_action                      [string]  What to do when peptides are found that do not follow a unified set of rules (since search engines sometimes differ in their interpretation of them).   (accepted: warn, error, remove) [default: warn]
  --IL_equivalent                         [boolean] Should isoleucine and leucine be treated interchangeably when mapping search engine hits to the database? Default: true [default: true]

PSM re-scoring (general)
  --skip_rescoring                        [boolean] Skip PSM rescoring steps for specific cases, such as studying pure search engine results and search engine ranks
  --ms2features_enable                    [boolean] Whether to enable MS2-based features for Percolator using deep learning predictors such as ms2pip, alphapeptdeep, and DeepLC.
  --ms2features_snr                       [boolean] Whether to add signal-to-noise ratio features for identification rescoring in percolator
  --ms2features_fine_tuning               [boolean] Whether to fine tune AlphaPeptDeep for identification rescoring
  --fine_tuning_sample_run                [integer] The number of sample ms run [default: 1]
  --ms2features_range                     [string]  MS2Features processing range: independent run, Sample or whole experiments  (accepted: independent_run, by_sample, by_project) [default: independent_run]
  --ms2features_model                     [string]  Deep learning model for MS2 prediction. Choose based on your feature generator: alphapeptdeep uses 'generic' model, ms2pip uses instrument-specific models.  (accepted: generic, HCD2021, HCD2019, CID, iTRAQ, iTRAQphospho, TMT, TTOF5600, HCDch2,
CIDch2, Immuno-HCD, CID-TMT, timsTOF2023, timsTOF2024) [default: generic]
  --ms2features_model_dir                 [string]  The path of ms2 prediction model files. Providing model file to avoid repeated download and slow internet connection
  --ms2features_generators                [string]  Feature generators for MS2-based rescoring. Comma-separated list of tools to use. [default: deeplc,alphapeptdeep]
  --ms2features_calibration               [number]  Percentage of PSMs used for DeepLC calibration set (retention time prediction) [default: 0.15]
  --ms2features_tolerance                 [number]  Fragment mass tolerance for MS2Features (fallback when SDRF values are incompatible) [default: 0.05]
  --ms2features_tolerance_unit            [string]  Unit for ms2features_tolerance (Da or ppm)  (accepted: Da, ppm) [default: Da]
  --ms2features_force                     [boolean] Force use of specified MS2 model without validation or best model selection
  --ms2features_modloss                   [boolean] Enable modification loss ion features for MS2 prediction
  --ms2features_best                      [boolean] Automatically find and use the best available MS2 model for your data [default: true]
  --force_transfer_learning               [boolean] Force save fine-tuning model [default: false]
  --epoch_to_train_ms2                    [integer] The number of fine-tuning epoch [default: 20]
  --transfer_learning_test_ratio          [number]  The proportion of test data used for comparing fine-tuned models with pre-trained models [default: 0.3]
  --ms2features_debug                     [boolean] Enable debug logging for quantms-rescoring tool.
  --run_fdr_cutoff                        [number]  FDR cutoff on PSM level (or peptide level; see Percolator options) *per run* before going into feature finding, map alignment and inference. This can be seen as a pre-filter. See  [default: 0.1]

PSM re-scoring (Percolator)
  --fdr_level                             [string]  Calculate FDR on PSM ('psm_level_fdrs') or peptide level ('peptide_level_fdrs')?  (accepted: peptide_level_fdrs, psm_level_fdrs) [default: psm_level_fdrs]
  --train_FDR                             [number]  The FDR cutoff to be used during training of the SVM. [default: 0.05]
  --test_FDR                              [number]  The FDR cutoff to be used during testing of the SVM. [default: 0.05]
  --subset_max_train                      [integer] Only train an SVM on a subset of PSMs, and use the resulting score vector to evaluate the other PSMs. Recommended when analyzing huge numbers (>1 million) of PSMs. When set to 0, all PSMs are used for training as normal. This is a runtime vs. quality
tradeoff. Default: 300,000 [default: 300000]
  --description_correct_features          [integer] Use additional features whose values are learnt by correct entries. See help text. Default: 0 = none

Consensus ID
  --consensusid_algorithm                 [string]  How to combine the probabilities from the single search engines: best, combine using a sequence similarity-matrix (PEPMatrix), combine using shared ion count of peptides (PEPIons). See help for further info.  (accepted: best, PEPMatrix, PEPIons)
[default: best]
  --consensusid_considered_top_hits       [integer] Only use the top N hits per search engine and spectrum for combination. Default: 0 = all
  --min_consensus_support                 [number]  A threshold for the ratio of occurrence/similarity scores of a peptide in other runs, to be reported. See help.

Isobaric analyzer
  --quant_activation_method               [string]  Operate only on MSn scans where any of its precursors features a certain activation method. Set to empty to disable.  (accepted: HCD, CID, ETD, ECD) [default: HCD]
  --reporter_mass_shift                   [number]  Allowed shift (left to right) in Th from the expected position [default: 0.002]
  --min_precursor_intensity               [number]  Minimum intensity of the precursor to be extracted [default: 1.0]
  --min_precursor_purity                  [number]  Minimum fraction of the total intensity. 0.0:1.0 [default: 0.0]
  --min_reporter_intensity                [number]  Minimum intensity of the individual reporter ions to be extracted. [default: 0.0]
  --precursor_isotope_deviation           [number]  Maximum allowed deviation (in ppm) between theoretical and observed isotopic peaks of the precursor peak [default: 10.0]
  --isotope_correction                    [boolean] Enable isotope correction (highly recommended)
  --plex_corr_matrix_file                 [string]  Path to the correction matrix file for isobaric labelling, defaults are in assets folder
  --iso_normalization                     [boolean] Enable normalization of the channel intensities
  --reference_channel                     [string]  The reference channel, e.g. for calculating ratios. [default: 126]

Protein inference
  --protein_inference_method              [string]  The inference method to use. 'aggregation' (default) or 'bayesian'.  (accepted: aggregation, bayesian) [default: aggregation]
  --protein_score                         [string]  [Ignored in Bayesian] How to aggregate scores of peptides matching to the same protein  (accepted: best, product, sum) [default: best]
  --use_shared_peptides                   [boolean] [Ignored in Bayesian] Also use shared peptides during score aggregation to protein level [default: true]
  --min_peptides_per_protein              [integer] [Ignored in Bayesian] Minimum number of peptides needed for a protein identification [default: 1]
  --top_PSMs                              [integer] Consider only the top X PSMs per spectrum to find the best PSM per peptide. 0 considers all. [default: 1]
  --protein_level_fdr_cutoff              [number]  The experiment-wide protein (group)-level FDR cutoff. Default: 0.01 [default: 0.01]
  --picked_fdr                            [boolean] Use picked protein FDRs [default: true]
  --psm_level_fdr_cutoff                  [number]  The experiment-wide PSM-level FDR cutoff. Default: 0.01 [default: 0.01]

Protein Quantification (DDA)
  --labelling_type                        [string]  Specify the labelling method that was used. Will be ignored if SDRF was given but is mandatory otherwise  (accepted: label free sample, itraq4plex, itraq8plex, tmt6plex, tmt10plex, tmt11plex, tmt16plex)
  --top                                   [integer] Calculate protein abundance from this number of proteotypic peptides (most abundant first; '0' for all, Default 3) [default: 3]
  --average                               [string]  Averaging method used to compute protein abundances from peptide abundances.  (accepted: median, mean, weighted_mean, sum) [default: median]
  --best_charge_and_fraction              [boolean] Distinguish between fraction and charge states of a peptide. (default: 'false')
  --ratios                                [boolean] Add the log2 ratios of the abundance values to the output.
  --normalize                             [boolean] Scale peptide abundances so that medians of all samples are equal.(Default false)
  --fix_peptides                          [boolean] Use the same peptides for protein quantification across all samples.(Default false)
  --include_all                           [boolean] Include results for proteins with fewer proteotypic peptide than indicated by top. [default: true]
  --protein_quant                         [string]  Quantify proteins based on:

* 'unique_peptides' = use peptides mapping to single proteins or a group of indistinguishable proteins (according to the set of experimentally identified peptides)
* 'strictly_unique_peptides' (only LFQ) = use peptides
mapping to a unique single protein only
* 'shared_peptides' = use shared peptides, too, but only greedily for its best group (by inference score and nr. of peptides)  (accepted: unique_peptides, strictly_unique_peptides, shared_peptides) [default:
unique_peptides]
  --export_mztab                          [boolean] Export the results in mzTab format. [default: true]

Protein Quantification (LFQ)
  --quantification_method                 [string]  Choose between feature-based quantification based on integrated MS1 signals ('feature_intensity'; default) or spectral counting of PSMs ('spectral_counting'). **WARNING:** 'spectral_counting' is not compatible with our MSstats step yet. MSstats will
therefore be disabled automatically with that choice.  (accepted: feature_intensity, spectral_counting) [default: feature_intensity]
  --mass_recalibration                    [boolean] Recalibrates masses based on precursor mass deviations to correct for instrument biases. (default: 'false')
  --targeted_only                         [boolean] Only looks for quantifiable features at locations with an identified spectrum. Set to false to include unidentified features so they can be linked and matched to identified ones (= match between runs). (default: 'true') [default: true]
  --feature_with_id_min_score             [number]  The minimum probability (e.g.: 0.25) an identified (=id targeted) feature must have to be kept for alignment and linking (0=no filter). [default: 0.1]
  --feature_without_id_min_score          [number]  The minimum probability (e.g.: 0.75) an unidentified feature must have to be kept for alignment and linking (0=no filter). [default: 0.75]
  --lfq_intensity_threshold               [number]  The minimum intensity for a feature to be considered for quantification. (default: '1000') [default: 1000.0]
  --alignment_order                       [string]  The order in which maps are aligned. Star = all vs. the reference with most IDs (default). TreeGuided = an alignment tree is calculated first based on similarity measures of the IDs in the maps.  (accepted: star, treeguided) [default: star]
  --quantify_decoys                       [boolean] Also quantify decoys? (Usually only needed for Triqler post-processing output with [`--add_triqler_output`](#add_triqler_output), where it is auto-enabled)

DIA-NN
  --mass_acc_automatic                    [boolean] Choosing the MS2 mass accuracy setting automatically [default: true]
  --scan_window_automatic                 [boolean] Choosing scan_window setting automatically [default: true]
  --scan_window                           [integer] Set the scan window radius to a specific value [default: 8]
  --performance_mode                      [boolean] Set Low RAM & High Speed Mode for DIANN, including min-corr, corr-diff, and time-corr-only three parameters [default: true]
  --quick_mass_acc                        [boolean] when choosing the MS2 mass accuracy setting automatically, DIA-NN will use a fast heuristical algorithm instead of IDs number optimisation [default: true]
  --pg_level                              [number]  Controls the protein inference mode  (accepted: 0, 1, 2) [default: 2]
  --species_genes                         [boolean] Instructs DIA-NN to add the organism identifier to the gene names
  --diann_speclib                         [string]  The spectral library to use for DIA-NN
  --diann_report_decoys                   [boolean] Save decoy PSMs to the main .parquet report for DIA-NN 2.0.*
  --diann_export_xic                      [boolean] instructs DIA-NN to extract MS1/fragment chromatograms for identified precursors within X seconds from the elution apex, with X set to 10s if not provided;equivalent to the 'XICs' option in the GUI
  --skip_preliminary_analysis             [boolean] Skip the preliminary analysis step, thus use the passed spectral library as-is insted of generating a local concensus library.
  --empirical_assembly_log                [string]  The log file for the empirical assembly, Only used if `--skip_preliminary_analysis` is set to `true` and `--diann_speclib` is passed. If passed, will use that log file to carry out the DIA-NN search, instead of running a preliminary search.
  --diann_normalize                       [boolean] Enable cross-run normalization between runs by diann. [default: true]
  --random_preanalysis                    [boolean] Enable random selection of spectrum files to generate empirical library.
  --random_preanalysis_seed               [integer] Set the random seed for the random selection of spectrum files to generate the empirical library. [default: 42]
  --enable_diann_mztab                    [boolean] Export the DIA-NN and DIA results to mzTab [default: false]

Statistical post-processing
  --skip_post_msstats                     [boolean] Skip MSstats/MSstatsTMT for statistical post-processing?
  --ref_condition                         [string]  Experimental: Instead of all pairwise contrasts (default), uses the given condition name/number (corresponding to your experimental design) as a reference and creates pairwise contrasts against it.
  --contrasts                             [string]  Experimental: Allows full control over contrasts by specifying a set of contrasts in a semicolon separated list of R-compatible contrasts with the condition names/numbers as variables (e.g. `1-2;1-3;2-3`). Overwrites
[`--ref_condition`](#ref_condition). [default: pairwise]
  --msstats_threshold                     [number]  The threshold value for differential expressed proteins in MSstats plots based on adjusted p-value [default: 0.05]
  --add_triqler_output                    [boolean] Also create an output in Triqler's format for an alternative manual post-processing with that tool
  --msstatslfq_feature_subset_protein     [string]  Which features to use for quantification per protein: 'top3' or 'highQuality' which removes outliers only  (accepted: top3, highQuality) [default: top3]
  --msstatslfq_quant_summary_method       [string]  which summary method to use: 'TMP' (Tukey's median polish) or 'linear' (linear mixed model)  (accepted: TMP, linear) [default: TMP]
  --msstats_remove_one_feat_prot          [boolean] Omit proteins with only one quantified feature? [default: true]
  --msstatslfq_removeFewMeasurements      [boolean] Keep features with only one or two measurements across runs? [default: true]
  --msstatsiso_useunique_peptide          [boolean] Use unique peptide for each protein [default: true]
  --msstatsiso_rmpsm_withfewmea_withinrun [boolean] Remove the features that have 1 or 2 measurements within each run [default: true]
  --msstatsiso_summaryformultiple_psm     [string]  select the feature with the largest summmation or maximal value  (accepted: sum, max) [default: sum]
  --msstatsiso_summarization_method       [string]  summarization methods to protein-level can be perfomed  (accepted: msstats, MedianPolish, Median, LogSum) [default: msstats]
  --msstatsiso_global_norm                [boolean] Reference channel based normalization between MS runs on protein level data? [default: true]
  --msstatsiso_remove_norm_channel        [boolean] Remove 'Norm' channels from protein level data [default: true]
  --msstatsiso_reference_normalization    [boolean] Reference channel based normalization between MS runs on protein level data [default: true]
  --msstats_plot_profile_qc               [boolean] Export MSstats profile QC plots including all proteins

Quality control
  --enable_pmultiqc                       [boolean] Enable generation of pmultiqc report? default: true [default: true]
  --pmultiqc_idxml_skip                   [boolean] Skip idXML files (do not generate search engine scores) in pmultiqc report? default: 'true' [default: true]
  --contaminant_string                    [string]  Contaminant affix string for pmultiqc report. This parameter maps to --contaminant_affix in pmultiqc. default: 'CONT' [default: CONT]

Generic options
  --skip_table_plots                      [boolean]         Skip protein/peptide table plots with pmultiqc for large dataset.
  --multiqc_methods_description           [string]          Custom MultiQC yaml file containing HTML including a methods description.
  --help                                  [boolean, string] Display the help message.
  --help_full                             [boolean]         Display the full detailed help message.
  --show_hidden                           [boolean]         Display hidden parameters in the help message (only works when --help or --help_full are provided).

 !! Hiding 50 param(s), use the `--showHidden` parameter to show them !!
------------------------------------------------------
```

</details>

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

Now it is safe to delete the `work` folder.

```bash
rm -r results/PXD040621 work
```

You could also delete the results folder in case you need some more storage

```
rm -r results/PXD040621
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
rm -r /vscode/extensionsCache/
rm -r /vscode/serverCache/
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
