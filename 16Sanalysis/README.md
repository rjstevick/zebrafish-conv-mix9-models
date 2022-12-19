# 16Sanalysis
This folder contains scripts used to perform QIIME2 analysis on the 16S rRNA amplicon data, the QIIME2 output artifacts, and the Rmd script to reproduce the figures in the manuscript. Software used: QIIME2 v2020.4 and R v3.4.1

## [scripts/](scripts/)
- [Mix9Conv16Samplicon_metadata.txt](scripts/Mix9Conv16Samplicon_metadata.txt) - metadata file formatted for QIIME2 import during `qiime feature-table summarize`
- [Mix9ConvControls16Samplicon_analysis.Rmd](scripts/Mix9ConvControls16Samplicon_analysis.Rmd) and knitted outputs of analysis in R
- [00_qiime2script_16Sv6.sh](scripts/00_qiime2script_16Sv6.sh) - main QIIME2 script with all bash commands

## [qiime2output/](qiime2output/)
QIIME2 artifacts and processed files  

## [figures/](figures/)
Output from Rmd file in scripts/  
