#!/bin/bash
#SBATCH --mem=20000
#SBATCH --mail-type=END,FAIL
#SBATCH --output="slurm-%x-%j.out"


# This script imports and QCs the data using DADA2
echo "QIIME2 bash script for V6 samples started running at: "; date

module load qiime2/2020.11.1
module list

# Put metadata file name here
METADATA="Mix9Conv16Samplicon_metadata.txt"
MANIFEST="Mix9Conv16Samplicon_manifest.txt"
CLASSIFIER="db/qiime2/2020.11/classifier_silva-138-99-V6.qza"

# Import data into QIIME
# Paired-end, based on sample-manifest.csv
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $MANIFEST --output-path paired-end-sequences.qza \
  --input-format PairedEndFastqManifestPhred33V2

# QC using dada2
qiime dada2 denoise-paired --verbose --i-demultiplexed-seqs paired-end-sequences.qza \
  --p-trunc-len-r 80  --p-trunc-len-f 80 \
  --p-trim-left-r 19 --p-trim-left-f 19 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --p-n-reads-learn 10000000

# Summarize feature table and sequences
qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file $METADATA
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

# This script assigns taxonomy based on the imported database (see X_qiime2gettaxonomydb.sh)

qiime feature-classifier classify-sklearn \
  --i-classifier $CLASSIFIER \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv
qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file $METADATA \
  --o-visualization taxa-bar-plots.qzv

# This script calculated phylogenetic trees for the data

# align and mask sequences
qiime alignment mafft \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza
qiime alignment mask \
  --i-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza

# calculate tree
qiime phylogeny fasttree \
  --i-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza
qiime phylogeny midpoint-root \
  --i-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv
qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file $METADATA \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file $METADATA \
  --m-metadata-column Timepoint \
  --o-visualization core-metrics-results/unweighted-unifrac-timepoint-significance.qzv \
  --p-pairwise

# This script calculates the rarefaction curve for the data
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 600 \
  --m-metadata-file $METADATA \
  --o-visualization alpha-rarefaction.qzv

echo "END $(date)"
