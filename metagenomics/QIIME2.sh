#!/bin/bash
source /opt/miniconda2/etc/profile.d/conda.sh
conda activate qiime2-2020.6

#our reads are demultiplexed paired-end data, need to be stored into .qza format to be able to analyse in QIIME2

dir="."
fastqd=${dir}/"subfastq"
greengenes="/home/public/16S/gg-13-8-99-515-806-nb-classifier.qza"
silva="/home/public/16S/silva-138-99-515-806-nb-classifier.qza"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $fastqd \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path demux-paired-end.qza

#Use dada2 which automatically does everything
#https://docs.qiime2.org/2020.6/tutorials/overview/#derep-denoise

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux-paired-end.qza \
  --p-trim-left-f 13 \
  --p-trim-left-r 13 \
  --p-trunc-len-f 240 \
  --p-trunc-len-r 200 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv

#generate a tree for phylogeny

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza


qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth 200 \
  --m-metadata-file metadata.txt \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file metadata.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file metadata.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Group \
  --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv \
  --p-pairwise

qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file metadata.txt \
  --p-custom-axes Day \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-days-since-experiment-start.qzv

#variable must be numeric, skipped

qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-max-depth 200 \
  --m-metadata-file metadata.txt \
  --o-visualization alpha-rarefaction.qzv

#you need the metadata which contains barcode information

#taxonomy generation

qiime feature-classifier classify-sklearn \
  --i-classifier $greengenes \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza

qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file metadata.txt \
  --o-visualization taxa-bar-plots.qzv

#make composition table

qiime composition add-pseudocount \
  --i-table table.qza \
  --o-composition-table comp-table.qza

qiime composition ancom \
  --i-table comp-table.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Group \
  --o-visualization ancom-subject.qzv

#can compare different taxnomy points

qiime taxa collapse \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --p-level 6 \
  --o-collapsed-table table-6.qza

qiime composition add-pseudocount \
  --i-table table-6.qza \
  --o-composition-table comp-table-6.qza

qiime composition ancom \
  --i-table comp-table-6.qza \
  --m-metadata-file metadata.txt \
  --m-metadata-column Group \
  --o-visualization 6-ancom-Group.qzv


#export

qiime tools export \
  --input-path table.qza \
  --output-path exported-feature-table

qiime tools export \
  --input-path unrooted-tree.qza \
  --output-path exported-tree

# ##using silva

# #taxonomy generation

# nohup qiime feature-classifier classify-sklearn \
#   --i-classifier $silva \
#   --i-reads rep-seqs.qza \
#   --o-classification taxonomy_silva.qza &

# nohup qiime metadata tabulate \
#   --m-input-file taxonomy_silva.qza \
#   --o-visualization taxonomy_silva.qzv &

# qiime taxa barplot \
#   --i-table table.qza \
#   --i-taxonomy taxonomy_silva.qza \
#   --m-metadata-file metadata.txt \
#   --o-visualization taxa-bar-plots_silva.qzv

# #can compare different taxnomy points

# qiime taxa collapse \
#   --i-table table.qza \
#   --i-taxonomy taxonomy_silva.qza \
#   --p-level 5 \
#   --o-collapsed-table table-5_silva.qza

# qiime composition add-pseudocount \
#   --i-table table-5_silva.qza \
#   --o-composition-table comp-table-5_silva.qza

# qiime composition ancom \
#   --i-table comp-table-5_silva.qza \
#   --m-metadata-file metadata.txt \
#   --m-metadata-column Group \
#   --o-visualization 5-ancom-Group_silva.qzv