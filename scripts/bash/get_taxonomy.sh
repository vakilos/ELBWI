#!/bin/bash

module load qiime/2-2019.7
source activate qiime2-2019.7
qiime feature-classifier classify-sklearn --i-classifier 16S-silva-341f-785r-classifier.qza --i-reads asv_seqs.qza --o-classification taxonomy_f.qza

qiime feature-classifier classify-sklearn --p-read-orientation reverse-complement --i-classifier 16S-silva-341f-785r-classifier.qza --i-reads asv_seqs.qza --o-classification taxonomy_r.qza
