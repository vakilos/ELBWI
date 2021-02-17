#!/bin/bash
module load seqmagick
source activate seqmagick
seqmagick convert --squeeze alignedOtus.fasta squeezed.fasta
conda deactivate
module unload seqmagick
