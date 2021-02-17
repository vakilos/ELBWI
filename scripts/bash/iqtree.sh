#!/bin/bash
module load iqtree
iqtree -nt 2 -s squeezed.fasta
mv squeezed.fasta.treefile asv_iq_tree.treefile

