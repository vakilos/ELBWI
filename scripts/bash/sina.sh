#!/bin/bash
gunzip SILVA_138.1_SSURef_opt.arb.gz
sina -i final_asv_seqs_clean.fasta -o alignedOtus.fasta --db SILVA_138.1_SSURef_opt.arb


