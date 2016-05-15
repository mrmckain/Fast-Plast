#!/bin/bash

#perl /home/mmckain/DASH_Phylogeny/GSS_Plastomes/sequence_based_ir_id.pl $1_afin.fa $1
#perl identify_coverage_variation_assembly-fix.pl testing_k.txt $1_afin.fa.coverage_20kmer.txt $1_afin.fa
blastn -query $1_regions_split.fsa -db $2 -evalue 1e-10 -outfmt 6 > $1.blastn
perl ~/bin/orientate_plastome.pl $1_regions_split.fsa $1.blastn $1
