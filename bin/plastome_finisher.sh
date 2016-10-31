#!/bin/bash

perl Fast-Plast/bin/sequence_based_ir_id.pl $1_afin.fa $1 2

<<<<<<< HEAD
blastn -query $1_regions_split.fsa -db  Fast-Plast/bin/position_genes.fsa -evalue 1e-10 -outfmt 6 > $1.split2.blastn

perl Fast-Plast/bin/orientate_plastome.pl $1_regions_split.fsa $1.split2.blastn $1
=======
blastn -query $1_regions_split2.fsa -db  Fast-Plast_master/bin/position_genes.fsa -evalue 1e-10 -outfmt 6 > $1.split2.blastn

perl Fast-Plast_master/bin/orientate_plastome.pl $1_regions_split2.fsa $1.split2.blastn $1
>>>>>>> c402b7d9cb99a60235fffcbdb315710dc0dbd4cd

mv $1_FULLCP.fsa $1_2_FULLCP.fsa

perl Fast-Plast/bin/sequence_based_ir_id.pl $1_afin.fa $1 3

<<<<<<< HEAD
blastn -query $1_regions_split.fsa -db Fast-Plast/bin/position_genes.fsa -evalue 1e-10 -outfmt 6 > $1.split3.blastn

perl Fast-Plast/bin/orientate_plastome.pl $1_regions_split.fsa $1.split3.blastn $1
=======
blastn -query $1_regions_split3.fsa -db Fast-Plast_master/bin/position_genes.fsa -evalue 1e-10 -outfmt 6 > $1.split3.blastn

perl Fast-Plast_master/bin/orientate_plastome.pl $1_regions_split3.fsa $1.split3.blastn $1
>>>>>>> c402b7d9cb99a60235fffcbdb315710dc0dbd4cd

mv $1_FULLCP.fsa $1_3_FULLCP.fsa


