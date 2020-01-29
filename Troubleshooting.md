<h1>Examples of fixing issues:</h1>

Look at Camelina_microcarpa_ 101_afin_iter2.fa found in 

If you run the script “sequence_based_ir_id.pl” from the Fast-Plast bin (see below), you can split an assembly into putative portions (large single copy, inverted repeat, and small single copy) of the chloroplast genome.

perl ~/Github/Fast-Plast/bin/sequence_based_ir_id.pl Camelina_microcarpa_ 101_afin_iter2.fa Camelina_microcarpa_101 0

Syntax:
perl sequence_based_ir_id.pl sequence_file_to_split name_base skip_parameter

The skip parameter is the number of the occurrence of “single copy” that gets remove. A “0” means none, “1” means the first, and so on. The reason this exists is that it is not always possible to ID the full IR sequence. If you have an assembly where the real identity of the assembly looks like this: partial_IR-SSC-full_IR-LSC-partial_IR, you will likely get a result from the script that looks like this: IR-SC-IR-SC-IR-SC-IR. The script only IDs single copy and IR, but you will notice there is an extra SC (single copy) in the middle of the full IR. This is region not covered by the two partial IRs, so it cannot be recognized as an inverted repeat. In this example, using a skip parameter of “2” would result in an accurate depiction of the composition. Fast-Plast does this for you by running this script multiple times using 0, 1, and 2. 

After running the script with skip parameter “0”, I get a fasta file with the following info:

Camelina_microcarpa_101_regions_split0.fsa:>ir_5295-31778
Camelina_microcarpa_101_regions_split0.fsa:>sc_31779-114049
Camelina_microcarpa_101_regions_split0.fsa:>ir_114050-140533
Camelina_microcarpa_101_regions_split0.fsa:>sc_140534-158353
Camelina_microcarpa_101_regions_split0.fsa:>ir_158354-184837

This represents: full_IR-LSC-full_IR-SSC-full_IR

This assembly is big enough that it loops around the entire plastomes and then extends in the IR for a third time. You will notice that the IR size is 26,483, the LSC size is 82,270, and the SSC is 17,819. These sizes are similar to Brassica juncea (LSC 83,286; IR: 26,212; SSC: 17,775) and are likely accurate.

To proceed, you can delete the entries :>ir_5295-31778 and >ir_158354-184837, leaving one IR copy and the two single copies. You then need to do the following:

blastn -query edited_afin_file -db Fast-Plast/bin/Angiosperm_Chloroplast_Genes.fsa -outfmt 6 -evalue 1e-10 -max_target_seqs 1000000 > base_name.blastn

*Make sure the blast database for Angiosperm_Chloroplast_Genes.fsa has been made.

With those results, run: 

perl Fast-Plast/bin/orientate_plastome_v.2.0.pl edited_afin_file base_name.blastn base_name

This will give you the files:

Base_name_FULLCP.fsa
Base_name_CP_pieces.fsa

The FULLCP file can then be used to run a coverage analysis with Fast-Plast to check the assembly. 
