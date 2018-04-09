Fast-Plast: Rapid de novo assembly and finishing for whole chloroplast genomes
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.973887.svg)](https://doi.org/10.5281/zenodo.973887)
=============
<b>Authors</b>: Michael R. McKain, <a href="https://github.com/afinit/afin">Mark Wilson</a><br>
</br>
Version 1.2.7<br>
</br>
<b>Contact</b>: https://github.com/mrmckain

<img src="https://github.com/mrmckain/Fast-Plast/blob/master/extras/Fast-Plast_Logo.png" width="256" alt="fast-plast_logo" align="middle">

<h1>Description</h1>

Fast-Plast is a pipeline that leverages existing and novel programs to quickly assemble, orient, and verify whole chloroplast genome sequences. For most datasets with sufficient data, Fast-Plast is able to produce a full-length de novo chloroplast genome assembly in approximately 30 minutes with no user mediation. In addition to a chloroplast sequence, Fast-Plast identifies chloroplast genes present in the final assembly.

Currently, Fast-Plast is written to accomodate Illumina data, although most data types could be used.

Fast-Plast uses a de novo assembly approach by combining the De Bruijn graph-based method of SPAdes with an iterative seed-based assembly implemented in afin to close gaps of contigs with low coverage. The pipeline then identifies regions from the quadripartite structure of the chloroplast genome, assigns identity, and orders them according to standard convention. A coverage analysis is then conducted to assess the quality of the final assembly. 


<h1>Dependencies</h1>

Fast-Plast requires a number of commonly used bioinformatics programs. We have included an installation script to help users properly prepare Fast-Plast for use.

* Perl 5.x which has been installed with threading enabled (see <a href="https://perlmaven.com/how-to-build-perl-from-source-code">How to build perl from source</a>, if needed).
* <a href="http://www.usadellab.org/cms/?page=trimmomatic">Trimmomatic</a> (initial read cleaning)
* <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a> (read reduction to only chloroplast-like reads)
* <a href="http://bioinf.spbau.ru/spades">SPAdes</a> (initial assembly)
* <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download">BLAST+</a> (multiple checks for gene content)
* <a href="https://github.com/nsoranzo/sspace_basic">SSPACE</a> (scaffolding when single contig not obtainable)
* <a href="https://sourceforge.net/projects/bowtie-bio/">Bowtie1</a> (required for SSPACE)

<h3>afin</h3>

* c++ complier with c++11 support and zlib.h.  zlib.h is a standard base library for most Unix systems but can be obtained <a href="http://www.zlib.net/">here</a>.

<h3>Coverage Analysis</h3>

* <a href="http://www.genome.umd.edu/jellyfish.html#Release">Jellyfish 2</a><br>
* R

Memory requirements will vary based on the size of your data set. Expect to use 1.5-2x the memory for the size of your reads files. If your data set is exceptionally large, we have found success in reducing the dataset to 5-10 million reads and running them through Fast-Plast.

<h1>Installation</h1>

Fast-Plast has been tested on Linux (CentOs 7) but should be compatible with any flavor that can handle dependencies. 

Clone the Github repository:

    git clone https://github.com/mrmckain/Fast-Plast.git

To install, run the INSTALL.pl script found in the Fast-Plast repository.

     perl INSTALL.pl

The installation script will walk you through installation.  If you already have the dependencies installed, you will need to give the full path to the executables. The installation script can also install all dependencies for you.  To do this, select "All" when prompted.  Dependencies will be installed in Fast-Plast/bin/. This script will also change paths in the control script, compile afin and jellyfish2, and unzip the default chloroplast genomes for mapping.

For advanced users, paths can be set directly in the fast-plast.pl file:

    ###directories
    my $FPROOT = "$FindBin::RealBin";
    my $AFIN_DIR = "$FPROOT/afin";
    my $COVERAGE_DIR = "$FPROOT/Coverage_Analysis";
    my $FPBIN = "$FPROOT/bin";
    my $TRIMMOMATIC;
    my $BOWTIE2;
    my $SPADES;
    my $BLAST;
    my $SSPACE;
    my $BOWTIE1;
    my $JELLYFISH;

Instructions for direct compilation of afin can be found <a href="https://github.com/mrmckain/Fast-Plast/tree/master/afin">here</a>. 

<h1>Input</h1>

<h3>Reads</h3>
Input files are in FASTQ format. The data are not expected to be adapter trimmed or quality filtered, though this will not impede the assembly. 

Fast-Plast was built for genome survey sequence (aka genome skimming or low-pass genome sequencing) data. Sequence capture data can be used but needs to be normalized first.  If your data are from a sequence capture experiment (aka target enrichments or anchored phylogenomics), we suggest using the normalization method packaged with <a href="https://github.com/trinityrnaseq/trinityrnaseq/wiki">Trinity</a>, <a href="http://ged.msu.edu/papers/2012-diginorm/">khmer</a>, or <a href="https://sourceforge.net/projects/bbmap/">bbnorm</a>. 

<h3>Bowtie Index</h3>

Fast-Plast is packaged with 1,021 whole chloroplast genomes from GenBank. These cover a wide range of diversity including marine algae, angiosperms, ferns, etc. To access these, use the order of your species in the --bowtie_index option. Fast-Plast will pull all members of that order and create a bowtie index. If your order is not present, Fast-Plast will use a representative sequence from all orders present. This option can be selected using "All" or "GenBank". "All" is the default.

Example:

    --bowtie_index Poales

This will use all available Poales plastomes in the data set (118).

Orders currently available in Fast-Plast (70):

	Alismatales		Dipsacales	Marchantiales
	Apiales			Ericales	Monomastigales
	Aquifoliales		Euglenales	Myrtales
	Araucariales		Eupodiscales	Nymphaeales
	Asparagales		Fabales		Orthotrichales
	Asterales		Fagales		Pinales
	Austrobaileyales	Fragilariales	Poales
	Bangiales		Fucales		Polypodiales
	Brassicales		Funariales	Proteales
	Bryopsidales		Garryales	Pyrenomonadales
	Buxales			Gentianales	Ranunculales
	Caryophyllales		Geraniales	Rosales
	Celastrales		Ginkgoales	Sapindales
	Chlamydomonadales	Gracilariales	Saxifragales
	Chloranthales		Hypnales	Solanales
	Chlorellales		Lamiales	Sphaeropleales
	Cornales		Laminariales	Takakiales
	Cucurbitales		Laurales	Ulvales
	Cupressales		Liliales	Vaucheriales
	Cyanidiales		Lycopodiales	Vitales
	Cyatheales		Magnoliales	Zingiberales
	Cycadales		Malpighiales	Zygnematales
	Desmidiales		Malvales	
	Dioscoreales		Mamiellales		


<h3>User Provided Bowtie Index</h3>

A user-made bowtie index can be provided using the --user_bowtie option. The full path to the index should be given. If this option is used, the --bowtie_index option will be ignored.

<h3>Name</h3>

The --name option should be used with each run. This simply gives a prefix to all files. The default is "Fast-Plast".

<h1>Output</h1>

Fast-Plast produces a number of files that allow the user to trace the steps of the pipeline. From the directory where Fast-Plast is called, three files and a new directory will be produced.  The directory will be named by the --name option. The three files include:

<b>name_Fast-Plast_Progress.log</b><br>
Gives time and results for each step in the pipeline. Information regarding paramters chosen based on reads (such as kmer size) and chloroplast gene content will be found here. 

<b>name_result_out.log</b><br>
Contains the STDOUT from all programs.

<b>name__results_error.log</b><br>
Contains the STDERR from all programs.

<h2>Directory Hierarchy</h2>
<b>name_Plastome_Summary.txt</b><br>
Provides information on the number of reads, reads mapped, assembly size, chloroplast region size, and average coverage of each region.

<h3>1_Trimmed_Reads</h3>
For paired-end data, three trimmed read files will be made. Files ending in *trimmed_P1.fq and *trimmed_P2.fq are still paired-end.  The files ending in *trimmed_UP.fq is single-end. If only single end files are used, then only the *trimmed_UP.fq file will be found.

<h3>2_Bowtie_Mapping</h3>
Fast-Plast created bowtie index files will be found in this directory. Reads that mapped to the bowtie index are in the files map_pair_hits.1.fq, map_pair_hits.2.fq, and map_hits.fq for paired-end and single-end respectively. The file name.sam is the standard bowtie2 mapping output but is not used.

<h3>3_Spades_Assembly</h3>
The directory "spades_iter1" will contain the SPAdes assembly and standard SPAdes output files.

<h3>4_Afin_Assembly</h3>
The file "filtered_spades_contigs.fsa" contains contigs from the SPAdes assemblies that fall within the range of minus one standard deviation of the weigthed mean coverage to plus 2.5 standard deviations. 

Output from afin is the files *_afin_iter0.fa, *_afin_iter1.fa, *_afin_iter2.fa, and *_afin.log.  The log file demonstrates the steps afin took in the extension and assembly process. 

Chloroplast_gene_composition_of_afin_contigs_nested_removed.txt contains information regarding the chloroplast gene content of the each contig found *_afin_iter2.fa.

If only single end data were used and more than one contig is found, Fast-Plast will quit here. The final output will be available in Final_Assembly.

<h4>Scaffolding</h4>
If more than one contig is found in the file afin assembly and paired-end reads were used, SSPACE will be invoked to attempt scaffolding of contigs.  These results will be found here. If more than one contig/scaffold is present in the final output, this will be the last step of the pipeline and results will be found in Final_Assembly.

<h3>5_Plastome_Finishing</h3>
If a single contig was found through either contig assembly or scaffolding, this directory will be created to contain files associated with identification of the large single copy, small single copy, and inverted repeats. If these are not present, this will be the last step of the pipeline and results of the single (unorientated) contig from the assembly steps will be found in Final_Assembly.

<h3>Final_Assembly</h3>
The final assembly will be found in this directory. If the pipeline was able to fully assemble and orientate the plastome, the files *_CP_pieces.fsa (plastome split into LSC, SSC, and IR), *_FULLCP.fsa (final assembly), and Chloroplast_gene_composition_of_final_contigs.txt (all chloroplast genes found in final assembly) will be present. The file Chloroplast_gene_composition_of_final_contigs.txt will always be made for the final assembly regardless of where Fast-Plast stops.

<h3>Coverage_Analysis</h3>
The coverage analysis option should also be used to ensure accurate assembly of the plastome. Multiple files associated with the coverage estimation process will be present in this directory. The three most important files are:

<b>name.coverage_25kmer.txt</b><br>
Contains 25-mer sequence, start position, and coverage across final assembly.

<b>name_coverage.pdf</b><br>
Graphical representation of name.coverage_25kmer.txt. Red circles indicate a coverage of 0 and potential assembly issue.

<b>name_problem_regions_plastid_assembly.txt</b><br>
Identified stretches of the assembly greater than 25 base pairs that have a coverage of 0. If this file is empty, the assembly is accepted.

Other files include the mapped reads from the Bowtie2 run (map_hits*) and those associated with Bowtie2 and Jellyfish.

<h3>4.5_Reassemble_Low_Coverage</h3>
If regions of low coverage are identified after the Coverage Analysis, then these regions are removed, the contig broken into pieces, and reassembled from the afin step. All of the reassembly steps (Afin Assembly and Plastome Finishing) will be conducted in this directory.

<h3>Final_Assembly_Fixed_Low_Coverage</h3>
The final reassembled plastome, chloroplast regions, and chloroplast gene recovery will be in this directory.

<h3>Coverage_Analysis_Reassembly</h3>
Coverage analysis and results of the reassembled plastome will be in this directory.

<h1>Usage</h1>

<h3>General Syntax</h3>

<code>fast-plast.pl [-1 <paired_end_file1> -2 <paired_end_file2> || --single <single_end_file>] -name <sample_name> [options]</code>


<h3>Example with Paired-End Data</h3>

	perl fast-plast.pl -1 /home/mmckain/Sequence_Vault/Washburn_Data/37_Urochloa_fusca_42940_RPGH_AGTCAA_L005_R1_001.fastq.gz -2 /home/mmckain/Sequence_Vault/Washburn_Data/37_Urochloa_fusca_42940_RPGH_AGTCAA_L005_R2_001.fastq.gz --name Urochloa_fusca-37 --bowtie_index Poales --coverage_analysis --clean light

In this example, one pair-end library is being used for assembly. The default adapters (NEB) are used for trimming, Poales species are used for the Bowtie2 index, the coverage analysis is invoked, and a light cleaning is done after completition.

<h3>Example with Single End Data</h3>

	perl fast-plast.pl --single /home/mmckain/Sequence_Vault/Andropogoneae_GSS/Chionachne_koenigii-TK057/K17_GTCCGC_L006_R1_001.fastq.gz --name Chionachne_koenigii-TK057 --adapters TruSeq --bowtie_index All --coverage_analysis

In this example, one single end library is being used for assembly. The TruSeq adapters are used for trimming, a representative from each order is used for the Bowtie2 index, the coverage analysis is invoked, and no cleaning is done.

<h3>Example with Mixed Libraries</h3>

	perl fast-plast.pl -1 /home/mmckain/Sequence_Vault/Andropogoneae_GSS/Monocymbium_ceresiiforme-TK203/TK203_GCTACGCT-AGAGTAGA_Both_R1.fastq.gz -2 /home/mmckain/Sequence_Vault/Andropogoneae_GSS/Monocymbium_ceresiiforme-TK203/TK203_GCTACGCT-AGAGTAGA_Both_R2.fastq.gz --single /home/mmckain/Sequence_Vault/Andropogoneae_GSS/Monocymbium_ceresiiforme-TK203/Monocymbium_ceresiiforme-TK203_GA_122949-34_S42_R1_001.fastq --name Monocymbium_ceresiiforme-TK203 --user_bowtie /home/mmckain/Andropogoneae_Plastomes/FINISHED/androplast --coverage_analysis --clean deep

In this example, one single end library and one paired-end library are being used for assembly. The default adapters (NEB) are used for trimming, a user-defined Bowtie2 index (base name given) is used for the Bowtie2 index, the coverage analysis is invoked, and a deep cleaning is done after completion.

<h3>Example with Coverage Analysis Only</h3>

    perl ~/bin/Fast-Plast/fast-plast.pl -1 /home/mmckain/Sequence_Vault/MSU_HiSeq4000_05032017/20170428_DNASeq_PE/20170428_DNASeq_PE/TK686R_S48_L005_R1_001.fastq.gz -2 /home/mmckain/Sequence_Vault/MSU_HiSeq4000_05032017/20170428_DNASeq_PE/20170428_DNASeq_PE/TK686R_S48_L005_R2_001.fastq.gz --name Schizachyrium_scoparium-TK686R --only_coverage /home/mmckain/DASH_Phylogeny/GSS_Plastomes/Schizachyrium_scoparium-TK686R/Schizachyrium_scoparium-TK686R/Final_Assembly/Schizachyrium_scoparium-TK686R_FULLCP.fsa --min_coverage 2 &

In this example, a paired end library is being used.  The --only_coverage option is used with the path to the fasta file of the chloroplast genome provided. A minimum coverage of 2 is used.



Definitions:

		-1 <filenames>		File with forward paired-end reads. Multiple files can be designated with a comma-delimited list.
					Read files should be in matching order with other paired end files.
		-2 <filenames>		File with reverse paired-end reads. Multiple files can be designated with a comma-delimited list.
					Read files should be in matching order with other paired end files.
		--single <filenames>	File with unpaired reads. Multiple files can be designated with a comma-delimited list.

		**PAIRED END AND SINGLE END FILES CAN BE PROVIDED SIMULTANEOUSLY.**

		-n <sample_name>	Name for current assembly. We suggest a species name/accession combination as Fast-Plast will use 
					this name as the FASTA ID in the final assembly. [Default = Fast-Plast]

		Advanced options:

		--subsample		Number of reads to subsample. Reads will be evenly pulled from all files. 
		--threads		Number of threads used by Fast-Plast.  [Default = 4]
		--min_coverage		Lowest acceptable coverage for 25-mer sliding window during coverage analysis. [Default = 0.25 * Average coverage]
		--adapters		[NEB|Nextera|TruSeq] Files of adapters used in making sequencing library. NEB, Nextera, and TruSeq options 
					available. Also accepts the path to a user created FASTA file of adapters.[Default = NEB]
		--bowtie_index		Taxonomic order of the sequenced species to pick references for bowtie2 indices. If the order is in the database, 
					then all available samples for that order will be used. If order does not exist in database or the terms "all" or "GenBank" are given, one 
					exemplar from each available order is used to build the Bowtie2 indices. Users may also specify multiple taxa separated by commas (","). Any of the following taxonomic levels is accepted for bowtie_index: genus, species epithet, tribe (if applicable), subfamily (if applicable), family, and order. [default="All"]
		--user_bowtie		User supplied bowtie2 indices. If this option is used, bowtie_index is ignored.
		--coverage_analysis 	Flag to run the coverage analysis of a final chloroplast assembly.[Recommended]
        	--skip          	Flag to skip trimming. Must include option "trim". [--skip trim]
		--only_coverage         Option allows user to run coverage analysis directly on a provided chloroplast genome. [requires: 
                                    read files, chloroplast genome sequence]
            --clean 		[light|deep] The "light" option will remove all bowtie indices, BLAST databases, SAM files,
					Jellyfish dumps, and Jellyfish kmer files. The "deep" option will remove all directories except for the 
					Final Assembly and Coverage Analysis directories. All files in the "light" option will also be removed. 
					Clean will only be invoked if a fully successful assembly is made. 
		--posgenes		User defined genes for identification of single copy/IR regions and orientation. Useful when major 
					rearrangments are present in user plastomes. This is a fairly advanced option; please contact if you are interested
					in using it.



## Changelog

* 07-Jan_2018 Fast-Plast v.1.2.7 <br>
    --added capacity for multiple taxa to be used for bowtie_index. These should be a comma separated list. Any of the following taxonomic levels is accepted for bowtie_index: genus, species epithet, tribe (if applicable), subfamily (if applicable), family, and order.
    <br>

* 01-June_2017  Fast-Plast v.1.2.6 <br>
    --skip parameter added. Allows users to skip the trimming step when "trim" option is include. Syntax: --skip trim
    <br>

* 21-May-2017 Fast-Plast v.1.2.5 <br>
    --only_coverage option added.  Allows users to provide a chloroplast genome and reads to run the coverage analysis portion of Fast-Plast.
    <br>
    --Coverage used by Fast-Plast is printed to log and summary files.

* 24-April-2017 Fast-Plast v.1.2.4 <br>
    --N's allowed in coverage and IR/SC identification from scaffolding. 
    <br>
    --Faster runs of afin after initial coverage analysis. 
    <br>
    --Minimum coverage for coverage analysis no longer defaulted as 0.  Now estimated as 0.25 of average coverage across plastome.
    <br>
    --General bug fixes accuracy.

* 31-March-2017 Fast-Plast v.1.2.3 <br>
	--Subsampling option added to allow for faster completion.

* 16-March-2017 Fast-Plast v.1.2.2 <br>
	--Nextera and TruSeq adapter options added from Trimmomatic adapter set.
	<br>
	--Cleaning option added to remove large files.
	<br>
	--More genes added to angiosperm chloroplast gene list.

* 10-March-2017 Fast-Plast v.1.2.0<br>
	--Added orignal code to identify low coverage areas, split assembly, and reassemble. This portion of the pipeline is still in alpha.<br>
	--Summary gives information on numbers of trimmed reads, reads mapped, size of plastome and its regions, and average coverage of regions.

* 07-March-2017	First Release Fast-Plast v.1.1.0




<h1>References</h1>

* Bankevich, A., S. Nurk, D. Antipov, A. A. Gurevich, M. Dvorkin, A. S. Kulikov, V. M. Lesin, S. I. Niklenko, S. Pham, A. D. Prjibelski, A. V. Pyshkin, A. V. Sirotkin, N. Vyahhi, G. Tesler, M.A. Alexkseyev, and P. A. Pevzner. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. <i>J. Comp. Biol.</i>, 19(5):455-477.
* Bolger, A. M., M. Lohse, and B. Usadel. 2014. Trimmomatic: A flexible trimmer for Illumina Sequence Data. <i>Bioinformatics</i>, btu170.
* Langmead, B. and S. Salzberg. 2012. Fast gapped-read alignment with Bowtie 2. <i>Nature Methods</i>, 9:357-359.
* Marçais, G. and C. Kingsford. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. <i>Bioinformatics</i>, 27:764-770.
