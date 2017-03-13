Fast-Plast: Rapid de novo assembly and finishing for whole chloroplast genomes
=============
<b>Authors</b>: Michael R. McKain, <a href="https://github.com/afinit/afin">Mark Wilson</a><br>
</br>
Version 1.2.1<br>
</br>
<b>Contact</b>: https://github.com/mrmckain

<img src="https://github.com/mrmckain/Fast-Plast/blob/master/extras/Fast-Plast_Logo.png" width="256" alt="fast-plast_logo" align="middle">

<h1>Description</h1>

Fast-Plast is a pipeline that leverages existing and novel programs to quickly assemble, orient, and verify whole chloroplast genome sequences. For most datasets with sufficient data, Fast-Plast is able to produce a full-length de novo chloroplast genome assembly in approximately 30 minutes with no user mediation. In addition to a chloroplast sequence, Fast-Plast identifies chloroplast genes present in the final assembly.

Currently, Fast-Plast is written to accomodate Illumina data, though most data types could be used.

Fast-Plast uses a de novo assembly approach by combining the de bruijn graph-based method of SPAdes with an iterative seed-based assembly implemented in afin to close gaps of contigs with low coverage. The pipeline then identifies regions from the quadripartite structure of the chloroplast genome, assigns identity, and orders them according to standard convention. A coverage analysis is then conducted to assess the quality of the final assembly. 


<h1>Dependencies</h1>

Fast-Plast requires a number of commonly used bioinformatics programs. We have included an installation script to help users properly prepare Fast-Plast for use.

* Perl 5+
* <a href="http://www.usadellab.org/cms/?page=trimmomatic">Trimmomatic</a> (initial read cleaning)
* <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a> (read reduction to only chloroplast-like reads)
* <a href="http://bioinf.spbau.ru/spades">SPAdes</a> (initial assembly)
* <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download">BLAST+</a> (multiple checks for gene content)
* <a href="https://github.com/nsoranzo/sspace_basic">SSPACE</a> (scaffolding when single contig not obtainable)
* <a href="https://sourceforge.net/projects/bowtie-bio/">Bowtie1</a> (required for SSPACE)

<h3>afin</h3>
* c++ complier with c++11 support and zlib.h.  zlib.h is a standard base library for most Unix systems but can be obtained <a href="http://www.zlib.net/">here</a>.

<h3>Coverage Analysis</h3>
* <a href="http://www.genome.umd.edu/jellyfish.html#Release">Jellyfish 2</a>
* R

Memory requirements will vary based on the size of your data set. Expect to use 1.5-2x the memory for the size of your reads files. If your data set is exceptionally large, we have found success in reducing the dataset to 50 million reads and running them through Fast-Plast.

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
Input files are in FASTQ format. The data is not expected to be adapter trimmed or quality filtered, though this will not impede the assembly. 

Fast-Plast was built for genome survey sequence (aka genome skimming or low-pass genome sequencing) data. Sequence capture data can be used but needs to be normalized first.  If your data is from a sequence capture experiment (aka target enrichments or anchored phylogenomics), we suggest using the normalization method packaged with <a href="https://github.com/trinityrnaseq/trinityrnaseq/wiki">Trinity</a>, <a href="http://ged.msu.edu/papers/2012-diginorm/">khmer</a>, or <a href="https://sourceforge.net/projects/bbmap/">bbnorm</a>. 

<h3>Bowtie Index</h3>

Fast-Plast is packaged with 1,020 whole chloroplast genomes from GenBank. These cover a wide range of diversity including marine algae, angiosperms, ferns, etc. To access these, use the order of your species in the --bowtie_index option. Fast-Plast will pull all members of that order and create a bowtie index. If your order is not present, Fast-Plast will use a representative sequence from all orders present. This option can be selected using "All" or "GenBank". "All" is the default.

Example:

    --bowtie_index Poales

This will use all available Poales plastomes in data set (117).

Orders currently available in Fast-Plast (68):

    Apiales				Ericales		Marchantiales
    Aquifoliales		Euglenales		Monomastigales
    Araucariales		Eupodiscales	Myrtales
    Asparagales			Fabales			Nymphaeales
    Asterales			Fagales			Orthotrichales
    Austrobaileyales	Fragilariales	Pinales
    Brassicales			Fucales			Poales
    Bryopsidales		Funariales		Polypodiales
    Buxales				Garryales		Proteales
    Caryophyllales		Gentianales		Pyrenomonadales
    Celastrales			Geraniales		Ranunculales
    Chlamydomonadales	Ginkgoales		Rosales
    Chloranthales		Gracilariales	Sapindales
    Chlorellales		Hypnales		Saxifragales
    Cornales			Lamiales		Solanales
    Cucurbitales		Laminariales	Sphaeropleales
    Cupressales			Laurales		Takakiales
    Cyanidiales			Liliales		Ulvales
    Cyatheales			Lycopodiales	Vaucheriales
    Cycadales			Magnoliales		Vitales
    Desmidiales			Malpighiales	Zingiberales
    Dioscoreales		Malvales		Zygnematales
    Dipsacales			Mamiellales	


<h3>User Provided Bowtie Index</h3>

A user-made bowtie index can be provided using the --user_bowtie option. The full path to the index should be given. If this option is used, the --bowtie_index option will be ignored.

<h3>Name</h3>

The --name option must be used with each run. No default is set. This simply gives a prefix to all files.

<h1>Output</h1>

Fast-Plast produces a number of files that allow the user to trace the steps of the pipeline. From the directory where Fast-Plast is called, three files and direcotry will be produced.  The directory will be named by the --name option. The three files include:

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
For paired-end data, four trimmed read files will be made. Files ending in *trimmed_P1.fq and *trimmed_P2.fq are still paired-end.  The files ending in *trimmed_UP.fq is single-end. If only single end files are used, then only the *trimmed_UP.fq file will be found.

<h3>2_Bowtie_Mapping</h3>
Fast-Plast created bowtie index files will be found in this directory. Reads that mapped to the bowtie index are in the files map_pair_hits.1.fq, map_pair_hits.2.fq, and map_hits.fq for paired-end and single-end respectively. The file name.sam is the standard bowtie2 mapping output but is not used.

<h3>3_Spades_Assembly</h3>
The directory "spades_iter1" will contain the SPAdes assembly and standard SPAdes output files.

<h3>4_Afin_Assembly</h3>
The file "filtered_spades_contigs.fsa" contains contigs from the SPAdes assemblies that fall within the range of minus one standard deviation of the weigthed mean coverage to plus 2.5 standard deviations. 

Output from afin is the files *_afin_iter0.fa, *_afin_iter1.fa, *_afin_iter2.fa, and *_afin.log.  The log file that shows the steps afin took in the extention and assembly process. 

Chloroplast_gene_composition_of_afin_contigs_nested_removed.txt contains information regarding the chloroplast gene content of the each contig found *_afin_iter2.fa.

If only single end data were used and more than one contig is found, Fast-Plast will quit here. The final output will be available in Final_Assembly.

<h4>Scaffolding</h4>
If more than one contig is found in the file afin assembly and paired-end reads were used, SSPACE will be invoked to attempt scaffolding of contigs.  These results will be found here. If more than one contig/scaffold is present in the final output, this will be the last step of the pipeline and results will be found in Final_Assembly.

<h3>5_Plastome_Finishing</h3>
If a single contig was found through either contig assembly or scaffolding, this directory will be created to contain files associated with identification of the large single copy, small single copy, and inverted repeats. If these are not present, this will be the last step of the pipeline and results of the single (unorientated) contig from the assembly steps will be found in Final_Assembly.

<h3>Final_Assembly</h3>
The final assembly will be found in this directory. If the pipeline was able to full assemble and orientate the plastome, the files *_CP_pieces.fsa (plastome split into LSC, SSC, and IR), *_FULLCP.fsa (final assembly), and Chloroplast_gene_composition_of_final_contigs.txt (all chloroplast genes found in final assembly) will be present. The file Chloroplast_gene_composition_of_final_contigs.txt will always be made for the final assembly regardless of where Fast-Plast stops.

<h3>Coverage_Analysis</h3>
The coverage analysis option should also be used to ensure accurate assemble of the plastome. Multiple files associated with the coverage estimation process will be present in this directory. The three most important files are:

<b>name.coverage_25kmer.txt</b><br>
Contains 25-mer sequence, start position, and coverage across final assembly.

<b>name_coverage.pdf</b><br>
Graphical representation of name.coverage_25kmer.txt. Red circles indicate a coverage of 0 and potential assembly issue.

<b>name_problem_regions_plastid_assembly.txt</b><br>
Identified stretches of the assembly greater than 25 base pairs that have a coverage of 0. If this file is empty, the assembly is accepted.

<h1>Usage</h1>

<code> fast-plast.pl [-1 <paired_end_file1> -2 <paired_end_file2> || -single <singe_end_file>] -name <sample_name> [options]  </code>

Definition:

        -1 <filenames>		File with forward paired-end reads. Multiple files can be designated with a comma-delimited list. 
							Read files should be in matching order with other paired end files.
		-2 <filenames>		File with reverse paired-end reads. Multiple files can be designated with a comma-delimited list. 
							Read files should be in matching order with other paired end files.
		-s <filenames>		File with unpaired reads. Multiple files can be designated with a comma-delimited list.

		**PAIRED END AND SINGLE END FILES CAN BE PROVIDED SIMULTAENOUSLY.**

		-n <sample_name>	Name for current assembly. We suggest a species name/accession combination as Fast-Plast will use 
							this name as the FASTA ID in the final assembly.

		Advanced options:

		--threads			Number of threads used by Fast-Plast.  [Default = 4]
		--adapters			Files of adapters used in making sequencing library. These should be in FASTA format. [Default = NEB-PE]
		--bowtie_index		Order for sample to draw references for mapping. If order exists, then all available samples for that 
							order will be used. If order does not exist in default set or the terms "all" or "GenBank" are given, 
							one exemplar from each available order is used to build the Bowtie2 indicies. [default="All"]
		--user_bowtie		User supplied bowtie2 indices. If this option is used, bowtie_index is ignored.
		--posgenes			User defined genes for identification of single copy/IR regions and orientation. Useful when major
							rearrangments are present in user plastomes.
		--coverage_analysis Flag to run the coverage analysis of a final chloroplast assembly.[Recommended]




## Changelog
* 07-March-2017	First Release Fast-Plast v.1.1.0

* 10-March-2017 Fast-Plast v.1.2.0
	--Added orignal code to identify low coverage areas, split assembly, and reassemble. This portion of the pipeline is still in alpha.<br>
	--Summary gives information on numbers of trimmed reads, reads mapped, size of plastome and its regions, and average coverage of regions.


<h1>References</h1>

* Bankevich, A., S. Nurk, D. Antipov, A. A. Gurevich, M. Dvorkin, A. S. Kulikov, V. M. Lesin, S. I. Niklenko, S. Pham, A. D. Prjibelski, A. V. Pyshkin, A. V. Sirotkin, N. Vyahhi, G. Tesler, M.A. Alexkseyev, and P. A. Pevzner. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. <i>J. Comp. Biol.</i>, 19(5):455-477.
* Bolger, A. M., M. Lohse, and B. Usadel. 2014. Trimmomatic: A flexible trimmer for Illumina Sequence Data. <i>Bioinformatics</i>, btu170.
* Langmead, B. and S. Salzberg. 2012. Fast gapped-read alignment with Bowtie 2. <i>Nature Methods</i>, 9:357-359.
* Marçais, G. and C. Kingsford. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. <i>Bioinformatics</i>, 27:764-770.
