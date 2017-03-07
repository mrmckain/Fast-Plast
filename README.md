Fast-Plast: Rapid de novo assembly and finishing for whole chloroplast genomes
=============
<b>Authors</b>: Michael R. McKain, <a href="https://github.com/afinit/afin">Mark Wilson</a><br>
</br>
Version 1.1.0<br>
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

<br>
<h3>afin</h3>
* c++ complier with c++11 support and zlib.h.  zlib.h is a standard base library for most Unix systems but can be obtained <a href="http://www.zlib.net/">here</a>.

<br>
<h3>Coverage Analysis</h3>
* <a href="http://www.genome.umd.edu/jellyfish.html#Release">Jellyfish 2</a>
* R

Memory requirements will vary based on the size of your data set. Expect to use 1.5-2x the memory for the size of your reads files. If your data set is exceptionally large, we have found success in reducing the dataset to 50 million reads and running them through Fast-Plast.

<h1>Installation</h1>

Fast-Plast has been tested on Linux (CentOs 7), but should be compatible with any flavor that can handle dependencies. 

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

Fast-Plast is packaged with 1,020 whole chloroplast genomes from GenBank. These cover a wide range of diversity including marine algae, angiosperms, ferns, etc. To access these, use the order of your species in the --bowtie_index option. Fast-Plast will pull all members of that order and create a bowtie index. If your order is not present, Fast-Plast will use a representative sequence from all orders present. This option can be selected using "All" or "GenBank".

Orders currently available in Fast-Plast:

<ol
	><li>Apiales</li
	><li>Aquifoliales</li
	><li>Araucariales
	><li>Asparagales
	><li>Asterales
	><li>Austrobaileyales
	><li>Brassicales
	><li>Bryopsidales
	><li>Buxales
	><li>Caryophyllales
	><li>Celastrales
	><li>Chlamydomonadales
	><li>Chloranthales
	><li>Chlorellales
	><li>Cornales
	><li>Cucurbitales
	><li>Cupressales
	><li>Cyanidiales
	><li>Cyatheales
	><li>Cycadales
	><li>Desmidiales
	><li>Dioscoreales</li
></ol> <ol
	><li>Dipsacales</li
	><li>Ericales</li
	><li>Euglenales</li
	><li>Eupodiscales</li
	><li>Fabales</li
	><li>Fagales</li
	><li>Fragilariales</li
	><li>Fucales</li
	><li>Funariales</li
	><li>Garryales</li
	><li>Gentianales</li
	><li>Geraniales</li
	><li>Ginkgoales</li
	><li>Gracilariales</li
	><li>Hypnales</li
	><li>Lamiales</li
	><li>Laminariales</li
	><li>Laurales</li
	><li>Liliales</li
	><li>Lycopodiales</li
	><li>Magnoliales</li
	><li>Malpighiales</li
></ol> <ol
	><li>Malvales</li
	><li>Mamiellales</li
	><li>Marchantiales</li
	><li>Monomastigales</li
	><li>Myrtales</li
	><li>Nymphaeales</li
	><li>Orthotrichales</li
	><li>Pinales</li
	><li>Poales</li
	><li>Polypodiales</li
	><li>Proteales</li
	><li>Pyrenomonadales</li
	><li>Ranunculales</li
	><li>Rosales</li
	><li>Sapindales</li
	><li>Saxifragales</li
	><li>Solanales</li
	><li>Sphaeropleales</li
	><li>Takakiales</li
	><li>Ulvales</li
	><li>Vaucheriales</li
	><li>Vitales</li
	><li>Zingiberales</li
	><li>Zygnematales</li
></ol>


<h4>Output</h4>

<b>Trimmed Reads</b>:
	For paired-end data, four trimmed read files will be made. Files ending in *trimmed_P1.fq and *trimmed_P2.fq are still paired-end.  The files ending in *trimmed_U1.fq and *trimmed_U2.fq are single-end.

<b>Chloroplast Mapped Reads</b>:
	Reads that mapped to the bowtie index are in the files map_pair_hits.1.fq, map_pair_hits.2.fq, and map_hits.fq for paired-end and single-end respectively. The file Taxon-ID.sam (see below) is the standard bowtie2 mapping output but is not used.

<b>SPAdes Directory</b>:
	The directory "spades_iter1" will contain the SPAdes assembly and standard SPAdes output files.

<b>Coverage Filtered SPAdes Assembly</b>:
	The file "filtered_spades_contigs.fsa" contains contigs from the SPAdes assemblies that have low coverage, relative to the average coverage of all contigs, removed.  This is to remove mitochondrial contamination.

<b>Afin Contigs</b>:
	Output from afin is in the file "Taxon-ID_afin.fa".  Often, this will be a single contig, but if more than one is present, either the coverage of the plastome is not high enough or there may be regions of the plastome with more than a single, dominant sequence. A log file that shows the steps afin took in the extention and assembly process if found in "Taxon-ID_afin.log".

<b>Regions Files</b>:
	Files ending in *regions_split2.fsa and *regions_split3.fsa will be created that Fast-Plast's best guesses for the LSC, IR, and SSC regions of the plastome.  These use the same sequence from the afin assembly, but use a different pattern for collapsing and removed putative single copy regions.

<b>Complete Chloroplast Genome</b>:
	Full chloroplast genomes for both 2 and 3 splits are completed and named  Taxon-ID_2_FULLCP.fsa and Taxon-ID_3_FULLCP.fsa. Users can check the size of these files to see which they think is correct for their plastome. Automation for this process is in development.


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
		--bowtie_index		Order for sample to draw references for mapping. If order exists, then all available samples for that order 					will be used. If order does not exist in default set or the terms "all" or "GenBank" are given, one 							exemplar from each available order is used to build the Bowtie2 indicies. [default="All"]
		--user_bowtie		User supplied bowtie2 indices. If this option is used, bowtie_index is ignored.
		--posgenes			User defined genes for identification of single copy/IR regions and orientation. Useful when major 								rearrangments are present in user plastomes.
		--coverage_analysis Flag to run the coverage analysis of a final chloroplast assembly.[Recommended]


<h3>Coverage Analysis</h3>

The coverage analysis pipeline uses a 25 base pair sliding window to look at coverage across your assembly.  A graphical representation of the coverage is provided.

<h4>Input</h4>

The chloroplast assembly provided by Fast-Plast and the trimmed reads from the initial step of the Fast-Plast pipeline.

<h4>Output</h4>

<b>Chloroplast Mapped Reads</b>:
	Reads that mapped to the bowtie index are in the files map_pair_hits.1.fq, map_pair_hits.2.fq, and map_hits.fq for paired-end and single-end respectively. The file plastid_coverage.sam (see below) is the standard bowtie2 mapping output but is not used.

<b>Jellyfish Output</b>:
	Binary output of kmer counts from Jellyfish is found in mer_counts.jf.  The FASTA format of counts for unique kmers is found in Accession_ID_dump.

<b>Coverage for Plastome</b>:
	Sliding window coverage for the chloroplast assembly is found in the file Accession_ID.coverage_25kmer.txt. This file is tab delimited with a sequence column, a position column, and a coverage column.

<b>Graph of Coverage</b>:
	Coverage is coveraged to a graph in R and presented as a pdf.  X-axis representes position across the chloroplast genome and the y-axis is the coverage.  If a dot in the graph is red, then this is a position with 0 coverage.

<b>Poorly Assembled Regions</b>:
	The coverage file is processed to identify all stretches of 25 or more that have a coverage of 0.  Users can alter this to any length or coverage they choose.  We default at 25 based on the 25-mer sliding window and a coverage of 0 showing no representation in the reads.  This file is Accession_ID_coverage_25kmer_problem_regions_plastid_assembly.txt


## Changelog
* 07-March-2017	First Release Fast-Plast v.1.1.0

<h1>References</h1>

Bankevich, A., S. Nurk, D. Antipov, A. A. Gurevich, M. Dvorkin, A. S. Kulikov, V. M. Lesin, S. I. Niklenko, S. Pham, A. D. Prjibelski, A. V. Pyshkin, A. V. Sirotkin, N. Vyahhi, G. Tesler, M.A. Alexkseyev, and P. A. Pevzner. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. <i>J. Comp. Biol.</i>, 19(5):455-477.<br>
Bolger, A. M., M. Lohse, and B. Usadel. 2014. Trimmomatic: A flexible trimmer for Illumina Sequence Data. <i>Bioinformatics</i>, btu170.<br>
Langmead, B. and S. Salzberg. 2012. Fast gapped-read alignment with Bowtie 2. <i>Nature Methods</i>, 9:357-359.<br>
Marçais, G. and C. Kingsford. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. <i>Bioinformatics</i>, 27:764-770.<br>
