Fast-Plast: Rapid de novo assembly and finishing for whole chloroplast genomes
=============
<b>Authors</b>: Michael R. McKain, <a href="https://github.com/afinit/afin">Mark Wilson</a><br>
</br>
Version 0.3 beta, August 12, 2016<br>
</br>
<b>Contact</b>: https://github.com/mrmckain
<h3>Description</h3>

Fast-Plast is a pipeline that leverages existing and novel programs to quickly assemble, orient, and verify whole chloroplast genome sequences. For most datasets with sufficient data, Fast-Plast is able to produce a full-length de novo chloroplast genome assembly in approximately 30 minutes with no user mediation. 

Currently, Fast-Plast is written to accomodate Illumina data, though most data types could be used with a few changes.

Fast-Plast uses a de novo assembly approach by combining the de bruijn graph-based method of SPAdes with an iterative seed-based assembly implemented in afin to close gaps of contigs with low coverage. The pipeline then identifies regions from the quadripartite structure of the chloroplast genome, assigns identity, and orders them according to standard convention. A coverage analysis is then conducted to assess the quality of the final assembly. 


<h4>Requirements</h4>

Fast-Plast requires <a href="http://www.usadellab.org/cms/?page=trimmomatic">Trimmomatic</a>, <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a>, <a href="http://bioinf.spbau.ru/spades">SPAdes</a>, and <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download">BLAST+</a>.

If you use the coverage analysis to verify the assembly, then <a href="http://www.genome.umd.edu/jellyfish.html#Release">Jellyfish 2</a> and R will be needed. We highly recommend the coverage analysis to check the Fast-Plast assembly. 

afin requires a c++ complier with c++11 support and zlib.h.  zlib.h is a standard base library for most *nix systems but can be obtained <a href="http://www.zlib.net/">here</a>.

Fast-Plast is coded to use 4 threads during the Trimmomatic, bowtie2, SPAdes, and afin steps. This can simply be changed by the user if this number is not available.

Memory requirements will vary based on the size of your data set. Expect to use 1.5-2x the memory for the size of your reads files. If your data set is exceptionally large, we have found success in reducing the dataset to 50 million reads and running them through Fast-Plast.

<h4>Installation</h4>

Download or clone this repository. afin needs to be compiled and paths to the various required components should be set.

We are developing a different control file that will allow users to change parameters throughout the pipeline and set the paths more easily.  Until then, follow these instructions.

<h5>afin</h5>
<ul>
	<li><code> cd afin </code></li>
	<li><code> make </code></li>
</ul>

<h5>Control File</h5>
The control file is run_fast-plast.sh and needs paths to the programs Trimmomatic, bowtie2, SPAdes, and the Fast-Plast repository.

The Trimmomatic command includes adapter trimming and needs a file.  We use the file included in the bin directory for adapters from the NEBNext DNA Ultra II library prep kit for Illumina. If this will not work for you data, it should be changed.

The bowtie2 command uses the bowtie index "Verdant" packaged with Fast-Plast.  This is a collection 320 whole chloroplast genomes from GenBank and the Verdant chloroplast database. These are all angiosperm chloroplast genomes. If your taxa are not angiosperms, we suggest using a data set that is phylogenetically closer to your samples.

<h5>Plastome Finisher</h5>
The plastome_finisher.sh scripts in the bin directory needs to be adjusted to add paths to blastn and the Fast-Plast repository.

The blastn search uses the "positional genes" database present in the Fast-Plast bin directory.  These include the genes psbA (assumed large single copy), rpl23 (assumed IR), and ndhF (assumed small single copy) all expected to be on the reverse strand.  If this does not follow with your samples, you can add a blast nucleotide database that includes genes from the LSC, IR, and SSC. The script orientate_plastome.pl expects the gene names psbA, rpl23, and ndhF to identify regions of the plastome. If you are not using these, then the names need to be changed in the script (or you can use these names in your blast database).  <b>We will update this for flexibility in the near future.</b>

<h4>Input</h4>

Input files are in FASTQ format and expected to be paired-end. The data is not expected to be adapter trimmed or quality filtered, though this will not impede the assembly. If you data is not paired-end, contact Michael through this repository and I will provide you with a verision that handles single-end data.  We will have this as an option in an update soon.

Fast-Plast was built for genome survey sequence (aka genome skimming or low-pass genome sequencing) data. Sequence capture data can be used, but need to be normalized first.  If your data is from a sequence capture experiment (aka target enrichments or anchored phylogenomics), contact Michael and I will be able to provide you with scripts for normalization.  We currently use the normalization method packaged with <a href="https://github.com/trinityrnaseq/trinityrnaseq/wiki">Trinity</a> but have not implemented this in the current version of Fast-Plast.

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


<h4>Usage</h4>

<code> run_fast-plast.sh fastq-file1 fastq-file2 Taxon-ID </code>

Definition:

          --fastq-file1        Fastq file for paired end read 1
          --fastq-file2        Fastq file for paired end read 2
          --Taxon-ID           Name for your sample


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


<h4>Usage</h4>

<code>Coverage_Analysis/plastome_coverage_paired-single.sh Chloroplast_Assembly Accession_ID Paired_End-1 Paired_End-2 Single_End</code>

Definition:

          --Chloroplast_Assembly        Assembly in fasta format
          --Accession_ID                Name for your sample
          --Paired_End-1                Fastq file for paired-end read 1
          --Paired_End-2                Fastq file for paired-end read 2
          --Single_End                  Comma delimited list of single-end fastq files of reads


The current version needs Jellyfish (v2+) and bowtie2, expected to be in PATH.  We usually use Trimmomatic trimmed paired end reads, so we use the two trimmed paired end files and both unpaired files that trimmomatic outputs.

<h3>References</h3>

Bankevich, A., S. Nurk, D. Antipov, A. A. Gurevich, M. Dvorkin, A. S. Kulikov, V. M. Lesin, S. I. Niklenko, S. Pham, A. D. Prjibelski, A. V. Pyshkin, A. V. Sirotkin, N. Vyahhi, G. Tesler, M.A. Alexkseyev, and P. A. Pevzner. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. <i>J. Comp. Biol.</i>, 19(5):455-477.<br>
Bolger, A. M., M. Lohse, and B. Usadel. 2014. Trimmomatic: A flexible trimmer for Illumina Sequence Data. <i>Bioinformatics</i>, btu170.<br>
Langmead, B. and S. Salzberg. 2012. Fast gapped-read alignment with Bowtie 2. <i>Nature Methods</i>, 9:357-359.<br>
Marçais, G. and C. Kingsford. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. <i>Bioinformatics</i>, 27:764-770.<br>
