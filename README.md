Fast-Plast: Rapid de novo assembly and finishing for whole chloroplast genomes
=============
<b>Authors</b>: Michael R. McKain, <a href="https://bitbucket.org/afinit/afin">Mark Wilson</a><br>
</br>
Version 1.0, August 12, 2016<br>
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
	<li><code> cd afinit-afin-c25d552842f5 </code></li>
	<li><code> make </code></li>
</ul>

<h5>Control File</h5>
The control file is run_fast-plast.sh and needs paths to the programs Trimmomatic, bowtie2, SPAdes, and the Fast-Plast repository.

The Trimmomatic command includes adapter trimming and needs a file.  We use the file included in the bin directory for adapters from the NEBNext DNA Ultra II library prep kit for Illumina. If this will not work for you data, it should be changed.

The bowtie2 command uses the bowtie index "Verdant" packaged with Fast-Plast.  This is a collection 320 whole chloroplast genomes from GenBank and the Verdant chloroplast database. These are all angiosperm chloroplast genomes. If you taxa are not angiosperms, we suggest using a data set that is phylogenetically closer to your samples.

<h5>Plastome Finisher</h5>
The plastome_finisher.sh scripts in the bin directory needs to be adjusted to add paths to blastn and the Fast-Plast repository.

<h4>Input</h4>

<h4>References</h4>

Bankevich, A., S. Nurk, D. Antipov, A. A. Gurevich, M. Dvorkin, A. S. Kulikov, V. M. Lesin, S. I. Niklenko, S. Pham, A. D. Prjibelski, A. V. Pyshkin, A. V. Sirotkin, N. Vyahhi, G. Tesler, M.A. Alexkseyev, and P. A. Pevzner. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. <i>J. Comp. Biol.</i>, 19(5):455-477.<br>
Bolger, A. M., M. Lohse, and B. Usadel. 2014. Trimmomatic: A flexible trimmer for Illumina Sequence Data. <i>Bioinformatics</i>, btu170.<br>
Langmead, B. and S. Salzberg. 2012. Fast gapped-read alignment with Bowtie 2. <i>Nature Methods</i>, 9:357-359.<br>
Marçais, G. and C. Kingsford. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. <i>Bioinformatics</i>, 27:764-770.<br>