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

Fast-Plast requires <a href="http://www.usadellab.org/cms/?page=trimmomatic">Trimmomatic</a>, <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a>, and <a href="http://bioinf.spbau.ru/spades">SPAdes</a>.

If you use the coverage analysis to verify the assembly, then <a href="http://www.genome.umd.edu/jellyfish.html#Release">Jellyfish 2</a> will be needed. We highly recommend the coverage analysis to check the Fast-Plast assembly. 

Fast-Plast is coded to use 4 threads during the Trimmomatic, bowtie2, SPAdes, and afin steps. This can simply be changed by the user if this number is not available.



<h4>Input</h4>

<h4>References</h4>

Bankevich, A., S. Nurk, D. Antipov, A. A. Gurevich, M. Dvorkin, A. S. Kulikov, V. M. Lesin, S. I. Niklenko, S. Pham, A. D. Prjibelski, A. V. Pyshkin, A. V. Sirotkin, N. Vyahhi, G. Tesler, M.A. Alexkseyev, and P. A. Pevzner. 2012. SPAdes: a new genome assembly algorithm and its applications to single-cell sequencing. <i>J. Comp. Biol.</i>, 19(5):455-477.<br>
Bolger, A. M., M. Lohse, and B. Usadel. 2014. Trimmomatic: A flexible trimmer for Illumina Sequence Data. <i>Bioinformatics</i>, btu170.<br>
Langmead, B. and S. Salzberg. 2012. Fast gapped-read alignment with Bowtie 2. <i>Nature Methods</i>, 9:357-359.<br>
Marçais, G. and C. Kingsford. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. <i>Bioinformatics</i>, 27:764-770.<br>