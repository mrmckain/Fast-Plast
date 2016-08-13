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

TreeIO from BioPerl is required to use PUG.

<h4>Input</h4>