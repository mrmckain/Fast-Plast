Fast-Plast: Rapid de novo assembly and finishing for whole chloroplast genomes
=============
<b>Authors</b>: Michael R. McKain<br>
</br>
				<a href="https://bitbucket.org/afinit/afin">Mark Wilson</a><br>
</br>
Version 1.0, August 12, 2016<br>
</br>
<br></br>
<b>Contact</b>: https://github.com/mrmckain
<h3>Description</h3>

PUG is designed to test the phylogenetic placment of a polyploid event over mulitple gene trees using a known species tree and a set of putative paralogs.
Despite the name, PUG works with genomic or transcriptomic data. 

Paralogs are queried in gene trees. Their placement, and the relative placement of taxa in the gene tree, are compared to the given species tree and if they are found to match a node of the species tree, the pair is counted towards supporting the polyploid event origin at that node. A single node may be represented in a gene tree by multiple paralog pairs, but that node is counted only once for the gene tree, if the "unique" option is used.  The "all pairs" option will all multiple mappings to the same node.

<h4>Requirements</h4>

TreeIO from BioPerl is required to use PUG.

<h4>Input</h4>