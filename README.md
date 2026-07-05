Fast-Plast: Rapid de novo assembly and finishing for whole chloroplast genomes
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.973887.svg)](https://doi.org/10.5281/zenodo.973887)
=============
<b>Authors</b>: Michael R. McKain, <a href="https://github.com/afinit/afin">Mark Wilson</a><br>
</br>
Version 1.3.0<br>
</br>
<b>Contact</b>: https://github.com/mrmckain

<img src="https://github.com/mrmckain/Fast-Plast/blob/master/extras/Fast-Plast_Logo.png" width="256" alt="fast-plast_logo" align="middle">

<h1>Description</h1>

Fast-Plast is a pipeline that leverages existing and novel programs to quickly assemble, orient, and verify whole chloroplast genome sequences. For most datasets with sufficient data, Fast-Plast is able to produce a full-length de novo chloroplast genome assembly in less than 20 minutes with no user mediation. In addition to a chloroplast sequence, Fast-Plast identifies chloroplast genes present in the final assembly.

Currently, Fast-Plast is written to accomodate Illumina data, although most data types could be used.

Fast-Plast uses a de novo assembly approach by combining the De Bruijn graph-based method of SPAdes with an iterative seed-based assembly implemented in afin to close gaps of contigs with low coverage. When more than one contig remains, the pipeline uses RagTag to order and orient the contigs against a closely related reference plastome and then re-extends with afin to close the joins. The pipeline then identifies regions from the quadripartite structure of the chloroplast genome, assigns identity, and orders them according to standard convention. A coverage analysis is then conducted to assess the quality of the final assembly.

<h1>What's new in 1.3.0</h1>

Version 1.3.0 modernizes the pipeline's dependencies and installation and moves the reference database out of the repository. In brief:

* Read trimming now uses <b>fastp</b> in place of Trimmomatic.
* Gzipped reads are streamed through <b>pigz</b> (falling back to gzip), which is faster and handles multi-member gzip transparently.
* Scaffolding now uses <b>RagTag</b> (reference-guided ordering/orienting) in place of SSPACE, removing the SSPACE and Bowtie1 dependencies. afin then re-extends the ordered scaffold to close gaps.
* Tools are resolved from your <b>PATH</b> at runtime (with optional per-tool environment overrides), so Fast-Plast installs cleanly from conda and no longer bakes absolute paths into the driver script.
* The bundled GenBank plastome reference set is now hosted on <b>Zenodo</b> and fetched on demand, rather than shipped in the repository. Reference headers are bare accessions, with taxonomy stored in a companion metadata table.

See the <a href="#changelog">Changelog</a> for details.

<h1>Dependencies</h1>

Fast-Plast requires a number of commonly used bioinformatics programs. The simplest way to obtain them is via conda (see Installation); they can also be installed by hand as long as each is available on your `PATH`.

* Perl 5.x which has been installed with threading enabled (see <a href="https://perlmaven.com/how-to-build-perl-from-source-code">How to build perl from source</a>, if needed).
* <a href="https://github.com/OpenGene/fastp">fastp</a> (initial read cleaning)
* <a href="https://zlib.net/pigz/">pigz</a> (parallel decompression of gzipped reads; gzip is used as a fallback)
* <a href="http://bowtie-bio.sourceforge.net/bowtie2/index.shtml">Bowtie2</a> (read reduction to only chloroplast-like reads)
* <a href="http://bioinf.spbau.ru/spades">SPAdes</a> (initial assembly)
* <a href="https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download">BLAST+</a> (multiple checks for gene content)
* <a href="https://github.com/malonge/RagTag">RagTag</a> (reference-guided scaffolding when a single contig is not obtainable)
* <a href="https://github.com/lh3/minimap2">minimap2</a> (required by RagTag)
* <a href="https://www.r-project.org/">R</a> (required to plot coverage)

<h3>afin</h3>

* c++ complier with c++11 support and zlib.h.  zlib.h is a standard base library for most Unix systems but can be obtained <a href="http://www.zlib.net/">here</a>.

<h3>Coverage Analysis</h3>

* <a href="http://www.genome.umd.edu/jellyfish.html#Release">Jellyfish 2</a><br>
* R

Memory requirements will vary based on the size of your data set. Expect to use 1.5-2x the memory for the size of your reads files. If your data set is exceptionally large, we have found success in reducing the dataset to 5-10 million reads and running them through Fast-Plast.

<h1>Installation</h1>

Fast-Plast has been tested on Linux but should be compatible with any flavor that can handle dependencies.

Clone the Github repository:

    git clone https://github.com/mrmckain/Fast-Plast.git

<h3>Dependencies via conda (recommended)</h3>

All external dependencies are available through conda (bioconda). From the cloned repository, an environment file is provided:

    conda env create -f environment.yml
    conda activate fast-plast

This installs fastp, pigz, Bowtie2, SPAdes, BLAST+, RagTag, minimap2, Jellyfish 2, and R. With the environment active, every tool is on your `PATH` and Fast-Plast will find it automatically.

<h3>Compiling afin</h3>

afin is the compiled C++ core and must be built once:

    cd afin
    make
    cd ..

Instructions for direct compilation of afin can also be found <a href="https://github.com/mrmckain/Fast-Plast/tree/master/afin">here</a>.

<h3>Reference database</h3>

The GenBank plastome reference set is no longer shipped in the repository (it is too large). It is hosted on Zenodo and installed with the provided fetch script:

    bin/fetch_plastome_db.sh

This downloads the database archive, verifies it against a checksum, and installs `GenBank_Plastomes` and `GenBank_Plastomes.metadata.tsv` into `bin/`, where Fast-Plast expects them. Run it once after cloning. See <a href="#reference-database">Reference Database</a> below for details.

<h3>Tool resolution and overrides</h3>

Fast-Plast locates each external tool on your `PATH` at runtime. If a tool is installed somewhere unusual, or you want to pin a specific build, you can override any tool without editing the script by exporting an environment variable, for example:

    export FP_FASTP=/opt/fastp/bin/fastp
    export FP_BLAST=/opt/blast/bin/blastn
    export FP_SPADES=/opt/spades/bin/spades.py

Corresponding variables exist for the other tools (`FP_BOWTIE2`, `FP_JELLYFISH`, `FP_PIGZ`, `FP_AFIN`, and so on). The location of the reference data can be overridden with `FASTPLAST_SHARE`.

<h1>Input</h1>

<h3>Reads</h3>
Input files are in FASTQ format. The data are not expected to be adapter trimmed or quality filtered, though this will not impede the assembly.

Fast-Plast was built for genome survey sequence (aka genome skimming or low-pass genome sequencing) data. Sequence capture data can be used but needs to be normalized first.  If your data are from a sequence capture experiment (aka target enrichments or anchored phylogenomics), we suggest using the normalization method packaged with <a href="https://github.com/trinityrnaseq/trinityrnaseq/wiki">Trinity</a>, <a href="http://ged.msu.edu/papers/2012-diginorm/">khmer</a>, or <a href="https://sourceforge.net/projects/bbmap/">bbnorm</a>.

<h3>Bowtie Index</h3>

Fast-Plast is packaged (via the Zenodo download above) with a large set of whole chloroplast genomes from GenBank, spanning a wide range of diversity including marine algae, angiosperms, ferns, and more. To choose references for read mapping, use the `--bowtie_index` option with a taxon of interest; Fast-Plast will pull all matching members and build a Bowtie2 index. If your taxon is not present, or you pass "All" or "GenBank", Fast-Plast uses one representative sequence from each order present. "All" is the default.

Example:

    --bowtie_index Poales

This will use all available Poales plastomes in the data set.

Any of the following taxonomic levels is accepted for `--bowtie_index`: genus, species epithet, tribe (if applicable), subfamily (if applicable), family, and order. Multiple taxa may be given as a comma-separated list. Matching is case-insensitive.

The taxa and accessions currently available are enumerated in the companion metadata table, `bin/GenBank_Plastomes.metadata.tsv` (installed by the fetch script). This table is the authoritative, versioned record of database contents; consult it (or the Zenodo record) for the current set of orders, families, and genera rather than a static list.

<h3>User Provided Bowtie Index</h3>

A user-made bowtie index can be provided using the --user_bowtie option. The full path to the index should be given. If this option is used, the --bowtie_index option will be ignored.

<h3>Name</h3>

The --name option should be used with each run. This simply gives a prefix to all files. The default is "Fast-Plast".

<h1 id="reference-database">Reference Database</h1>

The reference plastome database is distributed separately from the code, on Zenodo, and consists of two files:

* `GenBank_Plastomes` — a FASTA of whole chloroplast genomes whose headers are bare, unique accessions (e.g. `>NC_000932.1`). Short accession headers are compact and avoid BLAST identifier-length limits.
* `GenBank_Plastomes.metadata.tsv` — a tab-separated table mapping each accession to its taxonomy (`accession`, `genus`, `species`, `tribe`, `subfamily`, `family`, `order`, `source`). Ranks that are genuinely not assigned are recorded as `NA`.

Fast-Plast reads taxonomy from this metadata table for both per-order representative sampling and `--bowtie_index` taxon selection. (For backward compatibility, if the metadata table is absent, Fast-Plast falls back to reading taxonomy from legacy rich-format headers.)

<h3>Installing / updating the database</h3>

    bin/fetch_plastome_db.sh

fetches the current release, verifies its checksum, and installs both files into `bin/`. To pin or change the release, edit the version, URL, and checksum at the top of the script.

<h3>Maintainer tools</h3>

Two scripts in `bin/` are used to build a new database release and are not needed for normal runs:

* `normalize_plastome_db.pl` — converts a legacy rich-header FASTA into the accession-header FASTA plus the metadata table.
* `repair_taxonomy.py` — re-fetches taxonomy for accessions with missing ranks directly from NCBI (accession → TaxId → lineage) using Biopython. Requires network access and an NCBI-registered email.

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
For paired-end data, three trimmed read files will be made. Files ending in *trimmed_P1.fq and *trimmed_P2.fq are still paired-end.  The files ending in *trimmed_UP.fq is single-end. If only single end files are used, then only the *trimmed_UP.fq file will be found. Per-library fastp reports (JSON/HTML/log) are also written here.

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
If more than one contig is found in the afin assembly and paired-end reads were used, RagTag is invoked to order and orient the contigs against a reference plastome, and afin then re-extends the ordered scaffold to close the joins. These results are found in the `Scaffolding/` subdirectory, including the chosen reference (`scaffold_reference.fsa`) and the RagTag output (`ragtag_out/`). By default the reference is selected automatically as the best BLAST match among the bundled plastomes; a specific reference can be supplied with `--scaffold_reference` (see Usage). If more than one contig/scaffold remains in the final output, this is the last step of the pipeline and results are found in Final_Assembly.

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

	perl fast-plast.pl -1 sample_R1_001.fastq.gz -2 sample_R2_001.fastq.gz --name Urochloa_fusca-37 --bowtie_index Poales --coverage_analysis --clean light

In this example, one pair-end library is being used for assembly. The default adapters (NEB) are used for trimming, Poales species are used for the Bowtie2 index, the coverage analysis is invoked, and a light cleaning is done after completition.

<h3>Example with Single End Data</h3>

	perl fast-plast.pl --single sample_R1_001.fastq.gz --name Chionachne_koenigii-TK057 --adapters TruSeq --bowtie_index All --coverage_analysis

In this example, one single end library is being used for assembly. The TruSeq adapters are used for trimming, a representative from each order is used for the Bowtie2 index, the coverage analysis is invoked, and no cleaning is done.

<h3>Example with Mixed Libraries</h3>

	perl fast-plast.pl -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz --single sample_SE_R1_001.fastq --name Monocymbium_ceresiiforme-TK203 --user_bowtie /path/to/androplast --coverage_analysis --clean deep

In this example, one single end library and one paired-end library are being used for assembly. The default adapters (NEB) are used for trimming, a user-defined Bowtie2 index (base name given) is used for the Bowtie2 index, the coverage analysis is invoked, and a deep cleaning is done after completion.

<h3>Example with Coverage Analysis Only</h3>

    perl fast-plast.pl -1 sample_R1_001.fastq.gz -2 sample_R2_001.fastq.gz --name Schizachyrium_scoparium-TK686R --only_coverage /path/to/Schizachyrium_scoparium-TK686R_FULLCP.fsa --min_coverage 2 &

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

		--min_length_trim	Acceptable minimum length for reads after trimming for adapters and quality. (Default = 140]
		--subsample		Number of reads to subsample. Reads will be evenly pulled from all files. 
		--threads		Number of threads used by Fast-Plast.  [Default = 4]
		--min_coverage		Lowest acceptable coverage for 25-mer sliding window during coverage analysis. [Default = 0.25 * Average coverage]
  		--min_filter_spades	Minimum coverage allowed for SPAdes contig to be passed to afin. Only recommended to change if default is not working. [Default is one standard deviation of weigthed average length of contigs.]
		--adapters		[NEB|Nextera|TruSeq] Adapters used in making the sequencing library. NEB, Nextera, and TruSeq options 
					available. Also accepts the path to a user created FASTA file of adapters. [Default = NEB]
		--bowtie_index		Taxon used to pick references for the bowtie2 index. If the taxon is in the database, all matching samples are used. 
					If it is not present, or the terms "all" or "GenBank" are given, one exemplar from each available order is used. 
					Multiple taxa may be given as a comma-separated list. Accepted levels: genus, species epithet, tribe (if applicable), 
					subfamily (if applicable), family, and order. [default="All"]
		--user_bowtie		User supplied bowtie2 indices. If this option is used, bowtie_index is ignored.
		--scaffold_reference	Path to a reference plastome (FASTA) to use for RagTag scaffolding, overriding automatic reference selection. 
					Useful when a close relative is known. May also be set with the FP_REFERENCE environment variable.
		--coverage_analysis 	Flag to run the coverage analysis of a final chloroplast assembly. [Recommended]
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


<h1 id="changelog">Changelog</h1>

* 2026 Fast-Plast v.1.3.0 <br>
    --Read trimming switched from Trimmomatic to fastp; `--min_length_trim` now sets fastp's minimum length. <br>
    --Gzipped reads streamed through pigz (fallback gzip); multi-member gzip handled correctly. <br>
    --Scaffolding switched from SSPACE to RagTag (reference-guided ordering/orienting), followed by afin re-extension to close joins. SSPACE and Bowtie1 dependencies removed. <br>
    --Added `--scaffold_reference` (and `FP_REFERENCE`) to supply a scaffolding reference directly. <br>
    --External tools are now resolved from `PATH` at runtime, with optional `FP_*` environment overrides, in place of paths baked in at install. Conda/bioconda installation supported (`environment.yml`). <br>
    --Reference plastome database moved to Zenodo and fetched with `bin/fetch_plastome_db.sh`. Database headers are now bare accessions with taxonomy in a companion metadata table (`GenBank_Plastomes.metadata.tsv`); `--bowtie_index` selection reads taxonomy from this table (with legacy header fallback). Maintainer tools `normalize_plastome_db.pl` and `repair_taxonomy.py` added. <br>

* 14-June_2022 Fast-Plast v.1.2.9 <br>
    --Minimum length of trim reads option added to pass to Trimmomatic.
   <br>
    --Afin reduced to one iteration. Afin now uses mapped reads to complete plastome instead of all reads.
   <br>
	
* 28-June_2018 Fast-Plast v.1.2.8 <br>
    --Output directory of BLAST database for gene identification changed to user working directory.
    <br>

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
* Chen, S., Y. Zhou, Y. Chen, and J. Gu. 2018. fastp: an ultra-fast all-in-one FASTQ preprocessor. <i>Bioinformatics</i>, 34(17):i884-i890.
* Langmead, B. and S. Salzberg. 2012. Fast gapped-read alignment with Bowtie 2. <i>Nature Methods</i>, 9:357-359.
* Marçais, G. and C. Kingsford. 2011. A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. <i>Bioinformatics</i>, 27:764-770.
* Alonge, M., L. Lebeigle, M. Kirsche, K. Jenike, S. Ou, S. Aganezov, X. Wang, Z. Lippman, M. Schatz, and S. Soyk. 2022. Automated assembly scaffolding using RagTag elevates a new tomato system for high-throughput genome editing. <i>Genome Biology</i>, 23:258.
* Li, H. 2018. Minimap2: pairwise alignment for nucleotide sequences. <i>Bioinformatics</i>, 34(18):3094-3100.
