# README #

After long nights and careful consideration I've come to the conclusion that this project needs some documentation, so here it is.

## Purpose ##
`afin` is a program designed to supplement assembly of plastomes through de novo, seed based microassembly.

### Usage ###


    Usage: ./afin -c contigsfile(s) -r readsfile(s) [-o outfile] [-m sort_char] [-s sub_len]
            [-l search_loops] [-i min_cov] [-p min_overlap] [-t max_threads]
            [-d initial_trim] [-e max_missed] [-f stop_ext] [-g mismatch] [-x extend_len]
            [--silent] [--no_log] [--no_fusion] [--verbose] [--print_fused]

         ./afin -h [--help]

      -c, contigsfiles       Space (or comma) separated list of files containing contigs
      -r, readsfiles         Space (or comma) separated list of files containing reads
      -o, outfile            [default: afin_out] Output will be printed to the outfile specified, with a .fa extension for the contigs and .log extension for the logfile
      -m, sort_char          [default:   4] Sorts the reads by the first max_sort_char characters
      -s, sub_len            [default: 100] Will focus on the current last contig_sub_len characters of the contig in each search                                                 
      -l, search_loops       [default:  10] Will search against each contig a maximum of max_search_loops times before comparing them
      -i, min_cov            [default:   3] Will stop adding bases once the coverage falls below min_cov
      -p, min_overlap        [default:  20] Only those reads overlapping the contig by at least min_overlap bases will be returned in each search                                  
      -t, max_threads        [default:   4] Will only run max_threads threads at a time
      -d, initial_trim       [default:   0] Length to trim off the beginning and end of each contig at the start of the program
      -e, max_missed         [default:   5] Maximum allowable mismatched bases for each read when checking troubled contig fusions
      -f, stop_ext           [default:  .5] During extension, if the percentage of reads remaining after cleaning is below stop_ext, do not extend here
      -g, mismatch           [default:  .1] maximum percentage of mismatches allowed when fusing two contigs
      -x, extend_len         [default:  40] Will add a max of extend_len bases each search loop
      --silent               Suppress screen output
      --no_log               Suppress log file creation
      --no_fusion            Only extend, no attempt will be made to fuse contigs
      --verbose              Output additional information to logfile and/or screen (except if output to that location is suppressed)
      --print_fused          Print to file (_fused.fasta) fused contigs just before fusion, for inspecting the fusion locations

* contigsfiles  -- `/path/contigs.fa /path2/contigs.fasta,/path3/contigs.fsa` This is a list of files containing the contigs produced by an assembler of your choice. Typically this will be just one file, and testing indicated that the best input for `afin`, was output from the SPAdes assembler
* readsfiles  --  Same format as for contigsfiles. Input can be in plain text format or gzip'd
* outfile  --  File prefix for output files
* sort_char  --  Number of characters to sort and search by when organizing and referencing reads. A higher number here increases sort time and decreases extension time
* sub_len  --  Length of end of the contig to match against while extending
* search_loops  --  Iterations of extension and fusion
* min_cov  --  Minimum number of reads at a base to count during extension
* min_overlap  --  Minimum bases needed to overlap for a fusion to occur
* max_threads  --  Maximum number of threads capable of running during the extension process. The fusion process is not threaded
* initial_trim  --  Number of bases to trim off the beginning and end of contigs before processing
* max_missed  --  Rarely used. Used to set a max level for mismatched bases while attempting to check if the end of one contig is errant during the fusion process
* stop_ext  --  This number represents the fraction of reads required to remain after cleaning to continue with the extension. The purpose here is to avoid extending into regions of the genome where there are two clear paths that could be followed
* mismatch  --  Maximum percentage of mismatched bases to be allowed for the fusion of two contigs to occur
* extend_len  --  Maximum number of bases to extend a contig by during the extension process


### Installation ###

* Download the code base by cloning this git repository
* `cd ./afin`
* `make`

To ensure execution from any directory:
* `PATH=$PATH:/full/path/to/afin/`
-or-
* `ln -s /full/path/to/afin/afin ~/${dir_of_executables}`

### Requirements ###
* c++ compiler with c++11 support
* zlib.h installed (This is part of the base installation for most *nix systems )
    * If an error related to zlib is thrown during compilation, the library can be obtained [here](http://www.zlib.net/) or talk to your system administrator

afin has been tested on various linux distros and multiple Apple systems, using `gcc` and `clang` to compile.

### Overview of Approach ###
The popular approach to de novo assembly is to use de Bruijn graphs to create networks of kmers and then build contigs from strongest paths through these networks. Often this method is not able to assemble a plastome completely. By focusing on the ends of the contigs, the data can be used to draw connections between the contigs. By matching reads to the end of a contig and extending based on the results of that matching process, two things are accomplished that distinguish this method from typical de Bruijn graph assemblers. The entire information from each read is maintained, and information that may have been too uncertain when looking at the genome as a whole, can be used to support the fusion of two contigs. `afin` accomplishes this by breaking the process into two parts, Fusion and Extension, and alternating between the two.

#### Fusion ####
The Fusion process is approached by finding the best match score for fusing any two contigs. This score is calculated from overlap length and number of bases that do not match, or mismatches. If this best score falls below the mismatch_threshold, the two contigs are fused together. Any differences are resolved by dropping the end of each contig and only using inner half of the overlap from each contig.

#### Extension ####
The Extension process takes a piece of each end from each contig, then finds all reads that have an exact match to the end of this piece. The list of reads is cleaned by finding the reads most consistent with each other. The consensus bases are then added to the contig until a maximum extension length, or until the number of reads supporting the extension fall below the minimum coverage level at that base.


### Questions? Comments? Concerns?  ###

* For any of these concerning afin, contact at mark *dot* wilson *dot* 3121 *at* gmail
* Else, please contact authorities in the appropriate specialty
