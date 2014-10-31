# README #

After long nights and careful consideration I've come to the conclusion that this project needs some documentation. 

## Purpose ##
afin is a program designed to supplement assembly of plastomes through seed based microassembly.

### Usage ###


  Usage: ./afin -c contigsfile(s) -r readsfile(s) [-o outfile] [-m sort_char] [-s sub_len]
            [-l search_loops] [-i min_cov] [-p min_overlap] [-t max_threads]
            [-d initial_trim] [-e max_missed] [-g mismatch] [-x extend_len]
            [--silent] [--no_log] [--verbose] [--print_fused]

         ./afin -h [--help]

    -c, contigsfiles       Space (or comma) separated list of files containing contigs
    -r, readsfiles         Space (or comma) separated list of files containing reads
    -o, outfile            Output will be printed to the outfile specified with a .fasta extension
    -m, sort_char          [default:   4] Sorts the reads by the first max_sort_char characters
    -s, sub_len            [default: 100] Will focus on the current last contig_sub_len characters of the contig in each search                                                 
    -l, search_loops       [default:  10] Will search against each contig a maximum of max_search_loops times before comparing them
    -i, min_cov            [default:   3] Will stop adding bp's once the coverage falls below min_cov
    -p, min_overlap        [default:  20] Only those reads overlapping the contig by at least min_overlap bp's will be returned in each search                                  
    -t, max_threads        [default:   4] Will only run max_threads threads at a time
    -d, initial_trim       [default:   0] Length to trim off the beginning and end of each contig at the start of the program
    -e, max_missed         [default:   5] Maximum allowable mismatched bp's for each read
    -g, mismatch           [default:  .1] maximum percentage of mismatches allowed when fusing two contigs
    -x, extend_len         [default:  40] Will add a max of extend_len bp's each search loop
    --silent               Suppress screen output
    --no_log               Suppress log file creation
    --verbose              Output additional information to logfile and/or screen (except if output to that location is suppressed)
    --print_fused          Print to file (_fused.fasta) fused contigs just before fusion, for inspecting the fusion locations




### Installation ###

* Download the code base by cloning this git repository
* `cd ./afin`
* `make`

To ensure execution from any directory:
* `PATH=$PATH:/full/path/to/afin/`
-or-
* `ln -s /full/path/to/afin/afin ~/${dir_of_executables}`

### Questions? Comments? Concerns?  ###

* For any of these concerning afin, contact at benine.3121 at gmail
* Else, please contact authorities in the appropriate specialty
