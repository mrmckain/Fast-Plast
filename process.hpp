// $Author: benine $
// $Date$
// $Log$
// Contains the Process class for afin

#ifndef PROCESS_H
#define PROCESS_H

#include <unordered_map>
#include "contig.hpp"
#include "mismatch.hpp"

using namespace std;

extern vector<string> readlist;
extern vector<long int> rc_reflist;

// contains hash of ranges in the readlist corresponding to keys made up of the first max_sort_char characters in the read
extern unordered_map<string, tuple<long,long,long,long>> read_range;
extern int max_search_loops;
extern int contig_sub_len;
extern int extend_len;
extern int max_sort_char;
extern int min_cov_init;
extern int min_overlap;
extern int max_threads;
extern int trim_length;
extern int tip_length;
extern int end_depth;
extern int tip_depth;
extern int initial_trim;
extern int max_missed;
extern int bp_added_init;
extern int mismatch_threshold;

/////////////////////////////////////////////\
// Process Class: ////////////////////////////>
/////////////////////////////////////////////
//
// create object here containing both Read and Contig objects
// Process object, contains all classes, methods, data, and data references necessary for processing the contigs
// There will be only one Process object needed per iteration of this program
class Process{
  private:
    string logfile;
    time_t timer;
    fstream log_fs;

  public:
    vector<Contig> contigs;
    vector<Contig> contigs_fused;
    string outfile;
    string readsfiles;
    string contigsfiles;
    
    Process();

    // sorts the reads by the first max_sort_char characters
    void sort_reads();

    // Produces list of references to the readlist sorted based on the first max_sort_char characters of the reverse_compliment of the
    //  referenced read
    //  Must be done after the readlist is sorted as it contains the locations in the readlist of the referenced read
    //  Uses an insertion sort to build the list
    void sort_rc();

    // creates a hash table within the process object that contains the ranges corresponding to equivalent first max_sort_char characters in the reads
    // This increases the efficiency of searching
    void create_read_range();

    // parses the cov value from the contig_id and passes the result back as a double
    double parse_cov( string contig_id );

    // cycles through each contig and parses out the first section of the id
    void parse_ids();

    // check coverage of each contig, calculate the average coverage, then remove into a separate data structure any contigs that have more than 2xAvg coverage
    void contig_cov();

    // put reads from readfile into readlist
    void add_reads();

    // put contigs from contfile into contlist
    void add_contigs();

    // return contig from contigs_fused with index contig_ind
    string get_contig_fused( int contig_ind );
    
    // return contig with index contig_ind
    string get_contig( int contig_ind );

    // prints results to fasta file with outfile prefix and additional information is printed to a text based file with outfile prefix
    void print_to_outfile();

    // return a string of the current time in the program
    string get_time();

    // initalize logfile
    void logfile_init();

    // prints notes to file as the program progresses
    void print_to_logfile( string note );

    // creates id of fused contigs
    string get_fused_id( string contig1_id, string contig2_id );

    // complete contig_fusion process
    void contig_fusion_log( Mismatch fusion );

    // complete contig_fusion process
    void commit_fusion( string fused, string fused_id, int index_i, int index_j, int bp_added_fr, int bp_added_rr );

    // tally mismatches in substrings passed and return score in the form of misatches per length
    double mismatch_score( string contig_sub1, string contig_sub2 );

    // check overlap section for mismatches
    Mismatch overlap_check( string contig_a, string contig_b, int overlap, int end_i, int end_j );

    // sort the match_list for easier 
    vector<Mismatch> sort_matches( vector<Mismatch> match_list );

    // cleans match_list from conflicting matches
    void clean_matches( vector<Mismatch> &match_list );

    // create fused contig string
    string build_fusion_string( string contig_a, string contig_b, int overlap );

    // returns overlap length based on the indexes passed
    int get_overlap( int i, int j, int orientation );

    // remove duplicates from contig remove list
    void dedup_list( vector<int> &list );

    // remove fused contigs from contigs list
    void process_removals( vector<int> remove_list );

    // fuse contigs wherever possible
    void contig_fusion();

    // Initializes data structures and turns over control to run_manager()
    void start_run();
    
    // Manages run 
    void run_manager();

    // closes logfile
    void close_log();
};

///////////////////////
// End Process Class //
///////////////////////

#endif
