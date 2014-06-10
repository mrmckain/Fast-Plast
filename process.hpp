// $Author: benine $
// $Date$
// $Log$
// Contains the Process class for afin

#ifndef PROCESS_H
#define PROCESS_H

#include <unordered_map>
#include "contig.hpp"

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

/////////////////////////////////////////////\
// Process Class: ////////////////////////////>
/////////////////////////////////////////////
//
// create object here containing both Read and Contig objects
// Process object, contains all classes, methods, data, and data references necessary for processing the contigs
// There will be only one Process object needed per iteration of this program
class Process{
  private:
  public:
    vector<Contig> contigs;
    
    Process(){}

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

    // put reads from readfile into readlist
    void add_reads( string filename );

    // put contigs from contfile into contlist
    void add_contigs( string filename );

    // return contig with index contig_ind
    string get_contig( int contig_ind );
};

///////////////////////
// End Process Class //
///////////////////////

#endif
