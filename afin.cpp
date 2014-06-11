// $Author: benine $
// $Date: 2014/04/18 19:09:09 $
// $Log: afin.cpp,v $
// Revision 1.4  2014/04/18 19:09:09  benine
// Cleaned up some and found jesus
//
// Revision 1.3  2014/04/15 01:15:40  benine
// Process class added
// Cleaned up processing
// Next step is add functionality for testing various segments of the contig being processed at the time
//
// Revision 1.2  2014/04/09 23:31:14  benine
// Process is built but throws compile errors due to nested classes attempting to access members of parent class
// Reorganizing classes now
//
// Revision 1.1  2014/04/09 23:24:43  benine
// Initial revision
//
// Revision 1.2  2014/04/02 02:06:11  benine
// Cleaned up some of the testing bits and added some more comments
//
// Revision 1.1  2014/04/01 17:23:50  benine
// Initial revision
//
// $Header: /home/benine/code/git/bb_afin/RCS/afin.cpp,v 1.4 2014/04/18 19:09:09 benine Exp benine $
//
// Coded and compiled using c++11 standard

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <tuple>
#include <utility>
#include <unordered_map>
#include <thread>
#include "contig.hpp"
#include "print_time.hpp"
#include "process.hpp"
#include "read.hpp"
#include "queue.hpp"

using namespace std;

// TASK:: look into whether coverage should be based on number of similar bases rather than total bases.. maybe originally based on similar and when coverage drops switch and make a note
// TASK:: add output file section and command line option
// TASK:: output information about where contigs are combined and where contigs have had two or more options due to duplicate regions
// TASK:: remove unnecessary variables that have been replaced by global variables 
// TASK:: add print usage function
// TASK:: add parsing for multiple read/contig files
// TASK:: add long options
// TASK:: create separate methods for fasta and fastq files
// TASK:: create methods for detecting fasta vs fastq 
// TASK:: Add processing for IUPAC DNA ambiguity codes
// TASK:: Add processing for differences in reads( ie, create new contig objects for differing sets of matches, add method for splitting matchlist between two new contig objects ), determine which contig is correct
// TASK:: Add threading capability

// PRINT USAGE FUNCTION
void print_usage( string prog ){
  cout << "Usage: " << prog << " -c [contigfile(s)] -r [readfile(s)] [-m " << endl;
}

/////////////////////////////////////////////////////////////////////////////////\
// BEGIN MAIN FUNCTION ///////////////////////////////////////////////////////////>
/////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv ){
  //////////////////////////////////
  // Process Command Line Options //
  //////////////////////////////////
  max_search_loops = 10;
  contig_sub_len = 100;
  extend_len = 80;
  max_sort_char = 4;
  min_cov_init = 5;
  min_overlap = 20;
  max_threads = 6;
  int c;
  Process process;
  
  // prevent output to stderr if erroneous option is found
  opterr = 0;

  // get all options that have been provided on the command line
  while (( c = getopt (argc, argv, "hr:c:s:l:x:m:i:p:t:" )) != -1 ) {
    switch( c ) {
      case 'h':
        cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        break;
      // max_sort_char option
      case 'm':
        max_sort_char = atoi(optarg);
        break;
      // contig_sub_len option
      case 's':
        contig_sub_len = atoi(optarg);
        break;
      // extend_len option
      case 'x':
        extend_len = atoi(optarg);
        break;
      // max_sort_char option
      case 'l':
        max_search_loops = atoi(optarg);
        break;
      // min_cov_init option
      case 'i':
        min_cov_init = atoi(optarg);
        break;
      // min_overlap option
      case 'p':
        min_overlap = atoi(optarg);
        break;
      // max_threads option
      case 't':
        max_threads = atoi(optarg);
      // readfile option
      case 'r':
        cout << "readfile: " << optarg << endl;
        cout << "add_reads time: ";
        print_time();
        process.add_reads( optarg );
        print_time();
        break;
      // contig file option
      case 'c':
        cout << "contigfile: " << optarg << endl;
        cout << "add_contigs time: ";
        print_time();
        process.add_contigs( optarg );
        print_time();
        break;
      case '?':
        if ( optopt == 'r' ){
          fprintf( stderr, "%s: Error: Option -r requires an argument. ", argv[0] );
          cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        }
        else if ( optopt == 'c' ){
          fprintf( stderr, "%s: Error: Option -c requires an argument. ", argv[0] );
          cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        }
        else if ( isprint( optopt )){
          fprintf( stderr, "%s: Error: Unknown option -%c. \n", argv[0], optopt );
          cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        }
        else{
          fprintf( stderr, "%s: Error: Unknown option character %x.\n", argv[0], optopt );
          cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        }
        return 1;
      default:
        abort();
    }
  }
  
  // output starting option values
  cout << "OPTION VALUES" << endl;
  cout << "contig_sub_len: " << contig_sub_len << endl;
  cout << "extend_len: " << extend_len << endl;
  cout << "max_search_loops: " << max_search_loops << endl;
  cout << "max_sort_char: " << max_sort_char << endl;
  cout << "min_cov_init: " << min_cov_init << endl;
  cout << "min_overlap: " << min_overlap << endl;
  cout << "max_threads: " << max_threads << endl;

  while ( optind < argc ) {
    cout << argv[optind] << endl;
    optind++;
  }
  /////////////////
  // End Options //
  /////////////////



  /////////////////////////
  // Test Process Class ///
  /////////////////////////

  cout << "THIS IS A PROCESS CLASS TEST! THIS IS ONLY A TEST!" << endl;

  cout << "sort_reads start: ";
  print_time();
  process.sort_reads();
  cout << "sort_rc start: ";
  print_time();
  process.sort_rc();
  cout << "create_reads_range start: ";
  print_time();
  process.create_read_range();

  cout << "extend start: ";
  print_time();

  cout << "contig: " << process.get_contig(0) << endl;

  process.contigs[0].extend();

  cout << "contex: " << process.get_contig(0) << endl;

  cout << "THIS HAS BEEN A PROCESS CLASS TEST. THANK YOU FOR YOUR PATIENCE!" << endl;


  //////////////
  // End Test //
  //////////////
  
  ///////////////////////
  // Test Thread Queue //
  ///////////////////////

  cout << "This is testing thread methods with a queue and a max thread count" << endl;

  // create thread array with max_thread entries
  thread t[max_threads];





  return 0;
}

//////////////
// END MAIN //
//////////////
