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
#include <functional>
#include "contig.hpp"
#include "print_time.hpp"
#include "process.hpp"
#include "read.hpp"
#include "afin_util.hpp"
#include "queue.tcc"

using namespace std;

// TASK:: Check limits of readlist size and other limits at all points of the program
// TASK:: Discuss calculation for cov_avg and number of times to add contig to contigs_2x.. Determine the best way to know how many times a contig is duplicated within the genome.. Or should this be genome specific? Not portable this way
// TASK:: Clean up functions, break long functions into smaller ones and eliminate unused functions
// TASK:: Remove add_reads() and add_contigs() from the options section and place them in the Process constructor
// TASK:: expand Process::print_to_outfile() to include contigs_fused vector
// TASK:: develop Process::print_to_logfile( string )
// TASK:: Remove using line from each file and add std:: where necessary
// TASK:: break Contig::create_extension() into multiple files 
// TASK:: test multiple input files
//
// TASK:: Throws out_of_range error when min_overlap < 4.... ? Not that it should ever be that low
//
// TASK:: align contigs? Remove mismatched bp's at the end of contigs?
// TASK:: add contigs together
// TASK:: make considerations for splits in possibility (eg IR boundaries, RPL23 gene copy)
//
// TASK:: output information about where contigs are combined and where contigs have had two or more options due to duplicate regions
// TASK:: add long options
// TASK:: Add processing for IUPAC DNA ambiguity codes
// TASK:: Add processing for differences in reads( ie, create new contig objects for differing sets of matches, add method for splitting matchlist between two new contig objects ), determine which contig is correct

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
  Queue<int> qu;
  
  // prevent output to stderr if erroneous option is found
  opterr = 0;

  // get all options that have been provided on the command line
  while (( c = getopt (argc, argv, "hr:c:o:s:l:x:m:i:p:t:" )) != -1 ) {
    switch( c ) {
      case 'h':
        print_usage( argv[0] );
        exit(0);
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
        break;
      // outputfile option
      case 'o':
        cout << "output file: " << optarg << endl;
        process.outfile = optarg;
        print_time();
        break;
      // readfile option
      case 'r':
        process.add_reads( optarg );
        break;
      // contig file option
      case 'c':
        process.add_contigs( optarg );
        break;
      case '?':
        if ( optopt == 'r' ){
          fprintf( stderr, "%s: Error: Option -r requires an argument. ", argv[0] );
          print_usage( argv[0] );
        }
        else if ( optopt == 'o' ){
          fprintf( stderr, "%s: Error: Option -o requires an argument. ", argv[0] );
          print_usage( argv[0] );
        }
        else if ( optopt == 'c' ){
          fprintf( stderr, "%s: Error: Option -c requires an argument. ", argv[0] );
          print_usage( argv[0] );
        }
        else if ( isprint( optopt )){
          fprintf( stderr, "%s: Error: Unknown option -%c. \n", argv[0], optopt );
          print_usage( argv[0] );
        }
        else{
          fprintf( stderr, "%s: Error: Unknown option character %x.\n", argv[0], optopt );
          print_usage( argv[0] );
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

  
  ///////////////////////
  // Test Thread Queue //
  ///////////////////////

  int length[process.contigs.size()];
  for( int i=0; i<process.contigs.size(); i++ ){
    length[i] = process.get_contig(i).length();
  }
  
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

  cout << "This is testing thread methods with a queue and a max thread count" << endl;

  // create thread array with max_thread entries
  thread t[max_threads];

  // initialize threads
  for( int i=0; i<max_threads; i++ ){
    t[i] = thread( thread_worker, ref(process.contigs), ref(qu), i );
  }

  // push each thread onto queue
  for( int i=0; i<process.contigs.size(); i++ ){
    qu.push( i );
  }

  // push stop signals onto queue for each thread
  for( int i=0; i<max_threads; i++ ){
    qu.push( -1 );
  }

  // join threads
  for( int i=0; i<max_threads; i++ ){
    t[i].join();
  }

  // print out each contig and the number of bp added
  for( int i=0; i<process.contigs.size(); i++ ){
    cout << "Contig[" << i << "]: " << process.get_contig(i) << endl;
    cout << "\tbp added: " << process.get_contig(i).length() - length[i] << endl;
  }
  
  cout << "exit time: ";
  print_time();

  process.print_to_outfile();
  
  return 0;
}

//////////////
// END MAIN //
//////////////
