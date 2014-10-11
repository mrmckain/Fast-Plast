// $Author: benine $
//
// Coded and compiled using c++11 standard

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <functional>
#include "contig.hpp"
#include "print_time.hpp"
#include "process.hpp"
#include "read.hpp"
#include "afin_util.hpp"

using namespace std;

// TASK:: check logic on check_fusion_support()
// TASK:: test changing the create_read_range() to use the actual positions instead of +1
// TASK:: remove find_start() function from contig as it seems to be unused.. clear out other unused functions
// TASK:: change method of finding 2x coverage region to being more active.. use the read matching percentages.. when there's a 50/50 split, don't extenda
// TASK:: make sure that 2x cov contig ends retain this characterization when fused to other contigs
// TASK:: create installer with ability to test for presence of zlib? and/or install zlib
//          give installer capability to install to default directory or accept input directory to install executable to.. or just leave it in the base directory of the code
// TASK:: add signal handling
// TASK:: add support for gzipped files
// TASK:: give option to suppress output to screen
// TASK:: Review what goes into the logfile vs what gets printed to the screen
// TASK:: Look into whether N's and misses should be counted together? Or how they should relate to each other
// TASK:: Exact matching in the contig_fusion_support() may be too conservative
// TASK:: Review check_cov() and use a coverage that incorporates only the matching bp at that position
// TASK:: make global variables and options for the variables that need it
// TASK:: Change names of any functions that no longer title their function
// TASK:: Write up documentation explaining each option, its purpose, and why the default is set the way it is
// TASK:: Check limits of readlist size and other limits at all points of the program
// TASK:: Clean up functions, break long functions into smaller ones and eliminate unused functions
// TASK:: Remove using line from each file and add std:: where necessary
//
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
  max_threads = 4;
  initial_trim = 20;
  max_missed = 5;
  mismatch_threshold = 0.1;
  test_run = false;

  int c;
  bool quit_flag = false;
  Process process;
  
  // prevent output to stderr if erroneous option is found
  opterr = 0;

  // get all options that have been provided on the command line
  while (( c = getopt (argc, argv, "hr:c:o:s:l:x:m:i:p:t:a:b:d:e:f:g:z" )) != -1 ) {
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
      // initial_trim option
      case 'd':
        initial_trim = atoi(optarg);
        break;
      // max_missed option
      case 'e':
        max_missed = atoi(optarg);
        break;
      // mismatch_threshold option
      case 'g':
        mismatch_threshold = atof(optarg);
        break;
      // outputfile option
      case 'o':
        cout << "output file: " << optarg << endl;
        process.outfile = optarg;
        print_time();
        break;
      // readfile option
      case 'r':
        optind--;
        // loop through each file
        while ( optind < argc && argv[optind][0] != '-' ) { 
          if( process.readsfiles != "" ){
            process.readsfiles.append( "," );
          }
          process.readsfiles.append( argv[optind] );
          optind++;
        }   
        break;
      // contig file option
      case 'c':
        optind--;
        // loop through each file
        while ( optind < argc && argv[optind][0] != '-' ) { 
          if( process.contigsfiles != "" ){
            process.contigsfiles.append( "," );
          }
          process.contigsfiles.append( argv[optind] );
          optind++;
        }   
        break;
      // test_run option
      case 'z':
        test_run = true;
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
  
  /////////////////
  // End Options //
  /////////////////

  //////////////////
  // Input Errors //
  //////////////////
  if( process.contigsfiles == "" ){
    fprintf( stderr, "%s: Error: contigs_file(s) must be provided.\n", argv[0] );
    quit_flag = true;
  }

  if( process.readsfiles == "" ){
    fprintf( stderr, "%s: Error: reads_file(s) must be provided.\n", argv[0] );
    quit_flag = true;
  }

  if( min_overlap < max_sort_char ){
    fprintf( stderr, "%s: Error: min_overlap must be >= max_sort_char.\n", argv[0] );
    quit_flag = true;
  }

  if( quit_flag ){
    exit(1);
  }

  ////////////////
  // End Errors //
  ////////////////

  
  ///////////////////////
  // Test Thread Queue //
  ///////////////////////

  // start run
  process.start_run();

  

  process.print_to_outfile();
  
  return 0;
}

//////////////
// END MAIN //
//////////////
