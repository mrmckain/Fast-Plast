// benine
//
// Coded and compiled using c++11 standard

#include "process.hpp"
#include "log.hpp"
#include <getopt.h>
#include <iostream>

using namespace std;

// TASK:: include class for handling fastq or fasta.. return each read/contig with next call or something like it
// TASK:: consider storing things such as the more permanent and larger objects on the heap
// TASK:: remove unnecessary include statements
// TASK:: Print info about extension section?
// TASK:: Print fused_fasta file on option
// TASK:: change method of finding 2x coverage region to being more active.. use the read matching percentages.. when there's a 50/50 split, don't extenda
// TASK:: make sure that 2x cov contig ends retain this characterization when fused to other contigs
// TASK:: create installer with ability to test for presence of zlib? and/or install zlib
//          give installer capability to install to default directory or accept input directory to install executable to.. or just leave it in the base directory of the code
// TASK:: add signal handling
// TASK:: Review what goes into the logfile vs what gets printed to the screen
// TASK:: Change names of any functions that no longer title their function
// TASK:: Write up documentation explaining each option, its purpose, and why the default is set the way it is

// Usage function
void print_usage( string prog ){
  cout << "Usage: " << prog << " -c contigfile(s) -r readfile(s) [-o outfile] [-m max_sort_char] [-s contig_sub_len]" << endl;
  cout << "          [-l max_search_loops] [-i min_cov] [-p min_overlap] [-t max_threads]" << endl;
  cout << "       " << prog << " -h" << endl;
  cout << endl;
  cout << "  -c contigfile          Comma separated list of files containing contigs" << endl;
  cout << "  -r readfile            Comma separated list of files containing reads" << endl;
  cout << "  -o outfile             Output will be printed to the outfile specified with a .fasta extension" << endl;
  cout << "  -m sort_char           [default:   4] Sorts the reads by the first max_sort_char characters" << endl;
  cout << "  -s sub_len             [default: 100] Will focus on the current last contig_sub_len characters of the contig in each search" << endl;
  cout << "  -l search_loops        [default:  10] Will search against each contig a maximum of max_search_loops times before comparing them" << endl;
  cout << "  -i min_cov             [default:   3] Will stop adding bp's once the coverage falls below min_cov" << endl;
  cout << "  -p min_overlap         [default:  20] Only those reads overlapping the contig by at least min_overlap bp's will be returned in each search" << endl;
  cout << "  -t max_threads         [default:   4] Will only run max_threads threads at a time" << endl;
  cout << "  -d initial_trim        [default:  20] Length to trim off the beginning and end of each contig at the start of the program" << endl;
  cout << "  -e max_missed          [default:   5] Maximum allowable mismatched bp's for each read" << endl;
  cout << "  -g mismatch_threshold  [default:  .1] maximum percentage of mismatches allowed when fusing two contigs" << endl;
  cout << "  -x extend_len          [default:  40] Will add a max of 80 bp's each search loop" << endl << endl;
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
  extend_len = 40;
  max_sort_char = 4;
  min_cov = 3;
  min_overlap = 20;
  max_threads = 4;
  initial_trim = 20;
  max_missed = 5;
  mismatch_threshold = 0.1;
  test_run = false;
  print_fused = 0;
  screen_output = 1;
  log_output = 1;
  verbose = 0;

  int c;
  bool quit_flag = false;
  Process process;
  
  // initialize global Log object
  Log::Inst();
  
  // prevent output to stderr if erroneous option is found
  opterr = 0;
  
  int option_index = 0;

  static struct option long_options[] =
  {
    {"silent",  no_argument,  &screen_output, 0},
    {"no_log",  no_argument,  &log_output,  0},
    {"verbose",  no_argument,  &verbose,  1},
    {"print_fused", no_argument,  &print_fused, 1},
    {"help",  no_argument,  0,  'h'},
    {"contigfile",    required_argument,  0,  'c'},
    {"readfile",      required_argument,  0,  'r'},
    {"outfile",       required_argument,  0,  'o'},
    {"sort_char",     required_argument,  0,  'm'},
    {"sub_len",       required_argument,  0,  's'},
    {"search_loops",  required_argument,  0,  'l'},
    {"min_cov",       required_argument,  0,  'i'},
    {"min_overlap",   required_argument,  0,  'p'},
    {"max_threads",   required_argument,  0,  't'},
    {"initial_trim",  required_argument,  0,  'd'},
    {"max_missed",    required_argument,  0,  'e'},
    {"mismatch",      required_argument,  0,  'g'},
    {"extend_len",    required_argument,  0,  'x'},
    {"test_run",      no_argument,  0,  'z'},
    {0, 0, 0, 0}
  };
  
  // get all options that have been provided on the command line
  while (( c = getopt_long(argc, argv, "hr:c:o:s:l:x:m:i:p:t:d:e:g:z", long_options, &option_index )) != -1 ) {
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
      // min_cov option
      case 'i':
        min_cov = atoi(optarg);
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
        Log::Inst()->log_it( string("output file: ") + optarg );
        process.outfile = optarg;
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

  // start run
  process.start_run();

  return 0;
}

//////////////
// END MAIN //
//////////////
