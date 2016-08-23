// afinit
//
// Coded and compiled using c++11 standard

#include "process.hpp"
#include "log.hpp"
#include <getopt.h>
#include <iostream>

// Usage function
void print_usage( std::string prog ){
  std::cout << "Usage: " << prog << " -c contigsfile(s) -r readsfile(s) [-o outfile] [-m sort_char] [-s sub_len]" << std::endl;
  std::cout << "          [-l search_loops] [-i min_cov] [-p min_overlap] [-t max_threads]" << std::endl;
  std::cout << "          [-d initial_trim] [-e max_missed] [-f stop_ext] [-g mismatch] [-x extend_len]" << std::endl;
  std::cout << "          [--silent] [--no_log] [--no_fusion] [--verbose] [--print_fused]" << std::endl << std::endl;
  std::cout << "       " << prog << " -h [--help]" << std::endl << std::endl;
  std::cout << std::endl;
  std::cout << "  -c,--contigsfiles       Space (or comma) separated list of files containing contigs" << std::endl;
  std::cout << "  -r,--readsfiles         Space (or comma) separated list of files containing reads" << std::endl;
  std::cout << "  -o,--outfile            Output will be printed to the outfile specified, with a .fa extension for the contigs and .log extension for the logfile" << std::endl;
  std::cout << "  -m,--sort_char          [default:   4] Sorts the reads by the first max_sort_char characters" << std::endl;
  std::cout << "  -s,--sub_len            [default: 100] Will focus on the current last contig_sub_len characters of the contig in each search" << std::endl;
  std::cout << "  -l,--search_loops       [default:  10] Will search against each contig a maximum of max_search_loops times to attempt extension" << std::endl;
  std::cout << "  -i,--min_cov            [default:   3] Will stop adding bp's once the coverage falls below min_cov" << std::endl;
  std::cout << "  -p,--min_overlap        [default:  20] Only those reads overlapping the contig by at least min_overlap bp's will be returned in each search" << std::endl;
  std::cout << "  -t,--max_threads        [default:   4] Will only run max_threads threads at a time" << std::endl;
  std::cout << "  -d,--initial_trim       [default:   0] Length to trim off the beginning and end of each contig at the start of the program" << std::endl;
  std::cout << "  -e,--max_missed         [default:   5] Maximum allowable mismatched bp's for each read when checking troubled contig fusions" << std::endl;
  std::cout << "  -f,--stop_ext           [default:  .5] During extension, if the percentage of reads remaining after cleaning is below stop_ext, do not extend here" << std::endl;
  std::cout << "  -g,--mismatch           [default:  .1] maximum percentage of mismatches allowed when fusing two contigs" << std::endl;
  std::cout << "  -x,--extend_len         [default:  40] Will add a max of extend_len bp's each search loop" << std::endl;
  std::cout << "  --silent                Suppress screen output" << std::endl;
  std::cout << "  --no_log                Suppress log file creation" << std::endl;
  std::cout << "  --no_fusion             Only extend, no attempt will be made to fuse contigs" << std::endl;
  std::cout << "  --verbose               Output additional information to logfile and/or screen (except if output to that location is suppressed)" << std::endl;
  std::cout << "  --print_fused           Print to file (_fused.fasta) fused contigs just before fusion, for inspecting the fusion locations" << std::endl << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////\
// BEGIN MAIN FUNCTION ///////////////////////////////////////////////////////////>
/////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv ){
  //////////////////////////////////
  // Process Command Line Options //
  //////////////////////////////////
  max_threads = 4;
  test_run = false;
  print_fused = 0;
  screen_output = 1;
  log_output = 1;
  verbose = 0;
  no_fusion = 0;

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
    {"silent",        no_argument,  &screen_output, 0},
    {"no_log",        no_argument,  &log_output,  0},
    {"no_fusion",     no_argument,  &no_fusion,  1},
    {"verbose",       no_argument,  &verbose,  1},
    {"print_fused",   no_argument,  &print_fused, 1},
    {"contigsfiles",  required_argument,  0,  'c'},
    {"readsfiles",    required_argument,  0,  'r'},
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
    {"stop_ext",      required_argument,  0,  'f'},
    {"help",          no_argument,        0,  'h'},
    {"test_run",      no_argument,        0,  'z'},
    {0, 0, 0, 0}
  };

  // get all options that have been provided on the command line
  while (( c = getopt_long(argc, argv, "hr:c:o:s:l:x:m:i:p:t:d:e:f:g:z", long_options, &option_index )) != -1 ) {
    switch( c ) {
      case 0:
        /* If this option set a flag, do nothing else now. */
        if (long_options[option_index].flag != 0)
          break;
        printf ("option %s", long_options[option_index].name);
        if (optarg)
          printf (" with arg %s", optarg);
        printf ("\n");
        break;
      case 'h':
        print_usage( argv[0] );
        exit(0);
        break;
      // max_sort_char option
      case 'm':
        process.iterable_opts["max_sort_char"] = optarg;
        break;
      // contig_sub_len option
      case 's':
        process.iterable_opts["contig_sub_len"] = optarg;
        break;
      // extend_len option
      case 'x':
        process.iterable_opts["extend_len"] = optarg;
        break;
      // max_search_loops option
      case 'l':
        process.iterable_opts["max_search_loops"] = optarg;
        break;
      // min_cov option
      case 'i':
        process.iterable_opts["min_cov"] = optarg;
        break;
      // min_overlap option
      case 'p':
        process.iterable_opts["min_overlap"] = optarg;
        break;
      // max_threads option
      case 't':
        max_threads = std::stoi(optarg);
        break;
      // initial_trim option
      case 'd':
        initial_trim = std::stoi(optarg);
        break;
      // max_missed option
      case 'e':
        process.iterable_opts["max_missed"] = optarg;
        break;
      // stop_ext option
      case 'f':
        process.iterable_opts["stop_ext"] = optarg;
        break;
      // mismatch_threshold option
      case 'g':
        process.iterable_opts["mismatch_threshold"] = optarg;
        break;
      // outputfile option
      case 'o':
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
    fprintf( stderr, "%s: Error: contigsfile(s) must be provided.\n", argv[0] );
    quit_flag = true;
  }

  if( process.readsfiles == "" ){
    fprintf( stderr, "%s: Error: readsfile(s) must be provided.\n", argv[0] );
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

  Log::Inst()->close_log();
  return 0;
}

//////////////
// END MAIN //
//////////////
