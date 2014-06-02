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
#include "contig.hpp"
#include "print_time.hpp"
#include "process.hpp"
#include "read.hpp"

using namespace std;

// TASK:: split this into several header files
// TASK:: consider revcomp hash list that points to readlist positions
// TASK:: develop search based on sorted readlist
// TASK:: match up min overlap better with max extension? Extending doesn't add too much time but searching adds a great deal of time
// TASK:: create separate methods for fasta and fastq files
// TASK:: create methods for detecting fasta vs fastq 
// TASK:: Add processing for IUPAC DNA ambiguity codes
// TASK:: Add processing for differences in reads( ie, create new contig objects for differing sets of matches, add method for splitting matchlist between two new contig objects ), determine which contig is correct
// TASK:: Check to see if it would be beneficial to leave the offthefront() function out of find_part(), maybe include it when determining if two contigs connect, but probably not
// TASK:: Add threading capability
// TASK:: Put the following functions in class


/////////////////////////////////////////////////////////////////////////////////\
// BEGIN MAIN FUNCTION ///////////////////////////////////////////////////////////>
/////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv ){
  //////////////////////////////////
  // Process Command Line Options //
  //////////////////////////////////
  int c;
  Process process;
  
  // prevent output to stderr if erroneous option is found
  opterr = 0;

  // get all options that have been provided on the command line
  while (( c = getopt (argc, argv, "hr:c:" )) != -1 ) {
    switch( c ) {
      case 'h':
        cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        break;
      case 'r':
        cout << "readfile: " << optarg << endl;
        cout << "add_reads time: ";
        print_time();
        process.add_reads( optarg );
        print_time();
        break;
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


  cout << "contig: " << process.get_contig(0) << endl;

  process.contigs[0].extend( 4, 80, 100 );

  cout << "contex: " << process.get_contig(0) << endl;


  cout << "THIS HAS BEEN A PROCESS CLASS TEST. THANK YOU FOR YOUR PATIENCE!" << endl;

  //////////////
  // End Test //
  //////////////


  return 0;
}

//////////////
// END MAIN //
//////////////
