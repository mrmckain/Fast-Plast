// $Author: benine $
// $Date$
// $Log$
// Contains the Process class for afin

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <unordered_map>
#include <tuple>
#include "print_time.hpp"
#include "revcomp.hpp"
#include "process.hpp"
#include "contig.hpp"

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

#ifndef MAX_SORT
#define MAX_SORT 4
#endif

////////////////////////////////////
//////// PROCESS DEFINITIONS ///////
////////////////////////////////////
vector<string> readlist;
vector<long int> rc_reflist;
unordered_map<string, tuple<long,long,long,long>> read_range;

// uses a bubblesort to sort the read list based on the first MAX_SORT characters of each read
void Process::sort_reads(){
  for( long int i=1; i<readlist.size(); i++ ){
    for( long int j=readlist.size()-1; j>=i; j-- ){
      if( readlist[j].compare( 0, MAX_SORT, readlist[j-1] ) < 0 ){
        string temp = readlist[j];
        readlist[j] = readlist[j-1];
        readlist[j-1] = temp;
      }
    }
  }
}

// Produces list of references to the readlist sorted based on the first MAX_SORT characters of the reverse_compliment of the
//  referenced read
//  Must be done after the readlist is sorted as it contains the locations in the readlist of the referenced read
//  Uses an insertion sort to build the list
void Process::sort_rc(){
  vector<long int>::iterator iter;
  iter = rc_reflist.begin();

  rc_reflist.insert( iter, 0 );
  for( long int i=1; i<readlist.size(); i++ ){
    for( long int j=0; i<rc_reflist.size(); j++ ){
      string temp = readlist[i].substr( readlist[i].length() - MAX_SORT, MAX_SORT );
      if( readlist[rc_reflist[j]].compare( 0, MAX_SORT, revcomp( temp ) ) < 0 ){
        rc_reflist.insert( iter+j, i );
        break;
      }
    }
  }
}

// creates a hash table within the process object that contains the ranges corresponding to equivalent first MAX_SORT characters in the reads
// This increases the efficiency of searching
void Process::create_read_range(){
  string current( readlist[0].substr( 0, MAX_SORT ) );
  long curr_start = 0;
  // cycle through the readlist
  for( long int i=1; i<readlist.size(); i++ ){
    if( readlist[i].compare( 0, MAX_SORT, current ) != 0 ){
      // insert values into hash table, the ordered pair reflecting the range is increased by 1 to differentiate between a false search  
      read_range.insert( { current, make_tuple( curr_start+1, i, 0, 0 )});
      curr_start = i;
      current = readlist[i].substr( 0, MAX_SORT );
    }
  }
  
  // add the last entry
  read_range.insert( { current, make_tuple( curr_start+1, readlist.size(), 0, 0 )});

  string rc( readlist[0].substr( readlist[0].length() - MAX_SORT, MAX_SORT ) );
  current = revcomp( rc );
  curr_start = 0;

  // cycle through the rc_reflist to include these ranges in the hash as well
  for( long int i=1; i<rc_reflist.size(); i++ ){
    rc = readlist[i].substr( readlist[i].length() - MAX_SORT, MAX_SORT );
    rc = revcomp( rc );
    if( rc.compare( current ) != 0 ){
      // check if current already exists in hash
      if( get<0>( read_range[current] ) == 0 ){
        read_range.insert( { current, make_tuple( 0, 0, curr_start+1, i )});
      }
      else{
        get<2>( read_range[current] ) = curr_start+1;
        get<3>( read_range[current] ) = i;
      }
      curr_start = i;
      current = rc;
    }
  }
  
  // add the last entry
  if( get<0>( read_range[current] ) == 0 ){
    read_range.insert( { current, make_tuple( 0, 0, curr_start+1, rc_reflist.size() )});
  }
  else{
    get<2>( read_range[current] ) = curr_start+1;
    get<3>( read_range[current] ) = rc_reflist.size();
  }
}

// put reads from readfile into readlist
void Process::add_reads( string filename ){
  string buffer("");
  string line("");
  int line_count = 0;
  
  // open read file
  ifstream read( filename );
  
  // read in fastq reads
  while( getline( read, line )){
    line_count++;
    if( line[0] == '@' ){
      if( getline( read, line )){
        line_count++;
        readlist.push_back( line );
        if( getline( read, line )){
          line_count++;
          if( line[0] == '+'){
            if( getline( read, line )){
              line_count++;
              continue;
            }
            else{
              fprintf( stderr, "Error reading fastq file. Line missing. Line: %d\n", line_count );
            }
          }
          else{
            fprintf( stderr, "Error reading fastq file. '+' expected at this line. Line: %d\n", line_count );
          }
        }
        else{
          fprintf( stderr, "Error reading fastq file. Line missing. Line: %d\n", line_count );
        }
      }
      else{
        fprintf( stderr, "Error reading fastq file. Line missing. Line: %d\n", line_count );
      }
    }
    else{  
      fprintf( stderr, "Error reading fastq file. '@' expected at the beginning of this line. Line: %d\n", line_count );
    }
  }


  // close read file
  read.close();

  // insert last line into readlist
  if( buffer.length() != 0 ){
    readlist.push_back( buffer );
  }
}


// put contigs from contfile into contlist
void Process::add_contigs( string filename ){
  string buffer("");
  string line("");
  
  // open contig file
  ifstream cont( filename );

  // read in contig objects
  while( getline( cont, line ) ){
    if( line[0] == '>' && buffer.length() != 0 ){
      //cout << buffer << endl;
      contigs.push_back( Contig( buffer, 3 ));
      buffer = "";
    }
    else if ( line[0] == '>' ) {
    }
    else{
      buffer += line;
    }
  }
  
  // insert last line into contigs list
  if( buffer.length() != 0 ){
    contigs.push_back( Contig( buffer, 3 ) );
  }
  
  // close contig file
  cont.close();
}

// return contig with index contig_ind
string Process::get_contig( int contig_ind ){
  return contigs[ contig_ind ].getContig();
}

//////////////////////////////
// END DEFINITIONS ///////////
//////////////////////////////
