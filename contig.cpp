// $Author: benine $
// $Date$
// $Log$
// Contains the Contig class for afin

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <unordered_map>
#include <tuple>
#include "read.hpp"
#include "print_time.hpp"
#include "contig.hpp"
#include "process.hpp"
#include "revcomp.hpp"

using namespace std;

// TASK:: consider revcomp hash list that points to readlist positions.. Possibly store revcomp in same hash list? perhaps not.......
// TASK:: develop search based on sorted readlist
// TASK:: match up min overlap better with max extension? Extending doesn't add too much time but searching adds a great deal of time
// TASK:: create separate methods for fasta and fastq files
// TASK:: create methods for detecting fasta vs fastq 
// TASK:: Add processing for IUPAC DNA ambiguity codes
// TASK:: Add processing for differences in reads( ie, create new contig objects for differing sets of matches, add method for splitting matchlist between two new contig objects ), determine which contig is correct
// TASK:: Check to see if it would be beneficial to leave the offthefront() function out of find_part(), maybe include it when determining if two contigs connect, but probably not
// TASK:: Add threading capability
// TASK:: Put the following functions in class

// TASK:: Function to process alignment of matches
//    -> Find first starting place
//      -> based on number of reads that match the region searched for, only process portions where x% of coverage of matches
//    -> assess count of each bp at each position
//      ->keep track of misses in each read
//    -> possibility :::-->>> maintain multiple read lists when there is a dispute over the proper character at a position
//      -> It is possible that two branches comeback together
//      -> Rate value of each branch using weighted values 
//    -> take most frequently found character, potentially using weighted values involving how many misses in that read
//
//    -> Every x characters examined, if all looks good, add those characters to the contig and search reads again, adding to the current matchlist
//
//    -> Analyze read headings and assign weights of some sort based on coverage of each read if possible
//    -> Analyze contig headings and make decisions based on coverage of each contig
//

#ifndef MAX_SORT
#define MAX_SORT 4
#endif

////////// Contig FUNCTIONS ////////////
Contig::Contig( string str, int cov ) : min_cov( cov ), contig( str ){
  min_overlap = 20;
}

Contig::Contig( string str ) : contig( str ){
  min_overlap = 20;
  min_cov = 10;
}   

// adds a read to the read list
void Contig::push_match( string read, int pos ){
  matchlist.push_back(Read( read, pos ));
}

void Contig::push_match( string read, int pos, bool revcomp ){
  matchlist.push_back(Read( read, pos, revcomp ));
}

string Contig::getContig(){
  return contig;
}

int Contig::getListSize(){
  return matchlist.size();
}

string Contig::getRead_s( int i ){
  return matchlist[i].getRead();
}

int Contig::getStart( int i ){
  return matchlist[i].getStart();
}

Read Contig::getRead( int i ){
  return matchlist[i];
}

// clear matchlist to make room for new matches
void Contig::clear_matches(){
  matchlist.clear();
}

void Contig::check_pos( int pos ){
  // TASK:: should this be a class? return an object containing how many matches of each and where from? Prolly not
      
}

// finds the initial point at which the matches meet the min_cov requirement
int Contig::find_start(){
  long int cov;
  int next_char;
  first_read = contig.length();
  
  // Find the position of the first matching read on the contig
  for( int i=0; i<matchlist.size(); i++ ){
    if( -getRead(i).getStart() < first_read ){
      first_read = -getRead(i).getStart();
    }
  }

  // process each position checking coverage of the matching reads
  for( int j=first_read; j<contig.length(); j++ ){
    cov = 0;
    
    // get coverage at position j
    cov = check_cov( j );
    
    // if coverage is at the minimum level, return the current position
    if( cov >= min_cov ){
      return( j );
    }
  }
  
  return(contig.length() + 1);
}

// return the read range corresponding to the first MAX_SORT characters of the read_seg passed
tuple<long,long,long,long> Contig::get_read_range( string read_seg ){
  return read_range[read_seg];
}

// determines where the read passed matches the contig if at all
void Contig::match_contig(){
  // loop through possible substrings of contig_sub to check for matches in readlist
  for( int i=1; i<contig.length()-min_overlap; i++ ){
    string contig_sub( contig.substr( i, MAX_SORT ));
    tuple<long,long,long,long> range = get_read_range( contig_sub.substr( 0, MAX_SORT ) );

    // check reads in range against contig 
    if( get<0>(range) != 0 ){
      for( int j=get<0>(range)-1; j<get<1>(range); j++ ){
        // check if the current read matches the contig from this point
        if( readlist[j].compare( 0, contig_sub.length() - 1, contig_sub ) == 0 ){
  	      push_match( readlist[j], -i );
        }
      }

      // check reverse complements
      for( int j=get<2>(range)-1; j<get<3>(range); j++ ){
        // check if the current read matches the contig from this point
        string rc = readlist[rc_reflist[j]];
        rc = revcomp( rc );

        if( rc.compare( 0, contig_sub.length() - 1, contig_sub ) == 0 ){
  	      push_match( rc, -i, true );
        }
      }
    }
  }
}

/// checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
string Contig::check_match( int len ){
  int start = contig.length();

  matchlist.clear();
  match_contig();

  cout << "check_match1 start: " << start << "  contig: " << contig << endl;
  // create reference string for numeric based additions to the extension string
  string ATCGstr( "ATCG" );
  int i;
  string extension( "" );

  // loop len times processing 1 basepair at a time
  for( int j=0; j<len; j++ ){
    int ATCG[] = { 0,0,0,0 };
    int max = 0;
    int avg = 0;
    
    // check coverage at current position
    if( check_cov( start + j ) < min_cov ){
      break;
    }

    
    for( i=0; i<matchlist.size(); i++ ){
      int next_char = matchlist[i].getPos( start+j );

      if ( next_char == 'A' ) {
        ATCG[0]++;
      }
      else if ( next_char == 'T' ) {
        ATCG[1]++;
      }
      else if ( next_char == 'C' ) {
        ATCG[2]++;
      }
      else if ( next_char == 'G' ) {
        ATCG[3]++;
      }
    }

    // find character with greatest appearance
    for( i=1; i<4; i++ ){
      if( ATCG[i] > ATCG[max] ) {
        max = i;
      }
    }

    // add next base
    extension.append( ATCGstr.substr( max, 1 ) );
  }
 
  return extension;
}

// check_match( int ) with a default len provided
string Contig::check_match(){
  return check_match( 20 );
}

// extend() performs loops iterations of check_match with length len of each extension, at each iteration the extension is added to contig, and uses sublen characters from the end of the contig
void Contig::extend( int loops, int len, int sublen ){
  string extension("");
  
  // if sublen=0 use entire contig, else use the end of the contig len characters long and create new contig object to process and get extension from
  if( sublen > 0 ){
    for( int i=0; i<loops; i++ ){
      printf( "extend loop#: %2.d  time: ", i );
      print_time();
      Contig contig_sub( contig.substr( contig.length() - ( sublen + 1 ) ), 3 );
   
      // get extension through check_match
      extension = contig_sub.check_match( len );
      cout << "extension:" << extension << endl;
      if( extension.length() == 0 ){
        break;
      }

      contig.append( extension );
      printf( "extend loop#: %2.d  time: ", i );
    }
  }
  else {
    for( int i=0; i<loops; i++ ){
      extension = check_match( len );
      cout << extension << endl;
      if( extension.length() == 0 ){
        break;
      }
      
      contig.append( extension );
    }
  }
}

// checks the coverage of matches at the given positions, returns the coverage
long Contig::check_cov( long pos ){
  long cov = 0;

  for( long i=0; i<matchlist.size(); i++ ){
    if( matchlist[i].getPos( pos ) != -1 ){
      cov++;
    }
  }
  return cov;
}
