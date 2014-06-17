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
#include <stdexcept>
#include "read.hpp"
#include "print_time.hpp"
#include "process.hpp"
#include "contig.hpp"
#include "revcomp.hpp"

using namespace std;

// return the read range corresponding to the first max_sort_char characters of the read_seg passed
tuple<long,long,long,long> get_read_range( string read_seg ){
  try{
    return read_range.at(read_seg);
  }
  catch( const out_of_range& e ){
    return make_tuple( -1,-1,-1,-1 );
  }
}

////////// Contig FUNCTIONS ////////////
Contig::Contig( string str, int cov ) : min_cov( cov ), contig( str ){
}

Contig::Contig( string str ) : contig( str ){
  min_cov = min_cov_init;
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

// determines where the read passed matches the contig if at all for off the front matches
void Contig::match_contig_fr(){
  // loop through possible substrings of contig_sub to check for matches in readlist
  for( int i=contig.length()-2; i>=min_overlap; i-- ){
    string contig_sub_rc( contig.substr( 0, i ));
    contig_sub_rc = revcomp( contig_sub_rc );
    tuple<long,long,long,long> range = get_read_range( contig_sub_rc.substr( 0, max_sort_char ) );

    if( get<0>(range) != -1 ){
      // check reads in range against contig 
      if( get<0>(range) != 0 ){
        for( int j=get<0>(range)-1; j<get<1>(range); j++ ){
          // check if the current read matches the contig from this point
          if( readlist[j].compare( 0, contig_sub_rc.length(), contig_sub_rc ) == 0 ){
            string read_rc = revcomp( readlist[j] );
            push_match( read_rc, i-contig.length(), true );
          }
        }
      }

      // check if there are revcomp reads to be checked
      if( get<2>(range) != 0 ){
        // check reverse complements
        for( int j=get<2>(range)-1; j<get<3>(range); j++ ){
          // check if the current read matches the contig from this point
          string read_ = readlist[rc_reflist[j]];
          string rc = revcomp( read_ );

          if( rc.compare( 0, contig_sub_rc.length() - 1, contig_sub_rc ) == 0 ){
            push_match( read_, i-contig.length() );
          }
        }
      }
    }
  }
}

// determines where the read passed matches the contig if at all for off the back matches
void Contig::match_contig_rr(){
  // loop through possible substrings of contig_sub to check for matches in readlist
  for( int i=1; i<contig.length()-min_overlap; i++ ){
    string contig_sub( contig.substr( i, contig.length()-1 ));
    tuple<long,long,long,long> range = get_read_range( contig_sub.substr( 0, max_sort_char ) );

    if( get<0>(range) != -1 ){
      // check reads in range against contig 
      if( get<0>(range) != 0 ){
        for( int j=get<0>(range)-1; j<get<1>(range); j++ ){
          // check if the current read matches the contig from this point
          if( readlist[j].compare( 0, contig_sub.length(), contig_sub ) == 0 ){
            push_match( readlist[j], i );
          }
        }
      }

      // check if there are revcomp reads to be checked
      if( get<2>(range) != 0 ){
        // check reverse complements
        for( int j=get<2>(range)-1; j<get<3>(range); j++ ){
          // check if the current read matches the contig from this point
          string rc = readlist[rc_reflist[j]];
          rc = revcomp( rc );

          if( rc.compare( 0, contig_sub.length() - 1, contig_sub ) == 0 ){
            push_match( rc, i, true );
          }
        }
      }
    }
  }
}

/// checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
/// used for off the front matching
string Contig::check_match_fr( int len ){

  matchlist.clear();
  match_contig_fr();
  int start = -1;

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
    if( check_cov( start - j ) < min_cov ){
      break;
    }

    
    for( i=0; i<matchlist.size(); i++ ){
      int next_char = matchlist[i].getPos( start-j );

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
    extension.insert( 0, ATCGstr.substr( max, 1 ) );
  }
 
  return extension;
}

/// checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
/// used for off the back matching
string Contig::check_match_rr( int len ){

  matchlist.clear();
  match_contig_rr();
  int start = contig.length();

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

// extend performs loops iterations of check_match with length extend_len of each extension, at each iteration the extension is added to contig, and uses contig_sub_len characters from the front or back of the contig, which end is determined by the boolean value of back
void Contig::extend( bool back ){
  string extension("");
  string contig_sub_str("");
  
  // use the end of the contig len characters long and create new contig object to process and get extension from
  for( int i=0; i<max_search_loops; i++ ){
    printf( "extend loop#: %2.d  time: ", i );
    print_time();
    if( back ){
      contig_sub_str = contig.substr( contig.length() - ( contig_sub_len ) );
    }
    else{
      contig_sub_str = contig.substr( 0, contig_sub_len );
    }
      
    Contig contig_sub( contig_sub_str, min_cov_init );
 
    // get extension through check_match
    if( back ){
      extension = contig_sub.check_match_rr( extend_len );
    }
    else{
      extension = contig_sub.check_match_fr( extend_len );
    }

    cout << "extension:" << extension << endl;
    if( extension.length() == 0 ){
      break;
    }

    if( back ){
      contig.append( extension );
    }
    else{
      contig.insert( 0, extension );
    }
    printf( "extend loop#: %2.d  time: ", i );
    print_time();
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
