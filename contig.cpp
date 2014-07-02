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
Contig::Contig( string str, string id, int cov, int init_added_fr, int init_added_rr ) : min_cov( cov ), contig( str ), contig_id(id), bp_added_fr(init_added_fr), bp_added_rr(init_added_rr){}

Contig::Contig( string str, string id, int cov, int bp_added_init ) : min_cov( cov ), contig( str ), contig_id(id){
  bp_added_fr = bp_added_init;
  bp_added_rr = bp_added_init;
}

Contig::Contig( string str, string id, int cov ) : min_cov( cov ), contig( str ), contig_id(id){
  bp_added_fr = 0;
  bp_added_rr = 0;
}

Contig::Contig( string str, string id ) : contig( str ), contig_id(id){
  min_cov = min_cov_init;
  bp_added_fr = 0;
  bp_added_rr = 0;
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

// set contig_id
void Contig::set_contig_id( string new_contig_id ){
  contig_id = new_contig_id;
}

// return contig_id
string Contig::get_contig_id(){
  return contig_id;
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

// returns bp_added_fr
int Contig::get_bp_added_fr(){
  return bp_added_fr;
}

// returns bp_added_rr
int Contig::get_bp_added_rr(){
  return bp_added_rr;
}

// resets bp_added variables to 0
int Contig::reset_bp_added(){
  bp_added_rr = 0;
  bp_added_fr = 0;
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

// checks each matched read against the contig one bp at a time and filters out poorly aligning reads
string Contig::check_match( ){

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
string Contig::create_extension( int len, bool back ){
  // contains multiplier for position calculation
  int pos_mult = -1;
  int start = -1;

  // create reference string for numeric based additions to the extension string
  string ATCGstr( "ATCG" );
  string extension( "" );
  matchlist.clear();

  // vector of int vectors to hold values of nucleotides at each position, the 5th member of the in vector will be the max number for that position
  vector< vector< int > > ATCG;
  
  // get matches
  if( back ){
    match_contig_rr();
    start = contig.length();
    pos_mult = 1;
  }
  else{
    match_contig_fr();
  }

  // return if no matches are found
  if( matchlist.size() == 0 ){
    return "";
  }

  // create missed bp's vector to keep track of how many bp's each read contains that are below the max at that position
  vector<int> missed_bp( matchlist.size(), 0 );
  int missed_bp_tot = 0;
  int missed_bp_avg = 0;
 
  /////////////////////////////////////
  // STEP 1: First Pass Over Matches //
  // loop through bp's for initial pass to filter out errant matches
  for( int j=0; j<len; j++ ){
    // initialize temporary vector 
    vector<int> ATCG_curr = { 0,0,0,0,0 };
    
    // loop through matches to get count of each nucleotide present at the current position
    for( int i=0; i<matchlist.size(); i++ ){
      int next_char = matchlist[i].getPos( start+(pos_mult*j) );

      if ( next_char == 'A' ) {
        ATCG_curr[0]++;
      }
      else if ( next_char == 'T' ) {
        ATCG_curr[1]++;
      }
      else if ( next_char == 'C' ) {
        ATCG_curr[2]++;
      }
      else if ( next_char == 'G' ) {
        ATCG_curr[3]++;
      }
    }

    // determine maximum number of any nucleotide 
    for( int i=0; i<4; i++ ){
      if( ATCG_curr[i] > ATCG_curr[4] ){
        ATCG_curr[4] = ATCG_curr[i];
      }
    }

    // check coverage level and break if below min
    if( ATCG_curr[4] < min_cov ){
      len = j;
      break;
    }

    // initialize next bp to 0 for each nucleotide
    ATCG.push_back( ATCG_curr );
  }
  
  ///////////////////////////////////////////
  // STEP 2: Mark Errant Matches For Death //
  // loop through bp's to determine which reads, if any, should be eliminated from the matchlist
  // if nucleotide of read is < max for that position, count it against the read, otherwise don't count it
  for( int j=0; j<len; j++ ){
    for( int i=0; i<matchlist.size(); i++ ){
      int next_char = matchlist[i].getPos( start+(pos_mult*j));
      switch( next_char ){
        case 'A':
          if( ATCG[j][0] < ATCG[j][4] ){
            missed_bp[i]++;
            missed_bp_tot++;
          }
          break;
        case 'T':
          if( ATCG[j][1] < ATCG[j][4] ){
            missed_bp[i]++;
            missed_bp_tot++;
          }
          break;
        case 'C':
          if( ATCG[j][2] < ATCG[j][4] ){
            missed_bp[i]++;
            missed_bp_tot++;
          }
          break;
        case 'G':
          if( ATCG[j][3] < ATCG[j][4] ){
            missed_bp[i]++;
            missed_bp_tot++;
          }
          break;
        default:
          break;
      }
    }
  }

  if( matchlist.size() != 0 ){
    // calculate avg missed_bp's
    missed_bp_avg = missed_bp_tot / matchlist.size() + 1;
  }
  else{
    return "";
  }

  ///////////////////////////////////
  // STEP 3: Remove Errant Matches //
  // loop through missed_bp list to eliminate any matches that have a missed level greater than the avg 
  for( int i=0; i<missed_bp.size(); i++ ){
    if( missed_bp[i] > missed_bp_avg ){
      matchlist.erase( matchlist.begin() + i );
      missed_bp.erase( missed_bp.begin() + i );
      i--;
    }
  }
    
  //////////////////////////////
  // STEP 4: Create Extension //
  cout << "create_extension1 start: " << start << "  contig: " << contig << endl;
  // loop len times processing 1 basepair at a time
  for( int j=0; j<len; j++ ){
    vector<int> ATCG_curr = { 0,0,0,0 };
    int max = 0;
    int avg = 0;
    
    // check coverage at current position
    if( check_cov( start+(pos_mult*j) ) < min_cov ){
      break;
    }
    
    for( int i=0; i<matchlist.size(); i++ ){
      int next_char = matchlist[i].getPos( start+(pos_mult*j) );

      if ( next_char == 'A' ) {
        ATCG_curr[0]++;
      }
      else if ( next_char == 'T' ) {
        ATCG_curr[1]++;
      }
      else if ( next_char == 'C' ) {
        ATCG_curr[2]++;
      }
      else if ( next_char == 'G' ) {
        ATCG_curr[3]++;
      }
    }

    // find character with greatest appearance
    for( int i=1; i<4; i++ ){
      if( ATCG_curr[i] > ATCG_curr[max] ) {
        max = i;
      }
    }

    // add next base
    if( back ){
      extension.append( ATCGstr.substr( max, 1 ) );
      bp_added_rr++;
    }
    else{
      extension.insert( 0, ATCGstr.substr( max, 1 ) );
      bp_added_fr++;
    }
  }
 
  return extension;
}

// extend performs loops iterations of create_extension with length extend_len of each extension, at each iteration the extension is added to contig, and uses contig_sub_len characters from the front or back of the contig, which end is determined by the boolean value of back
void Contig::extend( bool back ){
  string extension("");
  string contig_sub_str("");
 
  // skip over any contigs that present at least double coverage
  if( contig_id.compare( 0,3,"2x_" ) == 0 ){
    return;
  }

  // get extension through create_extension
  if( back ){
    contig_sub_str = contig.substr( contig.length() - ( contig_sub_len ) );
  }
  else{
    contig_sub_str = contig.substr( 0, contig_sub_len );
  }

  Contig contig_sub( contig_sub_str, "temp", min_cov_init );
  extension = contig_sub.create_extension( extend_len, back );
  
  //cout << "extension:" << extension << endl;
  if( extension.length() == 0 ){
    return;
  }

  if( back ){
    contig.append( extension );
    bp_added_rr += contig_sub.get_bp_added_rr();
      
    // print messages to logfile about current actions
    //Process::print_to_logfile( contig_sub.get_bp_added_rr() + " basepairs were added to the back of " + contig_id ); 
  }
  else{
    contig.insert( 0, extension );
    bp_added_fr += contig_sub.get_bp_added_fr();
      
    // print messages to logfile about current actions
    //Process::print_to_logfile( contig_sub.get_bp_added_fr() + " basepairs were added to the front of " + contig_id ); 
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
