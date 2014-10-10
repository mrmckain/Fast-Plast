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
#include "afin_util.hpp"
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
Contig::Contig( string str, string id, double cov, int min_cov ) : cov(cov), min_cov(min_cov), contig(str), contig_id(id){
  doub_cov = false;
}

Contig::Contig( string str, string id ) : contig(str), contig_id(id){
  cov = 0;
  min_cov = min_cov_init;
  doub_cov = false;
}   

// adds a read to the read list
void Contig::push_match( string read, int pos ){
  matchlist.push_back(Read( read, pos ));
}

void Contig::push_match( string read, int pos, bool revcomp ){
  matchlist.push_back(Read( read, pos, revcomp ));
}

string Contig::get_contig(){
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

// return cov
double Contig::get_cov(){
  return cov;
}

// clear matchlist to make room for new matches
void Contig::clear_matches(){
  matchlist.clear();
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

// returns value associated with index of ATCG[] or -1 if not a member of ATCG
int Contig::get_ATCG_value( int ATCG_char ){
  if ( ATCG_char == 'A' ) {
    return 0;
  }
  else if ( ATCG_char == 'T' ) {
    return 1;
  }
  else if ( ATCG_char == 'C' ) {
    return 2;
  }
  else if ( ATCG_char == 'G' ) {
    return 3;
  }

  return -1;
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
          if( readlist[j].compare( 0, contig_sub_rc.length(), contig_sub_rc ) ){
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
          string rc = readlist[rc_reflist[j]];
          rc = revcomp( rc );

          if( rc.compare( 0, contig_sub_rc.length(), contig_sub_rc ) == 0 ){
            push_match( rc, i-contig.length() );
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
          if( readlist[j].compare( 0, contig_sub.length(), contig_sub ) ){
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

          if( rc.compare( 0, contig_sub.length(), contig_sub ) == 0 ){
            push_match( rc, i, true );
          }
        }
      }
    }
  }
}

// first step of create_extension: determine bp count and max represented bp at each position
void Contig::extension_bp_count( vector<vector<int>> &ATCG, int start, int pos_mult, int &len ){
  for( int i=0; i<len; i++ ){
    // initialize temporary vector 
    vector<int> ATCG_curr = { 0,0,0,0,0 };
    
    // loop through matches to get count of each nucleotide present at the current position
    for( int j=0; j<matchlist.size(); j++ ){
      int next_char = matchlist[j].getPos( start+(pos_mult*i) );

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
    for( int j=0; j<4; j++ ){
      if( ATCG_curr[j] > ATCG_curr[4] ){
        ATCG_curr[4] = ATCG_curr[j];
      }
    }

    // check coverage level and break if below min
    if( ATCG_curr[4] < min_cov ){
      len = i;
      break;
    }

    // initialize next bp to 0 for each nucleotide
    ATCG.push_back( ATCG_curr );
  }
}

// second step of the create_extension process: count missed bp's per read, or in other words, the bp's represented below the max for that position
void Contig::extension_missed_count( vector<vector<int>> &ATCG, vector<int> &missed_bp, int &missed_bp_tot, int start, int pos_mult, int len ){
  for( int i=0; i<len; i++ ){
    for( int j=0; j<matchlist.size(); j++ ){
      int next_char = matchlist[j].getPos( start+(pos_mult*i));
      switch( next_char ){
        case 'A':
          if( ATCG[i][0] < ATCG[i][4] ){
            missed_bp[j]++;
            missed_bp_tot++;
          }
          break;
        case 'T':
          if( ATCG[i][1] < ATCG[i][4] ){
            missed_bp[j]++;
            missed_bp_tot++;
          }
          break;
        case 'C':
          if( ATCG[i][2] < ATCG[i][4] ){
            missed_bp[j]++;
            missed_bp_tot++;
          }
          break;
        case 'G':
          if( ATCG[i][3] < ATCG[i][4] ){
            missed_bp[j]++;
            missed_bp_tot++;
          }
          break;
        default:
          break;
      }
    }
  }
}

// third step in create_extension(): removal of reads that have errors over the threshold
void Contig::extension_error_removal( vector<int> &missed_bp, int missed_bp_avg ){
  for( int i=0; i<missed_bp.size(); i++ ){
    if( missed_bp[i] > missed_bp_avg ){
      matchlist.erase( matchlist.begin() + i );
      missed_bp.erase( missed_bp.begin() + i );
      i--;
    }
  }
}

// fourth step in create_extension: build extension string
string Contig::extension_build_string( int start, int pos_mult, int len, bool back ){
  string ATCGstr( "ATCG" );
  string extension = "";

  // loop len times processing 1 basepair at a time
  for( int i=0; i<len; i++ ){
    vector<int> ATCG_curr = { 0,0,0,0 };
    int max = 0;
    int avg = 0;
    
    // check coverage at current position
    if( check_cov( start+(pos_mult*i) ) < min_cov ){
      break;
    }
    
    for( int j=0; j<matchlist.size(); j++ ){
      int next_char = matchlist[j].getPos( start+(pos_mult*i) );

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
    for( int j=1; j<4; j++ ){
      if( ATCG_curr[j] > ATCG_curr[max] ) {
        max = j;
      }
    }

    // add next base
    if( back ){
      extension.append( ATCGstr.substr( max, 1 ) );
    }
    else{
      extension.insert( 0, ATCGstr.substr( max, 1 ) );
    }
  }

  return extension;
}

/// checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
/// used for off the front matching
string Contig::create_extension( int len, bool back ){
  // contains multiplier for position calculation
  int pos_mult = -1;
  int start = -1;
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
  // loop through bp's for initial pass to tally up missed nucleotides
  extension_bp_count( ATCG, start, pos_mult, len );
  
  ///////////////////////////////////////////
  // STEP 2: Mark Errant Reads For Removal //
  // loop through bp's to determine which reads, if any, should be eliminated from the matchlist
  // if nucleotide of read is < max for that position, count it against the read, otherwise don't count it
  extension_missed_count( ATCG, missed_bp, missed_bp_tot, start, pos_mult, len );

  if( matchlist.size() != 0 ){
    // calculate avg missed_bp's
    missed_bp_avg = missed_bp_tot / matchlist.size() + 1;
  }
  else{
    return "";
  }

  /////////////////////////////////
  // STEP 3: Remove Errant Reads //
  // loop through missed_bp list to eliminate any matches that have a missed level greater than the avg 
  extension_error_removal( missed_bp, missed_bp_avg );
  
  //////////////////////////////
  // STEP 4: Create Extension //
  return extension_build_string( start, pos_mult, len, back );
}

// extend performs loops iterations of create_extension with length extend_len of each extension, at each iteration the extension is added to contig, and uses contig_sub_len characters from the front or back of the contig, which end is determined by the boolean value of back
void Contig::extend( bool back ){
  string extension("");
  string contig_sub_str("");
 
  // skip over any contigs that present at least double coverage
  if( doub_cov ){
    return;
  }

  // get extension through create_extension
  if( back ){
    contig_sub_str = contig.substr( contig.length() - ( contig_sub_len ) );
  }
  else{
    contig_sub_str = contig.substr( 0, contig_sub_len );
  }

  Contig contig_sub( contig_sub_str, "temp" );
  extension = contig_sub.create_extension( extend_len, back );
  
  //cout << "extension:" << extension << endl;
  if( extension.length() == 0 ){
    return;
  }

  if( back ){
    contig.append( extension );
  }
  else{
    contig.insert( 0, extension );
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

// set doub_cov var
void Contig::set_doub_cov( bool doub_cov ){
  this->doub_cov = doub_cov;
}

// contig_fusion: Attempt to support fusion in case of possibly poorly constructed end.. returns new score from section in question
//    ::> contig object is the second while contig_ref is the first and the extension is being made off the front of the object
int Contig::check_fusion_support( string contig_ref ){
  string support_string( "" );
  int score = 0;
  int pos = contig_ref.length() - 1;
  int start = -1;

  // create reference string for numeric based additions to the extension string
  string ATCGstr( "ATCG" );
  matchlist.clear();

  match_contig_fr();

  // return if matches found is less than min_cov
  if( matchlist.size() < min_cov ){
    return 1.0;
  }

  // create missed bp's vector to keep track of how many bp's each read contains that are below the max at that position
  vector<int> mismatch( matchlist.size(), 0 );
  vector<int> ambiguous_bp( matchlist.size(), 0 );
  int mismatch_tot = 0;
  int mismatch_avg = 0;
  int cov = min_cov;
  int cmp_len = 0;

  // count mismatch's in each read match
  for( int i=0; i<contig_ref.length(); i++ ){
    int next_char_ref = contig_ref[pos-i];
    for( int j=0; j<matchlist.size(); j++ ){
      int next_char = matchlist[j].getPos( start-i );
      if( next_char == 'N' ){
        ambiguous_bp[j]++;
      }
      // check if next_char matches or exists in the current read
      else if( next_char != next_char_ref && next_char != -1 ){
        mismatch[j]++;
      }
    }
  }

  // remove reads with more than the max_missed
  extension_error_removal( mismatch, max_missed );

  // remove reads with more than 3 N's
  extension_error_removal( ambiguous_bp, 3 );
  
  // if matchlist is smaller than min_cov, bail, return false
  if( matchlist.size() < min_cov ){
    return 1.0;
  }

  // get string built from remaining matches.
  support_string = extension_build_string( start, -1, contig_ref.length(), false );
  cmp_len = support_string.length();

  // use the smaller of contig_ref.length() and min_overlap/2 as min_length to compare support_string with matching contig
  if( cmp_len < contig_ref.length() && cmp_len < (min_overlap/2 + min_overlap%2) ){
    return 1.0;
  }

  // get smallest of string lengths
  if( cmp_len > contig_ref.length() ){
    cmp_len = contig_ref.length();
  }

  score = mismatch_score( support_string.substr( support_string.length() - cmp_len), contig_ref.substr( contig_ref.length() - cmp_len ) );
  
  return score;
}
