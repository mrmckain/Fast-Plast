

#include "extension.hpp"
#include "process.hpp"

Extension::Extension( Readlist *reads, int len ) : reads(reads), len(len), matches(Match(reads)){
  start = -1;
  missed_bp_tot = 0;
  missed_bp_avg = 0;
  back = false;
  contig = "";
  exten_seq = "";
  pos_mult = -1;
}

Extension::Extension( Readlist *reads, int len, std::string contig ) : reads(reads), len(len), matches(Match(reads,contig)){
  start = -1;
  missed_bp_tot = 0;
  missed_bp_avg = 0;
  back = false;
  contig = "";
  exten_seq = "";
  pos_mult = -1;
}

// set the value of the missed_bp vector
void Extension::set_missed_bp( std::vector<int> missed_bp ){
  this->missed_bp = missed_bp;
}

// set the value of missed_bp_avg
void Extension::set_missed_avg( int missed_avg ){
  missed_bp_avg = missed_avg;
}

// first step of get_extension: determine bp count and max represented bp at each position
void Extension::bp_count(){
  for( int i=0; i<len; i++ ){
    // initialize temporary vector
    std::vector<int> ATCG_curr = { 0,0,0,0,0 };
    int next_char = 0;

    // loop through matches to get count of each nucleotide present at the current position
    for( int j=0; j<matches.get_matchlist_size(); j++ ){
      next_char = matches.get_pos( j, start+(pos_mult*i) );

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

// second step of the get_extension process: count missed bp's per read, or in other words, the bp's represented below the max for that position
void Extension::missed_count(){
  for( int i=0; i<len; i++ ){
    int next_char = 0;
    for( int j=0; j<matches.get_matchlist_size(); j++ ){
      next_char = matches.get_pos( j, start+(pos_mult*i) );
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

// third step in get_extension(): removal of reads that have errors over the threshold
bool Extension::error_removal(){
  int start_size = (int)matches.get_matchlist_size();
  for( int i=0; i<missed_bp.size(); i++ ){
    if( missed_bp[i] > missed_bp_avg ){
      matches.remove_match( i );
      missed_bp.erase( missed_bp.begin() + i );
      i--;
    }
  }
  if( start_size > 5 && matches.get_matchlist_size() < start_size * stop_ext ){
    return false;
  }

  return true;
}

// fourth step in get_extension: build extension sequence
void Extension::build_string(){
  std::string ATCGstr( "ATCG" );

  // loop len times processing 1 basepair at a time
  for( int i=0; i<len; i++ ){
    std::vector<int> ATCG_curr = { 0,0,0,0 };
    int max = 0;
    int avg = 0;

    // check coverage at current position
    if( matches.check_cov( start+(pos_mult*i) ) < min_cov ){
      break;
    }

    int next_char = 0;

    for( int j=0; j<matches.get_matchlist_size(); j++ ){
      next_char = matches.get_pos( j, start+(pos_mult*i) );

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
      exten_seq.append( ATCGstr.substr( max, 1 ) );
    }
    else{
      exten_seq.insert( 0, ATCGstr.substr( max, 1 ) );
    }
  }
}

/// checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned
/// used for off the front matching
std::string Extension::get_extension( std::string contig, bool back ){
  // contains multiplier for position calculation
  exten_seq = "";
  missed_bp_tot = 0;
  missed_bp_avg = 0;
  this->back = back;
  this->contig = contig;

  // reset parameters of Match object
  matches.set_seq( contig, back );
  matches.start_match();

  // get matches
  if( back ){
    start = contig.length();
    pos_mult = 1;
  }
  else{
    pos_mult = -1;
    start = -1;
  }

  // return if no matches are found
  if( matches.get_matchlist_size() == 0 ){
    // reset lists
    ATCG.clear();
    missed_bp.clear();
    matches.clearlist();
    return "";
  }

  // create missed bp's vector to keep track of how many bp's each read contains that are below the max at that position
  missed_bp.resize( matches.get_matchlist_size(), 0 );

  /////////////////////////////////////
  // STEP 1: First Pass Over Matches //
  // loop through bp's for initial pass to tally up missed nucleotides
  bp_count();

  ///////////////////////////////////////////
  // STEP 2: Mark Errant Reads For Removal //
  // loop through bp's to determine which reads, if any, should be eliminated from the matchlist
  // if nucleotide of read is < max for that position, count it against the read, otherwise don't count it
  missed_count();

  if( matches.get_matchlist_size() != 0 ){
    // calculate avg missed_bp's
    missed_bp_avg = missed_bp_tot / matches.get_matchlist_size() + 1;
  }
  else{
    // reset lists
    ATCG.clear();
    missed_bp.clear();
    matches.clearlist();
    return "";
  }

  /////////////////////////////////
  // STEP 3: Remove Errant Reads //
  // loop through missed_bp list to eliminate any matches that have a missed level greater than the avg
  if ( !error_removal() ){
    // reset lists
    ATCG.clear();
    missed_bp.clear();
    matches.clearlist();
    return "";
  }

  //////////////////////////////
  // STEP 4: Build Extension //
  build_string();

  // reset lists
  ATCG.clear();
  missed_bp.clear();
  matches.clearlist();
  return exten_seq;
}

// simply returns the built extension string
std::string Extension::get_extension(){
  return exten_seq;
}
