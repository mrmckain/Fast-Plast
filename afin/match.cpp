// object for finding and storing matches to the contig in which the object is stored

#include "match.hpp"
#include "process.hpp"
#include "revcomp.hpp"

Match::Match( Readlist *reads ) : reads(reads) {
  back = false;
  sequence = "";
}

Match::Match( Readlist *reads, std::string sequence ) : reads(reads), sequence(sequence){
  back = false;
}

// set sequence and back variables
void Match::set_seq( std::string sequence, bool back ){
  this->back = back;
  this->sequence = sequence;
}

// return matchlist length
size_t Match::get_matchlist_size(){
  return matchlist.size();
}

// return char in matchlist read at pos
char Match::get_pos( int read, int pos ){
  return matchlist[read].get_pos(pos, back);
}

// adds a read to the read list
void Match::push_match( std::string read, int pos ){
  matchlist.push_back(Read( read, pos ));
}

void Match::push_match( std::string read, int pos, bool rev ){
  matchlist.push_back(Read( read, pos, rev ));
}

// perform match
void Match::start_match(){
  if( back ){
    match_contig_rr();
  }
  else{
    match_contig_fr();
  }
}

// determines where the read passed matches the contig if at all for off the front matches
void Match::match_contig_fr(){
  // loop through possible substrings of sequence to check for matches in readlist
  for( int i=sequence.length()-2; i>=min_overlap; i-- ){
    std::string seq_sub_rc( sequence.substr( 0, i ));
    seq_sub_rc = revcomp( seq_sub_rc );
    std::tuple<long,long,long,long> range = reads->get_read_range( seq_sub_rc.substr( 0, max_sort_char ) );

    if( std::get<0>(range) != -1 ){
      // check reads in range against sequence
      if( std::get<0>(range) != 0 ){
        for( int j=std::get<0>(range)-1; j<std::get<1>(range); j++ ){
          // check if the current read matches the sequence from this point
          if( reads->get_read(j).compare( 0, seq_sub_rc.length(), seq_sub_rc ) == 0 ){
            std::string read_rc = revcomp( reads->get_read(j) );
            push_match( read_rc, i, true );
          }
        }
      }

      // check if there are revcomp reads to be checked
      if( std::get<2>(range) != 0 ){
        // check reverse complements
        for( int j=std::get<2>(range)-1; j<std::get<3>(range); j++ ){
          // check if the current read matches the sequence from this point
          std::string rc = reads->get_rc_read(j);

          if( rc.compare( 0, seq_sub_rc.length(), seq_sub_rc ) == 0 ){
            push_match( rc, i );
          }
        }
      }
    }
  }
}

// determines where the read passed matches the contig if at all for off the back matches
void Match::match_contig_rr(){
  // loop through possible substrings of sequence to check for matches in readlist
  for( int i=1; i<sequence.length()-min_overlap; i++ ){
    std::string seq_sub( sequence.substr( i, sequence.length()-1 ));
    std::tuple<long,long,long,long> range = reads->get_read_range( seq_sub.substr( 0, max_sort_char ) );

    if( std::get<0>(range) != -1 ){
      // check reads in range against sequence
      if( std::get<0>(range) != 0 ){
        for( int j=std::get<0>(range)-1; j<std::get<1>(range); j++ ){
          // check if the current read matches the sequence from this point
          if( reads->get_read(j).compare( 0, seq_sub.length(), seq_sub ) == 0 ){
            push_match( reads->get_read(j), i );
          }
        }
      }

      // check if there are revcomp reads to be checked
      if( std::get<2>(range) != 0 ){
        // check reverse complements
        for( int j=std::get<2>(range)-1; j<std::get<3>(range); j++ ){
          // check if the current read matches the sequence from this point
          std::string rc = reads->get_rc_read(j);

          if( rc.compare( 0, seq_sub.length(), seq_sub ) == 0 ){
            push_match( rc, i, true );
          }
        }
      }
    }
  }
}

// checks the coverage of matches at the given positions, returns the coverage
int Match::check_cov( int pos ){
  int cov = 0;

  for( int i=0; i<matchlist.size(); i++ ){
    if( matchlist[i].get_pos( pos, back ) != -1 ){
      cov++;
    }
  }
  return cov;
}

// remove match at index ind
void Match::remove_match( int ind ){
  matchlist.erase( matchlist.begin() + ind );
}

// clear match list
void Match::clearlist(){
  matchlist.clear();
}
