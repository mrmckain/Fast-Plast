// $Author: benine $
// $Date$
// $Log$
// This contains the read class for afin

#include "read.hpp"
#include <cstdlib>
#include <unistd.h>

//////// READ FUNCTIONS //////////////
// constructor.. default revcomp will be false
Read::Read( std::string read, int match, bool rev ){
  this->read = read;
  start = match;
  misses = 0;
  this->rev = rev;

}

Read::Read( std::string read, int match ){
  this->read = read;
  start = match;
  misses = 0;
  rev = false;
}

// returns character at pos where pos is the nth postion of the match
char Read::get_pos( int pos, bool back ){
  if( back ){
    int read_pos = pos - start;
    if( read_pos >=0 && read_pos < read.length() ){
      return read[ read_pos ];
    }
  }
  else{
    int len = read.length();
    int read_pos = len - start + pos;
    if( read_pos >= 0 && read_pos < len ){
      return read[ read_pos ];
    }
  }

  return -1;
}
