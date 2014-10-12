// $Author: benine $
// $Date$
// $Log$
// This contains the read class for afin

#include <cstdlib>
#include <string>
#include <unistd.h>
#include "read.hpp"

using namespace std;

//////// READ FUNCTIONS //////////////
// constructor.. default revcomp will be false
Read::Read( string& read, int match, bool revcomp ){
  this->read = read;
  start = match;
  misses = 0;
  this->revcomp = revcomp;

}

Read::Read( string& read, int match ){
  this->read = read;
  start = match;
  misses = 0;
  revcomp = false;
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
