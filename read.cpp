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
char Read::getPos( int pos ){
  if( pos+start >=0 && pos+start < read.length() ){
    return read[ pos + start ];
  }
  return -1; 
}

int Read::getStart(){
  return start;
}

string Read::getRead(){
  return read;
}
