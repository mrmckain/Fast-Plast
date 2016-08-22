// $Author: benine $
// $Date$
// $Log$
// This contains the revcomp which will return a reverse compliment of the string passed

#include <string>
#include <cstdio>
#include "revcomp.hpp"

// return the reverse compliment of str
std::string revcomp( std::string str ){
  std::string rc(str);
  int len = str.length();
  int i;

  for( i=0; i<len; i++ ){
    switch( str[i] ){
      case 'A':
        rc[ len - i - 1 ] = 'T';
        break;
      case 'T':
        rc[ len - i - 1 ] = 'A';
        break;
      case 'G':
        rc[ len - i - 1 ] = 'C';
        break;
      case 'C':
        rc[ len - i - 1 ] = 'G';
        break;
      case 'N':
        rc[ len - i - 1 ] = 'N';
        break;
      case 't':
        rc[ len - i - 1 ] = 'A';
        break;
      case 'a':
        rc[ len - i - 1 ] = 'T';
        break;
      case 'c':
        rc[ len - i - 1 ] = 'G';
        break;
      case 'g':
        rc[ len - i - 1 ] = 'C';
        break;
      default:
        printf( "Incorrect character present in string.\n" );
    }
  }
  return rc;
}
