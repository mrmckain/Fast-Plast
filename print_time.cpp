// $Author: benine $
// $Date$
// $Log$
// Contains the print_time function to be used for diagnostic timekeeping

#include <iostream>
#include <ctime>
#include "print_time.hpp"

using namespace std;

void print_time(){
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  cout << now->tm_hour << ':' 
        << now->tm_min << ':'
        <<  now->tm_sec
        << endl;
}
