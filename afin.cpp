// $Author: benine $
// $Date: 2014/04/18 19:09:09 $
// $Log: afin.cpp,v $
// Revision 1.4  2014/04/18 19:09:09  benine
// Cleaned up some and found jesus
//
// Revision 1.3  2014/04/15 01:15:40  benine
// Process class added
// Cleaned up processing
// Next step is add functionality for testing various segments of the contig being processed at the time
//
// Revision 1.2  2014/04/09 23:31:14  benine
// Process is built but throws compile errors due to nested classes attempting to access members of parent class
// Reorganizing classes now
//
// Revision 1.1  2014/04/09 23:24:43  benine
// Initial revision
//
// Revision 1.2  2014/04/02 02:06:11  benine
// Cleaned up some of the testing bits and added some more comments
//
// Revision 1.1  2014/04/01 17:23:50  benine
// Initial revision
//
// $Header: /home/benine/code/git/bb_afin/RCS/afin.cpp,v 1.4 2014/04/18 19:09:09 benine Exp benine $
//
// Coded and compiled using c++11 standard

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// TASK:: Add processing for IUPAC DNA ambiguity codes
// TASK:: Add processing for differences in reads( ie, create new contig objects for differing sets of matches, add method for splitting matchlist between two new contig objects ), determine which contig is correct
// TASK:: Check to see if it would be beneficial to leave the offthefront() function out of find_part(), maybe include it when determining if two contigs connect, but probably not
// TASK:: Add threading capability
// TASK:: Put the following functions in class
// pass by reference rc and put in it the reverse compliment of str
void revcomp( std::string& str, std::string& rc ){
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
}

// return the reverse compliment of str
string revcomp( std::string& str ){
  string rc(str);
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

// returns position if the search string aligns with a portion off the front of the base 
int offthefront( string base, string sub, int min ){
  int length = sub.length();
  for( int i=1; i<length-min; i++ ){
    if( base.substr( 0, length-i ) ==  sub.substr( i, length-i ) ){
      return( i );
    }
  }
  return 0;
}

// returns position if the search string aligns with a portion off the back of the base
int offtheback( string base, string sub, int min ){
  int length = base.length();
  int sublen = sub.length();

  // compare base and sub to see if the end of base lines up with the beginning of sub anywhere with a minimum of min characters in common
  for( int i=length - sublen + 1; i<length-min; i++ ){
    if( base.substr( i, length-i-2 ) == sub.substr( 0, length-i-2 ) ){
      return( -i );
    }
  }
  return 0;
}

//////////////////////////////////////\
// Read Class:  ///////////////////////>
//////////////////////////////////////
//  Contains one read and information about that read in reference to the piece of the contig it matched to
class Read{
  private:
    int start;      // position of read where match (ie piece of contig) starts
                    // for instance:
                    //    match: aaaattta
                    //    read:  tttatt
                    //    start: -4
    int misses;    // counts the number of characters that don't match other reads
    string read;    // contains the read or rev compliment, whichever is found
    bool revcomp;   // indicates if the read is a reverse compliment of the matched string

  public:
    // constructor.. default revcomp will be false
    Read( string& read, int match, bool revcomp );

    Read( string& read, int match );
    
    // returns character at pos where pos is the nth postion of the match
    char getPos( int pos );

    int getStart();

    string getRead();
};

////////////////////
// End Read Class //
////////////////////



/////////////////////////////////////////////////////\
// Contig Class: //////////////////////////////////////>
//////////////////////////////////////////////////////
// Contains contigs and holds matches to those contigs during processing
// Contains methods for extending contigs
class Contig_c{
  private:
    vector<Read> matchlist;
    string contig;
    int first_read; // indicates the postition in contig where the first matching read begins    
    int min_cov; // contains the minimum coverage variable

  public:
    Contig_c( string str, int cov );

    Contig_c( string str );

    // adds a read to the read list
    void push_read( string read, int pos );

    void push_read( string read, int pos, bool revcomp );

    string getContig();

    int getListSize();

    string getRead_s( int i );

    int getStart( int i );

    Read getRead( int i );

    // clear matchlist to make room for new matches
    void clear_matches();


    // searches for read in contig with at least min characters in common, returns position in **pos or sets *pos to NULL if sub is not found anywhere
		// Checks reverse compliment of string if original fails
    // possible issues: sub may have less characters than min
		int find_part( string read, string contig_sub, int min );
		
		// overload find_part() to have a default of min = 10
		int find_part( string match, string contig_sub );

    // should this be a class? return an object containing how many matches of each and where from? Prolly not
		void check_pos( int pos );

    // finds the initial point at which the matches meet the min_cov requirement
    int find_start();

    // find reads from readlist that match the contig, the portion of the contig is indicated by sub_len, 0 indicates use of the full contig
    void find_match();

    // checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
    string check_match( int len );

    // check_match( int ) with a default len provided
    string check_match();

    // extend() performs loops iterations of check_match with length len of each extension, at each iteration the extension is added to contig, and uses sublen characters from the end of the contig
    void extend( int loops, int len, int sublen );

    // checks the coverage of matches at the given positions, returns the coverage
    int check_cov( int pos );
};
//////////////////////
// End Contig Class //
//////////////////////
    
    
/////////////////////////////////////////////\
// Process Class: ////////////////////////////>
/////////////////////////////////////////////
//
// create object here containing both Read and Contig_c objects
// Process object, contains all classes, methods, data, and data references necessary for processing the contigs
// There will be only one Process object needed per iteration of this program
class Process{
  private:
  public:
    vector<Contig_c> contigs;
    static vector<string> readlist;

    Process(){}

    // put reads from readfile into readlist
    void add_reads( string filename ){
		  string buffer("");
      string line("");
		  
      // open read file
      ifstream read( filename );
      
      // read in reads
		  while( getline( read, line )){
		    if( line[0] == '>' && buffer.length() != 0 ){
		      readlist.push_back( buffer );
		      buffer = "";
		    } 
		    else if( line[0] == '>' ){
		    }
		    else {
		      buffer += line;
		    }
      }

      // close read file
      read.close();

      // insert last line into readlist
		  if( buffer.length() != 0 ){
		    readlist.push_back( buffer );
		  }
    }

    // put contigs from contfile into contlist
    void add_contigs( string filename ){
		  string buffer("");
      string line("");
		  
      // open contig file
      ifstream cont( filename );

      // read in contig objects
		  while( getline( cont, line ) ){
		    if( line.length()>0 && line[0] != '>' && buffer.length() != 0 ){
		      contigs.push_back( Contig_c( buffer, 3 ));
		      buffer = "";
		    }
		    else {
          if( line.back() == '>' )
		        line.pop_back();
		      buffer += line;
		    }
		  }
      
      // insert last line into contigs list
		  if( buffer.length() != 0 ){
		    contigs.push_back( Contig_c( buffer, 3 ) );
		  }
		  
      // close contig file
      cont.close();
    }

    // return contig with index contig_ind
    string get_contig( int contig_ind ){
      return contigs[ contig_ind ].getContig();
    }
};

///////////////////////
// End Process Class //
///////////////////////







/////////////////////////////////////////////////////////////////////////////////\
// BEGIN MAIN FUNCTION ///////////////////////////////////////////////////////////>
/////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv ){
  //////////////////////////////////
  // Process Command Line Options //
  //////////////////////////////////
  int c;
  Process process;
  
  
  // prevent output to stderr if erroneous option is found
  opterr = 0;

  // get all options that have been provided on the command line
  while (( c = getopt (argc, argv, "hr:c:" )) != -1 ) {
    switch( c ) {
      case 'h':
        cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        break;
      case 'r':
        cout << "readfile: " << optarg << endl;
        process.add_reads( optarg );
        break;
      case 'c':
        cout << "contigfile: " << optarg << endl;
        process.add_contigs( optarg );
        break;
      case '?':
        if ( optopt == 'r' ){
          fprintf( stderr, "%s: Error: Option -r requires an argument. ", argv[0] );
          cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        }
        else if ( optopt == 'c' ){
          fprintf( stderr, "%s: Error: Option -c requires an argument. ", argv[0] );
          cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        }
        else if ( isprint( optopt )){
          fprintf( stderr, "%s: Error: Unknown option -%c. \n", argv[0], optopt );
          cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        }
        else{
          fprintf( stderr, "%s: Error: Unknown option character %x.\n", argv[0], optopt );
          cout << "Usage: " << argv[0] << " -c [contigfile(s)] -r [readfile(s)]\n";
        }
        return 1;
      default:
        abort();
    }
  }
  
  while ( optind < argc ) {
    cout << argv[optind] << endl;
    optind++;
  }
  /////////////////
  // End Options //
  /////////////////



  /////////////////////////
  // Test Process Class ///
  /////////////////////////

  cout << "THIS IS A PROCESS CLASS TEST! THIS IS ONLY A TEST!" << endl;


  cout << "contig: " << process.get_contig(0) << endl;

  process.contigs[0].extend( 4, 10, 100 );

  cout << "contex: " << process.get_contig(0) << endl;


  cout << "THIS HAS BEEN A PROCESS CLASS TEST. THANK YOU FOR YOUR PATIENCE!" << endl;

  //////////////
  // End Test //
  //////////////


  return 0;
}

//////////////
// END MAIN //
//////////////
    






///////////////////////////////////////////////////////////
//////// CLASS FUNCTION DEFINITIONS ///////////////////////
///////////////////////////////////////////////////////////


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

////////// CONTIG_C FUNCTIONS ////////////
Contig_c::Contig_c( string str, int cov ) : min_cov( cov ), contig( str ){}

Contig_c::Contig_c( string str ) : contig( str ){
  min_cov = 10;
}   

// adds a read to the read list
void Contig_c::push_read( string read, int pos ){
  matchlist.push_back(Read( read, pos ));
}

void Contig_c::push_read( string read, int pos, bool revcomp ){
  matchlist.push_back(Read( read, pos, revcomp ));
}

string Contig_c::getContig(){
  return contig;
}

int Contig_c::getListSize(){
  return matchlist.size();
}

string Contig_c::getRead_s( int i ){
  return matchlist[i].getRead();
}

int Contig_c::getStart( int i ){
  return matchlist[i].getStart();
}

Read Contig_c::getRead( int i ){
  return matchlist[i];
}

// clear matchlist to make room for new matches
void Contig_c::clear_matches(){
  matchlist.clear();
}

// searches for read in contig with at least min characters in common, returns position in **pos or sets *pos to NULL if sub is not found anywhere
// Checks reverse compliment of string if original fails
// possible issues: sub may have less characters than min
int Contig_c::find_part( string match, string contig_sub, int min ){
  int pos;
  string rc("");  
  
  // check for normal matches
  if( (pos = -contig_sub.find( match )) != 1 ){
  	    push_read( match, pos );
    return 1;
  }
  if( (pos = offthefront( contig_sub, match, min )) != 0){
  	    push_read( match, pos );
    return 1;
  }
  if( (pos = offtheback( contig_sub, match, min )) != 0){
  	    push_read( match, pos );
    return 1;
  }
	
  rc = revcomp( match );
  // check for reverse compliment matches
  if( (pos = -contig_sub.find( rc )) != 1 ){
  	    push_read( rc, pos, true );
    return 1;
  }
  if( (pos = offthefront( contig_sub, rc, min )) != 0){
  	    push_read( rc, pos, true );
    return 1;
  }
  if( (pos = offtheback( contig_sub, rc, min )) != 0){
  	    push_read( rc, pos, true );
    return 1;
  }

  return 0;
}

// overload find_part() to have a default of min = 10
int Contig_c::find_part( string match, string contig_sub ){
  return find_part( match, contig_sub, 10 );
}

// TASK:: Function to process alignment of matches
//    -> Find first starting place
//      -> based on number of reads that match the region searched for, only process portions where x% of coverage of matches
//    -> assess count of each bp at each position
//      ->keep track of misses in each read
//    -> possibility :::-->>> maintain multiple read lists when there is a dispute over the proper character at a position
//      -> It is possible that two branches comeback together
//      -> Rate value of each branch using weighted values 
//    -> take most frequently found character, potentially using weighted values involving how many misses in that read
//
//    -> Every x characters examined, if all looks good, add those characters to the contig and search reads again, adding to the current matchlist
//
//    -> Analyze read headings and assign weights of some sort based on coverage of each read if possible
//    -> Analyze contig headings and make decisions based on coverage of each contig
//

void Contig_c::check_pos( int pos ){
  // TASK:: should this be a class? return an object containing how many matches of each and where from? Prolly not
      
}

// finds the initial point at which the matches meet the min_cov requirement
int Contig_c::find_start(){
  int cov;
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

// find reads from readlist that match the contig
void Contig_c::find_match(){
  // cycle through readlist to check for matching reads
  for( int j=0; j<Process::readlist.size(); j++ ){
    find_part( Process::readlist[j], contig );
  }
}

// checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
string Contig_c::check_match( int len ){
  int start = contig.length();

  matchlist.clear();
  find_match();

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

// check_match( int ) with a default len provided
string Contig_c::check_match(){
  return check_match( 20 );
}

// extend() performs loops iterations of check_match with length len of each extension, at each iteration the extension is added to contig, and uses sublen characters from the end of the contig
void Contig_c::extend( int loops, int len, int sublen ){
  string extension("");
  
  // if sublen=0 use entire contig, else use the end of the contig len characters long and create new contig object to process and get extension from
  if( sublen > 0 ){
    for( int i=0; i<loops; i++ ){
      Contig_c contig_sub( contig.substr( contig.length() - ( sublen + 1 ) ), 3 );
   
      // get extension through check_match
      extension = contig_sub.check_match( len );
      contig.append( extension );
    }
  }
  else {
    for( int i=0; i<loops; i++ ){
      extension = check_match( len );
      contig.append( extension );
    }
  }
}

// checks the coverage of matches at the given positions, returns the coverage
int Contig_c::check_cov( int pos ){
  int cov = 0;

  for( int i=0; i<matchlist.size(); i++ ){
    if( matchlist[i].getPos( pos ) != -1 ){
      cov++;
    }
  }
  return cov;
}

//////// PROCESS DEFINITIONS ///////
vector<string> Process::readlist;

//////////////////////////////
// END DEFINITIONS ///////////
//////////////////////////////
