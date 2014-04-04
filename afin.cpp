// $Author: benine $
// $Date: 2014/04/02 02:06:11 $
// $Log: afin.cpp,v $
// Revision 1.2  2014/04/02 02:06:11  benine
// Cleaned up some of the testing bits and added some more comments
//
// Revision 1.1  2014/04/01 17:23:50  benine
// Initial revision
//
// $Header: /home/benine/code/git/afin/RCS/afin.cpp,v 1.2 2014/04/02 02:06:11 benine Exp benine $
//
// Coded and compiled using c++11 standard


#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

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
  for( int i=length - sublen + 1; i<length-min; i++ ){
    if( base.substr( i, length-i ) == sub.substr( 0, length-i ) ){
      return( -i );
    }
  }
  return 0;
}

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
    vector<ifstream*> readfile;
    vector<ifstream*> contfile;

    vector<Contig_c> contigs;
    vector<string> readlist;

    Process(){}

    // add read file inputs
    void add_read( string r_filename ){
      // add read file to readfile vector
      ifstream read( r_filename );
      contfile.push_back( &read );
    }

    // add contig file inputs
    void add_cont( string c_filename ){
      // add contig file to contfile vector
      ifstream cont( c_filename );
      contfile.push_back( &cont );
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
		    Read( string& read, int match, bool revcomp ){
		      this->read = read;
		      start = match;
		      misses = 0;
		      this->revcomp = revcomp;
		
		    }
		
		    Read( string& read, int match ){
		      this->read = read;
		      start = match;
		      misses = 0;
		      revcomp = false;
		    }
		    
		    // returns character at pos where pos is the nth postion of the match
		    char getPos( int pos ){
		      if( pos+start >=0 && pos+start < read.length() ){
		        return read[ pos + start ];
		      }
		      return -1; 
		    }
		
		    int getStart(){
		      return start;
		    }
		
		    string getRead(){
		      return read;
		    }
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
		    Contig_c( string str, int cov ) : min_cov( cov ), contig( str ){}
		
		    Contig_c( string str ) : contig( str ){
		      min_cov = 10;
		    }   
		
		    // adds a read to the read list
		    void push_read( string read, int pos ){
		      matchlist.push_back(Read( read, pos ));
		    }
		
		    void push_read( string read, int pos, bool revcomp ){
		      matchlist.push_back(Read( read, pos, revcomp ));
		    }
		
		    string getContig(){
		      return contig;
		    }
		
		    int getListSize(){
		      return matchlist.size();
		    }
		
		    string getRead_s( int i ){
		      return matchlist[i].getRead();
		    }
		
		    int getStart( int i ){
		      return matchlist[i].getStart();
		    }
		
		    Read getRead( int i ){
		      return matchlist[i];
		    }
		
		    // clear matchlist to make room for new matches
		    void clear_matches(){
		      matchlist.clear();
		    }
		
		
		    // searches for read in contig with at least min characters in common, returns position in **pos or sets *pos to NULL if sub is not found anywhere
				// Checks reverse compliment of string if original fails
		    // possible issues: sub may have less characters than min
				int find_part( string read, int min ){
				  int pos;
		      string rc("");  
				  
		      // check for normal matches
		      if( (pos = -contig.find( read )) != 1 ){
		  	    push_read( read, pos );
				    return 1;
				  }
				  if( (pos = offthefront( contig, read, min )) != 0){
		  	    push_read( read, pos );
				    return 1;
				  }
				  if( (pos = offtheback( contig, read, min )) != 0){
		  	    push_read( read, pos );
				    return 1;
				  }
			
		      rc = revcomp( read );
		      // check for reverse compliment matches
				  if( (pos = -contig.find( rc )) != 1 ){
		  	    push_read( rc, pos, true );
				    return 1;
				  }
				  if( (pos = offthefront( contig, rc, min )) != 0){
		  	    push_read( rc, pos, true );
				    return 1;
				  }
				  if( (pos = offtheback( contig, rc, min )) != 0){
		  	    push_read( rc, pos, true );
				    return 1;
				  }
		
				  return 0;
				}
				
				// overload find_part() to have a default of min = 10
				int find_part( string match ){
				  return find_part( match, 10 );
				}
		
		    // Function to process alignment of matches
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
				void check_pos( int pos ){
		      // should this be a class? return an object containing how many matches of each and where from? Prolly not
		      
		    }
		
		    // finds the initial point at which the matches meet the min_cov requirement
		    int find_start(){
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
            // check each read for a character at the current position
		        for( int i=0; i<matchlist.size(); i++ ){
		          next_char = matchlist[i].getPos( j );
		          if ( next_char != -1 ) {
		            cov++;
		          }
		        }
		        
		        // if coverage is at the minimum level, return the current position
            if( cov >= min_cov ){
		          return( j );
		        }
		      }
		      
		      return(contig.length() + 1);
		    }
		
		    // find reads from readlist that match the contig, the portion of the contig is indicated by sub_len, 0 indicates use of the full contig
		    void find_match( int sub_len ){
		      for( int j=0; j<readlist.size(); j++ ){
		        contig.find_part( readlist[j] );
		      }
		    }
		
		    // checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
		    string check_match( int len ){
		      int start = find_start();
		
		      matchlist.clear();
		      find_match( 0 );
		
		      // create reference string for numeric based additions to the extension string
		      string ATCGstr( "ATCG" );
		      int i;
		      string extension( "" );
		
		      // loop len times processing 1 basepair at a time
		      for( int j=0; j<len; j++ ){
		        int ATCG[] = { 0,0,0,0 };
		        int max = 0;
		        int avg = 0;
		        int cov = 0;
		        for( i=0; i<matchlist.size(); i++ ){
		          int next_char = matchlist[i].getPos( start+j );
		
		          if ( next_char == 'A' ) {
		            cov++;
		            ATCG[0]++;
		          }
		          else if ( next_char == 'T' ) {
		            ATCG[1]++;
		            cov++;
		          }
		          else if ( next_char == 'C' ) {
		            ATCG[2]++;
		            cov++;
		          }
		          else if ( next_char == 'G' ) {
		            ATCG[3]++;
		            cov++;
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
		
		        if( cov < min_cov ){
		          cout << "cov<min_cov\n";
		          break;
		        }
		      }
		      
		      return extension;
		    }
		
		    // check_match( int ) with a default len provided
		    string check_match(){
		      return check_match( 20 );
		    }
		
		
		    // extend() performs loops iterations of check_match with length len of each extension, at each iteration the extension is added to contig
		    void extend( int loops, int len ){
		      string extension("");
		      for( int i=0; i<loops; i++ ){
		        extension = check_match( len );
		        contig.append( extension );
		      }
		    }
		};
    //////////////////////
    // End Contig Class //
    //////////////////////
    
    
    
    // put reads from readfile into readlist
    void get_reads(){
		  string buffer("");
      string line("");
		  
      // read in reads
		  while( getline( *readfile[0], line )){
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
		  
      // insert last line into readlist
		  if( buffer.length() != 0 ){
		    readlist.push_back( buffer );
		  }
    }

    // put contigs from contfile into contlist
    void get_contigs(){
		  string buffer("");
      string line("");
		  
      // read in contig objects
		  while( getline( *contfile, line ) ){
		    if( line.length()>0 && line[0] != '>' && buffer.length() != 0 ){
		      contigs.push_back( Contig_c( buffer ));
		      buffer = "";
		    }
		    else {
		      line.pop_back();
		      buffer += line;
		    }
		  }
		  
      // put last line into contig list
      contigs.push_back( Contig_c( buffer, 3 ));
    }

    // close open read and contig files
    void close_files(){
      int i=0;
      for( i=0; i<contfile.length(); i++ ){
        *contfile[i].close();
      }

      for( i=0; i<readfile.length(); i++ ){
        *readfile[i].close();
      }
    }

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
  //vector<Contig_c> contigs;
  //vector<string> readlist;


  /////////////////////////
  // Test Process Class ///
  /////////////////////////

  cout << "THIS IS A PROCESS CLASS TEST! THIS IS ONLY A TEST!" << endl;


  Process process;
  process.add_cont( "contig.fasta" );
  process.add_read( "reads.fasta" );
  
  process.get_contigs();
  process.get_reads();

  process.close_files();

  cout << "contig: " << process.get_contig(0) << endl;

  process.contigs[0].extend( 4, 10 );

  cout << "contex: " << process.get_contig(0) << endl;


  cout << "THIS HAS BEEN A PROCESS CLASS TEST. THANK YOU FOR YOUR PATIENCE!" << endl;

  //////////////
  // End Test //
  //////////////



/*
  /// READ FILES INTO STORAGE
  // This section will read two fasta files, one containing the contig(s) to be expanded, the other containing the reads
  ifstream contfile( "contig.fasta" );
  ifstream readfile( "reads.fasta" );
  
  string buffer("");
  string line("");
  string extension("");


  // read in contigs
  while( getline( contfile, line ) ){
    if( line.length()>0 && line[0] != '>' && buffer.length() != 0 ){
      contigs.push_back( Contig_c( buffer ));
      buffer = "";
    }
    else {
      line.pop_back();
      buffer += line;
    }
  }
  contigs.push_back( Contig_c( buffer, 3 ));

  buffer = "";

  // read in reads
  while( getline( readfile, line )){
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
  
  if( buffer.length() != 0 ){
    readlist.push_back( buffer );
    buffer = "";
  }

  /// CHECK READS ///

  contfile.close();
  readfile.close();
  // END READ //

  
  // CHECK READS AGAINST CONTIGS //
  // This will check reads against contigs and store matching reads in corresponding Contig object 
  for( int i=0; i<contigs.size(); i++ ){
    for( int j=0; j<readlist.size(); j++ ){
      contigs[i].find_part( readlist[j] );
    }
  }
  // END CHECK //


  // DISPLAY MATCH LIST //
  // This will display any matches currently stored in the vector
  for( int i=0; i<contigs.size(); i++ ){
    for( int j=0; j<contigs[i].getListSize(); j++ ){
      cout << "FileRead: " << contigs[i].getRead( j ).getRead() << endl << "start: " << contigs[i].getRead(j).getStart() << "  pos(180): " << contigs[i].getRead(j).getPos(180) << endl;
    }
  }
  // END DISPLAY //

  // display results of check_match one time
  extension = contigs[0].check_match();
  cout << "EXTENSION:\t|" << extension << "|\n";

  cout << "Contig: " << contigs[0].getContig() << endl;

  contigs[0].extend( 4, 10 );

  cout << "Contex: " << contigs[0].getContig() << endl;
  */

  return 0;
}

//////////////
// END MAIN //
//////////////
