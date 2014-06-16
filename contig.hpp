// $Author: benine $
// $Date$
// $Log$
// Contains the Contig class for afin


#ifndef CONTIG_HPP
#define CONTIG_HPP

using namespace std;

#include "read.hpp"

// return the read range corresponding to the first MAX_SORT characters of the read passed
tuple<long,long,long,long> get_read_range( string read_seg );
    
/////////////////////////////////////////////////////\
// Contig Class: //////////////////////////////////////>
//////////////////////////////////////////////////////
// Contains contigs and holds matches to those contigs during processing
// Contains methods for extending contigs
class Contig{
  private:
    vector<Read> matchlist;
    string contig;
    int first_read; // indicates the postition in contig where the first matching read begins    
    int min_cov;

  public:
    Contig( string str, int cov );

    Contig( string str );

    // adds a read to the read list
    void push_match( string read, int pos );

    void push_match( string read, int pos, bool revcomp );

    string getContig();

    int getListSize();

    string getRead_s( int i );

    int getStart( int i );

    Read getRead( int i );

    // clear matchlist to make room for new matches
    void clear_matches();

    // should this be a class? return an object containing how many matches of each and where from? Prolly not
		void check_pos( int pos );

    // finds the initial point at which the matches meet the min_cov requirement
    int find_start();

    // determines where the read passed matches the contig if at all for off the front matches
    void match_contig_fr();

    // determines where the read passed matches the contig if at all for off the back matches
    void match_contig_rr();
    
    /// checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
    /// used for off the front matching
    string check_match_fr( int len );

    // checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
    /// used for off the back matching
    string check_match_rr( int len );

    // extend() performs loops iterations of check_match with length extend_len of each extension, at each iteration the extension is added to contig, and uses contig_sub_len characters from the front or back of the contig, which end is determined by the boolean value back
    void extend( bool back );

    // checks the coverage of matches at the given positions, returns the coverage
    long check_cov( long pos );
};
//////////////////////
// End Contig Class //
//////////////////////

#endif
