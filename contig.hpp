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
    string contig_id;
    int first_read; // indicates the postition in contig where the first matching read begins    
    int min_cov;

  public:
    Contig( string str, string id, int cov );

    Contig( string str, string id );

    // adds a read to the read list
    void push_match( string read, int pos );

    void push_match( string read, int pos, bool revcomp );

    // return contig
    string getContig();
    
    // return contig_id
    string get_contig_id();

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

    // checks each matched read against the contig one bp at a time and filters out poorly aligning reads
    string check_match();

    // determines where the read passed matches the contig if at all for off the front matches
    void match_contig_fr();

    // determines where the read passed matches the contig if at all for off the back matches
    void match_contig_rr();
    
    /// checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
    string create_extension( int len, bool back );

    // extend() performs loops iterations of create_extension with length extend_len of each extension, at each iteration the extension is added to contig, and uses contig_sub_len characters from the front or back of the contig, which end is determined by the boolean value back
    void extend( bool back );

    // checks the coverage of matches at the given positions, returns the coverage
    long check_cov( long pos );
};
//////////////////////
// End Contig Class //
//////////////////////

#endif
