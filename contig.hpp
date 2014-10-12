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
    double cov;
    bool doub_cov;

  public:
    Contig( string str, string id, double cov, int min_cov );
    
    Contig( string str, string id );

    // adds a read to the read list
    void push_match( string read, int pos );

    void push_match( string read, int pos, bool revcomp );

    // return contig
    string get_contig();

    // return contig_id
    string get_contig_id();

    // return cov
    double get_cov();

    // clear matchlist to make room for new matches
    void clear_matches();

    // determines where the read passed matches the contig if at all for off the front matches
    void match_contig_fr();

    // determines where the read passed matches the contig if at all for off the back matches
    void match_contig_rr();
    
    // first step of create_extension: determine bp count and max represented bp at each position
    void extension_bp_count( vector<vector<int>> &ATCG, int start, int pos_mult, int &len, bool back );

    // second step of the create_extension process: count missed bp's per read, or in other words, the bp's represented below the max for that position
    void extension_missed_count( vector<vector<int>> &ATCG, vector<int> &missed_bp, int &missed_bp_tot, int start, int pos_mult, int len, bool back );

    // third step in create_extension(): removal of reads that have errors over the threshold
    void extension_error_removal( vector<int> &missed_bp, int missed_bp_avg );

    // fourth step in create_extension: build extension string
    string extension_build_string( int start, int pos_mult, int len, bool back );

    /// checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
    string create_extension( int len, bool back );

    // extend() performs loops iterations of create_extension with length extend_len of each extension, at each iteration the extension is added to contig, and uses contig_sub_len characters from the front or back of the contig, which end is determined by the boolean value back
    void extend( bool back );

    // checks the coverage of matches at the given positions, returns the coverage
    long check_cov( long pos, bool back );

    // set doub_cov var
    void set_doub_cov( bool doub_cov );

    // contig_fusion: Attempt to support fusion in case of possibly poorly constructed end.. returns new score from section in question
    //    ::> contig object is the second while contig_ref is the first and the extension is being made off the front of the object
    int check_fusion_support( string contig_ref );
};
//////////////////////
// End Contig Class //
//////////////////////

#endif
