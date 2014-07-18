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
    int bp_added_fr;
    int bp_added_rr;

  public:
    Contig( string str, string id, int cov, int init_added_fr, int init_added_rr );
    
    Contig( string str, string id, int cov, int bp_added_init );
    
    Contig( string str, string id, int cov );

    Contig( string str, string id );

    // adds a read to the read list
    void push_match( string read, int pos );

    void push_match( string read, int pos, bool revcomp );

    // return contig
    string getContig();

    // set contig_id
    void set_contig_id( string new_contig_id );
    
    // return contig_id
    string get_contig_id();

    int getListSize();

    string getRead_s( int i );

    int getStart( int i );

    Read getRead( int i );

    // returns bp_added_fr
    int get_bp_added_fr();

    // returns bp_added_rr
    int get_bp_added_rr();

    // resets bp_added variables to 0
    int reset_bp_added();
    
    // clear matchlist to make room for new matches
    void clear_matches();

    // finds the initial point at which the matches meet the min_cov requirement
    int find_start();

    // returns value associated with index of ATCG[] or -1 if not a member of ATCG
    int get_ATCG_value( int ATCG_char );

    // determines where the read passed matches the contig if at all for off the front matches
    void match_contig_fr();

    // determines where the read passed matches the contig if at all for off the back matches
    void match_contig_rr();
    
    // first step of create_extension: determine bp count and max represented bp at each position
    void extension_bp_count( vector<vector<int>> &ATCG, int start, int pos_mult );

    // second step of the create_extension process: count missed bp's per read, or in other words, the bp's represented below the max for that position
    void extension_missed_count( vector<vector<int>> &ATCG, vector<int> &missed_bp, int &missed_bp_tot, int start, int pos_mult );

    // third step in create_extension(): removal of reads that have errors over the threshold
    void extension_error_removal( vector<int> &missed_bp, int missed_bp_avg );

    // fourth step in create_extension: build extension string
    string extension_build_string( int start, int pos_mult, bool back );

    /// checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
    string create_extension( bool back );

    // extend() performs loops iterations of create_extension with length extend_len of each extension, at each iteration the extension is added to contig, and uses contig_sub_len characters from the front or back of the contig, which end is determined by the boolean value back
    void extend( bool back );

    // checks the coverage of matches at the given positions, returns the coverage
    long check_cov( long pos );

    // contig_fusion: Attempt to support fusion in case of possibly poorly constructed end
    bool check_fusion_support( string contig_f, int pos, bool back );
};
//////////////////////
// End Contig Class //
//////////////////////

#endif
