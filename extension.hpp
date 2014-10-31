

#ifndef EXTENSION_HPP
#define EXTENSION_HPP

#include "match.hpp"
#include <vector>

class Extension{
  public:
    Match matches;
    Readlist *reads;
    std::vector<std::vector<int>> ATCG;
    std::vector<int> missed_bp;
    int start;
    bool back;
    std::string contig;
    std::string exten_seq;
    int pos_mult;
    int len;
    int missed_bp_tot;
    int missed_bp_avg;
    
    Extension( Readlist *reads, int len );
    Extension( Readlist *reads, int len, std::string contig );

    // set the value of the missed_bp std::vector
    void set_missed_bp( std::vector<int> missed_bp );

    // set the value of missed_bp_avg
    void set_missed_avg( int missed_avg );

    // first step of get_extension: determine bp count and max represented bp at each position
    void bp_count();

    // second step of the get_extension process: count missed bp's per read, or in other words, the bp's represented below the max for that position
    void missed_count();

    // third step in get_extension(): removal of reads that have errors over the threshold
    bool error_removal();

    // fourth step in get_extension: build extension string
    void build_string();

    // checks the matches against each other and the contig, compiles an extension of length len (or less if the length is limited by matches) that is returned 
    std::string get_extension( std::string contig, bool back );

    // simply returns the built extension string
    std::string get_extension();
};

#endif
