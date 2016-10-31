


#ifndef MATCH_HPP
#define MATCH_HPP

#include "readlist.hpp"
#include <string>
#include <vector>

class Match{
  private:
    Readlist *reads;
    std::vector<Read> matchlist;
    bool back;
    std::string sequence;
    
    // adds a read to the read list
    void push_match( std::string read, int pos );

    void push_match( std::string read, int pos, bool rev );

    // determines where the read passed matches the contig if at all for off the front matches
    void match_contig_fr();

    // determines where the read passed matches the contig if at all for off the back matches
    void match_contig_rr();  
  
  public:
    Match( Readlist *reads );

    Match( Readlist *reads, std::string sequence );

    // set sequence and back variables
    void set_seq( std::string sequence, bool back );

    // return matchlist length
    size_t get_matchlist_size();

    // return char in matchlist read at pos
    char get_pos( int read, int pos );
   
    // run the match for the current configuration
    void start_match();

    // checks the coverage of matches at the given positions, returns the coverage
    int check_cov( int pos );

    // remove match at index ind
    void remove_match( int ind );

    // clear match list
    void clearlist();
};

#endif
