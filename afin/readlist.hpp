// class for holding and processing the readlist and the reverse complement reference list (rc_reflist)

#ifndef READLIST_HPP
#define READLIST_HPP

#include <unordered_map>
#include <tuple>
#include <vector>
#include "read.hpp"

class Readlist{
  private:
    static std::vector<std::string> rlist;
    static std::vector<long int> rc_reflist;
    std::unordered_map<std::string, std::tuple<long,long,long,long>> read_range;
    
  public:
    Readlist( std::string readsfiles );

    static bool cmp_rc( const int ind1, const int ind2 );
    static bool cmp_read( const std::string str1, const std::string str2 );
    template<class Iter, typename Order>
    void merge_sort( Iter first, Iter last, Order order );
    // sorts the reads by the first max_sort_char characters
    void sort_reads();

    // Produces list of references to the readlist sorted based on the first max_sort_char characters of the reverse_compliment of the
    //  referenced read
    //  Must be done after the readlist is sorted as it contains the locations in the readlist of the referenced read
    //  Uses an insertion sort to build the list
    void sort_rc();

    // creates a hash table within the process object that contains the ranges corresponding to equivalent first max_sort_char characters in the reads
    // This increases the efficiency of searching
    void create_read_range();

    // put reads from readfile into readlist
    void add_reads( std::string readsfiles );

    // checks if string is a homopolymer
    bool homopolymer_check( std::string seq );

    // return the read range corresponding to the first MAX_SORT characters of the read passed
    std::tuple<long,long,long,long> get_read_range( std::string read_seg );

    // return read from rlist with index ind
    std::string get_read( int ind );

    // return revcomp std::string of rlist member referenced by rc_reflist[ind]
    std::string get_rc_read( int ind );
};

#endif
