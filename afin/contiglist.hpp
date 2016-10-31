

#ifndef CONTIGLIST_HPP
#define CONTIGLIST_HPP

#include "contig.hpp"
#include <vector>
#include <string>

class Contiglist{
  private:
    Readlist *reads;
    std::vector<Contig> contigs;
    std::vector<Contig> contigs_fused;
    std::string outfile;

  public:
    Contiglist( Readlist *reads, std::string contigsfiles, std::string outfile );

    // return contig at index ind
    Contig get_contig( int ind );

    // return contig at index ind
    Contig *get_contig_ref( int ind );

    // return conitg list size
    int get_list_size();

    // remove contig at index ind
    void remove_contig( int ind );

    // append contig to contigs
    void append_contig( int list_num, Contig cont );

    // cycles through each contig and parses out the first section of the id
    void parse_ids();

    // put contigs from contfile into contlist
    void add_contigs( std::string contigsfiles );

    // print contigs
    void output_contigs( int list_num, std::string file, std::string id_suffix );

    // prints results to fasta file with outfile prefix and additional information is printed to a text based file with outfile prefix
    void create_final_fasta( int current_iteration );
};

#endif
