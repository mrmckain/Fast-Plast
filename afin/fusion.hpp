

#ifndef FUSION_HPP
#define FUSION_HPP

#include "mismatch.hpp"
#include "contiglist.hpp"
#include <vector>
#include <string>

class Fusion{
  private:
    Contiglist *contigs;
    Readlist *reads;
    std::vector<Mismatch> match_list;
    std::vector<int> contig_remove_list;

  public:
    // include variable for counting fusions in each loop
    int fusions_completed;
    
    Fusion( Contiglist *contigs, Readlist *reads );

    // creates id of fused contigs
    std::string get_fused_id( std::string contig1_id, std::string contig2_id );

    // complete contig_fusion process
    void contig_fusion_log( Mismatch fus );

    // complete contig_fusion process
    void commit_fusion( std::string fused, std::string fused_id, int index_i, int index_j );

    // check overlap section for mismatches
    Mismatch overlap_check( std::string contig_a, std::string contig_b, int overlap, int end_i, int end_j );

    // create fused contig std::string
    std::string build_fusion_string( std::string contig_a, std::string contig_b, int overlap );

    // remove duplicates from contig remove list
    void dedup_removals();

    // sort index list for removing contigs
    void sort_removals();

    // remove fused contigs from contigs list
    void process_removals();

    // compile list of best mismatch scores between contigs that meet the mismatch threshold
    std::vector<Mismatch> get_mismatch_scores( bool first_run );

    // sort the match_list for easier
    void sort_matches();

    // cleans match_list from conflicting matches
    void clean_matches();

    // process the compiled list of fusions
    void process_fusions();

    // fuse contigs wherever possible
    void run_fusion( bool first_run );

    // contig_fusion: Attempt to support fusion in case of possibly poorly constructed end.. returns new score from section in question
    //    ::> contig object is the second while contig_ref is the first and the extension is being made off the front of the object
    double check_fusion_support( std::string contig, std::string contig_ref );

    // tally mismatches in substrings passed and return score in the form of misatches per length
    double mismatch_score( std::string contig_sub_a, std::string contig_sub_b );
};

#endif
