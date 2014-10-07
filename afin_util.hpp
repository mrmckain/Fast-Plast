// contains utility functions for afin

#ifndef AFIN_UTIL_H
#define AFIN_UTIL_H

#include "queue.tcc"
#include "contig.hpp"

using namespace std;

void thread_worker(vector<Contig>& contigs, Queue<int>& q, unsigned int id);

// tally mismatches in substrings passed and return score in the form of misatches per length
double mismatch_score( string contig_sub1, string contig_sub2 );

// checks if string is a homopolymer
bool homopolymer_check( string seq );

void print_usage( string prog );

#endif
