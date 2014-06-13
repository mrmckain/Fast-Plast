// contains utility functions for afin

#ifndef AFIN_UTIL_H
#define AFIN_UTIL_H

#include "queue.tcc"
#include "contig.hpp"

using namespace std;

void thread_worker(vector<Contig>& contigs, Queue<int>& q, unsigned int id);

void print_usage( string prog );

#endif
