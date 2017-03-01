// $Author: benine $
// $Date$
// $Log$
// Contains the Process class for afin

#ifndef PROCESS_H
#define PROCESS_H

#include <unordered_map>
#include "log.hpp"
#include "readlist.hpp"
#include "fusion.hpp"
#include "contiglist.hpp"
#include "mismatch.hpp"
#include "contig.hpp"
#include <queue>
#include <mutex>
#include <condition_variable>
#include <sstream>

// contains hash of ranges in the readlist corresponding to keys made up of the first max_sort_char characters in the read
extern int max_search_loops;
extern int contig_sub_len;
extern int extend_len;
extern int max_sort_char;
extern int min_cov;
extern int min_overlap;
extern int max_threads;
extern int initial_trim;
extern int max_missed;
extern double stop_ext;
extern bool test_run;
extern int print_fused;
extern int screen_output;
extern int log_output;
extern int verbose;
extern int no_fusion;
extern double mismatch_threshold;
extern std::mutex log_mut;

//////////////////////////////
// Queue Template Functions //
//////////////////////////////

template <typename T>
class Queue{
  private:
    std::queue<T> qu;
    std::mutex mtx;
    std::condition_variable cv;

  public:
    // pop front item and return value
    T pop();

    // pass item by reference and return front of queue through the reference
    void pop( T& item );

    // push item onto back of queue
    void push( const T& item );
};

// pop front item and return value
template< typename T >
T Queue<T>::pop(){
  std::unique_lock<std::mutex> mlock( mtx );
  while( qu.empty() ){
    cv.wait(mlock);
  }

  auto item = qu.front();
  qu.pop();
  return item;
}

// pass item by reference and return front of queue through the reference
template< typename T >
void Queue<T>::pop( T& item ){
  std::unique_lock<std::mutex> mlock( mtx );
  while( qu.empty() ){
    cv.wait(mlock);
  }

  item = qu.front();
  qu.pop();
}

// push item onto back of queue
template< typename T >
void Queue<T>::push( const T& item ){
  std::unique_lock<std::mutex> mlock( mtx );
  qu.push( item );
  mlock.unlock();
  cv.notify_one();
}

///////////////
// End Queue //
///////////////

/////////////////////////////////////////////\
// Process Class: ////////////////////////////>
/////////////////////////////////////////////
//
// Process object, contains all classes, methods, data, and data references necessary for processing the contigs
// There will be only one Process object needed per iteration of this program
class Process{
  private:
    Readlist *reads;
    Contiglist *contigs;
    Fusion *fuse;

    // vectors for each of the iterable options
    std::vector<int> max_sort_char_iter;
    std::vector<int> contig_sub_len_iter;
    std::vector<int> extend_len_iter;
    std::vector<int> max_search_loops_iter;
    std::vector<int> min_cov_iter;
    std::vector<int> min_overlap_iter;
    std::vector<int> initial_trim_iter;
    std::vector<int> max_missed_iter;
    std::vector<double> stop_ext_iter;
    std::vector<double> mismatch_threshold_iter;

  public:
    std::string outfile;
    std::string readsfiles;
    std::string contigsfiles;
    int max_iterations;
    std::vector<int> extension_count;

    // iterable_opts is built of 9 elements containing the 9 different options that can be used as iterating options in afin
    std::unordered_map<std::string,std::string> iterable_opts;
    //vector<vector<double>*> iterables; // create vector of pointers to filled iterable vectors only

    Process();
    ~Process();

    // initalize logfile
    void logfile_init();

    // print options values to logfile
    void logfile_print_options();

    // Parse option int
    void parse_option( std::string opt_key, std::vector<int>* iter_vect );

    // Parse option double
    void parse_option( std::string opt_key, std::vector<double>* iter_vect );

    // Populate iterable options vector
    void populate_iterables();

    // set iterables global values at each iteration
    void set_iterables( int i );

    // Initializes data structures and turns over control to run_manager()
    void start_run();

    // Manages run
    void run_manager(int current_iteration);

    // closes logfile
    void close_log();

    void thread_worker( Queue<int>& q, unsigned int id);
};

///////////////////////
// End Process Class //
///////////////////////

#endif
