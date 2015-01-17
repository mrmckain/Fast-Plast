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

using namespace std;

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
extern mutex log_mut;

//////////////////////////////
// Queue Template Functions //
//////////////////////////////

template <typename T>
class Queue{
  private:
    queue<T> qu;
    mutex mtx;
    condition_variable cv;
    
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
  unique_lock<mutex> mlock( mtx );
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
  unique_lock<mutex> mlock( mtx );
  while( qu.empty() ){
    cv.wait(mlock);
  }

  item = qu.front();
  qu.pop();
}

// push item onto back of queue
template< typename T >
void Queue<T>::push( const T& item ){
  unique_lock<mutex> mlock( mtx );
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
  
  public:
    string outfile;
    string readsfiles;
    string contigsfiles;
    
    Process();
    ~Process();

    // initalize logfile
    void logfile_init();

    // Initializes data structures and turns over control to run_manager()
    void start_run();
    
    // Manages run 
    void run_manager();

    // closes logfile
    void close_log();

    void thread_worker( Queue<int>& q, unsigned int id);
};

///////////////////////
// End Process Class //
///////////////////////

#endif
