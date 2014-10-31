// $Author: benine $
// $Date$
// $Log$
// Contains the Process class for afin

#include "process.hpp"
#include <thread>

using namespace std;

int max_search_loops;
int contig_sub_len;
int extend_len;
int max_sort_char;
int min_cov;
int min_overlap;
int max_threads;
int initial_trim;
int max_missed;
bool test_run;
int print_fused;
int screen_output;
int log_output;
int verbose;
double mismatch_threshold;


////////////////////////////////////
//////// PROCESS DEFINITIONS ///////
////////////////////////////////////

// Process constructor
Process::Process(){
  outfile = "afin_out";
  readsfiles = "";
  contigsfiles = "";
  reads = 0;
  contigs = 0;
  fuse = 0;
}

Process::~Process(){
  delete fuse;
  delete reads;
  delete contigs;
}

// initalize logfile
void Process::logfile_init(){
  Log::Inst()->open_log( outfile + ".log" );

  // output starting option values
  Log::Inst()->log_it( "OPTION VALUES" );
  Log::Inst()->log_it( "contig_sub_len: " + to_string(contig_sub_len) );
  Log::Inst()->log_it( "extend_len: " + to_string(extend_len) );
  Log::Inst()->log_it( "max_search_loops: " + to_string(max_search_loops) );
  Log::Inst()->log_it( "max_sort_char: " + to_string(max_sort_char) );
  Log::Inst()->log_it( "min_cov: " + to_string(min_cov) );
  Log::Inst()->log_it( "min_overlap: " + to_string(min_overlap) );
  Log::Inst()->log_it( "initial_trim: " + to_string(initial_trim) );
  Log::Inst()->log_it( "max_missed: " + to_string(max_missed) );
  Log::Inst()->log_it( "mismatch_threshold: " + to_string(mismatch_threshold) );
  Log::Inst()->log_it( "max_threads: " + to_string(max_threads) );
}

// Initializes data structures and turns over control to run_manager()
void Process::start_run(){
  // initialize logfile if logging enabled
  if( log_output || screen_output ){
    logfile_init();
  }

  // log output file
  Log::Inst()->log_it( string("output file: ") + outfile );
  
  // initialize objects
  reads = new Readlist( readsfiles ); 
  contigs = new Contiglist( reads, contigsfiles, outfile );
  fuse = new Fusion( contigs, reads );

  Log::Inst()->log_it( "End initialization phase" );
  
  // make initial attempt to fuse contigs  
  fuse->run_fusion();
  
  if( test_run )
    contigs->output_contigs( 0, outfile + ".fus", "mid" );

  run_manager();

  contigs->create_final_fasta();
}

// Manages run 
void Process::run_manager(){
  // create thread array with max_thread entries
  vector<thread> t;
  Queue<int> qu;

  // loop max search loops
  for( int j=0; j<max_search_loops; j++ ){
    Log::Inst()->log_it( "Begin Extensions" );
      
    // initialize threads
    for( int i=0; i<max_threads; i++ ){
      t.push_back(thread( &Process::thread_worker, this, ref(qu), i ));
    }

    if( test_run ){
      Log::Inst()->log_it( "contigs.size(): " + to_string(contigs->get_list_size()) + " max_threads: " + to_string(max_threads) );
    }

    // push each thread onto queue
    for( int i=0; i<contigs->get_list_size(); i++ ){
      qu.push( i );
    }

    // push stop signals onto queue for each thread
    for( int i=0; i<max_threads; i++ ){
      qu.push( -1 );
    }

    // join threads
    for( int i=0; i<max_threads; i++ ){
      t[i].join();
    }

    // remove threads from vector
    t.erase( t.begin(), t.begin()+max_threads );

    // removed for master branch until algorithm can be adjusted
    fuse->run_fusion();
    
    if( test_run ){
      contigs->output_contigs( 0, outfile + ".fus" + to_string(j), "mid" );
    }
  }
}

// Consume function which will act as the driver for an individual thread
void Process::thread_worker( Queue<int>& q, unsigned int id) {
  for (;;) {
    auto item = q.pop();
    if( item == -1 ){
      break;
    }
    else{
      contigs->get_contig_ref(item)->extend( false );
      contigs->get_contig_ref(item)->extend( true );
    }
  }
}

//////////////////////////
// END PROCESS ///////////
//////////////////////////
