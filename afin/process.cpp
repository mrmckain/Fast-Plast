// $Author: afinit $
// $Date$
// $Log$
// Contains the Process class for afin

#include "process.hpp"
#include <cstring>
#include <thread>

int max_search_loops;
int contig_sub_len;
int extend_len;
int max_sort_char;
int min_cov;
int min_overlap;
int max_threads;
int initial_trim;
int max_missed;
double stop_ext;
bool test_run;
int print_fused;
int screen_output;
int log_output;
int verbose;
int no_fusion;
double mismatch_threshold;
std::mutex log_mut;

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
  max_iterations = 1;

  // initialize unordered_map for iterable options
  iterable_opts["max_search_loops"] = "10";
  iterable_opts["contig_sub_len"] = "100";
  iterable_opts["extend_len"] = "40";
  iterable_opts["max_sort_char"] = "4";
  iterable_opts["min_cov"] = "3";
  iterable_opts["min_overlap"] = "20";
  iterable_opts["max_missed"] = "5";
  iterable_opts["stop_ext"] = "0.5";
  iterable_opts["mismatch_threshold"] = "0.1";

  // Declare vector variables for storing iterable_opts
  std::vector<int> max_search_loops_iter;
  std::vector<int> contig_sub_len_iter;
	std::vector<int> extend_len_iter;
	std::vector<int> max_sort_char_iter;
	std::vector<int> min_cov_iter;
	std::vector<int> min_overlap_iter;
	std::vector<int> max_missed_iter;
	std::vector<double> stop_ext_iter;
	std::vector<double> mismatch_threshold_iter;
}

Process::~Process(){
  delete fuse;
  delete reads;
  delete contigs;
}

// print options values to logfile
void Process::logfile_print_options(){
  // NOTE:: add capability to print iterable logfile options
  // output starting option values
  Log::Inst()->log_it( "OPTION VALUES" );
  Log::Inst()->log_it( "contig_sub_len: " + std::to_string(contig_sub_len) );
  Log::Inst()->log_it( "extend_len: " + std::to_string(extend_len) );
  Log::Inst()->log_it( "max_search_loops: " + std::to_string(max_search_loops) );
  Log::Inst()->log_it( "max_sort_char: " + std::to_string(max_sort_char) );
  Log::Inst()->log_it( "min_cov: " + std::to_string(min_cov) );
  Log::Inst()->log_it( "min_overlap: " + std::to_string(min_overlap) );
  Log::Inst()->log_it( "initial_trim: " + std::to_string(initial_trim) );
  Log::Inst()->log_it( "max_missed: " + std::to_string(max_missed) );
  Log::Inst()->log_it( "mismatch_threshold: " + std::to_string(mismatch_threshold) );
  Log::Inst()->log_it( "max_threads: " + std::to_string(max_threads) );
  Log::Inst()->log_it( "stop_ext: " + std::to_string(stop_ext) );
}

// initalize logfile
void Process::logfile_init(){
  Log::Inst()->open_log( outfile + ".log" );
}

// Parse option int
void Process::parse_option( std::string opt_key, std::vector<int>* iter_vect ){
  std::string opt = iterable_opts[opt_key];

  try{
		// check for commas in string
		if( ! opt.find(',') ){
      iter_vect->push_back(std::stoi(opt));
    }
		else{
			// split along commas if found
			std::stringstream ss(opt);
			std::string item;
			while( getline( ss, item, ',' )){
				iter_vect->push_back(std::stoi(item));
			}
		}
	}
	catch( std::exception const & e ){
		Log::Inst()->log_it( "Error: " + std::string(e.what(), strlen(e.what())) + " : Invalid values: " + opt + " For option: " + opt_key );
		exit(0);
	}

  // increase max_iterations if this option is longer than the current max_iterations
  if( iter_vect->size() > max_iterations )
    max_iterations = iter_vect->size();
}

// Parse option double
void Process::parse_option( std::string opt_key, std::vector<double>* iter_vect ){
  std::string opt = iterable_opts[opt_key];

	try{
		// check for commas in string
		if( ! opt.find(',') ){
			iter_vect->push_back(std::stof(opt));
		}
		else{
			// split along commas if found
			std::stringstream ss(opt);
			std::string item;
			while( getline( ss, item, ',' )){
				iter_vect->push_back(std::stof(item));
			}
		}
	}
	catch( std::exception const & e ){
		Log::Inst()->log_it( "Error: " + std::string(e.what(), strlen(e.what())) + " : Invalid values: " + opt + " For option: " + opt_key );
		exit(0);
	}

  // increase max_iterations if this option is longer than the current max_iterations
  if( iter_vect->size() > max_iterations )
    max_iterations = iter_vect->size();
}

// Populate iterable options vector
void Process::populate_iterables(){
  // process iterable options
  parse_option( "max_search_loops", &max_search_loops_iter );
  parse_option( "contig_sub_len", &contig_sub_len_iter );
  parse_option( "extend_len", &extend_len_iter );
  parse_option( "max_sort_char", &max_sort_char_iter );
  parse_option( "min_cov", &min_cov_iter );
  parse_option( "min_overlap", &min_overlap_iter );
  parse_option( "max_missed", &max_missed_iter );
  parse_option( "stop_ext", &stop_ext_iter );
  parse_option( "mismatch_threshold", &mismatch_threshold_iter );
}

// set iterables global values at each iteration
void Process::set_iterables( int i ){
  // max_search_loops_iter
  if( max_search_loops_iter.size() > i ){
  	max_search_loops = max_search_loops_iter[i];
  }
  else{
    max_search_loops = max_search_loops_iter.back();
  }

  // contig_sub_len_iter
  if( contig_sub_len_iter.size() > i ){
		contig_sub_len = contig_sub_len_iter[i];
  }
  else{
    contig_sub_len = contig_sub_len_iter.back();
  }

  // extend_len_iter
  if( extend_len_iter.size() > i ){
		extend_len = extend_len_iter[i];
  }
  else{
    extend_len = extend_len_iter.back();
  }

  // max_sort_char_iter
  if( max_sort_char_iter.size() > i ){
		max_sort_char = max_sort_char_iter[i];
  }
  else{
    max_sort_char = max_sort_char_iter.back();
  }

  // min_cov_iter
  if( min_cov_iter.size() > i ){
		min_cov = min_cov_iter[i];
  }
  else{
    min_cov = min_cov_iter.back();
  }

  // min_overlap_iter
  if( min_overlap_iter.size() > i ){
		min_overlap = min_overlap_iter[i];
  }
  else{
    min_overlap = min_overlap_iter.back();
  }

  // max_missed_iter
  if( max_missed_iter.size() > i ){
		max_missed = max_missed_iter[i];
  }
  else{
    max_missed = max_missed_iter.back();
  }

  // stop_ext_iter
  if( stop_ext_iter.size() > i ){
		stop_ext = stop_ext_iter[i];
  }
  else{
    stop_ext = stop_ext_iter.back();
  }

  // mismatch_threshold_iter
  if( mismatch_threshold_iter.size() > i ){
		mismatch_threshold = mismatch_threshold_iter[i];
  }
  else{
    mismatch_threshold = mismatch_threshold_iter.back();
  }
}

// Initializes data structures and turns over control to run_manager()
void Process::start_run(){
  // initialize logfile if logging enabled
  if( log_output || screen_output ){
    logfile_init();
  }

	// process iterable options
	populate_iterables();
  Log::Inst()->log_it( std::string("Max Iterations of afin: ") + std::to_string(max_iterations) );
  set_iterables(0);

  // import reads once
  reads = new Readlist( readsfiles );

  // initialize objects (read in contigs here)
  contigs = new Contiglist( reads, contigsfiles, outfile );
  fuse = new Fusion( contigs, reads );

  // loop through iterated options if present
	for( int i=0; i < max_iterations; i++ ){
    Log::Inst()->log_it( std::string("Begin Iteration of afin: ") + std::to_string(i+1) );

		// prevent printing of fused contigs
		if( no_fusion )
			print_fused = 0;

    // Set options for this iteration
    set_iterables(i);
    // print current options for current iteration
    logfile_print_options();

		// log output file
		Log::Inst()->log_it( std::string("output file: ") + outfile );

		Log::Inst()->log_it( "End initialization phase" );

		// make initial attempt to fuse contigs
		if( ! no_fusion )
			fuse->run_fusion( true );

		if( test_run )
			contigs->output_contigs( 0, outfile + ".fus", "mid" );

		run_manager(i);

    // write contigs to fasta
    contigs->create_final_fasta(i);

    // NOTE:: REMOVED TO allow single contigs to be further extended.. will remove this code later
    // // if there is only one contig left, no need to continue
    // if( contigs->get_list_size() == 1 ){
    //   Log::Inst()->log_it( "Assembled into 1 contig. Exiting" );
    //   break;
    // }
	}

  // close the log file dude
  if( log_output || screen_output ){
    Log::Inst()->close_log();
  }
}

// Manages run
void Process::run_manager(int current_iteration){
  // create thread array with max_thread entries
  std::vector<std::thread> t;
  Queue<int> qu;
  int extension_sum = 0;
  fuse->fusions_completed = 0;

  // loop max search loops
  for( int j=0; j<max_search_loops; j++ ){
    Log::Inst()->log_it( "Begin Extensions" );

    // initialize threads
    for( int i=0; i<max_threads; i++ ){
      t.push_back(std::thread( &Process::thread_worker, this, std::ref(qu), i ));
      extension_count.push_back(0);
    }

    if( test_run ){
      Log::Inst()->log_it( "contigs.size(): " + std::to_string(contigs->get_list_size()) + " max_threads: " + std::to_string(max_threads) );
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
      extension_sum += extension_count[i];
    }

    // remove threads from vector
    t.clear();
    extension_count.clear();

    // removed for master branch until algorithm can be adjusted
    if( ! no_fusion )
      fuse->run_fusion( false );

    // for test_runs, writes intermediate contigs to file to help troubleshoot
    if( test_run ){
      contigs->output_contigs( 0, outfile + ".fus" + std::to_string(j) + ".iter" + std::to_string(current_iteration), "mid" );
    }

    // if no extensions or fusions were made this round, no need to continue
    if( extension_sum == 0 && fuse->fusions_completed == 0 ){
      Log::Inst()->log_it( "Nothing to do here:" );
      Log::Inst()->log_it( std::string("  contigs:           ") + std::to_string(contigs->get_list_size()) );
      Log::Inst()->log_it( std::string("  extension_sum:     ") + std::to_string(extension_sum) );
      Log::Inst()->log_it( std::string("  fusions_completed: ") + std::to_string(fuse->fusions_completed) );
      return;
    }

    // reset extension_sum and fusions_completed
    extension_sum = 0;
    fuse->fusions_completed = 0;
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
      extension_count[id] += contigs->get_contig_ref(item)->extend( false );
      extension_count[id] += contigs->get_contig_ref(item)->extend( true );
    }
  }
}

//////////////////////////
// END PROCESS ///////////
//////////////////////////
