// $Author: afinit $
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
double stop_ext;
bool test_run;
int print_fused;
int screen_output;
int log_output;
int verbose;
int no_fusion;
double mismatch_threshold;
mutex log_mut;

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
  iterables_len = 0;
  
  // initialize unordered_map for iterable options
  iterable_opts["max_search_loops"] = "10";
  iterable_opts["contig_sub_len"] = "100";
  iterable_opts["extend_len"] = "40";
  iterable_opts["max_sort_char"] = "4";
  iterable_opts["min_cov"] = "3";
  iterable_opts["min_overlap"] = "20";
  iterable_opts["initial_trim"] = "0";
  iterable_opts["max_missed"] = "5";
  iterable_opts["stop_ext"] = "0.5";
  iterable_opts["mismatch_threshold"] = "0.1";
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
  Log::Inst()->log_it( "stop_ext: " + to_string(stop_ext) );
}

// Parse option int
void Process::parse_option( string opt_key, int* glob_var, vector<int>* iter_vect ){
  string opt = iterable_opts[opt_key];
  
  try{
		// check for commas in string
		if( ! opt.find(',') ){
			*glob_var = stoi( opt );
		}
		else{
			// split along commas if found
			stringstream ss(opt);
			string item;
			while( getline( ss, item, ',' )){
				iter_vect->push_back(stoi(item));
			}
		}
	}
	catch( exception const & e ){
		Log::Inst()->log_it( "Error: " + string(e.what(), strlen(e.what())) + " : Invalid values: " + opt + " For option: " + opt_key );
		exit(0);
	}
}

// Parse option double
void Process::parse_option( string opt_key, double* glob_var, vector<double>* iter_vect ){
  string opt = iterable_opts[opt_key];

	try{ 
		// check for commas in string
		if( ! opt.find(',') ){
			*glob_var = stof( opt );
		}
		else{
			// split along commas if found
			stringstream ss(opt);
			string item;
			while( getline( ss, item, ',' )){
				iter_vect->push_back(stof(item));
			}
		}
	}
	catch( exception const & e ){
		Log::Inst()->log_it( "Error: " + string(e.what(), strlen(e.what())) + " : Invalid values: " + opt + " For option: " + opt_key );
		exit(0);
	}
}

// Populate iterable options vector
void Process::populate_iterables(){
  // process iterable options
  parse_option( "max_search_loops", &max_search_loops, &max_search_loops_iter );
  parse_option( "contig_sub_len", &contig_sub_len, &contig_sub_len_iter );
  parse_option( "extend_len", &extend_len, &extend_len_iter );
  parse_option( "max_sort_char", &max_sort_char, &max_sort_char_iter );
  parse_option( "min_cov", &min_cov, &min_cov_iter );
  parse_option( "min_overlap", &min_overlap, &min_overlap_iter );
  parse_option( "initial_trim", &initial_trim, &initial_trim_iter );
  parse_option( "max_missed", &max_missed, &max_missed_iter );
  parse_option( "stop_ext", &stop_ext, &stop_ext_iter );
  parse_option( "mismatch_threshold", &mismatch_threshold, &mismatch_threshold_iter );

	// get all iterables lengths
	vector<int> it_lens;
  it_lens.push_back(max_search_loops_iter.size());
  it_lens.push_back(contig_sub_len_iter.size());
	it_lens.push_back(extend_len_iter.size());
	it_lens.push_back(max_sort_char_iter.size());
	it_lens.push_back(min_cov_iter.size());
	it_lens.push_back(min_overlap_iter.size());
	it_lens.push_back(initial_trim_iter.size());
	it_lens.push_back(max_missed_iter.size());
	it_lens.push_back(stop_ext_iter.size());
	it_lens.push_back(mismatch_threshold_iter.size());

  // set iterables_len
  iterables_len = *max_element(begin(it_lens),end(it_lens));
  
  // check to make sure all options with more than one iteration are equal in iterations
  for( int i=0; i<it_lens.size(); i++  ){
  	if( it_lens[i] != iterables_len || it_lens[i] != 0 ){
  		Log::Inst()->log_it( "Error: Iterated options must have some number of iterations");
  		exit(0);
  	}
  }
  
  if( iterables_len == 0 ){
  	iterables_len = 1;
  } 
}

// set iterables global values at each iteration
void Process::set_iterables( int i ){
	if( !max_search_loops_iter.empty() ){
  	max_search_loops = max_search_loops_iter[i];
  }
  if( !contig_sub_len_iter.empty() ){
		contig_sub_len = contig_sub_len_iter[i];
  }
  if( !extend_len_iter.empty() ){
		extend_len = extend_len_iter[i];
  }
  if( !max_sort_char_iter.empty() ){
		max_sort_char = max_sort_char_iter[i];
  }
  if( !min_cov_iter.empty() ){
		min_cov = min_cov_iter[i];
  }
  if( !min_overlap_iter.empty() ){
		min_overlap = min_overlap_iter[i];
  }
  if( !initial_trim_iter.empty() ){
		initial_trim = initial_trim_iter[i];
  }
  if( !max_missed_iter.empty() ){
		max_missed = max_missed_iter[i];
  }
  if( !stop_ext_iter.empty() ){
		stop_ext = stop_ext_iter[i];
  }
  if( !mismatch_threshold_iter.empty() ){
		mismatch_threshold = mismatch_threshold_iter[i];
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
	
  // import reads once  
  reads = new Readlist( readsfiles ); 
	
  // loop through iterated options if present
	for( int i=0; i < iterables_len; i++ ){	
		// prevent printing of fused contigs
		if( no_fusion )
			print_fused = 0;

		// log output file
		Log::Inst()->log_it( string("output file: ") + outfile );
	
		// initialize objects
		contigs = new Contiglist( reads, contigsfiles, outfile );
		fuse = new Fusion( contigs, reads );

		Log::Inst()->log_it( "End initialization phase" );
	
		// make initial attempt to fuse contigs  
		if( ! no_fusion )
			fuse->run_fusion( true );
	
		if( test_run )
			contigs->output_contigs( 0, outfile + ".fus", "mid" );

		run_manager();

    // write contigs to fasta
    contigs->create_final_fasta(i);
	}
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
    if( ! no_fusion )
      fuse->run_fusion( false );
    
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
