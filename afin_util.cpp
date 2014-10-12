// Worker function

#include <iostream>
#include "afin_util.hpp"
#include "process.hpp"
#include "contig.hpp"
#include "queue.tcc"

// Consume function which will act as the driver for an individual thread
void thread_worker(vector<Contig>& contigs, Queue<int>& q, unsigned int id) {
  for (;;) {
    auto item = q.pop();
    if( item == -1 ){
      break;
    }
    else{
      contigs[item].extend( false );
      contigs[item].extend( true );
    }
  }
}

// tally mismatches in substrings passed and return score in the form of misatches per length
double mismatch_score( string contig_sub_a, string contig_sub_b ){
  int mismatch = 0;

  // protect against division by zero
  if( contig_sub_a.length() == 0 ){
    return 1.0;
  }

  for( int i=0; i<contig_sub_a.length(); i++ ){
    if( contig_sub_a[i] != contig_sub_b[i] ){
      mismatch++;
    }
    // bail when too many mismatches have been tallied.. avoid excessive processing
    if( double(mismatch) /contig_sub_a.length() > mismatch_threshold ){
      return 1.0;
    }
  }

  double score = double(mismatch) / contig_sub_a.length();
  return score;
}

// checks if string is a homopolymer
bool homopolymer_check( string seq ){
  int i = 0;
  if( seq[0] == 'A' ){
    for( i=1; i<seq.length(); i++ ){
      if( seq[i] != 'A' ){
        break;
      }
    }
  }
  else if( seq[0] == 'T' ){
    for( i=1; i<seq.length(); i++ ){
      if( seq[i] != 'T' ){
        break;
      }
    }
  }
  else if( seq[0] == 'C' ){
    for( i=1; i<seq.length(); i++ ){
      if( seq[i] != 'C' ){
        break;
      }
    }
  }
  else if( seq[0] == 'G' ){
    for( i=1; i<seq.length(); i++ ){
      if( seq[i] != 'G' ){
        break;
      }
    }
  }
  if( i == seq.length() ){
    return true;
  }

  return false;
}

// PRINT USAGE FUNCTION
void print_usage( string prog ){
  cout << "Usage: " << prog << " -c contigfile(s) -r readfile(s) [-o outfile] [-m max_sort_char] [-s contig_sub_len]" << endl;
  cout << "          [-l max_search_loops] [-i min_cov_init] [-p min_overlap] [-t max_threads]" << endl;
  cout << "       " << prog << " -h" << endl;
  cout << endl;
  cout << "  -c contigfile          Comma separated list of files containing contigs" << endl;
  cout << "  -r readfile            Comma separated list of files containing reads" << endl;
  cout << "  -o outfile             Output will be printed to the outfile specified with a .fasta extension" << endl;
  cout << "  -m sort_char           [default:   4] Sorts the reads by the first max_sort_char characters" << endl;
  cout << "  -s sub_len             [default: 100] Will focus on the current last contig_sub_len characters of the contig in each search" << endl;
  cout << "  -l search_loops        [default:  10] Will search against each contig a maximum of max_search_loops times before comparing them" << endl;
  cout << "  -i min_cov             [default:   3] Will stop adding bp's once the coverage falls below min_cov_init" << endl;
  cout << "  -p min_overlap         [default:  20] Only those reads overlapping the contig by at least min_overlap bp's will be returned in each search" << endl;
  cout << "  -t max_threads         [default:   4] Will only run max_threads threads at a time" << endl;
  cout << "  -d initial_trim        [default:  20] Length to trim off the beginning and end of each contig at the start of the program" << endl;
  cout << "  -e max_missed          [default:   5] Maximum allowable mismatched bp's for each read" << endl;
  cout << "  -g mismatch_threshold  [default:  .1] maximum percentage of mismatches allowed when fusing two contigs" << endl;
  cout << "  -x extend_len          [default:  40] Will add a max of 80 bp's each search loop" << endl << endl;
}

