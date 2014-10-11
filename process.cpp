// $Author: benine $
// $Date$
// $Log$
// Contains the Process class for afin

#include <cstdio>
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unistd.h>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include <utility>
#include <thread>
#include "print_time.hpp"
#include "revcomp.hpp"
#include "contig.hpp"
#include "process.hpp"
#include "afin_util.hpp"
#include "mismatch.hpp"

using namespace std;

vector<string> readlist;
vector<long int> rc_reflist;
unordered_map<string, tuple<long,long,long,long>> read_range;
int max_search_loops;
int contig_sub_len;
int extend_len;
int max_sort_char;
int min_cov_init;
int min_overlap;
int max_threads;
int initial_trim;
int max_missed;
bool test_run;
double mismatch_threshold;

////////////////////////////////////
//////// PROCESS DEFINITIONS ///////
////////////////////////////////////

// compare function for sorting the reverse compliment reference list
bool cmp_rc( const int ind1, const int ind2 ){
  string str1 = readlist[ind1].substr( readlist[ind1].length()-max_sort_char, max_sort_char );
  str1 = revcomp( str1 );
  string str2 = readlist[ind2].substr( readlist[ind2].length()-max_sort_char, max_sort_char );
  str2 = revcomp( str2 );

  return (str1.compare( str2 ) < 0 ); 
}

// compare function for sorting the initial read list
bool cmp_read( const string str1, const string str2 ){
  return (str1.compare( 0, max_sort_char, str2.substr( 0, max_sort_char ) ) < 0 ); 
}

// mergesort template to be used with any datatype
template<class Iter, typename Order>
void merge_sort( Iter first, Iter last, Order order ){
  if (last - first > 1){
    Iter middle = first + (last - first)/2;
    merge_sort( first, middle, order );
    merge_sort( middle, last, order );
    inplace_merge( first, middle, last, order );
  }
}

// Process constructor
Process::Process(){
  outfile = "afin_out";
  readsfiles = "";
  contigsfiles = "";
}

// uses a mergesort to sort the read list based on the first max_sort_char characters of each read
void Process::sort_reads(){
  merge_sort( readlist.begin(), readlist.end(), cmp_read );
}

// Produces list of references to the readlist sorted based on the first max_sort_char characters of the reverse_compliment of the
//  referenced read
//  Must be done after the readlist is sorted as it contains the locations in the readlist of the referenced read
//  Uses a mergesort
void Process::sort_rc(){
  rc_reflist.resize( readlist.size() );

  for( int i=0; i<rc_reflist.size(); i++ ){
    rc_reflist[i]=i;
  }

  merge_sort( rc_reflist.begin(), rc_reflist.end(), cmp_rc );
}

// creates a hash table within the process object that contains the ranges corresponding to equivalent first max_sort_char characters in the reads
// This increases the efficiency of searching
void Process::create_read_range(){
  string current( readlist[0].substr( 0, max_sort_char ) );
  long curr_start = 0;
  // cycle through the readlist
  for( long int i=1; i<readlist.size(); i++ ){
    if( readlist[i].compare( 0, max_sort_char, current ) != 0 ){
      // insert values into hash table, the ordered pair reflecting the range is increased by 1 to differentiate between a false search  
      read_range.insert( { current, make_tuple( curr_start+1, i, 0, 0 )});
      curr_start = i;
      current = readlist[i].substr( 0, max_sort_char );
    }
  }
  
  // add the last entry
  read_range.insert( { current, make_tuple( curr_start+1, readlist.size(), 0, 0 )});

  string rc( readlist[0].substr( readlist[0].length() - max_sort_char, max_sort_char ) );
  current = revcomp( rc );
  curr_start = 0;

  // cycle through the rc_reflist to include these ranges in the hash as well
  for( long int i=1; i<rc_reflist.size(); i++ ){
    rc = readlist[i].substr( readlist[i].length() - max_sort_char, max_sort_char );
    rc = revcomp( rc );
    if( rc.compare( current ) != 0 ){
      // check if current already exists in hash
      if( get<0>( get_read_range(current) ) == -1 ){
        read_range.insert( { current, make_tuple( 0, 0, curr_start+1, i )});
      }
      else{
        get<2>( read_range[current] ) = curr_start+1;
        get<3>( read_range[current] ) = i;
      }
      curr_start = i;
      current = rc;
    }
  }
  
  // add the last entry
  if( get<0>( get_read_range(current) ) == -1 ){
    read_range.insert( { current, make_tuple( 0, 0, curr_start+1, rc_reflist.size() )});
  }
  else{
    get<2>( read_range[current] ) = curr_start+1;
    get<3>( read_range[current] ) = rc_reflist.size();
  }
}

// put reads from readfile into readlist
void Process::add_reads(){
  stringstream ss;
  ss.str( readsfiles );
  string filename;
  string buffer("");
  string line("");
  int line_count = 0;
    
  while( getline( ss, filename, ',' )){
    cout << "readfile: " << filename << endl;
    cout << "add_reads time: ";
    print_time();

    // open read file
    ifstream read( filename );
   
    // check what type of file it is
    if( getline( read, line ) ){
      if( line[0] == '@' ){
        // return to beginning of file
        read.seekg( 0, ios::beg );

        // read in fastq reads
        while( getline( read, line )){
          line_count++;
          if( line[0] == '@' ){
            if( getline( read, line )){
              line_count++;
              if( !homopolymer_check( line ) ){
                readlist.push_back( line );
              }
              if( getline( read, line )){
                line_count++;
                if( line[0] == '+'){
                  if( getline( read, line )){
                    line_count++;
                    continue;
                  }
                  else{
                    fprintf( stderr, "Error reading fastq file. Line missing. Line: %d\n", line_count );
                  }
                }
                else{
                  fprintf( stderr, "Error reading fastq file. '+' expected at this line. Line: %d\n", line_count );
                }
              }
              else{
                fprintf( stderr, "Error reading fastq file. Line missing. Line: %d\n", line_count );
              }
            }
            else{
              fprintf( stderr, "Error reading fastq file. Line missing. Line: %d\n", line_count );
            }
          }
          else{  
            fprintf( stderr, "Error reading fastq file. '@' expected at the beginning of this line. Line: %d\n", line_count );
          }
        }
      }
      else if( line[0] == '>' ){
        // return to beginning of file
        read.seekg( 0, ios::beg );
        
        // read in reads to vector from fasta file
        while( getline( read, line ) ){
          if( line[0] == '>' && buffer.length() != 0 ){
            //cout << buffer << endl;
            if( !homopolymer_check( buffer ) ){
              readlist.push_back( buffer );
            }
            buffer = "";
          }
          else if ( line[0] == '>' ) {
          }
          else{
            buffer += line;
          }
        }
      }
      else {
        fprintf( stderr, "Error: Unexpected file type. Needs to be fasta or fastq file for input." );
      }
    }

    // close read file
    read.close();
    print_time();
  }

  // insert last line into readlist
  if( buffer.length() != 0 ){
    if( !homopolymer_check( buffer ) ){
      readlist.push_back( buffer );
    }
  }
}

// parses the cov value from the contig_id and passes the result back as a double
// If cov is not in the header, return 1
double Process::parse_cov( string contig_id ){
  double cov = 1;
  size_t pos = contig_id.find( "cov_" ) + 4;
  if( pos != string::npos ){
    string contig_cov_str = contig_id.substr( pos, contig_id.length() - pos );
    pos = contig_cov_str.find( "_" );
    if( pos != string::npos ){
      contig_cov_str = contig_cov_str.substr( 0, pos );
    }
    cov = atof( contig_cov_str.c_str() );
  }

  return cov;
}

// check coverage of each contig, calculate the average coverage, then remove into a separate data structure any contigs that have more than 2xAvg coverage
void Process::contig_cov(){
  double cov_total = 0;
  int total_contigs = contigs.size();
 
  // protect against division by 0 and the end of the world
  if( total_contigs == 0 ){
    return;
  }

  for( int i=0; i<total_contigs; i++ ){
    // get contig id, parse out cov_##, add to the total of all cov values
    double cov = contigs[i].get_cov();
    cov_total += cov;
  }

  // find average of all cov values
  double cov_avg = cov_total / total_contigs;
  
  for( int i=0; i<total_contigs; i++ ){
    // get contig_id, parse out cov_## and compare this value to the avg
    double cov = contigs[i].get_cov();

    if( cov > cov_avg * 2.0 ){
      // set the value of doub_cov to false to indicate ignoring when extending contigs 
      contigs[i].set_doub_cov( true );
    }
  }
}

// cycles through each contig and parses out the first section of the id
void Process::parse_ids(){
  for( int i=0; i<contigs.size(); i++ ){
    string contig_id = contigs[i].get_contig_id(); 
    size_t pos = contig_id.find( "_", 5, 1 );

    if( pos != string::npos ){
      cout << "contig_id premod: " << contig_id << endl;
      contig_id = contig_id.substr( 0, pos );
      cout << "contig_id postmod: " << contig_id << endl;
    }
  }
}

// put contigs from contfile into contlist
void Process::add_contigs(){
  stringstream ss;
  ss.str( contigsfiles );
  string filename;
  string buffer("");
  string line("");
  string contig_id("");
    
  while( getline( ss, filename, ',' )){
    cout << "contigfile: " << filename << endl;
    cout << "add_contigs time: ";
    print_time();
  
    // open contig file
    ifstream cont( filename );

    // read in contig objects
    while( getline( cont, line ) ){
      if( line[0] == '>' && buffer.length() != 0 ){
        if( buffer.length() > 2*initial_trim + contig_sub_len ){
          buffer = buffer.substr( initial_trim, buffer.length() - 2*initial_trim );
          contigs.push_back( Contig( buffer, contig_id, parse_cov( contig_id ), min_cov_init ));
        }
        else if( buffer.length() > contig_sub_len ){
          int trim = (buffer.length() - contig_sub_len) / 2;
          buffer = buffer.substr( trim, buffer.length() - 2*trim );
          contigs.push_back( Contig( buffer, contig_id, parse_cov( contig_id ), min_cov_init ));
        }

        buffer = "";
        contig_id = line.substr(1);
      }
      else if ( line[0] == '>' ){
        contig_id = line.substr(1);
      }
      else{
        buffer += line;
      }
    }
    
    // close contig file
    cont.close();
 
  }

  // insert last line into contigs list
  if( buffer.length() != 0 ){
    contigs.push_back( Contig( buffer, contig_id, parse_cov( contig_id ), min_cov_init ) );
  }

  contig_cov();
  parse_ids();
}

// return contig from contigs_fused with index contig_ind
string Process::get_contig_fused( int contig_ind ){
  return contigs_fused[ contig_ind ].get_contig();
}

// return contig with index contig_ind
string Process::get_contig( int contig_ind ){
  return contigs[ contig_ind ].get_contig();
}

// print contigs
void Process::print_contigs_to_file( string file, string id_suffix ){
  // open outfile
  ofstream outfile_fp( file+".fasta");

  // print out each line to the 
  for( int i=0; i<contigs.size(); i++ ){
    outfile_fp << ">" << contigs[i].get_contig_id() << "_" << id_suffix << endl;
    outfile_fp << get_contig(i) << endl;
  }

  outfile_fp.close();
}

// prints results to fasta file with outfile prefix and additional information is printed to a text based file with outfile prefix
void Process::print_to_outfile(){
  // remove directories from outfile to form id_suffix if necessary
  size_t id_suffix_pos = outfile.find_last_of( "/" );
  string id_suffix = outfile;
  if( id_suffix_pos != string::npos ){
    id_suffix = id_suffix.substr( id_suffix_pos + 1 );
  }

  // print completed contigs to file
  print_contigs_to_file( outfile, id_suffix );

  // open outfile for contigs that have been fused and removed from the contigs vector
  ofstream fusedout_fp( outfile+"_fused.fasta");

  // print out each line to the 
  for( int i=0; i<contigs_fused.size(); i++ ){
    fusedout_fp << ">" << contigs_fused[i].get_contig_id() << "_" << id_suffix << endl;
    fusedout_fp << get_contig_fused(i) << endl;
  }

  fusedout_fp.close();
  close_log();
} 

// return a string of the current time in the program
string Process::get_time(){
  string curr_time = "";
  char buffer[100];

  int t_now = difftime( time(0), timer );   // get time now
  
  int t_hour = t_now / 3600;
  t_now = t_now % 3600;

  int t_min = t_now / 60;
  t_now = t_now % 60;

  // format time in xxx:xx:xx
  sprintf( buffer, "%d:%.2d:%.2d", t_hour, t_min, t_now);
  curr_time = buffer;

  return curr_time;
}


// initalize logfile
void Process::logfile_init(){
  logfile = outfile + ".log";
  log_fs.open( logfile, fstream::out | fstream::trunc );

  // output starting option values
  log_fs << "OPTION VALUES" << endl;
  log_fs << "contig_sub_len: " << contig_sub_len << endl;
  log_fs << "extend_len: " << extend_len << endl;
  log_fs << "max_search_loops: " << max_search_loops << endl;
  log_fs << "max_sort_char: " << max_sort_char << endl;
  log_fs << "min_cov_init: " << min_cov_init << endl;
  log_fs << "min_overlap: " << min_overlap << endl;
  log_fs << "initial_trim: " << initial_trim << endl;
  log_fs << "max_missed: " << max_missed << endl;
  log_fs << "mismatch_threshold: " << mismatch_threshold << endl;
  log_fs << "max_threads: " << max_threads << endl << endl;
}

// prints notes to file as the program progresses
void Process::print_to_logfile( string note ){
  log_fs << get_time() << "\t" << note << endl;
}

// creates id of fused contigs
string Process::get_fused_id( string contig1_id, string contig2_id ){
  if( contig1_id.length() >= 5 && contig1_id.compare( 0, 5, "fused" ) == 0 ){
    contig1_id = contig1_id.substr( 5 );
  }
  
  if( contig1_id.length() >= 5 && contig2_id.compare( 0, 5, "fused" ) == 0 ){
    contig2_id = contig2_id.substr( 5 );
  }
    
  return contig1_id+"_<>_"+contig2_id;
}

// complete contig_fusion process
void Process::contig_fusion_log( Mismatch fusion ){ 
  // print messages to logfile about current actions
  print_to_logfile( "Contig fused: " );
  print_to_logfile( "  Overlap length: " + to_string(fusion.get_length()) );
  print_to_logfile( "  Mismatch_score: " + to_string(fusion.get_score()) );
  print_to_logfile( "  Contig_i: " + contigs[fusion.get_index_i()].get_contig_id() );
  print_to_logfile( "  Contig_j: " + contigs[fusion.get_index_j()].get_contig_id() );
  print_to_logfile( "" );
}

// complete contig_fusion process
void Process::commit_fusion( string fused, string fused_id, int index_i, int index_j ){ 
  contigs_fused.push_back( contigs[index_i] );
  print_to_logfile( "Contig moved to fused file: " + contigs[index_i].get_contig_id() );

  contigs_fused.push_back( contigs[index_j] );
  print_to_logfile( "Contig moved to fused file: " + contigs[index_j].get_contig_id() );

  print_to_logfile( "Committing: " + fused_id );
  print_to_logfile( "\t" + fused );
  contigs.push_back( Contig( fused, fused_id, 1, min_cov_init ));
}

// check overlap section for mismatches
Mismatch Process::overlap_check( string contig_a, string contig_b, int overlap, int end_i, int end_j ){
  Mismatch mim;
  mim.set_end_i( end_i );
  mim.set_end_j( end_j );
  double score_f = 1.0;
  double score_r = 1.0;
  double score = 1.0;

  for( int i=overlap-1; i>=min_overlap; i-- ){
    score_f = mismatch_score( contig_a.substr( contig_a.length() - i, i/2 ), contig_b.substr( 0, i/2 ) );
    score_r = mismatch_score( contig_a.substr( contig_a.length() - (i-i/2) ), contig_b.substr( i/2, i-i/2 ) );
    score = (score_f + score_r)/2;

    // check if this overlap is the best so far
    if( score < mim.get_score() ){
      mim.set_score( score );
      mim.set_length( i );
    }
    // if the score for the end of contig_a is better than the threshold, check the end of contig_b to see if there could be support in the reads for that end going in a different direction
    else if( score_r <= mismatch_threshold ){
      Contig temp_cont( contig_b.substr(i/2, i-i/2), "temp" );
      score_f = temp_cont.check_fusion_support( contig_a.substr( contig_a.length() - i, i/2 ) );
      // make sure the support score is at least as good as the threshold
      if( score_f <= mismatch_threshold ){
        score = (score_f + score_r)/2;
        if( score < mim.get_score() ){
          mim.set_score( score );
          mim.set_length( i );
        }
      }
    }
    // if the score for the end of contig_b is better than the threshold, check the end of contig_a to see if there could be support in the reads for that end going in a different direction
    else if( score_f <= mismatch_threshold ){
      Contig temp_cont( contig_a.substr( contig_a.length() - (i-i/2) ), "temp" );
      score_r = temp_cont.check_fusion_support( contig_b.substr( i/2, i-i/2 ) );
      // make sure the support score is at least as good as the threshold
      if( score_r <= mismatch_threshold ){
        score = (score_r + score_f)/2;
        if( score < mim.get_score() ){
          mim.set_score( score );
          mim.set_length( i );
        }
      }
    }
  }

  return mim;
}

// sort the match_list for easier 
void Process::sort_matches( vector<Mismatch> &match_list ){
  vector<Mismatch> sorted_list;
  Mismatch temp_mim;

  // sort by overlap length
  for( int i=0; i<match_list.size(); i++ ){
    for( int j=i; j<match_list.size(); j++ ){
      if( match_list[i].get_length() < match_list[j].get_length() ){
        temp_mim = match_list[i];
        match_list[i] = match_list[j];
        match_list[j] = temp_mim;
      }
    }
  }

  // sort by score
  for( int i=0; i<match_list.size(); i++ ){
    for( int j=i; j<match_list.size(); j++ ){
      if( match_list[i].get_score() > match_list[j].get_score() ){
        temp_mim = match_list[i];
        match_list[i] = match_list[j];
        match_list[j] = temp_mim;
      }
    }
  }
}

// cleans match_list from conflicting matches
void Process::clean_matches( vector<Mismatch> &match_list ){
  for( int i=0; i<(int)match_list.size()-1; i++ ){
    int index_i = match_list[i].get_index_i();
    int end_i = match_list[i].get_end_i();
    int index_j = match_list[i].get_index_j();
    int end_j = match_list[i].get_end_j();
    for( int j=i+1; j<match_list.size(); j++ ){
      if( match_list[j].get_index_i() == index_i && match_list[j].get_end_i() == end_i ){
        match_list.erase( match_list.begin() + j );
        j--;
      }
      else if( match_list[j].get_index_j() == index_i && match_list[j].get_end_j() == end_i ){
        match_list.erase( match_list.begin() + j );
        j--;
      }
      else if( match_list[j].get_index_j() == index_j && match_list[j].get_end_j() == end_j ){
        match_list.erase( match_list.begin() + j );
        j--;
      }
      else if( match_list[j].get_index_i() == index_j && match_list[j].get_end_i() == end_j ){
        match_list.erase( match_list.begin() + j );
        j--;
      }
    }
  }
}

// create fused contig string
string Process::build_fusion_string( string contig_a, string contig_b, int overlap ){
  string fused( contig_a.substr( 0, contig_a.length() - (overlap/2) ) );
  fused.append( contig_b.substr( overlap/2 + overlap%2 ) );

  return fused;
}

// remove duplicates from contig remove list
void Process::dedup_list( vector<int> &list ){
  for( int i=0; i<(int)list.size()-1; i++ ){
    for( int j=i+1; j<list.size(); j++ ){
      if( list[i] == list[j] ){
        list.erase( list.begin() + j );
        j--;
      }
    }
  }
}

// sort index list for removing contigs
void Process::sort_removals( vector<int> &remove_list ){
  int plc_hld;
  for( int i=0; i<(int)remove_list.size()-1; i++ ){
    int low_idx = i;
    for( int j=i+1; j<(int)remove_list.size(); j++ ){
      if( remove_list[j] < remove_list[low_idx] ){
        low_idx = j;
      }
    }
    if( low_idx != i ){
      plc_hld = remove_list[i];
      remove_list[i] = remove_list[low_idx];
      remove_list[low_idx] = plc_hld;
    }
  }
}

// remove fused contigs from contigs list
void Process::process_removals( vector<int> remove_list ){
  dedup_list( remove_list );
  sort_removals( remove_list );
  for( int i=(int)remove_list.size()-1; i>=0; i-- ){
    contigs.erase( contigs.begin() + remove_list[i] );
  }
}


// fuse contigs wherever possible
void Process::contig_fusion(){
  // mismatch scores variables 
  // orientation is an integer 0-3 and corresponds as follows:
  //    0:  i to j
  //    1:  i to j_rev
  //    2:  j to i
  //    3:  j_rev to i
  vector< Mismatch > match_list;
  int id = 0;
  
  /////////////////////////
  // GET MISMATCH SCORES //
  /////////////////////////

  // loop through each contig to get the end of the contig
  for( int i=0; i<contigs.size(); i++ ){
    for( int j=i+1; j<contigs.size(); j++ ){
      string contig_i( get_contig( i ) );
      string contig_j( get_contig( j ) );
      string contig_j_rev( revcomp( contig_j ) );
      int overlap = extend_len * 2;
      Mismatch mim;

      //// Processing of rear end of the contig
      // orientation: i to j
      if( overlap > contig_i.length() ){
        overlap = contig_i.length();
      }
      if( overlap > contig_j.length() ){
        overlap = contig_j.length();
      }
      mim = overlap_check( contig_i, contig_j, overlap, 1, 0 );

      if( mim.get_score() <= mismatch_threshold ){
        mim.set_indices( i, j );
        match_list.push_back( mim );
      }

      // orientation: i to j_rev
      mim = overlap_check( contig_i, contig_j_rev, overlap, 1, 1 );

      if( mim.get_score() <= mismatch_threshold ){
        mim.set_indices( i, j );
        match_list.push_back( mim );
      }

      //// Processing of front end of the contig
      // orientation: j to i
      mim = overlap_check( contig_j, contig_i, overlap, 0, 1);

      if( mim.get_score() <= mismatch_threshold ){
        mim.set_indices( i, j );
        match_list.push_back( mim );
      }

      // orientation: j_rev to i
      mim = overlap_check( contig_j_rev, contig_i, overlap, 0, 0 );

      if( mim.get_score() <= mismatch_threshold ){
        mim.set_indices( i, j );
        match_list.push_back( mim );
      }
    }
  }

  ////////////////
  // END SCORES //
  ////////////////
  
  ///////////////////////////////
  // SORT AND CLEAN MATCH LIST //
  ///////////////////////////////

  sort_matches( match_list );
  clean_matches( match_list );

  ////////////////////////
  // END SORT AND CLEAN //
  ////////////////////////

  for( int i=0; i<match_list.size(); i++ ){
    print_to_logfile( "score: " + to_string( match_list[i].get_score() ) );
    print_to_logfile( "i: " + to_string( match_list[i].get_index_i() ) +  " j: " + to_string( match_list[i].get_index_j() ) );
  }
 
  //////////////////
  // FUSE CONTIGS //
  //////////////////
  vector<int> contig_remove_list;
  for( int i=0; i<match_list.size(); i++ ){
    // fuse contigs
    int index_i = match_list[i].get_index_i();
    int index_j = match_list[i].get_index_j();
    string contig_i = contigs[index_i].get_contig();
    string contig_j = contigs[index_j].get_contig();
    string fused_id("");
    string fused("");
    string rev_i("");
    string rev_j("");

    // find the reverse compliment of the contigs where necessary
    if( match_list[i].get_end_i() == 0 ){
      contig_i = revcomp( contig_i );
      rev_i = "_r";
    }
    
    if( match_list[i].get_end_j() == 1 ){
      contig_j = revcomp( contig_j );
      rev_j = "_r";
    }

    contig_fusion_log( match_list[i] );
    fused_id = get_fused_id( contigs[index_i].get_contig_id()+rev_i, contigs[index_j].get_contig_id()+rev_j );
    fused_id = "fused(" + fused_id + ")";
 
    // push each index onto the remove vector
    contig_remove_list.push_back( index_i );
    contig_remove_list.push_back( index_j );

    // create fused contig
    fused = build_fusion_string( contig_i, contig_j, match_list[i].get_length() );
            
    // commit fusion here
    commit_fusion( fused, fused_id, index_i, index_j );
    int new_index = (int)contigs.size()-1;

    // check for additional appearances of the current contigs in match_list
    for( int j=i+1; j<match_list.size(); j++ ){
      // change reference index and end variables where necessary
      if( index_i == match_list[j].get_index_i() ){
        match_list[j].set_index_i( new_index );
        match_list[j].set_end_i( 0 );
      }
      else if( index_i == match_list[j].get_index_j() ){
        match_list[j].set_index_j( new_index );
        match_list[j].set_end_j( 0 );
      }
      else if( index_j == match_list[j].get_index_i() ){
        match_list[j].set_index_i( new_index );
        match_list[j].set_end_i( 1 );
      }
      else if( index_j == match_list[j].get_index_j() ){
        match_list[j].set_index_j( new_index );
        match_list[j].set_end_j( 1 );
      }
    }
  }

  // remove fused contigs from the vector
  process_removals( contig_remove_list );

  //////////////
  // END FUSE //
  //////////////
}

// Initializes data structures and turns over control to run_manager()
void Process::start_run(){
  // initialize timer
  time( &timer );

  logfile_init();

  // initialize reads and contigs
  add_reads();
  add_contigs();

  cout << "sort_reads start: ";
  print_time();
  sort_reads();
  cout << "sort_rc start: ";
  print_time();
  sort_rc();
  cout << "create_reads_range start: ";
  print_time();
  create_read_range();
  cout << "End initialization phase: ";
  print_time();
  
  // make initial attempt to fuse contigs  
  // removed for master branch until algorithm can be adjusted
  contig_fusion();
  
  if( test_run ){
    print_contigs_to_file( outfile + ".fus", "mid" );
  }

  run_manager();
}

// Manages run 
void Process::run_manager(){
  // create thread array with max_thread entries
  vector<thread> t;
  Queue<int> qu;

  // loop max search loops
  for( int j=0; j<max_search_loops; j++ ){
      
    // initialize threads
    for( int i=0; i<max_threads; i++ ){
      cout << "Thread" << i << endl;
      t.push_back(thread( thread_worker, ref(contigs), ref(qu), i ));
    }

    cout << "contigs.size(): " << contigs.size() << " max_threads: " << max_threads << endl;

    // push each thread onto queue
    for( int i=0; i<contigs.size(); i++ ){
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
    contig_fusion();
    
    if( test_run ){
      print_contigs_to_file( outfile + ".fus" + to_string(j), "mid" );
    }
  }
}

// closes logfile
void Process::close_log(){
  log_fs << get_time() << "\tProcess terminated" << endl;

  log_fs.close();
}

//////////////////////////////
// END DEFINITIONS ///////////
//////////////////////////////
