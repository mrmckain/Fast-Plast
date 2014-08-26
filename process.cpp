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
#include "queue.tcc"
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
int trim_length;
int tip_length;
int end_depth;
int tip_depth;
int initial_trim;
int max_missed;
int bp_added_init;
int mismatch_threshold;

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
              readlist.push_back( line );
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
            readlist.push_back( buffer);
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
    readlist.push_back( buffer );
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
      // set the value of bp_added to -1 to indicate ignoring when extending contigs 
      contigs[i].set_bp_added( -1 );
      //contigs.push_back( contigs[i] );
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
          contigs.push_back( Contig( buffer, contig_id, parse_cov( contig_id ), min_cov_init, bp_added_init ));
        }
        else if( buffer.length() > contig_sub_len ){
          int trim = (buffer.length() - contig_sub_len) / 2;
          buffer = buffer.substr( trim, buffer.length() - 2*trim );
          contigs.push_back( Contig( buffer, contig_id, parse_cov( contig_id ), min_cov_init, bp_added_init ));
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
    contigs.push_back( Contig( buffer, contig_id, parse_cov( contig_id ), min_cov_init, bp_added_init ) );
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

// prints results to fasta file with outfile prefix and additional information is printed to a text based file with outfile prefix
void Process::print_to_outfile(){
  // remove directories from outfile to form id_suffix if necessary
  size_t id_suffix_pos = outfile.find_last_of( "/" );
  string id_suffix = outfile;
  if( id_suffix_pos != string::npos ){
    id_suffix = id_suffix.substr( id_suffix_pos + 1 );
  }

  // open outfile
  ofstream outfile_fp( outfile+".fasta");

  // print out each line to the 
  for( int i=0; i<contigs.size(); i++ ){
    outfile_fp << ">" << contigs[i].get_contig_id() << "_" << id_suffix << endl;
    outfile_fp << get_contig(i) << endl;
  }

  outfile_fp.close();

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

  int t_now = difftime( time(0), timer );   // get time now
  
  int t_hour = t_now / 3600;
  t_now = t_now % 3600;

  int t_min = t_now / 60;
  t_now = t_now % 60;

  curr_time = to_string( t_hour ) + ":" + to_string( t_min ) + ":" + to_string( t_now ); 
  
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
  log_fs << "trim_length: " << trim_length << endl;
  log_fs << "tip_length: " << tip_length << endl;
  log_fs << "end_depth: " << end_depth << endl;
  log_fs << "tip_depth: " << tip_depth << endl;
  log_fs << "initial_trim: " << initial_trim << endl;
  log_fs << "max_missed: " << max_missed << endl;
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
void Process::commit_fusion( string fused, string fused_id, vector<Mismatch> fusion_chain ){ 
  contigs_fused.push_back( contigs[fusion_chain[0].get_index_i()] );
  print_to_logfile( "Contig moved to fused file: " + contigs[fusion_chain[0].get_index_i()].get_contig_id() );

  for( int i=0; i<fusion_chain.size(); i++ ){
    contigs_fused.push_back( contigs[fusion_chain[i].get_index_j()] );
    print_to_logfile( "Contig moved to fused file: " + contigs[fusion_chain[i].get_index_j()].get_contig_id() );
  }
  
  contigs.push_back( Contig( fused, fused_id, 1, min_cov_init, contigs[fusion_chain[0].get_index_i()].get_bp_added_fr(), contigs[fusion_chain.back().get_index_j()].get_bp_added_rr() ));
}

// tally mismatches in substrings passed and return score in the form of misatches per length
double Process::mismatch_score( string contig_sub_a, string contig_sub_b ){
  int mismatch = 0;

  for( int i=0; i<contig_sub_a.length(); i++ ){
    if( contig_sub_a[i] != contig_sub_b[i] ){
      mismatch++;
    }
  }

  return double(mismatch) / contig_sub_a.length();
}

// check overlap section for mismatches
Mismatch Process::overlap_check( string contig_a, string contig_b, int overlap, int orientation ){
  Mismatch mim;
  mim.set_orientation( orientation );
  double score = 1.0;

  for( int i=min_overlap; i<overlap; i++ ){
    score = mismatch_score( contig_a.substr( contig_a.length() - i, i ), contig_b.substr( 0, i ) );

    // check if this overlap is the best so far
    if( score < mim.get_score() ){
      mim.set_score( score );
      mim.set_length( overlap );
    }
  }

  return mim;
}

// returns overlap length based on the indexes passed
int Process::get_overlap( int i, int j, int orientation ){
  int overlap = min_overlap;
  int bp_added_a, bp_added_b;
  int length_a, length_b;

  switch( orientation ){
    case 0:
      bp_added_a = contigs[i].get_bp_added_rr();
      bp_added_b = contigs[j].get_bp_added_fr();
      length_a = get_contig( i ).length();
      length_b = get_contig( j ).length();
      break;
    case 1:
      bp_added_a = contigs[i].get_bp_added_rr();
      bp_added_b = contigs[j].get_bp_added_rr();
      length_a = get_contig( i ).length();
      length_b = get_contig( j ).length();
      break;
    case 2:
      bp_added_a = contigs[j].get_bp_added_rr();
      bp_added_b = contigs[i].get_bp_added_fr();
      length_a = get_contig( j ).length();
      length_b = get_contig( i ).length();
      break;
    case 3:
      bp_added_a = contigs[j].get_bp_added_fr();
      bp_added_b = contigs[i].get_bp_added_fr();
      length_a = get_contig( j ).length();
      length_b = get_contig( i ).length();
      break;
    default:
      return 0;
      break;
  }

  if( bp_added_a < bp_added_b ){
    if( bp_added_b < length_a ){
      overlap = bp_added_b;
    }
    else{
      overlap = length_a - 1;
    }
  }
  else{
    if( bp_added_a < length_b ){
      overlap = bp_added_a;
    }
    else{
      overlap = length_b - 1;
    }
  }

  if( overlap < min_overlap ){
    overlap = min_overlap;
  }

  return overlap;
}

// create fused contig string
string Process::build_fusion_string( string contig_a, string contig_b, int overlap ){
  string fused( contig_a.substr( 0, contig_a.length() - (overlap/2) ) );
  fused.append( contig_b.substr( overlap/2 + overlap%2 ) );

  return fused;
}

// finds a chained list of mismatch objects, sorts them and orients them properly
vector< Mismatch > Process::find_chain( vector< Mismatch > &fusion_list ){
  vector< Mismatch > fusion_chain;
  fusion_chain.push_back( fusion_list[0] );
  int front_index = fusion_list[0].get_index_i();
  int rear_index = fusion_list[0].get_index_j();

  for( int i=fusion_list.size()-1; i>0; i--){
    // check orientation to mark if the match should be reversed to put index_i first
    if( fusion_list[i].get_orientation() < 2 ){
      // j<-i => contig
      if( fusion_list[i].get_index_i() == front_index ){
        fusion_list[i].set_rev(true);
        front_index = fusion_list[i].get_index_j();
        fusion_chain.insert( fusion_chain.begin(), fusion_list[i] );
        fusion_list.erase( fusion_list.begin() + i );
        continue;
      }
      // i->j => contig
      else if( fusion_list[i].get_index_j() == front_index ){
        front_index = fusion_list[i].get_index_i();
        fusion_chain.insert( fusion_chain.begin(), fusion_list[i] );
        fusion_list.erase( fusion_list.begin() + i );
        continue;
      }
      // contig => i->j
      else if( fusion_list[i].get_index_i() == rear_index ){
        rear_index = fusion_list[i].get_index_j();
        fusion_chain.push_back( fusion_list[i] );
        fusion_list.erase( fusion_list.begin() + i );
        continue;
      }
      // contig => j<-i
      else if( fusion_list[i].get_index_j() == rear_index ){
        fusion_list[i].set_rev(true);
        rear_index = fusion_list[i].get_index_i();
        fusion_chain.push_back( fusion_list[i] );
        fusion_list.erase( fusion_list.begin() + i );
        continue;
      }
    }
    else{
      // i->j => contig
      if( fusion_list[i].get_index_j() == front_index ){
        fusion_list[i].set_rev(true);
        front_index = fusion_list[i].get_index_i();
        fusion_chain.insert( fusion_chain.begin(), fusion_list[i] );
        fusion_list.erase( fusion_list.begin() + i );
        continue;
      }
      // j<-i => contig
      else if( fusion_list[i].get_index_i() == front_index ){
        front_index = fusion_list[i].get_index_j();
        fusion_chain.insert( fusion_chain.begin(), fusion_list[i] );
        fusion_list.erase( fusion_list.begin() + i );
        continue;
      }
      // contig => j<-i
      else if( fusion_list[i].get_index_j() == rear_index ){
        rear_index = fusion_list[i].get_index_i();
        fusion_chain.push_back( fusion_list[i] );
        fusion_list.erase( fusion_list.begin() + i );
        continue;
      }
      // contig => i->j
      else if( fusion_list[i].get_index_i() == rear_index ){
        fusion_list[i].set_rev(true);
        rear_index = fusion_list[i].get_index_j();
        fusion_chain.push_back( fusion_list[i] );
        fusion_list.erase( fusion_list.begin() + i );
        continue;
      }
    }
  }

  fusion_list.erase( fusion_list.begin() );
  return fusion_chain;
}

// remove fused contigs from contigs list
void Process::process_removals( vector<int> remove_list ){
  for( int i=remove_list.size() - 1; i>=0; i-- ){
    contigs.erase( contigs.begin() + i );
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
  vector< vector< Mismatch >> mismatch_scores;
  int id = 0;

  /////////////////////////
  // GET MISMATCH SCORES //
  /////////////////////////

  // loop through each contig to get the end of the contig
  for( int i=0; i<contigs.size(); i++ ){
    // vector to hold best mismatch scores for the current contig
    vector< Mismatch > mim_vec;
    
    // push two Mismatch objects onto the vector to account for each end of the contig
    mim_vec.push_back( Mismatch() );
    mim_vec.push_back( Mismatch() );

    // cycle through the remaining contigs
    for( int j=i+1; j<contigs.size(); i++ ){
      string contig_i( get_contig( i ) );
      string contig_j( get_contig( j ) );
      string contig_j_rev( revcomp( contig_j ) );

      // cycle through each possible orientation of the contigs
      for( int k=0; k<4; k++ ){
        Mismatch curr_mim;
        int overlap = get_overlap( i, j, k );

        switch( k ){
          case 0:
            curr_mim = overlap_check( contig_i, contig_j, overlap, k );
            
            // replace back end position of mim_vec if the current mismatch score is better than the previous best
            if( curr_mim.get_score() < mim_vec[0].get_score() ){
              curr_mim.set_id( id );
              id++;
              curr_mim.set_indices( i, j );
              mim_vec[0] = curr_mim;
            }
            break;
          case 1:
            curr_mim = overlap_check( contig_i, contig_j_rev, overlap, k );
            
            // replace back end position of mim_vec if the current mismatch score is better than the previous best
            if( curr_mim.get_score() < mim_vec[0].get_score() ){
              curr_mim.set_id( id );
              id++;
              curr_mim.set_indices( i, j );
              mim_vec[0] = curr_mim;
            }
            break;
          case 2:
            curr_mim = overlap_check( contig_j, contig_i, overlap, k );
            
            // replace front end position of mim_vec if the current mismatch score is better than the previous best
            if( curr_mim.get_score() < mim_vec[1].get_score() ){
              curr_mim.set_id( id );
              id++;
              curr_mim.set_indices( i, j );
              mim_vec[1] = curr_mim;
            }
            break;
          case 3:
            curr_mim = overlap_check( contig_j_rev, contig_i, overlap, k );
            
            // replace front end position of mim_vec if the current mismatch score is better than the previous best
            if( curr_mim.get_score() < mim_vec[1].get_score() ){
              curr_mim.set_id( id );
              id++;
              curr_mim.set_indices( i, j );
              mim_vec[1] = curr_mim;
            }
            break;
          default:
            continue;
            break;
        }
      }
    }

    // loop through previously processed contigs
    for( int j=0; j<i; j++ ){
      // check back end best score of contig_j
      if( mismatch_scores[j][0].get_index_j() == i ){
        if( mismatch_scores[j][0].get_orientation() == 0 ){
          if( mismatch_scores[j][0].get_score() < mim_vec[1].get_score() ){
            mim_vec[1] = mismatch_scores[j][0];
            mim_vec[1].set_indices( mismatch_scores[j][0].get_index_j(), mismatch_scores[j][0].get_index_i() );
            mim_vec[1].set_orientation( 2 );
          }
          else if( mismatch_scores[j][0].get_score() >= mim_vec[1].get_score() ){
            mismatch_scores[j][0] = mim_vec[1];
            mismatch_scores[j][0].set_indices( mim_vec[1].get_index_i(), mim_vec[1].get_index_j() );
            mismatch_scores[j][0].set_orientation( 2 );
          }
        }
        else if( mismatch_scores[j][0].get_orientation() == 1 ){
          if( mismatch_scores[j][0].get_score() < mim_vec[0].get_score() ){
            mim_vec[0] = mismatch_scores[j][0];
            mim_vec[0].set_indices( mismatch_scores[j][0].get_index_j(), mismatch_scores[j][0].get_index_i() );
          }
          else if( mismatch_scores[j][0].get_score() >= mim_vec[0].get_score() ){
            mismatch_scores[j][0] = mim_vec[0];
            mismatch_scores[j][0].set_indices( mim_vec[0].get_index_i(), mim_vec[0].get_index_j() );
          }
        }
      }
      
      // check front end best score of contig_j
      if( mismatch_scores[j][1].get_index_j() == i ){
        if( mismatch_scores[j][1].get_orientation() == 2 ){
          if( mismatch_scores[j][1].get_score() < mim_vec[0].get_score() ){
            mim_vec[0] = mismatch_scores[j][1];
            mim_vec[0].set_indices( mismatch_scores[j][1].get_index_j(), mismatch_scores[j][1].get_index_i() );
            mim_vec[0].set_orientation( 0 );
          }
          else if( mismatch_scores[j][1].get_score() >= mim_vec[0].get_score() ){
            mismatch_scores[j][1] = mim_vec[0];
            mismatch_scores[j][1].set_indices( mim_vec[0].get_index_i(), mim_vec[0].get_index_j() );
            mismatch_scores[j][1].set_orientation( 0 );
          }
        }
        else if( mismatch_scores[j][1].get_orientation() == 3 ){
          if( mismatch_scores[j][1].get_score() < mim_vec[1].get_score() ){
            mim_vec[1] = mismatch_scores[j][1];
            mim_vec[1].set_indices( mismatch_scores[j][1].get_index_j(), mismatch_scores[j][1].get_index_i() );
          }
          else if( mismatch_scores[j][1].get_score() >= mim_vec[1].get_score() ){
            mismatch_scores[j][1] = mim_vec[1];
            mismatch_scores[j][1].set_indices( mim_vec[1].get_index_i(), mim_vec[1].get_index_j() );
          }
        }
      }
    }

    // add mismatch scores for current contig to mismatch_scores
    mismatch_scores.push_back( mim_vec );
  }

  ////////////////
  // END SCORES //
  ////////////////

  
  /////////////////////////
  // CREATE FUSION QUEUE //
  /////////////////////////
  vector<Mismatch> fusion_list;

  // loop through mismatch_scores
  for( int i=0; i<mismatch_scores.size(); i++ ){
    int match_index = mismatch_scores[i][0].get_index_j();
    if( match_index > i && mismatch_scores[i][0].get_score() <= mismatch_threshold ){
      if( mismatch_scores[i][0].get_id() == mismatch_scores[match_index][1].get_id() || mismatch_scores[i][0].get_id() == mismatch_scores[match_index][0].get_id() ){
        fusion_list.push_back( mismatch_scores[i][0] );
      }
    }
    
    match_index = mismatch_scores[i][1].get_index_j();
    if( match_index > i && mismatch_scores[i][1].get_score() <= mismatch_threshold ){      
      if( mismatch_scores[i][1].get_id() == mismatch_scores[match_index][0].get_id() || mismatch_scores[i][1].get_id() == mismatch_scores[match_index][1].get_id() ){
        fusion_list.push_back( mismatch_scores[i][1] );
      }
    }
  }

  ////////////////
  // END CREATE //
  ////////////////

  //////////////////
  // FUSE CONTIGS //
  //////////////////
  vector< int > contig_remove_list;
////////////////////////////////////////////////////////// TASK:: Clean this up so that unnecessary storage of contig strings is not occurring  
  while( !fusion_list.empty() ){
    vector<Mismatch> fusion_chain = find_chain(fusion_list);
    int bp_added_fr = 0;
    int bp_added_rr = 0;
    string fused_id("");
    // contig_i will represent the fused contig from the previous iteration
    string contig_i("");
    string contig_j("");
    string fused("");
    int index_i = 0;
    int index_j = 0;

    // loop through the fusion_chain created by find_chain() in order to create one contig to add
    for( int i=0; i<fusion_chain.size(); i++){
      Mismatch current_fuse = fusion_chain[i];
      index_i = current_fuse.get_index_i();
      index_j = current_fuse.get_index_j();
      contig_i = contigs[index_i].get_contig();
      contig_j = contigs[index_j].get_contig();
      int orientation = current_fuse.get_orientation();

      // push each index onto the remove vector
      contig_remove_list.push_back( index_i );

      // log contig_fusion
      contig_fusion_log( current_fuse );

      // reverse the indices and revcomp the contigs if the fusion represented is reversed
      if( current_fuse.get_rev() ){
        index_i = current_fuse.get_index_j();
        index_j = current_fuse.get_index_i();
        contig_i = revcomp( contig_j );

        // to avoid unnecessary double revcomps
        if( orientation != 1 ){
          contig_j = revcomp( contig_i );
        }
      }
      else if( orientation == 1 ){
        contig_j = revcomp( contig_j );
      }

      // set fused_id and if this is the first iteration, set fused
      if( fused_id == "" ){
        fused_id = get_fused_id( contigs[index_i].get_contig_id(), contigs[index_j].get_contig_id() );
      }
      else{
        fused_id = get_fused_id( fused_id, contigs[index_j].get_contig_id() );
      }

      // create fused contig
      fused = build_fusion_string( contig_i, contig_j, current_fuse.get_length() );
            
      // log fusion during each iteration for ease of recording details such as mismatch_score
    }
      
    // push the last index onto the remove vector
    contig_remove_list.push_back( index_j );
    
    fused_id = "fused(" + fused_id + ")";

    // commit_fusion() here..
    commit_fusion( fused, fused_id, fusion_chain );

    // remove fused contigs from the vector
    process_removals( contig_remove_list );
  }

  //////////////
  // END FUSE //
  //////////////

  // reset bp_added variables for each contig
  for( int i=0; i<contigs.size(); i++ ){
    contigs[i].reset_bp_added();
  }
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

  run_manager();
}

// Manages run 
void Process::run_manager(){
  // create thread array with max_thread entries
  thread t[max_threads];
  Queue<int> qu;

  // loop max search loops
  for( int j=0; j<max_search_loops; j++ ){
      
    // initialize threads
    for( int i=0; i<max_threads; i++ ){
      cout << "Thread" << i << endl;
      t[i] = thread( thread_worker, ref(contigs), ref(qu), i );
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

    // removed for master branch until algorithm can be adjusted
    contig_fusion();
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
