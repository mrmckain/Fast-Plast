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
  return contigs_fused[ contig_ind ].getContig();
}

// return contig with index contig_ind
string Process::get_contig( int contig_ind ){
  return contigs[ contig_ind ].getContig();
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

// complete contig_fusion process
void Process::contig_fusion_wrapup( string fused, string fused_id, int index_i, int index_j, int bp_added_fr, int bp_added_rr, int total_missed, string overlap_seq ){ 
  // print messages to logfile about current actions
  print_to_logfile( "Fused contigs to form " + fused_id ); 
  print_to_logfile( "Total of " + to_string(total_missed) + " basepairs did not match in the overlap section: " + overlap_seq );
  print_to_logfile( "Contig moved to fused file: " + contigs[index_i].get_contig_id() );
  print_to_logfile( "Contig moved to fused file: " + contigs[index_j].get_contig_id() );
  print_to_logfile( "" );

  // put pre-fused contigs in contigs_fused vector and post-fused contig in contigs vector
  contigs_fused.push_back( contigs[index_i] );
  contigs_fused.push_back( contigs[index_j] );
  contigs.push_back( Contig( fused, fused_id, 1, min_cov_init, bp_added_fr, bp_added_rr ));
  
  // remove pre-fused vectors from contigs vector
  if( index_i>index_j ){
    contigs.erase( contigs.begin() + index_i );
    contigs.erase( contigs.begin() + index_j );
  }
  else{
    contigs.erase( contigs.begin() + index_j );
    contigs.erase( contigs.begin() + index_i );
  }
}

// creates id of fused contigs
string Process::get_fused_id( string contig1_id, string contig2_id ){
  cout << "1contig1_id: |" << contig1_id << "|" << endl;
  cout << "1contig2_id: |" << contig2_id << "|" << endl;
  if( contig1_id.length() >= 5 && contig1_id.compare( 0, 5, "fused" ) == 0 ){
    contig1_id = contig1_id.substr( 5 );
  }
  cout << "2contig1_id: |" << contig1_id << "|" << endl;
  cout << "2contig2_id: |" << contig2_id << "|" << endl;
  
  if( contig1_id.length() >= 5 && contig2_id.compare( 0, 5, "fused" ) == 0 ){
    contig2_id = contig2_id.substr( 5 );
  }
  cout << "3contig1_id: |" << contig1_id << "|" << endl;
  cout << "3contig2_id: |" << contig2_id << "|" << endl;
    
  return "fused("+contig1_id+"_<>_"+contig2_id+")";
}

// complete contig_fusion process
void Process::contig_fusion_wrapup( string fused, string fused_id, int index_i, int index_j, int bp_added_fr, int bp_added_rr ){ 
  // print messages to logfile about current actions
  print_to_logfile( "Fused contigs to form " + fused_id ); 
  print_to_logfile( "Contig moved to fused file: " + contigs[index_i].get_contig_id() );
  print_to_logfile( "Contig moved to fused file: " + contigs[index_j].get_contig_id() );
  print_to_logfile( "" );

  // put pre-fused contigs in contigs_fused vector and post-fused contig in contigs vector
  contigs_fused.push_back( contigs[index_i] );
  contigs_fused.push_back( contigs[index_j] );
  contigs.push_back( Contig( fused, fused_id, 1, min_cov_init, bp_added_fr, bp_added_rr ));
  
  // remove pre-fused vectors from contigs vector
  if( index_i>index_j ){
    contigs.erase( contigs.begin() + index_i );
    contigs.erase( contigs.begin() + index_j );
  }
  else{
    contigs.erase( contigs.begin() + index_j );
    contigs.erase( contigs.begin() + index_i );
  }
}

// process the ends of the contigs for fusion at the front end of the second contig
bool Process::contig_end_compare_fr( int index_i, int index_j, int pos, int bp_added_fr_i, string contig_i, string i_rev ){
  string contig_j = get_contig( index_j );
  // length of overlapping section
  int len = tip_depth + pos;
  //substring of contig_i to do comparisons with that represents the overlapping section from contig_i's point of view
  string contig_i_sub = contig_i.substr( contig_i.length() - len );

  int missed_i = 0;
  int missed_j = 0;
  int total_missed = 0;
  bool fuse_success = false;

  // tally misses in end of contig_j
  for( int k=0; k<pos; k++ ){
    if( contig_j.compare(k, 1, contig_i_sub.substr(k, 1) ) ){     // contig_i_sub starts at the beginning of contig_j here
      missed_j++;
    }
  }

  // tally misses in end of contig_i
  for( int k=0; k<trim_length; k++ ){
    if( contig_j.compare(k+pos+tip_length, 1, contig_i_sub.substr(k+pos+tip_length, 1) ) ){
      missed_i++;
    }
  }

  
  total_missed = missed_i + missed_j;

  if( (missed_i <= max_missed) && (missed_j <= max_missed) ){
    fuse_success = true;
  }
  else if( missed_i <= max_missed ){
    Contig contig_alt( contig_j.substr( pos ), "temp" );
    if( contig_alt.check_fusion_support( contig_i_sub, pos-1, false ) ){
      fuse_success = true;
    }
  }
  else if( missed_j <= max_missed ){
    Contig contig_alt( contig_i.substr( 0, contig_i.length()-trim_length ), "temp" );
    if( contig_alt.check_fusion_support( contig_j.substr(pos+tip_length), contig_i.length()-trim_length+1, true ) ){
      fuse_success = true;
    }
  }


  if( fuse_success ){
    // form fused contig and its id
    string fused( contig_i.substr( 0, contig_i.length() - tip_depth ) );
    fused.append( contig_j.substr( pos, contig_j.length()-pos ) );
    string fused_id = get_fused_id( contigs[index_i].get_contig_id()+i_rev, contigs[index_j].get_contig_id());

    if( total_missed == 0 ){
      contig_fusion_wrapup( fused, fused_id, index_i, index_j, bp_added_fr_i, contigs[index_j].get_bp_added_rr() );
    }
    else{
      contig_fusion_wrapup( fused, fused_id, index_i, index_j, bp_added_fr_i, contigs[index_j].get_bp_added_rr(), total_missed, 
          fused.substr( contig_i.length() - pos - tip_depth, pos + tip_depth ) );
    }
    return true;
  }
  return false;
}

// process the ends of the contigs for fusion at the rear end of the second contig
bool Process::contig_end_compare_rr( int index_i, int index_j, int pos, int bp_added_rr_i, string contig_i, string i_rev ){
  string contig_j = get_contig( index_j );
  int len = contig_j.length() - pos + trim_length;
  string contig_i_sub = contig_i.substr( 0, len );

  int missed_i = 0;
  int missed_j = 0;
  int total_missed = 0;
  bool fuse_success = false;

  // tally misses for contig_i end
  for( int k=0; k<trim_length; k++ ){
    if( contig_j.compare(pos-trim_length+k, 1, contig_i_sub.substr(k, 1) ) ){
      missed_i++;
    }
  }

  // tally misses for contig_j end
  for( int k=0; k<contig_j.length()-pos-tip_length; k++ ){
    if( contig_j.compare(k+pos+tip_length, 1, contig_i_sub.substr(k+tip_depth, 1) ) ){
      missed_j++;
    }
  }


  total_missed = missed_i + missed_j;

  if( (missed_i <= max_missed) && (missed_j <= max_missed) ){
    fuse_success = true;
  }
  else if ( missed_j <= max_missed ){
    Contig contig_alt( contig_i.substr( trim_length ), "temp" );
    if( contig_alt.check_fusion_support( contig_j, pos-1, false ) ){
      fuse_success = true;
    }
  }
  else if ( missed_i <= max_missed ){
    Contig contig_alt( contig_j.substr( 0, pos+tip_length ), "temp" );
    if( contig_alt.check_fusion_support( contig_i, tip_depth, true ) ){
      fuse_success = true;
    }
  }


  if( fuse_success ){
    string fused( contig_j );
    fused.append( contig_i.substr( trim_length, contig_i.length()-trim_length ) );
    string fused_id = get_fused_id( contigs[index_j].get_contig_id(), contigs[index_i].get_contig_id()+i_rev );

    if( total_missed == 0 ){
      contig_fusion_wrapup( fused, fused_id, index_i, index_j, contigs[index_j].get_bp_added_fr(), bp_added_rr_i );
    }
    else{
      contig_fusion_wrapup( fused, fused_id, index_i, index_j, contigs[index_j].get_bp_added_fr(), bp_added_rr_i, total_missed,
          fused.substr( pos - trim_length, contig_j.length() - pos ) );
    }
    return true;
  }
  return false;
}

// Compares the ends of the contigs with indices index_i and index_j which are related to the contigs from Process::contig_fusion()
// back indicates whether contig_i comes off the front or back of contig_j and changes the behavior of pos as follows:
//  1:  pos indicates the position in contig_j from which contig_i starts
//  0:  pos indicates the position in contig_j to which contig_i extends
// rev indicates whether contig_i is the reverse compliment of the original contig
// returns boolean value that indicates if a fusion was made
bool Process::contig_end_compare( int index_i, int index_j, int pos, bool back, bool rev ){
  string log_text = "";
  string i_rev( "" );
  string contig_i = get_contig( index_i );
  
  // set variables for creating new contig objects if necessary
  int bp_added_rr_i = contigs[index_i].get_bp_added_rr();
  int bp_added_fr_i = contigs[index_i].get_bp_added_fr();

  // set variables impacted by direction of contig_i
  if( rev ){
    bp_added_rr_i = contigs[index_i].get_bp_added_fr();
    bp_added_fr_i = contigs[index_i].get_bp_added_rr();

    contig_i = revcomp( contig_i );
    i_rev = "_rev_comp";
  }

  // back section
  if( back && contig_end_compare_rr( index_i, index_j, pos, bp_added_rr_i, contig_i, i_rev ) ){
    return true;
  }
  // front section
  else if( !back && contig_end_compare_fr( index_i, index_j, pos, bp_added_fr_i, contig_i, i_rev ) ){
    return true;
  }

  // unsuccessful fusion, return false
  return false;
}

// front end search for contig_fusion algorithm
bool Process::contig_fusion_front_search( string contig_i, string contig_j, string contig_tip, string contig_rev_tip, int index_i, int index_j ){
  // end_depth places the cursor at the depth of bp_added. Because this is at the front of contig_j, the matching segment extends further into contig_j and overlaps the old bp's that 
  //    were not covered in the last run
  end_depth = contigs[index_j].get_bp_added_fr();
  
  if( end_depth == -1 ){
    end_depth = bp_added_init;
  } 

  // limit to prevent out_of_bounds errors
  if( end_depth > contig_i.length() ){
    end_depth = contig_i.length();
  }

  if( end_depth > contig_j.length() ){
    end_depth = contig_j.length() - tip_length;
  }

  //// search the front of the contig_j for the contig_tips created above
  for( int k=end_depth; k>=0; k-- ){
    if( contig_j.compare( k, tip_length, contig_tip ) == 0 ){
      if( contig_end_compare( index_i, index_j, k, false, false ) ){
        return true;
      }
    }
    
    if( contig_j.compare( k, tip_length, contig_rev_tip ) == 0 ){
      if( contig_end_compare( index_i, index_j, k, false, true ) ){
        return true;
      }
    }
  }
  return false;
}

// back end search for contig_fusion algorithm
bool Process::contig_fusion_rear_search( string contig_i, string contig_j, string contig_tip, string contig_rev_tip, int index_i, int index_j ){
  // end_depth includes tip_length to cover the old bp's that are less than tip_length from the end from last run
  end_depth = contigs[index_j].get_bp_added_rr() + tip_length;

  // limit to prevent out_of_bounds errors
  if( end_depth > contig_i.length() ){
    end_depth = contig_i.length();
  }

  //// search the back of the contig_j for the contig_tips created above
  for( int k=contig_j.length()-end_depth; k<contig_j.length()-tip_length; k++ ){
    if( contig_j.compare( k, tip_length, contig_tip ) == 0 ){
      if( contig_end_compare( index_i, index_j, k, true, false ) ){
        return true;
      }
    }
    
    if( contig_j.compare( k, tip_length, contig_rev_tip ) == 0 ){
      if( contig_end_compare( index_i, index_j, k, true, true ) ){
        return true;
      }
    }
  }
  return false;
}

// fuse contigs wherever possible
void Process::contig_fusion(){
  // loop through each contig to get the end of the contig
  for( int i=0; i<contigs.size(); i++ ){
    string contig_i( get_contig( i ) );
    string contig_i_rev( revcomp( contig_i ) );

    // protect against high values of tip_length and trim_length
    if( tip_depth > contig_i.length() ){
      continue;
    }

    // front end tip of the contig and that of the reverse compliment of the contig
    string contig_tip_fr( contig_i.substr( trim_length, tip_length ) );
    string contig_rev_tip_fr( contig_i_rev.substr( trim_length, tip_length ) );

    // rear end tip of the contig and that of the reverse compliment of the contig
    string contig_tip_rr( contig_i.substr( get_contig(i).length()-tip_depth, tip_length ) );
    string contig_rev_tip_rr( contig_i_rev.substr( get_contig(i).length()-tip_depth, tip_length ) );

    // loop through contigs and use this contig as the base to compare the end of contig_i and contig_i_rev to
    for( int j=0; j<contigs.size(); j++ ){
      if( j != i ){
        string contig_j( get_contig( j ) );

        //cout << "1 index_i: " << i << " index_j: " << j << " contig_i.length(): " << contig_i.length() << " contig_j.length(): " << contig_j.length() << endl;
        //////////////////////
        // FRONT END SEARCH //
        //////////////////////
        if( contig_fusion_front_search( contig_i, contig_j, contig_tip_rr, contig_rev_tip_rr, i, j ) ){
          // if the two contigs match and are fused, adjustments must be made to the index markers so all contigs are checked
          if( i>j ){
            i-=2;
          }
          else{
            i--;
          }
          break;
        }
        //cout << "2 index_i: " << i << " index_j: " << j << endl;

        /////////////////////
        // BACK END SEARCH //
        /////////////////////
        if( contig_fusion_rear_search( contig_i, contig_j, contig_tip_fr, contig_rev_tip_fr, i, j ) ){
          // if the two contigs match and are fused, adjustments must be made to the index markers so all contigs are checked
          if( i>j ){
            i-=2;
          }
          else{
            i--;
          }
          break;
        }
        //cout << "3 index_i: " << i << " index_j: " << j << endl;
      }
    }
  }

  // reset bp_added variables to 0 for each contig
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
  cout << "End initialization phase";
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
