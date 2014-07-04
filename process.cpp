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
double Process::get_cov( string contig_id ){
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
    string contig_id = contigs[i].get_contig_id();
    double cov = get_cov( contig_id );
    cov_total += cov;
  }

  // find average of all cov values
  double cov_avg = cov_total / total_contigs;
  
  for( int i=0; i<total_contigs; i++ ){
    // get contig_id, parse out cov_## and compare this value to the avg
    string contig_id = contigs[i].get_contig_id();
    double cov = get_cov( contig_id );

    if( cov > cov_avg * 2.0 ){
      // Push an extra copy of the contig onto contigs and prepend "2x_" onto the contig_id
      contigs[i].set_contig_id( contig_id.insert( 0, "2x_" ) );
      //contigs.push_back( contigs[i] );
    }
  }
}

// put contigs from contfile into contlist
void Process::add_contigs(){
  int bp_added_init = 40;
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
        //cout << buffer << endl;
        contigs.push_back( Contig( buffer, contig_id, min_cov_init, bp_added_init ));
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
    contigs.push_back( Contig( buffer, contig_id, min_cov_init, bp_added_init ) );
  }

  contig_cov();
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
  log_fs << "max_threads: " << max_threads << endl << endl;
}

// prints notes to file as the program progresses
void Process::print_to_logfile( string note ){
  log_fs << get_time() << "\t" << note << endl;
}

// Compares the ends of the contigs with indices index_i and index_j which are related to the contigs from Process::contig_fusion()
// back indicates whether contig_i comes off the front or back of contig_j and changes the behavior of pos as follows:
//  1:  pos indicates the position in contig_j from which contig_i starts
//  0:  pos indicates the position in contig_j to which contig_i extends
// rev indicates whether contig_i is the reverse compliment of the original contig
// returns boolean value that indicates if a fusion was made
bool Process::contig_end_compare( int index_i, int index_j, int pos, bool back, bool rev ){
  string log_text = "";
  int max_missed = 10;
  string contig_j = get_contig( index_j );
  string contig_i = get_contig( index_i );
  string i_rev( "" );
  
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
  if( back ){
    pos -= trim_length;
    int len = contig_j.length() - pos;
    string contig_i_sub = contig_i.substr( 0, len );
    if( contig_j.compare( pos, len, contig_i_sub ) == 0 ){
      // form fused contig and its id
      string fused( contig_j );
      fused.append( contig_i.substr( len, contig_i.length()-len ) );
      string fused_id( "fused("+contigs[index_j].get_contig_id()+"_||_"+contigs[index_i].get_contig_id()+i_rev+")" );

      // print messages to logfile about current actions
      print_to_logfile( "Fused contigs to form " + fused_id ); 
      print_to_logfile( "Contig moved to fused file: " + contigs[index_i].get_contig_id() );
      print_to_logfile( "Contig moved to fused file: " + contigs[index_j].get_contig_id() );
      print_to_logfile( "" );

      // put pre-fused contigs in contigs_fused vector and post-fused contig in contigs vector
      contigs_fused.push_back( contigs[index_i] );
      contigs_fused.push_back( contigs[index_j] );
      contigs.push_back( Contig( fused, fused_id, min_cov_init, contigs[index_j].get_bp_added_fr(), bp_added_rr_i ));
      
      // remove pre-fused vectors from contigs vector
      if( index_i>index_j ){
        contigs.erase( contigs.begin() + index_i );
        contigs.erase( contigs.begin() + index_j );
      }
      else{
        contigs.erase( contigs.begin() + index_j );
        contigs.erase( contigs.begin() + index_i );
      }

      // successful fusion, return true
      return true;
    }
    // check to see if the unmatched portion is in the trim_length bp's at the end of either contig, if it is and the number of mismatched bp's is under max_missed, 
    //    then add contigs together ignoring trim_length section
    else{
      pos += trim_length;
      int overlap = len;
      len -= 2*trim_length;

      // check to make sure len>=0
      if( len>=0 ){
        if( contig_j.compare( pos, len, contig_i_sub.substr( trim_length, len ) ) == 0 ){
          string trim_section_i = contig_i.substr( 0, trim_length );
          string trim_section_j = contig_j.substr( contig_j.length()-trim_length, trim_length );
          int total_missed_i = 0;
          int total_missed_j = 0;

          // loop through bp's in trim_section for contig_i
          for( int k=0; k<trim_length; k++ ){
            // check trim_section_i
            if( trim_section_i.compare( k, 1, contig_j.substr( contig_j.length() - overlap + k, 1 ) ) != 0 ){
              total_missed_i++;
            }

            // check trim_section_j
            if( trim_section_j.compare( k, 1, contig_i.substr( overlap - trim_length + k, 1 ) ) != 0 ){
              total_missed_j++;
            }
          }

          // verify each trim section doesn't have too many misses
          if( total_missed_i <= max_missed && total_missed_j <= max_missed ){
            int total_missed = total_missed_i + total_missed_j;
            // form fused contig and its id
            string fused( contig_j.substr( 0, contig_j.length() - trim_length ) );
            fused.append( contig_i.substr( len + trim_length ) );
            string fused_id( "fused("+contigs[index_j].get_contig_id()+"_||_"+contigs[index_i].get_contig_id()+i_rev+")" );

            // print messages to logfile about current actions
            print_to_logfile( "Fused contigs to form " + fused_id ); 
            print_to_logfile( "Total of " + to_string(total_missed) + " basepairs did not match in the overlap section: " + contig_j.substr( contig_j.length() - overlap ) );
            print_to_logfile( "Contig moved to fused file: " + contigs[index_i].get_contig_id() );
            print_to_logfile( "Contig moved to fused file: " + contigs[index_j].get_contig_id() );
            print_to_logfile( "" );

            // put pre-fused contigs in contigs_fused vector and post-fused contig in contigs vector
            contigs_fused.push_back( contigs[index_i] );
            contigs_fused.push_back( contigs[index_j] );
            contigs.push_back( Contig( fused, fused_id, min_cov_init, contigs[index_j].get_bp_added_fr(), bp_added_rr_i ));
            
            // remove pre-fused vectors from contigs vector
            if( index_i>index_j ){
              contigs.erase( contigs.begin() + index_i );
              contigs.erase( contigs.begin() + index_j );
            }
            else{
              contigs.erase( contigs.begin() + index_j );
              contigs.erase( contigs.begin() + index_i );
            }

            // successful fusion, return true
            return true;
          }
        }
      }
    }
  }
  // front section
  else{
    int len = tip_depth + pos;
    string contig_i_sub = contig_i.substr( contig_i.length() - len );
    if( contig_j.compare( 0, len, contig_i_sub ) == 0 ){
      // form fused contig and its id
      string fused( contig_i );
      fused.append( contig_j.substr( len ) );
      string fused_id( "fused("+contigs[index_i].get_contig_id()+i_rev+"_||_"+contigs[index_j].get_contig_id()+")" );

      // print messages to logfile about current actions
      print_to_logfile( "Fused contigs to form " + fused_id ); 
      print_to_logfile( "Contig moved to fused file: " + contigs[index_i].get_contig_id() );
      print_to_logfile( "Contig moved to fused file: " + contigs[index_j].get_contig_id() );
      print_to_logfile( "" );

      // put pre-fused contigs in contigs_fused vector and post-fused contig in contigs vector
      contigs_fused.push_back( contigs[index_i] );
      contigs_fused.push_back( contigs[index_j] );
      contigs.push_back( Contig( fused, fused_id, min_cov_init, bp_added_fr_i, contigs[index_j].get_bp_added_rr() ));
      
      // remove pre-fused vectors from contigs vector
      if( index_i>index_j ){
        contigs.erase( contigs.begin() + index_i );
        contigs.erase( contigs.begin() + index_j );
      }
      else{
        contigs.erase( contigs.begin() + index_j );
        contigs.erase( contigs.begin() + index_i );
      }

      // successful fusion, return true
      return true;
    }
    // check to see if the unmatched portion is in the trim_length bp's at the end of either contig, if it is and the number of mismatched bp's is under max_missed, 
    //    then add contigs together ignoring trim_length section
    else{
      int overlap = len;
      len -= 2*trim_length;

      // check to make sure len>=0
      if( len>=0 ){
        if( contig_j.compare( trim_length, len, contig_i_sub.substr( trim_length, len ) ) == 0 ){
          string trim_section_i = contig_i.substr( contig_i.length()-trim_length, trim_length );
          string trim_section_j = contig_j.substr( 0, trim_length );
          int total_missed_i = 0;
          int total_missed_j = 0;

          // loop through bp's in trim_section for contig_i
          for( int k=0; k<trim_length; k++ ){
            // check trim_section_i
            if( trim_section_i.compare( k, 1, contig_j.substr( overlap - trim_length + k, 1 ) ) != 0 ){
              total_missed_i++;
            }

            // check trim_section_j
            if( trim_section_j.compare( k, 1, contig_i.substr( contig_i.length() - overlap + k, 1 ) ) != 0 ){
              total_missed_j++;
            }
          }

          // verify each trim section doesn't have too many misses
          if( total_missed_i <= max_missed && total_missed_j <= max_missed ){
            int total_missed = total_missed_i + total_missed_j;
            // form fused contig and its id
            string fused( contig_i.substr( 0, contig_i.length() - trim_length ) );
            fused.append( contig_j.substr( len + trim_length ) );
            string fused_id( "fused("+contigs[index_i].get_contig_id()+i_rev+"_||_"+contigs[index_j].get_contig_id()+")" );

            // print messages to logfile about current actions
            print_to_logfile( string( "Fused contigs to form " + fused_id ) ); 
            print_to_logfile( "Total of " + to_string(total_missed) + " basepairs did not match in the overlap section: " + contig_i.substr( contig_i.length() - overlap ) );
            print_to_logfile( "Contig moved to fused file: " + contigs[index_i].get_contig_id() );
            print_to_logfile( "Contig moved to fused file: " + contigs[index_j].get_contig_id() );
            print_to_logfile( "" );

            // put pre-fused contigs in contigs_fused vector and post-fused contig in contigs vector
            contigs_fused.push_back( contigs[index_i] );
            contigs_fused.push_back( contigs[index_j] );
            contigs.push_back( Contig( fused, fused_id, min_cov_init, bp_added_fr_i, contigs[index_j].get_bp_added_rr() ));
            
            // remove pre-fused vectors from contigs vector
            if( index_i>index_j ){
              contigs.erase( contigs.begin() + index_i );
              contigs.erase( contigs.begin() + index_j );
            }
            else{
              contigs.erase( contigs.begin() + index_j );
              contigs.erase( contigs.begin() + index_i );
            }

            // successful fusion, return true
            return true;
          }
        }
      }
    }
  }

  // unsuccessful fusion, return false
  return false;
}

// front end search for contig_fusion algorithm
bool Process::contig_fusion_front_search( string contig_i, string contig_j, string contig_tip, string contig_rev_tip, int index_i, int index_j ){
  // end_depth places the cursor at the depth of bp_added. Because this is at the front of contig_j, the matching segment extends further into contig_j and overlaps the old bp's that 
  //    were not covered in the last run
  end_depth = contigs[index_j].get_bp_added_fr();

  // limit to prevent out_of_bounds errors
  if( end_depth > contig_i.length() ){
    end_depth = contig_i.length();
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
  for( int k=contig_j.length()-end_depth; k<=contig_j.length()-tip_length; k++ ){
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

        //////////////////////
        // BACK END SEARCH //
        //////////////////////
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
  
  // make initial attempt to fuse contigs  
  // removed for master branch until algorithm can be adjusted
  //contig_fusion();

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
      t[i] = thread( thread_worker, ref(contigs), ref(qu), i );
    }

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
    //contig_fusion();
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
