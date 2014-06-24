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
#include "print_time.hpp"
#include "revcomp.hpp"
#include "contig.hpp"
#include "process.hpp"

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
  cout << "READLIST[0]: " << readlist[0] << "  READLIST[last]: " << readlist[readlist.size()-1] << endl;
  merge_sort( readlist.begin(), readlist.end(), cmp_read );
  cout << "READLIST[0]: " << readlist[0] << "  READLIST[last]: " << readlist[readlist.size()-1] << endl;
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
void Process::add_reads( string fnames ){
  stringstream ss;
  ss.str( fnames );
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


// put contigs from contfile into contlist
void Process::add_contigs( string fnames ){
  stringstream ss;
  ss.str( fnames );
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
        contigs.push_back( Contig( buffer, contig_id, min_cov_init ));
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
    contigs.push_back( Contig( buffer, contig_id, min_cov_init ) );
  }
}

// return contig with index contig_ind
string Process::get_contig( int contig_ind ){
  return contigs[ contig_ind ].getContig();
}

// prints results to fasta file with outfile prefix and additional information is printed to a text based file with outfile prefix
void Process::print_to_outfile(){
  // open outfile
  ofstream outfile_fp( outfile+".fasta");

  // print out each line to the 
  for( int i=0; i<contigs.size(); i++ ){
    outfile_fp << ">" << contigs[i].get_contig_id() << "_" << outfile << endl;
    outfile_fp << get_contig(i) << endl;
  }

  outfile_fp.close();
} 

// prints notes to file as the program progresses
void Process::print_to_logfile( string note ){

}

// Compares the ends of the contigs with indices index_i and index_j which are related to the contigs from Process::contig_fusion()
// back indicates whether contig_i comes off the front or back of contig_j and changes the behavior of pos as follows:
//  1:  pos indicates the position in contig_j from which contig_i starts
//  0:  pos indicates the position in contig_j to which contig_i extends
// rev indicates whether contig_i is the reverse compliment of the original contig
// returns boolean value that indicates if a fusion was made
bool Process::contig_end_compare( int index_i, int index_j, int pos, bool back, bool rev ){
  string contig_j = get_contig( index_j );
  string contig_i = get_contig( index_i );
  string i_rev( "" );

  if( rev ){
    contig_i = revcomp( contig_i );
    i_rev = "_rev_comp";
  }

  if( back ){
    int len = contig_j.length() - pos;
    string contig_i_sub = contig_i.substr( 0, len );
    if( contig_j.compare( pos, contig_j.length(), contig_i ) == 0 ){
      // form fused contig and its id
      string fused( contig_j );
      fused.append( contig_i.substr( len, contig_i.length()-len ) );
      string fused_id( "fused("+contigs[index_j].get_contig_id()+"_||_"+contigs[index_i].get_contig_id()+i_rev+")" );

      // put pre-fused contigs in contigs_fused vector and post-fused contig in contigs vector
      contigs_fused.push_back( contigs[index_i] );
      contigs_fused.push_back( contigs[index_j] );
      contigs.push_back( Contig( fused, fused_id, min_cov_init ));
      
      // remove pre-fused vectors from contigs vector
      if( index_i>index_j ){
        contigs.erase( contigs.begin() + index_i );
        contigs.erase( contigs.begin() + index_j );
      }
      else{
        contigs.erase( contigs.begin() + index_j );
        contigs.erase( contigs.begin() + index_i );
      }

      return true;

    }
  }
  else{
    int len = contig_i.length() - pos;
    string contig_i_sub = contig_i.substr( len, pos );
    if( contig_j.compare( 0, pos, contig_i ) == 0 ){
      // form fused contig and its id
      string fused( contig_i );
      fused.append( contig_j.substr( len, contig_j.length()-len ) );
      string fused_id( "fused("+contigs[index_i].get_contig_id()+i_rev+"_||_"+contigs[index_j].get_contig_id()+")" );

      // put pre-fused contigs in contigs_fused vector and post-fused contig in contigs vector
      contigs_fused.push_back( contigs[index_i] );
      contigs_fused.push_back( contigs[index_j] );
      contigs.push_back( Contig( fused, fused_id, min_cov_init ));
      
      // remove pre-fused vectors from contigs vector
      if( index_i>index_j ){
        contigs.erase( contigs.begin() + index_i );
        contigs.erase( contigs.begin() + index_j );
      }
      else{
        contigs.erase( contigs.begin() + index_j );
        contigs.erase( contigs.begin() + index_i );
      }

      return true;
    }
  }

  return false;
}

// fuse contigs wherever possible
void Process::contig_fusion(){
  int end_length = 20;
  int end_depth = 100;
  
  // loop through each contig to get the end of the contig
  for( int i=0; i<contigs.size(); i++ ){
    string contig_i( get_contig( i ) );
    string contig_i_rev( revcomp( contig_i ) );

    // front end tip of the contig and that of the reverse compliment of the contig
    string contig_tip_fr( contig_i.substr( 0, end_length ) );
    string contig_rev_tip_fr( contig_i_rev.substr( 0, end_length ) );

    // rear end tip of the contig and that of the reverse compliment of the contig
    string contig_tip_rr( contig_i.substr( get_contig(i).length()-end_length, end_length ) );
    string contig_rev_tip_rr( contig_i_rev.substr( get_contig(i).length()-end_length, end_length ) );

    // loop through contigs and use this contig as the base to compare the end of contig_i and contig_i_rev to
    for( int j=0; j<contigs.size(); j++ ){
      if( j != i ){
        string contig_j( get_contig( j ) );
        
        //////////////////////
        // FRONT END SEARCH //
        //////////////////////
        end_depth = contigs[j].get_bp_added_fr();

        // limit to prevent out_of_bounds errors
        if( end_depth > contig_i.length() ){
          end_depth = contig_i.length();
        }

        //// search the front of the contig_j for the contig_tips created above
        for( int k=end_depth; k>=end_length; k-- ){
          if( contig_j.compare( k, end_length, contig_tip_rr ) == 0 ){
            if( contig_end_compare( i, j, k, false, false ) ){
              // if the two contigs match and are fused, adjustments must be made to the index markers so all contigs are checked
              if( i>j ){
                j--;
                i-=2;
              }
              else{
                j-=2;
                i--;
              }
              continue;
            }
          }
          
          if( contig_j.compare( k, end_length, contig_rev_tip_rr ) == 0 ){
            if( contig_end_compare( i, j, k, false, true ) ){
              // if the two contigs match and are fused, adjustments must be made to the index markers so all contigs are checked
              if( i>j ){
                j--;
                i-=2;
              }
              else{
                j-=2;
                i--;
              }
            }
          }
        }
        
        //////////////////////
        // BACK END SEARCH //
        //////////////////////
        end_depth = contigs[j].get_bp_added_rr();

        // limit to prevent out_of_bounds errors
        if( end_depth > contig_i.length() ){
          end_depth = contig_i.length();
        }

        //// search the back of the contig_j for the contig_tips created above
        for( int k=contig_j.length()-end_depth; k<contig_j.length()-end_length; k++ ){
          if( contig_j.compare( k, end_length, contig_tip_fr ) == 0 ){
            if( contig_end_compare( i, j, k, true, false ) ){
              // if the two contigs match and are fused, adjustments must be made to the index markers so all contigs are checked
              if( i>j ){
                j--;
                i-=2;
              }
              else{
                j-=2;
                i--;
              }
              continue;
            }
          }
          
          if( contig_j.compare( k, end_length, contig_rev_tip_fr ) == 0 ){
            if( contig_end_compare( i, j, k, true, true ) ){
              // if the two contigs match and are fused, adjustments must be made to the index markers so all contigs are checked
              if( i>j ){
                j--;
                i-=2;
              }
              else{
                j-=2;
                i--;
              }
            }
          }
        }
      }
    }
  }
}

//////////////////////////////
// END DEFINITIONS ///////////
//////////////////////////////
