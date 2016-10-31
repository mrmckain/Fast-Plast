// class for holding and processing the readlist

#include "readlist.hpp"
#include "process.hpp"
#include "revcomp.hpp"
#include "gzip.hpp"
#include <sstream>
#include <algorithm>

std::vector<std::string> Readlist::rlist;
std::vector<long int> Readlist::rc_reflist;

Readlist::Readlist( std::string readsfiles ){
  Log::Inst()->log_it( "Begin add_reads()" );
  add_reads( readsfiles );

  Log::Inst()->log_it( "Begin sort_reads()" );
  sort_reads();

  Log::Inst()->log_it( "Begin sort_rc()" );
  sort_rc();

  Log::Inst()->log_it( "Begin create_reads_range()" );
  create_read_range();
}

// compare function for sorting the reverse compliment reference list
bool Readlist::cmp_rc( const int ind1, const int ind2 ){
  std::string str1 = rlist[ind1].substr( rlist[ind1].length()-max_sort_char, max_sort_char );
  str1 = revcomp( str1 );
  std::string str2 = rlist[ind2].substr( rlist[ind2].length()-max_sort_char, max_sort_char );
  str2 = revcomp( str2 );

  return (str1.compare( str2 ) < 0 );
}

// compare function for sorting the initial read list
bool Readlist::cmp_read( const std::string str1, const std::string str2 ){
  return (str1.compare( 0, max_sort_char, str2.substr( 0, max_sort_char ) ) < 0 );
}

// mergesort template to be used with any datatype
template<class Iter, typename Order>
void Readlist::merge_sort( Iter first, Iter last, Order order ){
  if (last - first > 1){
    Iter middle = first + (last - first)/2;
    merge_sort( first, middle, order );
    merge_sort( middle, last, order );
    inplace_merge( first, middle, last, order );
  }
}

// uses a mergesort to sort the read list based on the first max_sort_char characters of each read
void Readlist::sort_reads(){
  merge_sort( rlist.begin(), rlist.end(), cmp_read );
}

// Produces list of references to the rlist sorted based on the first max_sort_char characters of the reverse_compliment of the
//  referenced read
//  Must be done after the rlist is sorted as it contains the locations in the rlist of the referenced read
//  Uses a mergesort
void Readlist::sort_rc(){
  rc_reflist.resize( rlist.size() );

  // initialize the rc_reflist with indices equal to their own as a beginning point
  for( int i=0; i<rc_reflist.size(); i++ ){
    rc_reflist[i]=i;
  }

  // sort
  merge_sort( rc_reflist.begin(), rc_reflist.end(), cmp_rc );
}

// creates a hash table within the process object that contains the ranges corresponding to equivalent first max_sort_char characters in the reads
// This increases the efficiency of searching
void Readlist::create_read_range(){
  std::string current( rlist[0].substr( 0, max_sort_char ) );
  long curr_start = 0;
  // cycle through the rlist
  for( long int i=1; i<rlist.size(); i++ ){
    if( rlist[i].compare( 0, max_sort_char, current ) != 0 ){
      // insert values into hash table, the ordered pair reflecting the range is increased by 1 to differentiate between a false search
      read_range.insert( { current, std::make_tuple( curr_start+1, i, 0, 0 )});
      curr_start = i;
      current = rlist[i].substr( 0, max_sort_char );
    }
  }

  // add the last entry
  read_range.insert( { current, std::make_tuple( curr_start+1, rlist.size(), 0, 0 )});

  std::string rc( rlist[0].substr( rlist[0].length() - max_sort_char, max_sort_char ) );
  current = revcomp( rc );
  curr_start = 0;

  // cycle through the rc_reflist to include these ranges in the hash as well
  for( long int i=1; i<rc_reflist.size(); i++ ){
    rc = rlist[i].substr( rlist[i].length() - max_sort_char, max_sort_char );
    rc = revcomp( rc );
    if( rc.compare( current ) != 0 ){
      // check if current already exists in hash
      if( std::get<0>( get_read_range(current) ) == -1 ){
        read_range.insert( { current, std::make_tuple( 0, 0, curr_start+1, i )});
      }
      else{
        std::get<2>( read_range[current] ) = curr_start+1;
        std::get<3>( read_range[current] ) = i;
      }
      curr_start = i;
      current = rc;
    }
  }

  // add the last entry
  if( std::get<0>( get_read_range(current) ) == -1 ){
    read_range.insert( { current, std::make_tuple( 0, 0, curr_start+1, rc_reflist.size() )});
  }
  else{
    std::get<2>( read_range[current] ) = curr_start+1;
    std::get<3>( read_range[current] ) = rc_reflist.size();
  }
}

// put reads from readfile into rlist
void Readlist::add_reads( std::string readsfiles ){
  std::stringstream ss;
  ss.str( readsfiles );
  std::string filename;
  std::string buffer("");
  std::string line("");
  char log_buff[1000];
  int line_count = 0;
  FileIO* read = NULL;

  while( getline( ss, filename, ',' )){
    Log::Inst()->log_it( "readfile: " + filename );

    if( filename.rfind( ".gz" ) == filename.length() - 3 ){
      read = new Gzip( filename );
    }
    else{
      read = new IO_Wrapper( filename );
    }

    // check what type of file it is
    if( (line = read->getline()) != "" ){
      if( line[0] == '@' ){
        // read in fastq reads
        do{
          line_count++;
          if( line[0] == '@' ){
            if( (line = read->getline()) != "" ){
              line_count++;
              if( !homopolymer_check( line ) ){
                if( line.length() > max_sort_char ){
                  rlist.push_back( line );
                }
              }
              if( ( line = read->getline()) != "" ){
                line_count++;
                if( line[0] == '+'){
                  if( ( line = read->getline()) != "" ){
                    line_count++;
                    continue;
                  }
                  else{
                    sprintf( log_buff, "Error reading fastq file. Scoring line missing. Line_num: %d Line: %s", line_count, line.c_str() );
                    Log::Inst()->log_it( log_buff );
                    fprintf( stderr, "Error reading fastq file. Scoring line missing. Line_num: %d Line: %s\n", line_count, line.c_str() );
                    break;
                  }
                }
                else{
                  sprintf( log_buff, "Error reading fastq file. '+' expected at this line. Line_num: %d Line: %s", line_count, line.c_str() );
                  Log::Inst()->log_it( log_buff );
                  fprintf( stderr, "Error reading fastq file. '+' expected at this line. Line_num: %d Line: %s\n", line_count, line.c_str() );
                  break;
                }
              }
              else{
                sprintf( log_buff, "Error reading fastq file. Line missing. Should start with '+'. Line_num: %d Line: %s", line_count, line.c_str() );
                Log::Inst()->log_it( log_buff );
                fprintf( stderr, "Error reading fastq file. Line missing. Should start with '+'. Line_num: %d Line: %s\n", line_count, line.c_str() );
                break;
              }
            }
            else{
              sprintf( log_buff, "Error reading fastq file. Sequence line missing. Line_num: %d Line: %s", line_count, line.c_str() );
              Log::Inst()->log_it( log_buff );
              fprintf( stderr, "Error reading fastq file. Sequence line missing. Line_num: %d Line: %s\n", line_count, line.c_str() );
              break;
            }
          }
          else{
            sprintf( log_buff, "Error reading fastq file. '@' expected at the beginning of this line. Line_num: %d Line: %s", line_count, line.c_str() );
            Log::Inst()->log_it( log_buff );
            fprintf( stderr, "Error reading fastq file. '@' expected at the beginning of this line. Line_num: %d Line: %s\n", line_count, line.c_str() );
            break;
          }
        } while( ( line = read->getline()) != "" );
      }
      else if( line[0] == '>' ){
        // read in reads to vector from fasta file
        do{
          if( line[0] == '>' && buffer.length() != 0 ){
            if( !homopolymer_check( buffer ) ){
              if( buffer.length() > max_sort_char ){
                rlist.push_back( buffer );
              }
            }
            buffer = "";
          }
          else if ( line[0] == '>' ) {
          }
          else{
            buffer += line;
          }
        } while( ( line = read->getline()) != "" );
      }
      else {
        Log::Inst()->log_it( "Error: Unexpected file type. Needs to be fasta or fastq file for input." );
        fprintf( stderr, "Error: Unexpected file type. Needs to be fasta or fastq file for input.\n" );
        break;
      }
    }
  }

  // insert last line into rlist
  if( buffer.length() != 0 ){
    if( !homopolymer_check( buffer ) ){
      if( buffer.length() > max_sort_char ){
        rlist.push_back( buffer );
      }
    }
  }

  if( rlist.size() == 0 ){
    Log::Inst()->log_it( "Error: No reads found in files supplied" );
    fprintf( stderr, "Error: No reads found in files supplied\n" );
    exit(1);
  }
}

// checks if string is a homopolymer
bool Readlist::homopolymer_check( std::string seq ){
  if( seq.compare( std::string(seq.length(), seq[0]) ) == 0 ){
  	return true;
  }
  return false;
}

// return the read range corresponding to the first max_sort_char characters of the read_seg passed
std::tuple<long,long,long,long> Readlist::get_read_range( std::string read_seg ){
  try{
    return read_range.at(read_seg);
  }
  catch( const std::out_of_range& e ){
    return std::make_tuple( -1,-1,-1,-1 );
  }
}

// return read from rlist with index ind
std::string Readlist::get_read( int ind ){
  return rlist[ind];
}

// return revcomp string of rlist member referenced by rc_reflist[ind]
std::string Readlist::get_rc_read( int ind ){
  return revcomp( rlist[rc_reflist[ind]] );
}
