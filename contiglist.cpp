// Support for building and maintaing the contigs list

#include "contiglist.hpp"
#include "process.hpp"
#include <sstream>

using namespace std;

Contiglist::Contiglist( Readlist *reads, string contigsfiles, string outfile ) : reads(reads), outfile(outfile){
  Log::Inst()->log_it( "Begin add_contigs()" );
  add_contigs( contigsfiles );
}

// return contig at index ind
Contig Contiglist::get_contig( int ind ){
  return contigs[ind];
}

// return contig at index ind
Contig *Contiglist::get_contig_ref( int ind ){
  return (Contig*)&contigs[ind];
}

// return conitg list size
int Contiglist::get_list_size(){
  return (int)contigs.size();
}

// remove contig at index ind
void Contiglist::remove_contig( int ind ){
  contigs.erase( contigs.begin() + ind );
}

// append contig to contigs
void Contiglist::append_contig( int list_num, Contig cont ){
  if( list_num ){
    contigs_fused.push_back( cont );
  }
  else{
    contigs.push_back( cont );
  }
}

// parses the cov value from the contig_id and passes the result back as a double
// If cov is not in the header, return 1
double Contiglist::parse_cov( string contig_id ){
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
void Contiglist::contig_cov(){
  double cov_total = 0;
  double cov = 0;
  int total_contigs = contigs.size();
 
  // protect against division by 0 and the end of the world
  if( total_contigs == 0 ){
    return;
  }

  for( int i=0; i<total_contigs; i++ ){
    // get contig id, parse out cov_##, add to the total of all cov values
    cov = parse_cov( contigs[i].get_contig_id() );
    contigs[i].set_cov( cov );
    cov_total += cov;
  }

  // find average of all cov values
  double cov_avg = cov_total / total_contigs;
  
  for( int i=0; i<total_contigs; i++ ){
    // get contig_id, parse out cov_## and compare this value to the avg
    cov = contigs[i].get_cov();

    if( cov > cov_avg * 2.0 ){
      // set the value of doub_cov to false to indicate ignoring when extending contigs 
      contigs[i].set_doub_cov( true );
    }
  }
}

// cycles through each contig and parses out the first section of the id
void Contiglist::parse_ids(){
  for( int i=0; i<contigs.size(); i++ ){
    string contig_id = contigs[i].get_contig_id(); 
    size_t pos = contig_id.find( "_", 5, 1 );

    if( pos != string::npos ){
      contig_id = contig_id.substr( 0, pos );
    }
  }
}

// put contigs from contfile into contlist
void Contiglist::add_contigs( string contigsfiles ){
  stringstream ss;
  ss.str( contigsfiles );
  string filename;
  string buffer("");
  string line("");
  string contig_id("");
  
  while( getline( ss, filename, ',' )){
    Log::Inst()->log_it( "contigsfile: " + filename );
 
    // open contig file
    ifstream cont( filename );

    // read in contig objects
    while( getline( cont, line ) ){
      if( line[0] == '>' && buffer.length() != 0 ){
        if( buffer.length() > 2*initial_trim + contig_sub_len ){
          buffer = buffer.substr( initial_trim, buffer.length() - 2*initial_trim );
          contigs.push_back( Contig( reads, buffer, contig_id ));
        }
        else if( buffer.length() > contig_sub_len ){
          int trim = (buffer.length() - contig_sub_len) / 2;
          buffer = buffer.substr( trim, buffer.length() - 2*trim );
          contigs.push_back( Contig( reads, buffer, contig_id ));
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
    contigs.push_back( Contig( reads, buffer, contig_id ) );
  }

  contig_cov();
  parse_ids();
}

// print contigs
void Contiglist::output_contigs( int list_num, string file, string id_suffix ){
  vector<Contig> clist = contigs;  

  if( list_num )
    clist = contigs_fused;
  
  // open outfile
  ofstream outfile_fp( file+".fa");

  // print out each line to the
  for( int i=0; i<clist.size(); i++ ){
    outfile_fp << ">" << clist[i].get_contig_id() << "_" << id_suffix << endl;
    outfile_fp << clist[i].get_sequence() << endl;
  }

  outfile_fp.close();
}

// prints results to fasta file with outfile prefix and additional information is printed to a text based file with outfile prefix
void Contiglist::create_final_fasta(){
  // remove directories from outfile to form id_suffix if necessary
  size_t id_suffix_pos = outfile.find_last_of( "/" );
  string id_suffix = outfile;
  if( id_suffix_pos != string::npos ){
    id_suffix = id_suffix.substr( id_suffix_pos + 1 );
  }

  // print completed contigs to file
  output_contigs( 0, outfile, id_suffix );

  if( print_fused ){
    output_contigs( 1, outfile+"_fused", id_suffix );
  }

  if( log_output || screen_output ){
    Log::Inst()->close_log();
  }
} 
