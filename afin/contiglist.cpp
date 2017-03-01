// Support for building and maintaing the contigs list

#include "contiglist.hpp"
#include "process.hpp"
#include <sstream>
#include <iostream>

Contiglist::Contiglist( Readlist *reads, std::string contigsfiles, std::string outfile ) : reads(reads), outfile(outfile){
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

// cycles through each contig and parses out the first section of the id
void Contiglist::parse_ids(){
  for( int i=0; i<contigs.size(); i++ ){
    std::string contig_id = contigs[i].get_contig_id();
    size_t pos = contig_id.find( "_", 5, 1 );

    if( pos != std::string::npos ){
      contig_id = contig_id.substr( 0, pos );
    }
  }
}

// put contigs from contfile into contlist
void Contiglist::add_contigs( std::string contigsfiles ){
  std::stringstream ss;
  ss.str( contigsfiles );
  std::string filename;
  std::string buffer("");
  std::string line("");
  std::string contig_id("");

  while( getline( ss, filename, ',' )){
    Log::Inst()->log_it( "contigsfile: " + filename );

    // open contig file
    std::ifstream cont( filename );

    // read in contig objects
    while( getline( cont, line ) ){
      if( line[0] == '>' && buffer.length() != 0 ){
        if( buffer.length() > 2*initial_trim + contig_sub_len ){
          buffer = buffer.substr( initial_trim, buffer.length() - 2*initial_trim );
          contigs.push_back( Contig( reads, buffer, contig_id ));
          if( verbose ){
            Log::Inst()->log_it( "Added contig: " + contig_id );
          }
        }
        else if( buffer.length() > contig_sub_len ){
          int trim = (buffer.length() - contig_sub_len) / 2;
          buffer = buffer.substr( trim, buffer.length() - 2*trim );
          contigs.push_back( Contig( reads, buffer, contig_id ));
          if( verbose ){
            Log::Inst()->log_it( "Added contig: " + contig_id );
          }
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
}

// print contigs
void Contiglist::output_contigs( int list_num, std::string file, std::string id_suffix ){
  std::vector<Contig> clist = contigs;

  std::cout << "output_contigs" << std::endl;
  if( list_num )
    clist = contigs_fused;

  // open outfile
  std::ofstream outfile_fp( file+".fa");

  // print out each line to the
  for( int i=0; i<clist.size(); i++ ){
    std::string seq = clist[i].get_sequence();
    outfile_fp << ">" << "contig" << i << "_len_" << seq.length() << "_" << id_suffix << std::endl;
    outfile_fp << seq << std::endl;
  }

  outfile_fp.close();
}

// prints results to fasta file with outfile prefix and additional information is printed to a text based file with outfile prefix
void Contiglist::create_final_fasta( int current_iteration ){
  // remove directories from outfile to form id_suffix if necessary
  size_t id_suffix_pos = outfile.find_last_of( "/" );
  std::string id_suffix = outfile;
  if( id_suffix_pos != std::string::npos ){
    id_suffix = id_suffix.substr( id_suffix_pos + 1 );
  }

  std::cout << "create_final_fasta" << std::endl;
  // print completed contigs to file
  output_contigs( 0, outfile + "_iter" + std::to_string(current_iteration), id_suffix );

  if( print_fused ){
    output_contigs( 1, outfile + "_iter" + std::to_string(current_iteration) + "_fused", id_suffix );
  }
}
