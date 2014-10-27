// $Author: benine $
// $Date$
// $Log$
// Contains the Contig class for afin

#include "contig.hpp"
#include "process.hpp"

using namespace std;

////////// Contig FUNCTIONS ////////////
Contig::Contig( Readlist *reads, string str, string id ) : reads(reads), contig(str), contig_id(id){
  cov = 0;
  doub_cov = false;
  extension = new Extension( reads, extend_len );
}

Contig::~Contig(){
  delete extension;
}

// return contig_id
string Contig::get_sequence(){
  return contig;
}

string Contig::get_contig_id(){
  return contig_id;
}

// return cov
double Contig::get_cov(){
  return cov;
}

// set cov
void Contig::set_cov( int cov ){
  this->cov = cov;
}

// extend performs loops iterations of get_extension with length extend_len of each extension, at each iteration the extension is added to contig, and uses contig_sub_len characters from the front or back of the contig, which end is determined by the boolean value of back
void Contig::extend( bool back ){
  string exten_seq("");
  string contig_sub("");
 
  // skip over any contigs that present at least double coverage
  if( doub_cov ){
    return;
  }

  // get extension sequence through get_extension
  if( back ){
    contig_sub = contig.substr( contig.length() - ( contig_sub_len ) );
  }
  else{
    contig_sub = contig.substr( 0, contig_sub_len );
  }

  exten_seq = extension->get_extension( contig_sub, back );
  
  if( exten_seq.length() == 0 ){
    return;
  }

  if( back ){
    contig.append( exten_seq );
  }
  else{
    contig.insert( 0, exten_seq );
  }
}

// set doub_cov var
void Contig::set_doub_cov( bool doub_cov ){
  this->doub_cov = doub_cov;
}
