// $Author: benine $
// $Date$
// $Log$
// Contains the Contig class for afin

#include "contig.hpp"
#include "process.hpp"

////////// Contig FUNCTIONS ////////////
Contig::Contig(){
  reads = 0;
  extension = 0;
  contig = "";
  contig_id = "";
};

Contig::Contig( Readlist *reads, std::string str, std::string id ) : reads(reads), contig(str), contig_id(id){
  extension = new Extension( reads, extend_len );
}

Contig::~Contig(){
  delete extension;
}

// copy, move, =
Contig::Contig( const Contig& rhs ) : reads(rhs.reads), contig(rhs.contig), contig_id(rhs.contig_id){
  extension = new Extension( reads, extend_len );
  extension->matches = rhs.extension->matches;
  extension->ATCG = rhs.extension->ATCG;
  extension->missed_bp = rhs.extension->missed_bp;
  extension->missed_bp_tot = rhs.extension->missed_bp_tot;
  extension->missed_bp_avg = rhs.extension->missed_bp_avg;
  extension->start = rhs.extension->start;
  extension->back = rhs.extension->back;
  extension->contig = rhs.extension->contig;
  extension->exten_seq = rhs.extension->exten_seq;
  extension->pos_mult = rhs.extension->pos_mult;
}

Contig::Contig( Contig&& rhs ) : Contig(){
  swap( *this, rhs );
}

Contig& Contig::operator=( Contig rhs ){
  swap( *this, rhs );
  return *this;
}

void swap( Contig& data1, Contig& data2 ){
  using std::swap;

  swap( data1.extension, data2.extension );
  swap( data1.reads, data2.reads );
  swap( data1.contig, data2.contig );
  swap( data1.contig_id, data2.contig_id );
}

// return contig_id
std::string Contig::get_sequence(){
  return contig;
}

std::string Contig::get_contig_id(){
  return contig_id;
}

// extend performs loops iterations of get_extension with length extend_len of each extension, at each iteration the extension is added to contig, and uses contig_sub_len characters from the front or back of the contig, which end is determined by the boolean value of back
int Contig::extend( bool back ){
  std::string exten_seq("");
  std::string contig_sub("");
  std::string contig_end("");

  // get extension sequence through get_extension
  if( back ){
    contig_sub = contig.substr( contig.length() - ( contig_sub_len ) );
    contig_end = " rear ";
  }
  else{
    contig_sub = contig.substr( 0, contig_sub_len );
    contig_end = " front ";
  }

  exten_seq = extension->get_extension( contig_sub, back );

  if( verbose ){
    std::lock_guard<std::mutex> lk_g(log_mut);
    Log::Inst()->log_it( contig_id + contig_end + "extension: " + exten_seq );
  }

  if( exten_seq.length() > 0 ){
    if( back ){
      contig.append( exten_seq );
    }
    else{
      contig.insert( 0, exten_seq );
    }
  }

  return exten_seq.length();
}
