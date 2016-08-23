// $Author: benine $
// $Date$
// $Log$
// Contains the Contig class for afin


#ifndef CONTIG_HPP
#define CONTIG_HPP

#include "readlist.hpp"
#include "extension.hpp"
#include <vector>
#include <string>

/////////////////////////////////////////////////////\
// Contig Class: //////////////////////////////////////>
//////////////////////////////////////////////////////
// Contains contigs and holds matches to those contigs during processing
// Contains methods for extending contigs
class Contig{
  private:
    Readlist *reads;
    Extension *extension;
    std::string contig;
    std::string contig_id;

  public:
    Contig();
    Contig( Readlist *reads, std::string str, std::string id );
    ~Contig();

    // copy, move, =
    Contig( const Contig& rhs );
    Contig( Contig&& rhs );
    Contig& operator=( Contig rhs );
    friend void swap( Contig& data1, Contig& data2 );

    // return contig
    std::string get_sequence();

    // return contig_id
    std::string get_contig_id();

    // extend() performs loops iterations of create_extension with length extend_len of each extension, at each iteration the extension is added to contig, and uses contig_sub_len characters from the front or back of the contig, which end is determined by the boolean value back
    int extend( bool back );
};
//////////////////////
// End Contig Class //
//////////////////////

#endif
