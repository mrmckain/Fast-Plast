// $Author: benine $
// $Date$
// $Log$
// Contains the mismatch class for afin

#ifndef MISMATCH_H
#define MISMATCH_H

//////////////////////////////////////////////\
// Mismatch Class: ////////////////////////////>
//////////////////////////////////////////////
//
// Mismatch object, contains all classes, methods, data, and data references necessary for processing mismatches for contig fusion
// There will be only one Process object needed per iteration of this program
class Mismatch{
  private:
    double score;
    int length;
    int index_i;
    int index_j;
    int end_i;
    int end_j;
    
  public:
    Mismatch();

    Mismatch( double score, int length, int index_i, int index_j, int end_i, int end_j );

    // set mismatch score
    void set_score( double score );

    // set length
    void set_length( int length );

    // set index_i
    void set_index_i( int index );

    // set index_j
    void set_index_j( int index );

    // set index
    void set_indices( int index_i, int index_j );

    // set end_i
    void set_end_i( int end_i );

    // set end_j
    void set_end_j( int end_j );

    // return mismatch score
    double get_score();

    // return length
    int get_length();

    // return index i
    int get_index_i();

    // return index j
    int get_index_j();

    // return end_i
    int get_end_i();

    // return end_j
    int get_end_j();
};

#endif
