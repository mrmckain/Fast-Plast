// $Author: benine $
// $Date$
// $Log$
// Contains the mismatch class for afin

#ifndef MISMATCH_H
#define MISMATCH_H

using namespace std;

//////////////////////////////////////////////\
// Mismatch Class: ////////////////////////////>
//////////////////////////////////////////////
//
// Mismatch object, contains all classes, methods, data, and data references necessary for processing mismatches for contig fusion
// There will be only one Process object needed per iteration of this program
class Mismatch{
  private:
    int id;
    double score;
    int length;
    
    // orientation is an integer 0-3 and corresponds as follows:
    //    0:  i to j
    //    1:  i to j_rev
    //    2:  j to i
    //    3:  j_rev to i
    int orientation;
    int match_index_i;
    int match_index_j;
    bool rev;
    
  public:
    Mismatch();

    Mismatch( int id, double score, int length, int orientation, int match_index_i, int match_index_j );

    Mismatch( int id, int orientation, int match_index_i, int match_index_j );

    // set id
    void set_id( int id );

    // set mismatch score
    void set_score( double score );

    // set length
    void set_length( int length );

    // set orientation
    void set_orientation( int orientation );

    // set index
    void set_indices( int match_index_i, int match_index_j );

    // set rev
    void set_rev( bool rev );

    // return id
    int get_id();

    // return mismatch score
    double get_score();

    // return length
    int get_length();

    // return orientation
    int get_orientation();

    // return index i
    int get_index_i();

    // return index j
    int get_index_j();
    
    // return rev
    bool get_rev();
};

#endif
