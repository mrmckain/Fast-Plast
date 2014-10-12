// $Author: benine $
// $Date$
// $Log$
// This contains the read class for afin

#ifndef READ_HPP
#define READ_HPP

using namespace std;

//////////////////////////////////////\
// Read Class:  ///////////////////////>
//////////////////////////////////////
//  Contains one read and information about that read in reference to the piece of the contig it matched to
class Read{
  private:
    int start;      // position of read where match starts in reference to the contig it matches
                    // for instance:
                    //    match: aaaattta
                    //    read:  tttatt
                    //    start: 4
    int misses;    // counts the number of characters that don't match other reads
    string read;    // contains the read or rev compliment, whichever is found
    bool revcomp;   // indicates if the read is a reverse compliment of the matched string

  public:
    // constructor.. default revcomp will be false
    Read( string& read, int match, bool revcomp );

    Read( string& read, int match );
    
    // returns character at pos where pos is the nth postion of the match
    char get_pos( int pos, bool back );
};

////////////////////
// End Read Class //
////////////////////

#endif
