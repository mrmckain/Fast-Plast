// gzip handling
//

#ifndef GZIP_HPP
#define GZIP_HPP

using namespace std;

// abstract class to provide seemless use of getline() to access either plain text or gzip files
class FileIO{
  public:
    virtual string getline() =0; 
};

// wrapper class for ifstream object to provide inheritance to FileIO
class IO_Wrapper : public FileIO{
  private:
    ifstream *fs;

  public:
    IO_Wrapper( string filename );

    // access string getline to return next line in ifstream object
    string getline();
};

class Gzip : public FileIO{
  private:
    const int BUFF_SIZE;
    int status;
    string buffer;
    z_stream zs; 
    FILE* fp;
    string filename;
    
    // append next set of data to buffer
    void fill_buffer();

  public:
    Gzip( string filename );

    /* report a zlib or i/o error */
    void zerr();

    // return the next line in the gzip file
    string getline();
};

#endif
