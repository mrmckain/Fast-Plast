// gzip handling

#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include <zlib.h>
#include "gzip.hpp"
#include "log.hpp"

using namespace std;

IO_Wrapper::IO_Wrapper( string filename ){
  fs = new ifstream( filename );
}

// access string getline to return next line in ifstream object
string IO_Wrapper::getline(){
  string line("");
  std::getline( *fs, line );
  return line;
}

Gzip::Gzip( string filename ): BUFF_SIZE(16384), buffer(""), filename(filename){
  fp = fopen( filename.c_str(), "rb" );
  
  /* allocate inflate state */
  zs.zalloc = Z_NULL;
  zs.zfree = Z_NULL;
  zs.opaque = Z_NULL;
  zs.avail_in = 0;
  zs.avail_out = BUFF_SIZE;
  zs.next_in = Z_NULL;
  status = inflateInit2(&zs, 32+MAX_WBITS);

  if (status != Z_OK)
    zerr();
}

/* report a zlib or i/o error */
void Gzip::zerr(){
  char log_buff[1000];

  switch (status) {
    case Z_ERRNO:
      fprintf( stderr, "zlib: error reading file: %s\n", filename.c_str() );
      sprintf( log_buff, "zlib: error reading file: %s", filename.c_str() );
      Log::Inst()->log_it( log_buff );
      break;
    case Z_STREAM_ERROR:
      fprintf( stderr, "zlib: invalid compression level: %s\n", filename.c_str() );
      sprintf( log_buff, "zlib: invalid compression level: %s", filename.c_str() );
      Log::Inst()->log_it( log_buff );
      break;
    case Z_DATA_ERROR:
      fprintf( stderr, "zlib: invalid or incomplete deflate data: %s\n", filename.c_str() );
      sprintf( log_buff, "zlib: invalid or incomplete deflate data: %s", filename.c_str() );
      Log::Inst()->log_it( log_buff );
      break;
    case Z_MEM_ERROR:
      fprintf( stderr, "zlib: out of memory: %s\n", filename.c_str() );
      sprintf( log_buff, "zlib: out of memory: %s", filename.c_str() );
      Log::Inst()->log_it( log_buff );
      break;
    case Z_VERSION_ERROR:
      fprintf( stderr, "zlib: version mismatch: %s\n", filename.c_str() );
      sprintf( log_buff, "zlib: version mismatch: %s", filename.c_str() );
      Log::Inst()->log_it( log_buff );
  }
}

// append next set of data to buffer
void Gzip::fill_buffer(){
  unsigned char in[BUFF_SIZE];
  unsigned char out[BUFF_SIZE];
  int next_pos = 0;
  int new_data = 0;

  cout << "avail_out: " << zs.avail_out << endl;
  cout << "avail_in: " << zs.avail_in << endl;

  // read in data from gzipd file
  new_data = fread(in, 1, BUFF_SIZE-zs.avail_in, fp);
  zs.avail_in += new_data;

  cout << "ferror status: " << status << endl;

  // handle read errors
  if (ferror(fp)) {
    (void)inflateEnd(&zs);
    status = Z_ERRNO;
    return;
  }

  // handle data error
  if (zs.avail_in == 0){
    status = Z_DATA_ERROR;
    cout << "DATA_ERROR" << endl;
    return;
  }

  if( new_data == zs.avail_in ){
    memcpy( zs.next_in + zs.avail_in - new_data, &in, new_data );
  }
  else{
    zs.avail_in += new_data;


  cout << "predecompress status: " << status << endl;

  // prep z_stream members
  zs.avail_out = BUFF_SIZE;
  zs.next_out = out;

  // decompress data
  status = inflate(&zs, Z_NO_FLUSH);
  cout << "decompress status: " << status << endl;
  cout << "avail_out: " << zs.avail_out << endl;
  cout << "avail_in: " << zs.avail_in << endl;

  // check for errors
  switch (status) {
    case Z_NEED_DICT:
      status = Z_DATA_ERROR;     /* and fall through */
    case Z_DATA_ERROR:
    case Z_MEM_ERROR:
      (void)inflateEnd(&zs);
      cout << "ERROR" << endl;
      return;
  }
  
  stringstream s;

  s << out;

  // append data to buffer
  buffer.append( s.str() );
}

// return the next line in the gzip file
string Gzip::getline(){
  size_t next_pos = 0;
  string line( "" );

  // check if buffer has a full line in it, if not, fill_buffer
  while((next_pos = buffer.find( '\n' )) == string::npos && status == Z_OK ){
    cout << "fill_buffer()" << endl;
    fill_buffer();
    cout << "status: " << status << endl;
  }

  // if last line in file, store results for return
  if( status == Z_STREAM_END && next_pos == string::npos ){
    line = buffer;
    buffer = "";

    (void)inflateEnd(&zs);
  }
  // protect against errors
  else if( status != Z_OK && status != Z_STREAM_END ){
    (void)inflateEnd(&zs);
    zerr();
    return "";
  }
  // parse out next line in buffer
  else{
    line = buffer.substr( 0, next_pos );
    buffer.erase( 0, next_pos+1 );
  }

  return line;
}
