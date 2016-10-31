// gzip handling

#include "gzip.hpp"
#include <string>
#include <cstring>
#include <zlib.h>
#include "log.hpp"

IO_Wrapper::IO_Wrapper( std::string filename ){
  fs = new std::ifstream( filename );
}

// access string getline to return next line in ifstream object
std::string IO_Wrapper::getline(){
  std::string line("");
  std::getline( *fs, line );
  return line;
}

Gzip::Gzip( std::string filename ): BUFF_SIZE(16384), buffer(""), filename(filename){
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
  unsigned int char_in = 0;

  // read in data from gzipd file
  zs.avail_in = fread(in, 1, BUFF_SIZE, fp);

  // handle read errors
  if (ferror(fp)) {
    (void)inflateEnd(&zs);
    status = Z_ERRNO;
    return;
  }

  // handle data error
  if (zs.avail_in == 0){
    status = Z_DATA_ERROR;
    return;
  }

  // set next_in pointer to data from fread
  zs.next_in = in;
  do{
    zs.avail_out = BUFF_SIZE;
    zs.next_out = out;

    // decompress data
    status = inflate(&zs, Z_NO_FLUSH);

    // check for errors
    switch (status) {
      case Z_BUF_ERROR:
        status = Z_OK;
        return;
      case Z_NEED_DICT:
        status = Z_DATA_ERROR;     /* and fall through */
      case Z_DATA_ERROR:
      case Z_MEM_ERROR:
        (void)inflateEnd(&zs);
        zerr();
        return;
    }

    int have = BUFF_SIZE - zs.avail_out;
    char out_str[have+1];
    strncpy( out_str, (char*)out, have );
    out_str[have] = '\0';

    // append data to buffer
    buffer.append( out_str );
  }while( zs.avail_out == 0 );
}

// return the next line in the gzip file
std::string Gzip::getline(){
  size_t next_pos = 0;
  std::string line( "" );

  // check if buffer has a full line in it, if not, fill_buffer
  while((next_pos = buffer.find( '\n' )) == std::string::npos && status == Z_OK ){
    fill_buffer();
  }

  // if last line in file, store results for return
  if( status == Z_STREAM_END && next_pos == std::string::npos ){
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
