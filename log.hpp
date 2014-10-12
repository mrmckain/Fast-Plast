// Log class definitions for afin

#ifndef LOG_HPP
#define LOG_HPP

#include <fstream>
#include <ctime>
#include <vector>
#include "process.hpp"

class Log{
  private:
    Log();
    static Log* p_Inst;
    static time_t timer;
    fstream log_fs;
    string get_time();

  public:
    static Log* Inst();
    void open_log( string logfile );
    void close_log();
    void log_it( string note );
};

#endif
