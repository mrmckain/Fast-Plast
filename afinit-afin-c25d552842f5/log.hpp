// Log class definitions for afin

#ifndef LOG_HPP
#define LOG_HPP

#include "process.hpp"
#include <string>
#include <fstream>

class Log{
  private:
    Log();
    static Log* p_Inst;
    static time_t timer;
    std::fstream log_fs;
    std::string get_time();

  public:
    static Log* Inst();
    void open_log( std::string logfile );
    void close_log();
    void log_it( std::string note );
};

#endif
