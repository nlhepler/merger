
#include "bamfile.hpp"

#ifndef ARGPARSE_H
#define ARGPARSE_H

// program name
#define EXEC "merge"

// argument defaults
#define DEFAULT_MIN_OVERLAP 0
#define DEFAULT_MIN_READS 100
#define DEFAULT_TOL_GAPS true
#define DEFAULT_TOL_AMBIGS true
#define DEFAULT_BEGIN 0
#define DEFAULT_END 0
#define DEFAULT_STRIDE 25
#define DEFAULT_WINDOW_SIZE 125

class args_t
{
public:
    bamfile_t * bamin;
    bamfile_t * bamout;
    int min_overlap;
    int min_reads;
    int begin;
    int end;
    int window_size;
    int stride;
    bool tol_gaps;
    bool tol_ambigs;

    args_t( int, const char ** );
    ~args_t();
private:
    void parse_bamfile( const char *, const char * );
    void parse_minreads( const char * );
    void parse_begin( const char * );
    void parse_end( const char * );
    void parse_windowsize( const char * );
    void parse_stride( const char * );
};

#endif // ARGPARSE_H
