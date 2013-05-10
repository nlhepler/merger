
#include "bamfile.hpp"

#ifndef ARGPARSE_H
#define ARGPARSE_H

// program name
#define EXEC "merger"

// argument defaults
#define DEFAULT_MIN_OVERLAP 100
#define DEFAULT_MIN_READS 5
#define DEFAULT_TOL_GAPS true
#define DEFAULT_TOL_AMBIGS true

class args_t
{
public:
    bamfile::bamfile_t * bamin;
    bamfile::bamfile_t * bamout;
    bamfile::bamfile_t * bamdiscard;
    int min_overlap;
    int min_reads;
    bool tol_gaps;
    bool tol_ambigs;

    args_t( int, const char ** );
    ~args_t();
private:
    void parse_bamfile( const char *, const char * );
    void parse_bamdiscard( const char * );
    void parse_minoverlap( const char * );
    void parse_minreads( const char * );
    void parse_tolgaps();
    void parse_tolambigs();
};

#endif // ARGPARSE_H
