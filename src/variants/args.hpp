
#include "bamfile.hpp"

#ifndef ARGPARSE_H
#define ARGPARSE_H

// program name
#define EXEC "counter"

#define DEFAULT_CUTOFF 0.01

class args_t
{
public:
    bamfile::bamfile_t * bamin;
    double cutoff;

    args_t( int, const char ** );
    ~args_t();
private:
    void parse_bamfile( const char * );
    void parse_cutoff( const char * );
};

#endif // ARGPARSE_H
