
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "argparse.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

const char usage[] =
    "usage: " EXEC " [-h] "
    "[-o MIN_OVERLAP] "
    "[-r MIN_READS] "
    "[-g] [-a] "
    "(-B BAM_IN BAM_OUT)\n";

const char help_msg[] =
    "filter sequencing data using some simple heuristics\n"
    "\n"
    "required arguments:\n"
    "  -B BAM_IN BAM_OUT        BAM input and output files, respectively\n"
    "\n"
    "optional arguments:\n"
    "  -h, --help               show this help message and exit\n"
    "  -o MIN_OVERLAP           minimum overlap between two reads to merge them (default="
                                TO_STR( DEFAULT_MIN_OVERLAP ) ")\n"
    "  -r MIN_READS             minimum number of contributing reads to report a cluster (default="
                                TO_STR( DEFAULT_MIN_READS ) ")\n"
    "  -g                       don't tolerate gaps\n"
    "  -a                       don't tolerate ambigs\n";

inline
void help()
{
    fprintf( stderr, "%s\n%s", usage, help_msg );
    exit( 1 );
}

#define ERROR( msg, args... ) \
{ \
    fprintf( stderr, "%s" EXEC ": error: " msg "\n", usage , ##args ); \
    exit( 1 ); \
}

// args_t ----------------------------------------------------------------------------------------------------------- //

args_t::args_t( int argc, const char * argv[] ) :
    bamin( NULL ),
    bamout( NULL ),
    min_overlap( DEFAULT_MIN_OVERLAP ),
    min_reads( DEFAULT_MIN_READS ),
    tol_gaps( DEFAULT_TOL_GAPS ),
    tol_ambigs( DEFAULT_TOL_AMBIGS )
{
    int i;

    // skip arg[0], it's just the program name
    for ( i = 1; i < argc; ++i ) {
        const char * arg = argv[i];

        if ( arg[0] == '-' && arg[1] == '-' ) {
            if ( !strcmp( &arg[2], "help" ) ) help();
#if 0
            else if ( !strcmp( &arg[2], "fasta" ) ) parse_fasta( argv[++i] );
            else if ( !strcmp( &arg[2], "fastq" ) ) parse_fastq( argv[++i] );
            else if ( !strcmp( &arg[2], "output" ) ) parse_output( argv[++i] );
            else if ( !strcmp( &arg[2], "length" ) ) parse_length( argv[++i] );
#endif
            else
                ERROR( "unknown argument: %s", arg );
        }
        else if ( arg[0] == '-' ) {
            if ( !strcmp( &arg[1], "h" ) ) help();
            else if ( !strcmp( &arg[1], "B" ) ) {
                parse_bamfile( argv[i+1], argv[i+2] );
                i += 2;
            }
            else if ( !strcmp( &arg[1], "o" ) ) parse_minoverlap( argv[++i] );
            else if ( !strcmp( &arg[1], "r" ) ) parse_minreads( argv[++i] );
            else if ( !strcmp( &arg[1], "g" ) ) parse_tolgaps();
            else if ( !strcmp( &arg[1], "a" ) ) parse_tolambigs();
            else
                ERROR( "unknown argument: %s", arg );
        }
        else
            ERROR( "unknown argument: %s", arg );
    }

    if ( !bamin || !bamout )
        ERROR( "missing required argument -B BAM_IN BAM_OUT" );
}

args_t::~args_t()
{
    delete bamin;
    delete bamout;
}

void args_t::parse_bamfile( const char * input, const char * output )
{
    bamin = new bamfile_t( input, READ );
    bamout = new bamfile_t( output, WRITE );
}

void args_t::parse_minoverlap( const char * str )
{
    min_overlap = atoi( str );

    if ( min_overlap < 1 )
        ERROR( "minimum overlap must be an integer greater than 0, had: %s", str );
}

void args_t::parse_minreads( const char * str )
{
    min_reads = atoi( str );

    if ( min_reads < 1 )
        ERROR( "minimum reads must be an integer greater than 0, had: %s", str );
}

void args_t::parse_tolgaps()
{
    tol_gaps = false;
}

void args_t::parse_tolambigs()
{
    tol_ambigs = false;
}
