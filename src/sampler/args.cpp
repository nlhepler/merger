
/* argument parsing ------------------------------------------------------------------------------------------------- */

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "args.hpp"

// some crazy shit for stringifying preprocessor directives
#define STRIFY(x) #x
#define TO_STR(x) STRIFY(x)

const char usage[] =
    "usage: " EXEC " [-h] "
    "[-b BEGIN] "
    "[-e END] "
    "[-r MIN_READS] "
    "[-s STRIDE] "
    "[-w WINDOW_SIZE] "
    "(-B BAM_IN BAM_OUT)\n";

const char help_msg[] =
    "filter sequencing data using some simple heuristics\n"
    "\n"
    "required arguments:\n"
    "  -B BAM_IN BAM_OUT        BAM input and output files, respectively\n"
    "\n"
    "optional arguments:\n"
    "  -h, --help               show this help message and exit\n"
    "  -b BEGIN                 begin at position BEGIN in reference\n"
    "  -e END                   end at position END in reference\n"
    "  -r MIN_READS             minimum number of contributing reads to report a sample (default="
                                TO_STR( DEFAULT_MIN_READS ) ")\n"
    "  -s STRIDE                walk by STRIDE positions at a time (default="
                                TO_STR( DEFAULT_STRIDE ) ")\n"                
    "  -w WINDOW_SIZE           use windows of size WINDOW_SIZE (default="
                                TO_STR( DEFAULT_WINDOW_SIZE ) ")\n";

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
    begin( DEFAULT_BEGIN ),
    end( DEFAULT_END ),
    window_size( DEFAULT_WINDOW_SIZE ),
    stride( DEFAULT_STRIDE ),
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
            else if ( !strcmp( &arg[1], "b" ) ) parse_begin( argv[ ++i ] );
            else if ( !strcmp( &arg[1], "e" ) ) parse_end( argv[ ++i ] );
            else if ( !strcmp( &arg[1], "r" ) ) parse_minreads( argv[++i] );
            else if ( !strcmp( &arg[1], "s" ) ) parse_stride( argv[ ++i ] );
            else if ( !strcmp( &arg[1], "w" ) ) parse_windowsize( argv[ ++i ] );
            else
                ERROR( "unknown argument: %s", arg );
        }
        else
            ERROR( "unknown argument: %s", arg );
    }

    if ( !bamin || !bamout )
        ERROR( "missing required argument -B BAM_IN BAM_OUT" );

    if ( !end ) {
        if ( !bamin->hdr->target_len )
            ERROR( "no reference length information available" );

        end = bamin->hdr->target_len[ 0 ];
    }
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

void args_t::parse_minreads( const char * str )
{
    min_reads = atoi( str );

    if ( min_reads < 1 )
        ERROR( "minimum reads must be an integer greater than 0, had: %s", str );
}

void args_t::parse_begin( const char * str )
{
    begin = atoi( str );

    if ( begin < 0 )
        ERROR( "begin position must be a positive integer, had: %s", str );
}

void args_t::parse_end( const char * str )
{
    end = atoi( str );

    if ( end < 0 )
        ERROR( "end position must be a positive integer, had: %s", str );
}

void args_t::parse_stride( const char * str )
{
    stride = atoi( str );

    if ( stride < 1 )
        ERROR( "stride must be an integer greater than 0, had: %s", str );
}

void args_t::parse_windowsize( const char * str )
{
    window_size = atoi( str );

    if ( window_size < 1 )
        ERROR( "window size must be an integer greater than 0, had: %s", str );
}
