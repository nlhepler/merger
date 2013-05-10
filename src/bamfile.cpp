
#include <cstdlib>
#include <iostream>
#include <vector>

#include "aligned.hpp"
#include "bamfile.hpp"


using std::cerr;
using std::endl;
using std::vector;

using aligned::aligned_t;


#define likely( x ) __builtin_expect( ( x ), 1 )
#define unlikely( x ) __builtin_expect( ( x ), 0 )


namespace bamfile
{
    bamfile_t::bamfile_t( const char * path, bam_mode_t mode, bool index ) :
        fp( NULL ),
        buf( NULL ),
        idx( NULL ),
        hdr( NULL )
    {
        if ( strcmp( path, "-" ) ) {
            fp = bam_open( path, ( mode == READ ) ? "r" : "w" );
            if ( index && mode == READ )
                idx = bam_index_load( path );
        }
        else if ( mode == READ )
            fp = bam_dopen( fileno( stdin ), "r" );
        else
            fp = bam_dopen( fileno( stdout ), "w" );

        if ( fp == NULL ) {
            cerr << "failed to open BAM file: " << path << endl;
            exit( 1 );
        }

        hdr = ( mode == READ ) ? bam_header_read( fp ) : bam_header_init();
        buf = bam_init1();
        zero = bam_tell( fp );
    }

    bamfile_t::~bamfile_t()
    {
        if ( buf != NULL )
            bam_destroy1( buf );
        if ( hdr != NULL )
            bam_header_destroy( hdr );
        if ( idx != NULL )
            bam_index_destroy( idx );
        if ( fp != NULL )
            bam_close( fp );
    }


    typedef struct {
        int begin;
        int end;
        vector< aligned_t > & reads;
    } fetch_t;


    static int fetch_func( const bam1_t * const bam, void * tmp )
    {
        fetch_t * data = reinterpret_cast< fetch_t * >( tmp );
        
        data->reads.push_back( aligned_t( bam ) );

        return 0;
    }


    void bamfile_t::fetch( vector< aligned_t > & reads, int begin, int end, int tid )
    {
        fetch_t data = { .begin = begin, .end = end, .reads = reads };

        if ( !idx ) {
            cerr << "BAM index not found" << endl;
            exit( 1 );
        }

        bam_fetch( fp, idx, tid, begin, end, reinterpret_cast< void * >( &data ), fetch_func );
    }


    bool bamfile_t::next( bam1_t * const aln )
    {
        if ( fp->is_write )
            return false;

        if ( bam_read1( fp, aln ) >= 0 )
            return true;

        return false;
    }


    bool bamfile_t::write_header( const bam_header_t * hdr_ )
    {
        if ( !fp->is_write )
            return false;

        bam_header_write( fp, hdr_ ? hdr_ : hdr );

        return true;
    }


    bool bamfile_t::seek0()
    {
        return bam_seek( fp, zero, SEEK_SET ) == 0;
    }


    bool bamfile_t::write( const bam1_t * const bam )
    {
        if ( !bam_validate1( NULL, bam ) ) {
            cerr << "record failed validation" << endl;
            exit( 1 );
        }

        bam_write1( fp, bam );

        return true;
    }
}
