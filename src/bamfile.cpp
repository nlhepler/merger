
#include <cstdlib>
#include <iostream>
#include <vector>

#include "merge.hpp"
#include "bamfile.hpp"

using std::cerr;
using std::endl;
using std::vector;


#define likely( x ) __builtin_expect( ( x ), 1 )
#define unlikely( x ) __builtin_expect( ( x ), 0 )


inline
uint32_t cigval( int op, int nop )
{
    return ( BAM_CIGAR_MASK & op ) | ( nop << BAM_CIGAR_SHIFT );
}


inline
void calloc_data( bam1_t * const buf )
{
    kroundup32( buf->m_data );
    buf->data = reinterpret_cast< uint8_t * >( calloc( buf->m_data, sizeof( uint8_t ) ) );

    if ( !buf->data ) {
        cerr << "memory allocation failure" << endl;
        exit( 1 );
    }
}


inline
void realloc_data( bam1_t * const buf )
{
    kroundup32( buf->m_data );
    uint8_t * data = reinterpret_cast< uint8_t * >( realloc( buf->data, buf->m_data * sizeof( uint8_t ) ) );

    if ( !data ) {
        cerr << "memory allocation error" << endl;
        exit( 1 );
    }

    buf->data = data;
}


bamfile_t::bamfile_t( const char * path, bam_mode_t mode ) :
    fp( NULL ),
    buf( NULL ),
    idx( NULL ),
    hdr( NULL )
{
    if ( strcmp( path, "-" ) ) {
        fp = bam_open( path, ( mode == READ ) ? "r" : "w" );
        if ( mode == READ )
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

inline
bool bam2aln( const bam1_t * const buf, aligned_t & aln, int begin = 0, int end = 0 )
{
    const uint32_t * cigar = bam1_cigar( buf );
    int i, j, col = buf->core.pos, rpos = col, aln_idx = 0, seq_idx = 0, ins = 0;

    aln.len = buf->core.l_qseq;
    aln.data = reinterpret_cast< pos_t * >( malloc( aln.len * sizeof( pos_t ) ) );
    aln.ncontrib = 1;

    if ( !aln.data )
        return false;

    // we'll be incrementing before first use
    --col;

    for ( i = 0; i < buf->core.n_cigar; ++i ) {
        const int op = cigar[ i ] & BAM_CIGAR_MASK;
        const int nop = cigar[ i ] >> BAM_CIGAR_SHIFT;

        rpos += nop;

        if ( op == BAM_CDEL ) {
            col += nop;
            ins = 0;
            continue;
        }

        for ( j = 0; j < nop; ++j ) {
            if ( op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF ) {
                ++col;
                ins = 0;
            }
            else if ( op == BAM_CINS )
                ++ins;

            if ( col < begin ) {
                ++seq_idx;
                continue;
            }
            else if ( end && col > end )
                goto done;

            aln.data[ aln_idx ].col = col;
            aln.data[ aln_idx ].ins = ins;
            aln.data[ aln_idx ].cov = 1;
            aln.data[ aln_idx ].nuc = bam1_seqi( bam1_seq( buf ), seq_idx );

            ++aln_idx;
            ++seq_idx;
        }
    }

done:
    if ( aln_idx )
        aln.lpos = aln.data[ 0 ].col;

    aln.len = aln_idx;
    aln.rpos = rpos;

    /* this is probably not strictly necessary
    bam_destroy1( buf );
    buf = bam_init1();
    */

    return true;
}


typedef struct {
    int begin;
    int end;
    vector< aligned_t > & reads;
} fetch_t;


static int fetch_func( const bam1_t * b, void * tmp )
{
    fetch_t * data = reinterpret_cast< fetch_t * >( tmp );
    aligned_t read = { .data = NULL };

    if ( bam2aln( b, read, data->begin, data->end ) )
        data->reads.push_back( read );

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


bool bamfile_t::next( aligned_t & aln )
{
    if ( fp->is_write )
        return false;

    if ( bam_read1( fp, buf ) >= 0 ) {
        return bam2aln( buf, aln );
    }

    return false;
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


bool bamfile_t::write(
    const char * const qname,
    aligned_t & aln,
    const int begin,
    const int end
    )
{
    if ( !fp->is_write || !aln.len )
        return false;

    int i = 0, j = 0;

    for( ; aln.data[ i ].col < begin && i < aln.len; ++i );

    const int beg_idx = i;
    int next_col = aln.data[ i ].col + 1;
    int op = aln.data[ i ].ins ? BAM_CINS : BAM_CMATCH;
    int nop = 1;

    buf->core.tid = 0;
    buf->core.pos = aln.lpos;
    buf->core.qual = 0xFF;
    buf->core.l_qname = 1 + strlen( qname );
    buf->m_data = buf->core.l_qname + 2*aln.len;
    // there are no mates
    buf->core.mtid = -1;
    buf->core.mpos = -1;
    buf->core.isize = 0;

    calloc_data( buf );

    strncpy( reinterpret_cast< char * >( buf->data ), qname, buf->core.l_qname );

    for ( ++i, j = 0; i < aln.len; ++i ) {
        const int col = aln.data[ i ].col;
        const int op_ = aln.data[ i ].ins ? BAM_CINS : BAM_CMATCH;

        if ( end && col > end )
            break;

        if ( col > next_col ) {
            if ( unlikely( op == BAM_CDEL ) )
                nop += col - next_col;
            else {
                bam1_cigar( buf )[ j++ ] = cigval( op, nop );
                nop = col - next_col;
                op = BAM_CDEL;
            }
        }

        if ( op_ == op )
            ++nop;
        else {
            bam1_cigar( buf )[ j++ ] = cigval( op, nop );
            nop = 1;
            op = op_;
        }

        next_col = col + 1;

        if ( reinterpret_cast< uint8_t * >( bam1_cigar( buf ) + j ) >=
                buf->data + buf->m_data ) {
            ++buf->m_data;
            realloc_data( buf );
        }
    }

    const int end_idx = i;
    buf->core.l_qseq = end_idx - beg_idx;
    bam1_cigar( buf )[ j++ ] = cigval( op, nop );
    buf->core.n_cigar = j;

    // qname-cigar-seq-qual-aux
    buf->l_aux = 0;
    buf->data_len = (
        buf->core.l_qname +
        4*buf->core.n_cigar +
        ( buf->core.l_qseq + 1 ) / 2 +
        buf->core.l_qseq
        );
    buf->m_data = buf->data_len;
    realloc_data( buf );

    for ( i = beg_idx, j = 0; i < end_idx; ++i, ++j ) {
        bam1_seq_seti( bam1_seq( buf ), j, aln.data[ i ].nuc );
    }

    // magic value to say "no quality scores present"
    bam1_qual( buf )[ 0 ] = 0xFF;
    buf->core.bin = bam_reg2bin( buf->core.pos, bam_calend( &buf->core, bam1_cigar( buf ) ) );

    if ( !bam_validate1( NULL, buf ) ) {
        cerr << "record failed validation" << endl;
        exit( 1 );
    }

    bam_write1( fp, buf );

    free( buf->data );
    memset( buf, 0, sizeof( bam1_t ) );

    return true;
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
