
#include <cstdlib>
#include <iostream>

#include "merge.hpp"
#include "bamfile.hpp"

using std::cerr;
using std::endl;


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
    buf->data = reinterpret_cast< uint8_t * >( realloc( buf->data, buf->m_data * sizeof( uint8_t ) ) );

    if ( !buf->data ) {
        cerr << "memory allocation error" << endl;
        exit( 1 );
    }
}


bamfile_t::bamfile_t( const char * path, bam_mode_t mode ) :
    fp( NULL ),
    buf( NULL ),
    hdr( NULL )
{
    if ( strcmp( path, "-" ) )
        fp = bam_open( path, ( mode == READ ) ? "r" : "w" );
    else if ( mode == READ )
        fp = bam_dopen( fileno( stdin ), "r" );
    else
        fp = bam_dopen( fileno( stdout ), "w" );

    if ( fp != NULL ) {
        hdr = ( mode == READ ) ? bam_header_read( fp ) : bam_header_init();
        buf = bam_init1();
    }
}

bamfile_t::~bamfile_t()
{
    if ( buf != NULL )
        bam_destroy1( buf );
    if ( hdr != NULL )
        bam_header_destroy( hdr );
    if ( fp != NULL )
        bam_close( fp );
}

bool bamfile_t::next( aligned_t & aln )
{
    if ( fp->is_write )
        return false;

    if ( bam_read1( fp, buf ) >= 0 ) {
        const uint32_t * cigar = bam1_cigar( buf );
        int i, j, col = buf->core.pos, rpos = col, idx = 0, ins = 0;

        aln.len = buf->core.l_qseq;
        aln.lpos = col;
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

            for ( j = 0; j < nop; ++j ) {
                if ( op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF ) {
                    col += 1;
                    ins = 0;
                }
                else if ( op == BAM_CINS )
                    ins += 1;
                else if ( op == BAM_CDEL ) {
                    col += 1;
                    ins = 0;
                    continue;
                }

                aln.data[ idx ].col = col;
                aln.data[ idx ].ins = ins;
                aln.data[ idx ].cov = 1;
                aln.data[ idx ].nuc = bam1_seqi( bam1_seq( buf ), idx );

                ++idx;
            }
        }

        aln.rpos = rpos;

        bam_destroy1( buf );
        buf = bam_init1();

        return true;
    }

    return false;
}

bool bamfile_t::write_header( const bam_header_t * hdr_ )
{
    if ( !fp->is_write )
        return false;

    bam_header_write( fp, hdr_ ? hdr_ : hdr );

    return true;
}

bool bamfile_t::write( const char * const qname, aligned_t & aln )
{
    if ( !fp->is_write || !aln.len )
        return false;

    int i, j = 0;
    int next_col = aln.data[ 0 ].col + 1;
    int op = aln.data[ 0 ].ins ? BAM_CINS : BAM_CMATCH;
    int nop = 1;

    buf->core.tid = 0;
    buf->core.pos = aln.lpos;
    buf->core.qual = 0xFF;
    buf->core.l_qname = 1 + strlen( qname );
    buf->m_data = buf->core.l_qname + 2*aln.len;
    buf->core.l_qseq = aln.len;
    // there are no mates
    buf->core.mtid = -1;
    buf->core.mpos = -1;
    buf->core.isize = 0;

    calloc_data( buf );

    strncpy( reinterpret_cast< char * >( buf->data ), qname, buf->core.l_qname );

    for ( i = 1, j = 0; i < aln.len; ++i ) {
        const int col = aln.data[ i ].col;
        const int op_ = aln.data[ i ].ins ? BAM_CINS : BAM_CMATCH;

        if ( col > next_col ) {
            if ( op == BAM_CDEL )
                nop += col - next_col;
            else {
                bam1_cigar( buf )[ j++ ] = ( BAM_CIGAR_MASK & op ) | ( nop << BAM_CIGAR_SHIFT );
                nop = col - next_col;
                op = BAM_CDEL;
            }
        }

        next_col = col + 1;

        if ( op_ == op )
            nop += 1;
        else {
            bam1_cigar( buf )[ j++ ] = ( BAM_CIGAR_MASK & op ) | ( nop << BAM_CIGAR_SHIFT );
            nop = 1;
            op = op_;
        }

        if ( reinterpret_cast< uint8_t * >( bam1_cigar( buf ) + j ) >=
                buf->data + buf->m_data ) {
            ++buf->m_data;
            realloc_data( buf );
        }
    }

    bam1_cigar( buf )[ j++ ] = ( BAM_CIGAR_MASK & op ) | ( nop << BAM_CIGAR_SHIFT );

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

    for ( i = 0; i < aln.len; ++i ) {
        bam1_seq_seti( bam1_seq( buf ), i, aln.data[ i ].nuc );
    }

    // magic value to say "no quality scores present"
    bam1_qual( buf )[ 0 ] = 0xFF;
    buf->core.bin = bam_reg2bin( buf->core.pos, bam_calend( &buf->core, bam1_cigar( buf ) ) );

    if ( !bam_validate1( NULL, buf ) ) {
        cerr << "record failed validation" << endl;
        exit( 1 );
    }

    bam_write1( fp, buf );

    memset( buf, 0, sizeof( bam1_t ) );

    return false;
}
