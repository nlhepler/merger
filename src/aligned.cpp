
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "bam.h"

#include "aligned.hpp"
#include "util.hpp"


using std::cerr;
using std::endl;
using std::make_pair;
using std::pair;
using std::string;
using std::vector;

using util::bits2nuc;


namespace aligned
{
    pos_t::pos_t(
            const int col,
            const op_t op,
            const int cov
            ) :
        col( col ),
        cov( cov ),
        op( op )
    {
    }


    pos_t::pos_t(
            const int col,
            const op_t op,
            const vector< pair< char, char > > data,
            const int cov
            ) :
        vector< pair< char, char > >( data ),
        col( col ),
        cov( cov ),
        op( op )
    {
    }


    void pos_t::get_qual( char * qual ) const
    {
        pos_t::const_iterator it;

        for ( it = begin(); it != end(); ++it )
            *( qual++ ) = it->second;
    }


    void pos_t::get_seq( char * str ) const
    {
        pos_t::const_iterator it;

        for ( it = begin(); it != end(); ++it )
            *( str++ ) = it->first;
    }


    void pos_t::get_seq( string & str ) const
    {
        pos_t::const_iterator it;

        str.clear();

        for ( it = begin(); it != end(); ++it )
            str.push_back( bits2nuc( it->first ) );
    }


    void pos_t::get_seq( vector< char > & vec ) const
    {
        pos_t::const_iterator it;

        vec.clear();

        for ( it = begin(); it != end(); ++it )
            vec.push_back( it->first );
    }


    inline
    uint32_t cigval( const op_t op, const int nop )
    {
        return ( BAM_CIGAR_MASK & uint32_t( op ) ) | ( nop << BAM_CIGAR_SHIFT );
    }


    aligned_t::aligned_t() :
        vector< pos_t >(),
        tid( 0 ),
        qual( 0xFF ),
        flag( 0 ),
        mtid( -1 ),
        mpos( -1 ),
        isize( 0 ),
        ncontrib( 1 )
    {
    }

    aligned_t::aligned_t( const bam1_t * const bam ) :
        vector< pos_t >(),
        tid( bam->core.tid ),
        qual( bam->core.qual ),
        flag( bam->core.flag ),
        mtid( bam->core.mtid ),
        mpos( bam->core.mpos ),
        isize( bam->core.isize ),
        name( string( bam1_qname( bam ) ) ),
        ncontrib( 0 )
    {
        int idx = 0, col = bam->core.pos - 1;
        const bool has_quals = bam1_qual( bam )[ 0 ] != 0xFF;

        for ( int i = 0; i < bam->core.n_cigar; ++i ) {
            const int nop = bam1_cigar( bam )[ i ] >> BAM_CIGAR_SHIFT;
            const int op = bam1_cigar( bam )[ i ] & BAM_CIGAR_MASK;

            if ( op == BAM_CDEL ) {
                col += nop;
            }
            else if ( op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF ) {
                for ( int j = 0; j < nop; ++j, ++idx ) {
                    pos_t pos( ++col, op_t( op ) );

                    pos.push_back(
                        make_pair(
                            bam1_seqi( bam1_seq( bam ), idx ),
                            has_quals ? bam1_qual( bam )[ idx ] : 0xFF
                            )
                        );

                    push_back( pos );
                }
            }
            else if ( op == BAM_CINS ) {
                pos_t pos( col, op_t( op ) );

                for ( int j = 0; j < nop; ++j, ++idx )
                    pos.push_back(
                        make_pair(
                            bam1_seqi( bam1_seq( bam ), idx ),
                            has_quals ? bam1_qual( bam )[ idx ] : 0xFF
                            )
                        );

                push_back( pos );
            }
            else {
                cerr << "unhandled CIGAR operation encountered" << endl;
                col += nop;
            }
        }
    }


    int aligned_t::lpos() const
    {
        aligned_t::const_iterator it = begin();
        return ( it == end() ) ? -1: it->col;
    }


    int aligned_t::rpos() const
    {
        aligned_t::const_reverse_iterator it = rbegin();
        return ( it == rend() ) ? -1 : it->col;
    }


    inline
    aligned_t::const_iterator prev( aligned_t::const_iterator & it )
    {
        return ( --it )++;
    }


    bool
    aligned_t::to_bam( bam1_t * const bam ) const
    {
        aligned_t::const_iterator it;
        int n_cigar, seq_len, cig_idx, seq_idx, nop;
        op_t op;

        if ( !size() )
            return false;

        memset( bam, '\0', sizeof( bam1_t ) );

        bam->core.pos = ( lpos() < 0 ) ? 0 : lpos();
        bam->core.tid = tid;

        if ( ncontrib ) {
            bam->core.qual = 0xFF;
            bam->core.flag = 0;
            bam->core.mtid = -1;
            bam->core.mpos = -1;
            bam->core.isize = 0;
        }
        else {
            bam->core.qual = qual;
            bam->core.flag = flag;
            bam->core.mtid = mtid;
            bam->core.mpos = mpos;
            bam->core.isize = isize;
        }

        it = begin();
        n_cigar = 0;
        seq_len = 0;
        op = it->op;

        for ( ++it; it != end(); ++it ) {
            if ( prev( it )->col + 1 < it->col ) {
                ++n_cigar; // for the last op
                ++n_cigar; // for being a deletion
                op = it->op;
            }
            else if ( it->op != op ) {
                ++n_cigar;
                op = it->op;
            }
            seq_len += it->size();
        }

        ++n_cigar; // for the last op

        bam->core.l_qname = name.size() + 1;
        bam->core.l_qseq = seq_len;
        bam->core.n_cigar = n_cigar;
        bam->l_aux = 0;
        bam->data_len = (
                bam->core.l_qname +
                4 * bam->core.n_cigar +
                ( bam->core.l_qseq + 1 ) / 2 +
                bam->core.l_qseq +
                bam->l_aux
                );
        bam->m_data = bam->data_len;
        kroundup32( bam->m_data );
        bam->data = reinterpret_cast< uint8_t * >( calloc( bam->m_data, sizeof( uint8_t ) ) );

        if ( !bam->data )
            goto error;

        memcpy( bam1_qname( bam ), name.c_str(), bam->core.l_qname );

        cig_idx = 0;
        seq_idx = 0;

        it = begin();
        nop = 1;
        op = it->op;

        for ( ++it; it != end(); ++it ) {
            vector< pair< char, char > >::const_iterator pit;

            if ( prev( it )->col + 1 < it->col ) {
                bam1_cigar( bam )[ cig_idx++ ] = cigval( op, nop );
                bam1_cigar( bam )[ cig_idx++ ] = cigval( DEL, it->col - prev( it )->col - 1 );
                nop = 1;
                op = it->op;
            }
            else if ( it->op != op ) {
                bam1_cigar( bam )[ cig_idx++ ] = cigval( op, nop );
                nop = 1;
                op = it->op;
            }
            else
                ++nop;

            // I want to use get_seq and get_qual here,
            // but I can't because of the seq is bit-packed
            for ( pit = it->begin(); pit != it->end(); ++pit ) {
                bam1_seq_seti( bam1_seq( bam ), seq_idx, pit->first );
                bam1_qual( bam )[ seq_idx++ ] = pit->second;
            }
        }

        bam1_cigar( bam )[ cig_idx++ ] = cigval( op, nop );

        if ( cig_idx != n_cigar ) {
            cerr << "something broke while assembling the CIGAR string" << endl;
            goto error;
        }

        if ( !bam_validate1( NULL, bam ) ) {
            cerr << "BAM record failed validation" << endl;
            goto error;
        }

        return true;

error:
        free( bam->data );
        memset( bam, '\0', sizeof( bam1_t ) );

        return false;
    }

    vector< pos_t > aligned_t::to_vector() const
    {
        aligned_t::const_iterator it;
        vector< pos_t > vec;

        for ( it = begin(); it != end(); ++it ) {
            vec.push_back( *it );
        }

        return vec;
    }
}
