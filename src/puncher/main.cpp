
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <list>
#include <utility>
#include <vector>

#include "bam.h"

#include "aligned.h"
#include "args.hpp"
#include "bamfile.hpp"
#include "coverage.hpp"
#include "math.hpp"
#include "rateclass_em.hpp"
#include "util.hpp"


using std::cerr;
using std::endl;
using std::exp;
using std::list;
using std::log;
using std::make_pair;
using std::map;
using std::pair;
using std::vector;


typedef list< col_t >::const_iterator cov_citer;
typedef map< vector< char >, int >::const_iterator obs_citer;
typedef map< vector< char >, int >::iterator obs_iter;
typedef vector< col_t >::const_iterator var_citer;
typedef vector< triple_t >::const_iterator read_citer;


inline
uint32_t cigval( int op, int nop )
{
    return ( BAM_CIGAR_MASK & op ) | ( nop << BAM_CIGAR_SHIFT );
}


inline
void realloc_data( bam1_t * const bam )
{
    kroundup32( bam->m_data );
    uint8_t * data = reinterpret_cast< uint8_t * >( realloc( bam->data, bam->m_data * sizeof( uint8_t ) ) );

    if ( !data ) {
        cerr << "memory allocation error" << endl;
        exit( 1 );
    }

    bam->data = data;
}


bam1_t * punchout_read(
        const bam1_t * const in_bam,
        const vector< col_t > & variants,
        const vector< triple_t > & read
        )
{
    read_citer rit = read.begin();
    var_citer vit = variants.begin();
    vector< pair< bool, unsigned > > keep;
    vector< int > cols;

    bam1_t * const out_bam = bam_init1();

    if ( !out_bam ) {
        cerr << "memory allocation error" << endl;
        exit( 1 );
    }

    out_bam->core = in_bam->core;
    out_bam->core.l_qseq = 0;

    // build a vector of which positions we're keeping,
    // and how long our new seq will be
    for ( ; rit != read.end() && vit != variants.end(); ) {

        /*
        cerr << "r: ( " << rit->col << ", " << rit->ins << ", ";
        for ( unsigned j = 0; j < rit->elem.size(); ++j )
            cerr << bits2nuc( rit->elem[ j ] );
        cerr << " )" << endl;
        cerr << "v: ( " << vit->col << ", " << vit->ins << " )" << endl;
        */

        if ( rit->col == vit->col && rit->ins == vit->ins ) {
            obs_citer it = vit->obs.find( rit->elem );

            if ( it == vit->obs.end() ) {
                cerr << "unknown variant observed, which is weird...( 1 )" << endl;
                exit( 1 );
            }

            if ( it->second ) {
                keep.push_back( make_pair( true, it->first.size() ) );
                cols.push_back( rit->col );
                out_bam->core.l_qseq += it->first.size();
            }
            else
                keep.push_back( make_pair( false, it->first.size() ) );

            ++rit;
            ++vit;
        }
        // matching column in variants, but insertion in read
        else if ( rit->col == vit->col && rit->ins ) {
            ++vit;
        }
        // matching column in read, but insertion in variants:
        // we got ahead in read somehow, which should never happen
        else if ( rit->col == vit->col && vit->ins ) {
            cerr << "unknown variant observed, which is weird...( 2 )" << endl;
            exit( 1 );
        }
        // read is behind variants for some reason, which should never happen
        else if ( rit->col < vit->col ) {
            cerr << "unknown variant observed, which is weird...( 3 )" << endl;
            exit( 1 );
        }
        else // ( vit->col < rit->col )
            ++vit;
    }

    cerr << "seqlen: " << out_bam->core.l_qseq << endl;

    out_bam->data_len = in_bam->data_len;
    out_bam->m_data = out_bam->data_len;
    out_bam->data = NULL;

    realloc_data( out_bam );

    // copy the read name
    memcpy( bam1_qname( out_bam ), bam1_qname( in_bam ), in_bam->core.l_qname );

    // copy the cigar string
    {
        const uint32_t * cigar = bam1_cigar( in_bam );
        int in_idx = 0, out_idx = 0;

        for ( int i = 0; i < in_bam->core.n_cigar; ++i ) {
            const int op = cigar[ i ] & BAM_CIGAR_MASK;
            const int nop = cigar[ i ] >> BAM_CIGAR_SHIFT;

            if ( op == BAM_CDEL )
                bam1_cigar( out_bam )[ out_idx++ ] = cigar[ i ];
            else if ( op == BAM_CINS ) {
                if ( keep[ in_idx ].first )
                    bam1_cigar( out_bam )[ out_idx++ ] = cigval( op, keep[ in_idx ].second );
                else
                    bam1_cigar( out_bam )[ out_idx++ ] = cigval( BAM_CDEL, keep[ in_idx ].second );
                ++in_idx;
            }
            else if ( op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF ) {
                bool k = keep[ in_idx ].first;
                int n = keep[ in_idx ].second;

                ++in_idx;

                for ( int j = n; j < nop; ++in_idx ) {
                    bool k_ = keep[ in_idx ].first;
                    if ( k_ == k )
                        n += keep[ in_idx ].second;
                    else {
                        bam1_cigar( out_bam )[ out_idx++ ] = cigval( k ? op : BAM_CDEL, n );
                        k = k_;
                        n = keep[ in_idx ].second;
                    }
                    j += keep[ in_idx ].second;
                }

                bam1_cigar( out_bam )[ out_idx++ ] = cigval( k ? op : BAM_CDEL, n );
            }

            if ( out_bam->core.l_qname + out_idx >= out_bam->m_data ) {
                ++out_bam->m_data;
                realloc_data( out_bam );
            }
        }

        out_bam->core.n_cigar = out_idx;
    }

    out_bam->l_aux = in_bam->l_aux;

    out_bam->data_len = (
            out_bam->core.l_qname +
            out_bam->core.n_cigar +
            ( out_bam->core.l_qseq + 1 ) / 2 +
            out_bam->core.l_qseq +
            out_bam->l_aux
            );

    out_bam->m_data = out_bam->data_len;

    realloc_data( out_bam );

    memcpy( bam1_aux( out_bam ), bam1_aux( in_bam ), in_bam->l_aux );

    // copy the sequence and quality scores
    {
        int in_idx = 0, out_idx = 0;

        for ( unsigned i = 0; i < keep.size(); ++i ) {
            if ( keep[ i ].first )
                for ( unsigned k = 0; k < keep[ i ].second; ++k ) {
                    bam1_seq_seti( bam1_seq( out_bam ), out_idx, bam1_seqi( bam1_seq( in_bam ), in_idx ) );
                    ++in_idx;
                    ++out_idx;
            }
            else
                in_idx += keep[ i ].second;
        }

        if ( bam1_qual( in_bam )[ 0 ] == 0xFF )
            bam1_qual( out_bam )[ 0 ] = 0xFF;
        else {
            in_idx = 0, out_idx = 0;

            for ( unsigned i = 0; keep.size(); ++i ) {
                if ( keep[ i ].first )
                    for ( unsigned k = 0; k < keep[ i ].second; ++k )
                        bam1_qual( out_bam )[ out_idx++ ] = bam1_qual( in_bam )[ in_idx++ ];
                else
                    in_idx += keep[ i ].second;
            }
        }
    }

    if ( cols.size() )
        out_bam->core.pos = cols[ 0 ];

    out_bam->core.bin = bam_reg2bin( out_bam->core.pos, bam_calend( &out_bam->core, bam1_cigar( out_bam ) ) );

    if ( bam_cigar2qlen( &out_bam->core, bam1_cigar( out_bam ) ) != out_bam->core.l_qseq ) {
        cerr << "cig2qlen: " << bam_cigar2qlen( &out_bam->core, bam1_cigar( out_bam ) ) << endl;
        cerr << "l_qseq:   " << out_bam->core.l_qseq << endl;
        cerr << "invalid CIGAR string for sequence length" << endl;
        exit( 1 );
    }

    if ( !bam_validate1( NULL, out_bam ) ) {
        cerr << "record failed validation" << endl;
        exit( 1 );
    }

    return out_bam;
}


double prob_above_bg( const double lg_bg, const double lg_invbg, const int cov, const int k )
{
    double rv = exp( cov * lg_bg );

    for ( int i = 1; i < k; ++i )
        rv += exp( lg_choose( cov, i ) + i * lg_bg + ( cov - i ) * lg_invbg );

    return rv;
}


int main( int argc, const char * argv[] )
{
    args_t args = args_t( argc, argv );
    cov_citer cit;
    list< col_t > coverage;
    vector< col_t > variants;
    vector< pair< int, int > > data;

    // accumulate the data at each position in a linked list
    {
        bam1_t * in_bam = bam_init1();

        while ( args.bamin->next( in_bam ) ) {
            vector< triple_t > read;

            bam2vec( in_bam, read );

            include( coverage, read );
        }

        for ( cit = coverage.begin(); cit != coverage.end(); ++cit ) {
            obs_citer it = cit->obs.begin();
            int cov = 0, maj;

            if ( it == cit->obs.end() )
                continue;

            maj = it->second;
            cov += maj;

            for ( ++it; it != cit->obs.end(); ++it ) {
                if ( it->second > maj )
                    maj = it->second;
                cov += it->second;
            }

            data.push_back( make_pair( cov, maj ) );
        }

        bam_destroy1( in_bam );
    }

    // learn a joint multi-binomial model for the mutation rate classes
    {
        double lg_L, aicc, bg, lg_bg, lg_invbg;
        rateclass_t rc( data );
        vector< pair< double, double > > params;

        rc( lg_L, aicc, params );

        bg = params[ 0 ].second;
        lg_bg = log( bg );
        lg_invbg = log( 1.0 - bg );

        params_json_dump( stderr, lg_L, aicc, params );

        // determine which variants are above background and those which are not
        for ( cit = coverage.begin(); cit != coverage.end(); ++cit ) {
            col_t col = *cit;
            int cov = 0;

            for ( obs_citer it = col.obs.begin(); it != col.obs.end(); ++it )
                cov += it->second;

            for ( obs_iter it = col.obs.begin(); it != col.obs.end(); ++it ) {
                if ( prob_above_bg( lg_bg, lg_invbg, cov, it->second ) < args.cutoff )
                    it->second = 1;
                else
                    it->second = 0;
            }

            variants.push_back( col );
        }
    }

    // write out the input reads, but only with "real" variants this time
    {
        bam1_t * const in_bam = bam_init1();

        if ( !args.bamin->seek0() ) {
            cerr << "unable to seek( 0 )" << endl;
            exit( 1 );
        }

        if ( !args.bamout->write_header( args.bamin->hdr ) ) {
            cerr << "error writing out BAM header" << endl;
            exit( 1 );
        }

        while ( args.bamin->next( in_bam ) ) {
            vector< triple_t > read;

            bam2vec( in_bam, read );

            bam1_t * out_bam = punchout_read( in_bam, variants, read );

            if ( !args.bamout->write( out_bam ) ) {
                cerr << "error writing out read" << endl;
                exit( 1 );
            }

            bam_destroy1( out_bam );
        }

        bam_destroy1( in_bam );
    }

    return 0;
}
