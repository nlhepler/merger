
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <list>
#include <utility>
#include <vector>

#include "bam.h"

#include "aligned.hpp"
#include "args.hpp"
#include "bamfile.hpp"
#include "coverage.hpp"
#include "math.hpp"
#include "rateclass.hpp"
#include "util.hpp"


using std::cerr;
using std::cout;
using std::endl;
using std::exp;
using std::list;
using std::log;
using std::make_pair;
using std::map;
using std::pair;
using std::vector;

using aligned::INS;
using aligned::aligned_t;
using coverage::cov_t;
using coverage::coverage_t;
using coverage::elem_t;
using math::prob_background;
using rateclass::params_json_dump;
using rateclass::rateclass_t;
using util::bits2nuc;


typedef list< cov_t >::const_iterator cov_citer;
typedef list< cov_t >::iterator cov_iter;
typedef map< elem_t, int >::const_iterator obs_citer;
typedef map< elem_t, int >::iterator obs_iter;
typedef vector< cov_t >::const_iterator var_citer;


inline
uint32_t cigval( const int op, const int nop )
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
        const vector< cov_t > & variants,
        const aligned_t & read
        )
{
    aligned_t::const_iterator rit = read.begin();
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

#if 0
        cerr << "r: ( " << rit->col << ", " << rit->op == INS << ", ";
        for ( unsigned j = 0; j < rit->elem.size(); ++j )
            cerr << bits2nuc( rit->elem[ j ] );
        cerr << " )" << endl;
        cerr << "v: ( " << vit->col << ", " << vit->op == INS << " )" << endl;
#endif

        if ( rit->col == vit->col && rit->op == vit->op ) {
            elem_t elem;
            obs_citer it;

            rit->get_seq( elem );
            it = vit->obs.find( elem );

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
        else if ( rit->col == vit->col && rit->op == INS ) {
            ++vit;
        }
        // matching column in read, but insertion in variants:
        // we got ahead in read somehow, which should never happen
        else if ( rit->col == vit->col && vit->op == INS ) {
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

    // cerr << "seqlen: " << out_bam->core.l_qseq << endl;

    out_bam->data_len = in_bam->data_len;
    out_bam->m_data = out_bam->data_len;
    out_bam->data = NULL;

    realloc_data( out_bam );

    // copy the read name
    memcpy( bam1_qname( out_bam ), bam1_qname( in_bam ), in_bam->core.l_qname );

    // update the cigar string
    {
        const uint32_t * cigar = bam1_cigar( in_bam );
        int in_idx = 0, out_idx = 0;
        vector< int > cigar_;

        for ( int i = 0; i < in_bam->core.n_cigar; ++i ) {
            const int op = cigar[ i ] & BAM_CIGAR_MASK;
            const int nop = cigar[ i ] >> BAM_CIGAR_SHIFT;

            if ( op == BAM_CDEL )
                for ( int j = 0; j < nop; ++j )
                    cigar_.push_back( BAM_CDEL );
            else if ( op == BAM_CINS ) {
                if ( keep[ in_idx ].second != unsigned( nop ) ) {
                    cerr << "assumptions violated ( 1 )" << endl;
                    exit( 1 );
                }
                if ( keep[ in_idx ].first )
                    for ( int j = 0; j < nop; ++j )
                        cigar_.push_back( op );
                else
                    for ( int j = 0; j < nop; ++j )
                        cigar_.push_back( BAM_CDEL );
                ++in_idx;
            }
            else if ( op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF ) {
                for ( int j = 0; j < nop; ++j, ++in_idx ) {
                    if ( keep[ in_idx ].second != 1 ) {
                        cerr << "assumptions violated ( 2 )" << endl;
                        exit( 1 );
                    }
                    if ( keep[ in_idx ].first )
                        cigar_.push_back( op );
                    else
                        cigar_.push_back( BAM_CDEL );
                }
            }
        }

        if ( cigar_.size() ) {
            int op = cigar_[ 0 ];
            int nop = 1;
    
            for ( unsigned i = 1; i < cigar_.size(); ++i ) {
                if ( cigar_[ i ] == op )
                    ++nop;
                else {
                    bam1_cigar( out_bam )[ out_idx++ ] = cigval( op, nop );
                    op = cigar_[ i ];
                    nop = 1;
                }
            
                if ( out_bam->core.l_qname + out_idx >= out_bam->m_data ) {
                    ++out_bam->m_data;
                    realloc_data( out_bam );
                }
            }

            bam1_cigar( out_bam )[ out_idx++ ] = cigval( op, nop );
        }

        out_bam->core.n_cigar = out_idx;
    }

    out_bam->l_aux = in_bam->l_aux;

    out_bam->data_len = (
            out_bam->core.l_qname +
            4 * out_bam->core.n_cigar +
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


int main( int argc, const char * argv[] )
{
    args_t args = args_t( argc, argv );
    coverage_t coverage;
    vector< cov_t > variants;
    vector< pair< int, int > > data;

    // accumulate the data at each position in a linked list
    {
        cov_citer cit;
        bam1_t * in_bam = bam_init1();

        while ( args.bamin->next( in_bam ) ) {
            aligned_t read( in_bam );
            coverage.include( read );
        }

        for ( cit = coverage.begin(); cit != coverage.end(); ++cit ) {
            int cov = 0;

            for ( obs_citer it = cit->obs.begin(); it != cit->obs.end(); ++it )
                cov += it->second;
           
            for ( obs_citer it = cit->obs.begin(); it != cit->obs.end(); ++it )
                if ( it->second )
                    data.push_back( make_pair( cov, it->second ) );

#if 0
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
#endif
        }

        bam_destroy1( in_bam );
    }

    // learn a joint multi-binomial model for the mutation rate classes
    {
        cov_iter cit;
        double lg_L, aicc, bg, lg_bg, lg_invbg;
        rateclass_t rc( data );
        vector< pair< double, double > > params;

        rc( lg_L, aicc, params );

        bg = params[ 0 ].second;
        lg_bg = log( bg );
        lg_invbg = log( 1.0 - bg );

        params_json_dump( stderr, lg_L, aicc, params );

        // cerr << "background: " << bg << endl;

        // determine which variants are above background and those which are not
        for ( cit = coverage.begin(); cit != coverage.end(); ++cit ) {
            if ( cit->op == INS )
                continue;

            int cov = 0;

            for ( obs_citer it = cit->obs.begin(); it != cit->obs.end(); ++it )
                cov += it->second;

            for ( obs_iter it = cit->obs.begin(); it != cit->obs.end(); ++it ) {
                const double p = prob_background( lg_bg, lg_invbg, cov, it->second );
                if ( p < args.cutoff ) {
                    cout << cit->col << "\t" << cov << "\t" << it->second;
                    for ( unsigned i = 0; i < it->first.size(); ++i )
                        cout << bits2nuc( it->first[ i ] );
                    cout << ":" << p << endl;
                    it->second = 1;
                }
                else {
                    it->second = 0;
                }
            }

#if 0
            variants.push_back( *cit );
#endif
        }
    }

    return 0;

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
            aligned_t read( in_bam );

            bam1_t * const out_bam = punchout_read( in_bam, variants, read );

            if ( !out_bam->core.l_qseq )
                continue;

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
