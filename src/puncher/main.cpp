
#include <cmath>
#include <cstdlib>
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


bam1_t * punchout_read(
        const bam1_t * const read,
        const vector< col_t > & variants,
        const vector< triple_t > & read_
        )
{
    read_citer rit = read_.begin();
    var_citer vit = variants.begin();
    vector< pair< bool, unsigned > > keep;
    vector< int > cols;

    bam1_t * const read__ = bam_init1();

    if ( !read__ ) {
        cerr << "memory allocation error" << endl;
        exit( 1 );
    }

    read__->core = read->core;
    read__->core.l_qseq = 0;

    // build a vector of which positions we're keeping,
    // and how long our new cigar and seq will be
    for ( ; rit != read_.end() && vit != variants.end(); ) {
        if ( rit->col == vit->col && rit->ins == vit->ins ) {
            obs_citer it = vit->obs.find( rit->elem );

            if ( it == vit->obs.end() ) {
                cerr << "unknown variant observed, which is weird...( 1 )" << endl;
                exit( 1 );
            }

            if ( it->second ) {
                keep.push_back( make_pair( true, it->first.size() ) );
                cols.push_back( rit->col );
                read__->core.l_qseq += it->first.size();
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

    read__->data_len = read__->data_len;
    read__->m_data = read__->data_len;
    kroundup32( read__->m_data );
    read__->data = reinterpret_cast< uint8_t * >( malloc( read__->m_data * sizeof( uint8_t ) ) );

    if ( !read__->data ) {
        cerr << "memory allocation error" << endl;
        exit( 1 );
    }

    // copy the read name
    memcpy( read->data, read__->data, read->core.l_qname );

    // copy the cigar string
    {
        const uint32_t * cigar = bam1_cigar( read );
        int in_idx = 0, out_idx = 0;

        for ( int i = 0; i < read->core.n_cigar; ++i ) {
            const int op = cigar[ i ] & BAM_CIGAR_MASK;

            if ( op == BAM_CDEL ) {
                bam1_cigar( read__ )[ out_idx++ ] = cigar[ i ];
            }
            else if ( op == BAM_CINS || op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF ) {
                if ( keep[ in_idx++ ].first )
                    bam1_cigar( read__ )[ out_idx++ ] = cigar[ i ];
            }

            if ( out_idx >= read__->m_data ) {
                ++read__->m_data;
                kroundup32( read__->m_data );
                read__->data = reinterpret_cast< uint8_t * >( realloc( read__->data, read__->m_data * sizeof( uint8_t ) ) );
                if ( !read__->data ) {
                    cerr << "memory allocation error" << endl;
                    exit( 1 );
                }
            }
        }

        read__->core.n_cigar = out_idx;
    }

    // copy the sequence and quality scores
    {
        int in_idx = 0, out_idx = 0;

        for ( unsigned i = 0; i < keep.size(); ++i ) {
            if ( keep[ i ].first )
                for ( unsigned k = 0; k < keep[ i ].second; ++k ) {
                    bam1_seq_seti( bam1_seq( read__ ), out_idx, bam1_seqi( bam1_seq( read ), in_idx ) );
                    ++in_idx;
                    ++out_idx;
            }
            else
                in_idx += keep[ i ].second;
        }

        if ( bam1_qual( read )[ 0 ] == 0xFF )
            bam1_qual( read__ )[ 0 ] = 0xFF;
        else {
            in_idx = 0, out_idx = 0;
            for ( unsigned i = 0; keep.size(); ++i ) {
                if ( keep[ i ].first )
                    for ( unsigned k = 0; k < keep[ i ].second; ++k )
                        bam1_qual( read__ )[ out_idx++ ] = bam1_qual( read )[ in_idx++ ];
                else
                    in_idx += keep[ i ].second;
            }
        }
    }

    if ( cols.size() )
        read__->core.pos = cols[ 0 ];

    memcpy( bam1_aux( read ), bam1_aux( read__ ), read->l_aux );

    read__->core.bin = bam_reg2bin( read__->core.pos, bam_calend( &read__->core, bam1_cigar( read__ ) ) );

    if ( bam_cigar2qlen( &read__->core, bam1_cigar( read__ ) ) != read__->core.l_qseq ) {
        cerr << "invalid CIGAR string for sequence length" << endl;
        exit( 1 );
    }

    if ( !bam_validate1( NULL, read__ ) ) {
        cerr << "record failed validation" << endl;
        exit( 1 );
    }

    return read__;
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
        bam1_t * read = bam_init1();

        while ( args.bamin->next( read ) ) {
            vector< triple_t > read_;

            bam2vec( read, read_ );
            include( coverage, read_ );
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

        bam_destroy1( read );
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
        bam1_t * read = bam_init1();

        args.bamin->seek( 0 );

        if ( !args.bamout->write_header( args.bamin->hdr ) ) {
            cerr << "error writing out BAM header" << endl;
            exit( 1 );
        }

        while ( args.bamin->next( read ) ) {
            vector< triple_t > read_;

            bam2vec( read, read_ );

            bam1_t * read__ = punchout_read( read, variants, read_ );

            if ( !args.bamout->write( read__ ) ) {
                cerr << "error writing out read" << endl;
                exit( 1 );
            }

            bam_destroy1( read__ );
        }

        bam_destroy1( read );
    }

    return 0;
}
