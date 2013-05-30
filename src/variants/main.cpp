
#include <cstdio>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>
#include <utility>

#include "bam.h"

#include "aligned.hpp"
#include "args.hpp"
#include "bamfile.hpp"
#include "coverage.hpp"
#include "math.hpp"
#include "rateclass.hpp"
#include "util.hpp"


using std::cout;
using std::endl;
using std::exp;
using std::list;
using std::log;
using std::make_pair;
using std::map;
using std::pair;
using std::string;
using std::vector;

using aligned::MATCH;
using aligned::aligned_t;
using coverage::coverage_t;
using coverage::elem_t;
using math::prob_background;
using math::weighted_harmonic_mean;
using rateclass::params_json_dump;
using rateclass::rateclass_t;
using util::bits2nuc;


int main( int argc, const char * argv[] )
{
    args_t args = args_t( argc, argv );
    coverage_t::const_iterator cit;
    coverage_t coverage;
    vector< pair< int, int > > data;
    bam1_t * const in_bam = bam_init1();

    while ( args.bamin->next( in_bam ) ) {
        aligned_t read( in_bam );
        coverage.include( read );
    }

    bam_destroy1( in_bam );
    
    for ( cit = coverage.begin(); cit != coverage.end(); ++cit ) {
        map< elem_t, int >::const_iterator it;

        if ( cit->op != MATCH )
            continue;

        it = cit->obs.begin();

        if ( it == cit->obs.end() )
            continue;

        int cov = it->second;
        int max = it->second;

        for ( ++it; it != cit->obs.end(); ++it ) {
            cov += it->second;
            if ( it->second > max )
                max = it->second;
        }

        for ( it = cit->obs.begin(); it != cit->obs.end(); ++it )
            if ( it->second && it->second != max )
                data.push_back( make_pair( cov, cov - it->second ) );
    }

    rateclass_t rc( data, 3 );
    double lg_L, aicc;
    vector< pair< double, double > > params;

    rc( lg_L, aicc, params );

    const double bg = weighted_harmonic_mean( params );
    const double lg_bg = log( bg );
    const double lg_invbg = log( 1.0 - bg );

    params_json_dump( stderr, lg_L, aicc, params, bg );

    for ( cit = coverage.begin(); cit != coverage.end(); ++cit ) {
        map< elem_t, int >::const_iterator it;

        if ( cit->op != MATCH )
            continue;

        it = cit->obs.begin();

        if ( it == cit->obs.end() )
            continue;

        int cov = it->second;
        int max = it->second;

        for ( ++it; it != cit->obs.end(); ++it ) {
            cov += it->second;
            if ( it->second > max )
                max = it->second;
        }

        string css;

        for ( it = cit->obs.begin(); it != cit->obs.end(); ++it )
            if ( it->second == max ) {
                string elem;
                it->first.get_seq( elem );
                css.append( elem );
                css.push_back( '/' );
            }

        // erase the trailing slash, in a compatible way
        css.erase( --css.end() );

        for ( it = cit->obs.begin(); it != cit->obs.end(); ++it ) {
            string elem;

            if ( !it->second || it->second == max )
                continue;

            const double prob = prob_background( lg_bg, lg_invbg, cov, it->second );

            if ( prob >= args.cutoff )
                continue;

            fprintf( stdout, "%d\t%s\t%d\t", cit->col + 1, css.c_str(), cov );

            it->first.get_seq( elem );
            fprintf( stdout, "%s:%d:%.3e\n", elem.c_str(), it->second, prob );
        }

        fflush( stdout );
    }

    return 0;
}
