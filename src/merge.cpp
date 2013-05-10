
#include <algorithm>
#include <cstdio>
#include <iostream>
#include <utility>

#include "args.hpp"
#include "bamfile.hpp"
#include "merge.hpp"
#include "util.hpp"

using std::cerr;
using std::endl;
using std::make_pair;
using std::pair;
using std::sort;
using std::vector;

using aligned::INS;
using aligned::aligned_t;
using aligned::pos_t;
using bamfile::bamfile_t;
using util::bits2nuc;


#ifdef MAX
#undef MAX
#endif

#define MAX( a, b ) ( ( (a) > (b) ) ? (a) : (b) )

#ifdef MIN
#undef MIN
#endif

#define MIN( a, b ) ( ( (a) < (b) ) ? (a) : (b) )


namespace merge
{
    bool ncontrib_cmp( const aligned_t & x, const aligned_t & y )
    {
        return x.ncontrib > y.ncontrib;
    }


    bool merge_two(
        const aligned_t & x,
        const aligned_t & y,
        const unsigned min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps,
        aligned_t & merged
        )
    {
        aligned_t::const_iterator a = x.begin(), b = y.begin();
        unsigned overlap = 0;

        merged.clear();

        if ( a == x.end() || b == y.end() )
            return false;

        if ( a->col < b->col )
            for ( ; a->col < b->col && a != x.end(); merged.push_back( *( a++ ) ) );
        else if ( a->col > b->col )
            for ( ; b->col < a->col && b != x.end(); merged.push_back( *( b++ ) ) );

        while ( a != x.end() || b != y.end() ) {
            if ( a->col < b->col ) {
                if ( !tol_gaps )
                    goto abort;
                merged.push_back( *( a++ ) );
            }
            else if ( b->col < a->col ) {
                if ( !tol_gaps )
                    goto abort;
                merged.push_back( *( b++ ) );
            }
            else if ( a->op == b->op ) {
                vector< pair< char, char > >::const_iterator ap = a->data.begin();
                vector< pair< char, char > >::const_iterator bp = b->data.begin();
                vector< pair< char, char > > data;

                for ( ; ap != a->data.end() && bp != b->data.end(); ++ap, ++bp ) {
                    if ( ap->first == bp->first || ( tol_ambigs && ap->first & bp->first ) ) {
                        if ( ap->first == bp->first )
                            ++overlap;

                        data.push_back(
                            make_pair(
                                MIN( ap->first, bp->first ),
                                MAX( ap->second, bp->second )
                                )
                            );
                    }
                    else
                        goto abort;
                }

                merged.push_back(
                    pos_t(
                        a->col,
                        a->op,
                        data,
                        a->cov + b->cov
                        )
                    );

                ++a;
                ++b;
            }
            else if ( a->op == INS ) {
                if ( !tol_gaps )
                    goto abort;
                merged.push_back( *( b++ ) );
            }
            else { // b->op == INS
                if ( !tol_gaps )
                    goto abort;
                merged.push_back( *( a++ ) );
            }
        }

        if ( overlap < min_overlap )
            goto abort;

        if ( a != x.end() )
            for ( ; a != x.end(); merged.push_back( *( a++ ) ) );
        else if ( b != y.end() )
            for ( ; b != y.end(); merged.push_back( *( b++ ) ) );

        if ( !x.ncontrib )
            ++merged.ncontrib;
        else
            merged.ncontrib += x.ncontrib;

        if ( !y.ncontrib )
            ++merged.ncontrib;
        else
            merged.ncontrib += y.ncontrib;

        return true;

abort:
        merged.clear();
        
        return false;
    }


    void merge_clusters(
        const unsigned nread,
        const unsigned min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps,
        vector< aligned_t > & clusters
        )
    {
        vector< aligned_t >::iterator a, b;
        bool repeat;

        do {
            repeat = false;

            sort( clusters.begin(), clusters.end(), ncontrib_cmp );

            for ( a = clusters.begin(); a != clusters.end(); ++a ) {
begin:
                for ( b = a + 1; b != clusters.end(); ++b ) {
                    aligned_t merged;
                    if ( merge_two( *a, *b, min_overlap, tol_ambigs, tol_gaps, merged ) ) {
                        // replace i and remove j
                        *a = merged;
                        clusters.erase( b );
                        fprintf( stderr, "\rprocessed: %9u reads (%6lu clusters)", nread, clusters.size() );
                        repeat = true;
                        goto begin;
                    }
                }
            }
        } while ( repeat );
    }


    bool aln_cmp( const aligned_t & x, const aligned_t & y )
    {
        if ( x.lpos() < y.lpos() )
            return true;
        else if ( y.lpos() < x.lpos() )
            return false;
        else if ( x.size() >= y.size() )
            return true;
        return false;
    }


    int merge_reads(
        bamfile_t & bamfile,
        const unsigned min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps,
        const bool discard,
        vector< aligned_t > & clusters
        )
    {
        vector< aligned_t > discarded;
        vector< aligned_t >::iterator cluster;
        unsigned merge_size = MERGE_SIZE, nread = 0;
        bam1_t * const bam = bam_init1();

        if ( !bam )
            goto error;

        for ( ; bamfile.next( bam ); ++nread ) {
            aligned_t read( bam );

            // if ( nread % 100 == 0 )
                fprintf( stderr, "\rprocessed: %9u reads (%6lu clusters)", nread, clusters.size() );
                fflush( stderr );
            /*
            if ( nread % 1000 == 0 )
                sort( clusters.begin(), clusters.end(), ncontrib_cmp );
            */

            if ( read.size() < min_overlap ) {
                if ( !discard )
                    discarded.push_back( read );
                continue;
            }

            for ( cluster = clusters.begin(); cluster != clusters.end(); ++cluster ) {
                aligned_t merged;
                if ( merge_two( *cluster, read, min_overlap, tol_ambigs, tol_gaps, merged ) ) {
                    *cluster = merged;
                    if ( clusters.size() > merge_size ) {
                        merge_clusters( nread, min_overlap, tol_ambigs, tol_gaps, clusters );
                        merge_size *= 2;
                    }
                    goto next;
                }
            }
            clusters.push_back( read );
next:;
        }

        merge_clusters( nread, min_overlap, tol_ambigs, tol_gaps, clusters );

        fprintf( stderr, "\rprocessed: %9u reads (%6lu clusters)\n", nread, clusters.size() );

        sort( clusters.begin(), clusters.end(), aln_cmp );

        if ( !discard ) {
            clusters.reserve( clusters.size() + discarded.size() );
            clusters.insert( clusters.end(), discarded.begin(), discarded.end() );
        }

        bam_destroy1( bam );

        return clusters.size();

    error:
        bam_destroy1( bam );

        return -1;
    }
}
