
#include <algorithm>
#include <vector>

#include "aligned.h"
#include "args.hpp"
#include "bamfile.hpp"
#include "merge.hpp"

using std::sort;
using std::vector;


#define MERGE_SIZE 128


bool aln_cmp( const aligned_t & x, const aligned_t & y )
{
    if ( x.lpos < y.lpos )
        return true;
    else if ( y.lpos < x.lpos )
        return false;
    else if ( x.len >= y.len )
        return true;
    return false;
}

// merge ------------------------------------------------------------------------------------------------------------ //

int merge(
    args_t & args,
    vector< aligned_t > & clusters
    )
{
    vector< aligned_t > discard;
    vector< aligned_t >::iterator cluster;
    size_t merge_size = MERGE_SIZE;
    int nclusters, nread = 0;
    aligned_t read = { .data = NULL };

    while ( args.bamin->next( read ) ) {

        if ( nread % 100 == 0 )
            fprintf( stderr, "\rprocessed: %9i reads (%6lu clusters)", nread, clusters.size() );

        /*
        if ( nread % 1000 == 0 )
            sort( clusters.begin(), clusters.end(), ncontrib_cmp );
        */

        if ( read.len < args.min_overlap ) {
            if ( args.bamdiscard )
                discard.push_back( read );
            goto next;
        }

        for ( cluster = clusters.begin(); cluster != clusters.end(); ++cluster ) {
            aligned_t merged = { .data = NULL };
            res_t res = merge_two( read, *cluster, args, merged );
            if ( res == SUCCESS ) {
                // destroy our read and cluster
                aligned_destroy( read );
                aligned_destroy( *cluster );
                *cluster = merged;
                if ( clusters.size() > merge_size ) {
                    if ( merge_clusters( nread, args, clusters ) < 0 )
                        goto error;
                    merge_size *= 2;
                }
                goto next;
            }
            else if ( res == ERROR ) {
                aligned_destroy( merged );
                goto error;
            }
        }
        clusters.push_back( read );
next:
        ++nread;
    }

    nclusters = merge_clusters( nread, args, clusters );

    fprintf( stderr, "\rprocessed: %9d reads (%6lu clusters)\n", nread, clusters.size() );

    if ( nclusters < 0 )
        goto error;

    if ( args.bamdiscard ) {
        clusters.reserve( clusters.size() + discard.size() );
        clusters.insert( clusters.end(), discard.begin(), discard.end() );
    }

    sort( clusters.begin(), clusters.end(), aln_cmp );

    return nclusters;

error:
    aligned_destroy( read );

    for ( cluster = clusters.begin(); cluster != clusters.end(); ++cluster )
        aligned_destroy( *cluster );

    return -1;
}

// main ------------------------------------------------------------------------------------------------------------- //

int main( int argc, const char * argv[] )
{
    args_t args = args_t( argc, argv );
    vector< aligned_t > clusters;
    vector< aligned_t >::iterator cluster;
    size_t i, j;

    merge( args, clusters );

    args.bamout->write_header( args.bamin->hdr );

    if ( args.bamdiscard )
        args.bamdiscard->write_header( args.bamin->hdr );

    for ( cluster = clusters.begin(), i = 0, j = 0; cluster != clusters.end(); ++cluster ) {
        if ( cluster->ncontrib >= args.min_reads ) {
            char qname[ 256 ];
            snprintf( qname, 256, "cluster%lu_%dr", i++, cluster->ncontrib );
            args.bamout->write( qname, *cluster );
        }
        else if ( args.bamdiscard ) {
            char qname[ 256 ];
            snprintf( qname, 256, "cluster%lu_%dr", j++, cluster->ncontrib );
            args.bamdiscard->write( qname, *cluster );
        }
    }

    return 0;
}
