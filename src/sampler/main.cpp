
#include <algorithm>
#include <vector>

#include "aligned.h"
#include "args.hpp"
#include "bamfile.hpp"
#include "merge.hpp"


using std::sort;
using std::vector;

using merge::ERROR;
using merge::SUCCESS;
using merge::aligned_destroy;
using merge::merge_clusters;
using merge::merge_two;
using merge::ncontrib_cmp;
using merge::res_t;


#define MERGE_SIZE 128
#define MIN( a, b ) ( ( (a) < (b) ) ? (a) : (b) )


int merge_(
    args_t & args,
    vector< aligned_t > & reads,
    vector< aligned_t > & clusters
    )
{
    vector< aligned_t >::iterator read, cluster;
    size_t merge_size = MERGE_SIZE, nread = 0;

    for ( read = reads.begin(); read != reads.end(); ++read ) {
        
        if ( nread % 100 == 0 )
            fprintf( stderr, "\rprocessed: %9lu reads (%6lu clusters)", nread, clusters.size() );

        for ( cluster = clusters.begin(); cluster != clusters.end(); ++cluster ) {
            aligned_t merged = { .data = NULL };
            res_t res = merge_two( *read, *cluster, args, merged );
            if ( res == SUCCESS ) {
                // destroy our read and cluster
                aligned_destroy( *read );
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
        clusters.push_back( *read );
next:
        ++nread;
    }

    if ( merge_clusters( nread, args, clusters ) < 0 )
        goto error;

    fprintf( stderr, "\rprocessed: %9lu reads (%6lu clusters)\n", nread, clusters.size() );

    sort( clusters.begin(), clusters.end(), ncontrib_cmp );

    return clusters.size();

error:
    for ( ; read != reads.end(); ++read )
        aligned_destroy( *read );

    return 0;
}

// main ------------------------------------------------------------------------------------------------------------- //

int main( int argc, const char * argv[] )
{
    args_t args = args_t( argc, argv );
    size_t k = 0;

    args.bamout->write_header( args.bamin->hdr );

    for ( int i = args.begin; i < args.end; i += args.stride ) {
        int j = MIN( i + args.window_size, args.end );
        vector< aligned_t > reads;
        vector< aligned_t > clusters;
        vector< aligned_t >::iterator cluster;

        args.bamin->fetch( reads, i, j );

        if ( merge_( args, reads, clusters ) )
            for ( cluster = clusters.begin();
                    cluster->ncontrib >= args.min_reads && cluster != clusters.end();
                    ++cluster ) {
                char qname[ 256 ];
                snprintf( qname, 256, "cluster%lu_%ib_%ie_%dr", k++, i, j, cluster->ncontrib );
                args.bamout->write( qname, *cluster, i, j );
            }

        for ( cluster = clusters.begin(); cluster != clusters.end(); ++cluster )
            aligned_destroy( *cluster );
    }

    return 0;
}
