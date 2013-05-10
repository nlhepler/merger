
#include <algorithm>
#include <vector>

#include "bam.h"

#include "aligned.hpp"
#include "args.hpp"
#include "bamfile.hpp"
#include "merge.hpp"

using std::sort;
using std::vector;

using aligned::aligned_t;

using merge::merge_reads;


int main( int argc, const char * argv[] )
{
    args_t args = args_t( argc, argv );
    vector< aligned_t > clusters;
    vector< aligned_t >::iterator cluster;
    bam1_t * const bam = bam_init1();
    unsigned i, j;

    if ( !bam )
        goto error;

    if ( merge_reads(
                *args.bamin,
                args.min_overlap,
                args.tol_ambigs,
                args.tol_gaps,
                bool( args.bamdiscard ),
                clusters
                ) < 0 )
        goto error;

    args.bamout->write_header( args.bamin->hdr );

    if ( args.bamdiscard )
        args.bamdiscard->write_header( args.bamin->hdr );

    for ( cluster = clusters.begin(), i = 0, j = 0; cluster != clusters.end(); ++cluster ) {
        if ( cluster->ncontrib >= args.min_reads ) {
            char name[ 256 ];
            snprintf( name, 256, "cluster%u_%dr", i++, cluster->ncontrib );
            cluster->name += name;
            
            if ( !cluster->to_bam( bam ) )
                goto error;
            
            args.bamout->write( bam );
        }
        else if ( args.bamdiscard ) {
            char name[ 256 ];
            snprintf( name, 256, "cluster%u_%dr", j++, cluster->ncontrib );
            cluster->name += name;

            if ( !cluster->to_bam( bam ) )
                goto error;
            
            args.bamdiscard->write( bam );
        }
    }

    bam_destroy1( bam );

    return 0;

error:
    bam_destroy1( bam );

    return -1;
}
