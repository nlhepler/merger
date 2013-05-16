
#include <algorithm>
#include <iostream>
#include <vector>

#include "bam.h"

#include "aligned.hpp"
#include "args.hpp"
#include "bamfile.hpp"
#include "merge.hpp"

using std::cerr;
using std::endl;
using std::sort;
using std::vector;

using aligned::aligned_t;

using merge::merge_reads;


int main( int argc, const char * argv[] )
{
    args_t args = args_t( argc, argv );
    bam1_t * const bam = bam_init1();
    unsigned i, j;
    
    vector< aligned_t >::iterator cluster;
    vector< aligned_t > clusters = merge_reads(
        *args.bamin,
        args.min_overlap,
        args.tol_ambigs,
        args.tol_gaps,
        bool( args.bamdiscard )
        );

    if ( !bam ) {
        cerr << "memory allocation error" << endl; 
        goto error;
    }
    else if ( !clusters.size() ) {
        cerr << "no clusters found" << endl;
        goto error;
    }

    args.bamout->write_header( args.bamin->hdr );

    if ( args.bamdiscard )
        args.bamdiscard->write_header( args.bamin->hdr );

    for ( cluster = clusters.begin(), i = 0, j = 0; cluster != clusters.end(); ++cluster ) {
        if ( cluster->ncontrib >= args.min_reads ) {
            char name[ 256 ];
            snprintf( name, 256, "cluster%u_%dr", i++, cluster->ncontrib );
            cluster->name += name;
            
            if ( !cluster->to_bam( bam ) ) {
                cerr << "error converting to BAM format" << endl;
                goto error;
            }
            
            if ( !args.bamout->write( bam ) ) {
                cerr << "error writing to BAM_OUT" << endl;
                goto error;
            }
        }
        else if ( args.bamdiscard ) {
            char name[ 256 ];
            snprintf( name, 256, "cluster%u_%dr", j++, cluster->ncontrib );
            cluster->name += name;

            if ( !cluster->to_bam( bam ) ) {
                cerr << "error converting to BAM format" << endl;
                goto error;
            }
            
            if ( !args.bamdiscard->write( bam ) ) {
                cerr << "error writing to BAM_DISCARD" << endl;
                goto error;
            }
        }
    }

    bam_destroy1( bam );

    return 0;

error:
    bam_destroy1( bam );

    return -1;
}
