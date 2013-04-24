
#include <vector>

#include "aligned.h"
#include "argparse.hpp"
#include "bamfile.hpp"
#include "merge.hpp"

using std::vector;

// main ------------------------------------------------------------------------------------------------------------- //

int main( int argc, const char * argv[] )
{
    args_t args = args_t( argc, argv );
    vector< aligned_t > clusters;
    int i;

    merge( args, clusters );

    args.bamout->write_header( args.bamin->hdr );

    for ( i = 0; i < clusters.size(); ++i ) {
        if ( clusters[ i ].ncontrib >= args.min_reads ) {
            char qname[ 256 ];
            snprintf( qname, 256, "cluster%d_%dr", i, clusters[ i ].ncontrib );
            args.bamout->write( qname, clusters[ i ] );
        }
    }

    return 0;
}
