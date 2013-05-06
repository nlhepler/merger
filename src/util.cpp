
#include <vector>

#include "bam.h"

#include "util.hpp"


using std::vector;

void bam2vec( const bam1_t * const bam, vector< triple_t > & data )
{
    const uint32_t * cigar = bam1_cigar( bam );
    int col = bam->core.pos, idx = 0;
    
    for ( int i = 0; i < bam->core.n_cigar; ++i ) {
        const int op = cigar[ i ] & BAM_CIGAR_MASK;
        const int nop = cigar[ i ] >> BAM_CIGAR_SHIFT;

        if ( op == BAM_CDEL ) {
            col += nop;
            continue;
        }
        else if ( op == BAM_CINS ) {
            triple_t tri;
            
            for ( int j = 0; j < nop; ++j ) {
                tri.elem.push_back( bam1_seqi( bam1_seq( bam ), idx ) );
                ++idx;
            }

            tri.col = col;
            tri.ins = true;

            data.push_back( tri );

            col += nop;
        }
        else if ( op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF ) {
            for ( int j = 0; j < nop; ++j ) {
                triple_t tri;

                tri.elem.push_back( bam1_seqi( bam1_seq( bam ), idx ) );
                ++idx;

                tri.col = col;
                tri.ins = false;

                data.push_back( tri );

                col += 1;
            }
        }
    }
}
