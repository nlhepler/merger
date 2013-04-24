
#include <vector>

#include "argparse.hpp"
#include "aligned.h"

using std::vector;

typedef struct {
    int min_overlap;
    int min_reads;
    int tol_gaps;
    int tol_ambigs;
} opts_t;

char nuc2bits( const char nuc );

char bits2nuc( const char bits );

int merge(
    args_t & args,
    vector< aligned_t > & clusters
    );

void aligned_destroy( aligned_t & read );
