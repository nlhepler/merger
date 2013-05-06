
#include <vector>

#include "args.hpp"
#include "aligned.h"

enum res_t { SUCCESS, FAILURE, ERROR };

bool ncontrib_cmp( const aligned_t & x, const aligned_t & y );

res_t merge_two(
    const aligned_t & xs,
    const aligned_t & ys,
    const args_t & args,
    aligned_t & merged
    );

int merge_clusters(
    const int nread,
    const args_t & args,
    std::vector< aligned_t > & clusters
    );

void aligned_destroy( aligned_t & read );
