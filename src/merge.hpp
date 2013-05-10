
#include <vector>

#include "aligned.hpp"
#include "bamfile.hpp"


#define MERGE_SIZE 128


namespace merge
{
    bool ncontrib_cmp(
        const aligned::aligned_t & x,
        const aligned::aligned_t & y
        );

    bool merge_two(
        const aligned::aligned_t & x,
        const aligned::aligned_t & y,
        const unsigned min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps,
        aligned::aligned_t & merged
        );

    void merge_clusters(
        const unsigned nread,
        const unsigned min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps,
        std::vector< aligned::aligned_t > & clusters
        );

    int merge_reads(
        bamfile::bamfile_t & bamfile,
        const unsigned min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps,
        const bool discard,
        std::vector< aligned::aligned_t > & clusters
        );
}
