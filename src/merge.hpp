
#include <vector>

#include "aligned.hpp"
#include "bamfile.hpp"


#define MERGE_SIZE 128


namespace merge
{
    class nuc_t
    {
    public:
        int col;
        int cov;
        aligned::op_t op;
        char nuc;
        char qual;

        nuc_t(
            const int col,
            const aligned::op_t op,
            const char nuc,
            const char qual,
            const int cov = 1
            );
    };

    class cluster_t : public std::vector< nuc_t >
    {
    public:
        int ncontrib;

        cluster_t();
        cluster_t( const aligned::aligned_t & seq );
        
        int lpos() const;
        int rpos() const;

        aligned::aligned_t to_aligned() const;
        cluster_t merge(
            const cluster_t & other,
            const int min_overlap,
            const bool tol_ambigs,
            const bool tol_gaps
            ) const;
    };

    bool ncontrib_cmp(
        const cluster_t & x,
        const cluster_t & y
        );

    void merge_clusters(
        const unsigned nread,
        const int min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps,
        std::vector< cluster_t > & clusters
        );

    std::vector< aligned::aligned_t > merge_reads(
        bamfile::bamfile_t & bamfile,
        const int min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps,
        const bool discard
        );
}
