
#include <vector>

#include "bam.h"

#include "aligned.hpp"

#ifndef BAMFILE_H
#define BAMFILE_H

namespace bamfile
{
    enum bam_mode_t { READ, WRITE };

    class bamfile_t
    {
    private:
        bamFile fp;
        bam_index_t * idx;
        long zero;

    public:
        bam_header_t * hdr;

        bamfile_t( const char * path, bam_mode_t mode = READ, bool index = false );
        ~bamfile_t();
        bool next( bam1_t * const bam );
        void fetch(
                std::vector< aligned::aligned_t > & reads,
                const int begin,
                const int end,
                const int tid = 0
                );
        bool seek0();
        bool write_header( const bam_header_t * hdr_ = NULL );
        bool write( const bam1_t * const aln );
    };
}

#endif // BAMFILE_H
