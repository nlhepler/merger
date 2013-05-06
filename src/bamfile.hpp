
#include <vector>

#include "bam.h"

#include "aligned.h"

#ifndef BAMFILE_H
#define BAMFILE_H

enum bam_mode_t { READ, WRITE };

class bamfile_t
{
private:
    bamFile fp;
    bam1_t * buf;
    bam_index_t * idx;
    long zero;

public:
    bam_header_t * hdr;

    bamfile_t( const char * path, bam_mode_t mode = READ );
    ~bamfile_t();
    bool next( aligned_t & aln );
    bool next( bam1_t * const aln );
    void fetch( std::vector< aligned_t > & reads, int begin, int end, int tid = 0 );
    bool seek0();
    bool write_header( const bam_header_t * hdr_ = NULL );
    bool write( const char * const qname, aligned_t & aln, int begin = 0, int end = 0 );
    bool write( const bam1_t * const aln );
};

#endif // BAMFILE_H
