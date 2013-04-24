
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

public:
    bam_header_t * hdr;

    bamfile_t( const char * path, bam_mode_t mode = READ );
    ~bamfile_t();
    bool next( aligned_t & aln );
    bool write_header( const bam_header_t * hdr_ = NULL );
    bool write( const char * const qname, aligned_t & aln );
};

#endif // BAMFILE_H
