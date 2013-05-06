
#include <vector>

#include "bam.h"

#ifndef BAMITER_H
#define BAMITER_H

typedef struct {
    int col;
    bool ins;
    std::vector< char > elem;
} triple_t;

void bam2vec( const bam1_t * const bam, std::vector< triple_t > & data );

char nuc2bits( const char nuc );

char bits2nuc( const char bits );

#endif // BAMITER_H
