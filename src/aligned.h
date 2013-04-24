
#ifndef ALIGNED_H
#define ALIGNED_H

typedef struct {
    int col;
    int ins;
    int cov;
    char nuc;
} pos_t;

typedef struct {
    pos_t * data;
    int len;
    int lpos;
    int rpos;
    int ncontrib;
} aligned_t;

#endif // ALIGNED_H
