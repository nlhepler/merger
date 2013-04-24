
#include <cstdlib>
#include <algorithm>
#include <iostream>

#include "argparse.hpp"
#include "bamfile.hpp"
#include "merge.hpp"

using std::sort;
using std::cerr;
using std::endl;


#define MERGE_SIZE 128

#define MAX( a, b ) (((a) > (b)) ? (a) : (b))
#define MIN( a, b ) (((a) < (b)) ? (a) : (b))

enum cmp_t { LT, GT, EQ };
enum res_t { SUCCESS, FAILURE, ERROR };


char nuc2bits( const char nuc )
{
    switch ( nuc ) {
    case 'A': return 1;
    case 'C': return 2;
    case 'G': return 4;
    case 'T': return 8;
    case 'M': return 1 | 2;
    case 'R': return 1 | 4;
    case 'W': return 1 | 8;
    case 'S': return 2 | 4;
    case 'Y': return 2 | 8;
    case 'K': return 4 | 8;
    case 'V': return 1 | 2 | 4;
    case 'H': return 1 | 2 | 8;
    case 'D': return 1 | 4 | 8;
    case 'B': return 2 | 4 | 8;
    default:  return 1 | 2 | 4 | 8;
    }
}

char bits2nuc( const char bits )
{
    switch ( bits ) {
    case 1:           return 'A';
    case 2:           return 'C';
    case 4:           return 'G';
    case 8:           return 'T';
    case (1 | 2):     return 'M';
    case (1 | 4):     return 'R';
    case (1 | 8):     return 'W';
    case (2 | 4):     return 'S';
    case (2 | 8):     return 'Y';
    case (4 | 8):     return 'K';
    case (1 | 2 | 4): return 'V';
    case (1 | 2 | 8): return 'H';
    case (1 | 4 | 8): return 'D';
    case (2 | 4 | 8): return 'B';
    default:          return 'N';
    }
}

cmp_t pos_cmp( const pos_t & x, const pos_t & y )
{
    // test the column, then the insertion, then equal
    if ( x.col > y.col )
        return GT;
    else if ( y.col > x.col )
        return LT;
    else if ( x.ins > y.ins )
        return GT;
    else if ( y.ins < x.ins )
        return LT;
    else
        return EQ;
}


bool ncontrib_cmp( const aligned_t & x, const aligned_t & y )
{
    return x.ncontrib > y.ncontrib;
}


bool aln_cmp( const aligned_t & x, const aligned_t & y )
{
    if ( x.lpos < y.lpos )
        return true;
    else if ( y.lpos < x.lpos )
        return false;
    else if ( x.len >= y.len )
        return true;
    return false;
}

void cerr_triple( const pos_t & x, bool end=true )
{
    cerr << x.col << " " << x.ins << " " << bits2nuc( x.nuc );
    if ( end )
        cerr << endl;
    else
        cerr << ", ";
}


#ifdef DEBUG
#define ABORT( msg ) { \
    cerr << msg << ": "; \
    cerr_triple( xs.data[ xidx ], false ); \
    cerr_triple( ys.data[ yidx ] ); \
    return FAILURE; \
}
#else
#define ABORT( msg ) { return FAILURE; }
#endif


res_t merge_two(
    const aligned_t & xs,
    const aligned_t & ys,
    const args_t & args,
    aligned_t & merged
    )
{
    int overlap = 0;
    int xidx = 0;
    int yidx = 0;
    int midx = 0;
    int mlen = 0;
    int cmp = 0;

    if ( !xs.len || !ys.len )
        ABORT( "insufficient length" )

    // if there is absolutely no hope to reach min_overlap, skip
    if ( xs.rpos < ys.lpos + args.min_overlap &&
            ys.rpos < xs.lpos + args.min_overlap )
        ABORT( "no opportunity for sufficient overlap" )

    cmp = pos_cmp( xs.data[ xidx ], ys.data[ yidx ] );

    // cerr_triple( xs.data[ xidx ], false );
    // cerr_triple( xs.data[ yidx ] );

    // disregard overhangs
    if ( cmp == LT ) {
        for ( ; cmp == LT && xidx + 1 < xs.len; ) {
            ++xidx;
            cmp = pos_cmp( xs.data[ xidx ], ys.data[ yidx ] );
        }
        // if its not a match, it's a gap, then backup one in xs
        if ( cmp == GT && !args.tol_gaps )
            ABORT( "no gaps in ys" )
        // mlen is now the # of xs already visited
        mlen += xidx;
    }
    else if ( cmp == GT ) {
        for ( ; cmp == GT && yidx + 1 < ys.len; ) {
            ++yidx;
            cmp = pos_cmp( xs.data[ xidx ], ys.data[ yidx ] );
        }
        // if its not a match, it's a gap, then backup one in ys
        if ( cmp == GT && !args.tol_gaps )
            ABORT( "no gaps in xs" )
        // mlen is now the # of ys already visited
        mlen += yidx;
    }

    // compute the amount of overlap
    for ( ; xidx < xs.len && yidx < ys.len; ++mlen ) {
        cmp = pos_cmp( xs.data[ xidx ], ys.data[ yidx ] );
        if ( cmp == LT ) {
            if ( !args.tol_gaps )
                ABORT( "no gaps in xs" )
            ++xidx;
        }
        else if ( cmp == GT ) {
            if ( !args.tol_gaps )
                ABORT( "no gaps in ys" )
            ++yidx;
        }
        // if the nucleotides match, move ahead
        else if ( xs.data[ xidx ].nuc == ys.data[ yidx ].nuc ||
                    ( args.tol_ambigs && xs.data[ xidx ].nuc & ys.data[ yidx ].nuc ) ) {
                ++overlap;
                ++xidx;
                ++yidx;
        }
        // nucleotides do not match, abort early
        else
            ABORT( "mismatch" )
    }

    if ( overlap < args.min_overlap )
        ABORT( "insufficient overlap" )

    // get the remainder of either xs or ys, whichever remains
    if ( xidx < xs.len )
        mlen += xs.len - xidx;
    else if ( yidx < ys.len )
        mlen += ys.len - yidx;

    merged.len = mlen;
    merged.data = reinterpret_cast< pos_t * >( malloc( merged.len * sizeof( pos_t ) ) );
    merged.lpos = MIN( xs.lpos, ys.lpos );
    merged.rpos = MAX( xs.rpos, ys.rpos );
    merged.ncontrib = xs.ncontrib + ys.ncontrib;

    if ( !merged.data )
        ABORT( "memory allocation error" )

    for ( xidx = 0, yidx = 0; xidx < xs.len && yidx < ys.len; ) {
        cmp = pos_cmp( xs.data[ xidx ], ys.data[ yidx ] );
        if ( cmp == LT )
            merged.data[ midx++ ] = xs.data[ xidx++ ];
        else if ( cmp == GT )
            merged.data[ midx++ ] = ys.data[ yidx++ ];
        else {
            merged.data[ midx ] = xs.data[ xidx ];
            merged.data[ midx ].cov += ys.data[ yidx ].cov;
            merged.data[ midx ].nuc = MIN( xs.data[ xidx ].nuc, ys.data[ yidx ].nuc );
            ++xidx;
            ++yidx;
            ++midx;
        }
    }

    if ( xidx < xs.len )
        for ( ; xidx < xs.len; ++xidx )
            merged.data[ midx++ ] = xs.data[ xidx ];
    else if ( yidx < ys.len )
        for ( ; yidx < ys.len; ++yidx )
            merged.data[ midx++ ] = ys.data[ yidx ];

#ifndef DEBUG
    if ( midx < mlen )
        cerr << "error: failed to fill 'merged' data" << endl;
    else if ( midx > mlen )
        cerr << "error: overfilled 'merged' data" << endl;
#endif

    return SUCCESS;
}

int merge_clusters(
    const args_t & args,
    vector< aligned_t > & clusters
    )
{
    vector< aligned_t >::iterator a, b;
    int nclusters;

begin:
    sort( clusters.begin(), clusters.end(), ncontrib_cmp );

    for ( nclusters = 0, a = clusters.begin(); a != clusters.end(); ++a ) {
        for ( b = a + 1; b != clusters.end(); ++b ) {
            aligned_t merged = { .data = NULL };
            res_t res = merge_two( *a, *b, args, merged );
            if ( res == SUCCESS ) {
                aligned_destroy( *a );
                aligned_destroy( *b );
                // replace i and remove j
                *a = merged;
                clusters.erase( b );
                goto begin;
            }
            else if ( res == ERROR )
                return -1;
        }
        if ( a->ncontrib >= args.min_reads )
            ++nclusters;
    }

    return nclusters;
}

int merge(
    args_t & args,
    vector< aligned_t > & clusters
    )
{
    vector< aligned_t >::iterator cluster;
    size_t merge_size = MERGE_SIZE;
    int nclusters, nread = 0;
    aligned_t read = { .data = NULL };

    while ( args.bamin->next( read ) ) {
        if ( nread % 100 == 0 )
            fprintf( stderr, "\rprocessed: %9i reads (%6lu clusters)", nread, clusters.size() );
        for ( cluster = clusters.begin(); cluster != clusters.end(); ++cluster ) {
            aligned_t merged = { .data = NULL };
            res_t res = merge_two( read, *cluster, args, merged );
            if ( res == SUCCESS ) {
                // destroy our read and cluster
                aligned_destroy( read );
                aligned_destroy( *cluster );
                *cluster = merged;
                if ( clusters.size() > merge_size ) {
                    if ( merge_clusters( args, clusters ) < 0 )
                        goto error;
                    merge_size *= 2;
                }
                goto next;
            }
            else if ( res == ERROR )
                goto error;
        }
        clusters.push_back( read );
next:
        ++nread;
    }

    nclusters = merge_clusters( args, clusters );

    fprintf( stderr, "\rprocessed: %9d reads (%6lu clusters)", nread, clusters.size() );

    if ( nclusters < 0 )
        goto error;

    sort( clusters.begin(), clusters.end(), aln_cmp );

    return nclusters;

error:
    aligned_destroy( read );

    for ( cluster = clusters.begin(); cluster != clusters.end(); ++cluster )
        aligned_destroy( *cluster );

    return -1;
}

void aligned_destroy( aligned_t & read )
{
    free( read.data );
    read.data = NULL;
    read.len = 0;
}
