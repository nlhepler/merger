
#include <algorithm>
#include <cstdio>
#include <utility>

#include "args.hpp"
#include "bamfile.hpp"
#include "merge.hpp"
#include "util.hpp"

using std::make_pair;
using std::pair;
using std::sort;
using std::vector;

using aligned::INS;
using aligned::aligned_t;
using aligned::op_t;
using aligned::pos_t;
using bamfile::bamfile_t;
using util::bits2nuc;


#ifdef MAX
#undef MAX
#endif

#define MAX( a, b ) ( ( (a) > (b) ) ? (a) : (b) )

#ifdef MIN
#undef MIN
#endif

#define MIN( a, b ) ( ( (a) < (b) ) ? (a) : (b) )


namespace merge
{
    nuc_t::nuc_t( const int col, const op_t op, const char nuc, const char qual, const int cov ) :
        col( col ),
        cov( cov ),
        op( op ),
        nuc( nuc ),
        qual( qual )
    {
    }


    cluster_t::cluster_t() :
        ncontrib( 0 )
    {
    }


    cluster_t::cluster_t( const aligned_t & seq ) :
        ncontrib( seq.ncontrib ? seq.ncontrib : 1 )
    {
        aligned_t::const_iterator it;
        vector< pair< char, char > >::const_iterator jt;

        for ( it = seq.begin(); it != seq.end(); ++it )
            for ( jt = it->data.begin(); jt != it->data.end(); ++jt )
                push_back( nuc_t( it->col, it->op, jt->first, jt->second ) );
    }


    int cluster_t::lpos() const
    {
        cluster_t::const_iterator it = begin();
        return ( it == end() ) ? -1: it->col;
    }


    int cluster_t::rpos() const
    {
        if ( lpos() < 0 )
            return -1;
        return lpos() + size();
    }


    inline
    int mean( vector< int > & data )
    {
        vector< int >::const_iterator it;
        int sum = 0;

        if ( !data.size() )
            return 0;

        if ( data.size() == 1 )
            return data[ 0 ];

        for ( it = data.begin(); it != data.end(); ++it )
            sum += *it;

        return int( ( 0.5 + sum ) / data.size() );
    }


    aligned_t cluster_t::to_aligned() const
    {
        aligned_t cluster;
        cluster_t::const_iterator it = begin();
        int col = it->col;
        op_t op = it->op;
        vector< int > cov( 1, it->cov );
        vector< pair< char, char > > data( 1, make_pair( it->nuc, it->qual ) );

        for ( ++it; it != end(); ++it ) {
            if ( it->col == col && it->op == op ) {
                cov.push_back( it->cov );
                data.push_back( make_pair( it->nuc, it->qual ) );
            }
            else {
                cluster.push_back( pos_t( col, op, data, mean( cov ) ) );
                col = it->col;
                op = it->op;
                cov.clear();
                cov.push_back( it->cov );
                data.clear();
                data.push_back( make_pair( it->nuc, it->qual ) );
            }
        }

        cluster.push_back( pos_t( col, op, data, mean( cov ) ) );

        cluster.ncontrib = ncontrib;

        return cluster;
    }


    bool ncontrib_cmp( const cluster_t & x, const cluster_t & y )
    {
        return x.ncontrib > y.ncontrib;
    }


    cluster_t cluster_t::merge(
        const cluster_t & other,
        const int min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps
        ) const
    {
        cluster_t::const_iterator i = begin(), j = other.begin();
        cluster_t m;
        int overlap = 0;

        if ( i == end() || j == other.end() )
            return m;

        if ( rpos() < other.lpos() + min_overlap && other.rpos() < lpos() + min_overlap )
            return m;

        if ( i->col < j->col )
            for ( ; i->col < j->col && i != end(); m.push_back( *( i++ ) ) );
        else if ( i->col > j->col )
            for ( ; j->col < i->col && j != other.end(); m.push_back( *( j++ ) ) );

        while ( i != end() && j != other.end() ) {
            if ( i->col < j->col ) { // && i->col != INS ) {
                if ( !tol_gaps )
                    goto abort;
                m.push_back( *( i++ ) );
            }
            else if ( j->col < i->col ) { // && j->col != INS ) {
                if ( !tol_gaps )
                    goto abort;
                m.push_back( *( j++ ) );
            }
            else if ( i->op == j->op ) {
                if ( i->nuc == j->nuc || ( tol_ambigs && i->nuc & j->nuc ) ) {
                    if ( i->nuc == j->nuc )
                        ++overlap;
                    m.push_back(
                        nuc_t(
                            i->col, i->op,
                            MIN( i->nuc, j->nuc ), MAX( i->qual, j->qual ),
                            i->cov + j->cov
                            )
                        );
                    ++i;
                    ++j;
                }
                else
                    goto abort;
            }
            else if ( i->op == INS ) {
                if ( !tol_gaps )
                    goto abort;
                m.push_back( *( j++ ) );
            }
            else { // if ( j->op == INS ) {
                if ( !tol_gaps )
                    goto abort;
                m.push_back( *( i++ ) );
            }
#if 0
            else
                goto abort;
#endif
        }

        if ( overlap < min_overlap )
            goto abort;

        if ( i != end() )
            for ( ; i != end(); m.push_back( *( i++ ) ) );
        else if ( j != other.end() )
            for ( ; j != other.end(); m.push_back( *( j++ ) ) );

        m.ncontrib = ncontrib + other.ncontrib;

        return m;

abort:
        m.clear();

        return m;
    }


    void merge_clusters(
        const unsigned nread,
        const int min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps,
        vector< cluster_t > & clusters
        )
    {
        bool repeat;

        do {
            repeat = false;

            fprintf( stderr, "\rprocessed: %9u reads (%6lu clusters)", nread, clusters.size() );
            fflush( stderr );

            sort( clusters.begin(), clusters.end(), ncontrib_cmp );

            for ( int i = 0; i < int( clusters.size() ); ++i ) {
                bool stop = false;
                #pragma omp parallel for
                for ( int j = i + 1; j < int( clusters.size() ); ++j ) {
                    if ( stop )
                        continue;
                    
                    cluster_t merged = clusters[ i ].merge( clusters[ j ], min_overlap, tol_ambigs, tol_gaps );
                    if ( merged.size() ) {
                        #pragma omp critical
                        if ( !stop ) {
                            // replace i and remove j
                            clusters[ i ] = merged;
                            clusters.erase( clusters.begin() + j );
                            repeat = stop = true;
                            #pragma omp flush( stop )
                        }
                    }
                }
            }
        } while ( repeat );
    }


    bool aln_cmp( const cluster_t & x, const cluster_t & y )
    {
        if ( x.lpos() < y.lpos() )
            return true;
        else if ( y.lpos() < x.lpos() )
            return false;
        else if ( x.size() >= y.size() )
            return true;
        return false;
    }


    vector< aligned_t > merge_reads(
        bamfile_t & bamfile,
        const int min_overlap,
        const bool tol_ambigs,
        const bool tol_gaps,
        const bool discard
        )
    {
        vector< cluster_t >::iterator cluster;
        vector< cluster_t > clusters;
        vector< aligned_t > discards, rv;
        bam1_t * const bam = bam_init1();
        unsigned merge_size = MERGE_SIZE, nread = 1;

        if ( !bam )
            goto error;

        for ( ; bamfile.next( bam ); ++nread ) {
            aligned_t orig( bam );
            cluster_t read( orig );
            bool stop = false;

            /*
            if ( nread % 1000 == 0 )
                sort( clusters.begin(), clusters.end(), ncontrib_cmp );
            */

            if ( read.size() < unsigned( min_overlap ) ) {
                if ( !discard )
                    discards.push_back( orig );
                continue;
            }

            #pragma omp parallel for
            for ( int i = 0; i < int( clusters.size() ); ++i ) {
                if ( stop )
                    continue;

                cluster_t merged = clusters[ i ].merge( read, min_overlap, tol_ambigs, tol_gaps );

                if ( merged.size() ) {
                    #pragma omp critical
                    if ( !stop ) {
                        clusters[ i ] = merged;
                        stop = true;
                        #pragma omp flush( stop )
                    }
                }
            }

            if ( !stop )
                clusters.push_back( read );

            if ( clusters.size() >= merge_size ) {
                merge_clusters( nread, min_overlap, tol_ambigs, tol_gaps, clusters );
                merge_size *= 2;
            }

            if ( nread % 100 == 0 ) {
                fprintf( stderr, "\rprocessed: %9u reads (%6lu clusters)", nread, clusters.size() );
                fflush( stderr );
            }
        }

        merge_clusters( nread, min_overlap, tol_ambigs, tol_gaps, clusters );

        fprintf( stderr, "\rprocessed: %9u reads (%6lu clusters)\n", nread, clusters.size() );
        fflush( stderr );

        sort( clusters.begin(), clusters.end(), aln_cmp );

        rv.reserve( clusters.size() + discards.size() );

        for ( cluster = clusters.begin(); cluster != clusters.end(); ++cluster )
            rv.push_back( cluster->to_aligned() );

        if ( !discard )
            clusters.insert( clusters.end(), discards.begin(), discards.end() );

        bam_destroy1( bam );

        return rv;

    error:
        rv.clear();

        bam_destroy1( bam );

        return rv;
    }
}
