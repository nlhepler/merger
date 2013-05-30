
#include <vector>
#include <string>
#include <utility>

#include "bam.h"


#ifndef ALIGNED_H
#define ALIGNED_H

namespace aligned
{
    enum op_t {
        DEL = BAM_CDEL,
        INS = BAM_CINS,
        MATCH = BAM_CMATCH,
        EQUAL = BAM_CEQUAL,
        DIFF = BAM_CDIFF
    };

    class pos_t : public std::vector< std::pair< char, char > >
    {
    public:
        int col;
        int cov;
        op_t op;

        pos_t(
            const int col,
            const op_t op,
            const int cov = 1
            );

        pos_t(
            const int col,
            const op_t op,
            std::vector< std::pair< char, char > > data,
            const int cov = 1
            );

        void get_qual( char * qual ) const;
        void get_seq( char * seq ) const;
        void get_seq( std::string & str ) const;
        void get_seq( std::vector< char > & vec ) const;
    };


    class aligned_t : public std::vector< pos_t >
    {
    private:
        int tid;
        int qual;
        int flag;
        int mtid;
        int mpos;
        int isize;

    public:
        std::string name;
        int ncontrib;
        // TODO, implement auxiliary data

        aligned_t();
        aligned_t( const bam1_t * const bam );

        int lpos() const;
        int rpos() const;

        bool to_bam( bam1_t * const bam ) const;
        std::vector< pos_t > to_vector() const;
    };
}

#endif // ALIGNED_H
