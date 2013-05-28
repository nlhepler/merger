
#include <list>
#include <map>
#include <vector>

#include "aligned.hpp"


#ifndef COVERAGE_H
#define COVERAGE_H

namespace coverage
{
    class elem_t : public std::vector< char >
    {
    public:
        void get_seq( std::string & str ) const;
    };

    class cov_t
    {
    public:
        int col;
        aligned::op_t op;
        std::map< elem_t, int > obs;

        cov_t( const int col, const aligned::op_t op );
    };

    class coverage_t : public std::list< cov_t >
    {
    public:
        void include( const aligned::aligned_t & read );
    };
}

#endif // COVERAGE_H
