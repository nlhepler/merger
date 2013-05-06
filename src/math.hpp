
#include <cmath>

#ifndef MATH_H
#define MATH_H

inline
double lg_choose( const int n, const int k )
{
    double rv = 0.0;

    for ( int i = 1; i <= k; ++i )
        rv += std::log( n - k + i ) - std::log( i );

    return rv;
}

#endif // MATH_H
