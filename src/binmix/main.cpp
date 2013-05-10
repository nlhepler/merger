
#include <cstdio>
#include <iostream>
#include <utility>
#include <vector>

#include "rateclass.hpp"
#include "math.hpp"


using std::cin;
using std::cout;
using std::endl;
using std::make_pair;
using std::pair;
using std::vector;

using math::weighted_harmonic_mean;
using rateclass::params_json_dump;
using rateclass::rateclass_t;


int main( int argc, char * argv[] )
{
    double aicc, lg_L;
    vector< pair< double, double > > params;
    vector< pair< int, int > > data;

    while ( cin.good() ) {
        int c, m;
        cin >> c >> m;
        data.push_back( make_pair( c, m ) );
    }

    rateclass_t rc( data );

    rc( lg_L, aicc, params );

    params_json_dump( stderr, lg_L, aicc, params, weighted_harmonic_mean( params ) );

    return 0;
}
