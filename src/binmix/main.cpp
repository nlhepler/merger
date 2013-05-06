
#include <cstdio>
#include <iostream>
#include <utility>
#include <vector>

#include "rateclass_em.hpp"

using std::cin;
using std::cout;
using std::endl;
using std::make_pair;
using std::pair;
using std::vector;

int main( int argc, char * argv[] )
{
    double aicc, lg_L;
    vector< pair< double, double > > params;
    vector< pair< double, double > >::iterator it;
    vector< pair< int, int > > data;

    while ( cin.good() ) {
        int c, m;
        cin >> c >> m;
        data.push_back( make_pair( c, m ) );
    }

    rateclass_t rc( data );

    rc( lg_L, aicc, params );

    params_json_dump( stderr, lg_L, aicc, params );

    return 0;
}
