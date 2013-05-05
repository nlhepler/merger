
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

    fprintf( stdout, "{\n  \"lg_L\":     % .3f\n  \"aicc\":     % .3f\n  \"rates\":   [ ", lg_L, aicc );

    for ( it = params.begin(); it != params.end(); ++it ) {
        if ( it == params.begin() )
            fprintf( stdout, "%.7f", it->second );
        else
            fprintf( stdout, ", %.7f", it->second );
    }

    cout << " ]," << endl;
    cout << "  \"weights\": [ ";

    for ( it = params.begin(); it != params.end(); ++it ) {
        if ( it == params.begin() )
            fprintf( stdout, "%.7f", it->first );
        else
            fprintf( stdout, ", %.7f", it->first );
    }

    cout << " ]" << endl << "}" << endl;

    return 0;
}
