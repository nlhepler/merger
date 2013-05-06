
#include <list>
#include <map>
#include <vector>

#include "aligned.h"
#include "coverage.hpp"
#include "util.hpp"

using std::list;
using std::map;
using std::vector;


typedef list< col_t >::iterator cov_iter;
typedef map< vector< char >, int >::iterator obs_iter;
typedef vector< triple_t >::const_iterator read_citer;


void include( list< col_t > & cov, const vector< triple_t > & read )
{
    cov_iter cit = cov.begin();
    read_citer rit = read.begin();

    for ( ; cit != cov.end() && rit != read.end(); ) {
        if ( cit->col == rit->col && cit->ins == rit->ins ) {
            obs_iter mit = cit->obs.find( rit->elem );

            if ( mit != cit->obs.end() )
                ++mit->second;
            else
                cit->obs[ rit->elem ] = 1;

            // increment iterators here
            ++rit;
            ++cit;
        }
        else if ( cit->col == rit->col && rit->ins ) {
            col_t newcol = { .col = rit->col, .ins = 1 };

            ++cit;

            // we cover this case above
            if ( cit->col == rit->col && cit->ins )
                continue;

            // it's a new insertion, so add it to our coverage
            newcol.obs[ rit->elem ] = 1;
            cov.insert( cit, newcol );

            // cit has already been incremented, but not rit
            ++rit;
        }
        else if ( cit->col == rit->col && cit->ins )
            ++cit;
        else if ( rit->col < cit->col ) {
            col_t newcol = { .col = rit->col, .ins = 0 };

            newcol.obs[ rit->elem ] = 1;
            cov.insert( cit, newcol );

            ++rit;
        }
        else // if ( cit->col < rit->col )
            ++cit;
    }
}
