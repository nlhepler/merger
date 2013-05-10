
#include <list>
#include <map>
#include <vector>

#include "aligned.hpp"
#include "coverage.hpp"


using std::list;
using std::map;
using std::string;
using std::vector;

using aligned::INS;
using aligned::aligned_t;


namespace coverage
{
    void elem_t::get_seq( string & str ) const
    {
        vector< char >::const_iterator it;

        str.clear();

        for ( it = begin(); it != end(); ++it )
            str.push_back( *it );
    }

    void coverage_t::include( const aligned_t & read )
    {
        iterator cit = begin();
        aligned_t::const_iterator rit = read.begin();
        elem_t elem;

        for ( ; cit != end() && rit != read.end(); ) {
            if ( cit->col == rit->col && cit->op == rit->op ) {
                rit->get_seq( elem );
                map< elem_t, int >::iterator mit = cit->obs.find( elem );

                if ( mit != cit->obs.end() )
                    ++mit->second;
                else
                    cit->obs[ elem ] = 1;

                // increment iterators here
                ++rit;
                ++cit;
            }
            else if ( cit->col == rit->col && rit->op == INS ) {
                cov_t cov = { .col = rit->col, .op = INS };

                ++cit;

                // we cover this case above
                if ( cit->col == rit->col && cit->op == INS )
                    continue;

                // it's a new insertion, so add it to our coverage
                rit->get_seq( elem );
                cov.obs[ elem ] = 1;
                insert( cit, cov );

                // cit has already been incremented, but not rit
                ++rit;
            }
            else if ( rit->col < cit->col ) {
                cov_t cov = { .col = rit->col, .op = rit->op };

                rit->get_seq( elem );
                cov.obs[ elem ] = 1;
                insert( cit, cov );

                ++rit;
            }
            // we cover these cases above
            else if ( cit->col == rit->col && cit->op == INS )
                ++cit;
            else // if ( cit->col < rit->col )
                ++cit;
        }

        for ( ; rit != read.end(); ++rit ) {
            cov_t cov = { .col = rit->col, .op = rit->op };

            rit->get_seq( elem );
            cov.obs[ elem ] = 1;
            insert( cit, cov );
        }
    }
}
