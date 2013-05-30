
#include <list>
#include <map>
#include <vector>

#include "aligned.hpp"
#include "coverage.hpp"
#include "util.hpp"


using std::list;
using std::map;
using std::string;
using std::vector;

using aligned::INS;
using aligned::aligned_t;
using aligned::op_t;
using util::bits2nuc;


namespace coverage
{
    void elem_t::get_seq( string & str ) const
    {
        elem_t::const_iterator it;

        str.clear();

        for ( it = begin(); it != end(); ++it )
            str.push_back( bits2nuc( *it ) );
    }

    cov_t::cov_t( const int col, const op_t op ) :
        col( col ),
        op( op )
    {
    }

    void coverage_t::include( const aligned_t & read )
    {
        iterator cit = begin();
        aligned_t::const_iterator rit = read.begin();
        for ( ; cit != end() && rit != read.end(); ) {
            elem_t elem;

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
                cov_t cov( rit->col, INS );
                
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
                cov_t cov( rit->col, rit->op );

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
            elem_t elem;

            cov_t cov( rit->col, rit->op );

            rit->get_seq( elem );
            cov.obs[ elem ] = 1;
            insert( cit, cov );
        }
    }
}
