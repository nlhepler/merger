
#include <utility>
#include <vector>

#ifndef RATECLASS_EM_H
#define RATECLASS_EM_H

void params_json_dump(
        FILE * const file,
        const double lg_L,
        const double aicc,
        const std::vector< std::pair< double, double > > & params
        );

class rateclass_t
{
private:
    const std::vector< std::pair< int, int > > & data;

public:
    rateclass_t( const std::vector< std::pair< int, int > > & data );
    void operator()(
        double & lg_L,
        double & aicc,
        std::vector< std::pair< double, double > > & params,
        const int nrestart = 50
        ) const;
};

#endif
