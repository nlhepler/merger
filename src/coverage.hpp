
#include <list>
#include <map>
#include <vector>

#include "util.hpp"

#ifndef COVERAGE_H
#define COVERAGE_H

typedef struct {
    int col;
    int ins;
    std::map< std::vector< char >, int > obs;
} col_t;

void include( std::list< col_t > & cov, const std::vector< triple_t > & read );

#endif // COVERAGE_H
