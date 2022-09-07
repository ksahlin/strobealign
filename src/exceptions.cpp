#include "exceptions.hpp"

const char* k_exception::what() const throw()
{
    return "k not in [8,32]";
}

const char* s_lt_k_exception::what() const throw()
{
    return "s is larger than k";
}

const char* k_minus_s_even_exception::what() const throw()
{
    return "k - s should be an even number to create canonical syncmers. Please set s to e.g. k-2, k-4, k-6, ...";
}

const char*  max_seed_length_exception::what() const throw()
{
    return "maximum seed length (-m <max_dist>) is larger than 255";
}
