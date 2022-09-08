#ifndef exceptions_hpp
#define exceptions_hpp

#include <stdexcept>

class BadMappingParameter: public std::runtime_error {
public:
    BadMappingParameter( const char* what_arg ) : std::runtime_error(what_arg) {};
};

#endif /* exceptions.hpp */
