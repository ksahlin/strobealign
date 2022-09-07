#ifndef exceptions_hpp
#define exceptions_hpp

#include <stdexcept>

class mapping_parameter_exception: public std::runtime_error {
public:
    mapping_parameter_exception( const char* what_arg ) : std::runtime_error(what_arg) {};
};

#endif /* exceptions.hpp */
