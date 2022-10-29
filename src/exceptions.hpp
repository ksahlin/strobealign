#ifndef exceptions_hpp
#define exceptions_hpp

#include <stdexcept>
#include <string>

class BadParameter: public std::runtime_error {
public:
    BadParameter( const char* what_arg ) : std::runtime_error(what_arg) {};
};

class InvalidFasta : public std::runtime_error {
public:
    InvalidFasta(std::string message) : runtime_error(message) { }
};

class InvalidIndexFile : public std::runtime_error {
public:
    InvalidIndexFile(std::string message) : runtime_error(message) { }
};

class InvalidFile : public std::runtime_error {
public:
    InvalidFile(std::string message) : runtime_error(message) { }
};

#endif
