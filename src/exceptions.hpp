//
// exceptions.hpp
//
// Created by Paul-Theodor Pyl on 2022-09-07
//


#ifndef exceptions_hpp
#define exceptions_hpp

#include <exception>

class k_exception: public std::exception {
    virtual const char* what() const throw();
};

class s_lt_k_exception: public std::exception {
    virtual const char* what() const throw();
};

class k_minus_s_even_exception: public std::exception {
    virtual const char* what() const throw();
};

class max_seed_length_exception: public std::exception {
    virtual const char* what() const throw();
};

#endif /* exceptions.hpp */
