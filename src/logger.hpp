#ifndef STROBEALIGN_LOGGER_HPP
#define STROBEALIGN_LOGGER_HPP

#include <ostream>
#include <iostream>


/*
Very simple logging.

Usage:

    logger = Logger::get();  // returns the logging singleton
    logger.set_level(LOG_INFO);
    logger.info() << "info message" << std::endl; // printed
    logger.debug() << "debug message" << std::endl; // not printed

*/


enum LOG_LEVELS {
    LOG_TRACE = 1,
    LOG_DEBUG = 2,
    LOG_INFO = 3,
    LOG_WARNING = 4,
    LOG_ERROR = 5,
};

class Logger;

class LogStream {
public:
    LogStream(int level, Logger& logger): level(level), logger(logger) { }
    template <typename T>
    LogStream& operator<<(const T& val);
    LogStream& operator<<(std::ostream& (*f)(std::ostream&));
private:
    int level;
    Logger& logger;
};


class Logger {
public:
    static Logger& get() {
        static Logger instance;
        return instance;
    }
    Logger(Logger const&) = delete;
    void operator=(Logger const&) = delete;

    void set_level(int level) { this->_level = level; }
    int level() const { return this->_level; }
    LogStream& trace() { return _trace; }
    LogStream& debug() { return _debug; }
    LogStream& info() { return _info; }
    LogStream& warning() { return _warning; }
    LogStream& error() { return _error; }

private:
    Logger()
        : _level(LOG_WARNING)
        , _os(std::cerr)
        , _trace(LogStream(LOG_TRACE, *this))
        , _debug(LogStream(LOG_DEBUG, *this))
        , _info(LogStream(LOG_INFO, *this))
        , _warning(LogStream(LOG_WARNING, *this))
        , _error(LogStream(LOG_ERROR, *this))
    { }
    int _level;
    std::ostream& _os;
    LogStream _trace;
    LogStream _debug;
    LogStream _info;
    LogStream _warning;
    LogStream _error;

    friend class LogStream;
};


template <typename T>
LogStream& LogStream::operator<<(const T& val) {
    if (level >= logger._level) {
        logger._os << val;
    }
    return *this;
}

// This overload is required for supporting std::endl
inline LogStream& LogStream::operator<<(std::ostream& (*f)(std::ostream&)) {
    if (level >= logger._level) {
        f(logger._os);
    }
    return *this;
}

#endif
