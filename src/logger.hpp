#ifndef LOGGER_HPP
#define LOGGER_HPP

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
    LOG_DEBUG = 1,
    LOG_INFO = 2,
    LOG_WARNING = 3,
    LOG_ERROR = 4,
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

    void set_level(int level) { this->level = level; }
    LogStream& debug() { return _debug; }
    LogStream& info() { return _info; }
    LogStream& warning() { return _warning; }
    LogStream& error() { return _error; }

private:
    Logger()
        : level(0)
        , _os(std::cerr)
        , _debug(LogStream(LOG_DEBUG, *this))
        , _info(LogStream(LOG_INFO, *this))
        , _warning(LogStream(LOG_WARNING, *this))
        , _error(LogStream(LOG_ERROR, *this))
    { }
    int level;
    std::ostream& _os;
    LogStream _debug;
    LogStream _info;
    LogStream _warning;
    LogStream _error;

    friend class LogStream;
};


template <typename T>
LogStream& LogStream::operator<<(const T& val) {
    if (level >= logger.level) {
        logger._os << val;
    }
    return *this;
}

// This overload is required for supporting std::endl
inline LogStream& LogStream::operator<<(std::ostream& (*f)(std::ostream&)) {
    if (level >= logger.level) {
        f(logger._os);
    }
    return *this;
}

#endif
