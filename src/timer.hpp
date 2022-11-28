#ifndef TIMER_HPP
#define TIMER_HPP

#include <chrono>

// A timer that automatically starts on construction
class Timer {
public:
    Timer() : start_time(std::chrono::high_resolution_clock::now()) {
    }

    // Return time elapsed since construction
    std::chrono::duration<double> duration() const {
        return std::chrono::high_resolution_clock::now() - start_time;
    }

    std::chrono::duration<double>::rep elapsed() const {
        return duration().count();
    }

private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time;
};

#endif
