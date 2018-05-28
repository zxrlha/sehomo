#ifndef _ATHL_TIMER_HPP
#define _ATHL_TIMER_HPP 11

#include <chrono>

///\brief a timer to measure time interval.
///
///There are three mode:
///
///1: start(), then call passed_time() to get time between current and start
///
///2: start(), then stop(), next(), stop(), next(), ..., stop(), call time() to get time between previous two stop or between start and the first stop
///
///3: start(), then stop(), stop(), ..., stop(), call time() to get time between start and previous stop.
///
///a timer can be used in many times, each time you should first call start(), and different times it can work on different modes.

class timer
{
public:
    ///default constructor
    timer();

    ///default copy constructor
    timer(const timer&) = default;

    ///default move constructor
    timer(timer&&) = default;

    ///\brief start the timer.
    ///
    ///the user must call this function
    void start();

    ///\brief stop the timer
    void stop();

    ///used in mode 2
    void next();

    ///get previous measured time
    double time();

    ///return time from start to now
    double passed_time();

    ///\brief return it is started or not.
    ///
    ///after the first time you call start(), it will be true, until you call clear
    bool is_started();

    ///clear this timer
    void clear();
protected:
    typedef std::chrono::steady_clock _t_clock;

    _t_clock::time_point _start;
    _t_clock::time_point _stop;
    double _time; ///<time cache
    bool _started;
};

#endif // _ATHL_TIMER_HPP
