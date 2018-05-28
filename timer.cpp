#include "timer.hpp"

timer::timer()
{
    clear();
}

void timer::start()
{
    _started = true;
    _start = _t_clock::now();
}

void timer::stop()
{
    _stop = _t_clock::now();
    _time = double(std::chrono::duration_cast<std::chrono::nanoseconds>(_stop - _start).count()) / 1e9;
}

void timer::next()
{
    _start = _stop;
}

double timer::time()
{
    return _time;
}

double timer::passed_time()
{
    return double(std::chrono::duration_cast<std::chrono::nanoseconds>(_t_clock::now() - _start).count()) / 1e9;
}

bool timer::is_started()
{
    return _started;
}

void timer::clear()
{
    _started = false;
    _time = 0;
}
