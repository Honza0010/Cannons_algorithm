#ifndef TimerRTH
#define TimerRTH

class TimerRT
{
public:
    TimerRT();
    void reset();
    void stop();
    void start();
    double getTime();

protected:
    double stop_time;
    double total_time;
    bool stop_state;

private:
#ifdef _WIN32
    double frequency;
    __int64 start_count;
    __int64 stop_count;
#endif
};

extern TimerRT defaultTimer;

#ifdef _WIN32
#include <Windows.h>
#else
#include <sys/time.h>
#endif

TimerRT defaultTimer;

TimerRT::TimerRT() {
#ifdef _WIN32
    LARGE_INTEGER freq;
    QueryPerformanceFrequency(&freq);
    frequency = double(freq.QuadPart);
#endif
    reset();
}

void TimerRT::reset() {
#ifdef _WIN32
    LARGE_INTEGER count;
    QueryPerformanceCounter(&count);
    start_count = count.QuadPart;
    total_time = 0.0;
    stop_state = false;
#else
    struct timeval tp;
    gettimeofday(&tp, NULL);
    stop_time = (double)tp.tv_sec + 1.0e-6 * (double)tp.tv_usec;
    total_time = 0.0;
    stop_state = false;
#endif
}

void TimerRT::stop() {
    if (!stop_state) {
#ifdef _WIN32
        LARGE_INTEGER count;
        QueryPerformanceCounter(&count);
        stop_count = count.QuadPart;
        total_time += (stop_count - start_count) / frequency;
#else
        struct timeval tp;
        gettimeofday(&tp, NULL);
        total_time += (double)tp.tv_sec + 1.0e-6 * (double)tp.tv_usec - stop_time;
#endif
        stop_state = true;
    }
}

void TimerRT::start() {
#ifdef _WIN32
    LARGE_INTEGER count;
    QueryPerformanceCounter(&count);
    start_count = count.QuadPart;
#else
    struct timeval tp;
    gettimeofday(&tp, NULL);
    stop_time = (double)tp.tv_sec + 1.0e-6 * (double)tp.tv_usec;
#endif
    stop_state = false;
}

double TimerRT::getTime() {
    stop();
    start();
    return total_time;
}

#endif // TimerRTH
