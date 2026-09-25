/* Monotonic wall-clock timer and cycle counter. */
#ifndef TV_TIMER_H
#define TV_TIMER_H

#include <stdint.h>

#ifdef _WIN32
#include <windows.h>
static inline double now_ms(void) {
    static LARGE_INTEGER f;
    static int init = 0;
    if (!init) { QueryPerformanceFrequency(&f); init = 1; }
    LARGE_INTEGER c;
    QueryPerformanceCounter(&c);
    return (double)c.QuadPart * 1000.0 / (double)f.QuadPart;
}
#else
#include <time.h>
static inline double now_ms(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e3 + (double)ts.tv_nsec / 1e6;
}
#endif

#if defined(__x86_64__) || defined(_M_X64)
#include <x86intrin.h>
static inline uint64_t cycles(void) { return __rdtsc(); }
#else
static inline uint64_t cycles(void) { return 0; }
#endif

#endif
