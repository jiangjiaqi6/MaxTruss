#ifndef TIMER
#define TIMER

#include <iostream>
#include <chrono>
#include <iomanip>
#include "tools.h"
#include <sys/time.h>

class timeCount{
private:
    uint64_t _startTime;
    uint64_t _endTime;
    uint64_t _lastTime;
    bool _timeState;

public:
//explicit transform
    explicit timeCount(){
        _startTime = _endTime = _lastTime = 0;
        _timeState = false;
    }
    uint64_t GetSysTime(){
        return std::chrono::duration_cast<std::chrono::microseconds>(std::chrono::high_resolution_clock::now().time_since_epoch()).count();
    }
    void StartTime(){
        _timeState = true;
        _startTime = _lastTime = GetSysTime();
    }
    void EndTime(){
        _timeState = false;
        _endTime = GetSysTime();
    }
    void UpdateLastTime()
    {
        _lastTime = GetSysTime();
    }
    void Gaptime(const char* str){
        uint64_t curTime = GetSysTime();
        std::cout << CYAN << "timeCount Info" << RESET;
        std::cout << std::setfill('-') << std::setw(35) << str << " : "<< \
        std::fixed << std::setprecision(4) <<  static_cast<double>(curTime - _lastTime)/ 1e6 << "s";
        std::cout << " Total Time Usage: " << static_cast<double>(curTime - _startTime) / 1e6 << "s" << std::endl;
        _lastTime = curTime;
        std::cout << std::setfill(' ');
    }
    double QueryTime(){
        return static_cast<double>(_endTime - _startTime) / 1e6;
    }

};
#endif