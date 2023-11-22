#ifndef UTIL_H
#define UTIL_H

#include <omp.h>
#include <vector>
#include <sys/stat.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <unistd.h>
#include <iostream>
#include <fstream> 

class uvList
{
public:
    uint32_t u;
    uint32_t v;
    uint32_t uvId;
    uvList (uint32_t _u, uint32_t _v, uint32_t _uvId) : u(_u), v(_v), uvId(_uvId) {} 
};

// get current process pid
inline int GetCurrentPid()
{
    return getpid();
}
 
// FIXME: can also get cpu and mem status from popen cmd
// the info line num in /proc/{pid}/status file
#define VMRSS_LINE 22
#define PROCESS_ITEM 14

// get specific process physical memeory occupation size by pid (MB)
inline float GetMemoryUsage(int pid)
{
    char file_name[64] = { 0 };
    FILE* fd;
    char line_buff[512] = { 0 };
    sprintf(file_name, "/proc/%d/status", pid);
 
    fd = fopen(file_name, "r");
    if (nullptr == fd)
        return 0;
 
    char name[64];
    int vmrss = 0;
    for (int i = 0; i < VMRSS_LINE - 1; i++)
        fgets(line_buff, sizeof(line_buff), fd);
 
    fgets(line_buff, sizeof(line_buff), fd);
    sscanf(line_buff, "%s %d", name, &vmrss);
    fclose(fd);
 
    // cnvert VmRSS from KB to MB
    return vmrss / 1024.0;
}

bool orderComp (const uint32_t& u, const uint32_t& v);

void parmemset(void *const ptr, int const c, uint32_t const len);


#endif