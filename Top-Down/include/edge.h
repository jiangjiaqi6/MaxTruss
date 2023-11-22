#ifndef EDGE_H
#define EDGE_H
#include <iostream>
#include <vector>
#include <cstring>
#include "util.h"

struct Edge {
    uint32_t u;
    uint32_t v;
} __attribute__ ((aligned (4)));

typedef struct vertex_eid{
    uint32_t v;
    uint64_t eid;
}vertex_eid;
typedef struct eid_eid{
    uint64_t first;
    uint64_t second;
    eid_eid(){
        first = 0;second = 0;
    }
    eid_eid(uint64_t f, uint64_t s){
        first = f;second = s;
    }
}eid_eid;

typedef struct offset{
    uint64_t start;
    uint64_t end;
    offset(){
        start = 0;end = 0;
    }
    offset(uint64_t f, uint64_t s){
        start = f;end = s;
    }
}offset;

typedef struct TEdge{
    uint32_t u;
    uint32_t v;
    uint64_t eid;
}TEdge;

typedef struct Edge_pos{
    uint32_t u;
    uint32_t v;
    uint64_t eid;
    uint64_t pos;
}Edge_pos;

typedef struct EdgeSup{
    uint32_t u;
    uint32_t v;
    uint32_t sup;
    uint64_t eid;
}EdgeSup;

typedef struct uvSup{
    uint32_t u;
    uint32_t v;
    uint32_t sup;
}uvSup;

typedef struct ver_eid_eid{
    uint32_t w;
    uint64_t first;  //corresponding to the vertex u
    uint32_t u_sup;
    uint64_t second; //corresponding to the vertex v
    uint32_t v_sup;
    
}ver_eid_eid;

typedef struct ver_two_eid{
    uint32_t w;
    uint64_t first;  //corresponding to the vertex u
    uint64_t second; //corresponding to the vertex v
    
}ver_two_eid;

typedef struct es{
    uint64_t eid;
    uint32_t sup;
    es(){
        eid = 0;sup = 0;
    }
    es(uint64_t _eid, uint32_t _sup){
        eid = _eid;sup = _sup;
    }
}es;

#define LEFTCHILD(x) 2 * x + 1
#define RIGHTCHILD(x) 2 * x + 2
#define PARENT(x) (x - 1) / 2


bool operator <(const Edge &e1, const Edge &e2);
void ReadBaseLine(const char * path,  std::vector<Edge> & Edges);
uint32_t LoadEdge(const char * const bptr, const uint32_t len, std::vector<Edge> & Edges);
uint32_t mysscanf(const char * const eptr, const char *& cptr);
void sscanfSkip(const char * const eptr, const char *& cptr);

#endif
