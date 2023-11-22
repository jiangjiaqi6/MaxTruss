#pragma once
#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>

#include <cassert>
#include <cstring>
#include <algorithm>
#include <set>
#include <map>
#include <tuple>
#include "edge.h"
#include "tools.h"
#include "util.h"
#include "fileLoader.h"
#include "LinearHeapTrussSSD.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_set>




class Graph
{
private:
    int method;
    uint32_t nodeNum;
    uint64_t edgeNum;
    uint32_t maxSup;
    uint32_t minSup;
    uint32_t maxCore;
    uint32_t MEM;
    uint32_t topdown_k;
    int memEdges;
    int last_s;
    int last_sup;
    uint64_t total_io;
    float memUsage;

    std::unordered_map<uint64_t,uint32_t> supports;
    std::vector<std::map<uint32_t,uint32_t>> pos;

    uint32_t *degree;
    uint32_t *degree_;
    int *prefix;
    int *bin;
    int *cntClass;
    bool *vis;
    uint64_t *edgeListBegPtr;
    uint64_t *edgeListBegPtrPlus;

    uint32_t Triangles;
    uint32_t nBins;
	int mp;
    TEdge *binEdge;

    std::ofstream fout;
    Clock graphClock_;

public:
    Graph(int _nodeNum, int _edgeNum):graphClock_("graph process")
    {
        fout.open("result/out.txt");
        log_info(graphClock_.Start());
        Triangles = 0;
        method = 1;
        maxSup = 0;
        maxCore = 0;
        minSup = _nodeNum;
        memEdges = 10000000;   //the number of edges loaded in dram
        MEM = 800000000;
        topdown_k = 0;
        total_io = 0;
        memUsage = 0;
        
        nodeNum = _nodeNum;
        edgeNum = _edgeNum;
        
        nBins = 0;
	    mp = 0;
        last_s = 0;
        last_sup = 0;

        degree = (uint32_t *)malloc(sizeof(uint32_t) * nodeNum);
        degree_ = (uint32_t *)malloc(sizeof(uint32_t) * nodeNum);
        bin = (int*)malloc(sizeof(int)*nodeNum);
        prefix = (int*)malloc(sizeof(int)*nodeNum);
        cntClass = (int*)malloc(sizeof(int)*nodeNum);
        edgeListBegPtr = (uint64_t *)malloc(sizeof(uint64_t) * (nodeNum + 1)); 
        #ifdef DegSort 
        edgeListBegPtrPlus = (uint64_t *)malloc(sizeof(uint64_t) * (nodeNum + 1));
        memset(edgeListBegPtrPlus, 0, sizeof(uint64_t) * (nodeNum + 1));  
        #endif
        vis = (bool *)calloc(nodeNum,sizeof(bool));      

        memset(edgeListBegPtr, 0, sizeof(uint64_t) * (nodeNum + 1));
        memset(degree, 0, sizeof(uint32_t) * nodeNum);
        memset(degree_, 0, sizeof(uint32_t) * nodeNum);
        memset(bin, 0, sizeof(int) * nodeNum);
        memset(prefix, 0, sizeof(int) * nodeNum);
        memset(cntClass, 0, sizeof(int) * nodeNum);

    }
    ~Graph()
    {
        fout.close();
        free(edgeListBegPtr);
        #ifdef DegSort
        free(edgeListBegPtrPlus);
        #endif
        free(degree);
        free(degree_);
        free(bin);
        free(prefix);
        free(cntClass);
        free(vis);
    }

    template<typename T>
    void loadInfo(T* nbr, uint32_t num, uint64_t pos, MyReadFile& fDat);
    void Initial(readFile &file);
    void printClass(uint32_t u, uint32_t v, uint32_t cls);

    
    void WCUpperBound(readFile &file);  //semi-external method
    void upperBounding(readFile &file); //partition method
    void sortUpperBound(readFile &file);
    
    int WCHIndex(uint32_t *sup_nbr, uint64_t *eid_nbr,int count, uint64_t mark);
    int UpperBoundHIndex(uint64_t *eid_nbr,uint32_t *sup_arr,uint32_t *binCount,uint64_t start, uint64_t end, uint64_t mark, MyReadFile &fSup);
    void triangleListExtendedGraph(readFile &file);
    void processTriangleInMem(uint64_t &pos,uint32_t end,uint64_t &index,uint32_t &tmp_node, uint64_t &tmp_edge, uint32_t &triangle, uint32_t &size, uint32_t iter,
uint32_t *nbr, uint32_t *nbr_global, uint64_t *eid, uint64_t *eid_global, uint32_t *deg_global, uint32_t *ver,uint32_t *ver_global, offset *off,std::unordered_map<uint32_t,uint32_t> &verToMap, 
std::unordered_map<uint32_t,uint32_t> &verToMapNew, std::unordered_map<uint32_t,bool> &isInVerMap, std::unordered_map<uint64_t,uint32_t> &eidInTri, 
bool *isInVer, MyReadFile &fDat, FILE *fo, MyReadFile &fEid, FILE *fo_eid, MyReadFile &fOff, 
bool &writeIntoSSD, bool flag, FILE *fSupp);
    void processUBInMem(uint64_t &pos,uint32_t end,uint64_t &index,uint32_t *nbr,uint64_t *eid,uint32_t *ver,uint32_t *binCount, uint32_t *sup_arr, offset *off,std::unordered_map<uint32_t,bool> &isInVerMap,
    MyReadFile &fDat,MyReadFile &fEid,MyReadFile &fSup,FILE *fUB);

    void topDownTrussDecom(readFile &file);
    void writeTwoVerFromEid(readFile &file);
    void extractSubgraph(readFile &file);
    void updateEdgeDram(uint32_t u, uint32_t v, uint64_t eid, uint32_t minsup,std::unordered_map<uint32_t,uint32_t> &verToMap);


};

#endif