#pragma once
#ifndef GRAPH_H
#define GRAPH_H

#include <fstream>
#include <cassert>
#include <algorithm>
#include <random>
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
    uint32_t nodeNum;
    uint64_t edgeNum, maxKtrussEdge;
    uint32_t maxSup, minSup, maxCore;
    uint32_t zero_edge;
    uint64_t total_io;
    int memEdges,maxKtruss;
    int last_s;
    int last_sup;
    float memUsage;

    uint32_t *degree;
    uint32_t *degree_;
    uint32_t *prefix;
    int *bin;
    int *cntClass;
    bool *vis;
    bool isDeleteEdgeByFunc;
    uint64_t *edgeListBegPtr;
    uint64_t *edgeListBegPtrPlus;

    unordered_set<int>** m_dynamicDel;
    bool* m_delBit;

    unordered_set<int>** m_dynamicIns;
    bool* m_insBit;

    uint32_t Triangles;
    uint32_t nBins;
	int mp;
    string saveKtruss;

    std::ofstream fout;
    Clock graphClock_;
    std::unordered_set<uint64_t> isInDelQue;
    std::unordered_map<uint32_t,uint32_t> subToGlobal;
    std::unordered_map<uint32_t,uint32_t> trussClass;

    // std::unordered_map<uint64_t,uint64_t> delTwoVerMapEid;

public:
    struct DynamicHeap{
        int size;
        int cap;
        int min_sup;
        es *arr;
        std::unordered_map<uint64_t,uint64_t> hash;
        DynamicHeap(int sz){
            cap = sz;
            size = 0;
            arr = (es *)malloc(sz * sizeof(es));
        }
        ~DynamicHeap(){
            free(arr);
        }

        uint32_t getSup(uint64_t eid){
            if(hash.find(eid) != hash.end()){
                int i = hash[eid];
                return arr[i].sup;
            }
            else{
                printf("getSup function error\n");
                return 0;
            }
        }

        bool find(uint64_t eid){
            if(hash.find(eid) != hash.end())
                return true;
            return false;
        }

        void push(es blob){
            int i = size++;
            if(size > cap) printf("error");
            while(i && blob.sup < arr[PARENT(i)].sup){
                arr[i] = arr[PARENT(i)];
                hash[arr[i].eid] = i;
                i = PARENT(i);
            }
            arr[i] = blob;
            hash[arr[i].eid] = i;
            // print();
        }

        void heapify(int i) {
            int smallest = (LEFTCHILD(i) < size && arr[LEFTCHILD(i)].sup < arr[i].sup) ? LEFTCHILD(i) : i ;
            if(RIGHTCHILD(i) < size && arr[RIGHTCHILD(i)].sup < arr[smallest].sup) {
                smallest = RIGHTCHILD(i);
            }
            if(smallest != i) {
                es temp = arr[i];
                arr[i] = arr[smallest];
                hash[arr[i].eid] = i;
                arr[smallest] = temp;
                hash[temp.eid] = smallest;
                heapify(smallest) ;
            }
        }

        void erase(uint64_t eid){
            if(hash.find(eid) != hash.end()){
                int i = hash[eid];
                if(i < size-1){
                    arr[i].sup = 0;
                    es blob = arr[i];
                    while(i && blob.sup < arr[PARENT(i)].sup){
                        arr[i] = arr[PARENT(i)];
                        hash[arr[i].eid] = i;
                        i = PARENT(i);
                    }

                    arr[i] = arr[--size];
                    hash[arr[i].eid] = i;
                    heapify(i);
                }
                else
                    size--;
                hash.erase(eid);
            }
        }

        void adjust(uint64_t eid,int v){
            int i = hash[eid];
            arr[i].sup -= v;
            es blob = arr[i];
            while(i && blob.sup < arr[PARENT(i)].sup){
                arr[i] = arr[PARENT(i)];
                hash[arr[i].eid] = i;
                i = PARENT(i);
            }
            arr[i] = blob;
            hash[arr[i].eid] = i;
        }

        es pop(){
            es ret = arr[0];
            hash.erase(ret.eid);
            arr[0] = arr[--size];
            hash[arr[0].eid] = 0;
            heapify(0);
            return ret;
        }
        void print(){
            printf("-------------------\nDynamic Heap Display\n");
            for(int i = 0; i < size; i++){
                printf("i: %d, eid: %lu, sup: %u\n",i,arr[i].eid,arr[i].sup);
            }
        }
    };

public:
    Graph(int _nodeNum, int _edgeNum):graphClock_("graph process")
    {
        fout.open("result/out.txt");
        log_info(graphClock_.Start());
        Triangles = 0;
        zero_edge = 0;
        total_io = 0;
        memUsage = 0;
        maxSup = 0, maxCore = 0, maxKtruss = 0;
        minSup = _nodeNum;        
        nodeNum = _nodeNum;
        edgeNum = _edgeNum;
        
        nBins = 0;
	    mp = 0;
        last_s = 0;
        last_sup = 0;
        isDeleteEdgeByFunc = false;

        degree = (uint32_t *)malloc(sizeof(uint32_t) * nodeNum);
        degree_ = (uint32_t *)malloc(sizeof(uint32_t) * nodeNum);
        bin = (int*)malloc(sizeof(int)*nodeNum);
        prefix = (uint32_t*)malloc(sizeof(uint32_t)*nodeNum);
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

    void Intersect(std::vector<vertex_eid>& A, std::vector<vertex_eid>& B, std::vector<eid_eid>& common);
    void printClass(uint32_t u, uint32_t v, uint32_t cls);
    void changeEdgeSup(MyReadFile& fOff, MyReadFile& fSup, uint64_t eid, uint32_t& sup);
    void CoreTrussDecomPlus(readFile &file);
    void KCore(readFile &file, uint32_t *coreNum, bool isDynamicload);
    void reconCoreGraph(uint32_t *coreNum, readFile &file, readFile &newFile);

    void filterGlobalIntoSub(bool *isInSubG, readFile &file, readFile &newFile, bool isDynamicload);

    void binaryAndIncremental(readFile &file, uint32_t start);
    uint32_t binary(readFile &file, uint32_t start, uint32_t end);
    bool existTrussPlus(readFile &file, uint32_t &mid, uint32_t &TrussEdge, uint32_t &Truss);
    bool existTrussLazyUpdate(readFile &file, uint32_t &mid, uint32_t &TrussEdge, uint32_t &Truss, ListLinearHeapTruss *linear_heap, DynamicHeap &dheap);
    bool bottomUpDecom(readFile &file, uint32_t &mid, uint32_t &TrussEdge, uint32_t &Truss);

    bool inducedGraph(int mid, readFile &file, int &TrussEdge);
    
    void updateEdgeSSDNew(uint32_t u, uint32_t v, uint64_t eid, uint32_t uv_sup, uint32_t minsup, uint32_t *sup_arr,MyReadFile &fEdgePos, MyReadFile &fSup, MyReadFile &fBinEdge, MyReadFile &fOff);
    
    void CountTriangleSSD(readFile &file, bool saveSupEdge);
    void CountTriangleSSDPlus(readFile &file, bool isOrder = true);
    void CountTriangleSSDByDegOrder(readFile &file, bool saveSupEdge);

    void sortSupport(readFile &file);
    void binSortSSD(readFile &file);
    int IntersectTriangleSSDByDegOrder(uint32_t u, uint32_t u_nbrNum,uint32_t *nbr_u, uint32_t v, uint32_t v_nbrNum, uint32_t *nbr_v, 
    MyReadFile &fDat, MyReadFile &fEid, std::vector<eid_eid> &common, std::unordered_map<uint64_t,uint32_t> &map_pos);


    void IntersectTrussNew(uint32_t u, uint32_t *nbr_u, uint64_t *eid_u,  uint32_t *sup_u,
uint32_t v, uint32_t *nbr_v, uint64_t *eid_v, uint32_t *sup_v, MyReadFile &fDat, MyReadFile &fEid, MyReadFile &fSup,
uint32_t sup, vector<eid_eid>& comm, bool flag);

    bool peelVertex(int mid, readFile &file, bool flag = true);
    void constructBin(readFile &file);
    void trussDecomLinearList(readFile &file);
    void trussDecomLazyUpdate(readFile &file);
    void writeTwoVerFromEid(readFile &file);
    void prepareStage(readFile &file);

    bool inducedGraphUnOrder(uint32_t &left,uint32_t &right,uint32_t &mid, readFile &file, uint32_t &TrussEdge, uint32_t &Truss, bool delOnSubg);
    void InitialUnOrder(readFile &file);
    void InitialUnOrderDegSort(readFile &file);
    int IntersectTriangleSSDByDegOrderInSecondCore(uint32_t u, uint32_t u_nbrNum, uint32_t *nbr_u, uint32_t v, uint32_t v_nbrNum, uint32_t *nbr_v, 
MyReadFile &fDat, MyReadFile &fEid, std::vector<eid_eid> &common, std::unordered_map<uint64_t,uint32_t> &map_pos, bool *isInSubG);
    bool deleteEdge(readFile &file, readFile &newFile,uint32_t last_mid, uint32_t &mid, uint32_t *prefix, uint32_t &TrussDege, uint32_t &Truss);
    
    int getLast_sup(){
        return last_sup;
    }
    uint32_t getNodeNum(){
        return nodeNum;
    }
    uint32_t getTrussEdge(int i){
        return cntClass[i];
    }
    void updateNbrInDeleteEdgeLazyUpdate(uint32_t u, uint32_t *nbr_u, uint64_t *eid_u,  uint32_t *sup_u,
uint32_t v, uint32_t *nbr_v, uint64_t *eid_v, uint32_t *sup_v,
MyReadFile &fDat, MyReadFile &fEid, MyReadFile &fSup,MyReadFile &fOff, MyReadFile &fPres, MyReadFile &fNexts,
uint32_t sup, queue<EdgeSup> &del_que, uint32_t mid, DynamicHeap &dheap, ListLinearHeapTruss *linear_heap);
    void updateNbrInExistTrussLazyUpdate(uint32_t u, uint32_t *nbr_u, uint64_t *eid_u,  uint32_t *sup_u,
uint32_t v, uint32_t *nbr_v, uint64_t *eid_v, uint32_t *sup_v,
MyReadFile &fDat, MyReadFile &fEid, MyReadFile &fSup,MyReadFile &fOff, MyReadFile &fPres, MyReadFile &fNexts,
uint32_t sup, DynamicHeap &dheap, ListLinearHeapTruss *linear_heap);

    bool deleteEdgeLazyUpdateTest(readFile &file, readFile &newFile, uint32_t last_mid, uint32_t &mid, uint32_t *prefix_, uint32_t &TrussEdge, uint32_t &Truss, ListLinearHeapTruss *linear_heap, DynamicHeap &dheap);

    void dynamicMaxTrussMaintenance(readFile &file);
    void dynamicMaxTrussInsertion_YLJ(readFile &file);
    void dynamicMaxTrussDeletion_YLJ(readFile &file);
    void generateRandomEdges(readFile &file, uint32_t &generate_edge, Edge *dynamic, bool *isInMaxTruss, unordered_map<uint32_t, uint32_t>& globalToSub);
    void generateRandomEdgesInsertion(readFile &file, uint32_t &generate_edge, Edge *dynamic);

    uint32_t selectNbr(readFile &file, uint32_t a, int index);
    void reconMaxTrussGraph(readFile &file, readFile &newFile, bool firstUse = true);
    void recoveryKMaxTrussSub(readFile &file, unordered_map<uint32_t,uint32_t>& globalToSub, bool *verMaxKTrussSet, string dir);
    void recoveryKMaxTruss(readFile &file, unordered_map<uint32_t,uint32_t>& globalToSub, bool *verMaxKTrussSet, string dir);    
    void delEdgeDynamic(Edge del_e, readFile &newFile, uint32_t maxK, uint64_t &newFileEdge,
     unordered_map<uint32_t, uint32_t>& subToGlobal, unordered_map<uint64_t, uint64_t>& delTwoVerMapEid);
    void loadNbrForDynamic(uint32_t u, uint32_t* nbr, uint32_t& degree, uint64_t pos, MyReadFile& fDat, bool oriGra = true);
    void loadSupForDynamic(uint32_t u, uint32_t* nbr, uint64_t* eid_,uint32_t& _degree, uint64_t pos, MyReadFile& fDat, MyReadFile& fEid,
unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> &sup_dram,
unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> &eid_dram);

    void insEdgeDynamic(Edge e, readFile &file, unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> &sup_dram,
    unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> &eid_dram,
bool *verMaxKTruSet, vector<uint64_t> &begPtr, vector<uint32_t> &degSub);
    void insEdgeDynamicNew(Edge e, readFile &file, unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> &sup_dram,
    unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> &eid_dram,
    unordered_map<uint32_t,uint32_t> &subToGlobal);

    void delEdgeDynamicYLJ(Edge e, readFile &file, bool *verMaxKTruSet, vector<uint64_t> &begPtr, vector<uint32_t> &degSub);
    void loadNbrAndSupDynamic(uint32_t u, uint32_t* nbr_, uint32_t* sup_, uint64_t* eid_, uint32_t& _degree, 
uint64_t begPtr, MyReadFile& fDat, MyReadFile& fSup, MyReadFile& fEid, 
unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> &sup_dram, 
unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> &eid_dram);
    void loadNbrAndSupDynamicNew(uint32_t u, uint32_t* nbr_, uint32_t* sup_, uint64_t* eid_, uint32_t& _degree, 
uint64_t begPtr, MyReadFile& fDat, MyReadFile& fSup, MyReadFile& fEid);
    void IntersectOperaInsDynamic(uint32_t u, uint32_t* nbr_u, uint32_t* sup_u, uint64_t* eid_u, uint32_t u_degree, 
uint32_t v, uint32_t* nbr_v, uint32_t* sup_v, uint64_t* eid_v, uint32_t v_degree,uint32_t &max_sup,std::vector<ver_eid_eid> &comm);
    void IntersectOperaDelDynamic(uint32_t u, uint32_t* nbr_u, uint32_t* sup_u, uint64_t* eid_u, uint32_t u_degree, 
uint32_t v, uint32_t* nbr_v, uint32_t* sup_v, uint64_t* eid_v, uint32_t v_degree,std::vector<ver_eid_eid> &comm);

    void UpdateCoreDynamic(readFile &file, uint32_t *coreNum, uint32_t u, uint32_t v);


};

#endif