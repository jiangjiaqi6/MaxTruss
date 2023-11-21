#ifndef FILELOADER
#define FILELOADER
#include <cstdlib>
#include <unistd.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/types.h>  
#include <sys/socket.h>  
#include <netinet/in.h>  
#include <arpa/inet.h>  
#include <unistd.h>
#include "edge.h"  
#include "clock.h"
#include "log.h"
#include "MyFile.h"
#include <unordered_map>
#include <unordered_set>
#include <functional>
#include <set>
#include <map>
#include <algorithm>
#include <queue>

class readFile{
public:
    std::string m_idx;
	std::string m_dat;
	std::string m_base;
	std::string m_info;
    std::string m_eid;
    std::string m_supp;
    std::string m_suppSort;
    std::string m_binEdge;
    std::string m_ePos;
    std::string m_eidToVer;
    std::string m_offset;
    uint64_t len,fd,write_io;
    char* addr;
    int memEdges,buffersize;
    uint32_t *m_vertexMap;
    int m_maxID;
    uint32_t verNum;
    uint32_t maxDeg;
    uint64_t edgeNum;
    std::vector<uint32_t> vertexId;
    Clock readFileClock_;

    std::unordered_map<uint32_t,uint32_t> verDegMap;
public:
    
    readFile(std::string base):readFileClock_("ReadFile"){
        m_base = base;
        m_dat = base+"graph.dat";
        m_idx = base+"graph.idx";
        m_info = base+"graph.info";
        m_eid = base+"graph.eid";
        m_supp = base+"graph.support";
        m_suppSort = base+"graph.supportSorted";
        m_binEdge = base+"graph.binEdge";
        m_ePos = base+"graph.ePosInBin";
        m_eidToVer = base+"graph.eidToVertex";
        m_offset = base+"graph.eidToOffset";
        len = 0, fd = -1, write_io = 0;
        addr = NULL;
        memEdges = 500000000;   //the number of edges loaded in dram
        m_maxID = 150000000;   
        buffersize = 100*1024*1024;
        m_vertexMap = new uint32_t[m_maxID];
        log_debug(readFileClock_.Start());
    }
    
    ~readFile(){
        if(m_vertexMap != NULL)
            delete[] m_vertexMap;
        // printf("readfile finish\n");
    }
    void change(std::string base){
        m_base = base;
        m_dat = base+"graph.dat";
        m_idx = base+"graph.idx";
        m_info = base+"graph.info";
        m_eid = base+"graph.eid";
        m_supp = base+"graph.support";
        m_suppSort = base+"graph.supportSorted";
        m_binEdge = base+"graph.binEdge";
        m_ePos = base+"graph.ePosInBin";
        m_eidToVer = base+"graph.eidToVertex";
        m_offset = base+"graph.eidToOffset";
        len = 0, fd = -1, write_io = 0;
    }

    struct stat statbuf;
    void loadFile(const char* path){
        fd = open(path,O_RDONLY);
        if (fd < 0)
        {
            puts("can not open file");
            exit(-1);
        }
        int ret = fstat(fd, &statbuf);
        if (ret < 0)
        {
            puts("can not fstat");
            exit(-1);
        }
        len = statbuf.st_size;
        addr = (char*)mmap(0, len, PROT_READ, MAP_SHARED, fd, 0);
        if (addr == (void *)-1)
        {
            puts("can not mmap");
            exit(-1);
        }

    }
    void release(){
        if(fd != -1)
        {
            close(fd);
            munmap(addr, len);
        }
    }
    void initialGraphInfo(){
        MyReadFile fInfo( m_info );
	    fInfo.fopen( BUFFERED );
        uint32_t n,deg;
        uint64_t m;
        fInfo.fread(&n,sizeof(uint32_t));
        fInfo.fread(&deg,sizeof(uint32_t));
        fInfo.fread(&m,sizeof(uint64_t));
        verNum = n, maxDeg = deg, edgeNum = m;
        // printf("n: %u, deg: %u, m: %u\n",n,deg,m);
        fInfo.fclose();   
    }
    uint64_t getLen()
    {
        return len;
    }

    char *getAddr()
    {
        return addr;
    }

    uint32_t getVerNum(){
        return verNum;
    }
    uint64_t getEdgeNum(){
        return edgeNum;
    }

    uint32_t getOriginVer(uint32_t v){
        return m_vertexMap[v];
    }

    uint32_t getVertexID(uint32_t u,int& num);

    void createDir(string &s);
    void LoadGraph();
    void moduleInSaveEdges(uint32_t &u,uint32_t &v,int &num,TEdge* edges,
        unsigned long  &size,unsigned long &es,int &tmpFile,string &sub_dir);

    void moduleInSaveEdgesUnOrder(uint32_t &u,uint32_t &v,int &num,TEdge* edges,
        unsigned long  &size,uint32_t es,int &tmpFile,string &sub_dir);
    // void saveTmpEdges(TEdge* edges,int size,int tmpFile);
    template<typename T> 
    void saveTmpEdges(T* edges, int size, char fileName[], std::function<bool(const T &, const T &)> cmp){

        std::sort(edges,edges+size,cmp);
        log_debug(readFileClock_.Count(fileName));

        FILE* fo = fopen(fileName,"wb");
        fwrite( edges, sizeof(T), size, fo );
        write_io += 1;
        fclose(fo);
    }
    void partitionGraphByVertex(int percentage);
    void partitionGraphByEdge(int percentage);

    bool static edgeCompare(const TEdge &e1, const TEdge &e2);
    void merge(int size, uint32_t &vertexNum, uint32_t &maxDegree, string &name, bool flag);
    void mergeByDegSort(int size, uint32_t &vertexNum, uint32_t &maxDegree, string &name, bool flag, bool isBigG);
    bool mergeFinished(TEdge* es, int size);
    int minPartition(TEdge* es, int size);
    // void mergeBySup(int size);

    template<typename T> 
    void mergeBySup(int size){
        FILE **frl = new FILE*[size];
        FILE *fEdgeSup = fopen(m_suppSort.c_str(),"wb");
        char filename[100];
        T tmp;

        struct cmp{
            bool operator()(pair<int,T> a, pair<int,T> b){
                return a.second.sup > b.second.sup;
            }
        };
        std::priority_queue<pair<int,T>,std::vector<pair<int,T>>, cmp> pq; 
        string sub_dir = m_base + "support_sort_edge";

        for(int i = 0; i < size; i++){
            sprintf(filename,"%s/edges_tmp_%d",sub_dir.c_str(),i);
            frl[i] = fopen(filename,"rb");
            fread(&tmp,sizeof(T),1,frl[i]);
            pq.push({i,tmp});
        }
        write_io += size;
        int count =  0;
        while(!pq.empty()){
            int index = pq.top().first;
            tmp = pq.top().second;
            pq.pop();
            // printf("count: %d, u: %u, v: %u, sup: %u\n",count,tmp.u,tmp.v,tmp.sup);
            fwrite(&tmp,sizeof(T),1,fEdgeSup);
            write_io += 1;
            count++;
            if(fread(&tmp,sizeof(T),1,frl[index])){
                pq.push({index,tmp});
            }
        }

        for (int i = 0; i < size ; ++i)
            fclose(frl[i]);
        fclose(fEdgeSup);
        delete[] frl;
    }

};


#endif