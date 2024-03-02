#include "include/graph.h"
#include "include/Time.h"
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <queue>
#include <algorithm>

inline void Graph::printClass(uint32_t u, uint32_t v, uint32_t cls) {
	++cntClass[cls];
	// fout << "(" << u << "," << v << "):" << cls << std::endl;
}
template<typename T>
void Graph::loadInfo(T* nbr, uint32_t num, uint64_t pos, MyReadFile& fDat){
    fDat.fseek(pos);
    fDat.fread(nbr,sizeof(T)*num);
}


//load O(n) array, such as, offset array and degree array
void Graph::Initial(readFile &file){
    MyReadFile fIdx( file.m_idx );
	fIdx.fopen( BUFFERED );
	MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    uint64_t tmpa = 0,tmpb = 0;
    uint32_t degreeTmp = 0;
    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    for (uint32_t i = 0; i < nodeNum; ++i){
		fIdx.fread(&tmpa,sizeof(uint64_t));
        edgeListBegPtr[i] = tmpa;
        #ifdef DegSort
        fIdx.fread(&tmpb,sizeof(uint64_t));
        edgeListBegPtrPlus[i] = tmpb;
        #endif
		fIdx.fread(&degreeTmp,sizeof(uint32_t));
        degree[i] = degreeTmp;
        #ifdef DegSort
        // printf("ver: %d, deg: %d, pos: %lu, posPlus: %lu\n",i, degreeTmp, tmpa,tmpb);
        // uint32_t num = degree[i]-(edgeListBegPtrPlus[i]-edgeListBegPtr[i])/sizeof(uint32_t);
        // if(num == 0)continue;
        // loadInfo(nbr_u,num,edgeListBegPtrPlus[i],fDat);
        // uint32_t u = i;
        // for(uint32_t j = 0; j < num; j++){
        //     uint32_t v = nbr_u[j];
        //     printf("u: %d, v: %d\n",u,v);
        // }
        #endif
  	}
    total_io += fIdx.get_total_io();
    free(nbr_u);
    fDat.fclose();
	fIdx.fclose();
}
void Graph::WCUpperBound(readFile &file){
    MyReadFile fSup( file.m_supp );
	fSup.fopen( NOBUFFER );
    MyReadFile fEid( file.m_eid );
	fEid.fopen( BUFFERED );
    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fBinEdge( file.m_binEdge );
	fBinEdge.fopen( NOBUFFER );
    MyReadFile fEdgePos( file.m_ePos );
	fEdgePos.fopen( NOBUFFER );

    MyReadFile fOff( file.m_offset );
	fOff.fopen( BUFFERED );

    uint32_t sup,cur_sup = uint32_t(-1);

    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* nbr_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    

    uint64_t* eid_u = (uint64_t *)malloc(sizeof(uint64_t) * file.maxDeg);
    uint64_t* eid_v = (uint64_t *)malloc(sizeof(uint64_t) * file.maxDeg);

    uint32_t* sup_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* sup_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    TEdge edge;
    uint32_t u,v;
    uint64_t eid_uv,uv_pos, offset, eid_copy;
    int max_upperBound = 0;
    for(uint32_t i = 0; i < nodeNum; i++){
        loadInfo(nbr_u,degree[i],edgeListBegPtr[i],fDat);
        uint32_t u = i;
        loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
        loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);
        for(uint32_t j = 0; j < degree[i]; j++){
            uint32_t v = nbr_u[j];
            if(v<u)continue;
            loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
            loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);
            int upper_u = WCHIndex(sup_u,eid_u,degree[i],eid_u[j]);
            int upper_v = WCHIndex(sup_v,eid_v,degree[v],eid_u[j]);
            int min_uv = std::min(upper_u,upper_v);
            max_upperBound = std::max(min_uv,max_upperBound);
            min_uv = min_uv > sup_u[j] ? sup_u[j] : min_uv;
        }
    }
    printf("upper bound: %d\n",max_upperBound);
    fSup.fclose();
    fEid.fclose();
    fDat.fclose();
    fEdgePos.fclose();
    fBinEdge.fclose();
    fOff.fclose();
    free(nbr_u);
    free(nbr_v);
    free(eid_u);
    free(eid_v);
    free(sup_u);
    free(sup_v);
}

int Graph::WCHIndex(uint32_t *sup_nbr, uint64_t *eid_nbr, int count, uint64_t mark){
    int *binCount = new int[maxSup+1]();
    for(int i = 0; i < count; i++){
        if(eid_nbr[i] == mark)continue;
        binCount[sup_nbr[i]]++;
    }
    int incident_max = 0;
    for(int i = maxSup; i >= 0; i--){
        incident_max += binCount[i];
        if(incident_max >= i)
            break;
    }
    
    delete[] binCount;
    return incident_max;
}

int Graph::UpperBoundHIndex(uint64_t *eid_nbr,uint32_t *sup_arr,uint32_t *binCount,uint64_t start, uint64_t end, uint64_t mark, MyReadFile &fSup){
    memset(binCount,0,sizeof(uint32_t)*(maxSup+1));
    uint32_t sup;
    for(int i = start; i < end; i++){
        if(eid_nbr[i] == mark)continue;
        // fSup.fseek(eid_nbr[i]*sizeof(uint32_t));
        // fSup.fread(&sup,sizeof(uint32_t));
        sup = sup_arr[i];
        if(sup > maxSup)
            printf("error\n");
        binCount[sup]++;
    }
    int incident_max = 0;
    for(int i = maxSup; i >= 0; i--){
        incident_max += binCount[i];
        if(incident_max >= i)
            break;
    }
    
    return incident_max;    
}

void Graph::processTriangleInMem(uint64_t &pos, uint32_t end, uint64_t &index,uint32_t &tmp_node, uint64_t &tmp_edge,uint32_t &triangle, uint32_t &size, uint32_t iter,
uint32_t *nbr, uint32_t *nbr_global, uint64_t *eid, uint64_t *eid_global, uint32_t *deg_global, uint32_t *ver,uint32_t *ver_global,offset *off, std::unordered_map<uint32_t,uint32_t> &verToMap, 
std::unordered_map<uint32_t,uint32_t> &verToMapNew, std::unordered_map<uint32_t,bool> &isInVerMap, std::unordered_map<uint64_t,uint32_t> &eidInTri, 
bool *isInVer, MyReadFile &fDat, FILE *fo, MyReadFile &fEid, FILE *fo_eid,  MyReadFile &fOff,
bool &writeIntoSSD, bool flag, FILE *fSupp)
{
    uint32_t nbr_num = (off[end-1].end-pos)/sizeof(uint32_t);
    loadInfo(nbr,nbr_num,pos,fDat);
    loadInfo(eid,nbr_num,pos*2,fEid);

    // operation in memory
    for(int i = 0; i < index; i++){
        uint32_t u = ver[i],tmp_node_copy = tmp_node;
        bool exist_u = false;

        for(int j = 0; j < (off[verToMap[u]].end - off[verToMap[u]].start) / sizeof(uint32_t); j++){
            uint32_t v = nbr[(off[verToMap[u]].start - pos) / sizeof(uint32_t) + j];
            uint64_t eid_uv = eid[(off[verToMap[u]].start - pos)*2 / sizeof(uint64_t) + j];
            
            if(isInVerMap.find(v) == isInVerMap.end())
            {
                if(!exist_u)
                {
                    ver_global[tmp_node] = u;
                    verToMapNew[u] = tmp_node;
                    exist_u = true;
                    tmp_node++;
                }
                if(size >= MEM){
                    fwrite( nbr_global, sizeof(uint32_t), size, fo );
                    fwrite( eid_global, sizeof(uint64_t), size, fo_eid );
                    size = 0;
                    writeIntoSSD = true;
                    total_io += 2;
                }
                tmp_edge++;
                nbr_global[size] = v;
                eid_global[size++] = eid_uv;
                deg_global[tmp_node_copy]++;
                continue;
            }
            if(v < u) continue;
            uint32_t sup = 0, sup_copy = 0;
            uint64_t ptr_u = (off[verToMap[u]].start - pos) / sizeof(uint32_t), upper_u = (off[verToMap[u]].end - pos) / sizeof(uint32_t);
            uint64_t ptr_v = (off[verToMap[v]].start - pos) / sizeof(uint32_t), upper_v = (off[verToMap[v]].end - pos) / sizeof(uint32_t);
            uint64_t eid_ptr_u = (off[verToMap[u]].start - pos)*2 / sizeof(uint64_t), eid_upper_u = (off[verToMap[u]].end - pos)*2 / sizeof(uint64_t);
            uint64_t eid_ptr_v = (off[verToMap[v]].start - pos)*2 / sizeof(uint64_t), eid_upper_v = (off[verToMap[v]].end - pos)*2 / sizeof(uint64_t);
            while(ptr_u < upper_u && ptr_v < upper_v){
                if(nbr[ptr_u] < nbr[ptr_v]) ptr_u++,eid_ptr_u++;
                else if(nbr[ptr_u] > nbr[ptr_v]) ptr_v++,eid_ptr_v++;
                else{
                    // if(!isInVer[nbr[ptr_u]] || (nbr[ptr_u] > v && nbr[ptr_v] > v))
                    if(isInVerMap.find(nbr[ptr_u]) == isInVerMap.end()){
                        sup++;
                        eidInTri[eid[eid_ptr_u]] += 1;
                        eidInTri[eid[eid_ptr_v]] += 1;
                    }  
                    else if((nbr[ptr_u] > v && nbr[ptr_v] > v))
                    {
                        sup++;
                    }
                    ptr_u++,eid_ptr_u++;
                    ptr_v++,eid_ptr_v++;
                    sup_copy++;
                }
            }

            triangle += sup;
            if(eidInTri.find(eid_uv) != eidInTri.end()){
                sup_copy += eidInTri[eid_uv];
                eidInTri.erase(eid_uv);
            }
            maxSup = std::max(sup_copy,maxSup);
            eid_eid tmp;
            fOff.fseek(eid_uv*sizeof(eid_eid));
            fOff.fread(&tmp,sizeof(eid_eid));
            fseek(fSupp,tmp.first*sizeof(uint32_t),SEEK_SET);
            fwrite(&sup_copy,sizeof(uint32_t),1,fSupp);
            fseek(fSupp,tmp.second*sizeof(uint32_t),SEEK_SET);
            fwrite(&sup_copy,sizeof(uint32_t),1,fSupp);
            total_io += 2;
            // printf("u: %u, v: %u, eid_uv: %u, sup: %u\n",u,v,eid_uv,sup_copy);
        }
    }
    pos = off[end].start;
    if(writeIntoSSD && flag){
        fwrite( nbr_global, sizeof(uint32_t), size, fo );
        fwrite( eid_global, sizeof(uint64_t), size, fo_eid );
        total_io += 2;
    }
    index = 0;
    memset(isInVer,false,sizeof(bool)*nodeNum);
    memset(ver,0,sizeof(uint32_t)*nodeNum);
    isInVerMap.clear();
}

void Graph::processUBInMem(uint64_t &pos,uint32_t end,uint64_t &index,uint32_t *nbr,uint64_t *eid,uint32_t *ver, uint32_t *sup_arr, uint32_t *binCount,
offset *off,std::unordered_map<uint32_t,bool> &isInVerMap,
    MyReadFile &fDat,MyReadFile &fEid,MyReadFile &fSup,FILE *fUB){
    uint32_t nbr_num = (off[end-1].end-pos)/sizeof(uint32_t);
    loadInfo(nbr,nbr_num,pos,fDat);
    loadInfo(sup_arr,nbr_num,pos,fSup);
    loadInfo(eid,nbr_num,pos*2,fEid);

    // operation in memory
    for(int i = 0; i < index; i++){
        uint32_t u = ver[i];
        bool exist_u = false;

        for(int j = 0; j < (off[u].end - off[u].start) / sizeof(uint32_t); j++){
            uint32_t v = nbr[(off[u].start - pos) / sizeof(uint32_t) + j];
            uint64_t eid_uv = eid[(off[u].start - pos)*2 / sizeof(uint64_t) + j];
            uint32_t sup = sup_arr[(off[u].start - pos) / sizeof(uint32_t) + j];
            // printf("u: %u, v: %u, eid: %lu, sup: %u\n",u,v,eid_uv,sup);
            uint32_t sup_copy = 0, min_uv = 0;
            
            if(isInVerMap.find(v) == isInVerMap.end())
            {
                if(u < v){
                    uint64_t ptr_u = (off[u].start - pos) / sizeof(uint32_t), upper_u = (off[u].end - pos) / sizeof(uint32_t);
                    int x_u = UpperBoundHIndex(eid,sup_arr,binCount,ptr_u,upper_u,eid_uv,fSup);
                    min_uv = x_u < sup ? x_u : sup;
                }
                else{
                    fseek(fUB,eid_uv*sizeof(uint32_t),SEEK_SET);
                    fread(&min_uv,sizeof(uint32_t),1,fUB);
                    total_io += 1;
                    // printf("min_uv: %u\n",min_uv);
                    uint64_t ptr_u = (off[u].start - pos) / sizeof(uint32_t), upper_u = (off[u].end - pos) / sizeof(uint32_t);
                    int x_u = UpperBoundHIndex(eid,sup_arr,binCount,ptr_u,upper_u,eid_uv,fSup);

                    if(x_u < min_uv){
                        min_uv = x_u;
                        topdown_k = std::max(min_uv,topdown_k);
                    }
                    else{
                        topdown_k = std::max(min_uv,topdown_k);
                        continue;
                    }
                }
                fseek(fUB,eid_uv*sizeof(uint32_t),SEEK_SET);
                fwrite(&min_uv,sizeof(uint32_t),1,fUB);
                total_io += 1;
			    // fout << "#eid " << eid_uv << "* min_uv: " << min_uv << std::endl;
                continue;
            }
            if(v < u) continue;
            uint64_t ptr_u = (off[u].start - pos) / sizeof(uint32_t), upper_u = (off[u].end - pos) / sizeof(uint32_t);
            int x_u = UpperBoundHIndex(eid,sup_arr,binCount,ptr_u,upper_u,eid_uv,fSup);

            uint64_t ptr_v = (off[v].start - pos) / sizeof(uint32_t), upper_v = (off[v].end - pos) / sizeof(uint32_t);
            int x_v = UpperBoundHIndex(eid,sup_arr,binCount,ptr_v,upper_v,eid_uv,fSup);

            min_uv = std::min(x_u,x_v);
            min_uv = min_uv < sup ? min_uv : sup;

            fseek(fUB,eid_uv*sizeof(uint32_t),SEEK_SET);
            fwrite(&min_uv,sizeof(uint32_t),1,fUB);
            total_io += 1;
            topdown_k = std::max(min_uv,topdown_k);
            if(sup >= 2207)
                fout << "#eid " << eid_uv << "- ub_uv: " << min_uv << std::endl;
        }
    }
    pos = off[end].start;
    index = 0;
    memset(ver,0,sizeof(uint32_t)*nodeNum);
    isInVerMap.clear();    
}

void Graph::upperBounding(readFile &file){
    log_debug(graphClock_.Count("upperBounding begin"));
    FILE* fUB = fopen(file.m_upperBound.c_str(),"w+");
    uint32_t *binCount = new uint32_t[maxSup+1]();
    MyReadFile fDat( file.m_dat );
    fDat.fopen( BUFFERED );

    MyReadFile fEid( file.m_eid );
    fEid.fopen( BUFFERED );

    MyReadFile fSupp( file.m_supp );
    fSupp.fopen( BUFFERED );

    uint32_t node_copy = nodeNum, tmp_node = 0, tmp_file = 0, size = 0, triangle = 0;
    uint64_t pos = 0, end = 0, index = 0, tmp_edge = edgeNum*2, last_tmp_edge = tmp_edge;
    uint32_t *nbr = new uint32_t[MEM]();
    uint32_t *sup_arr = new uint32_t[MEM]();
    uint64_t *eid = new uint64_t[MEM]();
    uint32_t *ver = new uint32_t[nodeNum]();

    std::unordered_map<uint32_t,bool> isInVerMap;
    bool inMemDo = false,writeIntoSSD = false,inPartition = true;
    uint32_t inside_iter = 0;

    offset *off = new offset[nodeNum];
    std::unordered_map<uint32_t,uint32_t> verToMap;
    std::unordered_map<uint32_t,uint32_t> verToMapNew;

    for(int i = 0; i < nodeNum; i++){
        off[i] = {edgeListBegPtr[i],edgeListBegPtr[i]+degree[i]*sizeof(uint32_t)};
        // printf("i: %u, start: %u, end: %u\n",i,off[i].start,off[i].end);
        verToMap[i] = i;
    }
    while(end < node_copy){
        if((off[end].end - pos)/sizeof(uint32_t) <= MEM){
            isInVerMap[end] = true;
            ver[index++] = end;
            end++;
            inMemDo = false;
        }
        else{
            // printf("inside_iter: %u\n",inside_iter);
            inMemDo = true,inPartition = true;
            inside_iter++;
            processUBInMem(pos,end,index,nbr,eid,ver,sup_arr,binCount,off,isInVerMap,fDat,fEid,fSupp,fUB);
        }
    }
    int current_pid = GetCurrentPid();
    float memory_usage = GetMemoryUsage(current_pid);
    memUsage = max(memory_usage,memUsage);
    if(!inMemDo){
        // printf("inside_iter: %u\n",inside_iter);
        processUBInMem(pos,end,index,nbr,eid,ver,sup_arr,binCount,off,isInVerMap,fDat,fEid,fSupp,fUB);
        inside_iter++;
    }

    delete[] nbr;
    delete[] sup_arr;
    delete[] eid;
    delete[] off;
    delete[] ver;
    delete[] binCount;
    fclose(fUB);
    total_io += fDat.get_total_io();
    total_io += fEid.get_total_io();
    total_io += fSupp.get_total_io();
    fDat.fclose();
    fEid.fclose();
    fSupp.fclose();
    log_debug(graphClock_.Count("topdown_k: %u, memUsage: %f",topdown_k,memUsage));
}


void Graph::triangleListExtendedGraph(readFile &file){
    MyReadFile fOff( file.m_offset );
    fOff.fopen( BUFFERED );

    #ifdef opt
    string name = "support_sort_edge";
    string sub_dir = file.m_base + name;
    file.createDir(sub_dir);
    uvSup *edges = new uvSup[file.memEdges];
    #else
    string name = "WC_triangle";
    string sub_dir = file.m_base + name;
    #endif
    file.createDir(sub_dir);
    char filename[100];
    char filename_eid[100];
    FILE* fSupp = fopen(file.m_supp.c_str(),"wb");

    uint32_t node_copy = nodeNum, tmp_node = 0, tmp_file = 0, size = 0, triangle = 0;
    uint64_t pos = 0, end = 0, index = 0, tmp_edge = edgeNum*2, last_tmp_edge = tmp_edge;
    uint32_t *nbr = new uint32_t[MEM]();
    uint32_t *nbr_global = new uint32_t[MEM]();

    uint64_t *eid = new uint64_t[MEM]();
    uint64_t *eid_global = new uint64_t[MEM]();

    uint32_t *ver = new uint32_t[nodeNum]();
    uint32_t *ver_global = new uint32_t[nodeNum]();
    uint32_t *deg_global = new uint32_t[nodeNum]();
    bool *isInVer = new bool[nodeNum](); 
    std::unordered_map<uint32_t,bool> isInVerMap;
    offset *off = new offset[nodeNum];
    std::unordered_map<uint32_t,uint32_t> verToMap;
    std::unordered_map<uint32_t,uint32_t> verToMapNew;
    std::unordered_map<uint64_t,uint32_t> eidInTri;
    for(int i = 0; i < nodeNum; i++){
        off[i] = {edgeListBegPtr[i],edgeListBegPtr[i]+degree[i]*sizeof(uint32_t)};
        // printf("i: %u, start: %u, end: %u\n",i,off[i].start,off[i].end);
        ver_global[i] = i;
        verToMap[i] = i;
    }
    uint32_t iter = 0, last_iter = uint32_t(-1);
    int sz = 0, tmpfile = 0;
    string filename_string;
    string filename_string_eid;
    bool inPartition = false; 

    char fileName[100];
    // ofstream out("mytext.txt");

    while(tmp_edge != 0){
        if(tmp_edge > MEM){
            verToMapNew.clear();
            log_info(graphClock_.Count("iter: %u, node_copy: %u, tmp_edge: %lu, MEM: %u",iter,node_copy,tmp_edge,MEM));
            if(iter == 0){
                sprintf(filename,"%s%s/edges_tmp_%d",file.m_base.c_str(),name.c_str(),iter);
                filename_string = file.m_dat;
                sprintf(filename_eid,"%s%s/edges_tmp_eid_%d",file.m_base.c_str(),name.c_str(),iter);
                filename_string_eid = file.m_eid;                
                last_iter = iter;
            }
            else{
                filename_string = file.m_base + name + "/edges_tmp_" + to_string(last_iter);
                sprintf(filename,"%s%s/edges_tmp_%d",file.m_base.c_str(),name.c_str(),iter);
                filename_string_eid = file.m_base + name + "/edges_tmp_eid_" + to_string(last_iter);
                sprintf(filename_eid,"%s%s/edges_tmp_eid_%d",file.m_base.c_str(),name.c_str(),iter);
                last_iter = iter;
            }
            MyReadFile fDat( filename_string );
            fDat.fopen( BUFFERED );
            FILE* fo = fopen(filename,"wb");

            MyReadFile fEid( filename_string_eid );
            fEid.fopen( BUFFERED );
            FILE* fo_eid = fopen(filename_eid,"wb");

            pos = 0, end = 0, index = 0;
            tmp_edge = 0, tmp_node = 0, size = 0;
            bool inMemDo = false,writeIntoSSD = false;
            uint32_t inside_iter = 0;
            while(end < node_copy){
                if((off[end].end - pos)/sizeof(uint32_t) <= MEM){
                    isInVer[ver_global[end]] = true;
                    isInVerMap[ver_global[end]] = true;
                    ver[index++] = ver_global[end];
                    end++;
                    inMemDo = false;
                }
                else{
                    inMemDo = true,inPartition = true;
                    inside_iter++;
                    processTriangleInMem(pos,end,index,tmp_node,tmp_edge,triangle,size,iter,nbr,nbr_global,eid,eid_global,deg_global,ver,ver_global,
                    off,verToMap,verToMapNew,isInVerMap,eidInTri,isInVer,fDat,fo,fEid,fo_eid,fOff,writeIntoSSD,false,fSupp);
                }
            }
            if(!inMemDo){
                inside_iter++;
                // log_info(graphClock_.Count("inMemDo is false, size: %u",size));
                processTriangleInMem(pos,end,index,tmp_node,tmp_edge,triangle,size,iter,nbr,nbr_global,eid,eid_global,deg_global,ver,ver_global,
                off,verToMap,verToMapNew,isInVerMap,eidInTri,isInVer,fDat,fo,fEid,fo_eid,fOff,writeIntoSSD,true,fSupp);
            }
            // log_info(graphClock_.Count("inside_iter: %u, triangle: %u, MEM: %u, tmp_edge: %lu",inside_iter,triangle,MEM,tmp_edge));
            if(tmp_edge == last_tmp_edge){  // this bug costs me a day long.
                MEM += MEM/10;
                delete[] nbr_global;
                delete[] nbr;
                nbr = new uint32_t[MEM]();
                nbr_global = new uint32_t[MEM]();
                delete[] eid_global;
                delete[] eid;
                eid = new uint64_t[MEM]();
                eid_global = new uint64_t[MEM]();                
                
            }
            fclose(fo);
            total_io += fDat.get_total_io();
            fDat.fclose();
            fclose(fo_eid);
            total_io += fEid.get_total_io();
            fEid.fclose();            
            uint64_t count = 0;
            for(int i = 0; i < tmp_node; i++)
            {
                off[i] = {count*sizeof(uint32_t),(count+deg_global[i])*sizeof(uint32_t)};
                count += deg_global[i];
            }
            node_copy = tmp_node;

            // verToMap.clear();
            // for(auto it = verToMapNew.begin(); it != verToMapNew.end(); it++)
            //     verToMap.insert({it->first,it->second});
            verToMap = verToMapNew;
            memset(deg_global,0,sizeof(uint32_t)*nodeNum);
            last_tmp_edge = tmp_edge;
            iter++;
        }
        else{
            // printf("avalible mem can process, triangle: %u, tmp_edge: %u\n",triangle,tmp_edge);
            if(inPartition){
                swap(nbr,nbr_global);
                swap(eid,eid_global);
            }
            else{
                MyReadFile fDat( file.m_dat );
                fDat.fopen( BUFFERED );
                loadInfo(nbr,edgeNum*2,pos,fDat);
                total_io += fDat.get_total_io();
                fDat.fclose();

                MyReadFile fEid( file.m_eid );
                fEid.fopen( BUFFERED );
                loadInfo(eid,edgeNum*2,pos*2,fEid);
                total_io += fEid.get_total_io();
                fEid.fclose();
                tmp_node = nodeNum;
            }

            for(int i = 0; i < tmp_node; i++){
                uint32_t u = ver_global[i];
                for(int j = off[verToMap[u]].start/sizeof(uint32_t); j < off[verToMap[u]].end/sizeof(uint32_t); j++){
                    uint32_t v = nbr[j];
                    uint64_t eid_uv = eid[j];
                    // printf("u: %u, v: %u, u_start: %u, u_end: %u\n",u,v,off[verToMap[u]].start,off[verToMap[u]].end);
                    assert(eid_uv < edgeNum);
                    if(v < u) continue;
                    uint32_t sup = 0,sup_copy = 0;
                    uint64_t ptr_u = off[verToMap[u]].start/sizeof(uint32_t), upper_u = off[verToMap[u]].end/sizeof(uint32_t);
                    // uint64_t eid_ptr_u = off[verToMap[u]].start*2/sizeof(uint64_t), eid_upper_u = off[verToMap[u]].end*2/sizeof(uint64_t);
                    uint64_t ptr_v = off[verToMap[v]].start/sizeof(uint32_t), upper_v = off[verToMap[v]].end/sizeof(uint32_t);
                    // uint64_t eid_ptr_v = off[verToMap[v]].start*2/sizeof(uint64_t), eid_upper_v = off[verToMap[v]].end*2/sizeof(uint64_t);
                    while(ptr_u < upper_u && ptr_v < upper_v){
                        if(nbr[ptr_u] < nbr[ptr_v]) ptr_u++;
                        else if(nbr[ptr_u] > nbr[ptr_v]) ptr_v++;
                        else{
                            if(nbr[ptr_u] > v && nbr[ptr_v] > v)
                                sup++;
                            ptr_u++;
                            ptr_v++;
                            sup_copy++;
                        }
                    }
                    if(eidInTri.find(eid_uv) != eidInTri.end()){
                        sup_copy += eidInTri[eid_uv];
                        eidInTri.erase(eid_uv);
                    }
                    maxSup = std::max(sup_copy,maxSup);
                    prefix[sup_copy]++;
                    #ifdef opt
                    edges[sz].u = u;
                    edges[sz].v = v;
                    edges[sz].sup = sup_copy; // 这里可以优化，减少写入?
                    sz++;                
                    if(sz >= file.memEdges){            
                        sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpfile);
                        file.saveTmpEdges<uvSup>(edges,sz,fileName,[](const uvSup & a, const uvSup & b) {
                            return a.sup < b.sup;
                            });
                        sz = 0;
                        ++tmpfile;
                    }
                    //count range of edges' support
                    
                    #else
                    eid_eid tmp;
                    fOff.fseek(eid_uv*sizeof(eid_eid));
                    fOff.fread(&tmp,sizeof(eid_eid));
                    fseek(fSupp,tmp.first*sizeof(uint32_t),SEEK_SET);
                    fwrite(&sup_copy,sizeof(uint32_t),1,fSupp);
                    fseek(fSupp,tmp.second*sizeof(uint32_t),SEEK_SET);
                    fwrite(&sup_copy,sizeof(uint32_t),1,fSupp);
                    #endif

                    total_io += 2;
                    triangle += sup;
                }
            }
            // out << "maxSup: "<< maxSup << endl;
            // for(int i = 0; i <= maxSup; i++)
            //     if(prefix[i] != 0)
            //         out << "sup: "<< i << ",number: " << prefix[i] << endl;
            // out <<"---------------------------"<< endl;
            // for(int i = 1; i <= maxSup; i++)
            //     prefix[i]+= prefix[i-1];
            // for(int i = 1; i <= maxSup; i++)
            //     if(prefix[i] != prefix[i-1])
            //         out << "sup: "<< i << ",total: " << prefix[i] << endl; 
            // out.close();
            #ifdef opt
            sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpfile);
            file.saveTmpEdges<uvSup>(edges,sz,fileName,[](const uvSup & a, const uvSup & b) {
                return a.sup < b.sup;
            });
            file.mergeBySup<uvSup>(tmpfile+1);
            // print support of all edges
            
            #endif
            break;
        }
    }
    int current_pid = GetCurrentPid();
    float memory_usage = GetMemoryUsage(current_pid);
    memUsage = max(memory_usage,memUsage);
    log_debug(graphClock_.Count("Triangles: %u, memUsage: %f",triangle,memUsage));

    if(nbr != NULL)
        delete[] nbr;
    if(nbr_global != NULL)
        delete[] nbr_global;
    if(isInVer != NULL)
        delete[] isInVer;
    if(ver != NULL)
        delete[] ver;
    if(ver_global != NULL)
        delete[] ver_global;
    if(deg_global != NULL)
        delete[] deg_global;
    delete[] eid;
    delete[] eid_global;    
    fclose(fSupp);
    total_io += fOff.get_total_io();
    fOff.fclose();
}

void Graph::writeTwoVerFromEid(readFile &file){
    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fEid( file.m_eid );
	fEid.fopen( BUFFERED );
    //write the two vertex of edge into fEidToV 
    FILE* fEidToV = fopen(file.m_eidToVer.c_str(),"wb");
    uint32_t sup,u_nbrNum;
    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    #ifdef DegSort
    for(uint32_t i = 0; i < nodeNum; i++){
        u_nbrNum = degree[i]-(edgeListBegPtrPlus[i]-edgeListBegPtr[i])/sizeof(uint32_t);
        if(u_nbrNum == 0)continue;
        loadInfo(nbr_u,u_nbrNum,edgeListBegPtrPlus[i],fDat);
        uint32_t u = i;
        for(uint32_t j = 0; j < u_nbrNum; j++){
            uint32_t v = nbr_u[j];
            uint64_t offset = edgeListBegPtrPlus[i]*2+j*sizeof(uint64_t);
    #else
    for(uint32_t i = 0; i < nodeNum; i++){
        loadInfo(nbr_u,degree[i],edgeListBegPtr[i],fDat);
        uint32_t u = i;
        
        for(uint32_t j = 0; j < degree[i]; j++){
            uint32_t v = nbr_u[j];
            uint64_t offset = edgeListBegPtr[i]*2+j*sizeof(uint64_t);
    #endif
            uint64_t eid_uv;
            fEid.fseek(offset);
            fEid.fread(&eid_uv,sizeof(uint64_t));
            Edge tmp = {u,v};
            fseek(fEidToV,eid_uv*sizeof(Edge),SEEK_SET);
            fwrite(&tmp,sizeof(Edge),1,fEidToV);
            total_io += 1;
        }
    }

    fclose(fEidToV);
    total_io += fDat.get_total_io();
    fDat.fclose();
    total_io += fEid.get_total_io();
    fEid.fclose();
    free(nbr_u);
    log_info(graphClock_.Count("finish writing eid to vertex"));
}
void Graph::updateEdgeDram(uint32_t u, uint32_t v, uint64_t eid, uint32_t minsup,std::unordered_map<uint32_t,uint32_t> &verToMap)
{
    if (u > v) std::swap(u,v);
    uint32_t sup=supports[eid];
	if (sup<=minsup) return;
	int p=pos[verToMap[u]][verToMap[v]];
	int posbin=bin[sup];
	TEdge se=binEdge[posbin];
	TEdge e={u,v,eid};
	if (p!=posbin) {
		pos[verToMap[u]][verToMap[v]]=posbin;
		pos[verToMap[se.u]][verToMap[se.v]]=p;
		binEdge[p]=se;
		binEdge[posbin]=e;
	}
	++bin[sup];
    supports[eid] -= 1;
}
void Graph::extractSubgraph(readFile &file)
{
    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fEid( file.m_eid );
	fEid.fopen( BUFFERED );
    MyReadFile fEidToV( file.m_eidToVer );
	fEidToV.fopen( BUFFERED );
    MyReadFile fUB( file.m_suppSort );
	fUB.fopen( BUFFERED );
    std::set<uint32_t> vertexSet;
    // bool *vertexArr = new bool[nodeNum]();
    std::unordered_map<uint64_t,bool> delEdge;
    uint32_t e_ub,u,v; //upper bound of each edge
    es e_eidsup;

    #ifdef opt
    uvSup e;
    topdown_k = 2207;
    #endif;

    bool find = false;
    uint32_t successFind = 0, Truss = 0, TrussEdge = 0;

while(!find && topdown_k > 0)
{
    vertexSet.clear();
    // printf("prefix[topdown_k-1]: %u\n",prefix[topdown_k-1]/sizeof(es));

    #ifdef opt
    fUB.fseek(prefix[topdown_k-1] * sizeof(uvSup));
    #else
    fUB.fseek(prefix[topdown_k-1] * sizeof(es));
    #endif;

    
    if((edgeNum - prefix[topdown_k-1]) < (topdown_k+1)*(topdown_k+2) / 2)
    {
        topdown_k--;
        continue;
    }
    for(uint64_t i = prefix[topdown_k-1]; i < edgeNum; i++){
        #ifdef opt
        fUB.fread(&e,sizeof(uvSup));
        #else
        fUB.fread(&e_eidsup,sizeof(es));
        uint64_t eid = e_eidsup.eid;
        Edge e;
        fEidToV.fseek(eid*sizeof(Edge));
        fEidToV.fread(&e,sizeof(Edge));
        #endif;

        u = e.u;         
        if(vertexSet.find(u) == vertexSet.end())
            vertexSet.insert(u);
        v = e.v;
        if(vertexSet.find(v) == vertexSet.end())
            vertexSet.insert(v);
        
    }
    uint32_t tmp_node = vertexSet.size(), index = 0;
    uint64_t tmp_edge = 0;
    std::unordered_map<uint32_t,uint32_t> verToMap;
    std::vector<uint32_t> ver;
    std::vector<uint64_t> off;
    off.push_back(0);
    // uint32_t *ver = new uint32_t[tmp_node]();
    // uint64_t *off = new uint64_t[tmp_node+1]();
    
    for(auto it = vertexSet.begin(); it != vertexSet.end(); it++){
        verToMap[*it] = index;
        ver.push_back(*it);
        off.push_back(off.back() + degree[*it]);
        // ver[index++] = *it;
        // off[index] = off[index-1] + degree[*it];
        // deg_new[*it] = degree[*it];
        index++;
    }
    // for(int i = 0; i < nodeNum; i++){
    //     if(!vertexArr[i]) continue;
    //     verToMap[i] = index;
    //     ver.push_back(i);
    //     off.push_back(off.back() + degree[i]);
    //     index++;
    // }

    std::vector<uint32_t> nbr;
    std::vector<uint32_t> eid;


    uint32_t *nbr_copy = new uint32_t[file.maxDeg]();
    uint64_t *eid_copy = new uint64_t[file.maxDeg]();

    std::unordered_map<uint32_t,std::vector<vertex_eid>> reverseEdge;

    for(int i = 0; i < tmp_node; i++){
        u = ver[i];
        // printf("u: %u, deg: %u\n",u,degree[u]);
        loadInfo(nbr_copy,degree[u],edgeListBegPtr[u],fDat);
        loadInfo(eid_copy,degree[u],edgeListBegPtr[u]*2,fEid);
        
        for(int j = 0; j < degree[u]; j++){
            if(vertexSet.find(nbr_copy[j]) == vertexSet.end()){
                if(verToMap.find(nbr_copy[j]) == verToMap.end()){
                    verToMap[nbr_copy[j]] = index;
                    ver.push_back(nbr_copy[j]);
                    index++;
                }
                reverseEdge[nbr_copy[j]].push_back({u,eid_copy[j]});
            }
            nbr.push_back(nbr_copy[j]);
            eid.push_back(eid_copy[j]);
        }
    }
    for(int i = tmp_node; i < index; i++){
        u = ver[i];
        for(int j = 0; j < reverseEdge[u].size(); j++){
            nbr.push_back(reverseEdge[u][j].v);
            eid.push_back(reverseEdge[u][j].eid);
        }
        off.push_back(off.back() + reverseEdge[u].size());
    }

    supports.clear();
    delEdge.clear();
    printf("ver size: %u, edge: %u, index: %u, topDown_k: %u\n",ver.size(),off[off.size()-1]/2,index,topdown_k);
    tmp_node = index;

    memset(bin,0,sizeof(uint32_t)*nodeNum);
    uint32_t mp = 0, nBins = 0;
    for(int i = 0; i < index; i++){
        u = ver[i];
        for(uint64_t j = off[verToMap[u]]; j < off[verToMap[u]+1]; j++){
            v = nbr[j];
            uint32_t sup = 0;
            uint64_t eid_uv = eid[j];            
            if(u > v) continue;

            uint64_t ptr_u = off[verToMap[u]], upper_u = off[verToMap[u]+1];
            uint64_t ptr_v = off[verToMap[v]], upper_v = off[verToMap[v]+1];
            while(ptr_u < upper_u && ptr_v < upper_v){
                if(nbr[ptr_u] < nbr[ptr_v])ptr_u++;
                else if(nbr[ptr_u] > nbr[ptr_v]) ptr_v++;
                else{
                    sup++;
                    ptr_u++;
                    ptr_v++;
                }
            }
            if(sup == 0){
                printClass(u,v,2);
                delEdge[eid_uv] = true;
                continue;
            }

            // printf("u: %u, v: %u,eid: %u, sup: %u\n",u,v,eid_uv,sup);
            nBins = std::max(sup,nBins);
            supports[eid_uv] = sup;
            bin[sup]++;
            mp++;
        }
    }
    nBins++;
	int count=0,s;
	for (int i=0; i<nBins; ++i) {
		int binsize=bin[i];
		bin[i]=count;
		count+=binsize;
        // printf("i: %u, bin[i]: %u\n",i,bin[i]);
	}   
	pos.clear(); pos.resize(tmp_node);
	for (int i=0; i<tmp_node; ++i) pos[i].clear();
    
    binEdge = new TEdge[mp];

    for(uint32_t i = 0; i < tmp_node; i++){
        uint32_t u = ver[i];
        for(uint32_t j = off[verToMap[u]]; j < off[verToMap[u]+1]; j++){
            uint64_t eid_uv = eid[j];
            if(delEdge.find(eid_uv) != delEdge.end())continue;
            uint32_t v = nbr[j];   
            if(u > v) continue;         
            uint32_t sup = supports[eid_uv];
            // printf("u: %u, v: %u,eid: %u, sup: %u, bin[sup]: %u\n",u,v,eid_uv,sup,bin[sup]);
            TEdge e = {u,v,eid_uv};
            binEdge[bin[sup]] = e;
            pos[verToMap[u]][verToMap[v]] = bin[sup]++;
        }
    }
    for(int i = nBins; i > 0; i--) bin[i] = bin[i-1];
    bin[0] = 0;

    uint32_t last_sup = uint32_t(-1);
    for(s = 0; s < mp; s++){
        uint32_t u = binEdge[s].u;
		uint32_t v = binEdge[s].v;
        uint64_t eid_uv = binEdge[s].eid; 
		uint32_t sup = supports[eid_uv];
        // printf("u: %u, v: %u, eid: %u, sup: %u, pos: %u\n",u,v,eid_uv,sup,pos[verToMap[u]][verToMap[v]]);
		printClass(u,v,sup+2);

        if(sup >= topdown_k){
            int current_pid = GetCurrentPid();
            float memory_usage = GetMemoryUsage(current_pid);
            memUsage = max(memory_usage,memUsage);
            if(successFind)
                find = true;
            successFind += 1;
            topdown_k += 3;
            Truss = sup+2, TrussEdge = mp-s;
            break;
        }
        /* optimization */
        std::vector<ver_two_eid> comm;
        uint64_t ptr_u = off[verToMap[u]], upper_u = off[verToMap[u]+1];
        uint64_t ptr_v = off[verToMap[v]], upper_v = off[verToMap[v]+1];
        while(ptr_u < upper_u && ptr_v < upper_v){
            if(nbr[ptr_u] < nbr[ptr_v] || delEdge.find(eid[ptr_u]) != delEdge.end())ptr_u++;
            else if(nbr[ptr_u] > nbr[ptr_v] || delEdge.find(eid[ptr_v]) != delEdge.end()) ptr_v++;
            else{
                comm.push_back({nbr[ptr_u],eid[ptr_u],eid[ptr_v]});
                ptr_u++;
                ptr_v++;
            }
        }
        for(int i = 0; i < comm.size(); i++){
            updateEdgeDram(u,comm[i].w,comm[i].first,sup,verToMap);
			updateEdgeDram(v,comm[i].w,comm[i].second,sup,verToMap);
        }
        delEdge[eid_uv] = true;
    }
    delete[] binEdge;
    delete[] nbr_copy;

    if(successFind && s == mp)
        find = true;

    topdown_k -= 2;
    // fout << "############################"<< std::endl;
    // for (int i=0;i<nodeNum; ++i)
	// 	if (cntClass[i]>0)
	// 		fout << "#edges in " << i << "-class: " << cntClass[i] << std::endl;
    memset(cntClass,0,sizeof(int)*nodeNum);
}
    printf("max truss: %u, edge: %u, io: %lu, memUsage: %f\n",Truss,TrussEdge,total_io,memUsage);

    total_io = total_io + fDat.get_total_io() + fUB.get_total_io() + fEidToV.get_total_io() + fEid.get_total_io();
    fEidToV.fclose();
    fUB.fclose();
    fDat.fclose();
    fEid.fclose();
}

void Graph::sortUpperBound(readFile &file){
    es* arr = new es[file.memEdges];
    MyReadFile fUB( file.m_upperBound );
	fUB.fopen( BUFFERED );
    uint32_t e_ub;
    uint32_t size = 0,tmpfile = 0;
    string sub_dir = file.m_base + "support_sort_edge";
    char fileName[100];
    file.createDir(sub_dir);
    memset(prefix,0,nodeNum);
    for(uint64_t i = 0; i < edgeNum; i++){
        fUB.fread(&e_ub,sizeof(uint32_t));
        arr[size] = {i,e_ub};
        size++;                
        if(size >= file.memEdges){            
            sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpfile);
            file.saveTmpEdges<es>(arr,size,fileName,[](const es & a, const es & b) {
                return a.sup < b.sup;
                });
            size = 0;
            ++tmpfile;
        }
        // printf("u: %u, v: %u, eid: %lu, sup: %u\n",u,v,eid_uv,sup_uv);
        prefix[e_ub]++;
    }
    for(int i = 1; i <= topdown_k; i++){
        prefix[i] += prefix[i-1];
        // printf("i: %d, pos: %u\n",i,prefix[i]);
    }
    sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpfile);
    file.saveTmpEdges<es>(arr,size,fileName,[](const es & a, const es & b) {
        return a.sup < b.sup;
    });
    file.mergeBySup<es>(tmpfile+1);
    total_io = total_io + file.write_io + fUB.get_total_io();
    fUB.fclose();
    delete[] arr;
}

void Graph::topDownTrussDecom(readFile &file){

    triangleListExtendedGraph(file);
    log_info(graphClock_.Count("triangleListExtendedGraph done"));
    upperBounding(file);
    log_info(graphClock_.Count("upperBounding done"));
    
    sortUpperBound(file);
    log_info(graphClock_.Count("sortUpperBound done"));
    writeTwoVerFromEid(file);
    log_info(graphClock_.Count("writeTwoVerFromEid done"));
    extractSubgraph(file);
}
