#include "include/graph.h"
#include "include/Time.h"
#include "include/zip_iterator.h"
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

void Graph::binSortSSD(readFile &file)
{
    MyReadFile fSup( file.m_supp );
	fSup.fopen( NOBUFFER );
    MyReadFile fEid( file.m_eid );
	fEid.fopen( BUFFERED );
    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fOff( file.m_offset );
	fOff.fopen( BUFFERED );

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
            uint32_t sup,sup_copy;
            fEid.fseek(offset);
            fEid.fread(&eid_uv,sizeof(uint64_t));
            if(MOVE(eid_uv)) //this edge is deleted, the most left bit is 1
                continue;
            // if(!(degree[v]<degree[u] || (degree[u]==degree[v] && v<u)))
            //     continue;
            if(u > v) continue;
            eid_eid tmp;
            fOff.fseek(eid_uv*sizeof(eid_eid));
            fOff.fread(&tmp,sizeof(eid_eid));

            fSup.fseek(tmp.first*sizeof(uint32_t));
            fSup.fread(&sup,sizeof(uint32_t));
            // fSup.fseek(tmp.second*sizeof(uint32_t));
            // fSup.fread(&sup_copy,sizeof(uint32_t));
            
            if(sup == 0){   
                printClass(u,v,2);
                sup = DELETE(sup);  
                fSup.fseek(tmp.first*sizeof(uint32_t));
                fSup.fwrite(&sup,sizeof(uint32_t));
                fSup.fseek(tmp.second*sizeof(uint32_t));
                fSup.fwrite(&sup,sizeof(uint32_t));
                continue;
            }
            ++mp;
			++bin[sup];
			nBins=std::max(sup,nBins); 
        }
    }
    int m=mp;
	++nBins;
	int count=0;
	for (int i=0; i<nBins; ++i) {
		int binsize=bin[i];
		bin[i]=count;
        // if(bin[i] != 0){
        //     printf("i: %d, bin[i]: %d\n",i,bin[i]);
        // }
		count+=binsize;
	}

    FILE* fBinEdge = fopen(file.m_binEdge.c_str(),"wb");
    FILE* fEdgePos = fopen(file.m_ePos.c_str(),"wb");


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
            uint32_t sup;
            fEid.fseek(offset);
            fEid.fread(&eid_uv,sizeof(uint64_t));
            
            // if(!(degree[v]<degree[u] || (degree[u]==degree[v] && v<u)))
            //     continue;
            if(u>v || MOVE(eid_uv))continue;

            eid_eid tmp;
            fOff.fseek(eid_uv*sizeof(eid_eid));
            fOff.fread(&tmp,sizeof(eid_eid));

            fSup.fseek(tmp.first*sizeof(uint32_t));
            fSup.fread(&sup,sizeof(uint32_t));
            if(MOVE(sup)) //this edge is deleted, the most left bit is 1
                continue;
            TEdge e = {u,v,eid_uv};
            int pos = bin[sup];

            fseek(fBinEdge,pos*sizeof(TEdge),SEEK_SET);
            fwrite(&e,sizeof(TEdge),1,fBinEdge);

            fseek(fEdgePos,eid_uv*sizeof(int),SEEK_SET);
            fwrite(&pos,sizeof(int),1,fEdgePos);

            bin[sup]++;
        }
    }
    for(int i = nBins; i > 0; i--) bin[i] = bin[i-1];
    bin[0] = 0;


    fSup.fclose();
    fEid.fclose();
    fDat.fclose();
    fOff.fclose();
    free(nbr_u);
    fclose(fBinEdge);
    fclose(fEdgePos);
}

void Graph::updateEdgeSSDNew(uint32_t u, uint32_t v, uint64_t eid, uint32_t ptr, uint32_t minsup, 
uint32_t *sup_arr, MyReadFile &fEdgePos, MyReadFile &fSup, MyReadFile &fBinEdge, MyReadFile &fOff){
    if (!(degree[v]<degree[u] || (degree[u]==degree[v] && v<u))) std::swap(u,v);
    TEdge se;
    uint32_t sup = sup_arr[ptr];
    int p;
	if (sup<=minsup) return;

    fEdgePos.fseek(eid*sizeof(int));
    fEdgePos.fread(&p,sizeof(int));
	int posbin=bin[sup];
    fBinEdge.fseek(posbin*sizeof(TEdge));
    fBinEdge.fread(&se,sizeof(TEdge));
	TEdge e={u,v,eid};
	if (p!=posbin) {
        fEdgePos.fseek(eid*sizeof(int));
        fEdgePos.fwrite(&posbin,sizeof(int));

        fEdgePos.fseek(se.eid*sizeof(int));
        fEdgePos.fwrite(&p,sizeof(int));

        fBinEdge.fseek(p*sizeof(TEdge));
        fBinEdge.fwrite(&se,sizeof(TEdge));

        fBinEdge.fseek(posbin*sizeof(TEdge));
        fBinEdge.fwrite(&e,sizeof(TEdge));
	}
	++bin[sup];
    sup--;
    eid_eid tmp;

    fOff.fseek(eid*sizeof(eid_eid));
    fOff.fread(&tmp,sizeof(eid_eid));

    fSup.fseek(tmp.second*sizeof(uint32_t));
    fSup.fwrite(&sup,sizeof(uint32_t));       

    fSup.fseek(tmp.first*sizeof(uint32_t));
    fSup.fwrite(&sup,sizeof(uint32_t));
    sup_arr[ptr] = sup;
}

void Graph::updateNbrInDeleteEdgeLazyUpdate(uint32_t u, uint32_t *nbr_u, uint64_t *eid_u,  uint32_t *sup_u,
uint32_t v, uint32_t *nbr_v, uint64_t *eid_v, uint32_t *sup_v,
MyReadFile &fDat, MyReadFile &fEid, MyReadFile &fSup,MyReadFile &fOff, MyReadFile &fPres, MyReadFile &fNexts,
uint32_t sup, queue<EdgeSup> &del_que, uint32_t mid, DynamicHeap &dheap, ListLinearHeapTruss *linear_heap)
{
    eid_eid tmp_;
    std::vector<eid_eid> comm;
    IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,sup,comm,true);
    for(int i = 0; i < comm.size(); i++){
        ui first_sup = sup_u[comm[i].first];
        ui second_sup = sup_v[comm[i].second];
        uint64_t eid_fir = eid_u[comm[i].first];
        uint64_t eid_sec = eid_v[comm[i].second];

        if(dheap.find(eid_fir)){
            if(dheap.getSup(eid_fir) <= mid){
                // if(isInDelQue.find(eid_fir) == isInDelQue.end()){
                //     del_que.push({u,nbr_u[comm[i].first],first_sup,eid_fir});
                //     isInDelQue.insert(eid_fir);
                // }                
            }
            else
                dheap.adjust(eid_fir,1);
        }
        else{
            if(first_sup > mid)
                dheap.push({eid_fir,first_sup - 1});
            else if(first_sup <= mid)
            {
                // if(isInDelQue.find(eid_fir) == isInDelQue.end())
                // {
                //     isInDelQue.insert(eid_fir);
                //     del_que.push({u,nbr_u[comm[i].first],first_sup,eid_fir});
                // }
            }    
        }

        if(dheap.find(eid_sec)){
            if(dheap.getSup(eid_sec) <= mid){

                // if(isInDelQue.find(eid_sec) == isInDelQue.end())
                // {
                //     isInDelQue.insert(eid_sec);
                //     del_que.push({v,nbr_v[comm[i].second],second_sup,eid_sec});
                // }
            }
            else
                dheap.adjust(eid_sec,1);
        }
        else{
            if(second_sup > mid)
                dheap.push({eid_sec,second_sup - 1});
            else if(second_sup <= mid)
            {
                // if(isInDelQue.find(eid_sec) == isInDelQue.end()){
                //     isInDelQue.insert(eid_sec);
                //     del_que.push({v,nbr_v[comm[i].second],second_sup,eid_sec});
                // }

            }    
        }   
    }
}

void Graph::IntersectTrussNew(uint32_t u, uint32_t *nbr_u, uint64_t *eid_u,  uint32_t *sup_u,
uint32_t v, uint32_t *nbr_v, uint64_t *eid_v, uint32_t *sup_v, MyReadFile &fDat, MyReadFile &fEid, MyReadFile &fSup,
uint32_t sup, vector<eid_eid>& comm, bool flag){
    uint32_t ptr_u = 0;
    uint32_t upper_u = degree[u];
    uint32_t ptr_v = 0;
    uint32_t upper_v = degree[v];
    uint64_t u_pos = 0,v_pos = 0;

    if((upper_u > 0 && upper_v > 0) && (nbr_u[upper_u-1] < nbr_v[0] || nbr_v[upper_v-1] < nbr_u[0]))
        return ;

    uint32_t count = 0;
    while(ptr_u < upper_u && ptr_v < upper_v){
        
        if(flag && count == sup) break;
        if(vis[nbr_u[ptr_u]] || flag && MOVE(sup_u[ptr_u]))
        {
            ptr_u++;
            continue;
        }

        if(vis[nbr_v[ptr_v]] || flag && MOVE(sup_v[ptr_v]))
        {
            ptr_v++;
            continue;
        }
        if(nbr_u[ptr_u] < nbr_v[ptr_v])
            ptr_u++;
        else if(nbr_u[ptr_u] > nbr_v[ptr_v])
            ptr_v++;
        else
        {
            // printf("ptr_u: %u, ptr_v: %u\n",ptr_u,ptr_v);
            comm.push_back({ptr_u,ptr_v});
            ptr_u++;
            ptr_v++;
            count++;
        }
    }
}


void Graph::Initial(readFile &file){
    MyReadFile fIdx( file.m_idx );
	fIdx.fopen( BUFFERED );
	MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    uint64_t tmpa = 0,tmpb = 0;
    uint32_t degreeTmp = 0, max_deg = 0;
    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    for (uint32_t i = 0; i < nodeNum; ++i){
		fIdx.fread(&tmpa,sizeof(uint64_t));
        edgeListBegPtr[i] = tmpa;
        #ifdef DegSort
        fIdx.fread(&tmpb,sizeof(uint64_t));
        edgeListBegPtrPlus[i] = tmpb;
        #endif
		fIdx.fread(&degreeTmp,sizeof(uint32_t));
        max_deg = std::max(max_deg,degreeTmp);
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
    file.maxDeg = max_deg;
    total_io += fIdx.get_total_io();
    free(nbr_u);
    fDat.fclose();
	fIdx.fclose();
}

void Graph::InitialUnOrderDegSort(readFile &file){
    MyReadFile fIdx( file.m_idx );
    // cout<<"InitialUnOrderDegSort:"<<file.m_idx<<endl;

	fIdx.fopen( BUFFERED );
    uint64_t tmp = 0,tmpb;
    uint32_t degreeTmp = 0;
    for (uint32_t u = 0; u < file.verNum; ++u){
        int i = file.vertexId[u];
		fIdx.fread(&tmp,sizeof(uint64_t));
        #ifdef DegSort
        fIdx.fread(&tmpb,sizeof(uint64_t));
        edgeListBegPtrPlus[i] = tmpb;
        #endif
		fIdx.fread(&degreeTmp,sizeof(uint32_t));
        degree[i] = degreeTmp;
        degree_[i] = degreeTmp;
        edgeListBegPtr[i] = tmp;
        // printf("ver: %d, deg: %d, pos: %ld\n",i, degreeTmp, tmp);
  	}
	fIdx.fclose();
    total_io += fIdx.get_total_io();
}
void Graph::InitialUnOrder(readFile &file){
    MyReadFile fIdx( file.m_idx );
	fIdx.fopen( BUFFERED );
	MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    uint64_t tmp = 0,tmpb;
    uint32_t degreeTmp = 0;
    for (uint32_t u = 0; u < file.verNum; ++u){
        int i = file.vertexId[u];
		fIdx.fread(&tmp,sizeof(uint64_t));
		fIdx.fread(&degreeTmp,sizeof(uint32_t));
        degree[i] = degreeTmp;
        degree_[i] = degreeTmp;
        edgeListBegPtr[i] = tmp;
        // printf("ver: %d, deg: %d, pos: %ld\n",i, degreeTmp, tmp);
  	}
    total_io += 2*file.verNum;
    fDat.fclose();
	fIdx.fclose();
}


template<typename T>
void Graph::loadInfo(T* nbr, uint32_t num, uint64_t pos, MyReadFile& fDat){
    fDat.fseek(pos);
    fDat.fread(nbr,sizeof(T)*num);
}


int Graph::IntersectTriangleSSDByDegOrder(uint32_t u, uint32_t u_nbrNum, uint32_t *nbr_u, uint32_t v, uint32_t v_nbrNum, uint32_t *nbr_v, 
MyReadFile &fDat, MyReadFile &fEid, std::vector<eid_eid> &common, std::unordered_map<uint64_t,uint32_t> &map_pos)
{
    int ptr_u = 0;
    int upper_u = u_nbrNum;
    int ptr_v = 0;
    int upper_v = v_nbrNum;

    if(v_nbrNum == 0 || nbr_u[upper_u-1] < nbr_v[0] || nbr_v[upper_v-1] < nbr_u[0])
        return 0;

    int count = 0;
    while(ptr_u < upper_u && ptr_v < upper_v){
        if(vis[nbr_u[ptr_u]] || nbr_u[ptr_u] < nbr_v[ptr_v])
            ptr_u++;
        else if(vis[nbr_v[ptr_v]] || nbr_u[ptr_u] > nbr_v[ptr_v])
            ptr_v++;
        else
        {
            uint64_t upos = edgeListBegPtrPlus[u]*2+ptr_u*sizeof(uint64_t);
            map_pos[upos]++;
            uint64_t vpos = edgeListBegPtrPlus[v]*2+ptr_v*sizeof(uint64_t);
            map_pos[vpos]++;
            uint64_t eid_uw,eid_vw;            
            ptr_u++;
            ptr_v++;
            count++;
        }
    }
    return count;
}


void Graph::KCore(readFile &file, uint32_t *coreNum, bool isDynamicload){
    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    uint32_t maxDeg = file.maxDeg, d = 0;
    int *bin = new int[maxDeg+2]();
    
    int debug_c = 0;
    for (int i = 0; i < nodeNum; i++){
        if(isDynamicload)
        {
            bin[degree_[i]]++;
            // if(degree_[i] != degree[i]) debug_c++;
        }
            
        else
            bin[degree[i]]++;
    }
    // printf("miss : %d\n",debug_c);
        
    int lastBin = bin[0], nowBin;
    bin[0] = 0;
    for (int i = 1; i <= maxDeg; i++)
    {
        nowBin = lastBin + bin[i - 1];
        lastBin = bin[i];
        bin[i] = nowBin;
    }

    uint32_t *vert = new uint32_t[nodeNum](), *pos = new uint32_t[nodeNum](), *tmpDeg = new uint32_t[nodeNum]();
    for (uint32_t i = 0; i < nodeNum; i++)
    {
        if(isDynamicload) 
            d = degree_[i];
        else 
            d = degree[i];
        pos[i] = bin[d];
        vert[bin[d]++] = i;
        tmpDeg[i] = d;
    }

    
    for (int i = maxDeg; i >= 1; i--)
        bin[i] = bin[i - 1];
    bin[0] = 0;

    uint32_t maxCore_ = 0, nbr, deg_u = 0;
    uint32_t *nbr_u = new uint32_t[maxDeg];
    for(int i = 0; i < nodeNum; i++){
        uint32_t id = vert[i], nbr, binFrontId;
        coreNum[id] = tmpDeg[id];
        // degeneracyOrder[i] = id;
        maxCore_ = std::max(maxCore_,coreNum[id]);  
        deg_u = degree[id];
        if(isDynamicload) 
            loadNbrForDynamic(id,nbr_u,deg_u,edgeListBegPtr[id],fDat);
        else
            loadInfo(nbr_u,degree[id],edgeListBegPtr[id],fDat);
        for(int j = 0; j < deg_u; j++){
            nbr = nbr_u[j];
            if (tmpDeg[nbr] > tmpDeg[id])
            {
                binFrontId = vert[bin[tmpDeg[nbr]]];
                if (binFrontId != nbr)
                {
                    pos[binFrontId] = pos[nbr];
                    pos[nbr] = bin[tmpDeg[nbr]];
                    vert[bin[tmpDeg[nbr]]] = nbr;
                    vert[pos[binFrontId]] = binFrontId;
                }
                bin[tmpDeg[nbr]]++;
                tmpDeg[nbr]--;
            }
        }
    }
    maxCore = maxCore_;
    log_debug(graphClock_.Count("maxCore: %u",maxCore_));
    // string ss = file.m_base + "graph.core";
    // FILE* f = fopen(ss.c_str(),"wb");
    // for(int i = 0; i < nodeNum; i++)
    //     fwrite(&coreNum[i],sizeof(uint32_t),1,f);
    // fclose(f);

    delete[] bin;
    delete[] tmpDeg;
    delete[] pos;
    delete[] vert;
    delete[] nbr_u;
    total_io += fDat.get_total_io();
    fDat.fclose();
}


void Graph::UpdateCoreDynamic(readFile &file, uint32_t *coreNum, uint32_t u, uint32_t v){
    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    bool *vis = new bool[nodeNum]();
    uint32_t *nbr_u = new uint32_t[nodeNum]();
    uint32_t *deg_copy = new uint32_t[nodeNum]();
    memcpy(deg_copy, degree_, sizeof(uint32_t)*nodeNum);
    queue<uint32_t> ver; //affected vertex
    uint32_t coreness = 0;
    if(coreNum[u] == coreNum[v]){  // update coreness of vertexs related with u and v
        ver.push(u);
        coreness = coreNum[u];
    }
    else if(coreNum[u] > coreNum[v]){  // update coreness of vertexs related with v
        ver.push(v);
        coreness = coreNum[v];
    }
    else{ // update coreness of vertexs related with u
        ver.push(u);
        coreness = coreNum[u];
    }

    struct cmp{
        bool operator()(VerCore a, VerCore b){
            return a.deg > b.deg;
        }
    };
    std::priority_queue<VerCore,std::vector<VerCore>, cmp> pq; 
    unordered_set<uint32_t> isInVer;
    while(!ver.empty()){
        uint32_t x = ver.front(), u_degree = 0;
        vis[x] = true, u_degree = degree[x];
        ver.pop();
        
        isInVer.insert(x);
        loadNbrForDynamic(x,nbr_u,u_degree,edgeListBegPtr[x],fDat);
        for(int i = 0; i < u_degree; i++){
            if(coreNum[nbr_u[i]] == coreness && vis[nbr_u[i]] == false)
                ver.push(nbr_u[i]), vis[nbr_u[i]] = true;
            else if(coreNum[nbr_u[i]] < coreness)
                deg_copy[x]--;
        }
        pq.push({x,deg_copy[x]});
    }
    uint32_t record = 0;
    while(!pq.empty()){
        VerCore tmp = pq.top();
        pq.pop();
        if(isInVer.find(tmp.u) == isInVer.end()) continue;
        record = max(record,tmp.deg);
        coreNum[tmp.u] = record;
        uint32_t u_degree = degree[tmp.u];
        loadNbrForDynamic(tmp.u,nbr_u,u_degree,edgeListBegPtr[tmp.u],fDat);
        for(int i = 0; i < u_degree; i++){
            if(coreNum[nbr_u[i]] == coreness && isInVer.find(nbr_u[i]) != isInVer.end())
                pq.push({nbr_u[i],--deg_copy[nbr_u[i]]});
        }
        isInVer.erase(tmp.u);
    }

    total_io += fDat.get_total_io();

    delete[] vis;
    delete[] nbr_u;
    delete[] deg_copy;
    fDat.fclose();
}

void Graph::reconCoreGraph(uint32_t *coreNum, readFile &file, readFile &newFile){
    char fileName[200];

    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    bool *isInCoreG = new bool[nodeNum]();
    uint32_t newNode = 0, newEdge = 0;
    std::unordered_map<uint32_t,uint32_t> map;
    for(int i = 0; i < nodeNum; i++)
        if(coreNum[i] >= maxCore){
            isInCoreG[i] = true;
            map[i] = newNode;
            newNode++;
        }
    
    string name = "core_edge_tmp";
    string sub_dir = newFile.m_base + name;
    newFile.createDir(sub_dir);
    unsigned long size = 0,es = 0;
	int num = 0,tmpFile = 0;
    uint32_t u,v,u_nbrNum;
    TEdge* edges = new TEdge[file.memEdges];
    uint32_t *nbr_u = new uint32_t[file.maxDeg];
    for(int i = 0; i < nodeNum; i++){
        if(!isInCoreG[i])continue;
        #ifdef DegSort
        u_nbrNum = degree[i]-(edgeListBegPtrPlus[i]-edgeListBegPtr[i])/sizeof(uint32_t);
        if(u_nbrNum == 0)continue;
        loadInfo(nbr_u,u_nbrNum,edgeListBegPtrPlus[i],fDat);
        for(int j = 0; j < u_nbrNum; j++){
            uint32_t v = nbr_u[j];
            if(!isInCoreG[v])continue;
            newFile.moduleInSaveEdgesUnOrder(map[i],map[v],num,edges,size,es++,tmpFile,sub_dir);
        }
        #else
        loadInfo(nbr_u,degree[i],edgeListBegPtr[i],fDat);
        for(int j = 0; j < degree[i]; j++){
            uint32_t v = nbr_u[j];
            if(!isInCoreG[v])continue;
            if(u > v) continue;
            newFile.moduleInSaveEdgesUnOrder(map[i],map[v],num,edges,size,es++,tmpFile,sub_dir);
        }
        #endif
    }
    sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpFile);
	newFile.saveTmpEdges<TEdge>(edges,size,fileName,[](const TEdge & a, const TEdge & b) {
        if(a.u < b.u)
            return true;
        if( a.u > b.u )
            return false;
        return a.v < b.v;
        });
    log_debug(graphClock_.Count("Save tmp edges done, load %ld edges.",es));
    uint32_t vertexNum = 0, maxDegree = 0;
    newFile.edgeNum = es;
    #ifdef DegSort
    newFile.mergeByDegSort(tmpFile+1, vertexNum, maxDegree, name, true,false);
    #else
	newFile.merge(tmpFile+1, vertexNum, maxDegree, name, true);
    #endif
    
    delete[] edges;
    delete[] isInCoreG;
    delete[] nbr_u;
    total_io += newFile.write_io;
    total_io += tmpFile;
    total_io += fDat.get_total_io();
    fDat.fclose();
}

void Graph::filterGlobalIntoSub(bool *isInSubG, readFile &file, readFile &newFile, bool isDynamicload){
    char fileName[200];
    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );    
    MyReadFile fEid( file.m_eid );
	fEid.fopen( BUFFERED );  
    string name = "sort_edge_tmp";
    string sub_dir = newFile.m_base + name;
    newFile.createDir(sub_dir);
    unsigned long size = 0,es = 0;
	int num = 0,tmpFile = 0;
    uint32_t u,v,u_nbrNum;
    TEdge* edges = new TEdge[file.memEdges];
    uint32_t *nbr_u = new uint32_t[file.maxDeg];
    uint64_t *eid_u = new uint64_t[file.maxDeg];

    uint32_t newNode = 0, newEdge = 0;
    std::unordered_map<uint32_t,uint32_t> map;
    subToGlobal.clear();

    for(int i = 0; i < nodeNum; i++)
        if(isInSubG[i]){
            map[i] = newNode;
            subToGlobal[newNode] = i;
            newNode++;
        }


    for(int i = 0; i < nodeNum; i++){
        if(!isInSubG[i])continue;
        u = i;
        #ifdef DegSort
        u_nbrNum = degree[i]-(edgeListBegPtrPlus[i]-edgeListBegPtr[i])/sizeof(uint32_t);
        if(u_nbrNum == 0)continue;
        loadInfo(nbr_u,u_nbrNum,edgeListBegPtrPlus[i],fDat);
        // loadInfo(eid_u,u_nbrNum,edgeListBegPtrPlus[i]*2,fEid);
        for(int j = 0; j < u_nbrNum; j++){
            uint32_t v = nbr_u[j];
            if(!isInSubG[v])continue;
            if(isDynamicload && m_delBit[u]){
                if(m_dynamicDel[u]->find(v) != m_dynamicDel[u]->end())
                    continue;
            }
            newFile.moduleInSaveEdgesUnOrder(map[u],map[v],num,edges,size,es++,tmpFile,sub_dir);
        }
        #else
        loadInfo(nbr_u,degree[i],edgeListBegPtr[i],fDat);
        // loadInfo(eid_u,degree[i],edgeListBegPtr[i]*2,fEid);
        for(int j = 0; j < degree[i]; j++){
            uint32_t v = nbr_u[j];
            if(!isInSubG[v])continue;
            if(u > v) continue;
            if(isDynamicload && m_delBit[u]){
                if(m_dynamicDel[u]->find(v) != m_dynamicDel[u]->end())
                    continue;
            }
            newFile.moduleInSaveEdgesUnOrder(map[u],map[v],num,edges,size,es++,tmpFile,sub_dir);
        }
        #endif
    }
    sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpFile);
	newFile.saveTmpEdges<TEdge>(edges,size,fileName,[](const TEdge & a, const TEdge & b) {
        if(a.u < b.u)
            return true;
        if( a.u > b.u )
            return false;
        return a.v < b.v;
        });
    es = num;
    log_debug(graphClock_.Count("Save tmp edges done, load %ld edges.",es));
    uint32_t vertexNum = 0, maxDegree = 0;
    newFile.edgeNum = es;
    #ifdef DegSort
    newFile.mergeByDegSort(tmpFile+1, vertexNum, maxDegree, name, true,false);
    #else
	newFile.merge(tmpFile+1, vertexNum, maxDegree, name, true);
    #endif

    total_io = total_io + fDat.get_total_io() + fEid.get_total_io() + newFile.write_io;
    
    delete[] edges;
    delete[] nbr_u;
    delete[] eid_u; 
    fDat.fclose();
    fEid.fclose();
}

void Graph::CoreTrussDecomPlus(readFile &file){
    uint32_t *coreNum = new uint32_t[nodeNum];
    KCore(file,coreNum,false);
    string dir = file.m_base + "kCoreInfo";
    file.createDir(dir);
    dir += "/";
    readFile newFile(dir);

    reconCoreGraph(coreNum,file,newFile);

    Graph maxCoreG(newFile.verNum,newFile.edgeNum);
    maxCoreG.Initial(newFile);
    
    maxCoreG.CountTriangleSSDByDegOrder(newFile,true);


    // secondCoreG.prepareStage(newFile);
    int lower_bound = maxCoreG.Triangles / (maxCoreG.edgeNum-maxCoreG.zero_edge);
    last_sup = maxCoreG.binary(newFile,lower_bound,maxCore-1);

    total_io += maxCoreG.total_io;

    #ifndef Maintenance
    if(maxCoreG.getNodeNum() == last_sup+2){  // this is a prunning case: the subgraph is a clique
        maxKtruss = last_sup, maxKtrussEdge = maxCoreG.edgeNum;
        log_info(graphClock_.Count("Trussness: %d, Edge: %d, io: %lu\n", last_sup+2, maxCoreG.getTrussEdge(last_sup+2),total_io));
        return;
    }
    #endif

    uint32_t left = last_sup, right = maxCore-1, capacity, mid = last_sup;
    if(mid == 0)mid = 1;
    uint32_t TrussEdge = 0, last_mid, Truss = 0;


    bool *isInSubG = new bool[nodeNum]();
    for(int i = 0; i < nodeNum; i++)
        if(coreNum[i] >= last_sup+1)
            isInSubG[i] = true;
    
    delete[] coreNum;


    dir = file.m_base + "subGraphInfo";

    file.createDir(dir);
    dir += "/";
    readFile subFile(dir);
    filterGlobalIntoSub(isInSubG,file,subFile,false);  

    delete[] isInSubG;

    Graph subG(subFile.verNum,subFile.edgeNum);
    subG.Initial(subFile);
    subG.CountTriangleSSDByDegOrder(subFile,true);
    capacity = subG.prefix[subG.maxSup] - subG.prefix[mid-1];

    int current_pid = GetCurrentPid();
    float memory_usage = GetMemoryUsage(current_pid);
    memUsage = max(memory_usage,memUsage);
    log_info(graphClock_.Count("lower_bound done,mid: %d, capacity: %d,prefix[maxSup]: %d, memUsage: %f, io: %lu", mid, capacity,subG.prefix[subG.maxSup],memUsage,total_io));

    if(capacity >= (mid+2)*(mid+1) / 2){
        uint64_t edge_num = 0;
        uint32_t node_num = 0; 

        MyReadFile fSup( subFile.m_suppSort );
        fSup.fopen( BUFFERED );
        // create director in order to fill new subgraph constituted by edges whose support greater than mid
        string dir = file.m_base + "graphInfoCopy";
        subFile.createDir(dir);
        string dir_name = file.m_base + "graphInfoCopy/" + to_string(mid);
        subFile.createDir(dir_name);
        dir_name += "/";

        readFile newFile(dir_name);
        unsigned long size = 0,es = 0,eid;
        uint32_t max_degree = 0;
        int num = 0, tmpFile = 0;
        uint32_t u,v;

        memset(newFile.m_vertexMap,-1,sizeof(int)*newFile.m_maxID);
        string name = "sort_edge_tmp";
        string sub_dir = newFile.m_base + name;
        newFile.createDir(sub_dir);
        TEdge* edges = new TEdge[file.memEdges];
        char fileName[150];


        EdgeSup tmp;
        for(int i = subG.prefix[mid-1] ; i < subG.prefix[subG.maxSup]; i++){
            fSup.fseek(i*sizeof(EdgeSup));
            fSup.fread(&tmp,sizeof(EdgeSup));
            u = tmp.u;
            v = tmp.v;
            eid = tmp.eid;
            newFile.moduleInSaveEdgesUnOrder(u,v,num,edges,size,eid,tmpFile,sub_dir);
        }
        total_io += fSup.get_total_io();

        fSup.fclose();

        sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpFile);
        newFile.saveTmpEdges<TEdge>(edges,size,fileName,[](const TEdge & a, const TEdge & b) {
                    if(a.u < b.u)
                        return true;
                    if( a.u > b.u )
                        return false;
                    return a.v < b.v;
                    });
        newFile.edgeNum = num;
        delete[] edges;

        // newFile.merge(tmpFile+1, node_num, max_degree, name, false);
        newFile.mergeByDegSort(tmpFile+1, node_num, max_degree, name, false,false);

        int current_pid = GetCurrentPid();
        float memory_usage = GetMemoryUsage(current_pid);
        memUsage = max(memory_usage,memUsage);

        log_info(graphClock_.Count("new subgraph vertex: %d, edge: %lu, memUsage: %f",newFile.verNum,newFile.edgeNum,memUsage));
        total_io += newFile.write_io;   
        
        Graph tmp_g(subG.nodeNum,subG.edgeNum);
        // tmp_g.InitialUnOrder(newFile);
        tmp_g.InitialUnOrderDegSort(newFile);


        /* prune optimization -- kcore */
        if(!tmp_g.peelVertex(mid,newFile,false)) 
            printf("not exist\n");
        
        tmp_g.CountTriangleSSDByDegOrder(newFile,true);
        log_debug(graphClock_.Count("mid: %d, num: %d",mid,tmp_g.prefix[tmp_g.maxSup]-tmp_g.prefix[mid-1]));

        /* prune optimization -- delete edges on original graph, so that avoid reconstructing subgraph */
        if(tmp_g.minSup >= mid)
        {
            mid = tmp_g.minSup;
            log_debug(graphClock_.Count("minSup >= mid, mid: %d",mid));
        }
        std::string path = file.m_base + "linearList/";
        file.createDir(path);
        ListLinearHeapTruss *linear_heap = new ListLinearHeapTruss(edgeNum,tmp_g.maxSup,path,newFile.m_supp);
        DynamicHeap dheap = DynamicHeap(100000000);
        if(tmp_g.existTrussLazyUpdate(newFile,mid,TrussEdge,Truss,linear_heap,dheap)){
            last_mid = mid;
            for(uint32_t i = mid+1; i <= tmp_g.maxSup; i++){
                isDeleteEdgeByFunc = true;
                if(!tmp_g.deleteEdgeLazyUpdateTest(file,newFile,last_mid,i,prefix,TrussEdge,Truss,linear_heap,dheap))
                {
                    Truss = last_mid;
                    break;
                }
                last_mid = i;
            }
        }
        delete linear_heap;
        total_io += tmp_g.total_io;
        maxKtruss = Truss;
        saveKtruss = dir_name;
        maxKtrussEdge = subG.edgeNum;

        #ifdef Maintenance
        string isDeleteEdgeByFunc_ss = file.m_base + "graph.isDeleteEdgeByFunc";
        FILE* fp=fopen(isDeleteEdgeByFunc_ss.c_str(),"wb");
        fwrite(&isDeleteEdgeByFunc,sizeof(bool),1,fp);
        fclose(fp);

        string fileVertexId_ss = file.m_base + "graph.fileVertexId";
        fp=fopen(fileVertexId_ss.c_str(),"wb");
        int len = newFile.verNum;
        fwrite(&len,sizeof(int),1,fp);
        for (uint32_t u = 0; u < newFile.verNum; ++u){
            file.vertexId.push_back(newFile.vertexId[u]);
            fwrite(&newFile.vertexId[u],sizeof(int),1,fp);
        }
        fclose(fp);


        //save information of subToGlobal datastructure
        string subToGlobal_ss = file.m_base +"graph.subToGlobal";
        fp=fopen(subToGlobal_ss.c_str(),"wb");
        for(auto it = subToGlobal.begin(); it != subToGlobal.end(); it++){
            Edge e = {it->first,it->second};
            fwrite(&e,sizeof(Edge),1,fp);
        }
        fclose(fp);

        string finalInfo_ss = file.m_base +"graph.final";
        FILE* fInfo=fopen(finalInfo_ss.c_str(),"wb");
        fwrite(&last_sup,sizeof(uint32_t),1,fInfo);
        fwrite(&maxKtruss,sizeof(uint32_t),1,fInfo);
        fclose(fInfo);
        #endif
        
        log_debug(graphClock_.Count("Trussness: %u, Edge: %u, io: %lu\n", Truss+2, TrussEdge,total_io));
    }

}

int Graph::IntersectTriangleSSDPlus(uint32_t u, uint32_t *nbr_u, uint32_t v, uint32_t *nbr_v)
{
    int ptr_u = 0;
    int upper_u = degree[u];
    int ptr_v = 0;
    int upper_v = degree[v];

    if(nbr_u[upper_u-1] < nbr_v[0] || nbr_v[upper_v-1] < nbr_u[0])
        return 0;

    int count = 0;
    while(ptr_u < upper_u && ptr_v < upper_v){
        if(vis[nbr_u[ptr_u]] || nbr_u[ptr_u] < nbr_v[ptr_v])
            ptr_u++;
        else if(vis[nbr_v[ptr_v]] || nbr_u[ptr_u] > nbr_v[ptr_v])
            ptr_v++;
        else
        {
            ptr_u++;
            ptr_v++;
            count++;
        }
    }
    return count;
}

void Graph::CountTriangleSSDPlus(readFile &file, bool isOrder){
    log_info(graphClock_.Count("Enter CountTriangleSSDPlus function"));
    file.write_io = 0;
	MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fEid( file.m_eid );
	fEid.fopen( BUFFERED );
    MyReadFile fOff( file.m_offset );
	fOff.fopen( BUFFERED );
    EdgeSup* edges = new EdgeSup[file.memEdges];
    FILE* fSupp = fopen(file.m_supp.c_str(),"wb");

    uint32_t triangle_count = 0, size = 0, tmpfile = 0, c = 0;
    uint32_t upper,tmpi,i;
    char fileName[200];

    string sub_dir = file.m_base + "support_sort_edge";
    file.createDir(sub_dir);

    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* nbr_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    if(isOrder)
        upper = nodeNum;
    else
        upper = file.verNum;
    for(uint32_t tmpi = 0; tmpi < upper; tmpi++){
        if(isOrder)
            i = tmpi;
        else
            i = file.vertexId[tmpi];
        if(vis[i])continue;

        loadInfo(nbr_u,degree[i],edgeListBegPtr[i],fDat);
        uint32_t u = i;
        for(uint32_t j = 0; j < degree[i]; j++){
            uint32_t v = nbr_u[j];
            if(vis[v])continue;
            // if(!(degree[v]<degree[u] || (degree[u]==degree[v] && v<u)))
            //     continue;
            if(u > v)
                continue;

            loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
            uint32_t sup = IntersectTriangleSSDPlus(u,nbr_u,v,nbr_v);
            if(sup == 0)continue;

            uint64_t offset = edgeListBegPtr[i]*2+j*sizeof(uint64_t);
            uint64_t eid_uv;
            fEid.fseek(offset);
            fEid.fread(&eid_uv,sizeof(uint64_t));

            maxSup = max(sup,maxSup);
            minSup = min(sup,minSup);
            edges[size].u = u;
            edges[size].v = v;
            edges[size].sup = sup; // 这里可以优化，减少写入, if sup == 0, continue? it is not necessary to write.
            edges[size].eid = eid_uv;
            size++;
            c++;
            if(size >= file.memEdges){            
                sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpfile);
                file.saveTmpEdges<EdgeSup>(edges,size,fileName,[](const EdgeSup & a, const EdgeSup & b) {
                    return a.sup < b.sup;
                    });
                size = 0;
                ++tmpfile;
            }
            //count range of edges' support
            prefix[sup]++;
            eid_eid tmp;
            fOff.fseek(eid_uv*sizeof(eid_eid));
            fOff.fread(&tmp,sizeof(eid_eid));

            fseek(fSupp,tmp.first*sizeof(uint32_t),SEEK_SET);
            fwrite(&sup,sizeof(uint32_t),1,fSupp);
            fseek(fSupp,tmp.second*sizeof(uint32_t),SEEK_SET);
            fwrite(&sup,sizeof(uint32_t),1,fSupp);
            total_io += 2;
            triangle_count += sup;
        }
    }
    sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpfile);
    file.saveTmpEdges<EdgeSup>(edges,size,fileName,[](const EdgeSup & a, const EdgeSup & b) {
        return a.sup < b.sup;
    });
    
    file.mergeBySup<EdgeSup>(tmpfile+1);
    int current_pid = GetCurrentPid();
    float memory_usage = GetMemoryUsage(current_pid);
    memUsage = max(memory_usage,memUsage);

    Triangles = triangle_count/3;
    log_debug(graphClock_.Count("Original edges: %u, Now edges: %u, memUsage: %f",file.edgeNum,c,memUsage));
    file.edgeNum = c;
    total_io += file.write_io;

    // print support of all edges
    sortSupport(file);

    total_io += fDat.get_total_io();
    total_io += fEid.get_total_io();
    total_io += fOff.get_total_io();

    fDat.fclose();
    fEid.fclose();
    fOff.fclose();
    fclose(fSupp); 
    free(nbr_u);
    free(nbr_v);
    delete []edges;
    log_debug(graphClock_.Count("Triangles: %d, MaxSupport: %u, MinSupport: %u",Triangles,maxSup,minSup));
}


void Graph::sortSupport(readFile &file){
    // calculate the prefix of edges' support
    // printf("sup 0 : %d\n",prefix[0]);
    for(int i = 1; i <= maxSup; i++){ 
        prefix[i] += prefix[i-1]; 
        // printf("sup %d : %d\n",i,prefix[i]-prefix[i-1]);
    }
}

uint32_t Graph::binary(readFile &file, uint32_t start, uint32_t end){
    uint32_t left = start, right = end, capacity, mid;
    uint32_t TrussEdge = 0, last_mid, Truss = 0;

    log_info(graphClock_.Count("binary func => start: %d, end: %d",start,end));

    while(left <= right){
        mid = (left+right)/2;
        last_mid = mid;
        capacity = prefix[maxSup]-prefix[mid-1];
        // log_info(graphClock_.Count("mid: %d, capacity: %d,prefix[maxSup]: %d", mid, capacity,prefix[maxSup]));
        if(capacity < (mid+2)*(mid+1) / 2)
        {
            right = mid-1;
            continue;
        }
        if(!inducedGraphUnOrder(left,right,mid,file,TrussEdge,Truss,true))
        {
            log_debug(graphClock_.Count("error mid: %d, last_mid: %d",mid,last_mid));
            if(mid != last_mid){
                break;
            }
            else
                right = mid-1;
        }
        else
            left = mid+1;
    }

    log_debug(graphClock_.Count("final last_sup: %d\n",mid));
    return mid;
}


void Graph::binaryImproved(readFile &file, uint32_t start){
    uint32_t left = start, right = maxSup, capacity = 0, mid = 0;
    uint32_t TrussEdge = 0, last_mid = 0, Truss = 0;
    while(left <= right){
        mid = (left+right)/2;
        last_mid = mid;
        capacity = prefix[maxSup]-prefix[mid-1];

        if(capacity < (mid+2)*(mid+1) / 2)
        {
            // log_info(graphClock_.Count("xxx mid: %u",mid));
            right = mid-1;
            continue;
        }

        /* ID of vertex is same as original graph */
        if(!inducedGraphUnOrder(left,right,mid,file,TrussEdge,Truss,false))
        {
            right = mid-1;
        }
        else
            left = mid+1;
    }
    total_io += file.write_io;
    log_info(graphClock_.Count("Trussness: %d, Edge: %d, io: %lu\n", Truss+2, TrussEdge, total_io));
}

void Graph::binaryAndIncremental(readFile &file, uint32_t start){
    uint32_t left = start, right = maxSup;
    uint32_t capacity, mid;
    uint32_t TrussEdge = 0, last_mid, Truss = 0;
    while(left <= right){
        mid = (left+right)/2;
        last_mid = mid;
        capacity = prefix[maxSup]-prefix[mid-1];
        printf("mid: %d, capacity: %d,prefix[maxSup]: %d\n", mid, capacity,prefix[maxSup]);
        if(capacity < (mid+2)*(mid+1) / 2)
        {
            right = mid-1;
            continue;
        }
        /* ID of vertex is the same as original graph */
        if(!inducedGraphUnOrder(left,right,mid,file,TrussEdge,Truss,true))
        {
            printf("error mid: %d, last_mid: %d\n",mid,last_mid);
            if(mid != last_mid){
                break;
            }
            else
                right = mid-1;
        }
        else
            left = mid+1;
    }

    capacity = prefix[maxSup]-prefix[mid-1];
    printf("mid: %d, capacity: %d,prefix[maxSup]: %d\n", mid, capacity,prefix[maxSup]);

    if(capacity >= (mid+2)*(mid+1) / 2){
        uint64_t edge_num = 0;
        uint32_t node_num = 0; 

        MyReadFile fSup( file.m_suppSort );
        fSup.fopen( BUFFERED );
        // create director in order to fill new subgraph constituted by edges whose support greater than mid
        string dir =  file.m_base + "graphInfoCopy/";
        file.createDir(dir);
        string dir_name = dir + to_string(mid);
        file.createDir(dir_name);
        dir_name += "/";
        readFile newFile(dir_name);
        unsigned long size = 0,es = 0,eid;
        uint32_t max_degree = 0;
        int num = 0, tmpFile = 0;
        uint32_t u,v;

        memset(newFile.m_vertexMap,-1,sizeof(int)*newFile.m_maxID);
        string name = "sort_edge_tmp";
        string sub_dir = newFile.m_base + name;
        newFile.createDir(sub_dir);
        TEdge* edges = new TEdge[file.memEdges];
        char fileName[150];


        EdgeSup tmp;
        for(int i = prefix[mid-1] ; i < prefix[maxSup]; i++){
            fSup.fseek(i*sizeof(EdgeSup));
            fSup.fread(&tmp,sizeof(EdgeSup));
            u = tmp.u;
            v = tmp.v;
            eid = tmp.eid;
            newFile.moduleInSaveEdgesUnOrder(u,v,num,edges,size,eid,tmpFile,sub_dir);
        }
        fSup.fclose();

        sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpFile);
        newFile.saveTmpEdges<TEdge>(edges,size,fileName,[](const TEdge & a, const TEdge & b) {
                    if(a.u < b.u)
                        return true;
                    if( a.u > b.u )
                        return false;
                    return a.v < b.v;
                    });
        // log_debug(graphClock_.Count("Save new tmp edges done, load %ld edges.",es));
        newFile.edgeNum = num;
        delete[] edges;
        newFile.merge(tmpFile+1, node_num, max_degree, name, false);
        log_debug(graphClock_.Count("new subgraph vertex: %d, edge: %lu\n",newFile.verNum,newFile.edgeNum));

        Graph tmp_g(nodeNum,edgeNum);
        tmp_g.InitialUnOrder(newFile);
        // tmp_g.InitialUnOrderDegSort(newFile);


        /* prune optimization -- kcore */
        if(!tmp_g.peelVertex(mid,newFile,false)) 
            printf("not exist\n");
        
        tmp_g.CountTriangleSSDPlus(newFile,true);
        
        log_debug(graphClock_.Count("mid: %d, num: %d",mid,tmp_g.prefix[tmp_g.maxSup]-tmp_g.prefix[mid-1]));

        /* prune optimization -- delete edges on original graph, so that avoid reconstructing subgraph */
        if(tmp_g.minSup >= mid)
        {
            mid = tmp_g.minSup;
            log_debug(graphClock_.Count("minSup >= mid, mid: %d",mid));
        }
        if(tmp_g.existTrussPlus(newFile,mid,TrussEdge,Truss)){
            last_mid = mid;
            for(uint32_t i = mid+1; i <= maxSup; i++){
                if(!tmp_g.deleteEdge(file,newFile,last_mid,i,prefix,TrussEdge,Truss))
                {
                    Truss = last_mid;
                    break;
                }
                last_mid = i;
            }
        }
    }
    printf("Trussness: %d, Edge: %d\n", Truss+2, TrussEdge);

}


void Graph::constructBin(readFile &file){
    int start = 0;
    int count=0,binsize;
    memset(bin,0,sizeof(uint32_t)*nodeNum);
	for (int i = 0; i <= maxSup; ++i) {
		binsize = prefix[i]-start;
        start = prefix[i];
		bin[i] = count;
		count += binsize;
	}
    bin[maxSup+1] = count;


    FILE* fBinEdge = fopen(file.m_binEdge.c_str(),"wb");
    FILE* fEdgePos = fopen(file.m_ePos.c_str(),"wb");
    
    MyReadFile fSup( file.m_suppSort );
	fSup.fopen( BUFFERED );

    EdgeSup tmp;
    TEdge e;
    for(uint64_t i = 0; i < file.edgeNum; i++){
        fSup.fseek(i*sizeof(EdgeSup));
        fSup.fread(&tmp,sizeof(EdgeSup));
        
        e = {tmp.u,tmp.v,tmp.eid};
        uint32_t sup = tmp.sup;

        int pos = bin[sup];
        // printf("i: %u,u: %u,v: %u,sup: %u,eid: %lu,pos: %u\n",i,tmp.u,tmp.v,tmp.sup,tmp.eid,pos);
        fseek(fBinEdge,pos*sizeof(TEdge),SEEK_SET);
        fwrite(&e,sizeof(TEdge),1,fBinEdge);

        fseek(fEdgePos,tmp.eid*sizeof(int),SEEK_SET);
        fwrite(&pos,sizeof(int),1,fEdgePos);
        bin[sup]++;
    }

    for(int i = maxSup; i > 0; i--) bin[i] = bin[i-1];
    bin[0] = 0;
    total_io += 3*file.edgeNum;


    fclose(fBinEdge);
    fclose(fEdgePos);
    fSup.fclose();
}


bool Graph::bottomUpDecom(readFile &file, uint32_t &mid, uint32_t &TrussEdge, uint32_t &Truss){
    uint32_t m = file.edgeNum;
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
    uint32_t sup;

    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* nbr_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    uint64_t* eid_u = (uint64_t *)malloc(sizeof(uint64_t) * file.maxDeg);
    uint64_t* eid_v = (uint64_t *)malloc(sizeof(uint64_t) * file.maxDeg);

    uint32_t* sup_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* sup_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    TEdge edge;
    uint32_t u,v;
    uint64_t eid_uv,uv_pos;
    int finish_sup = 0;

    for(int s = 0; s < m; s++){
        fBinEdge.fseek(s*sizeof(TEdge));
        fBinEdge.fread(&edge,sizeof(TEdge));
        u = edge.u;
		v = edge.v;
        eid_uv = edge.eid; 
        eid_eid tmp;
        fOff.fseek(eid_uv*sizeof(eid_eid));
        fOff.fread(&tmp,sizeof(eid_eid));
        
        fSup.fseek(tmp.first*sizeof(uint32_t));
        fSup.fread(&sup,sizeof(uint32_t));
        
        // fSup.fseek(eid_uv*sizeof(uint32_t));
        // fSup.fread(&sup,sizeof(uint32_t));

        if(m-s < (mid+2)*(mid+1) / 2){
            mid = sup;
            total_io = total_io + fSup.get_total_io() + fEid.get_total_io() + fDat.get_total_io() + fEdgePos.get_total_io() + fBinEdge.get_total_io();
            log_debug(graphClock_.Count("error last_sup: %d",sup));
            fSup.fclose();
            fEid.fclose();
            fDat.fclose();
            fOff.fclose();
            fEdgePos.fclose();
            fBinEdge.fclose();
            free(nbr_u);
            free(nbr_v);
            free(eid_u);
            free(eid_v);
            free(sup_u);
            free(sup_v);
            return false;
        }
            
        
        if(sup >= mid){
            TrussEdge = m-s;
            Truss = sup;
            last_s = s;
            last_sup = sup;
            total_io = total_io + fSup.get_total_io() + fEid.get_total_io() + fDat.get_total_io() + fEdgePos.get_total_io() + fBinEdge.get_total_io();

            log_debug(graphClock_.Count("success last_sup: %u, s: %d, TrussEdge: %u",sup,s,TrussEdge));
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
            return true;
        }
        
        loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
        loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
        loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);

        loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
        loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
        loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);

        /* optimization */
        std::vector<eid_eid> comm;
        IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,sup,comm,true);
        for(int i = 0; i < comm.size(); i++){
            // updateEdgeSSD(u,comm[i].w,comm[i].first,comm[i].u_sup,sup,fEdgePos,fSup,fBinEdge);
			// updateEdgeSSD(v,comm[i].w,comm[i].second,comm[i].v_sup,sup,fEdgePos,fSup,fBinEdge);
            updateEdgeSSDNew(u,nbr_u[comm[i].first],eid_u[comm[i].first],comm[i].first,sup,sup_u,fEdgePos,fSup,fBinEdge,fOff);
			updateEdgeSSDNew(v,nbr_u[comm[i].first],eid_v[comm[i].second],comm[i].second,sup,sup_v,fEdgePos,fSup,fBinEdge,fOff);
        }
        finish_sup = sup;
        sup = DELETE(sup);
        --degree_[u];
        --degree_[v];
        // fSup.fseek(eid_uv*sizeof(uint32_t));
        // fSup.fwrite(&sup,sizeof(uint32_t));
        fSup.fseek(tmp.first*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t));

        fSup.fseek(tmp.second*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t));

    }

    fSup.fclose();
    fEid.fclose();
    fDat.fclose();
    fEdgePos.fclose();
    fBinEdge.fclose();
    free(nbr_u);
    free(nbr_v);
    free(eid_u);
    free(eid_v);

    total_io = total_io + fSup.get_total_io() + fEid.get_total_io() + fDat.get_total_io() + fEdgePos.get_total_io() + fBinEdge.get_total_io();

    log_debug(graphClock_.Count("error finish sup: %d",finish_sup));

    return false;
}

bool Graph::existTrussPlus(readFile &file, uint32_t &mid, uint32_t &TrussEdge, uint32_t &Truss){

    constructBin(file);
    return bottomUpDecom(file,mid,TrussEdge,Truss);
}


bool Graph::peelVertex(int mid, readFile &file, bool flag){
    log_debug(graphClock_.Count("enter peelVertex function, mid: %d",mid));

    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );

    std::queue<int> q;
    bool *isIn = (bool *)calloc(nodeNum,sizeof(bool));
    int edgeNumber = file.edgeNum;
    int c = 0, maxDegree = 0;
    uint32_t upper, i;
    if(flag)
        upper = nodeNum;
    else
        upper = file.verNum;
    for(int u = 0; u < upper; u++)
    {
        if(flag)
            i = u;
        else
            i = file.vertexId[u];
        if(!vis[i] && degree_[i] <= mid)
        {
            q.push(i);
            isIn[i] = true;
        }             
    }
    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    while(!q.empty()){
        int ver = q.front();
        q.pop();
        vis[ver] = true;
        c++;
        loadInfo(nbr_u,degree[ver],edgeListBegPtr[ver],fDat);
        for(uint32_t j = 0; j < degree[ver]; j++){
            uint32_t v = nbr_u[j];
            if(!vis[v]){
                degree_[v]--;
                edgeNumber--;
                if(edgeNumber < (mid+2)*(mid+1) / 2)
                    return false;
                if(!isIn[v] && degree_[v] <= mid)
                {
                    q.push(v);
                    isIn[v] = true;
                }    
            }
        }
    }
    total_io += fDat.get_total_io();
    fDat.fclose();
    free(isIn);
    free(nbr_u);
    log_info(graphClock_.Count("edge :%ld reduce: %ld, vertex :%d reduce: %d\n",file.edgeNum,file.edgeNum-edgeNumber,file.verNum,c));
    return true;
}

void Graph::prepareStage(readFile &file){

    MyReadFile fSuppSort( file.m_suppSort );
	fSuppSort.fopen( BUFFERED );
    uint64_t sup_sum = 0;

    for(uint32_t i = 0; i < file.edgeNum; i++)
    {
        uvSup tmp;
        fSuppSort.fread(&tmp,sizeof(uvSup));
        uint32_t sup = tmp.sup;
        assert(sup <= nodeNum);
        sup_sum += sup;
        prefix[sup]++;
        maxSup = max(maxSup,sup);
        minSup = max(minSup,sup);
    }
    sortSupport(file);
    Triangles = sup_sum / 3;
    fSuppSort.fclose();
}

void Graph::CountTriangleSSDByDegOrder(readFile &file, bool saveSupEdge){
    int threshold;
    log_info(graphClock_.Count("CountTriangleSSDByDegOrder begin"));
	MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fEid( file.m_eid );
	fEid.fopen( BUFFERED );
    MyReadFile fOff( file.m_offset );
	fOff.fopen( BUFFERED );
    EdgeSup* edges; 
    if(saveSupEdge)
        edges = new EdgeSup[file.memEdges];

    FILE* fSupp = fopen(file.m_supp.c_str(),"wb");

    uint32_t upper,i,triangle_count = 0,c=0, size = 0, tmpfile = 0, u_nbrNum = 0, v_nbrNum = 0, mapSize = 0;
    char fileName[200];

    string sub_dir = file.m_base + "support_sort_edge";
    file.createDir(sub_dir);

    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* nbr_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    std::unordered_map<uint64_t,uint32_t> map_pos;

    log_debug(graphClock_.Count("file.vertexId.size %u",file.vertexId.size()));

    bool isVerInOri = false;
    upper = nodeNum;
    if(file.vertexId.size() != 0){
        isVerInOri = true;
        upper = file.verNum;
    }
    
    for(uint32_t tmpi = 0; tmpi < upper; tmpi++){
        i = tmpi;
        if(isVerInOri)
            i = file.vertexId[tmpi];
        if(vis[i]) continue;
        u_nbrNum = degree[i]-(edgeListBegPtrPlus[i]-edgeListBegPtr[i])/sizeof(uint32_t);
        if(u_nbrNum == 0)continue;
        loadInfo(nbr_u,u_nbrNum,edgeListBegPtrPlus[i],fDat);
        uint32_t u = i;
        
        for(uint32_t j = 0; j < u_nbrNum; j++){
            uint32_t v = nbr_u[j];
            if(vis[v]) continue;
            v_nbrNum = degree[v]-(edgeListBegPtrPlus[v]-edgeListBegPtr[v])/sizeof(uint32_t);

            loadInfo(nbr_v,v_nbrNum,edgeListBegPtrPlus[v],fDat);
            std::vector<eid_eid> common;
            assert(u_nbrNum <= file.maxDeg);
            assert(v_nbrNum <= file.maxDeg);
            uint32_t sup = IntersectTriangleSSDByDegOrder(u,u_nbrNum,nbr_u,v,v_nbrNum,nbr_v,fDat,fEid,common,map_pos);
            
            triangle_count += sup;
            uint64_t offset = edgeListBegPtrPlus[i]*2+j*sizeof(uint64_t);
            map_pos[offset] += sup;
            
            uint64_t eid_uv = 0; //?
            if(map_pos[offset] == 0){
                zero_edge++;
                map_pos.erase(offset);
                // eid_uv = DELETE(eid_uv);  //affect bin struct process of naive truss decom                
                // fEid.fseek(offset);
                // fEid.fwrite(&eid_uv,sizeof(uint64_t));
                continue;
            }
            
            uint32_t sup_uv = 0;
            
            fEid.fseek(offset);
            fEid.fread(&eid_uv,sizeof(uint64_t));
            // printf("u: %u, v: %u, sup: %u, offset: %lu, eid_uv: %lu\n",u,v,sup,offset,eid_uv);

            eid_eid tmp;
            fOff.fseek(eid_uv*sizeof(eid_eid));
            fOff.fread(&tmp,sizeof(eid_eid));

            sup_uv = map_pos[offset];
            assert(sup_uv<=nodeNum);

            maxSup = max(sup_uv,maxSup);
            minSup = min(minSup,sup_uv);
            fseek(fSupp,tmp.first*sizeof(uint32_t),SEEK_SET);
            fwrite(&sup_uv,sizeof(uint32_t),1,fSupp);
            fseek(fSupp,tmp.second*sizeof(uint32_t),SEEK_SET);
            fwrite(&sup_uv,sizeof(uint32_t),1,fSupp);
            total_io += 2;
            c++;
            
            mapSize = mapSize > map_pos.size() ? mapSize : map_pos.size();
            map_pos.erase(offset);
            if(saveSupEdge){
                edges[size].u = u;
                edges[size].v = v;
                edges[size].sup = sup_uv; 
                edges[size].eid = eid_uv;
                size++;                
                if(size >= file.memEdges){            
                    sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpfile);
                    file.saveTmpEdges<EdgeSup>(edges,size,fileName,[](const EdgeSup & a, const EdgeSup & b) {
                        return a.sup < b.sup;
                        });
                    size = 0;
                    ++tmpfile;
                }
                //count range of edges' support
                prefix[sup_uv]++;
            }
        }
    }
    
    if(saveSupEdge){
        // if(tmpfile)
        {
            sortSupport(file);
            sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpfile);
            file.saveTmpEdges<EdgeSup>(edges,size,fileName,[](const EdgeSup & a, const EdgeSup & b) {
                return a.sup < b.sup;
            });
            file.mergeBySup<EdgeSup>(tmpfile+1);     
                   
        }
        // else //need not to write edges into ssd
        // {
        //     threshold = triangle_count / edgeNum;
        //     std::sort(edges,edges+size,[](const EdgeSup & a, const EdgeSup & b) {
        //         return a.sup < b.sup;
        //     });
        //     int cumulate = 0;
        //     for(int i = 0; i <= maxSup; i++){
        //         if(i < threshold){
        //             cumulate += prefix[i];
        //         }
        //         else if(i > threshold){
        //             prefix[i] += prefix[i-1];
        //         }
        //     }
        //     FILE *fEdgeSup = fopen(file.m_suppSort.c_str(),"wb");
        //     c = size-cumulate;
        //     fwrite( edges+cumulate, sizeof(EdgeSup), size-cumulate, fEdgeSup );
        //     fclose(fEdgeSup);
        // }
    }
    int current_pid = GetCurrentPid();
    float memory_usage = GetMemoryUsage(current_pid);
    memUsage = max(memory_usage,memUsage);
    
    log_debug(graphClock_.Count("Original edges: %u, Now edges: %u, memUsage: %lf, minSup: %u, maxSup: %u",file.edgeNum,c,memUsage,minSup,maxSup));
    file.edgeNum = c;
    Triangles = triangle_count;
    total_io += fDat.get_total_io();
    total_io += fEid.get_total_io();
    total_io += fOff.get_total_io();

    fDat.fclose();
    fEid.fclose();
    fOff.fclose();
    fclose(fSupp); 
    free(nbr_u);
    free(nbr_v);
    if(saveSupEdge)
        delete []edges;
    // log_debug(graphClock_.Count("Triangles: %u, MaxSupport: %u, MinSupport: %u, MapSize: %u",Triangles,maxSup,minSup,mapSize));
}

bool Graph::existTrussLazyUpdate(readFile &file, uint32_t &mid, uint32_t &TrussEdge, uint32_t &Truss, ListLinearHeapTruss *linear_heap, DynamicHeap &dheap){
    MyReadFile fSupSort( file.m_suppSort );
	fSupSort.fopen( BUFFERED );
    FILE* filePres = fopen(linear_heap->m_pres.c_str(),"wb");
	FILE* fileNexts = fopen(linear_heap->m_nexts.c_str(),"wb");
    FILE* fileEidToV = fopen(file.m_eidToVer.c_str(),"wb");

    EdgeSup tmp_;
    Edge e;
    for(int i = 0; i < file.edgeNum; i++){
        fSupSort.fseek(i*sizeof(EdgeSup));
        fSupSort.fread(&tmp_,sizeof(EdgeSup));
        e = {tmp_.u,tmp_.v};
        fseek(fileEidToV,tmp_.eid*sizeof(Edge),SEEK_SET);
        fwrite(&e,sizeof(Edge),1,fileEidToV);
        uint32_t sup = tmp_.sup;
        linear_heap->insert(tmp_.eid,tmp_.sup,filePres,fileNexts);
    }
    total_io += edgeNum;
    total_io += fSupSort.get_total_io();
    total_io += linear_heap->total_io;
    linear_heap->total_io = 0;

    fclose(filePres);
    fclose(fileNexts);
    fclose(fileEidToV);
    fSupSort.fclose();
    log_info(graphClock_.Count("finish constructing linearList"));

    MyReadFile fPres( linear_heap->m_pres );
	fPres.fopen( NOBUFFER );
    MyReadFile fNexts( linear_heap->m_nexts );
	fNexts.fopen( NOBUFFER );
    MyReadFile fEidToVer( file.m_eidToVer );
	fEidToVer.fopen( BUFFERED );
    MyReadFile fSup( file.m_supp );
	fSup.fopen( NOBUFFER );
    MyReadFile fOff( file.m_offset );
	fOff.fopen( BUFFERED );


    
    ui u,v;
    ui last_sup = ui(-1);
    int max_size = 0;
    Edge tmp;

    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fEid( file.m_eid );
	fEid.fopen( BUFFERED );
    struct timeval start_time, end_time;
    double insetsect_interval = 0;
    double update_interval = 0;
    uint32_t sup,u_nbrNum;
    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* nbr_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    uint64_t* eid_u = (uint64_t *)malloc(sizeof(uint64_t) * file.maxDeg);
    uint64_t* eid_v = (uint64_t *)malloc(sizeof(uint64_t) * file.maxDeg); 

    uint32_t* sup_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* sup_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);



    for(ui64 i = 0; i < file.edgeNum; i++){
        gettimeofday(&start_time, NULL);
        ui64 eid;
        linear_heap->pop_min(eid,sup,fPres,fNexts);

        if(linear_heap->exist_num+1 < (mid+2)*(mid+1) / 2){
            mid = sup;
            log_debug(graphClock_.Count("error last_sup: %d",sup));
            log_debug(graphClock_.Count("heap max size: %u, intersect_time: %lf, update_time: %lf",max_size,insetsect_interval,update_interval));
            total_io = total_io+fDat.get_total_io()+fEid.get_total_io()+fPres.get_total_io()+fNexts.get_total_io()+fEidToVer.get_total_io()+fSup.get_total_io()+fOff.get_total_io();

            fSup.fclose();
            fEid.fclose();
            fDat.fclose();
            fPres.fclose();
            fNexts.fclose();
            fEidToVer.fclose();
            fOff.fclose();
            free(nbr_u);
            free(nbr_v);
            free(eid_u);
            free(eid_v);
            free(sup_u);
            free(sup_v);
            return false;
        }
            
        
        if(sup >= mid){
            Truss = sup;
            last_s = i;
            last_sup = sup;
            linear_heap->insert(eid,sup,fOff,fSup,fPres,fNexts);
            TrussEdge = linear_heap->exist_num;
            
            log_debug(graphClock_.Count("success last_sup: %u, i: %d, TrussEdge: %u, %u",sup,i,TrussEdge,linear_heap->exist_num));
            log_debug(graphClock_.Count("heap max size: %u, intersect_time: %lf, update_time: %lf",max_size,insetsect_interval,update_interval));
            total_io = total_io+fDat.get_total_io()+fEid.get_total_io()+fPres.get_total_io()+fNexts.get_total_io()+fEidToVer.get_total_io()+fSup.get_total_io()+fOff.get_total_io();
            fSup.fclose();
            fEid.fclose();
            fDat.fclose();
            fPres.fclose();
            fNexts.fclose();
            fEidToVer.fclose();
            fOff.fclose();
            free(nbr_u);
            free(nbr_v);
            free(eid_u);
            free(eid_v);
            free(sup_u);
            free(sup_v);
            return true;
        }
        fEidToVer.fseek(eid*sizeof(Edge));
        fEidToVer.fread(&tmp,sizeof(Edge));
        u = tmp.u;
        v = tmp.v;
        printClass(u,v,sup+2);

        // if(sup != last_sup)
        // {
        //     log_info(graphClock_.Count("i: %u, u: %u, v: %u, eid: %u, sup: %u, size: %d",i,u,v,eid,sup,dheap.size));
        //     last_sup = sup;
        //     // linear_heap->print();
        // }
        loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
        loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
        loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);
        
        loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
        loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
        loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);

        std::vector<eid_eid> comm;
        IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,sup,comm,true);

        gettimeofday(&end_time, NULL);
        insetsect_interval += (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec)/1000000.0;
        gettimeofday(&start_time, NULL);
        for(int i = 0; i < comm.size(); i++){
            ui first_sup = sup_u[comm[i].first];
            ui second_sup = sup_v[comm[i].second];
            uint64_t eid_fir = eid_u[comm[i].first];
            uint64_t eid_sec = eid_v[comm[i].second];

            if(dheap.find(eid_fir)){
                if(dheap.getSup(eid_fir) == sup+1){
                    ui tmp = first_sup - (sup+1);
                    linear_heap->decrement(first_sup,eid_fir,fOff,fSup,fPres,fNexts,tmp+1);  
                    dheap.erase(eid_fir);
                }
                else if(dheap.getSup(eid_fir) > sup+1){
                    dheap.adjust(eid_fir,1);
                }
                else if(dheap.getSup(eid_fir) == sup){
                    ui tmp = first_sup - sup;
                    linear_heap->decrement(first_sup,eid_fir,fOff,fSup,fPres,fNexts,tmp);  
                    dheap.erase(eid_fir);
                }
            }
            else{
                if(first_sup > sup){
                    if(first_sup > sup+1)
                        dheap.push({eid_fir,first_sup - 1});
                    else
                        linear_heap->decrement(first_sup,eid_fir,fOff,fSup,fPres,fNexts);
                } 
            }


            if(dheap.find(eid_sec)){
                if(dheap.getSup(eid_sec) == sup+1){
                    ui tmp = second_sup - (sup+1);
                    linear_heap->decrement(second_sup,eid_sec,fOff,fSup,fPres,fNexts,tmp+1);  
                    dheap.erase(eid_sec);
                }
                else if(dheap.getSup(eid_sec) > sup+1){
                    dheap.adjust(eid_sec,1);
                }
                else if(dheap.getSup(eid_sec) == sup){
                    ui tmp = second_sup - sup;
                    linear_heap->decrement(second_sup,eid_sec,fOff,fSup,fPres,fNexts,tmp);  
                    dheap.erase(eid_sec);
                }
            }
            else{
                if(second_sup > sup){
                    if(second_sup > sup+1)
                        dheap.push({eid_sec,second_sup - 1});
                    else
                        linear_heap->decrement(second_sup,eid_sec,fOff,fSup,fPres,fNexts);
                } 
            }
        }

        linear_heap->empty();
        if(dheap.size > 0 && linear_heap->get_minkey() >= dheap.arr[0].sup){
            es ret = dheap.pop();
            ui _sup = linear_heap->get_key(ret.eid,fSup,fOff);
            ui tmp = _sup - ret.sup;
            if(tmp != 0)
                linear_heap->decrement(_sup,ret.eid,fOff,fSup,fPres,fNexts,tmp); 
        }

        max_size = std::max(max_size,dheap.size);  

        sup = DELETE(sup);

        eid_eid tmp_;
        fOff.fseek(eid*sizeof(eid_eid));
        fOff.fread(&tmp_,sizeof(eid_eid)); 

        fSup.fseek(tmp_.first*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t)); 
        fSup.fseek(tmp_.second*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t)); 

        gettimeofday(&end_time, NULL);
        update_interval += (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec)/1000000.0;
    }
    int current_pid = GetCurrentPid();
    float memory_usage = GetMemoryUsage(current_pid);
    memUsage = max(memory_usage,memUsage);

    total_io = total_io+fDat.get_total_io()+fEid.get_total_io()+fPres.get_total_io()+fNexts.get_total_io()+fEidToVer.get_total_io()+fSup.get_total_io()+fOff.get_total_io();
    
    fDat.fclose();
    fEid.fclose();
    fPres.fclose();
    fNexts.fclose();
    fEidToVer.fclose();
    fSup.fclose();
    fOff.fclose();
    free(nbr_u);
    free(nbr_v);
    free(eid_u);
    free(eid_v);
    free(sup_u);
    free(sup_v);

    log_debug(graphClock_.Count("heap max size: %u, intersect_time: %lf, update_time: %lf, memUsage: %f\n",max_size,insetsect_interval,update_interval, memUsage));
    for (int i=0;i<nodeNum; ++i)
		if (cntClass[i]>0)
			fout << "#edges in " << i << "-class: " << cntClass[i] << std::endl;
    return false;
}

bool Graph::inducedGraphUnOrder(uint32_t &left, uint32_t &right, uint32_t &mid, readFile &file, uint32_t &TrussEdge, uint32_t &Truss, bool delOnSubg){
    uint64_t edge_num = 0;
	uint32_t node_num = 0; 

    MyReadFile fSup( file.m_suppSort );
	fSup.fopen( BUFFERED );
    // create director in order to fill new subgraph constituted by edges whose support greater than mid
    string dir = file.m_base + "graphInfoCopy";
    file.createDir(dir);
    string dir_name = file.m_base + "graphInfoCopy/" + to_string(mid);
    file.createDir(dir_name);
    dir_name += "/";
    readFile newFile(dir_name);
    unsigned long size = 0,es = 0,eid;
    int num = 0, tmpFile = 0;
    uint32_t u,v,max_degree = 0;

    // memset(newFile.m_vertexMap,-1,sizeof(int)*newFile.m_maxID);
    string name = "sort_edge_tmp";
    string sub_dir = newFile.m_base + name;
    newFile.createDir(sub_dir);
    TEdge* edges = new TEdge[file.memEdges];
    char fileName[150];
    
    for(int i = prefix[mid-1] ; i < prefix[maxSup]; i++){
        EdgeSup tmp;
        fSup.fseek(i*sizeof(EdgeSup));
        fSup.fread(&tmp,sizeof(EdgeSup));
        u = tmp.u;
        v = tmp.v;
        eid = tmp.eid;
        newFile.moduleInSaveEdgesUnOrder(u,v,num,edges,size,eid,tmpFile,sub_dir);

    }
    total_io += fSup.get_total_io();

    fSup.fclose();

    sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpFile);
	newFile.saveTmpEdges<TEdge>(edges,size,fileName,[](const TEdge & a, const TEdge & b) {
                if(a.u < b.u)
                    return true;
                if( a.u > b.u )
                    return false;
                return a.v < b.v;
                });
    // log_debug(graphClock_.Count("Save new tmp edges done, load %ld edges.",es));
    newFile.edgeNum = num;
    int current_pid = GetCurrentPid();
    float memory_usage = GetMemoryUsage(current_pid);
    memUsage = max(memory_usage,memUsage);
	delete[] edges;

	// newFile.merge(tmpFile+1, node_num, max_degree, name,false);
    newFile.mergeByDegSort(tmpFile+1, node_num, max_degree, name,false,false);
    log_debug(graphClock_.Count("new subgraph vertex: %d, edge: %lu\n",newFile.verNum,newFile.edgeNum));

    total_io += newFile.write_io;
    
    Graph tmp_g(nodeNum,edgeNum);
    // tmp_g.InitialUnOrder(newFile);
    tmp_g.InitialUnOrderDegSort(newFile);

    /* prune optimization -- kcore */
    if(!tmp_g.peelVertex(mid,newFile,false)) 
        return false;

    tmp_g.CountTriangleSSDByDegOrder(newFile,true);
    log_debug(graphClock_.Count("mid: %d, num: %d",mid,tmp_g.prefix[tmp_g.maxSup]-tmp_g.prefix[mid-1]));
    
    if(!delOnSubg){
        bool flag_res = tmp_g.existTrussPlus(newFile,mid,TrussEdge,Truss);
        total_io += tmp_g.total_io;
        return flag_res;
    }
    else{
    
    /* prune optimization -- delete edges on original graph, so that avoid reconstructing subgraph */
    if(tmp_g.minSup >= mid)
    {
        mid = tmp_g.minSup;
        log_debug(graphClock_.Count("minSup >= mid, mid: %d",mid));
    }

    #ifdef LazyUpdate
    string path = file.m_base + "linearList";
    file.createDir(path);
    path += "/";
    ListLinearHeapTruss *linear_heap = new ListLinearHeapTruss(edgeNum,maxSup,path,newFile.m_supp);
    DynamicHeap dheap = DynamicHeap(100000000);
    if(tmp_g.existTrussLazyUpdate(newFile,mid,TrussEdge,Truss,linear_heap,dheap)){

    #else
    if(tmp_g.existTrussPlus(newFile,mid,TrussEdge,Truss)){

    #endif
        for(int u = 0; u < newFile.verNum; u++ ){
            if(!vis[u] && degree_[u] <= mid)
                vis[u] = true;
        }
        left = mid+1;
        uint32_t capacity, last_mid = mid;
        while(left <= right){
            mid = (left+right)/2;
            capacity = tmp_g.prefix[tmp_g.maxSup]-tmp_g.prefix[mid-1];
            log_debug(graphClock_.Count("mid: %u, capacity: %u, left: %u, right: %u",mid,capacity,left,right));
            if(capacity < (mid+2)*(mid+1) / 2)
            {
                right = mid-1;
                continue;
            }
            log_debug(graphClock_.Count("last_mid: %u, prefix[last_mid-1]: %u, mid: %u, prefix[mid-1]: %u",last_mid,tmp_g.prefix[last_mid-1],mid,tmp_g.prefix[mid]));
            
            #ifdef LazyUpdate
            if(!tmp_g.deleteEdgeLazyUpdateTest(file,newFile,last_mid,mid,prefix,TrussEdge,Truss,linear_heap,dheap))
            #else
            if(!tmp_g.deleteEdge(file,newFile,last_mid,mid,prefix,TrussEdge,Truss))
            #endif
            {
                #ifdef LazyUpdate
                delete linear_heap;
                #endif
                total_io += tmp_g.total_io;
                return false;
            }
            else
            {
                left = mid+1;
                for(int u = 0; u < newFile.verNum; u++ ){
                    if(!vis[u] && degree_[u] <= mid)
                        vis[u] = true;
                }
            }
            last_mid = mid;
        }
    }
    else
    {
        #ifdef LazyUpdate
            delete linear_heap;
        #endif
        total_io += tmp_g.total_io;
        return false;
    }
    #ifdef LazyUpdate
        delete linear_heap;
    #endif    
    total_io += tmp_g.total_io;
    return true;
    }
}

void Graph::updateNbrInExistTrussLazyUpdate(uint32_t u, uint32_t *nbr_u, uint64_t *eid_u,  uint32_t *sup_u,
uint32_t v, uint32_t *nbr_v, uint64_t *eid_v, uint32_t *sup_v,
MyReadFile &fDat, MyReadFile &fEid, MyReadFile &fSup,MyReadFile &fOff, MyReadFile &fPres, MyReadFile &fNexts,
uint32_t sup, DynamicHeap &dheap, ListLinearHeapTruss *linear_heap){
    std::vector<eid_eid> comm;
    IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,sup,comm,true);
    // printf("u: %u, v: %u, eid: %lu, sup: %u, comm: %u\n",u,v,eid,sup,comm.size());
    for(int i = 0; i < comm.size(); i++){
        ui first_sup = sup_u[comm[i].first];
        ui second_sup = sup_v[comm[i].second];
        uint64_t eid_fir = eid_u[comm[i].first];
        uint64_t eid_sec = eid_v[comm[i].second];

        if(dheap.find(eid_fir)){
            if(dheap.getSup(eid_fir) == sup+1){
                ui tmp = first_sup - (sup+1);
                linear_heap->decrement(first_sup,eid_fir,fOff,fSup,fPres,fNexts,tmp+1);  
                dheap.erase(eid_fir);
            }
            else if(dheap.getSup(eid_fir) > sup+1){
                dheap.adjust(eid_fir,1);
            }
            else if(dheap.getSup(eid_fir) == sup){
                ui tmp = first_sup - sup;
                linear_heap->decrement(first_sup,eid_fir,fOff,fSup,fPres,fNexts,tmp);  
                dheap.erase(eid_fir);
            }
            else{
                printf("first eid: %lu, oriSup: %u, now: %u, sup: %u\n",eid_fir,
                first_sup,dheap.getSup(eid_fir),sup);
                printf("size: %u\n",dheap.size);
            }
        }
        else{
            if(first_sup > sup){
                if(first_sup > sup+1)
                    dheap.push({eid_fir,first_sup - 1});
                else
                    linear_heap->decrement(first_sup,eid_fir,fOff,fSup,fPres,fNexts);
            } 
        }


        if(dheap.find(eid_sec)){
            if(dheap.getSup(eid_sec) == sup+1){
                ui tmp = second_sup - (sup+1);
                linear_heap->decrement(second_sup,eid_sec,fOff,fSup,fPres,fNexts,tmp+1);  
                dheap.erase(eid_sec);
            }
            else if(dheap.getSup(eid_sec) > sup+1){
                dheap.adjust(eid_sec,1);
            }
            else if(dheap.getSup(eid_sec) == sup){
                ui tmp = second_sup - sup;
                linear_heap->decrement(second_sup,eid_sec,fOff,fSup,fPres,fNexts,tmp);  
                dheap.erase(eid_sec);
            }
            else{
                printf("second eid: %lu, oriSup: %u, now: %u, sup: %u\n",eid_sec,
                linear_heap->get_key(eid_sec,fSup),dheap.getSup(eid_sec),sup);
                printf("size: %u\n",dheap.size);
            }
        }
        else{
            if(second_sup > sup){
                if(second_sup > sup+1)
                    dheap.push({eid_sec,second_sup - 1});
                else
                    linear_heap->decrement(second_sup,eid_sec,fOff,fSup,fPres,fNexts);
            } 
        }
    }


    linear_heap->empty();
    if(dheap.size > 0 && linear_heap->get_minkey() >= dheap.arr[0].sup){
        es ret = dheap.pop();
        ui _sup = linear_heap->get_key(ret.eid,fSup,fOff);
        ui tmp = _sup - ret.sup;
        if(tmp != 0)
            linear_heap->decrement(_sup,ret.eid,fOff,fSup,fPres,fNexts,tmp); 
    }

}


bool Graph::deleteEdgeLazyUpdateTest(readFile &file, readFile &newFile,uint32_t last_mid, uint32_t &mid, uint32_t *prefix_, uint32_t &TrussEdge, uint32_t &Truss, ListLinearHeapTruss *linear_heap, DynamicHeap &dheap){
    log_debug(graphClock_.Count("mid: %d,prefix[mid-1]: %u,last_mid: %d,prefix[last_mid-1]:%u",mid,prefix[mid-1],last_mid,prefix[last_mid-1]));
    MyReadFile fSupSort( newFile.m_suppSort );
	fSupSort.fopen( BUFFERED );

    MyReadFile fSup( newFile.m_supp );
	fSup.fopen( NOBUFFER );
    
    MyReadFile fPres( linear_heap->m_pres );
	fPres.fopen( NOBUFFER );
    MyReadFile fNexts( linear_heap->m_nexts );
	fNexts.fopen( NOBUFFER );
    MyReadFile fEidToVer( newFile.m_eidToVer );
	fEidToVer.fopen( BUFFERED );

    MyReadFile fDat( newFile.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fEid( newFile.m_eid );
	fEid.fopen( BUFFERED );
    MyReadFile fOff( newFile.m_offset );
	fOff.fopen( BUFFERED );

    uint32_t u,v,sup,p,record_sup;
    uint64_t eid,eid_uv;
    EdgeSup tmp;
    TEdge se,e;
    std::queue<EdgeSup> del_que;
    int max_size = 0;

    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * newFile.maxDeg);
    uint32_t* nbr_v = (uint32_t *)malloc(sizeof(uint32_t) * newFile.maxDeg);

    uint32_t* sup_u = (uint32_t *)malloc(sizeof(uint32_t) * newFile.maxDeg);
    uint32_t* sup_v = (uint32_t *)malloc(sizeof(uint32_t) * newFile.maxDeg);

    uint64_t* eid_u = (uint64_t *)malloc(sizeof(uint64_t) * newFile.maxDeg);
    uint64_t* eid_v = (uint64_t *)malloc(sizeof(uint64_t) * newFile.maxDeg);


    log_info(graphClock_.Count("before exist edge: %lu\n",linear_heap->exist_num));
    uint32_t test_c = 0;


    isInDelQue.clear();

    while(dheap.size > 0 && dheap.arr[0].sup < last_mid){
        es ret = dheap.pop();
        ui _sup = linear_heap->get_key(ret.eid,fSup,fOff);
        Edge tmp;
        fEidToVer.fseek(ret.eid*sizeof(Edge));
        fEidToVer.fread(&tmp,sizeof(Edge));
        del_que.push({tmp.u,tmp.v,_sup,ret.eid});
        isInDelQue.insert(ret.eid);
    }
    
    
    uint32_t _sup, ori_sup, dheap_sup;
    for(int i = prefix[last_mid-1] ; i < prefix[mid-1]; i++){
        fSupSort.fseek(i*sizeof(EdgeSup));
        fSupSort.fread(&tmp,sizeof(EdgeSup));
        u = tmp.u;
        v = tmp.v;
        eid = tmp.eid;
        
        if(vis[u] || vis[v])
            continue;
        
        eid_eid tmp_;
        fOff.fseek(eid*sizeof(eid_eid));
        fOff.fread(&tmp_,sizeof(eid_eid)); 
        
        fSup.fseek(tmp_.first*sizeof(uint32_t));
        fSup.fread(&sup,sizeof(uint32_t));
        dheap_sup = sup;
        if(dheap.find(eid)){
            // log_info(graphClock_.Count("**********ori_sup: %u,now_sup: %u,tmp.sup: %u, eid: %u\n",sup,dheap.getSup(eid),tmp.sup,eid));
            dheap_sup = dheap.getSup(eid);
            dheap.erase(eid);
        }


        if(MOVE(sup))
            continue;

        linear_heap->remove(eid,sup,fPres,fNexts);
        //update neighbor
        loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
        loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);
        loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
        
        loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
        loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);
        loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
        /* optimization */

        updateNbrInDeleteEdgeLazyUpdate(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v, 
        fDat,fEid,fSup,fOff,fPres,fNexts,sup,del_que,last_mid, dheap, linear_heap);
        // log_info(graphClock_.Count("[1]delete sup : %u,dheap_sup: %u, eid: %lu",sup,dheap_sup,eid));
        
        --degree_[u];
        --degree_[v];
        sup = DELETE(sup);
        fSup.fseek(tmp_.first*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t));
        fSup.fseek(tmp_.second*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t)); 
        

    }
    // log_info(graphClock_.Count("after deleting, exist edge: %lu",linear_heap->exist_num));
    // log_info(graphClock_.Count("del_que size: %u",del_que.size()));

    linear_heap->empty();
    while(dheap.size > 0 && linear_heap->get_minkey() > dheap.arr[0].sup){
        es ret = dheap.pop();
        ui _sup = linear_heap->get_key(ret.eid,fSup,fOff);
        // log_info(graphClock_.Count("dheap.arr[0].sup: %u, linear_heap sup: %u",dheap.arr[0].sup,_sup));
        ui tmp = _sup - ret.sup;
        if(tmp != 0)
            linear_heap->decrement(_sup,ret.eid,fOff,fSup,fPres,fNexts,tmp); 
    }

    while(linear_heap->exist_num){
        Edge tmp;
        linear_heap->pop_min(eid,sup,fPres,fNexts);
        if(MOVE(sup)){
            log_info(graphClock_.Count("eid: %u is deleted",eid));
            continue;            
        }

        record_sup = sup;      
        #ifndef Maintenance
        if(linear_heap->exist_num+1 < (mid+2)*(mid+1) / 2){
            mid = sup;
            log_info(graphClock_.Count("deleteEdge func error last_sup: %d, exist_edge: %d\n",sup,linear_heap->exist_num));
            total_io = total_io+fSupSort.get_total_io()+fDat.get_total_io()+fEid.get_total_io()+fPres.get_total_io()+fNexts.get_total_io()+fEidToVer.get_total_io()+fSup.get_total_io()+fOff.get_total_io();

            fSupSort.fclose();
            fSup.fclose();
            fEid.fclose();
            fDat.fclose();
            fPres.fclose();
            fNexts.fclose();
            fEidToVer.fclose();
            fOff.fclose();
            free(nbr_u);
            free(nbr_v);
            free(eid_u);
            free(eid_v);
            free(sup_u);
            free(sup_v);
            return false;
        }
        #endif   
        
        if(sup >= mid){
            Truss = sup;
            last_sup = sup;
            linear_heap->insert(eid,sup,fOff,fSup,fPres,fNexts);
            TrussEdge = linear_heap->exist_num;
            total_io = total_io+fSupSort.get_total_io()+fDat.get_total_io()+fEid.get_total_io()+fPres.get_total_io()+fNexts.get_total_io()+fEidToVer.get_total_io()+fSup.get_total_io()+fOff.get_total_io();

            log_info(graphClock_.Count("deleteEdge func success last_sup: %d, exist_num: %u, test_c: %u\n",sup,linear_heap->exist_num,test_c));
            fSupSort.fclose();
            fSup.fclose();
            fEid.fclose();
            fDat.fclose();
            fPres.fclose();
            fNexts.fclose();
            fEidToVer.fclose();
            fOff.fclose();
            free(nbr_u);
            free(nbr_v);
            free(eid_u);
            free(eid_v);
            free(sup_u);
            free(sup_v);
            return true;
        }
        fEidToVer.fseek(eid*sizeof(Edge));
        fEidToVer.fread(&tmp,sizeof(Edge));
        u = tmp.u;
        v = tmp.v;
        
        loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
        loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
        loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);
        
        loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
        loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
        loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);

        updateNbrInExistTrussLazyUpdate(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v, 
        fDat,fEid,fSup,fOff,fPres,fNexts,sup,dheap,linear_heap);

        max_size = std::max(max_size,dheap.size);  
        sup = DELETE(sup);
        eid_eid tmp_;
        fOff.fseek(eid*sizeof(eid_eid));
        fOff.fread(&tmp_,sizeof(eid_eid)); 

        fSup.fseek(tmp_.first*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t)); 
        fSup.fseek(tmp_.second*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t)); 
        --degree_[u];
        --degree_[v];
    }
    if(dheap.size > 0){
        log_info(graphClock_.Count("size: %d, minimal sup: %u",dheap.size,dheap.arr[0].sup));
    }

    total_io = total_io+fSupSort.get_total_io()+fDat.get_total_io()+fEid.get_total_io()+fPres.get_total_io()+fNexts.get_total_io()+fEidToVer.get_total_io()+fSup.get_total_io()+fOff.get_total_io();

    fSupSort.fclose();
    fSup.fclose();
    fPres.fclose();
    fNexts.fclose();
    fEidToVer.fclose();
    fDat.fclose();
    fEid.fclose();
    fOff.fclose();
    free(nbr_u);
    free(nbr_v);
    free(eid_u);
    free(eid_v);
    free(sup_u);
    free(sup_v);
    mid = record_sup;
    return false;

}

bool Graph::deleteEdge(readFile &file, readFile &newFile,uint32_t last_mid, uint32_t &mid, uint32_t *prefix_, uint32_t &TrussEdge, uint32_t &Truss){
    log_info(graphClock_.Count("mid: %d,prefix[mid-1]: %u,last_mid: %d,prefix[last_mid-1]:%u",mid,prefix[mid-1],last_mid,prefix[last_mid-1]));
    MyReadFile fSupSort( newFile.m_suppSort );
	fSupSort.fopen( BUFFERED );

    MyReadFile fSup( newFile.m_supp );
	fSup.fopen( NOBUFFER );

    MyReadFile fBinEdge( newFile.m_binEdge );
	fBinEdge.fopen( NOBUFFER );
    MyReadFile fEdgePos( newFile.m_ePos );
	fEdgePos.fopen( NOBUFFER );
    
    MyReadFile fDat( newFile.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fEid( newFile.m_eid );
	fEid.fopen( BUFFERED );

    MyReadFile fOff( newFile.m_offset );
	fOff.fopen( BUFFERED );

    uint32_t u,v,sup,p,record_sup;
    uint64_t eid,eid_uv;
    EdgeSup tmp_edge;
    TEdge se,e;

    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * newFile.maxDeg);
    uint32_t* nbr_v = (uint32_t *)malloc(sizeof(uint32_t) * newFile.maxDeg);

    uint64_t* eid_u = (uint64_t *)malloc(sizeof(uint64_t) * newFile.maxDeg);
    uint64_t* eid_v = (uint64_t *)malloc(sizeof(uint64_t) * newFile.maxDeg);
    
    uint32_t* sup_u = (uint32_t *)malloc(sizeof(uint32_t) * newFile.maxDeg);
    uint32_t* sup_v = (uint32_t *)malloc(sizeof(uint32_t) * newFile.maxDeg);

    int fixedEdgeCount = 0;
    uint32_t test_c = 0;
    for(int i = prefix[last_mid-1] ; i < prefix[mid-1]; i++){
        fSupSort.fseek(i*sizeof(EdgeSup));
        fSupSort.fread(&tmp_edge,sizeof(EdgeSup));
        u = tmp_edge.u;
        v = tmp_edge.v;
        eid = tmp_edge.eid;
        

        if(vis[u] || vis[v])
            continue;
        eid_eid tmp;
        fOff.fseek(eid*sizeof(eid_eid));
        fOff.fread(&tmp,sizeof(eid_eid));
        
        fSup.fseek(tmp.first*sizeof(uint32_t));
        fSup.fread(&sup,sizeof(uint32_t));
        
        if(MOVE(sup))
            continue;
        
        //update neighbor
        int start = bin[last_sup], end, pre_start;
        loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
        loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
        
        loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
        loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
        
        loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);
        loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);


        /* optimization */
        std::vector<eid_eid> comm;
        IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,sup,comm,true);        
        for(int i = 0; i < comm.size(); i++){
            // if( (eid_u[comm[i].first] == mark || eid_v[comm[i].second] == mark))
            //     printf("eid: %u , neighbor %u, %u, <%u,%u>\n",eid,eid_u[comm[i].first],eid_v[comm[i].second],u,v);
            updateEdgeSSDNew(u,nbr_u[comm[i].first],eid_u[comm[i].first],comm[i].first,sup,sup_u,fEdgePos,fSup,fBinEdge,fOff);
			updateEdgeSSDNew(v,nbr_u[comm[i].first],eid_v[comm[i].second],comm[i].second,sup,sup_v,fEdgePos,fSup,fBinEdge,fOff);
        }
        end = bin[last_sup];
        if(sup >= mid)
            fixedEdgeCount++;
        --degree_[u];
        --degree_[v];

        // log_debug(graphClock_.Count("eid: %d,sup:%u",eid,sup));

        sup = DELETE(sup);
        fSup.fseek(tmp.first*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t));

        fSup.fseek(tmp.second*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t));

        while(end - start > 0)
        {
            pre_start = bin[last_sup];
            for(int s = start; s < end; s++){
                fBinEdge.fseek(s*sizeof(TEdge));
                fBinEdge.fread(&e,sizeof(TEdge));
                u = e.u;
                v = e.v;
                eid_uv = e.eid; 

                fOff.fseek(eid_uv*sizeof(eid_eid));
                fOff.fread(&tmp,sizeof(eid_eid));
                
                fSup.fseek(tmp.first*sizeof(uint32_t));
                fSup.fread(&sup,sizeof(uint32_t));

                
                loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
                loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
                
                loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
                loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
                
                loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);
                loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);

                std::vector<eid_eid> comm;
                IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,sup,comm,true);        
                for(int i = 0; i < comm.size(); i++){
                    updateEdgeSSDNew(u,nbr_u[comm[i].first],eid_u[comm[i].first],comm[i].first,sup,sup_u,fEdgePos,fSup,fBinEdge,fOff);
                    updateEdgeSSDNew(v,nbr_u[comm[i].first],eid_v[comm[i].second],comm[i].second,sup,sup_v,fEdgePos,fSup,fBinEdge,fOff);
                }
                log_debug(graphClock_.Count("eid: %d,sup:%u in while",eid,sup));
                sup = DELETE(sup);
                test_c++;
                --degree_[u];
                --degree_[v];

                fSup.fseek(tmp.first*sizeof(uint32_t));
                fSup.fwrite(&sup,sizeof(uint32_t));

                fSup.fseek(tmp.second*sizeof(uint32_t));
                fSup.fwrite(&sup,sizeof(uint32_t));
            }

            end = bin[last_sup];
            start = pre_start;
        }
    }

    for(int s = bin[last_sup]; s < newFile.edgeNum; s++){
        fBinEdge.fseek(s*sizeof(TEdge));
        fBinEdge.fread(&e,sizeof(TEdge));
        u = e.u;
        v = e.v;
        eid_uv = e.eid; 

        eid_eid tmp;
        fOff.fseek(eid_uv*sizeof(eid_eid));
        fOff.fread(&tmp,sizeof(eid_eid));
        
        fSup.fseek(tmp.first*sizeof(uint32_t));
        fSup.fread(&sup,sizeof(uint32_t));
        
        if(MOVE(sup))
            continue;
        record_sup = sup;
        if(newFile.edgeNum-s-fixedEdgeCount < (mid+1)*(mid+2)/2){
            mid = sup;
            log_debug(graphClock_.Count("fail in deleteEdge mid: %d, s: %d\n",mid,s));
            total_io = total_io + fSupSort.get_total_io() + fBinEdge.get_total_io() + fEdgePos.get_total_io() + fDat.get_total_io() + fEid.get_total_io() + fSup.get_total_io() + fOff.get_total_io(); 

            fSupSort.fclose();
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
            return false;
        }

        if(sup >= mid){
            TrussEdge = newFile.edgeNum-s-fixedEdgeCount;
            Truss = mid;
            last_s = s;
            last_sup = sup;
            log_debug(graphClock_.Count("success in deleteEdge TrussEdge: %d,last_sup: %d, s: %d\n",TrussEdge,last_sup,s));
            total_io = total_io + fSupSort.get_total_io() + fBinEdge.get_total_io() + fEdgePos.get_total_io() + fDat.get_total_io() + fEid.get_total_io() + fSup.get_total_io() + fOff.get_total_io(); 

            fSupSort.fclose();
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
            return true;
        }
        
        loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
        loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
        
        loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
        loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);

        loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);
        loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);
        

        std::vector<eid_eid> comm;
        IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,sup,comm,true);        
        for(int i = 0; i < comm.size(); i++){
            updateEdgeSSDNew(u,nbr_u[comm[i].first],eid_u[comm[i].first],comm[i].first,sup,sup_u,fEdgePos,fSup,fBinEdge,fOff);
			updateEdgeSSDNew(v,nbr_u[comm[i].first],eid_v[comm[i].second],comm[i].second,sup,sup_v,fEdgePos,fSup,fBinEdge,fOff);
        }
        --degree_[u];
        --degree_[v];
        // log_debug(graphClock_.Count("delete in for eid: %d,sup:%u",eid,sup));

        sup = DELETE(sup);

        fSup.fseek(tmp.first*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t));
        fSup.fseek(tmp.second*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t));

    }
    log_debug(graphClock_.Count("error record_sup in deleteEdge function: %d\n",record_sup));
    free(nbr_u);
    free(nbr_v);
    free(eid_u);
    free(eid_v);
    free(sup_u);
    free(sup_v);
    total_io = total_io + fSupSort.get_total_io() + fBinEdge.get_total_io() + fEdgePos.get_total_io() + fDat.get_total_io() + fEid.get_total_io() + fSup.get_total_io() + fOff.get_total_io(); 

    fSupSort.fclose();
    fSup.fclose();
    fBinEdge.fclose();
    fEdgePos.fclose();
    fDat.fclose();
    fEid.fclose();
    fOff.fclose();
    mid = record_sup;
    return false;
}


void Graph::trussDecomLazyUpdate(readFile &file){
    log_info(graphClock_.Count("trussDecomLazyUpdate func begin"));
	MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fEid( file.m_eid );
	fEid.fopen( BUFFERED );
    MyReadFile fOff( file.m_offset );
	fOff.fopen( BUFFERED );
    struct timeval start_time, end_time;
    double insetsect_interval = 0;
    double update_interval = 0;
    
    //write the two vertex of edge into fEidToV 
    // writeTwoVerFromEid(file);

    std::string path = file.m_base + "linearList/";
    file.createDir(path);
    ListLinearHeapTruss *linear_heap = new ListLinearHeapTruss(edgeNum,maxSup,path,file.m_supp);
    /* linear_heap->init(edgeNum,maxSup); */

    // linear_heap->init(edgeNum,fOff);

    MyReadFile fSupSort( file.m_suppSort );
	fSupSort.fopen( BUFFERED );
    MyReadFile fSup( file.m_supp );
	fSup.fopen( NOBUFFER );

    FILE* filePres = fopen(linear_heap->m_pres.c_str(),"wb");
	FILE* fileNexts = fopen(linear_heap->m_nexts.c_str(),"wb");
    FILE* fileEidToV = fopen(file.m_eidToVer.c_str(),"wb");

    EdgeSup tmp_;
    Edge e;
    uint64_t fileEdgeNum = 0;
    uint32_t sup;
    memset(vis,0,sizeof(bool)*nodeNum);
    for(int i = 0; i < file.edgeNum; i++){
        fSupSort.fseek(i*sizeof(EdgeSup));
        fSupSort.fread(&tmp_,sizeof(EdgeSup));
        
        e = {tmp_.u,tmp_.v};

        #ifdef Maintenance
        uint32_t x = min(tmp_.u,tmp_.v);
        uint32_t y = max(tmp_.u,tmp_.v);
        uint64_t _l = COMBINE(x,y);
        eid_eid eid_;
        
        fOff.fseek(tmp_.eid*sizeof(eid_eid));
        fOff.fread(&eid_,sizeof(eid_eid)); 

        fSup.fseek(eid_.first*sizeof(uint32_t));
        fSup.fread(&sup,sizeof(uint32_t)); 
        if(isInDelQue.find(_l) != isInDelQue.end())   {
            continue;
        }
        
        if(MOVE(sup))
            printf("xx\n");
        tmp_.sup = sup;
        #endif
        // if(tmp_.eid > edgeNum)
        // {
        //     printf("u: %u, v: %u, i: %d, sup: %u, eid: %lu\n",tmp_.u,tmp_.v,i,tmp_.sup,tmp_.eid);
        // }
        fseek(fileEidToV,tmp_.eid*sizeof(Edge),SEEK_SET);
        fwrite(&e,sizeof(Edge),1,fileEidToV);

        linear_heap->insert(tmp_.eid,tmp_.sup,filePres,fileNexts);
        fileEdgeNum++;
    }
    file.edgeNum = fileEdgeNum;
    total_io += file.edgeNum;
    total_io += fSupSort.get_total_io();
    total_io += linear_heap->total_io;
    linear_heap->total_io = 0;

    fclose(filePres);
    fclose(fileNexts);
    fclose(fileEidToV);
    fSupSort.fclose();

    log_info(graphClock_.Count("finish constructing linearList"));

    MyReadFile fPres( linear_heap->m_pres );
	fPres.fopen( NOBUFFER );
    MyReadFile fNexts( linear_heap->m_nexts );
	fNexts.fopen( NOBUFFER );
    MyReadFile fEidToVer( file.m_eidToVer );
	fEidToVer.fopen( BUFFERED );
    
    DynamicHeap dheap = DynamicHeap(1000000000);
    ui u,v;
    last_sup = ui(-1);
    int max_size = 0;
    Edge tmp;

    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* nbr_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    uint32_t* sup_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* sup_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    uint64_t* eid_u = (uint64_t *)malloc(sizeof(uint64_t) * file.maxDeg);
    uint64_t* eid_v = (uint64_t *)malloc(sizeof(uint64_t) * file.maxDeg);  

    for(ui64 i = 0; i < file.edgeNum; i++){    
        gettimeofday(&start_time, NULL);
        ui64 eid;
        linear_heap->pop_min(eid,sup,fPres,fNexts);
        if(sup == 0)continue;
        fEidToVer.fseek(eid*sizeof(Edge));
        fEidToVer.fread(&tmp,sizeof(Edge));
        u = tmp.u;
        v = tmp.v;
        printClass(u,v,sup+2);

        if(sup != last_sup)
        {
            // log_info(graphClock_.Count("i: %u, u: %u, v: %u, eid: %u, sup: %u, size: %d",i,u,v,eid,sup,dheap.size));
            last_sup = sup;
            // linear_heap->print();
        }
        loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
        loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
        loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);

        
        loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
        loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
        loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);

        std::vector<eid_eid> comm;
        IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,sup,comm,true);
        
        gettimeofday(&end_time, NULL);
        insetsect_interval += (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec)/1000000.0;
        gettimeofday(&start_time, NULL);

        
        for(int i = 0; i < comm.size(); i++){
            ui first_sup = sup_u[comm[i].first];
            ui second_sup = sup_v[comm[i].second];
            uint64_t eid_fir = eid_u[comm[i].first];
            uint64_t eid_sec = eid_v[comm[i].second];

            if(dheap.find(eid_fir)){
                if(dheap.getSup(eid_fir) == sup+1){
                    ui tmp = first_sup - (sup+1);
                    linear_heap->decrement(first_sup,eid_fir,fOff,fSup,fPres,fNexts,tmp+1);  
                    dheap.erase(eid_fir);
                }
                else if(dheap.getSup(eid_fir) > sup+1){
                    dheap.adjust(eid_fir,1);
                }
                else if(dheap.getSup(eid_fir) == sup){
                    ui tmp = first_sup - sup;
                    linear_heap->decrement(first_sup,eid_fir,fOff,fSup,fPres,fNexts,tmp);  
                    dheap.erase(eid_fir);
                }
                else{
                    printf("first eid: %lu, oriSup: %u, now: %u, sup: %u\n",eid_fir,
                    first_sup,dheap.getSup(eid_fir),sup);
                    printf("size: %u\n",dheap.size);
                }
            }
            else{
                if(first_sup > sup){
                    if(first_sup > sup+1)
                        dheap.push({eid_fir,first_sup - 1});
                    else
                        linear_heap->decrement(first_sup,eid_fir,fOff,fSup,fPres,fNexts);
                } 
            }


            if(dheap.find(eid_sec)){
                if(dheap.getSup(eid_sec) == sup+1){
                    ui tmp = second_sup - (sup+1);
                    linear_heap->decrement(second_sup,eid_sec,fOff,fSup,fPres,fNexts,tmp+1);  
                    dheap.erase(eid_sec);
                }
                else if(dheap.getSup(eid_sec) > sup+1){
                    dheap.adjust(eid_sec,1);
                }
                else if(dheap.getSup(eid_sec) == sup){
                    ui tmp = second_sup - sup;
                    linear_heap->decrement(second_sup,eid_sec,fOff,fSup,fPres,fNexts,tmp);  
                    dheap.erase(eid_sec);
                }
                else{
                    printf("second eid: %lu, oriSup: %u, now: %u, sup: %u\n",eid_sec,
                    linear_heap->get_key(eid_sec,fSup),dheap.getSup(eid_sec),sup);
                    printf("size: %u\n",dheap.size);
                }
            }
            else{
                if(second_sup > sup){
                    if(second_sup > sup+1)
                        dheap.push({eid_sec,second_sup - 1});
                    else
                        linear_heap->decrement(second_sup,eid_sec,fOff,fSup,fPres,fNexts);
                } 
            }
        }
        linear_heap->empty();
        if(dheap.size > 0 && linear_heap->get_minkey() >= dheap.arr[0].sup){
            es ret = dheap.pop();
            ui _sup = linear_heap->get_key(ret.eid,fSup,fOff);
            ui tmp = _sup - ret.sup;
            if(tmp != 0)
                linear_heap->decrement(_sup,ret.eid,fOff,fSup,fPres,fNexts,tmp); 
        }
        max_size = std::max(max_size,dheap.size); 
 
        
        sup = DELETE(sup);
        eid_eid tmp_;
        fOff.fseek(eid*sizeof(eid_eid));
        fOff.fread(&tmp_,sizeof(eid_eid)); 

        fSup.fseek(tmp_.first*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t)); 
        fSup.fseek(tmp_.second*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t)); 

        gettimeofday(&end_time, NULL);
        update_interval += (end_time.tv_sec - start_time.tv_sec) + (end_time.tv_usec - start_time.tv_usec)/1000000.0;
    }
    total_io = total_io + fDat.get_total_io() + fEid.get_total_io() + fPres.get_total_io() + fNexts.get_total_io() + fEidToVer.get_total_io() + fSup.get_total_io() + fOff.get_total_io(); 

    fDat.fclose();
    fEid.fclose();
    fPres.fclose();
    fNexts.fclose();
    fEidToVer.fclose();
    fSup.fclose();
    fOff.fclose();
    delete linear_heap;
    free(sup_u);
    free(sup_v);
    free(nbr_u);
    free(nbr_v);
    free(eid_u);
    free(eid_v);
    log_info(graphClock_.Count("heap max size: %u, intersect_time: %lf, update_time: %lf",max_size,insetsect_interval,update_interval));
    for (int i=0;i<nodeNum; ++i)
		if (cntClass[i]>0)
			fout << "#edges in " << i << "-class: " << cntClass[i] << std::endl;
    log_info(graphClock_.Count("trussDecomLinearList func finish"));
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
        }
    }

    fclose(fEidToV);
    fDat.fclose();
    fEid.fclose();
    free(nbr_u);
    log_info(graphClock_.Count("finish writing eid to vertex"));
}


void Graph::trussDecomLinearList(readFile &file){
    log_info(graphClock_.Count("trussDecomLinearList func begin"));
	MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fEid( file.m_eid );
	fEid.fopen( BUFFERED );
    MyReadFile fOff( file.m_offset );
	fOff.fopen( BUFFERED );
    writeTwoVerFromEid(file);

    std::string path = "/home/jjq/research/ktruss/TrussProject/max_ktruss/linearList/";
    ListLinearHeapTruss *linear_heap = new ListLinearHeapTruss(edgeNum,maxSup,path,file.m_supp);
    // linear_heap->init(edgeNum,maxSup);
    linear_heap->init(edgeNum,fOff);
    MyReadFile fPres( linear_heap->m_pres );
	fPres.fopen( NOBUFFER );
    MyReadFile fNexts( linear_heap->m_nexts );
	fNexts.fopen( NOBUFFER );
    MyReadFile fEidToVer( file.m_eidToVer );
	fEidToVer.fopen( BUFFERED );
    MyReadFile fSup( file.m_supp );
	fSup.fopen( NOBUFFER );

    ui u,v,sup;
    Edge tmp;

    uint32_t* nbr_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* nbr_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    uint32_t* sup_u = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);
    uint32_t* sup_v = (uint32_t *)malloc(sizeof(uint32_t) * file.maxDeg);

    uint64_t* eid_u = (uint64_t *)malloc(sizeof(uint64_t) * file.maxDeg);
    uint64_t* eid_v = (uint64_t *)malloc(sizeof(uint64_t) * file.maxDeg);    
    for(ui64 i = 0; i < edgeNum; i++){
        ui64 eid;
        linear_heap->pop_min(eid,sup,fPres,fNexts);
        if(sup == 0)continue;
        fEidToVer.fseek(eid*sizeof(Edge));
        fEidToVer.fread(&tmp,sizeof(Edge));
        u = tmp.u;
        v = tmp.v;
        printClass(u,v,sup+2);

        if(sup != last_sup)
        {
            printf("i: %u, u: %u, v: %u, eid: %u, sup: %u\n",i,u,v,eid,sup);
            last_sup = sup;
            // linear_heap->print();
        }
        loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
        loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
        loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);

        loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
        loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
        loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);

        std::vector<eid_eid> comm;
        IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,sup,comm,true);

        for(int i = 0; i < comm.size(); i++){
            ui first_sup = sup_u[comm[i].first];
            ui second_sup = sup_v[comm[i].second];
            if(first_sup > sup){
                linear_heap->decrement(first_sup,eid_u[comm[i].first],fOff,fSup,fPres,fNexts);
            }
            if(second_sup > sup){
                linear_heap->decrement(second_sup,eid_v[comm[i].second],fOff,fSup,fPres,fNexts);
            }                
        }

        sup = DELETE(sup);
        eid_eid tmp_;
        fOff.fseek(eid*sizeof(eid_eid));
        fOff.fread(&tmp_,sizeof(eid_eid));
        fSup.fseek(tmp_.first*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t));    
        fSup.fseek(tmp_.second*sizeof(uint32_t));
        fSup.fwrite(&sup,sizeof(uint32_t));  
    }


    fDat.fclose();
    fEid.fclose();
    fPres.fclose();
    fNexts.fclose();
    fEidToVer.fclose();
    fSup.fclose();
    fOff.fclose();
    free(nbr_u);
    free(nbr_v);
    free(eid_u);
    free(eid_v);
    free(sup_u);
    free(sup_v);
    for (int i=0;i<nodeNum; ++i)
		if (cntClass[i]>0)
			fout << "#edges in " << i << "-class: " << cntClass[i] << std::endl;
    log_info(graphClock_.Count("trussDecomLinearList func finish"));
}


/*
dynamic graph maintainence  
*/



//maintenance: deletion first, insertion second
void Graph::dynamicMaxTrussMaintenance(readFile &file){
    total_io = 0;
    m_dynamicDel = new unordered_set<int>*[nodeNum];
	memset(m_dynamicDel, 0, sizeof(unordered_set<int>*)*nodeNum);

    m_delBit = new bool[nodeNum]();
    memcpy(degree_, degree, sizeof(uint32_t)*nodeNum);

    m_insBit = new bool[nodeNum]();

    uint32_t *coreNum = new uint32_t[nodeNum]();   
    bool *isInSubG = new bool[nodeNum](); 
    KCore(file,coreNum,false);  //here should be improved 
    

    string dir = file.m_base + "subGraphInfo/";
    readFile subG(dir);
    uint32_t subG_nodeNum = 0, degreeTmp = 0;
    uint64_t subG_edgeNum = 0; 
    MyReadFile fSubGInfo( subG.m_info );
	fSubGInfo.fopen( BUFFERED );
	fSubGInfo.fread(&subG_nodeNum,sizeof(uint32_t));
    fSubGInfo.fread(&degreeTmp,sizeof(uint32_t));
    fSubGInfo.fread(&subG_edgeNum,sizeof(uint64_t));
    fSubGInfo.fclose();
    total_io += 1;

    maxKtrussEdge = subG_edgeNum;
    MyReadFile fSubToGlobal( file.m_base+"graph.subToGlobal" );
	fSubToGlobal.fopen( BUFFERED );
    for(int i = 0; i < subG_nodeNum; i++){
        Edge e;
        fSubToGlobal.fread(&e,sizeof(Edge));
        subToGlobal[e.u] = e.v;
    }
    fSubToGlobal.fclose();
    MyReadFile fFinal( file.m_base+"graph.final" );
	fFinal.fopen( BUFFERED );
    fFinal.fread(&last_sup,sizeof(uint32_t));
    fFinal.fread(&maxKtruss,sizeof(uint32_t));
    fFinal.fclose();

    for(int k = 0; k < nodeNum; k++)
        if(coreNum[k] >= maxKtruss+1)
            isInSubG[k] = true;


    unordered_map<uint64_t,uint64_t> delTwoVerMapEid;
    unordered_map<uint32_t, uint32_t> globalToSub;
    for(auto it = subToGlobal.begin(); it != subToGlobal.end(); it++)
        globalToSub[it->second] = it->first;

    uint32_t generate_edge = 1000;
    Edge *dynamic = new Edge[generate_edge];


    dir = file.m_base + "maxTrussGraph/";
    file.createDir(dir);
    readFile newFile(dir);
    reconMaxTrussGraph(file,newFile);
    

    bool *isInMaxTruss = new bool[subG_nodeNum]();
    for(int i = 0; i < newFile.vertexId.size(); i++)
        isInMaxTruss[newFile.vertexId[i]] = true;  //subgraph
    
    printf("size: %d, subG_node: %u, subG_edge: %lu, nodeNum: %u\n",newFile.vertexId.size(),subG_nodeNum,subG_edgeNum,nodeNum);
    generateRandomEdges(file,generate_edge,dynamic,isInMaxTruss,globalToSub); 


    // first deletion

    Graph *tmp_g;
    tmp_g = new Graph(subG_nodeNum,subG_edgeNum);
    tmp_g->InitialUnOrderDegSort(newFile);
    tmp_g->CountTriangleSSDByDegOrder(newFile,true);

    uint64_t ori_newFile_edge = newFile.edgeNum, last_edge = newFile.edgeNum, ori_edgeNum = newFile.edgeNum;
    last_s = 0;

    for(int i = 0; i < generate_edge; i++){
        uint32_t a = dynamic[i].u, b = dynamic[i].v;
        if (!m_delBit[a]) {
            m_dynamicDel[a] = new unordered_set<int>;
            m_delBit[a] = true;
        }
        if (!m_delBit[b]) {
            m_dynamicDel[b] = new unordered_set<int>;
            m_delBit[b] = true;
        }
        if(m_dynamicDel[a]->find(b) != m_dynamicDel[a]->end()) continue;
        degree_[a]--;    // same operation in existTrussLazyUpdate func   
        m_dynamicDel[a]->insert(b);
        degree_[b]--;
        m_dynamicDel[b]->insert(a);
        
        if(globalToSub.find(a) != globalToSub.end() && globalToSub.find(b) != globalToSub.end() && isInMaxTruss[globalToSub[a]] && isInMaxTruss[globalToSub[b]] ){
            Edge e = {globalToSub[a],globalToSub[b]};
            ori_newFile_edge = newFile.edgeNum;
            uint64_t newFileEdge = newFile.edgeNum;
            tmp_g->delEdgeDynamic(e,newFile,maxKtruss,newFileEdge,subToGlobal,delTwoVerMapEid);
            log_debug(graphClock_.Count("i: %u, maxKtruss: %u, ori_newFile_edge: %u, newFileEdgeNum: %lu, u: %u, v: %u",i,maxKtruss,newFileEdge,newFile.edgeNum,a,b));
            if(newFileEdge == 0 || newFileEdge < (maxKtruss+1)*(maxKtruss+2)/2){
                maxKtruss -= 1;
                tmp_g->maxKtruss = maxKtruss;
                log_debug(graphClock_.Count("newFile.edgeNum == 0, maxKtruss: %u, origin_edgeNum: %lu",maxKtruss,ori_newFile_edge));                
                if(maxKtruss < last_sup)
                {
                    bool isSame = true;
                    KCore(file,coreNum,true);  //here should be improved 
                    bool *isInSubGCopy = new bool[nodeNum]();
                    int debug_c = 0;
                    for(int k = 0; k < nodeNum; k++)
                        if(coreNum[k] >= maxKtruss+1)
                            isInSubGCopy[k] = true,debug_c++;
                    for(int k = 0; k < nodeNum; k++)
                        if(isInSubG[k] != isInSubGCopy[k]){
                            isSame = false;
                            isInSubG[k] = isInSubGCopy[k];
                        }
                    delete[] isInSubGCopy;

                    if(isSame){
                        continue;
                    }
                    if(debug_c == 0){
                        debug_c = 0;
                        for(int k = 0; k < nodeNum; k++)
                            if(coreNum[k] >= maxCore)
                                isInSubG[k] = true,debug_c++;
                        // printf("after debug_c: %u\n",debug_c);
                    }
                    
                    total_io += tmp_g->total_io;
                    delete tmp_g;
                    dir = file.m_base + "subMainGraphInfo/";
                    file.createDir(dir);
                    readFile subFile(dir);
                    filterGlobalIntoSub(isInSubG,file,subFile,true);   // there is bug, because the variable in subToGlobal will be changed   

                    Graph subG(subFile.verNum,subFile.edgeNum);
                    subG_edgeNum = subFile.edgeNum, subG_nodeNum = subFile.verNum;
                    log_debug(graphClock_.Count("subFile.verNum: %u, subFile.edgeNum: %lu",subFile.verNum,subFile.edgeNum));
                    subG.Initial(subFile);
                    subG.CountTriangleSSDByDegOrder(subFile,true);
                    subG.trussDecomLazyUpdate(subFile);
                    total_io += subG.total_io;
                    maxKtruss = subG.last_sup;
                    maxKtrussEdge = subFile.edgeNum;
                    isDeleteEdgeByFunc = true;
                    reconMaxTrussGraph(subFile,newFile,false);
                    ori_edgeNum = newFile.edgeNum;
                    tmp_g = new Graph(subG_nodeNum,subG_edgeNum);
                    tmp_g->InitialUnOrderDegSort(newFile);
                    tmp_g->CountTriangleSSDByDegOrder(newFile,true);

                    memset(isInMaxTruss,0,sizeof(bool)*nodeNum);
                    for(int i = 0; i < newFile.vertexId.size(); i++)
                        isInMaxTruss[newFile.vertexId[i]] = true;

                    globalToSub.clear();
                    for(auto it = subToGlobal.begin(); it != subToGlobal.end(); it++)
                        globalToSub[it->second] = it->first;
                }
                else{
                    total_io += tmp_g->total_io;
                    delete tmp_g;   
                    reconMaxTrussGraph(file,newFile);
                    tmp_g = new Graph(subG_nodeNum,subG_edgeNum);
                    tmp_g->InitialUnOrderDegSort(newFile);
                    tmp_g->CountTriangleSSDByDegOrder(newFile,true);
                    memset(isInMaxTruss,0,sizeof(bool)*nodeNum);
                    for(int i = 0; i < newFile.vertexId.size(); i++)
                        isInMaxTruss[newFile.vertexId[i]] = true;
                    ori_newFile_edge = newFile.edgeNum;
                }
            }
        }           
    }
    total_io += tmp_g->total_io;
    log_debug(graphClock_.Count("total_io: %u, maxKtruss: %u",total_io, maxKtruss));

    // insertion: second method based sup of each edge

    total_io = 0;
    tmp_g->total_io = 0;
    tmp_g->maxKtruss = maxKtruss;

    unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> sup_dram;
    unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> eid_dram;

    uint32_t *nbr_u = new uint32_t[file.maxDeg]();
    uint32_t *nbr_v = new uint32_t[file.maxDeg]();
    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );


    for(int i = generate_edge-1; i >=0; i--){
        uint32_t a = dynamic[i].u, b = dynamic[i].v;
        uint32_t oriCoreU = coreNum[a], oriCoreV = coreNum[b];
        if (m_delBit[a]){
            m_dynamicDel[a]->erase(b);
            if(m_dynamicDel[a]->size() == 0) {
                m_delBit[a] = false;
                delete m_dynamicDel[a];
            }
        }
        if (m_delBit[b]){
            m_dynamicDel[b]->erase(a);
            if(m_dynamicDel[b]->size() == 0) {
                m_delBit[b] = false;
                delete m_dynamicDel[b];
            }
        }
        degree_[a]++,degree_[b]++;
        if(degree_[a] >= last_sup && degree_[b] >= last_sup)
            UpdateCoreDynamic(file,coreNum,a,b);
        eid_dram[a][b] = subG_edgeNum++;
        if(globalToSub.find(a) != globalToSub.end() && globalToSub.find(b) != globalToSub.end() && isInMaxTruss[globalToSub[a]] && isInMaxTruss[globalToSub[b]]){
            log_debug(graphClock_.Count("a: %u, degree[a]: %u,globalToSub[a]: %u, b: %u, degree[b]: %u,globalToSub[b]: %u",a,degree[a],globalToSub[a],b,degree[b],globalToSub[b]));
            eid_dram[a][b] = delTwoVerMapEid[COMBINE(a,b)];
            Edge e = {globalToSub[a],globalToSub[b]};
            tmp_g->insEdgeDynamicNew(e,newFile,sup_dram,eid_dram,subToGlobal);
            maxKtruss = tmp_g->maxKtruss;
        }
        else if(coreNum[a] >= last_sup+1 && coreNum[b] >= last_sup+1){
            MyReadFile fSubDat( newFile.m_dat );
	        fSubDat.fopen( BUFFERED );
            vector<uint32_t> vec_u, vec_v;
            if(globalToSub.find(a) == globalToSub.end())
            {
                loadInfo(nbr_u,degree[a],edgeListBegPtr[a],fDat);
                for(int j = 0; j < degree[a]; j++){
                    if(coreNum[nbr_u[j]] < last_sup+1) continue;
                    vec_u.push_back(nbr_u[j]);
                }
            }
            else{
                uint32_t u_degree = tmp_g->degree[globalToSub[a]];
                loadNbrForDynamic(globalToSub[a],nbr_u,u_degree,tmp_g->edgeListBegPtr[globalToSub[a]],fSubDat,true);
                for(int j = 0; j < u_degree; j++){
                    vec_u.push_back(subToGlobal[nbr_u[j]]);
                }
            }

            if(globalToSub.find(b) == globalToSub.end())
            {
                loadInfo(nbr_u,degree[b],edgeListBegPtr[b],fDat);
                for(int j = 0; j < degree[b]; j++){
                    if(coreNum[nbr_u[j]] < last_sup+1) continue;
                    vec_u.push_back(nbr_u[j]);
                }
            }
            else{
                uint32_t u_degree = tmp_g->degree[globalToSub[b]];
                loadNbrForDynamic(globalToSub[b],nbr_u,u_degree,tmp_g->edgeListBegPtr[globalToSub[b]],fSubDat,true);
                for(int j = 0; j < u_degree; j++){
                    vec_v.push_back(subToGlobal[nbr_u[j]]);
                }
            }
            uint32_t u_pos = 0, v_pos = 0, count = 0;
            while(u_pos < vec_u.size() && v_pos < vec_v.size()){
                if(vec_u[u_pos] < vec_v[v_pos]) u_pos++;
                else if(vec_u[u_pos] > vec_v[v_pos]) v_pos++;
                else{
                    count++;
                    u_pos++, v_pos++;
                }
            }
            total_io = total_io + fSubDat.get_total_io() + tmp_g->total_io;
            fSubDat.fclose();
            if(count < last_sup) continue;
            delete tmp_g;

            memset(isInSubG, false, sizeof(bool)*nodeNum);
            for(int k = 0; k < nodeNum; k++)
                if(coreNum[k] >= last_sup+1)
                    isInSubG[k] = true;

            readFile subFile(dir);
            filterGlobalIntoSub(isInSubG,file,subFile,true); 
            globalToSub.clear();
            for(auto it = subToGlobal.begin(); it != subToGlobal.end(); it++)
                globalToSub[it->second] = it->first;

            Graph subG(subFile.verNum,subFile.edgeNum);
            subG_edgeNum = subFile.edgeNum, subG_nodeNum = subFile.verNum;
            printf("subFile.verNum: %u, subFile.edgeNum: %lu\n",subFile.verNum,subFile.edgeNum);
            subG.Initial(subFile);
            subG.CountTriangleSSDByDegOrder(subFile,true);
            subG.trussDecomLazyUpdate(subFile);
            last_sup = subG.last_sup;
            maxKtruss = subG.last_sup;
            maxKtrussEdge = subFile.edgeNum;

            reconMaxTrussGraph(subFile,newFile,false);
            tmp_g = new Graph(subG_nodeNum,subG_edgeNum);
            tmp_g->InitialUnOrderDegSort(newFile);
            tmp_g->CountTriangleSSDByDegOrder(newFile,true);

            total_io += subG.total_io;
            dynamic[i].u = globalToSub[dynamic[i].u], dynamic[i].v = globalToSub[dynamic[i].v];
            printf("core: ori_u %u, now %u, ori_v %u, now %u\n",oriCoreU,coreNum[a],oriCoreV,coreNum[b]);
            
            globalToSub.clear();
            for(auto it = subToGlobal.begin(); it != subToGlobal.end(); it++)
                globalToSub[it->second] = it->first;
            
            Edge e = {globalToSub[a],globalToSub[b]};
            tmp_g->insEdgeDynamicNew(e,newFile,sup_dram,eid_dram,subToGlobal);
        }

    }
    total_io += tmp_g->total_io;
    log_debug(graphClock_.Count("total_io: %lu, maxKtruss: %u",total_io,maxKtruss));


	if (m_dynamicDel != NULL) {
		for (int i = 0;i<nodeNum;++i) {
			if (m_delBit[i]) {
				delete m_dynamicDel[i];
			}
		}
		delete[] m_dynamicDel;
	}

    delete[] dynamic;
    delete[] m_delBit;
    delete[] m_insBit;
    delete[] coreNum;
    delete[] isInSubG;
    delete[] isInMaxTruss;
    delete[] nbr_u;
    if(tmp_g)
        delete tmp_g;
    fDat.fclose();
}


void Graph::delEdgeDynamic(Edge del_e, readFile &newFile, uint32_t maxK, uint64_t &newFileEdge,
 unordered_map<uint32_t, uint32_t>& subToGlobal, unordered_map<uint64_t, uint64_t>& delTwoVerMapEid){
    //judge  whether vertex u and v are both in maxKtruss graph
    uint32_t u = del_e.u, v = del_e.v;
    vector<uint32_t>::iterator iter_u, iter_v;
    
    iter_u = find(newFile.vertexId.begin(), newFile.vertexId.end(), u);
    iter_v = find(newFile.vertexId.begin(), newFile.vertexId.end(), v);
    if(iter_u == newFile.vertexId.end() || iter_v == newFile.vertexId.end())
        return;

    printf("<%u,%u> is in maxTruss: %u, u'deg: %u, v'deg: %u, subToGlobal[u]: %u, subToGlobal[v]: %u\n",u,v,maxK,degree[u],degree[v],subToGlobal[u],subToGlobal[v]);
    
    MyReadFile fDat(newFile.m_dat);
    fDat.fopen( BUFFERED );
    MyReadFile fEid(newFile.m_eid);
    fEid.fopen( BUFFERED );
    MyReadFile fOff(newFile.m_offset);
    fOff.fopen( BUFFERED );
    MyReadFile fSup(newFile.m_supp);
    fSup.fopen( NOBUFFER );
    MyReadFile fSupSort(newFile.m_suppSort);
    fSupSort.fopen( NOBUFFER );

    uint32_t* nbr_u = new uint32_t[newFile.maxDeg];
    uint32_t* nbr_v = new uint32_t[newFile.maxDeg];

    uint32_t* sup_u = new uint32_t[newFile.maxDeg];
    uint32_t* sup_v = new uint32_t[newFile.maxDeg];

    uint64_t* eid_u = new uint64_t[newFile.maxDeg];
    uint64_t* eid_v = new uint64_t[newFile.maxDeg];

    loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
    loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);
    loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
    
    loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
    loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);
    loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
    
    uint64_t eid_uv;
    uint32_t sup_uv, index;
    for(int i = 0; i < degree[u]; i++)
    {
        if(nbr_u[i] == v)
        {
            eid_uv = eid_u[i], sup_uv = sup_u[i], index = i;
            delTwoVerMapEid[COMBINE(subToGlobal[u],subToGlobal[v])] = eid_uv;    //two ver in global graph map to eid in sub graph
            printf("eid_uv: %lu, subU: %u, subV: %u, combine: %lu\n",eid_uv,subToGlobal[u],subToGlobal[v],COMBINE(subToGlobal[u],subToGlobal[v]));
            break;
        }
    }
    if(newFile.edgeNum == 0 || MOVE(sup_uv))
        return ;
    vector<eid_eid> comm;
    IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,maxK,comm,true);

    sup_uv = DELETE(sup_uv);
    // printf("sup_uv: %u\n",sup_uv);
    changeEdgeSup(fOff,fSup,eid_uv,sup_uv);

    if(newFile.edgeNum == 0)
        return;
    newFile.edgeNum--;
    newFileEdge = newFile.edgeNum;

    queue<EdgeSup> del_edge;
    unordered_set<uint64_t> isInDelEdge;

    bool *sE = new bool[nodeNum]();

    for(int i = 0; i < comm.size(); i++){
        uint64_t eid_fir = eid_u[comm[i].first];
        uint64_t eid_sec = eid_v[comm[i].second];

        uint32_t first_sup = sup_u[comm[i].first];
        uint32_t second_sup = sup_v[comm[i].second];
        
        first_sup--;
        if(first_sup < maxK){
            EdgeSup tmp_es;
            tmp_es.u = u, tmp_es.v = nbr_u[comm[i].first];
            tmp_es.eid = eid_fir, tmp_es.sup = first_sup;
            del_edge.push(tmp_es);
            isInDelEdge.insert(eid_fir);
        }
        changeEdgeSup(fOff,fSup,eid_fir,first_sup);

        second_sup--;
        if(second_sup < maxK){
            EdgeSup tmp_es;
            tmp_es.u = v, tmp_es.v = nbr_v[comm[i].second];
            tmp_es.eid = eid_sec, tmp_es.sup = second_sup;
            del_edge.push(tmp_es);
            isInDelEdge.insert(eid_sec);
        }

        changeEdgeSup(fOff,fSup,eid_sec,second_sup);

        // printf("first intersec u: %u, v: %u, w: %u, first_sup: %u, second_sup: %u\n",u,v,nbr_u[comm[i].first],first_sup,second_sup);
    }

    while(!del_edge.empty()){
        EdgeSup tmp_es = del_edge.front();
        del_edge.pop();
        eid_uv = tmp_es.eid, sup_uv = tmp_es.sup, u = tmp_es.u, v = tmp_es.v;
        vis[u] = true, vis[v] = true;
        // printf("tmp_es.u: %u, tmp_es.v: %u, sup: %u\n",u,v,sup_uv);
        loadInfo(nbr_u,degree[u],edgeListBegPtr[u],fDat);
        loadInfo(sup_u,degree[u],edgeListBegPtr[u],fSup);
        loadInfo(eid_u,degree[u],edgeListBegPtr[u]*2,fEid);
        
        loadInfo(nbr_v,degree[v],edgeListBegPtr[v],fDat);
        loadInfo(sup_v,degree[v],edgeListBegPtr[v],fSup);
        loadInfo(eid_v,degree[v],edgeListBegPtr[v]*2,fEid);
        std::vector<eid_eid> comm;
        // IntersectTrussNew(u,nbr_u,eid_u,sup_u,v,nbr_v,eid_v,sup_v,fDat,fEid,fSup,maxK,comm,true);
        uint32_t ptr_u = 0;
        uint32_t upper_u = degree[u];
        uint32_t ptr_v = 0;
        uint32_t upper_v = degree[v];
        uint64_t u_pos = 0,v_pos = 0;
        while(ptr_u < upper_u && ptr_v < upper_v){
            if(MOVE(sup_u[ptr_u]) || sup_u[ptr_u] < maxK)
            {
                ptr_u++;
                continue;
            }
            if(MOVE(sup_v[ptr_v]) || sup_v[ptr_v] < maxK)
            {
                ptr_v++;
                continue;
            }
            if(nbr_u[ptr_u] < nbr_v[ptr_v])
                ptr_u++;
            else if(nbr_u[ptr_u] > nbr_v[ptr_v])
                ptr_v++;
            else
            {
                comm.push_back({ptr_u,ptr_v});
                ptr_u++;
                ptr_v++;
            }
        }

        for(int i = 0; i < comm.size(); i++){
            uint64_t eid_fir = eid_u[comm[i].first];
            uint64_t eid_sec = eid_v[comm[i].second];

            uint32_t first_sup = sup_u[comm[i].first];
            uint32_t second_sup = sup_v[comm[i].second];
            if(first_sup == maxK && isInDelEdge.find(eid_fir) == isInDelEdge.end()){
                // first_sup--;
                EdgeSup tmp_es;
                tmp_es.u = u, tmp_es.v = nbr_u[comm[i].first];
                tmp_es.eid = eid_fir, tmp_es.sup = first_sup;
                del_edge.push(tmp_es);
                isInDelEdge.insert(eid_fir);
            }
            // else if(first_sup > maxK){
            //     fOff.fseek(eid_fir * sizeof(eid_eid));
            //     fOff.fread(&tmp_,sizeof(eid_eid));
            //     first_sup--;
            //     fSup.fseek(tmp_.first*sizeof(uint32_t));
            //     fSup.fwrite(&first_sup,sizeof(uint32_t));
            //     fSup.fseek(tmp_.second*sizeof(uint32_t));
            //     fSup.fwrite(&first_sup,sizeof(uint32_t));
            // }
            
            if(second_sup == maxK && isInDelEdge.find(eid_sec) == isInDelEdge.end()){
                // second_sup--;
                EdgeSup tmp_es;
                tmp_es.u = v, tmp_es.v = nbr_v[comm[i].second];
                tmp_es.eid = eid_sec, tmp_es.sup = second_sup;
                del_edge.push(tmp_es);
                isInDelEdge.insert(eid_sec);
            }
            // else if(second_sup > maxK){
            //     fOff.fseek(eid_sec * sizeof(eid_eid));
            //     fOff.fread(&tmp_,sizeof(eid_eid));
            //     second_sup--;
            //     fSup.fseek(tmp_.first*sizeof(uint32_t));
            //     fSup.fwrite(&second_sup,sizeof(uint32_t));
            //     fSup.fseek(tmp_.second*sizeof(uint32_t));
            //     fSup.fwrite(&second_sup,sizeof(uint32_t));
            // }
        }
        // fOff.fseek(eid_uv * sizeof(eid_eid));
        // fOff.fread(&tmp_,sizeof(eid_eid));
        isInDelEdge.erase(eid_uv);
        // sup_uv = DELETE(sup_uv);
        assert(degree_[u] > 0);
        assert(degree_[v] > 0);
        
        if(newFileEdge == 0)
            return;
        
        newFileEdge--;
        // fSup.fseek(tmp_.first*sizeof(uint32_t));
        // fSup.fwrite(&sup_uv,sizeof(uint32_t));
        // fSup.fseek(tmp_.second*sizeof(uint32_t));
        // fSup.fwrite(&sup_uv,sizeof(uint32_t));
    }
    total_io = total_io + fDat.get_total_io()+ fEid.get_total_io()+ fSup.get_total_io() + fSupSort.get_total_io() + fOff.get_total_io();

    delete[] sE;
    delete[] nbr_u;
    delete[] nbr_v;
    delete[] eid_u;
    delete[] eid_v;
    delete[] sup_u;
    delete[] sup_v;
    fDat.fclose();
    fEid.fclose();
    fSup.fclose();
    fSupSort.fclose();
    fOff.fclose();
}

void Graph::reconMaxTrussGraph(readFile &file, readFile &newFile, bool firstUse)
{
    unsigned long size = 0,es = 0,eid;
    int num = 0, tmpFile = 0;
    uint32_t u,v,max_degree = 0, node_num;

    memset(newFile.m_vertexMap,-1,sizeof(int)*newFile.m_maxID);
    string name = "sort_edge_tmp";
    string sub_dir = newFile.m_base + name;
    newFile.createDir(sub_dir);
    TEdge* edges = new TEdge[file.memEdges];
    char fileName[150];

    string dir;
    if(firstUse)
        dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
    else
        dir = file.m_base;
    // cout<< "dir: " << dir << endl;
    readFile maxTrussFile(dir);
    MyReadFile fSaveSup( maxTrussFile.m_supp );
    fSaveSup.fopen( BUFFERED );
    MyReadFile fSaveOff( maxTrussFile.m_offset );
    fSaveOff.fopen( BUFFERED );
    MyReadFile fSaveEidToVer( maxTrussFile.m_eidToVer );
    fSaveEidToVer.fopen( BUFFERED );
    uint32_t cnt = 0,sup;
    log_debug(graphClock_.Count("subG.edgeNum: %lu, edgeNum: %lu, maxKTruss: %u\n",maxKtrussEdge,edgeNum,maxKtruss));
    
    if(maxKtrussEdge == (maxKtruss+1)*(maxKtruss+2)/2)
        isDeleteEdgeByFunc = false;
    else
        isDeleteEdgeByFunc = true;
    for(uint32_t i = 0; i < maxKtrussEdge; i++){
        eid_eid tmp;
        fSaveOff.fseek(i*sizeof(eid_eid));
        fSaveOff.fread(&tmp,sizeof(eid_eid));
        if(tmp.first == 0 && tmp.second == 0)continue;

        fSaveSup.fseek(tmp.first*sizeof(uint32_t));
        fSaveSup.fread(&sup,sizeof(uint32_t));
        // printf("eid: %u, tmp.first: %lu, sup: %u\n",i,tmp.first,sup);
        if(isDeleteEdgeByFunc && MOVE(sup) && RESTORE(sup) >= maxKtruss){
            Edge tmp_edge;
            fSaveEidToVer.fseek(i*sizeof(Edge));
            fSaveEidToVer.fread(&tmp_edge,sizeof(Edge));
            u = tmp_edge.u, v = tmp_edge.v, eid = i;
            if(m_delBit[u] && m_dynamicDel[u]->find(v) != m_dynamicDel[u]->end()) continue;
            cnt++;
            newFile.moduleInSaveEdgesUnOrder(u,v,num,edges,size,eid,tmpFile,sub_dir);
        }
        if(!isDeleteEdgeByFunc && sup >= maxKtruss){
            cnt++;
            Edge tmp_edge;
            fSaveEidToVer.fseek(i*sizeof(Edge));
            fSaveEidToVer.fread(&tmp_edge,sizeof(Edge));
            u = tmp_edge.u, v = tmp_edge.v, eid = i;
            if(m_delBit[u] && m_dynamicDel[u]->find(v) != m_dynamicDel[u]->end()) continue;
            newFile.moduleInSaveEdgesUnOrder(u,v,num,edges,size,eid,tmpFile,sub_dir);
        }
    }


    sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpFile);
	newFile.saveTmpEdges<TEdge>(edges,size,fileName,[](const TEdge & a, const TEdge & b) {
                if(a.u < b.u)
                    return true;
                if( a.u > b.u )
                    return false;
                return a.v < b.v;
                });
    newFile.edgeNum = num;
	delete[] edges;

    newFile.mergeByDegSort(tmpFile+1, node_num, max_degree, name,false,false);
    log_debug(graphClock_.Count("new subgraph vertex: %d, edge: %lu",newFile.verNum,newFile.edgeNum));

    total_io = total_io + fSaveSup.get_total_io() + fSaveOff.get_total_io() + fSaveEidToVer.get_total_io() + newFile.write_io;
    printf("cnt: %u\n",cnt);
    fSaveSup.fclose();
    fSaveOff.fclose();
    fSaveEidToVer.fclose();
}

void Graph::recoveryKMaxTrussSub(readFile &file, unordered_map<uint32_t, uint32_t>& globalToSub, bool *verMaxKTrussSet, string dir)
{
    uint32_t u,v,eid;
    // string dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
    readFile maxTrussFile(dir);
    MyReadFile fSaveSup( maxTrussFile.m_supp );
    fSaveSup.fopen( NOBUFFER );
    MyReadFile fSaveOff( maxTrussFile.m_offset );
    fSaveOff.fopen( BUFFERED );
    MyReadFile fSaveEidToVer( maxTrussFile.m_eidToVer );
    fSaveEidToVer.fopen( BUFFERED );
    uint32_t cnt = 0,sup;
    // isDeleteEdgeByFunc = true;
    log_info(graphClock_.Count("subG.edgeNum: %lu, edgeNum: %lu, maxKTruss: %u\n",maxKtrussEdge,edgeNum,maxKtruss));
    ofstream out("mytext.txt");
    unordered_map<uint32_t,uint32_t> map;
    memset(verMaxKTrussSet,0,sizeof(bool)*nodeNum);
    int method = 0;
    for(uint32_t i = 0; i < maxKtrussEdge; i++){
        eid_eid tmp;
        fSaveOff.fseek(i*sizeof(eid_eid));
        fSaveOff.fread(&tmp,sizeof(eid_eid));
        if(tmp.first == 0 && tmp.second == 0)continue;

        fSaveSup.fseek(tmp.first*sizeof(uint32_t));
        fSaveSup.fread(&sup,sizeof(uint32_t));
        
        
        if(!isDeleteEdgeByFunc && sup >= maxKtruss){
            cnt++;
            Edge tmp_edge;
            fSaveEidToVer.fseek(i*sizeof(Edge));
            fSaveEidToVer.fread(&tmp_edge,sizeof(Edge));
            verMaxKTrussSet[tmp_edge.u] = true, verMaxKTrussSet[tmp_edge.v] = true;
            method += 0;
        }
        if(isDeleteEdgeByFunc && MOVE(sup) && isInDelQue.find(i) == isInDelQue.end()){
            sup = RESTORE(sup);
            if(sup == maxKtruss){
                Edge tmp_edge;
                fSaveEidToVer.fseek(i*sizeof(Edge));
                fSaveEidToVer.fread(&tmp_edge,sizeof(Edge));
                if(m_delBit[subToGlobal[tmp_edge.u]] && m_dynamicDel[subToGlobal[tmp_edge.u]]->find(subToGlobal[tmp_edge.v]) != m_dynamicDel[subToGlobal[tmp_edge.u]]->end())
                    continue;

                fSaveSup.fseek(tmp.first*sizeof(uint32_t));
                fSaveSup.fwrite(&sup,sizeof(uint32_t));
                fSaveSup.fseek(tmp.second*sizeof(uint32_t));
                fSaveSup.fwrite(&sup,sizeof(uint32_t));
                
                cnt++;
                method += 1;
                out<< "eid: " <<i <<", u: " << tmp_edge.u<< ", v: " << tmp_edge.v << ", sup: " << sup <<endl;
                map[tmp_edge.u]++;
                map[tmp_edge.v]++;
                // printf("eid: %u, u: u%, v: u%, sup: %u\n",i,tmp_edge.u,tmp_edge.v,sup);
                verMaxKTrussSet[tmp_edge.u] = true, verMaxKTrussSet[tmp_edge.v] = true;
            }
        }
    }
    total_io = total_io + fSaveSup.get_total_io() + fSaveOff.get_total_io() + fSaveEidToVer.get_total_io();

    // for(auto it = map.begin(); it != map.end(); it++)
    //     out<<it->first<<",deg: "<< it->second <<endl;

    printf("cnt: %u, method: %u\n",cnt,method);
    fSaveSup.fclose();
    fSaveOff.fclose();
    fSaveEidToVer.fclose();
}


void Graph::recoveryKMaxTruss(readFile &file, unordered_map<uint32_t, uint32_t>& globalToSub, bool *verMaxKTrussSet, string dir)
{
    uint32_t u,v;
    MyReadFile fSup( file.m_supp );
    fSup.fopen( NOBUFFER );
    MyReadFile fSupSort( file.m_suppSort );
    fSupSort.fopen( BUFFERED );
    MyReadFile fOff( file.m_offset );
    fOff.fopen( BUFFERED );
    MyReadFile fSaveEidToVer( file.m_eidToVer );
    fSaveEidToVer.fopen( BUFFERED );
    uint32_t cnt = 0,sup;
    isDeleteEdgeByFunc = true;
    printf("subG.edgeNum: %lu, edgeNum: %lu, maxKTruss: %u, file.edgeNum: %u\n",maxKtrussEdge,edgeNum,maxKtruss,file.edgeNum);
    // ofstream out("mytext.txt");
    unordered_map<uint32_t,uint32_t> map;
    memset(verMaxKTrussSet,0,sizeof(bool)*nodeNum);
    
    file.edgeNum += isInDelQue.size();

    for(int i = 0; i < file.edgeNum; i++){
        EdgeSup tmp_;
        fSupSort.fseek(i*sizeof(EdgeSup));
        fSupSort.fread(&tmp_,sizeof(EdgeSup));
        if(isInDelQue.find(tmp_.eid) != isInDelQue.end())   {
            continue;
        }
        eid_eid eid_;
        fOff.fseek(tmp_.eid*sizeof(eid_eid));
        fOff.fread(&eid_,sizeof(eid_eid)); 

        fSup.fseek(eid_.first*sizeof(uint32_t));
        fSup.fread(&sup,sizeof(uint32_t)); 

        if(MOVE(sup)){
            sup = RESTORE(sup);
            Edge tmp_edge;
            fSaveEidToVer.fseek(tmp_.eid*sizeof(Edge));
            fSaveEidToVer.fread(&tmp_edge,sizeof(Edge));
            // out<< "eid: " <<tmp_.eid <<", u: " << tmp_edge.u<< ", v: " << tmp_edge.v << ", sup: " << sup <<endl;
            if(sup == maxKtruss){
                fSup.fseek(eid_.first*sizeof(uint32_t));
                fSup.fwrite(&sup,sizeof(uint32_t));
                fSup.fseek(eid_.second*sizeof(uint32_t));
                fSup.fwrite(&sup,sizeof(uint32_t));
                
                cnt++;
                // out<< "eid: " <<i <<", u: " << tmp_edge.u<< ", v: " << tmp_edge.v << ", sup: " << sup <<endl;
                map[tmp_edge.u]++;
                map[tmp_edge.v]++;
                // printf("eid: %u, u: u%, v: u%, sup: %u\n",i,tmp_edge.u,tmp_edge.v,sup);
                verMaxKTrussSet[tmp_edge.u] = true, verMaxKTrussSet[tmp_edge.v] = true;
            }
        }

    }
    
    total_io = total_io + fSup.get_total_io() + fOff.get_total_io() + fSaveEidToVer.get_total_io();

    // for(auto it = map.begin(); it != map.end(); it++)
    //     out<<it->first<<",deg: "<< it->second <<endl;

    printf("cnt: %u\n",cnt);
    fSup.fclose();
    fSupSort.fclose();
    fOff.fclose();
    fSaveEidToVer.fclose();
}



void Graph::loadNbrAndSupDynamic(uint32_t u, uint32_t* nbr_, uint32_t* sup_, uint64_t* eid_, uint32_t& _degree, 
 uint64_t begPtr,  MyReadFile& fDat, MyReadFile& fSup, MyReadFile& fEid, 
unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> &sup_dram, 
unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> &eid_dram){
    uint32_t deg = _degree;
    loadNbrForDynamic(u, nbr_, _degree, begPtr, fDat, false);
    _degree = deg;
    loadSupForDynamic(u, sup_, eid_, _degree, begPtr, fSup, fEid, sup_dram, eid_dram);
    std::sort(
      make_zip_iterator(nbr_, sup_, eid_),
      make_zip_iterator(nbr_+_degree, sup_+_degree, eid_+_degree),
      [](const std::tuple<uint32_t, uint32_t, uint64_t>& v, const std::tuple<uint32_t, uint32_t, uint64_t>& w) {
        return std::get<0>(v) < std::get<0>(w);
      });
}

void Graph::loadNbrAndSupDynamicNew(uint32_t u, uint32_t* nbr_, uint32_t* sup_, uint64_t* eid_, uint32_t& _degree, 
 uint64_t begPtr,  MyReadFile& fDat, MyReadFile& fSup, MyReadFile& fEid){
    int purDegree = _degree;
	fDat.fseek(begPtr);
    fSup.fseek(begPtr);
    fEid.fseek(begPtr*2);
	// load all neighbors of vertex u
	_degree = 0;
	uint32_t nbru,supu;
    uint64_t eid;
	for (int i = 0; i < purDegree; ++i) {
		fDat.fread(&nbru, sizeof(int));
        fSup.fread(&supu, sizeof(int));
        fEid.fread(&eid, sizeof(uint64_t));
        if(MOVE(supu)) continue;
        nbr_[_degree] = nbru;
        sup_[_degree] = supu;
        eid_[_degree++] = eid;
	}
 }


void Graph::loadNbrForDynamic(uint32_t u, uint32_t* nbr, uint32_t& _degree, uint64_t pos, MyReadFile& fDat, bool oriGra) 
{
    uint32_t uu = u;
    loadInfo(nbr,_degree,pos,fDat);

	int purDegree = _degree;
	fDat.fseek(pos);
	// load all neighbors of vertex u
	_degree = 0;
	uint32_t t1,t;
    if(oriGra == false)
        u = subToGlobal[u];

    for(int i = 0; i < purDegree; i++)
    {
        t1 = nbr[i];
        if(oriGra == false)
            t = subToGlobal[t1];
		if (m_delBit[u] && m_dynamicDel[u]->find(t) != m_dynamicDel[u]->end()) {
            // printf("u: %u, v: %u in del\n",uu,t1);
			continue;
		}
        nbr[_degree++] = t1;
    }
    
    
	if (m_insBit[u]) {
		int addDegree = m_dynamicIns[u]->size();
        for (auto iter = m_dynamicIns[u]->begin(); iter != m_dynamicIns[u]->end(); iter++){
		// for (int i = 0;i<addDegree;++i) {
			t = *iter;
            if(oriGra == false)
                t= subToGlobal[t];
		    if (!m_delBit[u] || m_dynamicDel[u]->find(subToGlobal[t]) == m_dynamicDel[u]->end()) {
				nbr[_degree++] = t;
			}
		}
	}
}

void Graph::loadSupForDynamic(uint32_t u, uint32_t* nbr, uint64_t* eid_,uint32_t& _degree, uint64_t pos, MyReadFile& fDat, MyReadFile& fEid,
unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> &sup_dram,
unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> &eid_dram) 
{
	int purDegree = _degree;
	fDat.fseek(pos);
    fEid.fseek(pos*2);
	// load all neighbors of vertex u
	_degree = 0;
	int t;
    uint64_t tmp;
	for (int i = 0; i < purDegree; ++i) {
		fDat.fread(&t, sizeof(int));
        fEid.fread(&tmp, sizeof(uint64_t));
        if(MOVE(t)) continue;
        nbr[_degree] = t;
        eid_[_degree++] = tmp;
	}
    u = subToGlobal[u];
	if (m_insBit[u]) {
		int addDegree = m_dynamicIns[u]->size();
        for (auto iter = m_dynamicIns[u]->begin(); iter != m_dynamicIns[u]->end(); iter++){
		// for (int i = 0;i<addDegree;++i) {
			t = subToGlobal[*iter];
			if (!m_delBit[u] || m_dynamicDel[u]->find(t) == m_dynamicDel[u]->end()) {
				nbr[_degree] = sup_dram[u][t];
                eid_[_degree++] = eid_dram[u][t];
			}
		}
	}
}

uint32_t Graph::selectNbr(readFile &file, uint32_t a, int index){
    MyReadFile fIdx(file.m_idx);
	fIdx.fopen(BUFFERED);
	MyReadFile fDat(file.m_dat);
	fDat.fopen(BUFFERED);
	uint32_t* nbr = new uint32_t[file.maxDeg];
	uint32_t deg = degree[a];

    // #ifdef DegSort
	// fIdx.fseek(a*(sizeof(long)*2 + sizeof(uint32_t)));
    // #else
    // fIdx.fseek(a*(sizeof(long) + sizeof(uint32_t)));
    // #endif

	long pos, posPlus;
	// fIdx.fread(&pos, sizeof(long));
    // fIdx.fread(&posPlus, sizeof(long));
	// fIdx.fread(&deg, sizeof(uint32_t));
	// fDat.fseek(pos);

    pos = edgeListBegPtr[a];
    posPlus = edgeListBegPtrPlus[a];

    uint32_t u_nbrNum = deg-(posPlus-pos)/sizeof(uint32_t);
    if(u_nbrNum == 0 || index == u_nbrNum){
        delete[] nbr;
        return 0;
    }

    loadInfo(nbr,u_nbrNum,posPlus,fDat);
	uint32_t r = nbr[index];

	delete[] nbr;
	fDat.fclose();
	fIdx.fclose();
	return r;
} 

void Graph::generateRandomEdges(readFile &file, uint32_t &generate_edge, Edge *dynamic, bool *isInMaxTruss, unordered_map<uint32_t, uint32_t>& globalToSub)
{
    int c = 2;
    /* method 1*/
    // uint32_t num = generate_edge;
    // int sep = nodeNum / num;
    // cout << "sep: " << sep <<"  nodeNum: "<<nodeNum<< "  num: " << num << endl;
    // uint32_t curVertex = 0, cut = 0;
    // for (int i = 0; curVertex < nodeNum && i < num; i++) {
    //     uint32_t b = selectNbr(file,curVertex,0);
    //     if (b%num == 0) {    
    //         i -= 1;
    //         curVertex += 1;
    //     }
    //     else{
    //         ++cut;
    //         dynamic[i] = {curVertex,b};
    //         curVertex += sep;
    //     }
    // }
    // // num -= cut;
    // generate_edge = cut;

    /* method 2*/
    srand(time(0));
    /* Seed */ 
    std::random_device rd; 
    /* Random number generator */ 
    std::default_random_engine generator(rd()); 
    /* Distribution on which to apply the generator */ 
    std::uniform_int_distribution<long long unsigned> distribution(0,0xFFFFFFFFFFFFFFFF);

    unordered_set<uint32_t> **exist = new unordered_set<uint32_t>*[nodeNum];
	memset(exist, 0, sizeof(unordered_set<uint32_t>*)*nodeNum);
    bool *isInExist = new bool[nodeNum]();

    for(int i = 0; i < generate_edge-c; ){
        uint32_t start = distribution(generator) % nodeNum;
        // if(!isInMaxTruss[start])
        //     continue;
        if(!isInExist[start]){
            isInExist[start] = true;
            exist[start] = new unordered_set<uint32_t>;  
            uint32_t end = selectNbr(file,start,0);
            // if(end != 0 && isInMaxTruss[end]){
            if(end != 0){
                exist[start]->insert(end);
                dynamic[i] = {start,end};
                i++;
            } 
        }
        else{
            int index = exist[start]->size();
            uint32_t end = selectNbr(file,start,index);
            // if(end != 0 && end < nodeNum && isInMaxTruss[end]){
            if(end != 0 && end < nodeNum){
                exist[start]->insert(end);
                dynamic[i] = {start,end};
                i++;
            } 
        }
    }

    MyReadFile fDat(file.m_dat);
	fDat.fopen(BUFFERED);
    uint32_t *nbr = new uint32_t[file.maxDeg]();
    c = generate_edge - c;
    // printf("c_start: %u\n",c);
    for(auto it = globalToSub.begin(); it != globalToSub.end() && c < generate_edge; it++)
    {
        uint32_t i = it->first;
    // }
    // for(uint32_t i = 0; i < nodeNum && c < generate_edge; i++){
        if(!isInMaxTruss[globalToSub[i]]) continue;
        uint64_t pos = edgeListBegPtr[i];
        uint64_t posPlus = edgeListBegPtrPlus[i];
        uint32_t u_nbrNum = degree[i]-(posPlus-pos)/sizeof(uint32_t);
        if(u_nbrNum == 0) continue;
        loadInfo(nbr,u_nbrNum,posPlus,fDat);
        for(int j = 0; j < u_nbrNum && c < generate_edge; j += 3)
        {
            if(globalToSub.find(nbr[j]) == globalToSub.end()) continue;
            if(!isInMaxTruss[globalToSub[nbr[j]]]) continue;
            if(isInExist[i] == false || (isInExist[i] && exist[i]->find(nbr[j]) == exist[i]->end()))
            {
                dynamic[c++] = {i,nbr[j]};
            }
        }
    }
    fDat.fclose();
    delete[] nbr;
    delete[] isInExist;

    std::sort(dynamic,dynamic+generate_edge,[](const Edge & a, const Edge & b) {
        if(a.u == b.u)
            return a.v < b.v;
        return a.u < b.u;
    });
    
    
    for(int i = 0; i < nodeNum; i++){
        if(exist[i])
            delete exist[i];
    }
    delete[] exist;
    log_info(graphClock_.Count("select %d edges, c: %d\n", generate_edge,c));
}


void Graph::dynamicMaxTrussDeletion_YLJ(readFile &file){
    total_io = 0;
    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );

    m_dynamicDel = new unordered_set<int>*[nodeNum];
	memset(m_dynamicDel, 0, sizeof(unordered_set<int>*)*nodeNum);
    m_dynamicIns = new unordered_set<int>*[nodeNum];
	memset(m_dynamicIns, 0, sizeof(unordered_set<int>*)*nodeNum);

    m_delBit = new bool[nodeNum]();
    m_insBit = new bool[nodeNum]();

    memcpy(degree_, degree, sizeof(uint32_t)*nodeNum);
    memset(vis,0,nodeNum*sizeof(bool));

    unordered_map<uint32_t, uint32_t> globalToSub;

    uint32_t generate_edge = 1000;
    Edge *dynamic = new Edge[generate_edge];

    uint32_t *coreNum = new uint32_t[nodeNum]();  
    uint32_t *nbr_u = new uint32_t[nodeNum](); 
    uint32_t *nbr_v = new uint32_t[nodeNum](); 
    KCore(file,coreNum,false);

    bool *isInMaxTruss = new bool[nodeNum]();
    // for(int i = 0; i < file.vertexId.size(); i++)
    //     isInMaxTruss[file.vertexId[i]] = true;  //subgraph
    // printf("size: %d, nodeNum: %u\n",file.vertexId.size(),nodeNum);
    
    bool *verMaxKTruSet = new bool[nodeNum]();
    bool *isInSubG = new bool[nodeNum](); 

    string dir = file.m_base + "subGraphInfo/";
    readFile subG(dir);
    uint32_t subG_nodeNum = 0, subG_maxDeg = 0;
    uint64_t subG_edgeNum = 0;
    MyReadFile fSubGInfo( subG.m_info );
	fSubGInfo.fopen( NOBUFFER );
	fSubGInfo.fread(&subG_nodeNum,sizeof(uint32_t));
    fSubGInfo.fread(&subG_maxDeg,sizeof(uint32_t));
    fSubGInfo.fread(&subG_edgeNum,sizeof(uint64_t));

    MyReadFile fileisDeleteEdgeByFunc(file.m_base + "graph.isDeleteEdgeByFunc");
    fileisDeleteEdgeByFunc.fopen( NOBUFFER );
	fileisDeleteEdgeByFunc.fread(&isDeleteEdgeByFunc,sizeof(bool));
    fileisDeleteEdgeByFunc.fclose();

    int len = 0;
    MyReadFile fileVertexId(file.m_base + "graph.fileVertexId");
    fileVertexId.fopen( BUFFERED );
    fileVertexId.fread(&len,sizeof(int));
    file.vertexId.clear();
    for(int i = 0; i < len; i++){
        uint32_t ver;
        fileVertexId.fread(&ver,sizeof(uint32_t));
        file.vertexId.push_back(ver);
    }
    fileVertexId.fclose();

    maxKtrussEdge = subG_edgeNum;
    MyReadFile fSubToGlobal( file.m_base+"graph.subToGlobal" );
	fSubToGlobal.fopen( BUFFERED );
    for(int i = 0; i < subG_nodeNum; i++){
        Edge e;
        fSubToGlobal.fread(&e,sizeof(Edge));
        subToGlobal[e.u] = e.v;
    }
    fSubToGlobal.fclose();

    for(auto it = subToGlobal.begin(); it != subToGlobal.end(); it++)
        globalToSub[it->second] = it->first;

    MyReadFile fFinal( file.m_base+"graph.final" );
	fFinal.fopen( BUFFERED );
    fFinal.fread(&last_sup,sizeof(uint32_t));
    fFinal.fread(&maxKtruss,sizeof(uint32_t));
    fFinal.fclose();

    dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
    recoveryKMaxTrussSub(file,globalToSub,verMaxKTruSet,dir);
    generateRandomEdges(file,generate_edge,dynamic,verMaxKTruSet,globalToSub); 

    readFile maxTrussFile(dir);

    vector<uint32_t> degSub(nodeNum,0);
    vector<uint64_t> begPtr(nodeNum,0);

    MyReadFile fIdx( maxTrussFile.m_idx );
	fIdx.fopen( BUFFERED );

    uint64_t tmp = 0,tmpb;
    uint32_t degreeTmp = 0;

    printf("file.vertexId.size(): %u, globalToSub.size(): %u\n",file.vertexId.size(),globalToSub.size());
    
    for (uint32_t u = 0; u < file.vertexId.size(); ++u){
        int i = file.vertexId[u];
		fIdx.fread(&tmp,sizeof(uint64_t));
        #ifdef DegSort
        fIdx.fread(&tmpb,sizeof(uint64_t));
        #endif
		fIdx.fread(&degreeTmp,sizeof(uint32_t));
        begPtr[i] = tmp;
        degSub[i] = degreeTmp;
  	}
	fIdx.fclose();

    for(int i = 0; i < generate_edge; i++){
        uint32_t a = dynamic[i].u, b = dynamic[i].v;
        if (!m_delBit[a]) {
            m_dynamicDel[a] = new unordered_set<int>;
            m_delBit[a] = true;
        }
        if (!m_delBit[b]) {
            m_dynamicDel[b] = new unordered_set<int>;
            m_delBit[b] = true;
        }
        degree_[a]--,degree_[b]--;
        m_dynamicDel[a]->insert(b);
        m_dynamicDel[b]->insert(a);
        uint32_t oriCoreU = coreNum[a], oriCoreV = coreNum[b];

        
        if(globalToSub.find(a) != globalToSub.end() && globalToSub.find(b) != globalToSub.end()){
            // printf("a: %u, degree[a]: %u,globalToSub[a]: %u, b: %u, degree[b]: %u,globalToSub[b]: %u\n",a,degree[a],globalToSub[a],b,degree[b],globalToSub[b]);
            dynamic[i].u = globalToSub[a], dynamic[i].v = globalToSub[b];
            if(!verMaxKTruSet[globalToSub[a]] || !verMaxKTruSet[globalToSub[b]]) continue;
            uint32_t old_maxKtruss = maxKtruss;
            dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
            readFile maxTrussFileTmp(dir);

            delEdgeDynamicYLJ(dynamic[i],maxTrussFileTmp,verMaxKTruSet,begPtr,degSub);
            if(maxKtruss < old_maxKtruss && maxKtruss >= last_sup){
                dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
                recoveryKMaxTrussSub(file,globalToSub,verMaxKTruSet,dir);
            }
            else if(maxKtruss < last_sup){
                last_sup = maxKtruss;
                KCore(file,coreNum,true);
                memset(isInSubG,0,nodeNum*sizeof(bool));
                for(int i = 0; i < nodeNum; i++)
                    if(coreNum[i] >= maxKtruss+1)
                        isInSubG[i] = true;
                dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
                file.createDir(dir);
                readFile subFile(dir);
                filterGlobalIntoSub(isInSubG,file,subFile,true); 
                globalToSub.clear();
                for(auto it = subToGlobal.begin(); it != subToGlobal.end(); it++)
                    globalToSub[it->second] = it->first;
                fill(begPtr.begin(), begPtr.end(), 0);
                fill(degSub.begin(), degSub.end(), 0);
                Graph subG(subFile.verNum,subFile.edgeNum);
                subG_edgeNum = subFile.edgeNum, subG_nodeNum = subFile.verNum;
                printf("subFile.verNum: %u, subFile.edgeNum: %lu\n",subFile.verNum,subFile.edgeNum);
                subG.Initial(subFile);
                for(int k = 0; k < subG.nodeNum; k++){
                    begPtr[k] = subG.edgeListBegPtr[k];
                    degSub[k] = subG.degree[k];
                }
                subG.CountTriangleSSDByDegOrder(subFile,true);
                subG.trussDecomLazyUpdate(subFile);
                maxKtruss = subG.last_sup;
                maxKtrussEdge = subG.edgeNum;
                dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
                isDeleteEdgeByFunc = true;
                recoveryKMaxTrussSub(file,globalToSub,verMaxKTruSet,dir);
                total_io += subG.total_io;
            }

        }
    }
    total_io = total_io + fSubGInfo.get_total_io() + fDat.get_total_io();
    log_debug(graphClock_.Count("delete total_io: %lu\n",total_io));

    total_io = 0;

    unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> sup_dram;
    unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> eid_dram;

    for(int i = generate_edge-1; i >=0; i--){
        uint32_t a = dynamic[i].u, b = dynamic[i].v;
        if (m_delBit[a]){
            m_dynamicDel[a]->erase(b);
            if(m_dynamicDel[a]->size() == 0) {
                m_delBit[a] = false;
                delete m_dynamicDel[a];
            }
        }
        if (m_delBit[b]){
            m_dynamicDel[b]->erase(a);
            if(m_dynamicDel[b]->size() == 0) {
                m_delBit[b] = false;
                delete m_dynamicDel[b];
            }
        }
        if (!m_insBit[a]) {
            m_dynamicIns[a] = new unordered_set<int>;
            m_insBit[a] = true;
        }
        if (!m_insBit[b]) {
            m_dynamicIns[b] = new unordered_set<int>;
            m_insBit[b] = true;
        }
        degree_[a]++,degree_[b]++;
        m_dynamicIns[a]->insert(b);
        m_dynamicIns[b]->insert(a);
        uint32_t oriCoreU = coreNum[a], oriCoreV = coreNum[b];

        if(degree_[a] >= last_sup && degree_[b] >= last_sup)
            UpdateCoreDynamic(file,coreNum,a,b);
            // KCore(file,coreNum,true);
        eid_dram[a][b] = subG_edgeNum++;
        
        if(globalToSub.find(a) != globalToSub.end() && globalToSub.find(b) != globalToSub.end()){
            printf("a: %u, degree[a]: %u,globalToSub[a]: %u, b: %u, degree[b]: %u,globalToSub[b]: %u\n",a,degree[a],globalToSub[a],b,degree[b],globalToSub[b]);
            dynamic[i].u = globalToSub[a], dynamic[i].v = globalToSub[b];
            insEdgeDynamic(dynamic[i],maxTrussFile,sup_dram,eid_dram,verMaxKTruSet,begPtr,degSub);
        }
        else if(coreNum[dynamic[i].u] >= last_sup+1 && coreNum[dynamic[i].v] >= last_sup+1){
            vector<uint32_t> vec_u, vec_v;
            MyReadFile fDatSub( maxTrussFile.m_dat );
	        fDatSub.fopen( BUFFERED );

            if(globalToSub.find(dynamic[i].u) == globalToSub.end())
            {
                loadInfo(nbr_u,degree[dynamic[i].u],edgeListBegPtr[dynamic[i].u],fDat);
                for(int j = 0; j < degree[dynamic[i].u]; j++){
                    if(coreNum[nbr_u[j]] < last_sup+1) continue;
                    vec_u.push_back(nbr_u[j]);
                }
            }
            else{
                uint32_t u_degree = degSub[globalToSub[a]];
                loadNbrForDynamic(globalToSub[a],nbr_u,u_degree,begPtr[globalToSub[a]],fDatSub,true);
                for(int j = 0; j < u_degree; j++){
                    vec_u.push_back(subToGlobal[nbr_u[j]]);
                }
            }

            if(globalToSub.find(dynamic[i].v) == globalToSub.end())
            {
                loadInfo(nbr_u,degree[dynamic[i].v],edgeListBegPtr[dynamic[i].v],fDat);
                for(int j = 0; j < degree[dynamic[i].v]; j++){
                    if(coreNum[nbr_u[j]] < last_sup+1) continue;
                    vec_u.push_back(nbr_u[j]);
                }
            }
            else{
                uint32_t u_degree = degSub[globalToSub[b]];
                loadNbrForDynamic(globalToSub[b],nbr_u,u_degree,begPtr[globalToSub[b]],fDatSub,true);
                for(int j = 0; j < u_degree; j++){
                    vec_v.push_back(subToGlobal[nbr_u[j]]);
                }
            }
            uint32_t u_pos = 0, v_pos = 0, count = 0;
            while(u_pos < vec_u.size() && v_pos < vec_v.size()){
                if(vec_u[u_pos] < vec_v[v_pos]) u_pos++;
                else if(vec_u[u_pos] > vec_v[v_pos]) v_pos++;
                else{
                    count++;
                    u_pos++, v_pos++;
                }
            }
            total_io += fDatSub.get_total_io();
            fDatSub.fclose();
            if(count < last_sup) continue; 

            memset(isInSubG, false, sizeof(bool)*nodeNum);
            for(int k = 0; k < nodeNum; k++)
                if(coreNum[k] >= last_sup+1)
                    isInSubG[k] = true;
            dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
            readFile subFile(dir);
            filterGlobalIntoSub(isInSubG,file,subFile,true); 
            globalToSub.clear();
            for(auto it = subToGlobal.begin(); it != subToGlobal.end(); it++)
                globalToSub[it->second] = it->first;
            fill(begPtr.begin(), begPtr.end(), 0);
            fill(degSub.begin(), degSub.end(), 0);

            Graph subG(subFile.verNum,subFile.edgeNum);
            subG_edgeNum = subFile.edgeNum, subG_nodeNum = subFile.verNum;
            printf("subFile.verNum: %u, subFile.edgeNum: %lu\n",subFile.verNum,subFile.edgeNum);
            subG.Initial(subFile);
            for(int k = 0; k < subG.nodeNum; k++){
                begPtr[k] = subG.edgeListBegPtr[k];
                degSub[k] = subG.degree[k];
            }
            subG.CountTriangleSSDByDegOrder(subFile,true);
            subG.trussDecomLazyUpdate(subFile);
            last_sup = subG.last_sup;
            dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
            recoveryKMaxTrussSub(file,globalToSub,verMaxKTruSet,dir);
            maxTrussFile.change(dir);
            total_io += subG.total_io;
            dynamic[i].u = globalToSub[dynamic[i].u], dynamic[i].v = globalToSub[dynamic[i].v];
            printf("core: ori_u %u, now %u, ori_v %u, now %u\n",oriCoreU,coreNum[a],oriCoreV,coreNum[b]);
            insEdgeDynamic(dynamic[i],maxTrussFile,sup_dram,eid_dram,verMaxKTruSet,begPtr,degSub);
        }
    }

    total_io = total_io + fSubGInfo.get_total_io() + fDat.get_total_io();
    log_debug(graphClock_.Count("total_io: %lu\n",total_io));


    delete[] dynamic;
    delete[] m_insBit;
    delete[] m_delBit;
    delete[] m_dynamicDel;
    delete[] coreNum;
    delete[] nbr_u;
    delete[] nbr_v;
    fSubGInfo.fclose();
    fDat.fclose();
}

void Graph::delEdgeDynamicYLJ(Edge e, readFile &file, bool *verMaxKTruSet, vector<uint64_t> &begPtr, vector<uint32_t> &degSub){
    printf("u: %u, degsub[u]: %u, v: %u\n",e.u,degSub[e.u],e.v);
    MyReadFile fInfo(file.m_info);
    fInfo.fopen(BUFFERED);
    MyReadFile fIdx(file.m_idx);
	fIdx.fopen(BUFFERED);
	MyReadFile fDat(file.m_dat);
	fDat.fopen(BUFFERED); 
    MyReadFile fEid(file.m_eid);
	fEid.fopen(BUFFERED);
	MyReadFile fOff(file.m_offset);
    fOff.fopen(BUFFERED);
    MyReadFile fSup(file.m_supp);
	fSup.fopen(NOBUFFER); 

    uint32_t u = e.u, v = e.v, u_degree = 0, v_degree = 0, subGraV = 0, sup_tmp = 0;
    uint64_t pos = 0, posPlus = 0, eid = 0;
    uint32_t* nbr_u = new uint32_t[nodeNum]();
    uint32_t* nbr_v = new uint32_t[nodeNum]();

    uint64_t* eid_u = new uint64_t[nodeNum]();
    uint64_t* eid_v = new uint64_t[nodeNum]();

    uint32_t* sup_u = new uint32_t[nodeNum]();
    uint32_t* sup_v = new uint32_t[nodeNum]();
    
    loadInfo(nbr_u,degSub[u],begPtr[u],fDat);
    loadInfo(sup_u,degSub[u],begPtr[u],fSup);
    loadInfo(eid_u,degSub[u],begPtr[u]*2,fEid);
    bool find = false;
    for(int i = 0; i < degSub[u]; i++){
        printf("%u ",nbr_u[i]);
        if(nbr_u[i] == v){
            eid = eid_u[i];
            sup_tmp = sup_u[i];
            find = true;
            break;
        }
    }
    if(!find)
    {
        printf("error, not find\n");
        exit(0);
    }


    struct cmp{
        bool operator()(EdgeSup a, EdgeSup b){
            return a.sup > b.sup;
        }
    };
    std::priority_queue<EdgeSup,std::vector<EdgeSup>, cmp> pq; 

    eid_eid tmp;
    fOff.fseek(eid*sizeof(eid_eid));
    fOff.fread(&tmp,sizeof(eid_eid));
    if(sup_tmp != maxKtruss){
        printf("error, sup_tmp: %u, eid: %lu, degSub[u]: %u\n",sup_tmp,eid,degSub[u]);
        exit(0);
    } 
    sup_tmp = DELETE(sup_tmp);
    fSup.fseek(tmp.first*sizeof(uint32_t));
    fSup.fwrite(&sup_tmp,sizeof(uint32_t));
    fSup.fseek(tmp.second*sizeof(uint32_t));
    fSup.fwrite(&sup_tmp,sizeof(uint32_t));

    loadInfo(nbr_v,degSub[v],begPtr[v],fDat);
    loadInfo(sup_v,degSub[v],begPtr[v],fSup);
    loadInfo(eid_v,degSub[v],begPtr[v]*2,fEid);

    std::vector<ver_eid_eid> comm;
    IntersectOperaDelDynamic(u,nbr_u,sup_u,eid_u,degSub[u],v,nbr_v,sup_v,eid_v,degSub[v],comm);

    queue<TEdge> que;
    unordered_map<uint32_t,unordered_map<uint32_t,bool>> isInLk;
    unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> support;

    for(int i = 0; i < comm.size(); i++){
        if(comm[i].u_sup == maxKtruss){
            uint32_t x = min(u, comm[i].w);
            uint32_t y = max(u, comm[i].w);
            que.push({x,y,comm[i].first});
            isInLk[x][y] = true;
        }
        if(comm[i].v_sup == maxKtruss){
            uint32_t x = min(v, comm[i].w);
            uint32_t y = max(v, comm[i].w);
            que.push({x,y,comm[i].second});
            isInLk[x][y] = true;
        }
    }

    while(!que.empty()){
        TEdge tmpe = que.front();
        u = tmpe.u, v = tmpe.v, eid = tmpe.eid;
        que.pop();
        support[u][v] = 0;

        loadInfo(nbr_u,degSub[u],begPtr[u],fDat);
        loadInfo(sup_u,degSub[u],begPtr[u],fSup);
        loadInfo(eid_u,degSub[u],begPtr[u]*2,fEid);

        loadInfo(nbr_v,degSub[v],begPtr[v],fDat);
        loadInfo(sup_v,degSub[v],begPtr[v],fSup);
        loadInfo(eid_v,degSub[v],begPtr[v]*2,fEid);
        comm.clear();
        IntersectOperaDelDynamic(u,nbr_u,sup_u,eid_u,degSub[u],v,nbr_v,sup_v,eid_v,degSub[v],comm);

        for(int i = 0; i < comm.size(); i++){
            if(comm[i].u_sup < maxKtruss || comm[i].v_sup < maxKtruss) continue;
            support[u][v] += 1;
            // printf("u: %u, v: %u, w: %u, comm[i].u_sup: %u, comm[i].v_sup: %u\n",u,v,comm[i].w,comm[i].u_sup,comm[i].v_sup); 
            uint32_t x = min(u,comm[i].w);
            uint32_t y = max(u,comm[i].w);
            if(comm[i].u_sup == maxKtruss){
                if(isInLk.find(x) == isInLk.end() || isInLk[x].find(y) == isInLk[x].end()){
                    que.push({x,y,comm[i].first});
                    isInLk[x][y] = true;
                }
            }
            x = min(v,comm[i].w);
            y = max(v,comm[i].w);
            if(comm[i].v_sup == maxKtruss){
                if(isInLk.find(x) == isInLk.end() || isInLk[x].find(y) == isInLk[x].end()){
                    que.push({x,y,comm[i].second});
                    isInLk[x][y] = true;
                }
            }
        } 
        pq.push({u,v,support[u][v],eid});
    }
    log_debug(graphClock_.Count("after for! pq.size(): %d",pq.size()));

    while(!pq.empty()){
        EdgeSup es_tmp = pq.top();  pq.pop();
        u = es_tmp.u, v = es_tmp.v, eid = es_tmp.eid;
        uint32_t sup = es_tmp.sup;
        if(isInLk.find(u) == isInLk.end() || isInLk[u].find(v) == isInLk[u].end()) continue;
        if(sup >= maxKtruss) {
            pq.push({u,v,sup,eid});
            break;
        }

        loadInfo(nbr_u,degSub[u],begPtr[u],fDat);
        loadInfo(sup_u,degSub[u],begPtr[u],fSup);
        loadInfo(eid_u,degSub[u],begPtr[u]*2,fEid);

        loadInfo(nbr_v,degSub[v],begPtr[v],fDat);
        loadInfo(sup_v,degSub[v],begPtr[v],fSup);
        loadInfo(eid_v,degSub[v],begPtr[v]*2,fEid);
        comm.clear();
        IntersectOperaDelDynamic(u,nbr_u,sup_u,eid_u,degSub[u],v,nbr_v,sup_v,eid_v,degSub[v],comm);

        for(int i = 0; i < comm.size(); i++){
            if(comm[i].u_sup < maxKtruss || comm[i].v_sup < maxKtruss) continue;
            uint32_t x = min(u,comm[i].w);
            uint32_t y = max(u,comm[i].w);
            if(comm[i].u_sup == maxKtruss){
                if(isInLk.find(x) == isInLk.end() || isInLk[x].find(y) == isInLk[x].end())
                    continue;
            }
            uint32_t xx = min(v,comm[i].w);
            uint32_t yy = max(v,comm[i].w);
            if(comm[i].v_sup == maxKtruss){
                if(isInLk.find(xx) == isInLk.end() || isInLk[xx].find(yy) == isInLk[xx].end())
                    continue;  
            }

            support[x][y] -= 1;
            if(support[x][y] < maxKtruss && support[x][y] >= sup) pq.push({x,y,support[x][y],comm[i].first});
            support[xx][yy] -= 1;
            if(support[xx][yy] < maxKtruss && support[xx][yy] >= sup) pq.push({xx,yy,support[xx][yy],comm[i].second});
        }
        isInLk[u].erase(v);
        support[u].erase(v);
        if(isInLk[u].size() == 0)
        {
            isInLk.erase(u);
            support.erase(u);
        }
        sup = maxKtruss - 1;
        sup = DELETE(sup);
        changeEdgeSup(fOff,fSup,eid,sup);
    }
    log_debug(graphClock_.Count("after while! pq.size(): %d",pq.size()));
    memset(verMaxKTruSet,0,nodeNum*sizeof(bool));

    if(pq.size() == 0){
        maxKtruss--;
    }
    else{
        while(!pq.empty()){
            EdgeSup es_tmp = pq.top();
            pq.pop();
            u = es_tmp.u, v = es_tmp.v, eid = es_tmp.eid;
            if(es_tmp.sup >= maxKtruss)
                verMaxKTruSet[u] = true, verMaxKTruSet[v] = true;
        }
    }
    log_debug(graphClock_.Count("maxKTruss: %d, delete successfully!", maxKtruss));
    total_io = total_io + fInfo.get_total_io() + fIdx.get_total_io() + fDat.get_total_io() + fSup.get_total_io() + fOff.get_total_io() + fEid.get_total_io();
    delete[] nbr_u;
    delete[] nbr_v;

    delete[] eid_u;
    delete[] eid_v;

    delete[] sup_u;
    delete[] sup_v;
    fInfo.fclose();
    fIdx.fclose();
    fDat.fclose();
    fSup.fclose();
    fOff.fclose();
    fEid.fclose();
}

void Graph::dynamicMaxTrussInsertion_YLJ(readFile &file){
    total_io = 0;
    MyReadFile fDat( file.m_dat );
	fDat.fopen( BUFFERED );

    m_dynamicIns = new unordered_set<int>*[nodeNum];
	memset(m_dynamicIns, 0, sizeof(unordered_set<int>*)*nodeNum);

    m_dynamicDel = new unordered_set<int>*[nodeNum];
	memset(m_dynamicDel, 0, sizeof(unordered_set<int>*)*nodeNum);

    m_insBit = new bool[nodeNum]();
    m_delBit = new bool[nodeNum]();

    memcpy(degree_, degree, sizeof(uint32_t)*nodeNum);

    unordered_map<uint32_t, uint32_t> globalToSub;
    for(auto it = subToGlobal.begin(); it != subToGlobal.end(); it++)
        globalToSub[it->second] = it->first;

    uint32_t generate_edge = 1000;
    Edge *dynamic = new Edge[generate_edge];

    uint32_t *coreNum = new uint32_t[nodeNum]();  
    uint32_t *nbr_u = new uint32_t[nodeNum](); 
    uint32_t *nbr_v = new uint32_t[nodeNum](); 
    KCore(file,coreNum,false);

    
    bool *verMaxKTruSet = new bool[nodeNum]();
    bool *isInSubG = new bool[nodeNum](); 

    string dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
    recoveryKMaxTrussSub(file,globalToSub,verMaxKTruSet,dir);
    
    generateRandomEdgesInsertion(file,generate_edge,dynamic); 

    readFile maxTrussFile(dir);

    dir = file.m_base + "subGraphInfo/";
    readFile subG(dir);
    uint32_t subG_nodeNum = 0, subG_maxDeg = 0;
    uint64_t subG_edgeNum = 0;
    MyReadFile fSubGInfo( subG.m_info );
	fSubGInfo.fopen( NOBUFFER );
	fSubGInfo.fread(&subG_nodeNum,sizeof(uint32_t));
    fSubGInfo.fread(&subG_maxDeg,sizeof(uint32_t));
    fSubGInfo.fread(&subG_edgeNum,sizeof(uint64_t));
    
    vector<uint32_t> degSub(nodeNum,0);
    vector<uint64_t> begPtr(nodeNum,0);

    MyReadFile fIdx( maxTrussFile.m_idx );
	fIdx.fopen( BUFFERED );
    uint64_t tmp = 0,tmpb;
    uint32_t degreeTmp = 0;

    printf("file.vertexId.size(): %u, globalToSub.size(): %u\n",file.vertexId.size(),globalToSub.size());
    
    for (uint32_t u = 0; u < file.vertexId.size(); ++u){
        int i = file.vertexId[u];
		fIdx.fread(&tmp,sizeof(uint64_t));
        #ifdef DegSort
        fIdx.fread(&tmpb,sizeof(uint64_t));
        #endif
		fIdx.fread(&degreeTmp,sizeof(uint32_t));
        begPtr[i] = tmp;
        degSub[i] = degreeTmp;
        // printf("i: %u, tmp: %u\n",i,tmp);
  	}
	fIdx.fclose();


    unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> sup_dram;
    unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> eid_dram;

    for(int i = 0; i < generate_edge; i++){
        uint32_t a = dynamic[i].u, b = dynamic[i].v;
        if (!m_insBit[a]) {
            m_dynamicIns[a] = new unordered_set<int>;
            m_insBit[a] = true;
        }
        if (!m_insBit[b]) {
            m_dynamicIns[b] = new unordered_set<int>;
            m_insBit[b] = true;
        }
        degree_[a]++,degree_[b]++;
        m_dynamicIns[a]->insert(b);
        m_dynamicIns[b]->insert(a);
        uint32_t oriCoreU = coreNum[a], oriCoreV = coreNum[b];

        if(degree_[a] >= last_sup && degree_[b] >= last_sup)
            UpdateCoreDynamic(file,coreNum,a,b);
            // KCore(file,coreNum,true);
        eid_dram[a][b] = subG_edgeNum++;
        
        if(globalToSub.find(a) != globalToSub.end() && globalToSub.find(b) != globalToSub.end()){
            printf("a: %u, degree[a]: %u,globalToSub[a]: %u, b: %u, degree[b]: %u,globalToSub[b]: %u\n",a,degree[a],globalToSub[a],b,degree[b],globalToSub[b]);
            dynamic[i].u = globalToSub[a], dynamic[i].v = globalToSub[b];
            insEdgeDynamic(dynamic[i],maxTrussFile,sup_dram,eid_dram,verMaxKTruSet,begPtr,degSub);
        }
        else if(coreNum[dynamic[i].u] >= last_sup+1 && coreNum[dynamic[i].v] >= last_sup+1){

            vector<uint32_t> vec_u, vec_v;
            MyReadFile fDatSub( maxTrussFile.m_dat );
	        fDatSub.fopen( BUFFERED );

            if(globalToSub.find(dynamic[i].u) == globalToSub.end())
            {
                loadInfo(nbr_u,degree[dynamic[i].u],edgeListBegPtr[dynamic[i].u],fDat);
                for(int j = 0; j < degree[dynamic[i].u]; j++){
                    if(coreNum[nbr_u[j]] < last_sup+1) continue;
                    vec_u.push_back(nbr_u[j]);
                }
            }
            else{
                uint32_t u_degree = degSub[globalToSub[a]];
                loadNbrForDynamic(globalToSub[a],nbr_u,u_degree,begPtr[globalToSub[a]],fDatSub,true);
                for(int j = 0; j < u_degree; j++){
                    vec_u.push_back(subToGlobal[nbr_u[j]]);
                }
            }

            if(globalToSub.find(dynamic[i].v) == globalToSub.end())
            {
                loadInfo(nbr_u,degree[dynamic[i].v],edgeListBegPtr[dynamic[i].v],fDat);
                for(int j = 0; j < degree[dynamic[i].v]; j++){
                    if(coreNum[nbr_u[j]] < last_sup+1) continue;
                    vec_u.push_back(nbr_u[j]);
                }
            }
            else{
                uint32_t u_degree = degSub[globalToSub[b]];
                loadNbrForDynamic(globalToSub[b],nbr_u,u_degree,begPtr[globalToSub[b]],fDatSub,true);
                for(int j = 0; j < u_degree; j++){
                    vec_v.push_back(subToGlobal[nbr_u[j]]);
                }
            }
            uint32_t u_pos = 0, v_pos = 0, count = 0;
            while(u_pos < vec_u.size() && v_pos < vec_v.size()){
                if(vec_u[u_pos] < vec_v[v_pos]) u_pos++;
                else if(vec_u[u_pos] > vec_v[v_pos]) v_pos++;
                else{
                    count++;
                    u_pos++, v_pos++;
                }
            }
            total_io += fDatSub.get_total_io();
            fDatSub.fclose();
            if(count < last_sup) continue; 

            memset(isInSubG, false, sizeof(bool)*nodeNum);
            for(int k = 0; k < nodeNum; k++)
                if(coreNum[k] >= last_sup+1)
                    isInSubG[k] = true;
            dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
            readFile subFile(dir);
            filterGlobalIntoSub(isInSubG,file,subFile,true); 
            globalToSub.clear();
            for(auto it = subToGlobal.begin(); it != subToGlobal.end(); it++)
                globalToSub[it->second] = it->first;
            fill(begPtr.begin(), begPtr.end(), 0);
            fill(degSub.begin(), degSub.end(), 0);

            Graph subG(subFile.verNum,subFile.edgeNum);
            subG_edgeNum = subFile.edgeNum, subG_nodeNum = subFile.verNum;
            printf("subFile.verNum: %u, subFile.edgeNum: %lu\n",subFile.verNum,subFile.edgeNum);
            subG.Initial(subFile);
            for(int k = 0; k < subG.nodeNum; k++){
                begPtr[k] = subG.edgeListBegPtr[k];
                degSub[k] = subG.degree[k];
            }
            subG.CountTriangleSSDByDegOrder(subFile,true);
            subG.trussDecomLazyUpdate(subFile);
            last_sup = subG.last_sup;
            dir = file.m_base + "graphInfoCopy/"+ to_string(last_sup) + "/";
            recoveryKMaxTrussSub(file,globalToSub,verMaxKTruSet,dir);
            maxTrussFile.change(dir);
            total_io += subG.total_io;
            dynamic[i].u = globalToSub[dynamic[i].u], dynamic[i].v = globalToSub[dynamic[i].v];
            printf("core: ori_u %u, now %u, ori_v %u, now %u\n",oriCoreU,coreNum[a],oriCoreV,coreNum[b]);
            insEdgeDynamic(dynamic[i],maxTrussFile,sup_dram,eid_dram,verMaxKTruSet,begPtr,degSub);
        }
    }

    total_io = total_io + fSubGInfo.get_total_io() + fDat.get_total_io();
    log_debug(graphClock_.Count("total_io: %lu\n",total_io));


    delete[] dynamic;
    delete[] m_delBit;
    delete[] m_insBit;
    delete[] m_dynamicDel;
    delete[] m_dynamicIns;
    delete[] coreNum;
    delete[] nbr_u;
    delete[] nbr_v;
    fSubGInfo.fclose();
    fDat.fclose();
    
}

void Graph::IntersectOperaInsDynamic(uint32_t u, uint32_t* nbr_u, uint32_t* sup_u, uint64_t* eid_u, uint32_t u_degree, 
uint32_t v, uint32_t* nbr_v, uint32_t* sup_v, uint64_t* eid_v, uint32_t v_degree,uint32_t &max_sup,std::vector<ver_eid_eid> &comm){
    uint32_t u_pos = 0, v_pos = 0;    
    while(u_pos < u_degree && v_pos < v_degree){
        if(nbr_u[u_pos] < nbr_v[v_pos]) u_pos++;
        else if(nbr_u[u_pos] > nbr_v[v_pos]) v_pos++;
        else{
            max_sup = max(max_sup,sup_u[u_pos]);
            max_sup = max(max_sup,sup_v[v_pos]);
            ver_eid_eid tmp;
            tmp.w = nbr_u[u_pos], tmp.first = eid_u[u_pos], tmp.second = eid_v[v_pos];
            tmp.u_sup = sup_u[u_pos], tmp.v_sup = sup_v[v_pos];
            comm.push_back(tmp);
            u_pos++, v_pos++;
        }
    }
}

void Graph::IntersectOperaDelDynamic(uint32_t u, uint32_t* nbr_u, uint32_t* sup_u, uint64_t* eid_u, uint32_t u_degree, 
uint32_t v, uint32_t* nbr_v, uint32_t* sup_v, uint64_t* eid_v, uint32_t v_degree,std::vector<ver_eid_eid> &comm){
    uint32_t u_pos = 0, v_pos = 0;  
    if((v_degree > 0 && u_degree > 0) && (nbr_u[u_degree-1] < nbr_v[0] || nbr_v[v_degree-1] < nbr_u[0]))
        return ;  
    while(u_pos < u_degree && v_pos < v_degree){
        if(MOVE(sup_u[u_pos]))
        {
            u_pos++;
            continue;
        }
        if(MOVE(sup_v[v_pos]))
        {
            v_pos++;
            continue;
        }
        if(nbr_u[u_pos] < nbr_v[v_pos]) u_pos++;
        else if(nbr_u[u_pos] > nbr_v[v_pos]) v_pos++;
        else{
            ver_eid_eid tmp;
            tmp.w = nbr_u[u_pos], tmp.first = eid_u[u_pos], tmp.second = eid_v[v_pos];
            tmp.u_sup = sup_u[u_pos], tmp.v_sup = sup_v[v_pos];
            comm.push_back(tmp);
            u_pos++, v_pos++;
        }
    }
}

void Graph::insEdgeDynamic(Edge e, readFile &file, unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> &sup_dram,
unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> &eid_dram,
bool *verMaxKTruSet, vector<uint64_t> &begPtr, vector<uint32_t> &degSub) 
{
    ofstream out("mytext.txt");
    // printf("insert u: %u, v: %u\n",e.u,e.v);
    MyReadFile fInfo(file.m_info);
    fInfo.fopen(BUFFERED);
    MyReadFile fIdx(file.m_idx);
	fIdx.fopen(BUFFERED);
	MyReadFile fDat(file.m_dat);
	fDat.fopen(BUFFERED); 
    MyReadFile fEid(file.m_eid);
	fEid.fopen(BUFFERED);
	MyReadFile fOff(file.m_offset);
    fOff.fopen(BUFFERED);
    MyReadFile fSup(file.m_supp);
	fSup.fopen(NOBUFFER); 
    uint32_t u = e.u, v = e.v, u_degree = 0, v_degree = 0, subGraV = 0;
    uint32_t* nbr_u = new uint32_t[nodeNum]();
    uint32_t* nbr_v = new uint32_t[nodeNum]();

    uint64_t* eid_u = new uint64_t[nodeNum]();
    uint64_t* eid_v = new uint64_t[nodeNum]();

    uint32_t* sup_u = new uint32_t[nodeNum]();
    uint32_t* sup_v = new uint32_t[nodeNum]();
    uint64_t pos = 0, posPlus = 0, eid = eid_dram[subToGlobal[u]][subToGlobal[v]];
    u_degree = degSub[u], v_degree = degSub[v];
    sup_dram[subToGlobal[u]][subToGlobal[v]] = 0;

    eid_eid tmp;
    fOff.fseek(eid*sizeof(eid_eid));
    fOff.fread(&tmp,sizeof(eid_eid));
    fSup.fseek(tmp.first*sizeof(uint32_t));
    fSup.fwrite(&u,sizeof(uint32_t));
    fSup.fseek(tmp.second*sizeof(uint32_t));
    fSup.fwrite(&u,sizeof(uint32_t));

    loadNbrAndSupDynamic(u, nbr_u, sup_u, eid_u, u_degree, begPtr[u], fDat, fSup, fEid, sup_dram, eid_dram);
    loadNbrAndSupDynamic(v, nbr_v, sup_v, eid_v, v_degree, begPtr[v], fDat, fSup, fEid, sup_dram, eid_dram);
    
    uint32_t kmin = 0, kmax = 0, max_sup = 0;
    vector<ver_eid_eid> comm;
    IntersectOperaInsDynamic(u,nbr_u,sup_u,eid_u,u_degree,v,nbr_v,sup_v,eid_v,v_degree,max_sup,comm);

    kmin = max_sup;
    while(kmin){
        uint32_t count = 0;
        for(int i = 0; i < comm.size(); i++)
            if(comm[i].u_sup >= kmin && comm[i].v_sup >= kmin)
                count++;
        if(count >= kmin)
            break;
        kmin--;
    }

    kmax = max_sup+1;
    while(kmax){
        uint32_t count = 0;
        for(int i = 0; i < comm.size(); i++)
            if(comm[i].u_sup >= kmax-1 && comm[i].v_sup >= kmax-1)
                count++;
        if(count >= kmax)
            break;
        kmax--;
    }

    struct cmp{
        bool operator()(EdgeSup a, EdgeSup b){
            return a.sup > b.sup;
        }
    };
    std::priority_queue<EdgeSup,std::vector<EdgeSup>, cmp> pq; 

    
    uint32_t k2 = kmax-1;
    if(k2 == maxKtruss)
        memset(verMaxKTruSet,0,sizeof(bool)*nodeNum);
    queue<TEdge>** que = new std::queue<TEdge>*[k2+1];
    queue<TEdge> tmp_q;
    for(int i = 0; i <= k2; i++)
        que[i] = new queue<TEdge>();
    unordered_map<uint32_t,unordered_map<uint32_t,bool>> isInLk;
    unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> support;
    unordered_map<uint64_t, uint32_t> contain;
    uint64_t combine = COMBINE(u,v);
    que[k2]->push({u,v,eid});
    log_debug(graphClock_.Count("kmin: %u, kmax: %u, eid: %lu",kmin,kmax,eid));
    contain[combine] = k2; 
    sup_dram[subToGlobal[u]][subToGlobal[v]] = k2;
    
    fOff.fseek(eid*sizeof(eid_eid));
    fOff.fread(&tmp,sizeof(eid_eid));
    fSup.fseek(tmp.first*sizeof(uint32_t));
    fSup.fwrite(&k2,sizeof(uint32_t));
    fSup.fseek(tmp.second*sizeof(uint32_t));
    fSup.fwrite(&k2,sizeof(uint32_t));

    isInLk[u][v] = true;
    for(int i = 0; i < comm.size(); i++){
        if(comm[i].u_sup <= k2){
            uint32_t x = min(u,comm[i].w);
            uint32_t y = max(u,comm[i].w);
            que[comm[i].u_sup]->push({x,y,comm[i].first});
            isInLk[x][y] = true;
            uint64_t combine = COMBINE(x,y);
            contain[combine] = comm[i].u_sup;
        }
        if(comm[i].v_sup <= k2){
            uint32_t x = min(v,comm[i].w);
            uint32_t y = max(v,comm[i].w);
            que[comm[i].v_sup]->push({x,y,comm[i].second});
            isInLk[x][y] = true;
            uint64_t combine = COMBINE(x,y);
            contain[combine] = comm[i].v_sup;
        }
        // printf("u: %u, v: %u, w: %u, supu: %u, supv: %u\n",u,v,comm[i].w,comm[i].u_sup,comm[i].v_sup); 
    }
    uint32_t c = 0;
    for(uint32_t k = k2; k >=2; k--){
        if(que[k]->size() == 0) continue;
        while(!que[k]->empty()){
            TEdge tmpe = que[k]->front();
            u = tmpe.u, v = tmpe.v, eid = tmpe.eid;
            que[k]->pop();
            support[u][v] = 0;
            u_degree = degSub[u], v_degree = degSub[v];
            loadNbrAndSupDynamic(u, nbr_u, sup_u, eid_u, u_degree, begPtr[u], fDat, fSup,fEid, sup_dram, eid_dram);
            loadNbrAndSupDynamic(v, nbr_v, sup_v, eid_v, v_degree, begPtr[v], fDat, fSup,fEid, sup_dram, eid_dram);
            comm.clear();
            IntersectOperaInsDynamic(u,nbr_u,sup_u,eid_u,u_degree,v,nbr_v,sup_v,eid_v,v_degree,max_sup,comm);
            for(int i = 0; i < comm.size(); i++){
                if(comm[i].u_sup < k || comm[i].v_sup < k) continue;
                support[u][v] += 1;
                // printf("u: %u, v: %u, w: %u, comm[i].u_sup: %u, comm[i].v_sup: %u\n",u,v,comm[i].w,comm[i].u_sup,comm[i].v_sup); 
                uint32_t x = min(u,comm[i].w);
                uint32_t y = max(u,comm[i].w);
                if(comm[i].u_sup == k){
                    if(isInLk.find(x) == isInLk.end() || isInLk[x].find(y) == isInLk[x].end()){
                        que[k]->push({x,y,comm[i].first});
                        isInLk[x][y] = true;
                        uint64_t combine = COMBINE(x,y);
                        contain[combine] = comm[i].u_sup;
                    }
                }
                x = min(v,comm[i].w);
                y = max(v,comm[i].w);
                if(comm[i].v_sup == k){
                    if(isInLk.find(x) == isInLk.end() || isInLk[x].find(y) == isInLk[x].end()){
                        que[k]->push({x,y,comm[i].second});
                        isInLk[x][y] = true;
                        uint64_t combine = COMBINE(x,y);
                        contain[combine] = comm[i].v_sup;
                    }
                }
            } 
            // out<<"u: " << u <<", v: " << v <<", sup: " <<support[u][v] << ",que.size: "<<que[k]->size() <<", pq.size: "<<pq.size()<<endl;
            // if(c >= 643000)
            //     printf("u: %u, v: %u, eid: %lu, sup: %u, que[k].size: %u, pq.size: %u, c: %u\n",u,v,eid,support[u][v],que[k]->size(),pq.size(),c);
            pq.push({u,v,support[u][v],eid});
            c++;
        }
        log_debug(graphClock_.Count("after for!"));
        while(!pq.empty()){
            EdgeSup es_tmp = pq.top(); 
            u = es_tmp.u, v = es_tmp.v, eid = es_tmp.eid;
            uint32_t sup = es_tmp.sup;
            pq.pop();
            if(isInLk.find(u) == isInLk.end() || isInLk[u].find(v) == isInLk[u].end()) continue;
            if(sup > k) {
                pq.push({u,v,sup,eid});
                break;
            }
            
            u_degree = degSub[u], v_degree = degSub[v];
            loadNbrAndSupDynamic(u, nbr_u, sup_u, eid_u, u_degree, begPtr[u], fDat, fSup,fEid, sup_dram, eid_dram);
            loadNbrAndSupDynamic(v, nbr_v, sup_v, eid_v, v_degree, begPtr[v], fDat, fSup,fEid, sup_dram, eid_dram);

            comm.clear();
            IntersectOperaInsDynamic(u,nbr_u,sup_u,eid_u,u_degree,v,nbr_v,sup_v,eid_v,v_degree,max_sup,comm);
            for(int i = 0; i < comm.size(); i++){
                if(comm[i].u_sup < k || comm[i].v_sup < k) continue;
                uint32_t x = min(u,comm[i].w);
                uint32_t y = max(u,comm[i].w);
                if(comm[i].u_sup == k){
                    if(isInLk.find(x) == isInLk.end() || isInLk[x].find(y) == isInLk[x].end())
                        continue;
                }
                uint32_t xx = min(v,comm[i].w);
                uint32_t yy = max(v,comm[i].w);
                if(comm[i].v_sup == k){
                    if(isInLk.find(xx) == isInLk.end() || isInLk[xx].find(yy) == isInLk[xx].end())
                        continue;  
                }
                if(contain.find(COMBINE(x,y)) != contain.end() && contain[COMBINE(x,y)] == k)
                    support[x][y] -= 1, pq.push({x,y,support[x][y],comm[i].first});
                if(contain.find(COMBINE(xx,yy)) != contain.end() && contain[COMBINE(xx,yy)] == k)
                    support[xx][yy] -= 1, pq.push({xx,yy,support[xx][yy],comm[i].second});
            }
            // printf("u: %u, v: %u, eid: %lu, sup: %u, size: %u\n",u,v,eid,support[u][v],comm.size());   

            isInLk[u].erase(v);
            support[u].erase(v);
            contain.erase(COMBINE(u,v));
        }
        log_debug(graphClock_.Count("after first while!"));

        
        if(pq.size() != 0){
            minSup = INF;
            while(!pq.empty()){
                EdgeSup es_tmp = pq.top();
                pq.pop();
                u = es_tmp.u, v = es_tmp.v, eid = es_tmp.eid;
                if(es_tmp.sup == maxKtruss)
                    verMaxKTruSet[u] = true, verMaxKTruSet[v] = true;
                // printf("u: %u, v: %u, sup: %u\n",u,v,es_tmp.sup);
                uint32_t sup = es_tmp.sup;
                if(sup != 0)
                    minSup = min(sup,minSup);
                if(isInLk.find(u) == isInLk.end() || isInLk[u].find(v) == isInLk[u].end()) continue;
                uint32_t u_global = subToGlobal[u], v_global = subToGlobal[v];
                if(m_insBit[u_global] && m_dynamicIns[u_global]->find(v_global) != m_dynamicIns[u_global]->end()){
                    sup_dram[u_global][v_global] += 1;
                }
                else{
                    eid_eid tmp;
                    fOff.fseek(eid*sizeof(eid_eid));
                    fOff.fread(&tmp,sizeof(eid_eid));
                    fSup.fseek(tmp.first*sizeof(uint32_t));
                    fSup.fwrite(&sup,sizeof(uint32_t));
                    fSup.fseek(tmp.second*sizeof(uint32_t));
                    fSup.fwrite(&sup,sizeof(uint32_t));
                }
                isInLk[u].erase(v);
                support[u].erase(v);
                contain.erase(COMBINE(u,v));
                if(isInLk[u].size() == 0)
                {
                    isInLk.erase(u);
                    support.erase(u);
                }
            }
            maxKtruss = minSup;
            
        }
    }
    

    log_debug(graphClock_.Count("maxKTruss: %d, insert successfully!\n", maxKtruss));
    total_io = total_io + fInfo.get_total_io() + fIdx.get_total_io() + fDat.get_total_io() + fSup.get_total_io() + fOff.get_total_io() + fEid.get_total_io();

    delete[] nbr_u;
    delete[] nbr_v;

    delete[] eid_u;
    delete[] eid_v;

    delete[] sup_u;
    delete[] sup_v;
    fInfo.fclose();
    fIdx.fclose();
    fDat.fclose();
    fSup.fclose();
    fOff.fclose();
    fEid.fclose();
}


void Graph::insEdgeDynamicNew(Edge e, readFile &file, unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> &sup_dram,
unordered_map<uint32_t,unordered_map<uint32_t,uint64_t>> &eid_dram, unordered_map<uint32_t,uint32_t> &subToGlobal) 
{
    // printf("insert u: %u, v: %u\n",e.u,e.v);
    MyReadFile fInfo(file.m_info);
    fInfo.fopen(BUFFERED);
    MyReadFile fIdx(file.m_idx);
	fIdx.fopen(BUFFERED);
	MyReadFile fDat(file.m_dat);
	fDat.fopen(BUFFERED); 
    MyReadFile fEid(file.m_eid);
	fEid.fopen(BUFFERED);
	MyReadFile fOff(file.m_offset);
    fOff.fopen(BUFFERED);
    MyReadFile fSup(file.m_supp);
	fSup.fopen(NOBUFFER); 
    uint32_t u = e.u, v = e.v, u_degree = 0, v_degree = 0, subGraV = 0;
    uint32_t* nbr_u = new uint32_t[nodeNum]();
    uint32_t* nbr_v = new uint32_t[nodeNum]();

    uint64_t* eid_u = new uint64_t[nodeNum]();
    uint64_t* eid_v = new uint64_t[nodeNum]();

    uint32_t* sup_u = new uint32_t[nodeNum]();
    uint32_t* sup_v = new uint32_t[nodeNum]();
    uint64_t pos = 0, posPlus = 0, eid = eid_dram[subToGlobal[u]][subToGlobal[v]];
    u_degree = degree[u], v_degree = degree[v];
    sup_dram[subToGlobal[u]][subToGlobal[v]] = 0;

    memset(vis,0,sizeof(bool)*nodeNum);
    memset(degree_,0,sizeof(bool)*nodeNum);

    eid_eid tmp;
    fOff.fseek(eid*sizeof(eid_eid));
    fOff.fread(&tmp,sizeof(eid_eid));
    fSup.fseek(tmp.first*sizeof(uint32_t));
    fSup.fwrite(&u,sizeof(uint32_t));
    fSup.fseek(tmp.second*sizeof(uint32_t));
    fSup.fwrite(&u,sizeof(uint32_t));

    loadNbrAndSupDynamicNew(u, nbr_u, sup_u, eid_u, u_degree, edgeListBegPtr[u], fDat, fSup, fEid);
    loadNbrAndSupDynamicNew(v, nbr_v, sup_v, eid_v, v_degree, edgeListBegPtr[v], fDat, fSup, fEid);
    
    uint32_t kmin = 0, kmax = 0, max_sup = 0;
    vector<ver_eid_eid> comm;
    IntersectOperaInsDynamic(u,nbr_u,sup_u,eid_u,u_degree,v,nbr_v,sup_v,eid_v,v_degree,max_sup,comm);
    uint32_t k2 = comm.size();


    unordered_map<uint32_t,unordered_map<uint32_t,bool>> isInLk;
    unordered_map<uint32_t,unordered_map<uint32_t,uint32_t>> support;
    unordered_map<uint64_t, uint32_t> contain;
    uint64_t combine = COMBINE(u,v);
    log_debug(graphClock_.Count("maxKtruss: %u",maxKtruss));

    changeEdgeSup(fOff,fSup,eid,k2);
    
    isInLk[u][v] = true;
    queue<uint32_t> tranvers;
    uint64_t po_edge = 0;
    uint32_t po_node = 0;
    k2 = 0;
    for(int i = 0; i < comm.size(); i++){
        comm[i].u_sup += 1;
        changeEdgeSup(fOff,fSup,comm[i].first,comm[i].u_sup);

        comm[i].v_sup += 1;
        changeEdgeSup(fOff,fSup,comm[i].second,comm[i].v_sup);

        if(comm[i].u_sup > maxKtruss){
            uint32_t x = min(u,comm[i].w);
            uint32_t y = max(u,comm[i].w);
            if(!vis[x])
                tranvers.push(x);
            if(!vis[y])
                tranvers.push(y);
            vis[x] = true, vis[y] = true;
            degree_[x]++,degree_[y]++;
        }
        if(comm[i].v_sup > maxKtruss){
            uint32_t x = min(v,comm[i].w);
            uint32_t y = max(v,comm[i].w);
            
            if(!vis[x])
                tranvers.push(x);
            if(!vis[y])
                tranvers.push(y);
            vis[x] = true, vis[y] = true;
            degree_[x]++,degree_[y]++;
        }
        if(comm[i].v_sup > maxKtruss && comm[i].u_sup > maxKtruss)
            k2++;
        // printf("u: %u, v: %u, w: %u, supu: %u, supv: %u\n",u,v,comm[i].w,comm[i].u_sup,comm[i].v_sup); 
    }
    if(k2 == maxKtruss){
        log_debug(graphClock_.Count("maxKTruss: %d, insert successfully, not increase!", maxKtruss));
        total_io = total_io + fInfo.get_total_io() + fIdx.get_total_io() + fDat.get_total_io() + fSup.get_total_io() + fOff.get_total_io() + fEid.get_total_io();
        return;
    }
    
    uint32_t len = 100000000, index = 0;
    Edge *buf = new Edge[len];
    memset(buf,0,sizeof(Edge)*len);
    
    unordered_set<uint64_t> isInBuf;
    while(!tranvers.empty()){
        u = tranvers.front();
        tranvers.pop();
        u_degree = degree[u];
        loadNbrAndSupDynamicNew(u, nbr_u, sup_u, eid_u, u_degree, edgeListBegPtr[u], fDat, fSup, fEid);
        for(int i = 0; i < u_degree; i++)
        {
            v = nbr_u[i];
            if(sup_u[i] > maxKtruss){
                uint32_t x = min(u,v);
                uint32_t y = max(u,v);
                if(!vis[v])
                    tranvers.push(v);
                vis[v] = true;
                degree_[u]++,degree_[v]++;
                po_edge++;
            }
            else{
                uint32_t x = min(u,v), y = max(u,v);   
                uint64_t combine = COMBINE(x,y);
                if(isInBuf.find(combine) == isInBuf.end())
                    buf[index++] = {x,y}, isInBuf.insert(combine);
            }
        }
    }
    po_edge /= 2;
    for(int i = 0; i < nodeNum; i++)
        if(vis[i])
            po_node++;
    log_debug(graphClock_.Count("po_node: %u, po_edge: %u, index: %u",po_node,po_edge,index));
    if((maxKtruss+2)*(maxKtruss+3)/2 > po_edge)
        return;
    bool res = true;
    for(int j = 0; j < index; j++){
        // printf("index: %u, po_edge: %u, j: %u, u: %u, v: %u\n",index,po_edge,j,u,v);
        if((maxKtruss+2)*(maxKtruss+3)/2 > po_edge ) {
            res = false;
            break;
        }
        Edge tmp_edge = buf[j];
        u = tmp_edge.u, v = tmp_edge.v;
        uint64_t combine = COMBINE(u,v);

        u_degree = degree[u], v_degree = degree[v];
        loadNbrAndSupDynamicNew(u, nbr_u, sup_u, eid_u, u_degree, edgeListBegPtr[u], fDat, fSup, fEid);
        loadNbrAndSupDynamicNew(v, nbr_v, sup_v, eid_v, v_degree, edgeListBegPtr[v], fDat, fSup, fEid);
        vector<ver_eid_eid> comm;
        uint32_t u_pos = 0, v_pos = 0;    

        while(u_pos < u_degree && v_pos < v_degree){
            if(nbr_u[u_pos] < nbr_v[v_pos]) u_pos++;
            else if(nbr_u[u_pos] > nbr_v[v_pos]) v_pos++;
            else{
                if(sup_u[u_pos] <= maxKtruss && sup_v[v_pos] <= maxKtruss) {
                    u_pos++, v_pos++;
                    continue;
                }
                ver_eid_eid tmp;
                tmp.w = nbr_u[u_pos], tmp.first = eid_u[u_pos], tmp.second = eid_v[v_pos];
                tmp.u_sup = sup_u[u_pos], tmp.v_sup = sup_v[v_pos];
                comm.push_back(tmp);
                u_pos++, v_pos++;
            }
        }
        for(int i = 0; i < comm.size(); i++){  

            if(comm[i].u_sup > maxKtruss){
                if(contain.find(comm[i].first) == contain.end()) contain[comm[i].first] = comm[i].u_sup;
                comm[i].u_sup -= 1;
                changeEdgeSup(fOff,fSup,comm[i].first,comm[i].u_sup);
                if(comm[i].u_sup == maxKtruss) {
                    uint32_t x = min(u,comm[i].w), y = max(u,comm[i].w); 
                    uint64_t cb = COMBINE(x,y);
                    if(isInBuf.find(cb) == isInBuf.end()){
                        buf[index++] = {x,y};
                        isInBuf.insert(cb);
                    }
                    po_edge--;
                }
            }
            
            if(comm[i].v_sup > maxKtruss){
                if(contain.find(comm[i].second) == contain.end()) contain[comm[i].second] = comm[i].v_sup;
                comm[i].v_sup -= 1;
                changeEdgeSup(fOff,fSup,comm[i].second,comm[i].v_sup);
                if(comm[i].v_sup == maxKtruss) {
                    uint32_t x = min(v,comm[i].w), y = max(v,comm[i].w); 
                    uint64_t cb = COMBINE(x,y);
                    if(isInBuf.find(cb) == isInBuf.end()){
                        buf[index++] = {x,y};
                        isInBuf.insert(cb);
                    }
                    po_edge--;
                }
            }
        }
    }
    if(res)
        maxKtruss++;
    else
    {
        for(auto it = contain.begin(); it != contain.end(); it++){
            uint64_t eid_ = it->first;
            uint32_t sup_ = it->second;
            changeEdgeSup(fOff,fSup,eid_,sup_);
        }
        log_debug(graphClock_.Count("not success, contain.size(): %u",contain.size()));
    }

    log_debug(graphClock_.Count("maxKTruss: %d, po_edge : %u, insert successfully!\n", maxKtruss,po_edge));
    total_io = total_io + fInfo.get_total_io() + fIdx.get_total_io() + fDat.get_total_io() + fSup.get_total_io() + fOff.get_total_io() + fEid.get_total_io();

    delete[] nbr_u;
    delete[] nbr_v;

    delete[] eid_u;
    delete[] eid_v;

    delete[] sup_u;
    delete[] sup_v;
    fInfo.fclose();
    fIdx.fclose();
    fDat.fclose();
    fSup.fclose();
    fOff.fclose();
    fEid.fclose();
}

void Graph::changeEdgeSup(MyReadFile& fOff, MyReadFile& fSup, uint64_t eid, uint32_t& sup){
    eid_eid tmp;
    fOff.fseek(eid*sizeof(eid_eid));
    fOff.fread(&tmp,sizeof(eid_eid));
    fSup.fseek(tmp.first*sizeof(uint32_t));
    fSup.fwrite(&sup,sizeof(uint32_t));
    fSup.fseek(tmp.second*sizeof(uint32_t));
    fSup.fwrite(&sup,sizeof(uint32_t));
}

void Graph::generateRandomEdgesInsertion(readFile &file, uint32_t &generate_edge, Edge *dynamic)
{
    srand(time(0));
    /* Seed */ 
    std::random_device rd; 
    /* Random number generator */ 
    std::default_random_engine generator(rd()); 
    /* Distribution on which to apply the generator */ 
    std::uniform_int_distribution<long long unsigned> distribution(0,0xFFFFFFFFFFFFFFFF);

    unordered_set<uint32_t> exist;

    MyReadFile fIdx(file.m_idx);
	fIdx.fopen(BUFFERED);
	MyReadFile fDat(file.m_dat);
	fDat.fopen(BUFFERED);
	uint32_t* nbr = new uint32_t[file.maxDeg];
	uint32_t d;


    for(int i = 0; i < generate_edge; ){
        uint32_t start = distribution(generator) % nodeNum;
        if(exist.find(start) == exist.end()){
            #ifdef DegSort
            fIdx.fseek(start*(sizeof(long)*2 + sizeof(uint32_t)));
            #else
            fIdx.fseek(start*(sizeof(long) + sizeof(uint32_t)));
            #endif

            long pos, posPlus;
            fIdx.fread(&pos, sizeof(long));
            fIdx.fread(&posPlus, sizeof(long));
            fIdx.fread(&d, sizeof(uint32_t));
            fDat.fseek(pos);
            uint32_t u_nbrNum = d-(posPlus-pos)/sizeof(uint32_t);
            if(u_nbrNum == 0)
                continue;
            
            loadInfo(nbr,u_nbrNum,posPlus,fDat);
            uint32_t r = nbr[0], end = 0;

            for(uint32_t j = 1; j < u_nbrNum; j++){
                if(nbr[j] - r > 1){
                    end = r+1;
                    break;
                }
                r = nbr[j];
            }
            if(end == 0 && nbr[u_nbrNum-1] < nodeNum-1)
                end = nbr[u_nbrNum-1]+1;
            else if(end == 0 && nbr[u_nbrNum-1] == nodeNum-1)
                continue;
            dynamic[i] = {start,end};
            i++;
            exist.insert(start);
        }
    }

    std::sort(dynamic,dynamic+generate_edge,[](const Edge & a, const Edge & b) {
                return a.u < b.u;
            });
    delete[] nbr;
	fDat.fclose();
	fIdx.fclose();
}
