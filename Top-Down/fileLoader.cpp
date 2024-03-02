#include "include/edge.h"
#include "include/tools.h"
#include "include/fileLoader.h"
#include <cassert>


void readFile::createDir(string &sub_dir){

    if (0 != access(sub_dir.c_str(), 0)){
        int isCreate = mkdir(sub_dir.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
        if( !isCreate )
            printf("create path:%s\n",sub_dir.c_str());
        else
            printf("create path failed! error code : %s \n",sub_dir.c_str());
    }
}


//divide and conquer for sorting a extreme file
//ordered internal information for subfile
//merge many subfiles into an ordered extreme file 
void readFile::moduleInSaveEdges(uint32_t &u,uint32_t &v,int &num,TEdge* edges,
        unsigned long  &size,unsigned long &es,int &tmpFile,string &sub_dir){
    u = getVertexID(u,num);
    v = getVertexID(v,num);
    verDegMap[u]++;
    verDegMap[v]++;
    edges[size].u = u;
    edges[size].v = v;
    edges[size].eid = es;
    size++;

    edges[size].u = v;
    edges[size].v = u;
    edges[size].eid = es;
    size++;

    ++es;
    char fileName[150];

    if(size >= memEdges){         
        sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpFile);
        saveTmpEdges<TEdge>(edges,size,fileName,[](const TEdge & a, const TEdge & b) {
            if(a.u < b.u)
                return true;
            if( a.u > b.u )
                return false;
            return a.v < b.v;
            });
        size = 0;
        ++tmpFile;
    }
}


void readFile::moduleInSaveEdgesUnOrder(uint32_t &u,uint32_t &v,int &num,TEdge* edges,
        unsigned long  &size,uint64_t es,int &tmpFile,string &sub_dir){
    // u = getVertexID(u,num);
    // v = getVertexID(v,num);
    edges[size].u = u;
    edges[size].v = v;
    edges[size].eid = es;
    size++;

    edges[size].u = v;
    edges[size].v = u;
    edges[size].eid = es;
    size++;
    num++;
    char fileName[150];

    if(size >= memEdges){         
        sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpFile);
        saveTmpEdges<TEdge>(edges,size,fileName,[](const TEdge & a, const TEdge & b) {
            if(a.u < b.u)
                return true;
            if( a.u > b.u )
                return false;
            return a.v < b.v;
            });
        size = 0;
        ++tmpFile;
    }
}



//load graph data into binary file about graph

void readFile::LoadGraph(){
    char *bptr = getAddr();
    uint64_t len = getLen();
    uint32_t vertexNum = 0, maxDegree = 0;
    char fileName[200];

    const char *cptr = bptr;
    const char *const eptr = bptr + len;

    unsigned long size = 0,es = 0;
	int num = 0,tmpFile = 0;
    uint32_t u,v;
    memset(m_vertexMap,-1,sizeof(int)*m_maxID);
    TEdge* edges = new TEdge[memEdges];
    string name = "sort_edge_tmp";

    string sub_dir = m_base + name;
    createDir(sub_dir);

    while (cptr < eptr)
    {
        u = mysscanf(eptr, cptr);
        if (cptr >= eptr)
            break;
        v = mysscanf(eptr, cptr);
        if (cptr >= eptr)
            break;
        //  自环判断
        if (u == v)
            continue;
        
        moduleInSaveEdges(u,v,num,edges,size,es,tmpFile,sub_dir);
    }
    
    sprintf(fileName,"%s/edges_tmp_%d",sub_dir.c_str(),tmpFile);
	saveTmpEdges<TEdge>(edges,size,fileName,[](const TEdge & a, const TEdge & b) {
                if(a.u < b.u)
                    return true;
                if( a.u > b.u )
                    return false;
                return a.v < b.v;
                });
    log_debug(readFileClock_.Count("Save tmp edges done, load %ld edges.",es));
    edgeNum = es;
	delete[] edges;
    #ifdef DegSort
    mergeByDegSort(tmpFile+1, vertexNum, maxDegree, name, true);
    #else
	merge(tmpFile+1, vertexNum, maxDegree,name, true);
    #endif
    
}

// template<typename T> 
// void readFile::mergeBySup<T>(int size)
// {
//     FILE **frl = new FILE*[size];
//     FILE *fEdgeSup = fopen(m_suppSort.c_str(),"wb");
//     char filename[100];
//     EdgeSup tmp;

//     struct cmp{
//         bool operator()(pair<int,EdgeSup> a, pair<int,EdgeSup> b){
//             return a.second.sup > b.second.sup;
//         }
//     };

//     std::priority_queue<pair<int,EdgeSup>,std::vector<pair<int,EdgeSup>>, cmp> pq; 

//     string sub_dir = m_base + "support_sort_edge";

//     for(int i = 0; i < size; i++){
//         sprintf(filename,"%s/edges_tmp_%d",sub_dir.c_str(),i);
//         frl[i] = fopen(filename,"rb");
//         fread(&tmp,sizeof(EdgeSup),1,frl[i]);
//         pq.push({i,tmp});
//     }
//     int count =  0;
//     while(!pq.empty()){
//         int index = pq.top().first;
//         tmp = pq.top().second;
//         pq.pop();
//         // printf("count: %d, u: %u, v: %u, sup: %u\n",count,tmp.u,tmp.v,tmp.sup);
        
//         fwrite(&tmp,sizeof(EdgeSup),1,fEdgeSup);
//         count++;
//         if(fread(&tmp,sizeof(EdgeSup),1,frl[index])){
// 			pq.push({index,tmp});
// 		}
//     }

//     for (int i = 0; i < size ; ++i)
// 		fclose(frl[i]);
// 	fclose(fEdgeSup);
// 	delete[] frl;
// }



//merge data of sort_edge_tmp dir 
void readFile::merge(int size, uint32_t &vertexNum, uint32_t &maxDegree, string &name, bool flag){
    // if flag is true, graph vertex is reorder, otherwise, graph vertex is unorder
    FILE** frl = new FILE*[size];
    TEdge* es = new TEdge[size];
    FILE* fIdx = fopen(m_idx.c_str(),"wb");
    FILE* fEid = fopen(m_eid.c_str(),"wb");
	FILE* fDat = fopen(m_dat.c_str(),"wb");

    log_info(readFileClock_.Count("Start merge"));

    for(int i = 0; i < size; i++){
        char filename[100];
        sprintf(filename,"%s%s/edges_tmp_%d",m_base.c_str(),name.c_str(),i);
        frl[i] = fopen(filename,"rb");
        // get first edge of all edges_tmp files
        fread(&es[i],sizeof(TEdge),1,frl[i]);
    }
    int minIndex,previousA = -1,previousB = -1;
    uint32_t maxDegVer = 0;

	int degree = -1;
	int f = 0;
    long pos,last_pos;
    
    bool sig = false;
    while(mergeFinished(es,size)){
        minIndex = minPartition(es,size);
        uint32_t u = es[minIndex].u;
		uint32_t v = es[minIndex].v;
        uint64_t eid = es[minIndex].eid;        
		
		if(u != previousA){
			// u != previousA demonstrates that all previousA's neighbors have been writen
            if(!flag)//for unorder method
                vertexId.push_back(u);

			if(degree != -1){
				// write the vertex degree in .idx file
                fwrite(&degree,sizeof(int),1,fIdx);
                if(degree>maxDegree)
                {
                    maxDegree = degree;
                    maxDegVer = previousA;
                }
			}

			// write start position of neighbor of vertex u from .dat file to .idx file
			pos = ftell(fDat);
            last_pos = pos;
            sig = false;
			fwrite(&pos,sizeof(long),1,fIdx);
			degree = 1;
			fwrite(&v,sizeof(int),1,fDat);
            fwrite(&eid,sizeof(uint64_t),1,fEid);
            
		}
        else if(v != previousB){
			fwrite(&v,sizeof(int),1,fDat);  //save neighbor v of vertex u
            fwrite(&eid,sizeof(uint64_t),1,fEid);
			++degree;
		}
		// if u==previousA & v==previousB, ignore edge(u,v) cause it is same as previous one.
		previousA = u;
		previousB = v;

		// replace es[minIndex] by picking up the first edge from file edges_tmp_minIndex
		if(!fread(&es[minIndex],sizeof(TEdge),1,frl[minIndex])){
			es[minIndex].u = m_maxID;
		}
    }
    fwrite(&degree,sizeof(int),1,fIdx);
    if(degree>maxDegree)
    {
        maxDegree = degree;
        maxDegVer = previousA;
    }
    // maxDegree = degree > maxDegree ? degree : maxDegree;
    // maxDegVer = degree > maxDegree ? previousA : maxDegVer;

    log_info(readFileClock_.Count("Merge operation done"));
    
    // write the vertex num and max degree
    if(flag)
    {
        vertexNum = previousA+1;
        verNum = vertexNum;
    }
    else
        verNum = vertexId.size();
	    
    log_debug(readFileClock_.Count("vertex num: %d",vertexNum));
    log_debug(readFileClock_.Count("Max degree: %d",maxDegree));
    log_debug(readFileClock_.Count("Max degree vertex: %d",maxDegVer));

    maxDeg = maxDegree;
    

	FILE* fInfo = fopen(m_info.c_str(),"wb");
	fwrite(&verNum,sizeof(int),1,fInfo);
	fwrite(&maxDegree,sizeof(int),1,fInfo);
    fwrite(&edgeNum,sizeof(uint64_t),1,fInfo);
	fclose(fInfo);

	for (int i = 0; i < size ; ++i){
		fclose(frl[i]);
	}

	fclose(fIdx);
	fclose(fDat);
    fclose(fEid);

	delete[] frl;
	delete[] es;
}


template<typename T>
struct comp{
        bool operator()(pair<T,uint32_t> a, pair<T,uint32_t> b){
            return a.second > b.second;
        }
};

//merge data of sort_edge_tmp dir in degree increasing method 
void readFile::mergeByDegSort(int size, uint32_t &vertexNum, uint32_t &maxDegree,string &name, bool flag){
    // if flag is true, graph vertex is reorder, otherwise, graph vertex is unorder

    FILE** frl = new FILE*[size];
    TEdge* es = new TEdge[size];
    FILE* fIdx = fopen(m_idx.c_str(),"wb");
    FILE* fEid = fopen(m_eid.c_str(),"wb");
	FILE* fDat = fopen(m_dat.c_str(),"wb");
    FILE* fOff = fopen(m_offset.c_str(),"wb");

    log_info(readFileClock_.Count("Start merge"));    

    /*  vertex ID is not need to reorder by degree increasing order */

    for(int i = 0; i < size; i++){
        char filename[100];
        sprintf(filename,"%s%s/edges_tmp_%d",m_base.c_str(),name.c_str(),i);
        frl[i] = fopen(filename,"rb");
        // get first edge of all edges_tmp files
        fread(&es[i],sizeof(TEdge),1,frl[i]);
    }
    int minIndex,previousA = -1,previousB = -1;
    uint32_t maxDegVer = 0;

	int degree = -1, f = 0;
    long pos,last_pos;
    uint64_t off = 0;
    
    bool sig = false;
    std::unordered_map<uint64_t,uint64_t> eidToPos;
    std::unordered_map<uint64_t,uint64_t> reorderEid;
    uint64_t max_size = 0,index = 0, eidSize = 0;
    uint32_t buffersize = 100*1024*1024;
    uint32_t *buf = new uint32_t[buffersize]();
    uint64_t *buf_eid = new uint64_t[buffersize]();

    while(mergeFinished(es,size)){
        minIndex = minPartition(es,size);
        uint32_t u = es[minIndex].u;
		uint32_t v = es[minIndex].v;
        uint64_t eid = es[minIndex].eid;  
        uint64_t com = COMBINE(u,v);
        if((u == previousA && v == previousB)) {
            if(!fread(&es[minIndex],sizeof(TEdge),1,frl[minIndex]))
			    es[minIndex].u = m_maxID;
            continue; 
        }  
        if(eidToPos.find(com) != eidToPos.end()){
            eid = reorderEid[com];
            eid_eid tmp = {eidToPos[com],off};
            // random write costs a lot
            fseek(fOff,eid*sizeof(eid_eid),SEEK_SET);
            fwrite(&tmp,sizeof(eid_eid),1,fOff);
            eidToPos.erase(com);
            reorderEid.erase(com);            
        }   
        else{
            eidToPos[COMBINE(v,u)] = off;
            reorderEid[COMBINE(v,u)] = eidSize;
            eid = eidSize++;
        } 
        // printf("u: %u, v: %u, eid: %lu\n",u,v,eid);
        off++;
        max_size = std::max(eidToPos.size(),max_size);
		
		if(u != previousA){
			// u != previousA demonstrates that all previousA's neighbors have been writen
            if(!flag)//for unorder method
                vertexId.push_back(u);

			if(degree != -1){
				// write the vertex degree in .idx file
                if(!sig && last_pos == pos){
                    last_pos = pos+degree*sizeof(uint32_t);
                    fwrite(&last_pos,sizeof(long),1,fIdx);
                }
                fwrite(&degree,sizeof(int),1,fIdx);
                if(degree>maxDegree)
                {
                    maxDegree = degree;
                    maxDegVer = previousA;
                }
			}

			// write start position of neighbor of vertex u from .dat file to .idx file
			pos = ftell(fDat)+index*sizeof(uint32_t);
            last_pos = pos;
            sig = false;
			fwrite(&pos,sizeof(long),1,fIdx);
            if(!sig && v > u ){
                fwrite(&last_pos,sizeof(long),1,fIdx);
                sig = true;
            }
			degree = 1;

            if(index >= buffersize){
                fwrite(buf,sizeof(uint32_t),buffersize,fDat);
                fwrite(buf_eid,sizeof(uint64_t),buffersize,fEid);
                index = 0;
            }
            else
            {
                buf[index] = v;
                buf_eid[index++] = eid;
            }
			// fwrite(&v,sizeof(int),1,fDat);
            // fwrite(&eid,sizeof(uint64_t),1,fEid);
		}
        else if(v != previousB){
            if(!sig && v > previousA ){
                last_pos = pos+degree*sizeof(uint32_t);
                fwrite(&last_pos,sizeof(long),1,fIdx);
                sig = true;
            }
            if(index >= buffersize){
                fwrite(buf,sizeof(uint32_t),buffersize,fDat);
                fwrite(buf_eid,sizeof(uint64_t),buffersize,fEid);
                index = 0;
            }
            else
            {
                buf[index] = v;
                buf_eid[index++] = eid;
            }

			++degree;
		}
		// if u==previousA & v==previousB, ignore edge(u,v) cause it is same as previous one.
		previousA = u;
		previousB = v;

		// replace es[minIndex] by picking up the first edge from file edges_tmp_minIndex
		if(!fread(&es[minIndex],sizeof(TEdge),1,frl[minIndex])){
			es[minIndex].u = m_maxID;
		}
    }
    if(index)
    {
        fwrite(buf,sizeof(uint32_t),index,fDat);
        fwrite(buf_eid,sizeof(uint64_t),index,fEid);
    }
    if(!sig && last_pos == pos){
        last_pos = pos + degree*sizeof(uint32_t);
        fwrite(&last_pos,sizeof(long),1,fIdx);
    }
    fwrite(&degree,sizeof(int),1,fIdx);
    if(degree>maxDegree)
    {
        maxDegree = degree;
        maxDegVer = previousA;
    }
    log_info(readFileClock_.Count("Merge operation done"));
    
    // write the vertex num and max degree
    if(flag)
    {
        vertexNum = previousA+1;
        verNum = vertexNum;
    }
    else
        verNum = vertexId.size();
    
    edgeNum = off/2;
	    
    log_debug(readFileClock_.Count("vertex num: %d, edge num: %d",vertexNum, edgeNum));
    log_debug(readFileClock_.Count("Max degree: %d, Max degree vertex: %d",maxDegree, maxDegVer));
    log_debug(readFileClock_.Count("Max map size: %u",max_size));

    maxDeg = maxDegree;
    

	FILE* fInfo = fopen(m_info.c_str(),"wb");
	fwrite(&verNum,sizeof(int),1,fInfo);
	fwrite(&maxDegree,sizeof(int),1,fInfo);
    fwrite(&edgeNum,sizeof(uint64_t),1,fInfo);
	fclose(fInfo);

	for (int i = 0; i < size ; ++i){
		fclose(frl[i]);
	}

	fclose(fIdx);
	fclose(fDat);
    fclose(fEid);
    fclose(fOff);

	delete[] frl;
	delete[] es;
    delete[] buf;
    delete[] buf_eid;
}



// get minimal edge from edge list
int readFile::minPartition(TEdge* es,int size){
	int min = 0;
    for(int i = 1; i < size; i++){
        if(es[i].u < es[min].u){
            min = i;
        }else if(es[i].u > es[min].u){
            continue;
        }else if(es[i].v < es[min].v){
            min = i;
        }
    }
    return min;
}


bool readFile::mergeFinished(TEdge* es, int size){
    for(int i = 0 ; i < size; i++){
        if(es[i].u != m_maxID)
            return true;
    }
    return false;
}

uint32_t readFile::getVertexID(uint32_t u,int &num){
	if(m_vertexMap[u] == INF_U32){
		m_vertexMap[u] = num++;
	}
	return m_vertexMap[u];
}


bool readFile::edgeCompare(const TEdge &e1, const TEdge &e2){
	if(e1.u < e2.u){
		return true;
	}
	if( e1.u > e2.u ){
		return false;
	}
	return e1.v < e2.v;
}
// sort edges and save
// template<typename T> 
// void readFile::saveTmpEdges(T* edges, int size, char fileName[], std::function<bool(const T &, const T &)> cmp){
// 	// std::sort(edges,edges+size,[&](int i, int j) {
//     //     if(edges[i].u < edges[j].u){
//     //         return true;
//     //     }
//     //     if( edges[i].u > edges[j].u ){
//     //         return false;
//     //     }
//     //     return edges[i].v < edges[j].v;
//     // });
//     std::sort(edges,edges+size,cmp);
// 	// char fileName[200];
// 	// sprintf(fileName,"%ssort_edge_tmp/edges_tmp_%d",m_base.c_str(),tmpFile);

//     log_debug(readFileClock_.Count(fileName));

// 	FILE* fo = fopen(fileName,"wb");
// 	// for (int i = 0; i < size; ++i)
//     // {
// 	// 	fwrite( edges+i, sizeof(TEdge), 1, fo );
// 	// 	// printf("edge[%d,%d]\n",edges[i].a,edges[i].b );
// 	// }
//     fwrite( edges, sizeof(T), size, fo );
// 	fclose(fo);
// }



uint32_t LoadEdge(const char *const bptr, const uint32_t len, std::vector<Edge> &Edges)
{
    const char *cptr = bptr;
    const char *const eptr = bptr + len;
    Edges.clear();
    std::set<uint32_t> vertex;
    std::map<uint32_t,int> mp;

    while (cptr < eptr)
    {

        struct Edge curEdge;
        curEdge.u = mysscanf(eptr, cptr);
        if (cptr >= eptr)
            break;
        curEdge.v = mysscanf(eptr, cptr);
        if (cptr >= eptr)
            break;
        vertex.insert(curEdge.u);
        vertex.insert(curEdge.v);
        //  自环判断
        if (curEdge.u == curEdge.v)
            continue;
        Edges.push_back(curEdge);

        // sscanfSkip(eptr, cptr);
    }
    uint32_t c = 0;
    for(auto it = vertex.begin(); it != vertex.end(); it++)
    {
        mp[*it] = c++;
    }

    #pragma omp parallel for schedule(static)
    for(uint64_t i = 0; i < Edges.size(); i++){
        uint32_t u = mp[Edges[i].u];
        uint32_t v = mp[Edges[i].v];

        if(v < u)
            std::swap(u,v);
        Edges[i].u = u;
        Edges[i].v = v;
    }
    return vertex.size();
}

uint32_t mysscanf(const char *const eptr, const char *&cptr)
{
    uint32_t ret = 0;
    bool readF = false;
    while (cptr < eptr)
    {
        const char ch = *cptr;
        if (ch >= '0' && ch <= '9')
        {
            ret = ret * 10 + (ch - '0');
            readF = true;
        }
        else
        {
            if (readF)
                return ret;
        }
        cptr++;
    }
    if (readF)
        return ret;
    else
        return INF;
}

void sscanfSkip(const char *const eptr, const char *&cptr)
{
    while (*cptr != '1' && cptr < eptr)
        cptr++;
    cptr++;
}
