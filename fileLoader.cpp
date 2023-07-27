#include "include/edge.h"
#include "include/tools.h"
#include "include/fileLoader.h"
#include <cassert>


void readFile::createDir(string &sub_dir){

    if (0 != access(sub_dir.c_str(), 0)){
        int isCreate = mkdir(sub_dir.c_str(),S_IRUSR | S_IWUSR | S_IXUSR | S_IRWXG | S_IRWXO);
        // if( !isCreate )
        //     printf("create path:%s\n",sub_dir.c_str());
        // else
        //     printf("create path failed! error code : %s \n",sub_dir.c_str());
    }
}

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


void readFile::moduleInSaveEdgesUnOrder(uint32_t &u, uint32_t &v, int &num, TEdge* edges,
        unsigned long  &size,uint32_t es,int &tmpFile,string &sub_dir){
    // u = getVertexID(u,num);
    // v = getVertexID(v,num);
    // printf("u: %u, v: %u, eid: %u\n",u,v,es);
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
    mergeByDegSort(tmpFile+1, vertexNum, maxDegree, name, true, true);
    #else
	merge(tmpFile+1, vertexNum, maxDegree,name, true);
    #endif

}

/*
void readFile::mergeBySup(int size)
{
    FILE **frl = new FILE*[size];
    FILE *fEdgeSup = fopen(m_suppSort.c_str(),"wb");
    char filename[100];
    EdgeSup tmp;

    struct cmp{
        bool operator()(pair<int,EdgeSup> a, pair<int,EdgeSup> b){
            return a.second.sup > b.second.sup;
        }
    };
    std::priority_queue<pair<int,EdgeSup>,std::vector<pair<int,EdgeSup>>, cmp> pq; 

    string sub_dir = m_base + "support_sort_edge";

    for(int i = 0; i < size; i++){
        sprintf(filename,"%s/edges_tmp_%d",sub_dir.c_str(),i);
        frl[i] = fopen(filename,"rb");
        fread(&tmp,sizeof(EdgeSup),1,frl[i]);
        pq.push({i,tmp});
    }
    write_io += size;
    while(!pq.empty()){
        int index = pq.top().first;
        tmp = pq.top().second;
        pq.pop();
        // printf("count: %d, u: %u, v: %u, sup: %u\n",count,tmp.u,tmp.v,tmp.sup);
        
        fwrite(&tmp,sizeof(EdgeSup),1,fEdgeSup);
        write_io += 1;
        if(fread(&tmp,sizeof(EdgeSup),1,frl[index])){
			pq.push({index,tmp});
            write_io += 1;
		}
    }

    for (int i = 0; i < size ; ++i)
		fclose(frl[i]);
	fclose(fEdgeSup);
	delete[] frl;
}
*/



//merge data of sort_edge_tmp dir 
void readFile::merge(int size, uint32_t &vertexNum, uint32_t &maxDegree, string &name, bool flag){
    // if flag is true, graph vertex is reorder, otherwise, graph vertex is unorder
    FILE** frl = new FILE*[size];
    TEdge* es = new TEdge[size];
    FILE* fIdx = fopen(m_idx.c_str(),"wb");
    FILE* fEid = fopen(m_eid.c_str(),"wb");
	FILE* fDat = fopen(m_dat.c_str(),"wb");
    FILE* fOff = fopen(m_offset.c_str(),"wb");
    uint32_t *buf = new uint32_t[buffersize]();
    uint64_t *buf_eid = new uint64_t[buffersize]();

    log_info(readFileClock_.Count("Start merge"));

    for(int i = 0; i < size; i++){
        char filename[100];
        sprintf(filename,"%s%s/edges_tmp_%d",m_base.c_str(),name.c_str(),i);
        frl[i] = fopen(filename,"rb");
        // get first edge of all edges_tmp files
        fread(&es[i],sizeof(TEdge),1,frl[i]);
    }
    write_io += size;
    int minIndex,previousA = -1,previousB = -1;
    uint32_t maxDegVer = 0,index = 0;
    uint64_t max_size = 0;
	int degree = -1;
    long pos,last_pos;


    std::unordered_map<uint64_t,uint64_t> eidToPos;
    uint64_t off = 0;

    while(mergeFinished(es,size)){
        minIndex = minPartition(es,size);
        uint32_t u = es[minIndex].u;
		uint32_t v = es[minIndex].v;
        uint64_t eid = es[minIndex].eid;    
        if(eidToPos.find(eid) != eidToPos.end()){
            eid_eid tmp = {eidToPos[eid],off};
            // random write costs a lot
            fseek(fOff,eid*sizeof(eid_eid),SEEK_SET);
            fwrite(&tmp,sizeof(eid_eid),1,fOff);
            eidToPos.erase(eid);
        }   
        else{
            eidToPos[eid] = off;
        } 
        off++;
        max_size = std::max(eidToPos.size(),max_size);     
		
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
			pos = ftell(fDat)+index*sizeof(uint32_t);
            last_pos = pos;
			fwrite(&pos,sizeof(long),1,fIdx);
			degree = 1;

            //use buffer for nbr
            if(index >= buffersize){
                fwrite(buf,sizeof(uint32_t),buffersize,fDat);
                fwrite(buf_eid,sizeof(uint64_t),buffersize,fEid);
                index = 0;
            }
            buf[index] = v;
            buf_eid[index] = eid;
            index++;

            // //not use buffer for nbr
			// fwrite(&v,sizeof(int),1,fDat);
            // fwrite(&eid,sizeof(uint64_t),1,fEid);
            
		}
        else if(v != previousB){
            //use buffer for nbr
            if(index >= buffersize){
                fwrite(buf,sizeof(uint32_t),buffersize,fDat);
                fwrite(buf_eid,sizeof(uint64_t),buffersize,fEid);
                index = 0;
            }
            buf[index] = v;
            buf_eid[index] = eid;
            index++;
			++degree;

            // fwrite(&v,sizeof(int),1,fDat);  //save neighbor v of vertex u
            // fwrite(&eid,sizeof(uint64_t),1,fEid);
		}
		// if u==previousA & v==previousB, ignore edge(u,v) cause it is same as previous one.
		previousA = u;
		previousB = v;

		// replace es[minIndex] by picking up the first edge from file edges_tmp_minIndex
		if(!fread(&es[minIndex],sizeof(TEdge),1,frl[minIndex])){
			es[minIndex].u = m_maxID;
            write_io += 1;
		}
    }
    if(index)
    {
        fwrite(buf,sizeof(uint32_t),index,fDat);
        fwrite(buf_eid,sizeof(uint64_t),index,fEid);
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
    
    write_io += 5*edgeNum/buffersize;
    write_io += 3*verNum;
	    
    log_debug(readFileClock_.Count("vertex num: %d",vertexNum));
    log_debug(readFileClock_.Count("Max degree: %d",maxDegree));
    log_debug(readFileClock_.Count("Max degree vertex: %d",maxDegVer));
    log_debug(readFileClock_.Count("unordered map size: %d",max_size));

    maxDeg = maxDegree;
    

	FILE* fInfo = fopen(m_info.c_str(),"wb");
	fwrite(&verNum,sizeof(int),1,fInfo);
	fwrite(&maxDegree,sizeof(int),1,fInfo);
    fwrite(&edgeNum,sizeof(uint64_t),1,fInfo);
	fclose(fInfo);
    write_io += 3;

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


template<typename T>
struct comp{
        bool operator()(pair<T,uint32_t> a, pair<T,uint32_t> b){
            return a.second > b.second;
        }
};

//merge data of sort_edge_tmp dir in degree increasing method 
void readFile::mergeByDegSort(int size, uint32_t &vertexNum, uint32_t &maxDegree,
string &name, bool flag, bool isBigG){
    // if flag is true, graph vertex is reorder, otherwise, graph vertex is unorder
    // if isBigG is true, need not to write fEid and fOff datastructure

    // #ifdef semiBinary
    // delete[] m_vertexMap;
    // #endif

    FILE** frl = new FILE*[size];
    TEdge* es = new TEdge[size];
    FILE* fIdx = fopen(m_idx.c_str(),"wb");
	FILE* fDat = fopen(m_dat.c_str(),"wb");
    FILE* fEid;
    FILE* fOff;
    uint32_t buffersize = 10*1024*1024;
    uint32_t *buf = new uint32_t[buffersize]();
    uint64_t *buf_eid = new uint64_t[buffersize]();
    if(!isBigG){
        fEid = fopen(m_eid.c_str(),"wb");
        fOff = fopen(m_offset.c_str(),"wb");
    }
    log_info(readFileClock_.Count("Start merge sort by deg"));

    /*  vertex ID is not need to reorder by degree increasing order */
    for(int i = 0; i < size; i++){
        char filename[100];
        sprintf(filename,"%s%s/edges_tmp_%d",m_base.c_str(),name.c_str(),i);
        frl[i] = fopen(filename,"rb");
        // get first edge of all edges_tmp files
        fread(&es[i],sizeof(TEdge),1,frl[i]);
        
    }
    write_io += size;
    
    int minIndex,previousA = -1,previousB = -1;
    uint32_t maxDegVer = 0;
	int degree = -1, f = 0;
    long pos,last_pos;
    bool sig = false;
    
    uint64_t max_size = 0,index = 0;

    std::unordered_map<uint64_t,uint64_t> eidToPos;
    uint64_t off = 0;
    if(vertexId.size() != 0)
        vertexId.clear();

    while(mergeFinished(es,size)){
        minIndex = minPartition(es,size);
        uint32_t u = es[minIndex].u;
		uint32_t v = es[minIndex].v;
        uint64_t eid = es[minIndex].eid; 
        // if(vertexNum == 108)
        //     printf("u: %u, v: %u, eid: %lu, minIndex: %u, off: %lu\n",u,v,eid,minIndex,off); 
        // if(off == 33028782)
        // {
        //     for(int i = 0 ; i < size; i++){
        //         if(es[i].u != m_maxID)
        //         {
        //             printf("i: %d, %u\n",i,es[i].u);
        //             break;
        //         }
                    
        //     }
        // }
        if(u == previousA && v == previousB) {
            if(!fread(&es[minIndex],sizeof(TEdge),1,frl[minIndex]))
			    es[minIndex].u = m_maxID;
            continue; 
        }
        
        if(!isBigG){
            if(eidToPos.find(eid) != eidToPos.end()){
                eid_eid tmp = {eidToPos[eid],off};
                // random write costs a lot
                fseek(fOff,eid*sizeof(eid_eid),SEEK_SET);
                fwrite(&tmp,sizeof(eid_eid),1,fOff);
                write_io += 1;
                eidToPos.erase(eid);
            }   
            else{
                eidToPos[eid] = off;
            } 
            max_size = std::max(eidToPos.size(),max_size);
        }
        off++;
		
		if(u != previousA){ 			// u != previousA demonstrates that all previousA's neighbors have been writen
            if(!flag)//for unorder method
                vertexId.push_back(u);

			if(degree != -1){ // write the vertex degree in .idx file
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
                if(!isBigG)
                    fwrite(buf_eid,sizeof(uint64_t),buffersize,fEid);
                index = 0;
            }
            buf[index] = v;
            if(!isBigG)
                buf_eid[index] = eid;
            index++;
			// fwrite(&v,sizeof(int),1,fDat);
            // if(!isBigG)
            //     fwrite(&eid,sizeof(uint64_t),1,fEid);
		}
        else if(v != previousB){
            if(!sig && v > previousA ){
                last_pos = pos+degree*sizeof(uint32_t);
                fwrite(&last_pos,sizeof(long),1,fIdx);
                sig = true;
            }
            if(index >= buffersize){
                fwrite(buf,sizeof(uint32_t),buffersize,fDat);
                if(!isBigG)
                    fwrite(buf_eid,sizeof(uint64_t),buffersize,fEid);
                index = 0;
            }
            buf[index] = v;
            if(!isBigG)
                buf_eid[index] = eid;
            index++;
			// fwrite(&v,sizeof(int),1,fDat);  //save neighbor v of vertex u
            // if(!isBigG)
            //     fwrite(&eid,sizeof(uint64_t),1,fEid);
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
        if(!isBigG)
            fwrite(buf_eid,sizeof(uint64_t),index,fEid);
    }

    if(!sig && last_pos == pos){
        last_pos = pos + degree*sizeof(uint32_t);
        fwrite(&last_pos,sizeof(long),1,fIdx);
    }
    fwrite(&degree,sizeof(int),1,fIdx);
    edgeNum = off/2;

    write_io += 2*edgeNum/buffersize;
    if(!isBigG)
        write_io += 3*edgeNum/buffersize;

    if(degree>maxDegree)
    {
        maxDegree = degree;
        maxDegVer = previousA;
    }
    // log_info(readFileClock_.Count("Merge operation done"));
    
    // write the vertex num and max degree
    if(flag)
    {
        vertexNum = previousA+1;
        verNum = vertexNum;
    }
    else
        verNum = vertexId.size();
	write_io += 3*verNum;
    log_debug(readFileClock_.Count("vertex num: %d, edge num: %d",verNum, edgeNum));
    
    // log_debug(readFileClock_.Count("Max degree: %d",maxDegree));
    // log_debug(readFileClock_.Count("Max degree vertex: %d",maxDegVer));
    // log_debug(readFileClock_.Count("Max map size: %u",max_size));

    maxDeg = maxDegree;
    
	FILE* fInfo = fopen(m_info.c_str(),"wb");
	fwrite(&verNum,sizeof(int),1,fInfo);
	fwrite(&maxDegree,sizeof(int),1,fInfo);
    fwrite(&edgeNum,sizeof(uint64_t),1,fInfo);
	fclose(fInfo);

    write_io += 3;

	for (int i = 0; i < size ; ++i){
		fclose(frl[i]);
	}

	fclose(fIdx);
	fclose(fDat);
    if(!isBigG){
        fclose(fEid);
        fclose(fOff);
    }
    

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

void readFile::partitionGraphByVertex(int percentage)
{
    string base = m_base + "maintenance/";
    createDir(base);
    string subDir = base + to_string(percentage) +"/";
    createDir(subDir);
    string dat = subDir+"graph.dat";
    string idx = subDir+"graph.idx";
    string info = subDir+"graph.info";

    FILE* fileIdx = fopen(idx.c_str(),"wb");
	FILE* fileDat = fopen(dat.c_str(),"wb");
    FILE* fileInfo = fopen(info.c_str(),"wb");

    uint32_t buffersize = 10*1024*1024;
    uint32_t *buf = new uint32_t[buffersize]();

    uint32_t total_ver = 0;
    uint64_t total_edge = 0;

    uint32_t sub_ver = verNum * percentage / 10, newNode = 0, last_node = 0;
    double stride = 10.0 / percentage, start = 0;
    stride = round(stride*10)/10.0;
    unordered_set<uint32_t> existVer;
    unordered_map<uint32_t,uint32_t> verMap;
    for(uint32_t i = 0; i < sub_ver; i++){
        start = (10.0 * newNode) / percentage;
        uint32_t node = round(start);
        existVer.insert(node);
        if(last_node == node)
            printf("node: %u, newNode: %u, start: %f, %f\n",node,newNode,start,(10.0 * newNode) / percentage);
        last_node = node;
        verMap[node] = newNode++;  
    }
    total_ver = existVer.size();
    printf("verNum: %u, sub_ver: %u, total_ver: %u, stride: %f, newNode: %u\n",verNum,sub_ver,total_ver,stride,newNode);

    MyReadFile fDat( m_dat );
	fDat.fopen( BUFFERED );
    MyReadFile fIdx( m_idx );
	fIdx.fopen( BUFFERED );

    uint64_t *begPtr = new uint64_t[verNum]();
    uint32_t *degArr = new uint32_t[verNum]();
    uint64_t tmpa = 0,tmpb = 0;
    uint32_t degreeTmp = 0;
    for (uint32_t i = 0; i < verNum; ++i){
		fIdx.fread(&tmpa,sizeof(uint64_t));
        begPtr[i] = tmpa;
        #ifdef DegSort
        fIdx.fread(&tmpb,sizeof(uint64_t));
        #endif
		fIdx.fread(&degreeTmp,sizeof(uint32_t));
        degArr[i] = degreeTmp;
  	}
    uint32_t maxDegVer = 0, maxDegree = 0;
	int degree = -1, f = 0, previousA = -1,previousB = -1, index = 0;
    long pos,last_pos;
    bool sig = false; 

    uint32_t *nbr = new uint32_t[maxDeg]();
    start = 0;
    newNode = 0;
    for(uint32_t i = 0; i < sub_ver; i++){
        start = (10.0 * newNode) / percentage;
        uint32_t node = round(start);
        fDat.fseek(begPtr[node]);
        fDat.fread(nbr,sizeof(uint32_t)*degArr[node]);
        newNode++;
        bool secFor = false;
        for(uint32_t j = 0; j < degArr[node]; j++){
            if(existVer.find(nbr[j]) == existVer.end()) continue;
            secFor = true;
            uint32_t u = verMap[node], v = verMap[nbr[j]];
            // if(i >= 6856439)
            //     printf("u: %u, v: %u, deg: %u, pos: %lu, last_pos: %lu\n",u,v,degree,pos,last_pos);
            total_edge++;
            if(u != previousA){ 			// u != previousA demonstrates that all previousA's neighbors have been writen
                if(degree != -1){ // write the vertex degree in .idx file
                    if(!sig && last_pos == pos){
                        last_pos = pos+degree*sizeof(uint32_t);
                        fwrite(&last_pos,sizeof(long),1,fileIdx);
                    }
                    fwrite(&degree,sizeof(int),1,fileIdx);
                    if(degree>maxDegree)
                    {
                        maxDegree = degree;
                        maxDegVer = previousA;
                    }
                }

                // write start position of neighbor of vertex u from .dat file to .idx file
                pos = ftell(fileDat)+index*sizeof(uint32_t);
                last_pos = pos;
                sig = false;
                fwrite(&pos,sizeof(long),1,fileIdx);
                if(!sig && v > u ){
                    fwrite(&last_pos,sizeof(long),1,fileIdx);
                    sig = true;
                }
                degree = 1;

                if(index >= buffersize){
                    fwrite(buf,sizeof(uint32_t),buffersize,fileDat);
                    index = 0;
                }
                buf[index] = v;
                index++;
            }
            else if(v != previousB){
                if(!sig && v > previousA ){
                    last_pos = pos+degree*sizeof(uint32_t);
                    fwrite(&last_pos,sizeof(long),1,fileIdx);
                    sig = true;
                }
                if(index >= buffersize){
                    fwrite(buf,sizeof(uint32_t),buffersize,fileDat);
                    index = 0;
                }
                buf[index] = v;
                index++;
                ++degree;
            }
            previousA = u;
            previousB = v;
        }
        if(!secFor){
            // printf("xxx\n");
            if(!sig && last_pos == pos){
                last_pos = pos+degree*sizeof(uint32_t);
                fwrite(&last_pos,sizeof(long),1,fileIdx);
            }
            fwrite(&degree,sizeof(int),1,fileIdx);
            pos = ftell(fileDat)+index*sizeof(uint32_t);

            last_pos = pos;
            sig = false;
            fwrite(&pos,sizeof(long),1,fileIdx);
            degree = 0;
        }
    }
    if(index)
    {
        fwrite(buf,sizeof(uint32_t),index,fileDat);
    }

    if(!sig && last_pos == pos){
        last_pos = pos + degree*sizeof(uint32_t);
        fwrite(&last_pos,sizeof(long),1,fileIdx);
    }
    fwrite(&degree,sizeof(int),1,fileIdx);

    total_edge = total_edge/2;
    fwrite(&total_ver,sizeof(int),1,fileInfo);
	fwrite(&maxDegree,sizeof(int),1,fileInfo);
    fwrite(&total_edge,sizeof(uint64_t),1,fileInfo);

    printf("total_ver: %u, total_edge: %lu, maxDeg: %u\n",total_ver,total_edge,maxDegree);

    fDat.fclose();
    fIdx.fclose();
    fclose(fileIdx);
    fclose(fileDat);
    fclose(fileInfo);
    delete[] begPtr;
    delete[] degArr;
    delete[] nbr;
    delete[] buf;
}
void readFile::partitionGraphByEdge(int percentage){

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


