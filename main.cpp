#include <cassert>
#include <iomanip>
#include "include/Time.h"
#include "include/graph.h"
#include "include/fileLoader.h"
#include "include/util.h"



int main(int argc, char **argv){
    if(argc != 7 || strcmp(argv[1],"-f") != 0){
        printf("Usage: -f [data_file_path] -m [method]  -d [dir]\n");
        exit(1);
    }
    const char* filepath = argv[2];
    uint32_t method = atoi(argv[4]);
    // uint32_t persentate = atoi(argv[6]);
    char* name = argv[6];
    string str = name;


    Clock allClock("All");
    log_info(allClock.Start());
    string midResultBase = "/home/jjq/research/ktruss/midResult/"+str;
    
    readFile file(midResultBase+"/");
    file.createDir(midResultBase);
    file.loadFile(filepath);
    file.LoadGraph();
    log_info(allClock.Count("Load graph done"));

    Graph g(file.verNum,file.edgeNum);
    g.Initial(file);

    switch (method)
    {
    case 0:   //SemiBinary
        #ifdef DegSort
        g.CountTriangleSSDByDegOrder(file,true);
        #else
        g.CountTriangleSSD(file,true);
        #endif
        log_info(allClock.Count("Triangle count done"));

        g.binaryAndIncremental(file,1);
        log_info(allClock.Count("binaryAndIncremental done"));
        break;
    
    case 1:  // If -DLazyUpdate is not defined, it is the SemiGreedyCore, otherwise it is the SemiLazyUpdate.
        g.CoreTrussDecomPlus(file);
        log_info(allClock.Count("CoreTrussDecomPlus done"));
        #ifdef Maintenance
        g.dynamicMaxTrussMaintenance(file);  
        // g.dynamicMaxTrussInsertion_YLJ(file);
        // g.dynamicMaxTrussDeletion_YLJ(file);
        // log_info(allClock.Count("dynamic MaxTruss done"));
        #endif
        break;
    default:
        break;
    }

    return 0;

}