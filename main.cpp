#include <cassert>
#include <iomanip>
#include "include/Time.h"
#include "include/graph.h"
#include "include/fileLoader.h"
#include "include/util.h"



int main(int argc, char **argv){
    if(argc != 9 || strcmp(argv[1],"-f") != 0){
        printf("Usage: -f [data_file_path] -m [method] -s [name] -d [dir]\n");
        exit(1);
    }
    const char* filepath = argv[2];
    uint32_t method = atoi(argv[4]);
    uint32_t persentate = atoi(argv[6]);
    char* name = argv[8];
    string str = name;


    Clock allClock("All");
    log_info(allClock.Start());
    string midResultBase = "/home/jjq/research/ktruss/midResult/"+str;
    
    readFile file(midResultBase+"/");
    // file.createDir(midResultBase);
    file.loadFile(filepath);
    log_debug(allClock.Count("file len: %ld",file.getLen()));
    file.LoadGraph();
    log_info(allClock.Count("Load graph done"));


    Graph g(file.verNum,file.edgeNum);
    g.Initial(file);


    g.kCoreTrussDecomPlus(file);
    log_info(allClock.Count("kCoreTrussDecomPlus done"));
    g.dynamicMaxTruss(file);  
    log_info(allClock.Count("dynamic MaxTruss done"));

    return 0;

}