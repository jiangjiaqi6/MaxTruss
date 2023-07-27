# MaxTruss : Efficient algorithm for max Truss on Large-Scale Graph
## Notice for Compilation Options
### -DDegSort
```
All vertices on the graph are renumbered in ascending degree order, so the process of calculating triangles can be accelerated.
```

### -DLazyUpdate

```
In this function "secondCoreTrussDecomPlus", we need to use this -DLazyUpdate.
```

### graph_origin.cpp

```
The difference between graph.cpp and graph_origin.cpp is that the former implements optimal intersect function with sup_u and sup_v array when delete edges.The latter method is to find the intersection of two vertexs. Each time I/O queries whether the edge is deleted, which reduces the efficiency.
```


## Executable file
- main_1: binary search optimization
- main_3: binary search optimization and maxCore reduction
- main_6: binary search, maxCore reduction and lazyUpdate strategy
- maintenance_pre: perform complete lazyUpdate method, in order to make preparation for dynamic updates on graph 
- maintenance_3: dynamically delete and insert (XH method) edges
- maintenance_6: dynamically delete and insert edges
