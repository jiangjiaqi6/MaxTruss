# Top-Down: Truss Decomposition in Massive Networks

## run command
```
./main -f ../../ktruss-data/correct-youtube.csv -m 1
```
- -f parameter represents graph data location
- -m parameter represents method which to get max truss subgraph
    - m = 0, the implementation is similar with peel algorithm,iteratively delete edges whose sup is the most small
    - m = 1, the implementation of the closed-source TopDown algorithm in paper "Truss Decomposition in Massive Networks"
