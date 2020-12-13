/* Simple CUDA library for APSP problem
 *
 * Author: Matuesz Bojanowski
 *  Email: bojanowski.mateusz@gmail.com
 */

#ifndef _CUDA_APSP_
#define _CUDA_APSP_

#include "common.h"

// CONSTS for CUDA FW
#define BLOCK_SIZE 16
#define MAX_DISTANCE UBND_FLT

struct graphAPSPTopology {
    unsigned int nvertex; // number of vertex in graph
    std::unique_ptr<int[]> pred; // predecessors matrix
    std::unique_ptr<float[]> graph; // graph matrix

    /* Constructor for init fields */
    graphAPSPTopology(int nvertex) : nvertex(nvertex) {
        int size = nvertex * nvertex;
        pred = std::unique_ptr<int[]>(new int[size]);
        graph = std::unique_ptr<float[]>(new float[size]);
    }
 
};


/**
 * Naive implementation of Floyd Warshall algorithm in CUDA
 *
 * @param data: unique ptr to graph data with allocated fields on host
 */
void cudaNaiveFW(const std::unique_ptr<graphAPSPTopology>& dataHost);

/**
 * Blocked implementation of Floyd Warshall algorithm in CUDA
 *
 * @param data: unique ptr to graph data with allocated fields on host
 */
void cudaBlockedFW(const std::unique_ptr<graphAPSPTopology>& dataHost);


#endif /* _APSP_ */