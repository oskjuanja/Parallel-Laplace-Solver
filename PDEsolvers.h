#ifndef PDESOLVERS_H   /* Include guard */
#define PDESOLVERS_H
    #ifndef _REENTRANT 
    #define _REENTRANT 
    #endif 
    #include <pthread.h>
    #include <stdlib.h>
    #include <stdio.h>
    #include <string.h>
    #include <stdbool.h>
    #include <time.h>
    #include <sys/time.h>
    #include <unistd.h>
    #include <math.h>
    #define GRIDSIZE 100  // matrix size
    #define NUMITERS 1000   // number of iterations
    #define MAXWORKERS 8   // maximum number of workers
    #define NUMRUNS 5       // number of runs per computation

    extern int gridSize;
    extern int numIters;
    extern int numWorkers;     
    extern int numRuns;
    extern double maxErr;          // shared variable
    extern int numArrived;         // shared variable for the barrier

    extern pthread_mutex_t barrier;    // barrier lock
    extern pthread_mutex_t maxErrLock; // lock for shared variable maxErr
    extern pthread_cond_t go;          // cond var for barrier

    /* memory functions */
    void initMatrixMem(double ***mat, int size);
    void initMatrixVals(double ***mat, int size);
    void terminate(double ***mat);

    /* finite difference functions*/
    void restriction(long myid, int coarseSize, double ***fine, double ***coarse);
    void interpolation(long myid, int fineSize, double ***fine, double ***coarse);
    void interpUpdater(double ***fine, double ***coarse, int i, int j);

    /* misc */
    double read_timer();
    void Barrier(long id);
    void getWorkerRows(int myid, int size, int *start, int *end);
    int cmp(const void *a, const void *b);

#endif