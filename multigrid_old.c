#include "PDEsolvers.h"
#define SMOOTHNUM 2000

int currSize, cycle;
double **fine, **prevFine, **coarse, **prevCoarse;
// double **matHierarchy[3], **prevHierarchy[3]; // array of matrices. 0: fine, 1: coarser, 2: coarsest
double ***matHierarchy, ***prevHierarchy; // array of matrices. 0: fine, 1: coarser, 2: coarsest

void *MG(void *arg);
void Jacobi(int iters, long myid, double ***mat, double ***prevMat);
void VCycle2(long myid, int smoothIterNum, int iterNum);
void VCycle3(long myid, int smoothIterNum, int iterNum);
void parallel(int numRuns);
void serial(int numRuns);

int main(int argc, char *argv[]) {
    int i, j;
    double start_time, end_time;
    FILE *fp;

    // command line arguments
    gridSize = (argc > 1)? atoi(argv[1]) : GRIDSIZE; // smallest grid
    numIters = (argc > 2)? atoi(argv[2]) : NUMITERS;
    numWorkers = (argc > 3)? atoi(argv[3]) : MAXWORKERS;
    numRuns = (argc > 4)? atoi(argv[4]) : NUMRUNS;

    if (numWorkers > MAXWORKERS) {
        numWorkers = MAXWORKERS;
    }

    int cycle = 3;

    if (cycle == 2) {
        // Initialize matrix memory
        currSize = gridSize*2;
        initMatrixMem(&fine, gridSize*2);
        initMatrixMem(&coarse, gridSize);
        initMatrixVals(&fine, gridSize*2);
        initMatrixVals(&coarse, gridSize);
        
        // initilize "old array" memory
        initMatrixMem(&prevFine, gridSize*2);
        initMatrixMem(&prevCoarse, gridSize);
        initMatrixVals(&prevFine, gridSize*2);
        initMatrixVals(&prevCoarse, gridSize);
    } else {
        // Initialize matrix memory
        currSize = gridSize*2*2;
        matHierarchy = (double ***) malloc(sizeof(double**)*3); 
        prevHierarchy = (double ***) malloc(sizeof(double**)*3);
        
        initMatrixMem(&matHierarchy[0], gridSize*2*2);
        initMatrixMem(&matHierarchy[1], gridSize*2);
        initMatrixMem(&matHierarchy[2], gridSize);

        initMatrixVals(&matHierarchy[0], gridSize*2*2);
        initMatrixVals(&matHierarchy[1], gridSize*2);
        initMatrixVals(&matHierarchy[2], gridSize);
        
        // initilize "old array" memory
        initMatrixMem(&prevHierarchy[0], gridSize*2*2);
        initMatrixMem(&prevHierarchy[1], gridSize*2);
        initMatrixMem(&prevHierarchy[2], gridSize);

        initMatrixVals(&prevHierarchy[0], gridSize*2*2);
        initMatrixVals(&prevHierarchy[1], gridSize*2);
        initMatrixVals(&prevHierarchy[2], gridSize);

        // printf("h: %p\n 2h: %p\n 4h: %p\n", &matHierarchy[0],&matHierarchy[1],&matHierarchy[2]);
    }
    

    // // // declare worker array and attributes
    // pthread_t workers[numWorkers];
    // pthread_attr_t attr;

    // // initialize locks and cond vars
    // pthread_cond_init(&go, NULL);
    // pthread_mutex_init(&barrier, NULL);
    // pthread_mutex_init(&maxErrLock, NULL);

    // // set global thread attributes
    // pthread_attr_init(&attr);
    // pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

    // start_time = read_timer();
    // // create threads
    // for (long l = 0; l < numWorkers; l++) {
    //     pthread_create(&workers[l], &attr, MG, (void *) l);
    // }
    // // join threads whn they are done
    // for (long l = 0; l < numWorkers; l++) {
    //     pthread_join(workers[l], NULL);
    // }
    // end_time = read_timer();

    // // print parallel results
    // printf("max error: %f\n", maxErr);
    // printf("parallel execution time: %f\n", end_time - start_time);

    // parallel
    fp = fopen("./exec_time.data", "w+");
    fprintf(fp, "%d ", gridSize);
    fclose(fp);

    for (int i = 1; i <=4; i++) {
        numWorkers = i;
        parallel(numRuns);
    }
    // prevs implicitly re-initilized
    // serial(numRuns);

    // numWorkers = 1;
    // start_time = read_timer();
    // MG((void *) 0);
    // end_time = read_timer();

    // // print serial results
    // printf("max error: %f\n", maxErr);
    // printf("serial execution time: %f\n", end_time - start_time);

    // write to filedata.out
    if (cycle == 2) {
        fp = fopen("./filedata.out", "w+");
        for (i=0; i<currSize; i++) {
            for (j = 0; j<currSize; j++) {
                fprintf(fp, "%f ", fine[i][j]);
            }
            fputs("\n", fp); // new row
        }
        fclose(fp);

        // terminate matrices, freeing memory
        terminate(&fine);
        terminate(&coarse);
        terminate(&prevFine);
        terminate(&prevCoarse);
    } else {
        fp = fopen("./filedata.out", "w+");
        for (i=0; i<currSize; i++) {
            for (j = 0; j<currSize; j++) {
                fprintf(fp, "%f ", matHierarchy[0][i][j]);
            }
            fputs("\n", fp); // new row
        }
        fclose(fp);

        for (i = 0; i <= 2; i++) {
            terminate(&matHierarchy[i]);
            terminate(&prevHierarchy[i]);
        }

        terminate(matHierarchy);
    }


    const char *command = "python3 plotMat.py";
    system(command);

    command = "python3 plotTimes.py";
    system(command);
}


void parallel(int numRuns) {
    FILE *fp;
    double start_time, end_time;
    double mTime, mErr;
    double times[numRuns], err[numRuns];

    // declare worker array and attributes
    pthread_t workers[numWorkers];
    pthread_attr_t attr;

    // initialize locks and cond vars
    pthread_cond_init(&go, NULL);
    pthread_mutex_init(&barrier, NULL);
    pthread_mutex_init(&maxErrLock, NULL);

    // set global thread attributes
    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);

    for (int i = 0; i < numRuns; i++) {
        start_time = read_timer();
        // create threads
        for (long l = 0; l < numWorkers; l++) {
            pthread_create(&workers[l], &attr, MG, (void *) l);
        }
        // join threads whn they are done
        for (long l = 0; l < numWorkers; l++) {
            pthread_join(workers[l], NULL);
        }
        end_time = read_timer();
        // store current
        times[i] = end_time - start_time;
        err[i] = maxErr;

        // re-initialize prev
        if (cycle != 2) {
            // initMatrixVals(&prevHierarchy[0], gridSize*2*2);
            // initMatrixVals(&prevHierarchy[1], gridSize*2);
            // initMatrixVals(&prevHierarchy[2], gridSize);
        } else {
            initMatrixVals(&prevFine, gridSize*2);
            initMatrixVals(&prevCoarse, gridSize);
        }   
    }
    // sort
    qsort(times, numRuns, sizeof(double), cmp);
    qsort(err, numRuns, sizeof(double), cmp);

    if (numRuns % 2 == 1) {
        mTime = times[numRuns/2];
        mErr = err[numRuns/2];
    } else {
        mTime = (times[numRuns/2-1] + times[numRuns/2]) / 2;
        mErr = (err[numRuns/2-1] + err[numRuns/2]) / 2;
    }

    fp = fopen("/home/conc/Documents/Kurser 2023/VT/ID1217/project/exec_time.data", "a+");
    fprintf(fp, "%f ", mTime);
    fclose(fp);

    // print parallel results
    printf("median max error: %f\n", mErr);
    printf("# of workers: %d, execution median time (%d runs): %f\n\n", numWorkers, numRuns, mTime);
}


void serial(int numRuns) {
    double start_time, end_time;
    double mTime, mErr;
    double times[numRuns], err[numRuns];
    numWorkers = 1; // i.e. serial

    for (int i = 0; i < numRuns; i++) {
        start_time = read_timer();
        MG((void *) 0);
        end_time = read_timer();
        // store time and error for current simulation
        times[i] = end_time - start_time;
        err[i] = maxErr;
        // re-initialize prev
        // initMatrixVals(&prevHierarchy[0], gridSize*2*2);
        // initMatrixVals(&prevHierarchy[1], gridSize*2);
        // initMatrixVals(&prevHierarchy[2], gridSize);

        initMatrixVals(&prevFine, gridSize*2);
        initMatrixVals(&prevCoarse, gridSize);

    }

    // sort time and error arrays
    qsort(times, numRuns, sizeof(double), cmp);
    qsort(err, numRuns, sizeof(double), cmp);

    // find median
    if (numRuns % 2 == 1) {
        mTime = times[numRuns/2];
        mErr = err[numRuns/2];
    } else {
        mTime = (times[numRuns/2-1] + times[numRuns/2]) / 2;
        mErr = (err[numRuns/2-1] + err[numRuns/2]) / 2;
    }

    // print serial results
    printf("median max error: %f\n", mErr);
    printf("serial execution median time (%d runs): %f\n", numRuns, mTime);
}


// void *MG(void *arg) {
//     long myid = (long) arg;
//     int startRow, endRow;

//     int N = 10;
//     for (int i = 1; i <= N; i++) {
//         VCycle2(myid, SMOOTHNUM, numIters/N);
//     }

//     getWorkerRows(myid, gridSize*2, &startRow, &endRow);
//     // max error
//     for (int i = startRow; i <= endRow; i++) { // rows
//         for (int j = 0; j < gridSize; j++) { // columns
//             if (fabs(fine[i][j] - 1) > maxErr) {

//                 if (numWorkers != 1) {
//                     pthread_mutex_lock(&maxErrLock);
//                     if (fabs(fine[i][j] - 1) > maxErr) {
//                         maxErr = fabs(fine[i][j] - 1);
//                     }
//                     pthread_mutex_unlock(&maxErrLock);
//                 } 
//                 else {
//                     maxErr = fabs(fine[i][j] - 1);
//                 }
//             }
//         }
//     }
// }


void *MG(void *arg) {
    long myid = (long) arg;
    int startRow, endRow;
    int fineSize = currSize;

    if (cycle == 2) {
        VCycle2(myid, SMOOTHNUM, numIters);
    } else {
        VCycle3(myid, SMOOTHNUM, numIters);
    }
    
    getWorkerRows(myid, fineSize, &startRow, &endRow);
    // max error
    if (cycle != 2) {
        for (int i = startRow; i <= endRow; i++) { // rows
            for (int j = 0; j < fineSize; j++) { // columns
                if (fabs(matHierarchy[0][i][j] - 1) > maxErr) {
                    pthread_mutex_lock(&maxErrLock);
                    if (fabs(matHierarchy[0][i][j] - 1) > maxErr) {
                        maxErr = fabs(matHierarchy[0][i][j] - 1);
                    }
                    pthread_mutex_unlock(&maxErrLock);
                }
            }
        }
    } 
    else {
        for (int i = startRow; i <= endRow; i++) { // rows
            for (int j = 0; j < fineSize; j++) { // columns
                if (fabs(fine[i][j] - 1) > maxErr) {
                    pthread_mutex_lock(&maxErrLock);
                    if (fabs(fine[i][j] - 1) > maxErr) {
                        maxErr = fabs(fine[i][j] - 1);
                    }
                    pthread_mutex_unlock(&maxErrLock);
                }
            }
        }
    }
    
}


void VCycle2(long myid, int smoothIterNum, int iterNum) {
    int startRow, endRow;
    getWorkerRows(myid, currSize, &startRow, &endRow);

    // V-cycle
    // step 1: smooth (iterate using jacobi)
    Jacobi(smoothIterNum, myid, &fine, &prevFine);

    // step 2: restrict fine to coarse
    if (myid == 0) {
        currSize /= 2;
    }
    Barrier(myid);
    restriction(myid, currSize, &fine, &coarse);
    Barrier(myid);

    // step 3: solve (iterate using jacobi)
    Jacobi(iterNum, myid, &coarse, &prevCoarse);

    // step 4: interpolate coarse to fine
    if (myid == 0) {
        currSize *= 2;
    }
    Barrier(myid);
    interpolation(myid, currSize, &fine, &coarse);
    Barrier(myid);

    // step 5: smooth (iterate using jacobi)
    Jacobi(smoothIterNum, myid, &fine, &prevFine);
}


void VCycle3(long myid, int smoothIterNum, int iterNum) {
    int startRow, endRow;
    getWorkerRows(myid, currSize, &startRow, &endRow);

    // V-cycle
    // step 1: smooth (iterate using jacobi) (h)
    Jacobi(smoothIterNum, myid, &matHierarchy[0], &prevHierarchy[0]);

    // step 2: restrict fine to coarse (2h)
    if (myid == 0) {
        currSize /= 2;
    }
    Barrier(myid);
    restriction(myid, currSize, &matHierarchy[0], &matHierarchy[1]);
    Barrier(myid);

    // step 3: smooth again (iterate using jacobi) (2h)
    Jacobi(smoothIterNum, myid, &matHierarchy[1], &prevHierarchy[1]);

    // step 4: restrict again (4h)
    if (myid == 0) {
        currSize /= 2;
    }
    Barrier(myid);
    restriction(myid, currSize, &matHierarchy[1], &matHierarchy[2]);
    Barrier(myid);

    // step 5: solve (iterate using Jacobi) (4h)
    Jacobi(numIters, myid, &matHierarchy[2], &prevHierarchy[2]);

    // step 6: interpolate (2h)
    if (myid == 0) {
        currSize *= 2;
    }
    Barrier(myid);
    interpolation(myid, currSize, &matHierarchy[1], &matHierarchy[2]);
    Barrier(myid);

    // step 7: smooth (iterate using jacobi) (2h)
    Jacobi(smoothIterNum, myid, &matHierarchy[1], &prevHierarchy[1]);

    // step 8: interpolate again (h)
    if (myid == 0) {
        currSize *= 2;
    }
    Barrier(myid);
    interpolation(myid, currSize, &matHierarchy[0], &matHierarchy[1]);
    Barrier(myid);

    // step 9: smooth (iterate using jacobi) (h)
    Jacobi(smoothIterNum, myid, &matHierarchy[0], &prevHierarchy[0]);
}


void Jacobi(int iters, long myid, double ***mat, double ***prevMat) {
    double **tmp;
    int startRow, endRow; // worker rows

    getWorkerRows(myid, currSize, &startRow, &endRow);

    if (myid == 0) {
        startRow++;
    }
    if (myid == numWorkers-1) {
        endRow--;
    }
    
    for (int k = 0; k < iters; k++) { // outer loop
        // inner loops
        for (int i = startRow; i <= endRow; i++) { // rows
            for (int j = 1; j < currSize-1; j++) { // columns
                (*mat)[i][j] = ((*prevMat)[i-1][j] + (*prevMat)[i+1][j] + (*prevMat)[i][j-1] + (*prevMat)[i][j+1]) / 4; // FD
            }
        }

        // wait for all threads to be done before preceeding to the next time step
        Barrier(myid);
        // copy new matrix to old matrix then wait before continuing next time step
        if (myid == 0) {
            tmp = *prevMat;
            *prevMat = *mat;
            *mat = tmp;
        }
        Barrier(myid);
    }
}