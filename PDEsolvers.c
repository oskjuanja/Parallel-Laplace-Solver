#include "PDEsolvers.h"

int gridSize;
int numIters;
int numWorkers;
int numRuns;     
double maxErr = 0;          // shared variable
int numArrived = 0;         // shared variable for the barrier


pthread_mutex_t barrier;    // barrier lock
pthread_mutex_t maxErrLock; // lock for shared variable maxErr
pthread_cond_t go;          // cond var for barrier

/* 
restrict from fine to coars
Params: current size, pointer to fine matrix, pointer to coarse matrix
Returns: void
 */
void restriction(long myid, int coarseSize, double ***fine, double ***coarse) {
    int i, j; // fine grid coordinates
    int startRow, endRow;

    getWorkerRows(myid, coarseSize, &startRow, &endRow);

    for (int I = startRow; I < endRow; I+=2) {
        for (int J = 0; J < coarseSize; J++) {
            i = 2*I + 1;
            j = 2*J + 1;
            (*coarse)[I][J] = ((*fine)[i][j]) / 2 + 
                              ((*fine)[i-1][j] + (*fine)[i+1][j] + (*fine)[i][j-1] + (*fine)[i][j+1]) / 8;
        }
    }
}

/* 
interpolate from coarse to fine
Params: current size, pointer to fine matrix, pointer to coarse matrix
Returns: void
 */
void interpolation(long myid, int fineSize, double ***fine, double ***coarse) {
    int i,j; // fine grid coordinates
    int startRow, endRow;

    getWorkerRows(myid, fineSize, &startRow, &endRow);

    // printf("myid: %ld, start: %d, end: %d\n", myid, startRow, endRow);

    for (int i = startRow; i < endRow; i++) {
        for (int i = 0; j < fineSize; j++) {
            interpUpdater(fine, coarse, i, j);
        }
    }
}


void interpUpdater(double ***fine, double ***coarse, int i, int j) {
    int I, J; // coarse grid coordinates
    double res = 0.0;

    if (i % 2 == 0) { // even rows
        if (j % 2 == 0) { // even columns
            // coarse grid contributions will be weighted by 1/4
            if (i >= 2 && j >= 2) {
                I = ((i-1) - 1)/2; 
                J = ((j-1) - 1)/2;
                res += (*coarse)[I][J] / 4; // upper left

                J = ((j+1) - 1)/2;
                res += (*coarse)[I][J] / 4; // upper right

                I = ((i+1) - 1)/2;
                J = ((j+1) - 1)/2;
                res += (*coarse)[I][J] / 4; // lower right

                J = ((j-1) - 1)/2;
                res += (*coarse)[I][J] / 4; // lower left

                res /= 4; // average of the points
            }
            else if (i >= 2) {
                I = ((i-1) - 1)/2;
                res += (*coarse)[I][J] / 4; // upper right

                I = ((i+1) - 1)/2;
                J = ((j+1) - 1)/2;
                res += (*coarse)[I][J] / 4; // lower right

                res /= 2; // average
            } 
            else if (j >= 2) {
                I = ((i+1) - 1)/2;
                J = ((j+1) - 1)/2;
                res += (*coarse)[I][J] / 4; // lower right

                J = ((j-1) - 1)/2;
                res += (*coarse)[I][J] / 4; // lower left
            }
            else { // (0,0)
                I = ((i+1) - 1)/2;
                J = ((j+1) - 1)/2;
                res = (*coarse)[I][J] / 4; // lower right
            }
        }
        else { // odd columns
            // coarse grid contributions will be weighted by 1/2
            if (i >= 2) {
                I = ((i+1) - 1)/2;
                J = (j - 1)/2;
                res += (*coarse)[I][J] / 2; // lower

                I = ((i-1) - 1)/2;
                res += (*coarse)[I][J] / 2; // upper

                res /= 2; // average
            }
            else {
                I = ((i+1) - 1)/2;
                J = (j - 1)/2;
                res = (*coarse)[I][J] / 2; // lower
            }
        }
    }
    else { // odd rows
        if (j % 2 == 0) { // even columns
            // coarse grid contributions will be weighted by 1/2
            if (j >= 2) {
                I = (i - 1)/2;
                J = ((j-1) - 1)/2;
                res += (*coarse)[I][J] / 2; // left point

                J = ((j+1) - 1)/2;
                res += (*coarse)[I][J] / 2; // right point

                res /= 2; // average
            } 
            else {
                I = (i - 1)/2;
                J = ((j+1) - 1)/2;
                res = (*coarse)[I][J] / 2; // right point
            }
        }
        else { // odd columns
            // coarse grid contributions will be weighted by 1
            I = (i-1)/2;
            J = (j-1)/2;
            res = (*coarse)[I][J];
        }
    }

    (*fine)[i][j] = res;
}


void getWorkerRows(int myid, int size, int *start, int *end) {
    int startRow, endRow; // worker rows
    int rpw = (size / numWorkers);
    int slack = size % numWorkers;
    int offset;

    // partition worker rows
    startRow = myid * rpw;
    endRow = (myid + 1)*rpw - 1;

    if (myid != 0) {
        offset = (myid < slack) ? myid : slack;
        startRow += offset;
        endRow += offset;
    }
    if (myid + 1 <= slack) {
        endRow++;
    }

    *start = startRow;
    *end = endRow;
}


/* a reusable counter barrier */
// nicked this from first assignment
void Barrier(long myid) {
    // if (numWorkers == 1) {
    //     return;
    // }
    pthread_mutex_lock(&barrier);
    numArrived++;
    if (numArrived == numWorkers) 
    {
        numArrived = 0;
        pthread_cond_broadcast(&go);
    } else {
        pthread_cond_wait(&go, &barrier);
    }
    pthread_mutex_unlock(&barrier);
}


double read_timer() {
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}


void initMatrixMem(double ***mat, int size) {
    *mat = (double **) malloc(sizeof(double)*size*size); // allocate rows
    for (int i = 0; i < size; i++) {
        (*mat)[i] = (double *) malloc(sizeof(double)*size); // allocate columns in rows
        memset((*mat)[i], 0, sizeof(double)*size); // good practice
    }
    // printf("size: %d, size of mat: %ld\n", size, sizeof(*mat));
}


void initMatrixVals(double ***mat, int size) {
    // Initialize the boundary points of the grids to 1.0 and the interior points to 0.0
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            if (i == 0 || j == 0 || i == size-1 || j == size-1) {
            // if (i == 0 || i == size-1) {
                (*mat)[i][j] = 1.0;
            } else {
                (*mat)[i][j] = 0.0;
            }
        }
    }
}


void terminate(double ***mat) {
    // free memory for each word
    int size = sizeof(*mat)/sizeof(double);
    for (int i = 0; i < size; i++) {
        free((*mat)[i]); // free each row
    }
    // free rows
    free(*mat);
}


int cmp(const void *a, const void *b) {
    if (*(float*) a < *(float*) b) {
        return -1;
    } else if (*(float*) a > *(float*) b) {
        return 1;
    } else {
        return 0;
    }
}