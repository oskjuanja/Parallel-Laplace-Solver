#include "PDEsolvers.h"

double **matrix, **prev;    // shared matrices

void parallel(int numRuns);
void serial(int numRuns);
void *Jacobi(void *arg);

int main(int argc, char *argv[]) {
    int i, j;
    FILE *fp;

    // command line arguments
    gridSize = (argc > 1)? atoi(argv[1]) : GRIDSIZE;
    numIters = (argc > 2)? atoi(argv[2]) : NUMITERS;
    numWorkers = (argc > 3)? atoi(argv[3]) : MAXWORKERS;
    numRuns = (argc > 4)? atoi(argv[4]) : NUMRUNS;

    if (numWorkers > MAXWORKERS) {
        numWorkers = MAXWORKERS;
    }

    // Initialize matrix memory
    initMatrixMem(&matrix, gridSize);
    initMatrixVals(&matrix, gridSize);
    // initilize "old array" memory
    initMatrixMem(&prev, gridSize);
    initMatrixVals(&prev, gridSize);

    // parallel
    fp = fopen("/home/conc/Documents/Kurser 2023/VT/ID1217/project/exec_time.data", "w+");
    fprintf(fp, "%d ", gridSize);
    fclose(fp);
    
    for (i=1; i<=4; i++) {
        numWorkers = i;
        parallel(numRuns);
    }

    // write to filedata.out
    fp = fopen("/home/conc/Documents/Kurser 2023/VT/ID1217/project/filedata.out", "w+");
    for (i=0; i<gridSize; i++) {
        for (j = 0; j<gridSize; j++) {
            fprintf(fp, "%f ", matrix[i][j]);
        }
        fputs("\n", fp); // new row
    }
    fclose(fp);

    // terminate matrix, freeing memory
    terminate(&matrix);
    terminate(&prev);

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
            pthread_create(&workers[l], &attr, Jacobi, (void *) l);
        }
        // join threads whn they are done
        for (long l = 0; l < numWorkers; l++) {
            pthread_join(workers[l], NULL);
        }
        end_time = read_timer();

        times[i] = end_time - start_time;
        err[i] = maxErr;

        // re-initialize prev
        initMatrixVals(&prev, gridSize);
    }

    // sort time and error arrays
    qsort(times, numRuns, sizeof(double), cmp);
    qsort(err, numRuns, sizeof(double), cmp);

    if (numRuns % 2 == 1) {
        mTime = times[numRuns/2];
        mErr = err[numRuns/2];
    } else {
        mTime = (times[numRuns/2-1] + times[numRuns/2]) / 2;
        mErr = (err[numRuns/2-1] + err[numRuns/2]) / 2;
    }

    // print parallel results
    printf("median max error: %f\n", mErr);
    printf("# of workers: %d, execution median time (%d runs): %f\n\n", numWorkers, numRuns, mTime);

    fp = fopen("/home/conc/Documents/Kurser 2023/VT/ID1217/project/exec_time.data", "a+");
    fprintf(fp, "%f ", mTime);
    fclose(fp);
}


void serial(int numRuns) {
    double start_time, end_time;
    double mTime, mErr;
    double times[numRuns], err[numRuns];

    numWorkers = 1;

    for (int i = 0; i < numRuns; i++) {
        start_time = read_timer();
        Jacobi((void *) 0);
        end_time = read_timer();

        // re-initialize prev
        initMatrixVals(&prev, gridSize);
    }

    // sort time and error arrays
    qsort(times, numRuns, sizeof(double), cmp);
    qsort(err, numRuns, sizeof(double), cmp);

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


void *Jacobi(void *arg) {
    double **tmp;
    long myid = (long) arg;
    
    int startRow, endRow; // worker rows

    getWorkerRows(myid, gridSize, &startRow, &endRow);

    if (myid == 0) {
        startRow++;
    }
    if (myid == numWorkers-1) {
        endRow--;
    }
    
    for (int k = 0; k < numIters; k++) { // outer loop
        // inner loops
        for (int i = startRow; i <= endRow; i++) { // rows
            for (int j = 1; j < gridSize-1; j++) { // columns
                matrix[i][j] = (prev[i-1][j] + prev[i+1][j] + prev[i][j-1] + prev[i][j+1]) / 4; // FD
            }
        }

        // wait for all threads to be done before preceeding to the next time step
        Barrier(myid);
        // copy new matrix to old matrix then wait before continuing next time step
        if (myid == 0) {
            // memcpy(prev, matrix, sizeof(double)*gridSize*gridSize);
            tmp = prev;
            prev = matrix;
            matrix = tmp;
        }
        Barrier(myid);
    }

    // max error
    for (int i = startRow; i <= endRow; i++) { // rows
        for (int j = 0; j < gridSize; j++) { // columns
            if (fabs(matrix[i][j] - 1) > maxErr) {

                if (numWorkers != 1) {
                    pthread_mutex_lock(&maxErrLock);
                    if (fabs(matrix[i][j] - 1) > maxErr) {
                        maxErr = fabs(matrix[i][j] - 1);
                    }
                    pthread_mutex_unlock(&maxErrLock);
                    }
                else {
                    maxErr = fabs(matrix[i][j] - 1);
                }
            }
        }
    }
}
