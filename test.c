#include <stdio.h>
#include <stdlib.h>
#define FINESIZE 10
#define COARSESIZE 5

double fine[FINESIZE][FINESIZE];
double coarse[COARSESIZE][COARSESIZE];

void restriction();
void interpolation();
void init();

int main() {
    int i, j;

    FILE *fp;
    
    init();
    restriction();

    // write to filedata.out
    fp = fopen("/home/conc/Documents/Kurser 2023/VT/ID1217/project/test.txt", "w+");
    for (i=0; i<FINESIZE; i++) {
        for (j = 0; j<FINESIZE; j++) {
            fprintf(fp, "%f ", fine[i][j]);
        }
        fputs("\n", fp); // new row
    }
    fputs("\n", fp);
    for (int i = 0; i < COARSESIZE; i++) {
        for (int j = 0; j < COARSESIZE; j++) {
            fprintf(fp, "%f ", coarse[i][j]);
        }
        fputs("\n", fp); // new row
    }
    fclose(fp);

}


void restriction() {
    for (int i = 1; i < COARSESIZE-1; i++) {
        for (int j = 1; j < COARSESIZE-1; j++) {
            coarse[i][j] = (fine[i*2][j*2]) / 2 + (fine[i*2-1][j*2] + fine[i*2 + 1][j*2] + fine[i*2][j*2-1] + fine[i*2][j*2+1]) / 8;
        }
    }
}

void interpolation() {
    for (int i = 1; i < COARSESIZE-1; i++) {
        for (int j = 1; j < COARSESIZE-1; j++) {
            // 4 points need to be filled in 
        }
    }
}

void init() {
    int i, j;

    for (i = 0; i < FINESIZE; i++) {
        for (j = 0; j < FINESIZE; j++) {
            if (i == 0 || j == 0 || i == FINESIZE-1 || j == FINESIZE-1) {
                fine[i][j] = 1.0;
            } else {
                fine[i][j] = 0.0;
            }
        }
    }

    for (i = 0; i < COARSESIZE; i++) {
        for (j = 0; j < COARSESIZE; j++) {
            if (i == 0 || j == 0 || i == COARSESIZE-1 || j == COARSESIZE-1) {
                coarse[i][j] = 1.0;
            } else {
                coarse[i][j] = 0.0;
            }
        }
    }
}