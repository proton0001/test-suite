#include <stdio.h>
#include <time.h>
#include <assert.h>
#include <stdlib.h> 
#include <math.h> 
#include <sys/time.h>
#include <string.h>

#define STR_SIZE (256)
#define MAX_PD	(3.0e6)
/* required precision in degrees	*/
#define PRECISION	0.001
#define SPEC_HEAT_SI 1.75e6
#define K_SI 100
/* capacitance fitting factor	*/
#define FACTOR_CHIP	0.5


#define ITERATION 100
#define NUMLAYERS 4


#ifdef SMALL_DATASET
#define NUMCOLS 256
#else
#define NUMCOLS 512
#endif



/* chip parameters	*/
float t_chip = 0.0005;
float chip_height = 0.016;
float chip_width = 0.016; 
/* ambient temperature, assuming no package at all	*/
float amb_temp = 80.0;

void readinput(float * __restrict__ vect, int grid_rows, int grid_cols, int layers, char *__restrict__ file) {
    int i,j,k;
    FILE *fp;
    char str[STR_SIZE];
    float val;

    if( (fp  = fopen(file, "r" )) ==0 )
      fprintf(stderr, "The file was not opened\n");

    for (i=0; i <= grid_rows-1; i++) 
      for (j=0; j <= grid_cols-1; j++)
        for (k=0; k <= layers-1; k++)
          {
            if (fgets(str, STR_SIZE, fp) == NULL)
            {
                fprintf(stderr, "Error reading file\n");
                exit(1);
            }
            if (feof(fp))
            {
              fprintf(stderr, "not enough lines in file\n");
              exit(1);
            }
            if ((sscanf(str, "%f", &val) != 1))
            {
              fprintf(stderr, "invalid file format\n");
              exit(1);
            }
            vect[i*grid_cols+j+k*grid_rows*grid_cols] = val;
          }
    fclose(fp);	
}


void writeoutput(float *vect, int grid_rows, int grid_cols, int layers) {
    int i,j,k, index=0;
    for (i=0; i < grid_rows; i++) 
      for (j=0; j < grid_cols; j++)
        for (k=0; k < layers; k++) {
          printf("%d\t%.6f\n", index, vect[i*grid_cols+j+k*grid_rows*grid_cols]);
          index++;
        }
}



void computeTempCPU(float *__restrict__ pIn, float* __restrict__ tIn, float *__restrict__ tOut, 
        int nx, int ny, int nz, float Cap, 
        float Rx, float Ry, float Rz, 
        float dt, int numiter) 
{   float ce, cw, cn, cs, ct, cb, cc;
    float stepDivCap = dt / Cap;
    ce = cw =stepDivCap/ Rx;
    cn = cs =stepDivCap/ Ry;
    ct = cb =stepDivCap/ Rz;

    cc = 1.0 - (2.0*ce + 2.0*cn + 3.0*ct);

    int c,w,e,n,s,b,t;
    int x,y,z;
    int i = 0;
    do{
        for(z = 0; z < nz; z++)
            for(y = 0; y < ny; y++)
                for(x = 0; x < nx; x++)
                {
                    c = x + y * nx + z * nx * ny;

                    w = (x == 0) ? c      : c - 1;
                    e = (x == nx - 1) ? c : c + 1;
                    n = (y == 0) ? c      : c - nx;
                    s = (y == ny - 1) ? c : c + nx;
                    b = (z == 0) ? c      : c - nx * ny;
                    t = (z == nz - 1) ? c : c + nx * ny;


                    tOut[c] = tIn[c]*cc + tIn[n]*cn + tIn[s]*cs + tIn[e]*ce + tIn[w]*cw + tIn[t]*ct + tIn[b]*cb + (dt/Cap) * pIn[c] + ct*amb_temp;
                }
        float *temp = tIn;
        tIn = tOut;
        tOut = temp; 
        i++;
    }
    while(i < numiter);

}

int main(int argc, char** argv)
{
    if (argc != 3)
    {
       fprintf(stderr, "Usage: %s <powerFile> <initfile>\n", argv[0]);
       exit(1);
    }

    char *pfile, *tfile, *ofile;

    int iterations = ITERATION;
    tfile = argv[2];
    pfile = argv[1];
    int numCols = NUMCOLS;
    int numRows = NUMCOLS;
    int layers = NUMLAYERS;

    /* calculating parameters*/
    float dx = chip_height/numRows;
    float dy = chip_width/numCols;
    float dz = t_chip/layers;

    float Cap = FACTOR_CHIP * SPEC_HEAT_SI * t_chip * dx * dy;
    float Rx = dy / (2.0 * K_SI * t_chip * dx);
    float Ry = dx / (2.0 * K_SI * t_chip * dy);
    float Rz = dz / (K_SI * dx * dy);

    float max_slope = MAX_PD / (FACTOR_CHIP * t_chip * SPEC_HEAT_SI);
    float dt = PRECISION / max_slope;


    float *powerIn, *tempOut, *tempIn, *tempCopy;// *pCopy;
    //    float *d_powerIn, *d_tempIn, *d_tempOut;
    int size = numCols * numRows * layers;

    powerIn = (float*)calloc(size, sizeof(float));
    tempCopy = (float*)malloc(size * sizeof(float));
    tempIn = (float*)calloc(size,sizeof(float));
    tempOut = (float*)calloc(size, sizeof(float));
    //pCopy = (float*)calloc(size,sizeof(float));
    float* answer = (float*)calloc(size, sizeof(float));

    // outCopy = (float*)calloc(size, sizeof(float));
    readinput(powerIn,numRows, numCols, layers, pfile);
    readinput(tempIn, numRows, numCols, layers, tfile);

    memcpy(tempCopy,tempIn, size * sizeof(float));

    computeTempCPU(powerIn, tempCopy, answer, numCols, numRows, layers, Cap, Rx, Ry, Rz, dt,iterations);
    
    writeoutput(answer,numRows, numCols, layers);
    free(tempIn);
    free(tempOut);
    free(powerIn);
    return 0;
}	


