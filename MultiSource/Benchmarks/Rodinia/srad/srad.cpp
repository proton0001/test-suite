// srad.cpp : Defines the entry point for the console application.

/*
 * Modified by Pankaj Kukreja
 * Indian Institute of Technology, Hyderabad, India
 */

#define OUTPUT
#define OPEN
#define ITERATION
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define ROWS 2048
#define COLS 2048
#define Y1 0
#define Y2 127
#define X1 0
#define X2 127
#define LAMDBA 0.5

#ifdef SMALL_DATASET
#define ITER 2
#else
#define ITER 5
#endif


void random_matrix(float * I, int rows, int cols);

int main(int argc, char *argv[]) {
  int rows, cols, size_I, size_R, niter, iter, k;
  float * __restrict__ I, * __restrict__ J, q0sqr, sum, sum2, tmp, meanROI, varROI;
  float Jc, G2, L, num, den, qsqr;
  int * __restrict__ iN, * __restrict__ iS, * __restrict__ jE, * __restrict__ jW;
  float * __restrict__ dN, * __restrict__ dS, * __restrict__ dW, * __restrict__ dE;
  int r1, r2, c1, c2;
  float cN, cS, cW, cE;
  float * __restrict__ c, D;
  float lambda;
  int i, j;

  rows = ROWS; // number of rows in the domain
  cols = COLS; // number of cols in the domain
  if ((rows % 16 != 0) || (cols % 16 != 0)) {
    fprintf(stderr, "rows and cols must be multiples of 16\n");
    exit(1);
  }
  r1 = Y1;         // y1 position of the speckle
  r2 = Y2;         // y2 position of the speckle
  c1 = X1;         // x1 position of the speckle
  c2 = X2;         // x2 position of the speckle
  lambda = LAMDBA; // Lambda value
  niter = ITER;    // number of iterations

  size_I = cols * rows;
  size_R = (r2 - r1 + 1) * (c2 - c1 + 1);

  I = (float *)malloc(size_I * sizeof(float));
  J = (float *)malloc(size_I * sizeof(float));
  c = (float *)malloc(sizeof(float) * size_I);

  iN = (int *)malloc(sizeof(unsigned int *) * rows);
  iS = (int *)malloc(sizeof(unsigned int *) * rows);
  jW = (int *)malloc(sizeof(unsigned int *) * cols);
  jE = (int *)malloc(sizeof(unsigned int *) * cols);

  dN = (float *)malloc(sizeof(float) * size_I);
  dS = (float *)malloc(sizeof(float) * size_I);
  dW = (float *)malloc(sizeof(float) * size_I);
  dE = (float *)malloc(sizeof(float) * size_I);

  for (int i = 0; i < rows; i++) {
    iN[i] = i - 1;
    iS[i] = i + 1;
  }
  for (int j = 0; j < cols; j++) {
    jW[j] = j - 1;
    jE[j] = j + 1;
  }
  iN[0] = 0;
  iS[rows - 1] = rows - 1;
  jW[0] = 0;
  jE[cols - 1] = cols - 1;


  random_matrix(I, rows, cols);

  for (k = 0; k < size_I; k++) {
    J[k] = (float)exp(I[k]);
  }


#ifdef ITERATION
  for (iter = 0; iter < niter; iter++) {
#endif
    sum = 0;
    sum2 = 0;
    for (i = r1; i <= r2; i++) {
      for (j = c1; j <= c2; j++) {
        tmp = J[i * cols + j];
        sum += tmp;
        sum2 += tmp * tmp;
      }
    }
    meanROI = sum / size_R;
    varROI = (sum2 / size_R) - meanROI * meanROI;
    q0sqr = varROI / (meanROI * meanROI);

    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {

        k = i * cols + j;
        Jc = J[k];

        // directional derivates
        dN[k] = J[iN[i] * cols + j] - Jc;
        dS[k] = J[iS[i] * cols + j] - Jc;
        dW[k] = J[i * cols + jW[j]] - Jc;
        dE[k] = J[i * cols + jE[j]] - Jc;

        G2 = (dN[k] * dN[k] + dS[k] * dS[k] + dW[k] * dW[k] + dE[k] * dE[k]) /
             (Jc * Jc);

        L = (dN[k] + dS[k] + dW[k] + dE[k]) / Jc;

        num = (0.5 * G2) - ((1.0 / 16.0) * (L * L));
        den = 1 + (.25 * L);
        qsqr = num / (den * den);

        // diffusion coefficent (equ 33)
        den = (qsqr - q0sqr) / (q0sqr * (1 + q0sqr));
        c[k] = 1.0 / (1.0 + den);

        // saturate diffusion coefficent
        if (c[k] < 0) {
          c[k] = 0;
        } else if (c[k] > 1) {
          c[k] = 1;
        }
      }
    }
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {

        // current index
        k = i * cols + j;

        // diffusion coefficent
        cN = c[k];
        cS = c[iS[i] * cols + j];
        cW = c[k];
        cE = c[i * cols + jE[j]];

        // divergence (equ 58)
        D = cN * dN[k] + cS * dS[k] + cW * dW[k] + cE * dE[k];

        // image update (equ 61)
        J[k] = J[k] + 0.25 * lambda * D;
      }
    }
#ifdef ITERATION
  }
#endif

#ifdef OUTPUT
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      printf("%.5f\n", J[i * cols + j]);
    }
    printf("\n");
  }
#endif
  free(I);
  free(J);
  free(iN);
  free(iS);
  free(jW);
  free(jE);
  free(dN);
  free(dS);
  free(dW);
  free(dE);
  free(c);
  return 0;
}

void random_matrix(float *I, int rows, int cols) {
  srand(7);
  for (int i = 0; i < rows; i++) {
    for (int j = 0; j < cols; j++) {
      I[i * cols + j] = rand() / (float)RAND_MAX;
    }
  }
}
