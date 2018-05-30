/*
 ******************************************************************
 * HISTORY
 * 15-Oct-94  Jeff Shufelt (js), Carnegie Mellon University
 * Prepared for 15-681, Fall 1994.
 * Modified by Shuai Che
 * 28-May-2018: Modified by Pankaj Kukreja,
 * Indian Institute of Technology Hyderabad, India
 ******************************************************************
 */

#include "backprop.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// GAP Between array index, which we want to compare
// Otherwise very big reference output
#define GAP 100

#define ABS(x) (((x) > 0.0) ? (x) : (-(x)))

void bpnn_randomize_weights(float *w, int m, int n) {
  int i, j;
  for (i = 0; i <= m; i++) {
    for (j = 0; j <= n; j++) {
      w[i * n + i + j] = (float)rand() / RAND_MAX;
    }
  }
}

void bpnn_randomize_row(float *w, int m) {
  int i;
  for (i = 0; i <= m; i++) {
    w[i] = 0.1;
  }
}

void bpnn_zero_weights(float *w, int m, int n) {
  int i, j;
  for (i = 0; i <= m; i++) {
    for (j = 0; j <= n; j++) {
      w[i * n + i + j] = 0.0;
    }
  }
}

void bpnn_internal_create(int n_in, int n_hidden, int n_out) {

  input_n = n_in;
  hidden_n = n_hidden;
  output_n = n_out;

  input_units =
      /*alloc_1d_dbl(n_in + 1);*/ (float *)malloc(((n_in + 1) * sizeof(float)));
  if (input_units == NULL) {
    fprintf(stderr, "Couldn't allocate array of floats to input_units\n");
    exit (EXIT_FAILURE);
  }

  hidden_units = /*alloc_1d_dbl(n_hidden + 1); */ (float *)malloc(
      ((n_hidden + 1) * sizeof(float)));
  if (hidden_units == NULL) {
    fprintf(stderr, "Couldn't allocate array of floats to hidden_units\n");
    exit (EXIT_FAILURE);
  }

  output_units = /* alloc_1d_dbl(n_out + 1);*/ (float *)malloc(
      ((n_out + 1) * sizeof(float)));
  if (output_units == NULL) {
    fprintf(stderr, "Couldn't allocate array of floats to output_units\n");
    exit (EXIT_FAILURE);
  }

  hidden_delta = /*alloc_1d_dbl(n_hidden + 1);*/ (float *)malloc(
      ((n_hidden + 1) * sizeof(float)));
  if (hidden_delta == NULL) {
    fprintf(stderr, "Couldn't allocate array of floats to hidden_delta\n");
    exit (EXIT_FAILURE);
  }

  output_delta = /*alloc_1d_dbl(n_out + 1);*/ (float *)malloc(
      ((n_out + 1) * sizeof(float)));
  if (output_delta == NULL) {
    fprintf(stderr, "Couldn't allocate array of floats to output_delta\n");
    exit (EXIT_FAILURE);
  }

  target = /*alloc_1d_dbl(n_out + 1);*/ (float *)malloc(
      ((n_out + 1) * sizeof(float)));
  if (target == NULL) {
    fprintf(stderr, "Couldn't allocate array of floats to target\n");
    exit (EXIT_FAILURE);
  }

  input_weights = /*alloc_2d_dbl(n_in + 1, n_hidden + 1);*/ (float *)malloc(
      ((n_in + 1) * (n_hidden + 1) * sizeof(float)));
  if (input_weights == NULL) {
    fprintf(
        stderr,
        "ALLOC_2D_DBL: Couldn't allocate array of dbl ptrs input_weights \n");
    exit (EXIT_FAILURE);
  }

  hidden_weights = /*alloc_2d_dbl(n_hidden + 1, n_out + 1); */ (float *)malloc(
      ((n_hidden + 1) * (n_out + 1) * sizeof(float)));
  if (hidden_weights == NULL) {
    fprintf(stderr, " Couldn't allocate array of dbl ptrs to hidden_weights\n");
    exit (EXIT_FAILURE);
  }

  input_prev_weights = /* alloc_2d_dbl(n_in + 1, n_hidden + 1); */ (
      float *)malloc(((n_in + 1) * (n_hidden + 1) * sizeof(float)));
  if (input_prev_weights == NULL) {
    fprintf(stderr,
            " Couldn't allocate array of dbl ptrs to input_prev_weights\n");
    exit (EXIT_FAILURE);
  }

  hidden_prev_weights = /* alloc_2d_dbl(n_hidden + 1, n_out + 1); */ (
      float *)malloc(((n_hidden + 1) * (n_out + 1) * sizeof(float)));
  if (hidden_prev_weights == NULL) {
    fprintf(stderr,
            " Couldn't allocate array of dbl ptrs to hidden_prev_weights\n");
    exit (EXIT_FAILURE);
  }
}

void bpnn_free() {
  free(input_units);
  free(hidden_units);
  free(output_units);
  free(hidden_delta);
  free(output_delta);
  free(target);
  free(input_weights);
  free(input_prev_weights);
  free(hidden_weights);
  free(hidden_prev_weights);
}

/*** Creates a new fully-connected network from scratch,
     with the given numbers of input, hidden, and output units.
     Threshold units are automatically included.  All weights are
     randomly initialized.

     Space is also allocated for temporary storage (momentum weights,
     error computations, etc).
***/
void bpnn_create(int n_in, int n_hidden, int n_out) {

  bpnn_internal_create(n_in, n_hidden, n_out);

#ifdef INITZERO
  bpnn_zero_weights(input_weights, n_in, n_hidden);
#else
  bpnn_randomize_weights(input_weights, n_in, n_hidden);
#endif
  bpnn_randomize_weights(hidden_weights, n_hidden, n_out);
  bpnn_zero_weights(input_prev_weights, n_in, n_hidden);
  bpnn_zero_weights(hidden_prev_weights, n_hidden, n_out);
  bpnn_randomize_row(target, n_out);
}

// conn= [n1][n2] , l1[n1] , l2[n2]
void bpnn_layerforward(float *__restrict__ l1, float *__restrict__ l2,
                       float *__restrict__ conn, int n1, int n2) {
  float sum;
  int j, k;

  /*** Set up thresholding unit ***/
  l1[0] = 1.0;
  /*** For each unit in second layer ***/
  for (j = 1; j <= n2; j++) {
    /*** Compute weighted sum of its inputs ***/
    sum = 0.0;
    for (k = 0; k <= n1; k++) {
      sum += conn[k * n2 + k + j] * l1[k];
    }
    l2[j] = (1.0 / (1.0 + exp(-sum)));
  }
}

void bpnn_output_error(float *__restrict__ delta, float *__restrict__ target,
                       float *__restrict__ output, int nj,
                       float *__restrict__ err) {
  int j;
  float o, t, errsum;
  errsum = 0.0;
  for (j = 1; j <= nj; j++) {
    o = output[j];
    t = target[j];
    delta[j] = o * (1.0 - o) * (t - o);
    errsum += ABS(delta[j]);
  }
  *err = errsum;
}

// delta_h[nh], delta_o[no], who[nh][no], err = int, hidden[nh]
void bpnn_hidden_error(float *__restrict__ delta_h, int nh,
                       float *__restrict__ delta_o, int no,
                       float *__restrict__ who, float *__restrict__ hidden,
                       float *__restrict__ err) {
  int j, k;
  float h, sum, errsum;
  errsum = 0.0;
  for (j = 1; j <= nh; j++) {
    h = hidden[j];
    sum = 0.0;
    for (k = 1; k <= no; k++) {
      sum += delta_o[k] * who[j * no + j + k];
    }
    delta_h[j] = h * (1.0 - h) * sum;
    errsum += ABS(delta_h[j]);
  }
  *err = errsum;
}

// delta[ndelta], ly[nly], w[nly][ndelta], oldw[nly][ndelta]
void bpnn_adjust_weights(float *__restrict__ delta, int ndelta,
                         float *__restrict__ ly, int nly, float *__restrict__ w,
                         float *__restrict__ oldw) {
  float new_dw;
  int k, j;
  ly[0] = 1.0;

  for (j = 1; j <= ndelta; j++) {
    for (k = 0; k <= nly; k++) {
      new_dw =
          ((ETA * delta[j] * ly[k]) + (MOMENTUM * oldw[k * ndelta + k + j]));
      w[k * ndelta + k + j] += new_dw;
      oldw[k * ndelta + k + j] = new_dw;
    }
  }
}

void bpnn_train(float *__restrict__ eo, float *__restrict__ eh) {
  int in, hid, out;
  float out_err, hid_err;

  in = input_n;
  hid = hidden_n;
  out = output_n;

  /*** Feed forward input activations. ***/
  bpnn_layerforward(input_units, hidden_units, input_weights, in, hid);
  bpnn_layerforward(hidden_units, output_units, hidden_weights, hid, out);

  /*** Compute error on output and hidden units. ***/
  bpnn_output_error(output_delta, target, output_units, out, &out_err);
  bpnn_hidden_error(hidden_delta, hid, output_delta, out, hidden_weights,
                    hidden_units, &hid_err);

  *eo = out_err;
  *eh = hid_err;

  /*** Adjust input and hidden weights. ***/
  bpnn_adjust_weights(output_delta, out, hidden_units, hid, hidden_weights,
                      hidden_prev_weights);
  bpnn_adjust_weights(hidden_delta, hid, input_units, in, input_weights,
                      input_prev_weights);
}

void bpnn_dump() {
  int n1, n2, n3, i, j;

  n1 = input_n;
  n2 = hidden_n;
  n3 = output_n;
  fflush(stdout);

  for (i = 0; i <= n1; i++) {
    for (j = 0; j <= n2; j++) {
      if((i*n2+i+j) % GAP ==  0)
      {
        printf("%.6f\n", input_weights[i * n2 + i + j]);
      }
    }
  }


  for (i = 0; i <= n2; i++) {
    for (j = 0; j <= n3; j++) {
      if((i*n3+i+j) % GAP ==  0)
      {
        printf("%.6f\n", hidden_weights[i * n3 + i + j]);
      }
    }
  }
  return;
}