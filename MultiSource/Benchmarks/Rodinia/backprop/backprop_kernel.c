#include "backprop.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

int layer_size = 0;

////////////////////////////////////////////////////////////////////////////////

extern void bpnn_layerforward(float *__restrict__ l1, float *__restrict__ l2,
                              float *__restrict__ conn, int n1, int n2);

extern void bpnn_output_error(float *__restrict__ delta,
                              float *__restrict__ target,
                              float *__restrict__ output, int nj,
                              float *__restrict__ err);

extern void bpnn_hidden_error(float *__restrict__ delta_h, int nh,
                              float *__restrict__ delta_o, int no,
                              float *__restrict__ who,
                              float *__restrict__ hidden,
                              float *__restrict__ err);

extern void bpnn_adjust_weights(float *__restrict__ delta, int ndelta,
                                float *__restrict__ ly, int nly,
                                float *__restrict__ w,
                                float *__restrict__ oldw);

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main() {

  layer_size = LAYERSIZE;

  // Seed for random number generator
  srand(SEED);

  int i;
  float out_err, hid_err;
  bpnn_create(layer_size, 16, 1); // (16, 1 can not be changed)

  // printf("Input layer size : %d\n", layer_size);
  int k;
  k = 1;
  for (i = 0; i < layer_size; i++) {
    input_units[k] = (float)rand() / RAND_MAX;
    k++;
  }

  // entering the training kernel, only one iteration
  bpnn_train_kernel(&out_err, &hid_err);

  bpnn_dump();

  bpnn_free();
}

void bpnn_train_kernel() {
  int in, hid, out;
  float out_err, hid_err;

  in = input_n;
  hid = hidden_n;
  out = output_n;

  // printf("Performing CPU computation\n");
  bpnn_layerforward(input_units, hidden_units, input_weights, in, hid);
  bpnn_layerforward(hidden_units, output_units, hidden_weights, hid, out);
  bpnn_output_error(output_delta, target, output_units, out, &out_err);
  bpnn_hidden_error(hidden_delta, hid, output_delta, out, hidden_weights,
                    hidden_units, &hid_err);
  bpnn_adjust_weights(output_delta, out, hidden_units, hid, hidden_weights,
                      hidden_prev_weights);
  bpnn_adjust_weights(hidden_delta, hid, input_units, in, input_weights,
                      input_prev_weights);
  printf("Output Error: %.6f\nHidden Error: %.6f\n", out_err, hid_err);
}
