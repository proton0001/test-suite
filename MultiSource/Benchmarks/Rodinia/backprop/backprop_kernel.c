#include "backprop.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

int layer_size = 0;

void bpnn_train_kernel();
////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main() {

  layer_size = LAYERSIZE;

  // Seed for random number generator
  srand(SEED);

  int i;
  bpnn_create(layer_size, 16, 1); // (16, 1 can not be changed)

  int k;
  k = 1;
  for (i = 0; i < layer_size; i++) {
    input_units[k] = (float)rand() / RAND_MAX;
    k++;
  }

  // entering the training kernel, only one iteration
  bpnn_train_kernel();

  bpnn_dump();

  bpnn_free();
}

void bpnn_train_kernel() {
  int in, hid, out;
  float out_err, hid_err;

  in = input_n;
  hid = hidden_n;
  out = output_n;

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
