#ifndef _BACKPROP_H_
#define _BACKPROP_H_

#define BIGRND 0x7fffffff

#define SEED 7
#define ETA 0.3      // eta value
#define MOMENTUM 0.3 // momentum value

#ifdef SMALL_DATASET
#define LAYERSIZE 165536
#else
#define LAYERSIZE 555365
#endif

extern int layer_size;

int input_n;  /* number of input units */
int hidden_n; /* number of hidden units */
int output_n; /* number of output units */

float *__restrict__ input_units;  /* the input units */
float *__restrict__ hidden_units; /* the hidden units */
float *__restrict__ output_units; /* the output units */

float *__restrict__ hidden_delta; /* storage for hidden unit error */
float *__restrict__ output_delta; /* storage for output unit error */

float *__restrict__ target; /* storage for target vector */

float *__restrict__ input_weights;  /* weights from input to hidden layer */
float *__restrict__ hidden_weights; /* weights from hidden to output layer */

/*** The next two are for momentum ***/
float *__restrict__ input_prev_weights;  /* previous change on input to hidden
                                            wgt */
float *__restrict__ hidden_prev_weights; /* previous change on hidden to output
                                            wgt */

void bpnn_randomize_weights(float *w, int m, int n);
void bpnn_randomize_row(float *w, int m);
void bpnn_zero_weights(float *w, int m, int n);
void bpnn_internal_create(int n_in, int n_hidden, int n_out);
void bpnn_free();
void bpnn_create(int n_in, int n_hidden, int n_out);
void bpnn_layerforward(float *__restrict__ l1, float *__restrict__ l2,
                       float *__restrict__ conn, int n1, int n2);
void bpnn_output_error(float *__restrict__ delta, float *__restrict__ target,
                       float *__restrict__ output, int nj,
                       float *__restrict__ err);
void bpnn_hidden_error(float *__restrict__ delta_h, int nh,
                       float *__restrict__ delta_o, int no,
                       float *__restrict__ who, float *__restrict__ hidden,
                       float *__restrict__ err);
void bpnn_adjust_weights(float *__restrict__ delta, int ndelta,
                         float *__restrict__ ly, int nly, float *__restrict__ w,
                         float *__restrict__ oldw);
void bpnn_train(float *__restrict__ eo, float *__restrict__ eh);
void bpnn_dump();

#endif
