#ifndef _BACKPROP_H_
#define _BACKPROP_H_

#define BIGRND 0x7fffffff

#define SEED 7
#define ETA 0.3       //eta value
#define MOMENTUM 0.3  //momentum value


#ifdef SMALL_DATASET
#define LAYERSIZE 165536
#else
#define LAYERSIZE 555365 
#endif

extern int layer_size	;

int input_n;                  /* number of input units */
int hidden_n;                 /* number of hidden units */
int output_n;                 /* number of output units */

float *__restrict__ input_units;          /* the input units */
float *__restrict__ hidden_units;         /* the hidden units */
float *__restrict__ output_units;         /* the output units */

float *__restrict__ hidden_delta;         /* storage for hidden unit error */
float *__restrict__ output_delta;         /* storage for output unit error */

float *__restrict__ target;               /* storage for target vector */

float *__restrict__ input_weights;       /* weights from input to hidden layer */
float *__restrict__ hidden_weights;      /* weights from hidden to output layer */

                            /*** The next two are for momentum ***/
float *__restrict__ input_prev_weights;  /* previous change on input to hidden wgt */
float *__restrict__ hidden_prev_weights; /* previous change on hidden to output wgt */

void bpnn_create();
void bpnn_free();

void bpnn_train();
void bpnn_feedforward();

void bpnn_save();


#endif
