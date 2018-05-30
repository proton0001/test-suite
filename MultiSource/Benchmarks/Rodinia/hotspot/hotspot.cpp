/*
 * Modified by Pankaj Kukreja
 * Indian Institute of Technology Hyderabad, India
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <iostream>
#include <fstream>
#include <cstdlib>

// Gap in which array should be printed
#define GAP 25


#ifdef SMALL_DATASET
  #define IN_SIZE 128
  #define MULTIPLIER 4 
#else
  #define IN_SIZE 1024
  #define MULTIPLIER 2
#endif



#define OUT_SIZE IN_SIZE*MULTIPLIER
#define ITERATIONCOUNT 4

using namespace std;

#define BLOCK_SIZE 16
#define BLOCK_SIZE_C BLOCK_SIZE
#define BLOCK_SIZE_R BLOCK_SIZE
#define STR_SIZE 256

/* maximum power density possible (say 300W for a 10mm x 10mm chip) */
#define MAX_PD (3.0e6)

/* required precision in degrees    */
#define PRECISION 0.001
#define SPEC_HEAT_SI 1.75e6
#define K_SI 100

/* capacitance fitting factor   */
#define FACTOR_CHIP 0.5

typedef float FLOAT;

/* chip parameters  */
const FLOAT t_chip = 0.0005;
const FLOAT chip_height = 0.016;
const FLOAT chip_width = 0.016;

/* ambient temperature, assuming no package at all  */
const FLOAT amb_temp = 80.0;

/* Single iteration of the transient solver in the grid model.
 * advances the solution of the discretized difference equations
 * by one time step
 */
void single_iteration(FLOAT *__restrict__ result, FLOAT *__restrict__ temp,
                      FLOAT *__restrict__ power, int row, int col, FLOAT Cap_1,
                      FLOAT Rx_1, FLOAT Ry_1, FLOAT Rz_1, FLOAT step,
                      int num_iterations) {
{
    FLOAT *tmp;
    for (int i = 0; i < num_iterations; i++) {
      /*++++++++++++++++++++++++++++++*/
      FLOAT delta;
      int r, c;
      int chunk;
      int num_chunk = row * col / (BLOCK_SIZE_R * BLOCK_SIZE_C);
      int chunks_in_row = col / BLOCK_SIZE_C;
      int chunks_in_col = row / BLOCK_SIZE_R;

      for (chunk = 0; chunk < num_chunk; ++chunk) {
        int r_start = BLOCK_SIZE_R * (chunk / chunks_in_col);
        int c_start = BLOCK_SIZE_C * (chunk % chunks_in_row);
        int r_end = r_start + BLOCK_SIZE_R > row ? row : r_start + BLOCK_SIZE_R;
        int c_end = c_start + BLOCK_SIZE_C > col ? col : c_start + BLOCK_SIZE_C;

        if (r_start == 0 || c_start == 0 || r_end == row || c_end == col) {
          for (r = r_start; r < r_start + BLOCK_SIZE_R; ++r) {
            for (c = c_start; c < c_start + BLOCK_SIZE_C; ++c) {
              /* Corner 1 */
              if ((r == 0) && (c == 0)) {
                delta = (Cap_1) * (power[0] + (temp[1] - temp[0]) * Rx_1 +
                                   (temp[col] - temp[0]) * Ry_1 +
                                   (amb_temp - temp[0]) * Rz_1);
              } /* Corner 2 */
              else if ((r == 0) && (c == col - 1)) {
                delta = (Cap_1) * (power[c] + (temp[c - 1] - temp[c]) * Rx_1 +
                                   (temp[c + col] - temp[c]) * Ry_1 +
                                   (amb_temp - temp[c]) * Rz_1);
              } /* Corner 3 */
              else if ((r == row - 1) && (c == col - 1)) {
                delta = (Cap_1) *
                        (power[r * col + c] +
                         (temp[r * col + c - 1] - temp[r * col + c]) * Rx_1 +
                         (temp[(r - 1) * col + c] - temp[r * col + c]) * Ry_1 +
                         (amb_temp - temp[r * col + c]) * Rz_1);
              } /* Corner 4 */
              else if ((r == row - 1) && (c == 0)) {
                delta =
                    (Cap_1) * (power[r * col] +
                               (temp[r * col + 1] - temp[r * col]) * Rx_1 +
                               (temp[(r - 1) * col] - temp[r * col]) * Ry_1 +
                               (amb_temp - temp[r * col]) * Rz_1);
              } /* Edge 1 */
              else if (r == 0) {
                delta = (Cap_1) *
                        (power[c] +
                         (temp[c + 1] + temp[c - 1] - 2.0 * temp[c]) * Rx_1 +
                         (temp[col + c] - temp[c]) * Ry_1 +
                         (amb_temp - temp[c]) * Rz_1);
              } /* Edge 2 */
              else if (c == col - 1) {
                delta = (Cap_1) *
                        (power[r * col + c] +
                         (temp[(r + 1) * col + c] + temp[(r - 1) * col + c] -
                          2.0 * temp[r * col + c]) *
                             Ry_1 +
                         (temp[r * col + c - 1] - temp[r * col + c]) * Rx_1 +
                         (amb_temp - temp[r * col + c]) * Rz_1);
              } /* Edge 3 */
              else if (r == row - 1) {
                delta = (Cap_1) *
                        (power[r * col + c] +
                         (temp[r * col + c + 1] + temp[r * col + c - 1] -
                          2.0 * temp[r * col + c]) *
                             Rx_1 +
                         (temp[(r - 1) * col + c] - temp[r * col + c]) * Ry_1 +
                         (amb_temp - temp[r * col + c]) * Rz_1);
              } /* Edge 4 */
              else if (c == 0) {
                delta = (Cap_1) * (power[r * col] +
                                   (temp[(r + 1) * col] + temp[(r - 1) * col] -
                                    2.0 * temp[r * col]) *
                                       Ry_1 +
                                   (temp[r * col + 1] - temp[r * col]) * Rx_1 +
                                   (amb_temp - temp[r * col]) * Rz_1);
              }
              result[r * col + c] = temp[r * col + c] + delta;
            }
          }
          continue;
        }

        for (r = r_start; r < r_start + BLOCK_SIZE_R; ++r) {
          for (c = c_start; c < c_start + BLOCK_SIZE_C; ++c) {
            /* Update Temperatures */
            result[r * col + c] =
                temp[r * col + c] +
                (Cap_1 * (power[r * col + c] +
                          (temp[(r + 1) * col + c] + temp[(r - 1) * col + c] -
                           2.f * temp[r * col + c]) *
                              Ry_1 +
                          (temp[r * col + c + 1] + temp[r * col + c - 1] -
                           2.f * temp[r * col + c]) *
                              Rx_1 +
                          (amb_temp - temp[r * col + c]) * Rz_1));
          }
        }
      }
      /*++++++++++++++++++++++++++++++*/
      tmp = temp;
      temp = result;
      result = tmp;
    }
  }
}

/*
 * Transient solver driver routine: simply converts the heat
 * transfer differential equations to difference equations
 * and solves the difference equations by iterating
 */
void compute_tran_temp(FLOAT *__restrict__ result, int num_iterations,
                       FLOAT *__restrict__ temp, FLOAT *__restrict__ power,
                       int row, int col) {
  FLOAT grid_height = chip_height / row;
  FLOAT grid_width = chip_width / col;

  FLOAT Cap = FACTOR_CHIP * SPEC_HEAT_SI * t_chip * grid_width * grid_height;
  FLOAT Rx = grid_width / (2.0 * K_SI * t_chip * grid_height);
  FLOAT Ry = grid_height / (2.0 * K_SI * t_chip * grid_width);
  FLOAT Rz = t_chip / (K_SI * grid_height * grid_width);

  FLOAT max_slope = MAX_PD / (FACTOR_CHIP * t_chip * SPEC_HEAT_SI);
  FLOAT step = PRECISION / max_slope / 1000.0;

  FLOAT Rx_1 = 1.f / Rx;
  FLOAT Ry_1 = 1.f / Ry;
  FLOAT Rz_1 = 1.f / Rz;
  FLOAT Cap_1 = step / Cap;
  //
  {
    single_iteration(result, temp, power, row, col, Cap_1, Rx_1, Ry_1, Rz_1,
                     step, num_iterations);
  }
}

void writeoutput(FLOAT *vect, int grid_rows, int grid_cols) {
  int i, j, index = 0;
  for (i = 0; i < grid_rows; i++) {
    for (j = 0; j < grid_cols; j++) {
      if(index%GAP==0)
      {
        printf("%d\t%g\n", index, vect[i * grid_cols + j]);
      }
      index++;
    }
  }
}

void read_input(FLOAT * __restrict__ vect, char * __restrict__ file) {

  const int x = MULTIPLIER;

  double val;

  fstream fs;
  
  // copy values into larger array
  fs.open(file, ios::in );
  if ( !fs )
    std::cerr << "Failed to open input file.\n";

    for ( int row = 0; row < IN_SIZE; row++ )
    {
      for ( int col = 0; col < IN_SIZE; col++ )
      {
        fs >> val;
        for ( int rowOff = 0; rowOff < x; rowOff++ )
        {
          for ( int colOff = 0; colOff < x; colOff++ )
          {
            vect[(x * row + rowOff) * OUT_SIZE + (x * col + colOff)] = val;
          }
        }
      }
    }
  fs.close();
}

void usage(int argc, char **argv) {
  fprintf( stderr, "Usage: %s <init_file> <power_file>\n", argv[0]);
  fprintf(stderr, "\t<temp_file>  - name of the file containing the initial temperature values of each cell\n");
  fprintf(stderr, "\t<power_file> - name of the file containing the dissipated power values of each cell\n");
  exit(1);
}

int main(int argc, char **argv) {
  int grid_rows, grid_cols, sim_time;
  FLOAT * __restrict__ temp, * __restrict__ power, * __restrict__ result;
  char *__restrict__ tfile, *__restrict__ pfile;

  /* check validity of inputs */
  if (argc != 3) 
    usage(argc, argv);

  grid_cols = OUT_SIZE;
  grid_rows = OUT_SIZE;
  sim_time = ITERATIONCOUNT; 


  /* allocate memory for the temperature and power arrays */
  temp = (FLOAT *)calloc(grid_rows * grid_cols, sizeof(FLOAT));
  power = (FLOAT *)calloc(grid_rows * grid_cols, sizeof(FLOAT));
  result = (FLOAT *)calloc(grid_rows * grid_cols, sizeof(FLOAT));

  if (!temp || !power) fprintf(stderr, " unable to allocate memory\n");

  /* read initial temperatures and input power    */
  tfile = argv[1]; // temp_file
  pfile = argv[2]; // Power_file

  read_input(temp, tfile);
  read_input(power, pfile);


  // long long start_time = get_time();
  compute_tran_temp(result, sim_time, temp, power, grid_rows, grid_cols);
  // long long end_time = get_time();
  // printf("Total time: %.3f seconds\n", ((float) (end_time - start_time)) / (1000*1000));

  // Since we are swapping res and temp array in function, we need to check which was final result array
  writeoutput((1 & sim_time) ? result : temp, grid_rows, grid_cols);

  /* cleanup  */
  free(temp);
  free(power);
  return 0;
}

