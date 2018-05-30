#define LIMIT -999
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>


#define BLOCK_SIZE 16
#define PENALTY 10

#ifdef SMALL_DATASET
#define DIMENSION 2048 
#else
#define DIMENSION 4096 
#endif

////////////////////////////////////////////////////////////////////////////////
// declaration, forward
void runTest();



int maximum(int a, int b, int c)
{
	int k;
	if( a <= b )
		k = b;
	else 
    	k = a;

	if( k <=c )
    	return(c);
	else
    	return(k);
}

int blosum62[24][24] = {
{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
{-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
{-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
{ 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
};


////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv)
{
    runTest();
    return EXIT_SUCCESS;
}

void usage(int argc, char **argv)
{
	fprintf(stderr, "Usage: %s <max_rows/max_cols> <penalty>\n", argv[0]);
	fprintf(stderr, "\t<dimension>      - x and y dimensions\n");
	fprintf(stderr, "\t<penalty>        - penalty(positive integer)\n");
	exit(1);
}

void nw_kernel(int * __restrict__ input_itemsets, int * __restrict__ output_itemsets, int * __restrict__ referrence, int max_rows, int max_cols, int penalty)
{
    {
        for( int blk = 1; blk <= (max_cols-1)/BLOCK_SIZE; blk++ )
        {
            for( int b_index_x = 0; b_index_x < blk; ++b_index_x)
            {
                int b_index_y = blk - 1 - b_index_x;
                int input_itemsets_l[(BLOCK_SIZE + 1) *(BLOCK_SIZE+1)] __attribute__ ((aligned (64)));
                int reference_l[BLOCK_SIZE * BLOCK_SIZE] __attribute__ ((aligned (64)));

                // Copy referrence to local memory
                for ( int i = 0; i < BLOCK_SIZE; ++i )
                {
                    for ( int j = 0; j < BLOCK_SIZE; ++j)
                    {
                        reference_l[i*BLOCK_SIZE + j] = referrence[max_cols*(b_index_y*BLOCK_SIZE + i + 1) + b_index_x*BLOCK_SIZE +  j + 1];
                    }
                }

                // Copy input_itemsets to local memory
                for ( int i = 0; i < BLOCK_SIZE + 1; ++i )
                {
                    for ( int j = 0; j < BLOCK_SIZE + 1; ++j)
                    {
                        input_itemsets_l[i*(BLOCK_SIZE + 1) + j] = input_itemsets[max_cols*(b_index_y*BLOCK_SIZE + i) + b_index_x*BLOCK_SIZE +  j];
                    }
                }

                // Compute
                for ( int i = 1; i < BLOCK_SIZE + 1; ++i )
                {
                    for ( int j = 1; j < BLOCK_SIZE + 1; ++j)
                    {
                        input_itemsets_l[i*(BLOCK_SIZE + 1) + j] = maximum( input_itemsets_l[(i - 1)*(BLOCK_SIZE + 1) + j - 1] + reference_l[(i - 1)*BLOCK_SIZE + j - 1],
                                input_itemsets_l[i*(BLOCK_SIZE + 1) + j - 1] - penalty,
                                input_itemsets_l[(i - 1)*(BLOCK_SIZE + 1) + j] - penalty);
                    }
                }

                // Copy results to global memory
                for ( int i = 0; i < BLOCK_SIZE; ++i )
                {
                    for ( int j = 0; j < BLOCK_SIZE; ++j)
                    {
                        input_itemsets[max_cols*(b_index_y*BLOCK_SIZE + i + 1) + b_index_x*BLOCK_SIZE +  j + 1] = input_itemsets_l[(i + 1)*(BLOCK_SIZE+1) + j + 1];
                    }
                }
                
            }
        }    
            
        // printf("Processing bottom-right matrix\n");

        for ( int blk = 2; blk <= (max_cols-1)/BLOCK_SIZE; blk++ )
        {
            for( int b_index_x = blk - 1; b_index_x < (max_cols-1)/BLOCK_SIZE; ++b_index_x)
            {
                int b_index_y = (max_cols-1)/BLOCK_SIZE + blk - 2 - b_index_x;

                int input_itemsets_l[(BLOCK_SIZE + 1) *(BLOCK_SIZE+1)] __attribute__ ((aligned (64)));
                int reference_l[BLOCK_SIZE * BLOCK_SIZE] __attribute__ ((aligned (64)));
     
                // Copy referrence to local memory
                for ( int i = 0; i < BLOCK_SIZE; ++i )
                {
                    for ( int j = 0; j < BLOCK_SIZE; ++j)
                    {
                        reference_l[i*BLOCK_SIZE + j] = referrence[max_cols*(b_index_y*BLOCK_SIZE + i + 1) + b_index_x*BLOCK_SIZE +  j + 1];
                    }
                }

                // Copy input_itemsets to local memory
                for ( int i = 0; i < BLOCK_SIZE + 1; ++i )
                {
                    for ( int j = 0; j < BLOCK_SIZE + 1; ++j)
                    {
                        input_itemsets_l[i*(BLOCK_SIZE + 1) + j] = input_itemsets[max_cols*(b_index_y*BLOCK_SIZE + i) + b_index_x*BLOCK_SIZE +  j];
                    }
                }

                // Compute
                for ( int i = 1; i < BLOCK_SIZE + 1; ++i )
                {
                    for ( int j = 1; j < BLOCK_SIZE + 1; ++j)
                    {
                        input_itemsets_l[i*(BLOCK_SIZE + 1) + j] = maximum( input_itemsets_l[(i - 1)*(BLOCK_SIZE + 1) + j - 1] + reference_l[(i - 1)*BLOCK_SIZE + j - 1],
                                input_itemsets_l[i*(BLOCK_SIZE + 1) + j - 1] - penalty,
                                input_itemsets_l[(i - 1)*(BLOCK_SIZE + 1) + j] - penalty);
                    }
                }
                // Copy results to global memory
                for ( int i = 0; i < BLOCK_SIZE; ++i )
                {

                    for ( int j = 0; j < BLOCK_SIZE; ++j)
                    {
                        input_itemsets[max_cols*(b_index_y*BLOCK_SIZE + i + 1) + b_index_x*BLOCK_SIZE +  j + 1] = input_itemsets_l[(i + 1)*(BLOCK_SIZE+1) + j +1];
                    }
                }
            }
        }   
    }
}
////////////////////////////////////////////////////////////////////////////////
//! Run a simple test
////////////////////////////////////////////////////////////////////////////////
void runTest() 
{
    int max_rows, max_cols, penalty;
    int * __restrict__ input_itemsets, * __restrict__ output_itemsets, * __restrict__ referrence;
    
    // the lengths of the two sequences should be able to divided by 16.
    // And at current stage  max_rows needs to equal max_cols
    
    max_rows = DIMENSION;
    max_cols = DIMENSION;
    penalty = PENALTY;


    max_rows = max_rows + 1;
    max_cols = max_cols + 1;

    referrence = (int *)malloc(max_rows * max_cols * sizeof(int) );
    input_itemsets = (int *)malloc( max_rows * max_cols * sizeof(int) );
    output_itemsets = (int *)malloc( max_rows * max_cols * sizeof(int) );

    if (!(input_itemsets && output_itemsets && referrence))
        fprintf(stderr, "Error in allocating memory to arrays");
    
    srand(7);
    
    for (int i = 0 ; i < max_cols; i++){
        for (int j = 0 ; j < max_rows; j++){
            input_itemsets[i*max_cols+j] = 0;
        }
    }

    for( int i=1; i< max_rows ; i++){   
        input_itemsets[i*max_cols] = rand() % 10 + 1;
    }

    for( int j=1; j< max_cols ; j++){   
        input_itemsets[j] = rand() % 10 + 1;
    }

    for (int i = 1 ; i < max_cols; i++){
        for (int j = 1 ; j < max_rows; j++){
            referrence[i*max_cols+j] = blosum62[input_itemsets[i*max_cols]][input_itemsets[j]];
        }
    }

    for( int i = 1; i< max_rows ; i++)
        input_itemsets[i*max_cols] = -i * penalty;

    for( int j = 1; j< max_cols ; j++)
        input_itemsets[j] = -j * penalty;



    //Compute top-left matrix 
  

    nw_kernel( input_itemsets, output_itemsets, referrence,
        max_rows, max_cols, penalty );


#define TRACEBACK
#ifdef TRACEBACK
    bool flag = false;
    for (int i = max_rows - 2,  j = max_rows - 2; i>=0 && j>=0;){
        int nw, n, w, traceback;
        if ( i == max_rows - 2 && j == max_rows - 2 )
            printf("%d ", input_itemsets[ i * max_cols + j]); //print the first element
        if ( i == 0 && j == 0 )
            break;
        if ( i > 0 && j > 0 ){
            nw = input_itemsets[(i - 1) * max_cols + j - 1];
            w  = input_itemsets[ i * max_cols + j - 1 ];
            n  = input_itemsets[(i - 1) * max_cols + j];
        }
        else if ( i == 0 ){
            nw = n = LIMIT;
            w  = input_itemsets[ i * max_cols + j - 1 ];
        }
        else if ( j == 0 ){
            nw = w = LIMIT;
            n  = input_itemsets[(i - 1) * max_cols + j];
        }
        else{
        }

        //traceback = maximum(nw, w, n);
        int new_nw, new_w, new_n;
        new_nw = nw + referrence[i * max_cols + j];
        new_w = w - penalty;
        new_n = n - penalty;

        traceback = maximum(new_nw, new_w, new_n);
        if(traceback == new_nw)
            traceback = nw;
        if(traceback == new_w)
            traceback = w;
        if(traceback == new_n)
            traceback = n;

        printf("%d ", traceback);
        flag = true;

        if(traceback == nw )
        {i--; j--; continue;}

        else if(traceback == w )
        {j--; continue;}

        else if(traceback == n )
        {i--; continue;}

        else
            ;
    }
    if(flag) printf("\n");

#endif
    free(referrence);
    free(input_itemsets);
    free(output_itemsets);

}



