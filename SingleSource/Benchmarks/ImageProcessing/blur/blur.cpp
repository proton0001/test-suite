#include <iostream>
#include <cmath>
#include <stdlib.h>
using namespace std;

#ifdef SMALL_DATASET
#define HEIGHT 1920
#define WIDTH 1080
#else
#define HEIGHT 7680
#define WIDTH 4320
#endif




#define WINDOW 10
#define DUMPIMAGE(s) std::cout << s


void Blur(int height, int width, int **img2d);

int main(/*int argc, char *argv[]*/)
{
    int** img2d = (int**)malloc(HEIGHT*sizeof(int *));

    //  Initialize a random image  
    for (int i=0; i<HEIGHT; i++){
        img2d[i] = (int*)malloc(WIDTH*sizeof(int));
        for (int j=0; j<WIDTH; j++) {
            img2d[i][j] = (i*WIDTH*j +i+j)%256; //Any Random Arbitary Input Should Work
        }
    }

    Blur(HEIGHT, WIDTH, img2d);


    for (int i=0; i<HEIGHT; i++)
        free(img2d[i]);
    free(img2d);   

    return 0;
}

void Blur(int height, int width, int **img2d)
{
    
    int ** img2dblur = (int**)malloc(height*sizeof(int *));

    // Initializing output Image 
    for (int i=0; i<height; i++){
        img2dblur[i] = (int*)malloc(width*sizeof(int));
        for (int j=0; j<width; j++) {
            img2dblur[i][j] = 0;
        }
    }

    // Averaging over a window
    int sum = 0;
    int window_size = WINDOW;
    int offset = (window_size-1)/2;

    int max= -200, min=2000;
    for (int i=offset; i<height-offset; i++)
    {
        for (int j=offset; j<width-offset; j++)
        {
            // Computing Sum Over Window
            sum=0;
            for (int k= -1 * offset; k<offset; k++)
            {
                for (int l= -1 * offset; l<offset; l++)
                {
                    sum += img2d[i+k][j+l];
                }
            }
            // Averaging it
            img2dblur[i][j] = (sum)/(window_size*window_size);
            // Get Max and Min (to Normalize it between 0-255)
            if (img2dblur[i][j]>max)
                max = img2dblur[i][j];
            if (img2dblur[i][j]<min)
                min = img2dblur[i][j];
        }
    }

    //  To make it a bit sharper 
    int diff = max - min;
    for (int i=0; i<height; i++)
    {
        for (int j=0; j<width; j++)
        {
            float abc = (img2dblur[i][j]-min)/(diff*1.0);
            img2dblur[i][j] = abc* 255;
        }
    }

    // Print Image
    for (int i=0; i<height; i++) {
        for (int j=0; j<width; j++) {
            DUMPIMAGE(img2dblur[i][j]);
        }
        free(img2dblur[i]);
    }
    free(img2dblur);
}