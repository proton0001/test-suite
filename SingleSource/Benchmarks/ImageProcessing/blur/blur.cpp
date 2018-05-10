#include <iostream>
#include <cmath>
#include <stdlib.h>

#ifdef SMALL_DATASET
#define HEIGHT 1920
#define WIDTH 1080
#else
#define HEIGHT 7680
#define WIDTH 4320
#endif


#define WINDOW 10
#define PRINTIMAGE(s) std::cout << s




void init_image(int height, int width, int **image);
void kernel_box_blur(int height, int width, int **image);
void print_image(int height, int width, int **image);



int main(int argc, char *argv[])
{
    int ** image = (int**)malloc(HEIGHT*sizeof(int *));

    init_image(HEIGHT, WIDTH, image);
    kernel_box_blur(HEIGHT, WIDTH, image);

    print_image( HEIGHT,  WIDTH,  image);
    
    for (int i=0; i<HEIGHT; i++)
        free(image[i]);
    free(image);   

    return EXIT_SUCCESS;
}



void init_image(int height, int width, int **image)
{
    //  Initialize a random image  
    for (int i=0; i<HEIGHT; i++){
        image[i] = (int*)malloc(WIDTH*sizeof(int));
        for (int j=0; j<WIDTH; j++) {
            image[i][j] = (i*j+i+j)%256; //Any Random Arbitary Input Should Work
        }
    }
}

void kernel_box_blur(int height, int width, int **image)
{
    // Allocating memory for output image
    int ** img2dblur = (int**)malloc(height*sizeof(int *));

    // Initializing output Image 
    for (int i=0; i<height; i++){
        img2dblur[i] = (int*)malloc(width*sizeof(int));
        for (int j=0; j<width; j++) {
            img2dblur[i][j] = 0;
        }
    }

    int sum_in_window = 0;
    int window_size = WINDOW;
    int offset = (window_size-1)/2;
    int n   = WINDOW*WINDOW;
    int max = -200;
    int min = 2000;

    for (int i=offset; i<height-offset; i++)
    {
        for (int j=offset; j<width-offset; j++)
        {
            /* Computing sum of elements in window centered at i,j */
            sum_in_window=0;
            for (int k= -1 * offset; k<offset; k++)
            {
                for (int l= -1 * offset; l<offset; l++)
                {
                    sum_in_window += image[i+k][j+l];
                }
            }
            /* Averaging it */
            img2dblur[i][j] = (sum_in_window)/(n);
            /* Get Max and Min (to Scale it later between 0-255) */
            if (img2dblur[i][j]>max)
                max = img2dblur[i][j];
            if (img2dblur[i][j]<min)
                min = img2dblur[i][j];
        }
    }

    // Scale everything from 0-255
    int diff = max - min;

    // if max = min then image is constant all over hence no edges(0 pixels only)
    if(diff==0)
    {
        diff = 1;
    }

    for (int i=0; i<height; i++)
    {
        for (int j=0; j<width; j++)
        {
            float abc = (img2dblur[i][j]-min)/(diff*1.0);
            img2dblur[i][j] = abc* 255;
        }
    }


    for (int i=0; i<height; i++)
    {
        for (int j=0; j<width; j++)
        {
            image[i][j] =  img2dblur[i][j];
        }
    }


    // Clear allocated space for image
    for (int i=0; i<height; i++) {
        free(img2dblur[i]);
    }
    free(img2dblur);
}

void print_image(int height, int width, int **image)
{    
    // Print Image
    for (int i=0; i<height; i++) {
        for (int j=0; j<width; j++) {
            PRINTIMAGE(image[i][j]);
        }
    }
}