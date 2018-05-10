#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

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


void SobelEdgeDetection(int height, int width ,int **image);

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


    SobelEdgeDetection(HEIGHT, WIDTH, img2d);


    for (int i=0; i<HEIGHT; i++)
        free(img2d[i]);
    free(img2d); 

    return 0;

}

void SobelEdgeDetection(int height, int width ,int **img2d)
{
    

    // int img2d[height][width]; // Seg fault for large image due to this
    // std::cout << "Image Size = " << height << " x " << width << std::endl;
    int ** img2dhororg = (int **)malloc(height*sizeof(int *));
    int ** img2dverorg = (int **)malloc(height*sizeof(int *));
    int ** img2dmag = (int **)malloc(height*sizeof(int *));



    for (int i=0; i<height; i++){
        img2dhororg[i] = (int *)malloc(width*sizeof(int));
        for (int j=0; j<width; j++) {
            img2dhororg[i][j] = 0;
        }
    }

    for (int i=0; i<height; i++){
        img2dverorg[i] = (int *)malloc(width*sizeof(int));
        for (int j=0; j<width; j++) {
            img2dverorg[i][j] = 0;
        }
    }

    for (int i=0; i<height; i++){
        img2dmag[i] = (int *)malloc(width*sizeof(int));
        for (int j=0; j<width; j++) {
            img2dmag[i][j] = 0;
        }
    }

    // int img2dhororg[height][width];
    // int img2dverorg[height][width];
    // int img2dmag[height][width];




    ///horizontal
    int max=-200, min=2000;
    for (int i=1; i<height-1; i++){
        for (int j=1; j<width-1; j++) {
            int curr=img2d[i-1][j-1]+2*img2d[i-1][j]+img2d[i-1][j+1]-img2d[i+1][j-1]-2*img2d[i+1][j]-img2d[i+1][j+1];
            img2dhororg[i][j] = curr;
            if (curr>max)
                max = curr;
            if (curr<min)
                min = curr;
        }
    }


  ///vertical:
  max=-200; min=2000;

    for (int i=1; i<height-1; i++){
        for (int j=1; j<width-1; j++) {
            int curr=img2d[i-1][j-1]+2*img2d[i][j-1]+img2d[i+1][j-1]-img2d[i-1][j+1]-2*img2d[i][j+1]-img2d[i+1][j+1];
            img2dverorg[i][j] = curr;
            if (curr>max)
                max = curr;
            if (curr<min)
                min = curr;
        }
    }

  ///magnitude
  max=-200; min=2000;

    for (int i=0; i<height; i++) {
        for (int j=0; j<width; j++) {
            img2dmag[i][j] = sqrt(pow(img2dhororg[i][j], 2)+pow(img2dverorg[i][j], 2));
            if (img2dmag[i][j]>max)
                max = img2dmag[i][j];
            if (img2dmag[i][j]<min)
                min = img2dmag[i][j];
        }
    }

    int diff = max - min;
    for (int i=0; i<height; i++) {
        for (int j=0; j<width; j++){
            float abc = (img2dmag[i][j]-min)/(diff*1.0);
            img2dmag[i][j] = abc* 255;
        }
    }
    // Print Image
    for (int i=0; i<height; i++) {
        for (int j=0; j<width; j++) {
            DUMPIMAGE(img2dmag[i][j]);
        }
        free(img2dhororg[i]);
        free(img2dverorg[i]);
        free(img2dmag[i]);
    }

    free(img2dhororg);
    free(img2dverorg);
    free(img2dmag);

}