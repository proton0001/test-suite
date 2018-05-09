#include <stdio.h>
#include <stdlib.h>
#include <jpeglib.h>
#include <iostream>
#include <cmath>

using namespace std;
#define PI 3.14159265

static inline int read_JPEG_file(char * __restrict__ filename, int * __restrict__ width, int * __restrict__ height, int *__restrict__ channels, unsigned char * __restrict__ (image[]))
{
    FILE *infile;
    if ((infile = fopen(filename, "rb")) == NULL) {
        fprintf(stderr, "can't open %s\n", filename);
        return 0;
    }
    struct jpeg_error_mgr jerr;
    struct jpeg_decompress_struct cinfo;
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    jpeg_stdio_src(&cinfo, infile);
    (void) jpeg_read_header(&cinfo, TRUE);
    (void) jpeg_start_decompress(&cinfo);

    *width = cinfo.output_width, *height = cinfo.output_height;
    *channels = cinfo.num_components;
    // printf("width=%d height=%d c=%d\n", *width, *height, *channels);
    *image = (unsigned char *)malloc(*width * *height * *channels * sizeof(*image));
    JSAMPROW rowptr[1];
    int row_stride = *width * *channels;

    while (cinfo.output_scanline < cinfo.output_height) {
        rowptr[0] = *image + row_stride * cinfo.output_scanline;
        jpeg_read_scanlines(&cinfo, rowptr, 1);
    }
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    fclose(infile);
    return 1;
}


/** Writes the image in the specified file.
  NOTE: works with Grayscale or RGB modes only (based on number of channels)
 */
static inline void write_JPEG_file(char * __restrict__ filename, int width, int height, int channels, unsigned char image[], int quality)
{
    FILE *outfile;
    if ((outfile = fopen(filename, "wb")) == NULL)
    {
        fprintf(stderr, "can't open %s\n", filename);
        exit(1);
    }
    struct jpeg_error_mgr jerr;
    struct jpeg_compress_struct cinfo;
    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);
    jpeg_stdio_dest(&cinfo,outfile);

    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = channels;
    cinfo.in_color_space = channels == 1 ? JCS_GRAYSCALE : JCS_RGB;
    jpeg_set_defaults(&cinfo);
    jpeg_set_quality(&cinfo, quality, TRUE);

    jpeg_start_compress(&cinfo, TRUE);
    JSAMPROW rowptr[1];
    int row_stride = width * channels;
    while (cinfo.next_scanline < cinfo.image_height) {
        rowptr[0] = & image[cinfo.next_scanline * row_stride];
        jpeg_write_scanlines(&cinfo, rowptr, 1);
    }
    jpeg_finish_compress(&cinfo);

    fclose(outfile);
    jpeg_destroy_compress(&cinfo);
}

void SobelEdgeDetection(char * __restrict__ filename, int height, int width,int  channels,unsigned char image[]);



int main(/*int argc, char *argv[]*/)
{
    // if (argc != 1) {
    //     printf("Usage: ./ex1\n");
    //     return 1;
    // }
    char *infile = "Sample.jpeg";
    char *outfile = "Sobel_Sample.jpeg";
    unsigned char *image;
    int width, height, channels;
    read_JPEG_file(infile, &width, &height, &channels, &image);

    //  If RBG then convert it to GrayScale Image
    if (channels == 3) {
        unsigned char *image2;
        image2 = (unsigned char *)malloc(width * height * channels * sizeof(image));
        channels = 1;
        for (int i=0; i<height; i++) {
            for (int j=0; j<width; j++) {
                image2[i*width+j] = (image[i*width*3+j*3] + image[i*width*3+j*3+1] + image[i*width*3+j*3+2])/3;
            }
        }
        free(image);
        image = image2;
        write_JPEG_file("Sample_Grayscale.jpg", width, height, channels, image, 95);
    }
    // return 0;
    SobelEdgeDetection(outfile , height, width, channels, image);
    // RobertsEdgeDetection(argv[3], height, width, channels, image);
    free(image);
    std::cout << "Output = Sobel_Sample.jpeg" << std::endl;
    return 0;

}

void SobelEdgeDetection(char *filename, int height, int width,int  channels,unsigned char image[])
{
    

    // int img2d[height][width]; // Seg fault for large image due to this
    std::cout << "Image Size = " << height << " x " << width << std::endl;
    int ** img2d = malloc(height*sizeof(int *));
    int ** img2dhororg = malloc(height*sizeof(int *));
    int ** img2dverorg = malloc(height*sizeof(int *));
    int ** img2dmag = malloc(height*sizeof(int *));
    for (int i=0; i<height; i++){
        img2d[i] = malloc(width*sizeof(int));
        for (int j=0; j<width; j++) {
            img2d[i][j] = image[i*width+j];
        }
    }


    for (int i=0; i<height; i++){
        img2dhororg[i] = malloc(width*sizeof(int));
        for (int j=0; j<width; j++) {
            img2dhororg[i][j] = 0;
        }
    }

    for (int i=0; i<height; i++){
        img2dverorg[i] = malloc(width*sizeof(int));
        for (int j=0; j<width; j++) {
            img2dverorg[i][j] = 0;
        }
    }

    for (int i=0; i<height; i++){
        img2dmag[i] = malloc(width*sizeof(int));
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

    for (int i=0; i<height; i++) {
        for (int j=0; j<width; j++) {
            image[i*width+j]=img2dmag[i][j];
            std::cout << img2dmag[i][j];
        }
    }
    write_JPEG_file(filename, width, height, channels, image, 95);
}