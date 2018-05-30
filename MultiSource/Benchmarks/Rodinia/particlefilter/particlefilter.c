/**
 * @file particlefilter.c
 * @author Michael Trotter & Matt Goodrum
 * 
 * Modified by Pankaj Kukreja, Indian Institute of technology, Hyderabad
 */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include <string.h>
#define PI 3.1415926535897932

#ifdef SMALL_DATASET
#define DIMX 128
#define DIMY 128
#define NFRAME 10
#define NPART 10000
#else
#define DIMX 256
#define DIMY 256
#define NFRAME 20
#define NPART 10000
#endif


/**
@var M value for Linear Congruential Generator (LCG); use GCC's value
*/
long M = INT_MAX;
/**
@var A value for LCG
*/
int A = 1103515245;
/**
@var C value for LCG
*/
int C = 12345;



/*****************************

* Takes in a double and returns an integer that approximates to that double
* @return if the mantissa < .5 => return value < input value; else return value > input value
*/
double roundDouble(double value){
	int newValue = (int)(value);
	if(value - newValue < .5)
		return newValue;
	else
		return newValue++;
}
/**
* Set values of the 3D array to a newValue if that value is equal to the testValue
* @param testValue The value to be replaced
* @param newValue The value to replace testValue with
* @param array3D The image vector
* @param dimX The x dimension of the frame
* @param dimY The y dimension of the frame
* @param dimZ The number of frames
*/
void setIf(int testValue, int newValue, int * array3D, int dimX, int dimY, int dimZ){
	int x, y, z;
	for(x = 0; x < dimX; x++) {
		for(y = 0; y < dimY; y++) {
			for(z = 0; z < dimZ; z++) {
				if(array3D[x * dimY * dimZ + y * dimZ + z] == testValue)
					array3D[x * dimY * dimZ + y * dimZ + z] = newValue;
			}
		}
	}
}
/**
* Generates a uniformly distributed random number using the provided seed and GCC's settings for the Linear Congruential Generator (LCG)
* @see http://en.wikipedia.org/wiki/Linear_congruential_generator
* @note This function is thread-safe
* @param seed The seed array
* @param index The specific index of the seed to be advanced
* @return a uniformly distributed number [0, 1)
*/
double randu(int * seed, int index)
{
	int num = A*seed[index] + C;
	seed[index] = num % M;
	return fabs(seed[index]/((double) M));
}
/**
* Generates a normally distributed random number using the Box-Muller transformation
* @note This function is thread-safe
* @param seed The seed array
* @param index The specific index of the seed to be advanced
* @return a double representing random number generated using the Box-Muller algorithm
* @see http://en.wikipedia.org/wiki/Normal_distribution, section computing value for normal random distribution
*/
double randn(int * seed, int index){
	/*Box-Muller algorithm*/
	double u = randu(seed, index);
	double v = randu(seed, index);
	double cosine = cos(2*PI*v);
	double rt = -2*log(u);
	return sqrt(rt)*cosine;
}
/**
* Sets values of 3D matrix using randomly generated numbers from a normal distribution
* @param array3D The video to be modified
* @param dimX The x dimension of the frame
* @param dimY The y dimension of the frame
* @param dimZ The number of frames
* @param seed The seed array
*/
void addNoise(int * __restrict__ array3D, int dimX, int dimY, int dimZ, int * __restrict__ seed){
	int x, y, z;
	for(x = 0; x < dimX; x++){
		for(y = 0; y < dimY; y++){
			for(z = 0; z < dimZ; z++){
				array3D[x * dimY * dimZ + y * dimZ + z] = array3D[x * dimY * dimZ + y * dimZ + z] + (int)(5*randn(seed, 0));
			}
		}
	}
}
/**
* Dilates the target matrix using the radius as a guide
* @param matrix The reference matrix
* @param dimX The x dimension of the video
* @param dimY The y dimension of the video
* @param dimZ The z dimension of the video
* @param error The error radius to be dilated
* @param newMatrix The target matrix
*/
void imdilate_disk(int * __restrict__ matrix, int dimX, int dimY, int dimZ, int error, int * __restrict__ newMatrix)
{
	int x, y, z;
	for(z = 0; z < dimZ; z++){
		for(x = 0; x < dimX; x++){
			for(y = 0; y < dimY; y++){
				if(matrix[x*dimY*dimZ + y*dimZ + z] == 1){
					int startX = x - error;
					while(startX < 0)
						startX++;
					int startY = y - error;
					while(startY < 0)
						startY++;
					int endX = x + error;
					while(endX > dimX)
						endX--;
					int endY = y + error;
					while(endY > dimY)
						endY--;
					int i, j;
					for(i = startX; i < endX; i++){
						for(j = startY; j < endY; j++){
							double distance = sqrt( pow((double)(i-x),2) + pow((double)(j-y),2) );
							if(distance < error)
								newMatrix[i*dimY*dimZ + j*dimZ + z] = 1;
						}
					}
				}
			}
		}
	}
}

/**
* Fills a 2D array describing the offsets of the disk object
* @param se The disk object
* @param numOnes The number of ones in the disk
* @param neighbors The array that will contain the offsets
* @param radius The radius used for dilation
*/
void getneighbors(int * __restrict__ se, int numOnes, double * __restrict__ neighbors, int radius) {
	int x, y;
	int neighY = 0;
	int center = radius - 1;
	int diameter = radius*2 -1;
	for(x = 0; x < diameter; x++) {
		for(y = 0; y < diameter; y++) {
			if(se[x*diameter + y]) {
				neighbors[neighY*2] = (int)(y - center);
				neighbors[neighY*2 + 1] = (int)(x - center);
				neighY++;
			}
		}
	}
}
/**
* The synthetic video sequence we will work with here is composed of a
* single moving object, circular in shape (fixed radius)
* The motion here is a linear motion
* the foreground intensity and the backgrounf intensity is known
* the image is corrupted with zero mean Gaussian noise
* @param I The video itself
* @param IszX The x dimension of the video
* @param IszY The y dimension of the video
* @param Nfr The number of frames of the video
* @param seed The seed array used for number generation
*/
void videoSequence(int * __restrict__ I, int IszX, int IszY, int Nfr, int * __restrict__ seed){
	int k;
	int max_size = IszX*IszY*Nfr;
	/*get object centers*/
	int x0 = (int)roundDouble(IszY/2.0);
	int y0 = (int)roundDouble(IszX/2.0);
	I[x0 *IszY *Nfr + y0 * Nfr  + 0] = 1;
	
	/*move point*/
	int xk, yk, pos;
	for(k = 1; k < Nfr; k++){
		xk = abs(x0 + (k-1));
		yk = abs(y0 - 2*(k-1));
		pos = yk * IszY * Nfr + xk *Nfr + k;
		if(pos >= max_size)
			pos = 0;
		I[pos] = 1;
	}
	
	/*dilate matrix*/
	int * newMatrix = (int *)malloc(sizeof(int)*IszX*IszY*Nfr);
	imdilate_disk(I, IszX, IszY, Nfr, 5, newMatrix);
	
	int x, y;
	for(x = 0; x < IszX; x++) {
		for(y = 0; y < IszY; y++) {
			for(k = 0; k < Nfr; k++) {
				I[x*IszY*Nfr + y*Nfr + k] = newMatrix[x*IszY*Nfr + y*Nfr + k];
			}
		}
	}
	free(newMatrix);
	/*define background, add noise*/
	setIf(0, 100, I, IszX, IszY, Nfr);
	setIf(1, 228, I, IszX, IszY, Nfr);
	/*add noise*/
	addNoise(I, IszX, IszY, Nfr, seed);
}


/**
* Finds the first element in the CDF that is greater than or equal to the provided value and returns that index
* @note This function uses sequential search
* @param CDF The CDF
* @param lengthCDF The length of CDF
* @param value The value to be found
* @return The index of value in the CDF; if value is never found, returns the last index
*/
int findIndex(double * CDF, int lengthCDF, double value){
	int index = -1;
	int x;
	for(x = 0; x < lengthCDF; x++){
		if(CDF[x] >= value){
			index = x;
			break;
		}
	}
	if(index == -1){
		return lengthCDF-1;
	}
	return index;
}


/**
* The implementation of the particle filter using OpenMP for many frames
* @see http://openmp.org/wp/
* @note This function is designed to work with a video of several frames. In addition, it references a provided MATLAB function which takes the video, the objxy matrix and the x and y arrays as arguments and returns the likelihoods
* @param I The video to be run
* @param IszX The x dimension of the video
* @param IszY The y dimension of the video
* @param Nfr The number of frames
* @param seed The seed array used for random number generation
* @param Nparticles The number of particles to be used
*/
void particleFilter(int * __restrict__ I, int IszX, int IszY, int Nfr, int * __restrict__ seed, int Nparticles)
{
	int max_size = IszX*IszY*Nfr;
	//original particle centroid
	double xe = roundDouble(IszY/2.0);
	double ye = roundDouble(IszX/2.0);
	
	//expected object locations, compared to center
	int radius = 5;
	int diameter = radius*2 - 1;
	int * disk = (int *)malloc(diameter*diameter*sizeof(int));
	int _x, _y;
	for(_x = 0; _x < diameter; _x++){
		for(_y = 0; _y < diameter; _y++){
			double _distance = sqrt(pow((double)(_x-radius+1),2) + pow((double)(_y-radius+1),2));
			if(_distance < radius)
				disk[_x*diameter + _y] = 1;
		}
	}

	int countOnes = 0;
	int x, y;
	for(x = 0; x < diameter; x++){
		for(y = 0; y < diameter; y++){
			if(disk[x*diameter + y] == 1)
				countOnes++;
		}
	}
	double * objxy = (double *)malloc(countOnes*2*sizeof(double));
	getneighbors(disk, countOnes, objxy, radius);
	
	//initial weights are all equal (1/Nparticles)
	double * weights = (double *)malloc(sizeof(double)*Nparticles);
	for(x = 0; x < Nparticles; x++){
		weights[x] = 1/((double)(Nparticles));
	}
	//initial likelihood to 0.0
	double * likelihood = (double *)malloc(sizeof(double)*Nparticles);
	double * arrayX = (double *)malloc(sizeof(double)*Nparticles);
	double * arrayY = (double *)malloc(sizeof(double)*Nparticles);
	double * xj = (double *)malloc(sizeof(double)*Nparticles);
	double * yj = (double *)malloc(sizeof(double)*Nparticles);
	double * CDF = (double *)malloc(sizeof(double)*Nparticles);
	double * u = (double *)malloc(sizeof(double)*Nparticles);
	int * ind = (int*)malloc(sizeof(int)*countOnes*Nparticles);
	for(x = 0; x < Nparticles; x++){
		arrayX[x] = xe;
		arrayY[x] = ye;
	}
	int k;
	
	int indX, indY;
	for(k = 1; k < Nfr; k++){
		//apply motion model
		//draws sample from motion model (random walk). The only prior information
		//is that the object moves 2x as fast as in the y direction
		for(x = 0; x < Nparticles; x++){
			arrayX[x] += 1 + 5*randn(seed, x);
			arrayY[x] += -2 + 2*randn(seed, x);
		}
		//particle filter likelihood
		for(x = 0; x < Nparticles; x++){
			//compute the likelihood: remember our assumption is that you know
			// foreground and the background image intensity distribution.
			// Notice that we consider here a likelihood ratio, instead of
			// p(z|x). It is possible in this case. why? a hometask for you.		
			//calc ind
			for(y = 0; y < countOnes; y++){
				indX = roundDouble(arrayX[x]) + objxy[y*2 + 1];
				indY = roundDouble(arrayY[x]) + objxy[y*2];
				ind[x*countOnes + y] = fabs(indX*IszY*Nfr + indY*Nfr + k);
				if(ind[x*countOnes + y] >= max_size)
					ind[x*countOnes + y] = 0;
			}
			likelihood[x] = 0;
			for(y = 0; y < countOnes; y++)
				likelihood[x] += (pow((I[ind[x*countOnes + y]] - 100),2) - pow((I[ind[x*countOnes + y]]-228),2))/50.0;
			likelihood[x] = likelihood[x]/((double) countOnes);
		}
		// update & normalize weights
		// using equation (63) of Arulampalam Tutorial
		for(x = 0; x < Nparticles; x++){
			weights[x] = weights[x] * exp(likelihood[x]);
		}
		double sumWeights = 0;
		for(x = 0; x < Nparticles; x++){
			sumWeights += weights[x];
		}
		for(x = 0; x < Nparticles; x++){
			weights[x] = weights[x]/sumWeights;
		}
		xe = 0;
		ye = 0;
		// estimate the object location by expected values
		for(x = 0; x < Nparticles; x++){
			xe += arrayX[x] * weights[x];
			ye += arrayY[x] * weights[x];
		}
		double distance = sqrt( pow((double)(xe-(int)roundDouble(IszY/2.0)),2) + pow((double)(ye-(int)roundDouble(IszX/2.0)),2) );
		
		printf("%lf\n", xe);
		printf("%lf\n", ye);
		printf("%lf\n", distance);

		CDF[0] = weights[0];
		for(x = 1; x < Nparticles; x++){
			CDF[x] = weights[x] + CDF[x-1];
		}
		double u1 = (1/((double)(Nparticles)))*randu(seed, 0);
		for(x = 0; x < Nparticles; x++){
			u[x] = u1 + x/((double)(Nparticles));
		}
		int j, i;
		
		for(j = 0; j < Nparticles; j++){
			i = findIndex(CDF, Nparticles, u[j]);
			if(i == -1)
				i = Nparticles-1;
			xj[j] = arrayX[i];
			yj[j] = arrayY[i];
			
		}
		
		for(x = 0; x < Nparticles; x++){
			//reassign arrayX and arrayY
			arrayX[x] = xj[x];
			arrayY[x] = yj[x];
			weights[x] = 1/((double)(Nparticles));
		}
	}
	free(disk);
	free(objxy);
	free(weights);
	free(likelihood);
	free(xj);
	free(yj);
	free(arrayX);
	free(arrayY);
	free(CDF);
	free(u);
	free(ind);

}

int main(int argc, char * argv[])
{
	int IszX, IszY, Nfr, Nparticles;
	IszX = DIMX;
	IszY = DIMY;
	Nparticles = NPART;
	Nfr = NFRAME;
	
	//establish seed
	int * __restrict__ seed = (int *)malloc(sizeof(int)*Nparticles);
	int i;
	for(i = 0; i < Nparticles; i++)
		seed[i] = i*i;
	//malloc matrix
	int * __restrict__ I = (int *)malloc(sizeof(int)*IszX*IszY*Nfr);
	
	//call video sequence
	videoSequence(I, IszX, IszY, Nfr, seed);
	
	//call particle filter
	particleFilter(I, IszX, IszY, Nfr, seed, Nparticles);
	
	free(seed);
	free(I);
	return 0;
}
