#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>


#define MAX_ARGS 10
#define REC_LENGTH 49	// size of a record in db
#define REC_WINDOW 10	// number of records to read at a time
#define LATITUDE_POS 28	// location of latitude coordinates in input record
#define OPEN 10000	// initial value of nearest neighbors


#define NUM_VALUES 444444 // Aroung 3MB of input data is generated
#define LATITUDE 30
#define LONGITUDE 90



#ifdef SMALL_DATASET
#define NUMNEAREST 10
#else 
#define NUMNEAREST 50
#endif


int main(int argc, char* argv[]) 
{
	int    i=0,j=0, k=0, rec_count=0, done=0;
	float target_lat, target_long, tmp_lat=0, tmp_long=0;
	char   sandbox[REC_LENGTH * REC_WINDOW * 2];


	srand(10);
	

	k = NUMNEAREST;
	target_lat = LATITUDE;
	target_long = LONGITUDE;

	char * __restrict__ entry  = malloc(k*REC_LENGTH* sizeof(char));
	if(entry == NULL) {
		fprintf(stderr, "no room for neighbor:entry \n");
		exit(0);
	}

	double * __restrict__ dist = malloc(k*sizeof(double));
	if(dist == NULL) {
		fprintf(stderr, "no room for neighbor:dist\n");
		exit(0);
	}

	//Initialize list of nearest neighbors to very large dist
	for( j = 0 ; j < k ; j++ ) {
		dist[j] = OPEN;
	}

	float *z;
	z  = (float *) malloc(REC_WINDOW * sizeof(float));
	

	int year,month,date,hour,num,speed,press;
	float lat,lon;
	int hours[4] = {0,6,12,18};
	char *name;
	char names[21][10] = {"ALBERTO", "BERYL", "CHRIS","DEBBY","ERNESTO","FLORENCE","GORDON",
	    "HELENE","ISAAC","JOYCE","KIRK","LESLIE","MICHAEL","NADINE","OSCAR","PATTY","RAFAEL",
	    "SANDY","TONY","VALERIE","WILLIAM"};

	done=0; 
	char src[REC_LENGTH+1];
	while(done<=NUM_VALUES) {

		rec_count = 0;
		for(i=0;i<REC_WINDOW && done <= NUM_VALUES; i++)
	    {
	    	rec_count +=1;
			done+=1;
			year = 1950 + rand() % 55;
			month = 1 + rand() % 12;
			date = 1 + rand() % 28;
			hour = hours[rand()%4];
			num = 1 + rand() % 28;
			name = names[rand()%21];
			lat = ((float)(7 + rand() % 63)) + ((float) rand() / (float) 0x7fffffff);
			lon = ((float)(rand() % 358)) + ((float) rand() / (float) 0x7fffffff); 
			speed = 10+ rand() % 155;
			press = rand() % 900;

			sprintf(src, "%4d %2d %2d %2d %2d %-9s %5.1f %5.1f %4d %4d\n", year, month, date, hour, num, name, lat, lon, speed, press);
			if(i==0)
			{
				strcpy(sandbox, src);
			}
			else
			{
				strcat(sandbox, src);
			}
		}
		


        for (i = 0; i < rec_count; i++){
            float tmp_lat = atof(sandbox+(i * REC_LENGTH + LATITUDE_POS - 1));
            float tmp_long = atof(sandbox+(i * REC_LENGTH + LATITUDE_POS - 1)+5);
			z[i] = sqrt(( (tmp_lat-target_lat) * (tmp_lat-target_lat) )+( (tmp_long-target_long) * (tmp_long-target_long) ));
        }

		
        for( i = 0 ; i < rec_count ; i++ ) {
			float max_dist = -1;
			int max_idx = 0;
			// find a neighbor with greatest dist and take his spot if allowed!
			for( j = 0 ; j < k ; j++ ) {
				if( dist[j] > max_dist ) {
					max_dist = dist[j];
					max_idx = j;
				}
			}
			// compare each record with max value to find the nearest neighbor
			if( z[i] < dist[max_idx] ) {
				sandbox[(i+1)*REC_LENGTH-1] = '\0';
			  	strcpy((entry + max_idx*REC_LENGTH), (sandbox+i*REC_LENGTH));
			  	dist[max_idx] = z[i];
			}
		}
	}//End while loop
	
	for( j = 0 ; j < k ; j++ ) {
		if( !(dist[j] == OPEN))
			printf("%s\t%f\n", (entry + j*REC_LENGTH), dist[j]);
	}
}

