#include "binding.h"
#include "defines.h"
#include "stdio.h"
#include "time.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"

int main (int argc, char const *argv[])
{
	int i,n,icoef;

	int npoints;
	double coords[MAX_COORD],radius[MAX_ATOM],coef[20*MAX_ATOM];
	char buf[999];
	int check;
	
	if( argc < 2 ) {
		scanf("%s",buf);
		if( strncmp(buf,"NPOINTS",7) ) {
			printf("ERROR: expecting 'NPOINTS'\n");
			exit(-1);
		}
		n = scanf("%i",&npoints);
		if( 1 != n ) {
			printf("ERROR: reading npoints\n");
			exit(-1);
		}
		// printf("reading %i coords from stdin\n",npoints);
		n = scanf("%s",buf);
		if( 1 != n ) {
			printf("ERROR: reading 'COORDS'\n");
			exit(-1);
		}
		if( strncmp(buf,"COORDS",6) ) {
			printf("ERROR: expecting 'COORDS'\n");
			exit(-1);
		}
		double x,y,z,r;
		for(i=0;i<npoints;i++) {
			n = scanf("%lf %lf %lf %lf",&x,&y,&z,&r);
			if( 4 != n ) {
				printf("ERROR: reading x, y, z, r\n");
				exit(-1);
			}
			coords[3*i+0] = x;
			coords[3*i+1] = y;
			coords[3*i+2] = z;
			radius[i] = r;
			// printf("%i %i %f %f %f %f\n",i,n,x,y,z,r);
		}
		n = scanf("%s",buf);
		if( 1 != n ) {
			printf("ERROR: reading 'WEIGHTS'");
			exit(-1);
		}
		if( strncmp(buf,"WEIGHTS",7) ) {
			printf("ERROR: expecting 'WEIGHTS'\n");
			exit(-1);
		}
		for(i=0;i<npoints;i++) {
			for(icoef=0; icoef<20; icoef++) {
				double tmp;
				n = scanf("%lf",&tmp);
				if( 1 != n ) {
					printf("ERROR: reading coef %i %i\n",i,icoef);
					exit(-1);
				}
				coef[icoef*npoints+i] = tmp;
			}
		}
		n = scanf("%s",buf);
		if( 1 != n ) {
			printf("ERROR: reading npoints");
			exit(-1);
		}
		if( strncmp(buf,"END",3) ) {
			printf("ERROR: expecting 'END'\n");
			exit(-1);
		}
		// printf("points: \n",npoints);
		// for(i=0;i<npoints;i++) {
		// 	printf("%i p = %lf %lf %lf; r = %lfm\n",i,coords[3*i+0],coords[3*i+1],coords[3*i+2],radius[i]);
		// }
	} else {
		// printf("reading coords from %s\n",argv[1]);

		int j,npoint;
		double x,y,z,r;
		char line[MAX_LENGTH];
		char* val;
		FILE* fp;
	
		if ( (fp = fopen( argv[1] ,"r")) == NULL ) return;
	
		val=fgets(line, MAX_LENGTH, fp);
		sscanf(line,"%d",&npoint);
		val=fgets(line, MAX_LENGTH, fp);
		npoints = npoint;
	
		for (i=0;i<npoint;i++)
		{
			val=fgets(line, MAX_LENGTH, fp);
			sscanf(line,"%d %lf %lf %lf %lf",&j,&x,&y,&z,&r);
	//		printf("%d %d %f %f %f %f\n",i,j,x,y,z,r);
			coords[3*i+0] = x;
			coords[3*i+1] = y;
			coords[3*i+2] = z;
			radius[i] = r;
		}
	
		fclose(fp);
		// printf("Data read properly !!\n");
	}

	double *wt = malloc(20*npoints*sizeof(double));
	for( i = 0; i < 20*npoints; i++ ) wt[i] = 1.0;
	double *surf_out  = malloc((npoints+4)*sizeof(double));
	double *deriv_out = malloc((npoints+4)*sizeof(double));
	get_alpha20_surf_weighted_deriv( npoints, coords, radius, wt, surf_out, deriv_out );
	
	return 0;
}