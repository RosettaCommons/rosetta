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
	// This causes a stack overflow if you go past a certan size
	// double coords[MAX_COORD],radius[MAX_ATOM],coef[20*MAX_ATOM];
	double * coords = malloc((MAX_COORD)*sizeof(double));
	double * radius = malloc((MAX_ATOM)*sizeof(double));
	double * coef = malloc((20*MAX_ATOM)*sizeof(double));
	char buf[999];
	int check;
	int READ_WEIGHTS = 0;
	
	if( argc < 2 ) {
		printf("at least one arg is required! surf, surf_vol, surf_vol_deriv, ...\n");
		exit(-1);
	}
	
	if( !strncmp(argv[1],"alpha20_deriv",13) ) READ_WEIGHTS = 1;
	
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
	
	if( READ_WEIGHTS ) {
	
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
	
	}
	
	n = scanf("%s",buf);
	if( 1 != n ) {
		printf("ERROR: reading END");
		exit(-1);
	}
	if( strncmp(buf,"END",3) ) {
		printf("ERROR: expecting 'END'\n");
		exit(-1);
	}
	
	if(        !strncmp(argv[1],"alpha20_surf",99) ) {
		get_alpha20_surf( npoints, coords, radius );
	} else if( !strncmp(argv[1],"alpha20_deriv_surf",99) ) {
		get_alpha20_surf_weighted_deriv( npoints, coords, radius, coef );
	} else if( !strncmp(argv[1],"surf_vol",99) ) {
		double * surf_out = malloc((npoints+4)*sizeof(double));
		double *  vol_out = malloc((npoints+4)*sizeof(double));		
		get_surf_vol( npoints, coords, radius,  surf_out, vol_out );
		free(surf_out);
		free(vol_out);
	} else if( !strncmp(argv[1],"surf_vol_deriv",99) ) {
		double * surf_out = malloc((npoints+4)*sizeof(double));
		double *  vol_out = malloc((npoints+4)*sizeof(double));		
		double *dsurf_out = malloc((npoints+4)*3*sizeof(double));
		double * dvol_out = malloc((npoints+4)*3*sizeof(double));		
		get_surf_vol_deriv( npoints, coords, radius,  surf_out, vol_out, dsurf_out, dvol_out );
		free( surf_out);
		free(  vol_out);
		free(dsurf_out);
		free( dvol_out);
	} else if( !strncmp(argv[1],"cavballs",99) ) {
		get_cavballs( npoints, coords, radius );
	}
		
	free( coords );
	free( radius );
	free( coef );
	return 0;
}
