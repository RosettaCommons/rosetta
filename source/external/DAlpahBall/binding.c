#include "defines.h"
#include "stdio.h"
#include "time.h"
#include "string.h"
#include "stdlib.h"
#include "math.h"

	#define F77Name1(x) x##_   /* Need to add one underscore to Fortran program */

	#if defined intel
	#define F77Name2(x) x##_   /* Need to add one underscore to Fortran program for 
					Intel compiler */
	#else
	#define F77Name2(x) x##__   /* Need to add two underscores to Fortran program for GNU 
					f77 compiler */
	#endif

	extern struct{
		int npoints;
		int nvertex;
		char redinfo[MAX_POINT];
	} F77Name2(vertex_zone);

	extern struct{
		double coord[MAX_COORD];
		double radius[MAX_POINT];
		double coord4[MAX_POINT];
	} F77Name2(xyz_vertex);

	extern struct{
		int   depth[MAX_TETRA];
		char isvoid[MAX_TETRA];
	} F77Name2(voids);

	extern void clear_all_gmp_arrays(int *nvertex);
	extern void clear_alf_gmp();
	extern void F77Name2(adjust_nsphere)();
	extern void F77Name1(setup)();
	extern void F77Name1(regular3)(int *nredundant, int list_redundant[MAX_ATOM]);
	extern void F77Name2(remove_inf)();
	extern void F77Name1(alfcx)(double *alpha_val, int *nred, int list_redalpha[MAX_ATOM]);
	extern void F77Name2(readjust_nsphere)(int *nat, int *nred, int list_redalpha[MAX_ATOM]);
	extern void F77Name2(write_del)(char delaunay_file[FSIZE]);
	extern void F77Name2(write_alpha)(char alpha_file[FSIZE]);
	extern void F77Name2(surface_only)(double coef[MAX_ATOM], double *wsurf, double *surf, double surface [MAX_ATOM]);
	extern void F77Name2(surface_only_void)(double coef[MAX_ATOM], double *wsurf, double *surf, double surface [MAX_ATOM] );
	extern void F77Name2(surface_deriv_coord)(double coef[MAX_ATOM], double *wsurf, double *surf, double surface [MAX_ATOM], double deriv_surf [3*MAX_ATOM]);
	extern void F77Name2(volume_only)       ( double coef[MAX_ATOM], double *wsurf, double *wvol, double *surf, double *vol, double ballwsurf[MAX_ATOM], double ballwvol[MAX_ATOM] );

	extern void F77Name2(volume_deriv_coord)( double coef[MAX_ATOM], double *wsurf, double *wvol, 
			double *surf, double *vol, double ballwsurf[MAX_ATOM], double ballwvol[MAX_ATOM], 
			double dwsurf[3*MAX_ATOM], double dwvol[3*MAX_ATOM] );

	extern void F77Name2(cavballs)( double *x, double *y, double *z, double *r, int *ncav );

void get_surf_vol( int npoints, double *coords, double *radius, double *surf_out, double *vol_out ) {
	int i;

	double tmp_surf,tmp_wsurf,tmp_vol,tmp_wvol, alpha_val=0.0;
	int tmp_nredundant,tmp_nredalpha,tmp_nat_orig;
	double * tmp_coef = malloc((MAX_COORD)*sizeof(double));
	int * tmp_list_redalpha = malloc((MAX_COORD)*sizeof(int));
	int * tmp_list_redundant = malloc((MAX_COORD)*sizeof(int));
	
	for( i = 0; i<MAX_ATOM; i++ ) tmp_coef[i] = 1.0;
	
	F77Name2(vertex_zone).npoints = npoints;
	for( i=0; i < npoints; i++ ) {
		F77Name2(xyz_vertex).coord[3*i+0] = coords[3*i+0];
		F77Name2(xyz_vertex).coord[3*i+1] = coords[3*i+1];
		F77Name2(xyz_vertex).coord[3*i+2] = coords[3*i+2];
		F77Name2(xyz_vertex).radius[i] = radius[i];
	}		
		
	F77Name2(adjust_nsphere)();
	F77Name1(setup)();
	tmp_nat_orig = npoints;

	F77Name1(regular3)(&tmp_nredundant,tmp_list_redundant);	/* computes triangulation; stores tetrahedra*/
	
	F77Name1(alfcx)(&alpha_val, &tmp_nredalpha, tmp_list_redalpha); /* computes alpha shape, alpha=0 */

	F77Name2(readjust_nsphere)(&tmp_nat_orig, &tmp_nredalpha, tmp_list_redalpha);

	F77Name2(volume_only)(tmp_coef, &tmp_wsurf, &tmp_wvol, &tmp_surf, &tmp_vol, surf_out, vol_out );
	
	for(i=0;i<npoints;i++) {
		printf("%d %.20lf %.20lf\n",i+1,surf_out[i+4],vol_out[i+4]);
	}

	free( tmp_coef );
	free( tmp_list_redalpha );
	free( tmp_list_redundant );
	
	return;
}


void get_surf_vol_deriv( int npoints, double *coords, double *radius, 
              double *surf_out, double *vol_out, double *dsurf_out, double *dvol_out  ) 
{
	int i;
	double tmp_coef[MAX_ATOM],tmp_surf,tmp_wsurf,tmp_vol,tmp_wvol, alpha_val=0.0;
	int tmp_nredundant,tmp_nredalpha,tmp_nat_orig;
	int * tmp_list_redalpha = malloc((MAX_COORD)*sizeof(int));
	int * tmp_list_redundant = malloc((MAX_COORD)*sizeof(int));
	
	for( i = 0; i<MAX_ATOM; i++ ) tmp_coef[i] = 1.0;
	
	F77Name2(vertex_zone).npoints = npoints;
	for( i=0; i < npoints; i++ ) {
		F77Name2(xyz_vertex).coord[3*i+0] = coords[3*i+0];
		F77Name2(xyz_vertex).coord[3*i+1] = coords[3*i+1];
		F77Name2(xyz_vertex).coord[3*i+2] = coords[3*i+2];
		F77Name2(xyz_vertex).radius[i] = radius[i];
	}		
		
	F77Name2(adjust_nsphere)();
	F77Name1(setup)();
	tmp_nat_orig = npoints;

	F77Name1(regular3)(&tmp_nredundant,tmp_list_redundant);	/* computes triangulation; stores tetrahedra*/
	
	F77Name1(alfcx)(&alpha_val, &tmp_nredalpha, tmp_list_redalpha); /* computes alpha shape, alpha=0 */

	F77Name2(readjust_nsphere)(&tmp_nat_orig, &tmp_nredalpha, tmp_list_redalpha);

	F77Name2(volume_deriv_coord)(tmp_coef, &tmp_wsurf, &tmp_wvol, &tmp_surf, &tmp_vol, surf_out, vol_out, dsurf_out, dvol_out );
	
	for(i=0;i<npoints;i++) {
		printf("%d    %.20lf    %.20lf    %.20lf %.20lf %.20lf    %.20lf %.20lf %.20lf\n",
		   i+1,surf_out[i+4],vol_out[i+4], 
			dsurf_out[3*(i+4)+0],dsurf_out[3*(i+4)+1],dsurf_out[3*(i+4)+2],
   		 dvol_out[3*(i+4)+0], dvol_out[3*(i+4)+1], dvol_out[3*(i+4)+2] );
	}
	

	free( tmp_list_redalpha );
	free( tmp_list_redundant );
	
}

double alphas[20] = { 0.0000000,0.6167577,0.8835883,1.0960643,
                      1.2816093,1.4507516,1.6089340,1.7595542,
							 1.9051524,2.0480757,2.1910429,2.3378414,
 							 2.4943877,2.6704755,2.8827284,3.1596130,
							 3.5499700,4.1377324,5.0680891,6.5954530 };

double get_alpha20_surf( 
	int npoints, double *coords, double *radius
) {
	int i;
	int ialpha;
	double alpha_val=0;
	double coordat[3], radiusat, rad0, rad, x,y,z;

	double tmp_surf,tmp_wsurf;
	double void_surf,void_wsurf;
	int tmp_nredundant;
	int tmp_nredalpha;
	int tmp_nat_orig;	
	double * tmp_coef = malloc((MAX_COORD)*sizeof(double));
	int * tmp_list_redalpha = malloc((MAX_COORD)*sizeof(int));
	int * tmp_list_redundant = malloc((MAX_COORD)*sizeof(int));

	double *surf_out = malloc((npoints+4)*sizeof(double));

	// printf("\n\n    397   \n\n\n");
		
	F77Name2(vertex_zone).npoints = npoints;
	for( i=0; i < npoints; i++ ) {
		F77Name2(xyz_vertex).coord[3*i+0] = coords[3*i+0];
		F77Name2(xyz_vertex).coord[3*i+1] = coords[3*i+1];
		F77Name2(xyz_vertex).coord[3*i+2] = coords[3*i+2];
		F77Name2(xyz_vertex).radius[i] = radius[i];
	}		
		
	F77Name2(adjust_nsphere)();
	F77Name1(setup)();
	tmp_nat_orig = npoints;

	F77Name1(regular3)(&tmp_nredundant,tmp_list_redundant);	/* computes triangulation; stores tetrahedra*/

	for( ialpha = 0; ialpha < 20; ialpha++ ) {
	
		for (i=0;i<npoints;i++) tmp_coef[i+4] = 1.0;

		alpha_val = alphas[ialpha];// 6.5954530;
		for( i=0;i < npoints;i++ ) {
			rad0 = radius[i];
			rad  = sqrt( alpha_val*alpha_val + rad0*rad0 );
			F77Name2(xyz_vertex).radius[i+4] = rad;
			x = coords[3*i+0];
			y = coords[3*i+1];
			z = coords[3*i+2];
			F77Name2(xyz_vertex).coord4[i+4] = x*x + y*y + z*z - rad*rad;
		}
		alpha_val = 0;
	
		F77Name1(alfcx)(&alpha_val, &tmp_nredalpha, tmp_list_redalpha); /* computes alpha shape, alpha=0 */

		F77Name2(readjust_nsphere)(&tmp_nat_orig, &tmp_nredalpha, tmp_list_redalpha);

		F77Name2(surface_only)(tmp_coef, &tmp_wsurf , &tmp_surf , surf_out );
		for(i=0;i<npoints;i++) {
			if( surf_out[i+4] < 0.000001 && surf_out[i+4] > -0.000001 ) surf_out[i+4] = 0.0;
			printf("%d %d %.20lf\n",ialpha+1,i+1,surf_out[i+4]);
		}

	}

	free( tmp_coef );
	free( tmp_list_redalpha );
	free( tmp_list_redundant );
	
	free(surf_out);
	return( tmp_wsurf );

}


double get_alpha20_surf_weighted_deriv( 
	int npoints, double *coords, double *radius, double *weights
) {
	int i;
	int ialpha;
	double alpha_val=0;
	double coordat[3], radiusat, rad0, rad, x,y,z;

	double tmp_surf,tmp_wsurf;
	// double tmp_deriv_surf[3*MAX_ATOM];
	double void_surf,void_wsurf;
	// double void_deriv_surf[3*MAX_ATOM];
	int tmp_nredundant;
	int tmp_nredalpha;
	int tmp_nat_orig;	
	double * tmp_coef = malloc((MAX_COORD)*sizeof(double));
	int * tmp_list_redalpha = malloc((MAX_COORD)*sizeof(int));
	int * tmp_list_redundant = malloc((MAX_COORD)*sizeof(int));

	double * surf_out = malloc((npoints+4)  *sizeof(double));
	double *deriv_out = malloc((npoints+4)*3*sizeof(double));
		
	F77Name2(vertex_zone).npoints = npoints;
	for( i=0; i < npoints; i++ ) {
		F77Name2(xyz_vertex).coord[3*i+0] = coords[3*i+0];
		F77Name2(xyz_vertex).coord[3*i+1] = coords[3*i+1];
		F77Name2(xyz_vertex).coord[3*i+2] = coords[3*i+2];
		F77Name2(xyz_vertex).radius[i] = radius[i];
	}		
   
   double *deriv_x   = (double*)malloc(sizeof(double)*npoints);
   double *deriv_y   = (double*)malloc(sizeof(double)*npoints);
   double *deriv_z   = (double*)malloc(sizeof(double)*npoints);
   double *tot_score = (double*)malloc(sizeof(double)*npoints);
   for( i=0; i < npoints; i++ ) {
      tot_score[i] = 0.0;
      deriv_x[i]   = 0.0;
      deriv_y[i]   = 0.0;
      deriv_z[i]   = 0.0;
   }
		
	F77Name2(adjust_nsphere)();
	F77Name1(setup)();
	tmp_nat_orig = npoints;

	// for( i = 0; i < 16; i++ ) {
	// 	printf("%f %f %f %f\n",coords[i], F77Name2(xyz_vertex).coord[i], radius[i], F77Name2(xyz_vertex).radius[i] );
	// }		

	F77Name1(regular3)(&tmp_nredundant,tmp_list_redundant);	/* computes triangulation; stores tetrahedra*/
	
	for( ialpha = 0; ialpha < 20; ialpha++ ) {
	
		for (i=0;i<npoints;i++) tmp_coef[i+4] = weights[ialpha*npoints+i];

		alpha_val = alphas[ialpha];// 6.5954530;
		for( i=0;i < npoints;i++ ) {
			rad0 = radius[i];
			rad  = sqrt( alpha_val*alpha_val + rad0*rad0 );
			F77Name2(xyz_vertex).radius[i+4] = rad;
			x = coords[3*i+0];
			y = coords[3*i+1];
			z = coords[3*i+2];
			F77Name2(xyz_vertex).coord4[i+4] = x*x + y*y + z*z - rad*rad;
		}
		alpha_val = 0;
	
		F77Name1(alfcx)(&alpha_val, &tmp_nredalpha, tmp_list_redalpha); /* computes alpha shape, alpha=0 */
		// F77Name1(markvoids)();
		F77Name2(readjust_nsphere)(&tmp_nat_orig, &tmp_nredalpha, tmp_list_redalpha);

		F77Name2(surface_deriv_coord)(tmp_coef, &tmp_wsurf , &tmp_surf , surf_out, deriv_out );
		for(i=0;i<npoints;i++) {
         tot_score[i] += surf_out[i+4];
         deriv_x[i] += deriv_out[3*i+12];
         deriv_y[i] += deriv_out[3*i+13];
         deriv_z[i] += deriv_out[3*i+14];
		}

	}

	for(i=0;i<npoints;i++) {
	   printf("%i %.20lf %.20lf %.20lf %.20lf\n",i+1,tot_score[i],deriv_x[i],deriv_y[i],deriv_z[i]);
   }

	free( tmp_coef );
	free( tmp_list_redalpha );
	free( tmp_list_redundant );
   
	free( surf_out);
	free(deriv_out);
	return( tmp_wsurf );

}


void get_cavballs( int npoints, double *coords, double *radius ) {
	int i;
	double tmp_surf,tmp_wsurf,tmp_vol,tmp_wvol, alpha_val=0.0;
	int tmp_nredundant,tmp_nredalpha,tmp_nat_orig;
	double * tmp_coef = malloc((MAX_COORD)*sizeof(double));
	int * tmp_list_redalpha = malloc((MAX_COORD)*sizeof(int));
	int * tmp_list_redundant = malloc((MAX_COORD)*sizeof(int));
	
	for( i = 0; i<MAX_ATOM; i++ ) tmp_coef[i] = 1.0;
	
	F77Name2(vertex_zone).npoints = npoints;
	for( i=0; i < npoints; i++ ) {
		F77Name2(xyz_vertex).coord[3*i+0] = coords[3*i+0];
		F77Name2(xyz_vertex).coord[3*i+1] = coords[3*i+1];
		F77Name2(xyz_vertex).coord[3*i+2] = coords[3*i+2];
		F77Name2(xyz_vertex).radius[i] = radius[i];
	}		
		
	F77Name2(adjust_nsphere)();
	F77Name1(setup)();
	tmp_nat_orig = npoints;

	F77Name1(regular3)(&tmp_nredundant,tmp_list_redundant);	/* computes triangulation; stores tetrahedra*/
	
	F77Name1(alfcx)(&alpha_val, &tmp_nredalpha, tmp_list_redalpha); /* computes alpha shape, alpha=0 */

	double *x = malloc((MAX_TETRA)*sizeof(double));
	double *y = malloc((MAX_TETRA)*sizeof(double));
	double *z = malloc((MAX_TETRA)*sizeof(double));
	double *r = malloc((MAX_TETRA)*sizeof(double));
	int ncav = 0;
	cavballs_( x, y, z, r, &ncav );
	
	printf("%d\n",ncav);
	for( i = 0; i<ncav; i++ ) {
		printf("%i %f %f %f %f\n",i+1,x[i],y[i],z[i],r[i]);
	}
	
	free( tmp_coef );
	free( tmp_list_redalpha );
	free( tmp_list_redundant );
}


