
double get_alpha20_surf               ( int npoints, double *coords, double *radius );
double get_alpha20_surf_weighted_deriv( int npoints, double *coords, double *radius, double *weights );

void get_surf_vol( int npoints, double *coords, double *radius, double *surf_out, double *vol_out );

void get_surf_vol_deriv( int npoints, double *coords, double *radius, 
	double *surf_out, double *vol_out, double *dsurf_out, double *dvol_out  );

	
void get_cavballs( int npoints, double *coords, double *radius );

