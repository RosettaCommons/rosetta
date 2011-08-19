/*****************************************************************************
 * DAlphaVol.c  : Computes the volume and surface area of a union of balls
 *             using the Alpha Shape theory introduced by Herbert
 *             Edelsbrunner.
 *             Also computes the derivatives of the volume and surface
 *             areas with respect to the coordinates of the centers of
 *	       the balls, if needed.
 *
 * Information on the Alpha Shape method, and on the calculationof the derivatives 
 * can be found in:
 *             
 *		1. H. Edelsbrunner, "The Union of balls and its dual shape", Discrete 
 *		   Comput. Geom., 13:415-440 (1995)
 *
 *		2. H. Edelsbrunner and E.P. Mucke, "Three-dimensional alpha shapes",
 *                 ACM Trans. Graphics, 13:43-72 (1994)
 *
 *              3. J. Liang, H. Edelsbrunner, P. Fu, P.V. Sudhakar and S. Subramaniam,
 *                 "Analytical shape computation of macromolecules: I. Molecular area",
 *                 Proteins: Struct. Func. Genet., 33:1-17 (1998)
 *
 *		4. H. Edelsbrunner and P. Koehl, "The weighted volume derivative of a
 *		   space filling diagram", Proc. Natl. Acad. Sci. (USA), 100:2203-2208,
 *		   (2003)
 *
 *		5. R. Bryant, H. Edelsbrunner, P. Koehl and M. Levitt, "The area
 *		   derivative of a space-filling diagram", Discrete Comput. Geom.,
 *                 32:293-308 (2004)
 *
 *		6. H. Edelsbrunner and P. Koehl, "The geometry of biomolecular solvation",
 *		   Combinatorial and Discrete Geometry (MSRI Series), to appear (2005)
 *
 * Calculations in the program DAlphaVol are performed first in floating point; if the
 * results are inconclusive, the programs switches to arbitrary precision arithmetics,
 * using the package GMP. 
 *
 * In cases of degeneracies (i.e. a geometric predicate has a value of 0), the
 * program applies the procedure defined as "Simulation of Simplicity"
 *
 *		1. H. Edelsbrunner and E.P. Mucke, "Simulation of Simplicity: a technique
 *		   to cope with degenerate cases in geometric algorithms", ACM trans.
 *		   Graphics, 9:66-104 (1990)
 *
 * For historical reasons, the program is written as a mixture of Fortran and C.
 *
 * The whole suite of program is distributed under the GNU Lesser general Public License.
 *
 * Copyright (C) 2007 Patrice Koehl
 *
 */

#include "defines.h"
#include "stdio.h"
#include "time.h"
#include "string.h"
#include "stdlib.h"

/*****************************************************************************************
 * Naming convention between C and Fortran 
 *
 * Let us consider two Fortran subroutines : foo, and foo_with_underscore
 *
 * Case 1: the name of the Fortran subroutine (foo) does not contain an underscore (_)
 *	   In the C program, we must add ONE underscore to the Fortran name:
 *	   the C program refers to the subroutine as:
 *		foo_(...)
 *	   while the Fortran code writes it as:
 *		foo(...)
 *	   This is independent of the compiler pair (at least for gcc/f77 and 
 *	   Intel icc/ifort)
 *
 * Case 2: the name of the Fortran subroutine (foo_with_underscore) contains at least one 
 * underscore (_)
 *	   Treatment of case 2 is compiler dependent!
 *	   - The Intel compiler treats this case as if it was case 1, i.e. requires 
 *	     ONE underscore at the end of the Fortran name in the C program
 *	   - the gnu f77 however requires TWO underscores at the end of the Fortran name 
 *	     in the C program.
 *
 * I solve this by introducing two functions, FTName1 and FTName2, where FTName2 is 
 * compiler dependent
 */

#define F77Name1(x) x##_   /* Need to add one underscore to Fortran program */

#if defined intel
#define F77Name2(x) x##_   /* Need to add one underscore to Fortran program for 
				Intel compiler */
#else
#define F77Name2(x) x##__   /* Need to add two underscores to Fortran program for GNU 
				f77 compiler */
#endif

/**************************************************************************************
* local common variables
*/

static char input_file[FSIZE], output_file[FSIZE], alpha_file[FSIZE];
static char delaunay_file[FSIZE];

static int const YES = 1;
static int const NO = 0;
static int const ERROR = 0;
static int const SAFE = 1;
static int const ONE = 1;

static int AlphaFileName = 0;
static int DelaunayFileName = 0;

static int npoint;
static int nat_orig;
static int nredundant;
static int nredalpha;
static double vol,surf,wvol,wsurf;

static int list_redundant[MAX_ATOM];
static int list_redalpha[MAX_ATOM];

/* Variables for volume and surface calculation:
	coef:		: we compute a weighted sum of the individual surface 
			  and volume; coef is the array containing the weights
	volume, surface	: contribution of each atom to the total surface and volume of the molecule
	deriv_vol
	deriv_surf	: derivatives of surface and volume
								*/
static double coef[MAX_ATOM];
static double volume[MAX_ATOM], surface[MAX_ATOM];
static double deriv_vol[3*MAX_ATOM], deriv_surf[3*MAX_ATOM];
static double deriv_vol_dist[MAX_PAIR], deriv_surf_dist[MAX_PAIR];


/**************************************************************************************/
/* program parameters */

static char usage[] = "\n\
usage:    \n\
	%s DATA  -s <switch> -o <res_file> \n\
with:     \n\
	DATA ............... path name of input data file \n\
	-s <switch> ........ switch = 0: computes surface area and volume \n\
				      1: computes surface area, volume and derivatives wrt coord  \n\
				      2: computes surface area, volume and derivatives wrt internal distances  \n\
                                      3: compute Delaunay only                  \n\
                                      4: compute Delaunay + Dual Complex only   \n\
	-o <res_file> ...... path name of output result file (default DATA.vol) \n\
        -d <delaunay_file>.. path name of file for simplices of Delaunay (optional) \n\
        -a <alpha_file>..... path name of file for simplices of dual complex (optional) \n\
";

/**************************************************************************************/
/* Common structures between Fortran and C                                             */

/* Common blocks in Fortran are equivalent to data structures in C:

if in Fortran we have:

	real*8	a,b,c
	integer i,j,k

	common /values/ a,b,c,i,j,k

We use in C:

	extern struct {
		double a,b,c;
		int i,j,k;
	} values_

Some remarks:

	- the name of the common black (here "values") should be the same as the name of
	  the struct in C, with the addition of one (or two!) underscore

	- the underscore in the name of the structure follows the same rules as for
	  procedure (see above the definition of F77Name1 and F77Name2)

	- use "extern" if the common block / structure is first referenced in the Fortran
	  program; if the first reference is in the C program, remove extern

	- put variables in the same order in the structure and in the common block!!

For the Delaunay triangulation and the alpha complex, there are a series of variables that can be used indifferently in the C code and the Fortran code:

1) Variables that describe each tetrahedron of the (weighted) Delaunay triangulation:

	ntetra		: number of tetrehedrons
	tetra		: indices of the 4 vertices defining the tetrahedron
	tetra_neighbour	: the four neighbours of the tetrahedron. If the tetrahedron is (a,b,c,d),
			  the four neighbors are given in that order:
				tetraneighbor[i][0] is the neighbor sharing face b,c,d
				tetraneighbor[i][1] is the neighbor sharing face a,c,d
				tetraneighbor[i][2] is the neighbor sharing face a,b,d
				tetraneighbor[i][3] is the neighbor sharing face a,b,c
	tetra_nindex	: for each neighbor of the tetrahedron:
			  the neighbor of a,b,c,d sharing face b,c,d is a tetrahedron A,B,C,D
			  that in fact corresponds to a permutation of (b,c,d,e) where e is the
			  fourth vertex. tetra_nindex[][0] gives the position of e in (A,B,C,D)
			  (1, 2 3 or 4)
	tetra_link	: indices (in triangle list), of the 4 triangles (b,c,d), (a,c,d),
			  (a,b,d) and (a,b,c)
	tetra_status	: flag for each tetrahedron:
				1 if the tetrahedron is active
				0 if the tetrahedron is "dead" (i.e. does not belong to Delaunay
								triangulation)
	tetra_orient	: flag for orientation of the tetrahedron:
				1 if (a,b,c,d) in ccw order
				-1 otherwise
	tetra_hull	: flag for position of the tetrahedron:
				1 if (at least) one face is on the convex hull
				0 if tetrehedron is interior
	order		: Order of the tetrahedrons, to fit with a first depth search

2) Variables that describe each triangle of the (weighted) Delaunay triangulation:

	ntrig		: number of triangles
	trig		: indices of the 3 vertices defining the triangle
	trig_alp	: flag for each triangle:
				1 if the triangle is active
				0 if the triangle is "dead" (i.e. does not belong to Alpha complex)
	trig_link	: indices (in edge list), of the 3 edges (b,c), (a,c),
			  and (b,d)
	trig_hull	: flag for position of the triangle:
				1 if triangle is on the convex hull
				0 if triangle is interior
	trig_coef	: flag for each triangle (used for volume calculation only)
	
3) Variables that describe the edges of the (weighted) Delaunay triangulation:

	nedge		: number of edges
	edge		: indices of the 2 vertices defining the edge
	edge_hull	: flag for position of the edge:
				1 if both vertices are on convex hull
				0 otherwise
	edge_alp	: flag for each edge:
				1 if the edge is active
				0 if the edge is "dead" (i.e. does not belong to Alpha complex)

4) Variables that link vertices to tetrahedra

	vertex_nlink	: number of tetrahedra attached to each vertex
	vertex_link	: list of tetrahedra attached to each vertex
												*/
/* Now list all Fortran common blocks / C data structures that lies at the interface */

extern struct{
	int ntetra;
	int tetra[MAX_TETRA][4];
	int tetra_neighbor[MAX_TETRA][4];
	} F77Name2(tetra_zone);

extern struct{
	char tetra_info[MAX_TETRA];
	char tetra_nindex[MAX_TETRA];
	} F77Name2(tetra_stat);

extern struct{
	char tetra_edge[MAX_TETRA];
	} F77Name2(alp_zone);

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

/**************************************************************************************/
/* local procedures with common variables */

int command_line(char *argv[], int argc);
int read_data();
int alpha_setup();
int alpha_geom(int switchval);
void get_volume();
void get_deriv();
void print_volume(double vol, double surf);
void write_delaunay();
void write_alphacx();
extern void clear_all_gmp_arrays(int *nvertex);
extern void clear_alf_gmp();
extern void F77Name2(adjust_nsphere)();
extern void F77Name1(setup)();
extern void F77Name1(regular3)(int *nredundant, int list_redundant[MAX_ATOM]);
extern void F77Name2(remove_inf)();
extern void F77Name1(alfcx)(double *alpha_val, int *nred, int list_redalpha[MAX_ATOM]);
extern void F77Name2(readjust_nsphere)(int *nat, int *nred, int list_redalpha[MAX_ATOM]);

extern void F77Name2(volume_only)(double coef[MAX_ATOM], double *wsurf, double *wvol, double *surf, double *vol, 
	double surface [MAX_ATOM], double volume [MAX_ATOM]);
extern void F77Name2(volume_deriv_coord)(double coef[MAX_ATOM], double *wsurf, double *wvol, double *surf, double *vol, 
	double surface [MAX_ATOM], double volume [MAX_ATOM], double deriv_surf [3*MAX_ATOM], double deriv_vol[3*MAX_ATOM]);
extern void F77Name2(volume_deriv_dist)(double coef[MAX_ATOM], double *wsurf, double *wvol, double *surf, double *vol, 
	double surface [MAX_ATOM], double volume [MAX_ATOM], double deriv_surf [3*MAX_ATOM], double deriv_vol[3*MAX_ATOM],
	double deriv_surf_dist[MAX_PAIR], double deriv_vol_dist[MAX_PAIR]);

extern void F77Name2(write_del)(char delaunay_file[FSIZE]);
extern void F77Name2(write_alpha)(char alpha_file[FSIZE]);


/**************************************************************************************/
/* command line */

int command_line (char *List[], int ListLength)
{
	register int i;
	static int switchval;
	static int OutputFileName = NO;
	static int InputFileName = NO;

	static char OPTION;
	static char *cmd_name;

	cmd_name = List[0];

	for ( i=1; i<ListLength; i++) {
		if( *(List[i]) == '-' ) {
			OPTION = (*(List[i]+1));
			if(OPTION == 'o') {
				OutputFileName = YES;
				strcpy(output_file,List[i+1]);
				i++;
			}
			else if(OPTION == 'a') {
				AlphaFileName = YES;
				strcpy(alpha_file,List[i+1]);
				i++;
			}
			else if(OPTION == 'd') {
				DelaunayFileName = YES;
				strcpy(delaunay_file,List[i+1]);
				i++;
			}
			else if( OPTION == 's') {
				switchval = atoi(List[i+1]);
				i++;
			}
			else if( OPTION == 'h') {
				fprintf(stderr,"\n");
				fprintf (stderr, usage, cmd_name);
				exit(1);
			}
			else {
				fprintf(stderr,"\n");
				fprintf(stderr,"Allowed options: h, o, s, d and a. You used: %s\n",&OPTION);
				fprintf (stderr, usage, cmd_name);
				exit(1);
			}
		}
		else {
			InputFileName = YES;
			strcpy(input_file,List[i]);
		}
	}
	if(InputFileName == NO)
	{
		fprintf(stderr,"\n");
		fprintf(stderr,"You did not provide an input file name\n");
		fprintf (stderr, usage, cmd_name);
		exit(1);
	}

	if (OutputFileName == 0) sprintf(output_file,"%s.vol",input_file);

	return(switchval); /* return switch selection */
}
	

/***************************************************************************************/
/* read data points*/

int read_data()
{
	int i,j;
	double x,y,z,r;
	char line[MAX_LENGTH];
	char* val;
	FILE* fp;

	if ( (fp = fopen(input_file,"r")) == NULL ) return ERROR;

	val=fgets(line, MAX_LENGTH, fp);
	sscanf(line,"%d",&npoint);
	val=fgets(line, MAX_LENGTH, fp);
	F77Name2(vertex_zone).npoints = npoint;

	for (i=0;i<npoint;i++)
	{
		val=fgets(line, MAX_LENGTH, fp);
		sscanf(line,"%d %lf %lf %lf %lf",&j,&x,&y,&z,&r);
//		printf("%d %d %f %f %f %f\n",i,j,x,y,z,r);
		F77Name2(xyz_vertex).coord[3*i]=x;
		F77Name2(xyz_vertex).coord[3*i+1]=y;
		F77Name2(xyz_vertex).coord[3*i+2]=z;
		F77Name2(xyz_vertex).radius[i]=r;
	}

	fclose(fp);
	printf("Data read properly !!\n");
	return SAFE;
}

/***************************************************************************************/
/* wrapper to Fortran "setup" function that prepares for Volume computation */

int alpha_setup()
{
	int i;
	int ival;

	nat_orig = npoint;
	F77Name2(adjust_nsphere)();
	F77Name1(setup)();

	ival = 1;

	return(ival);
}

/***************************************************************************************/
/* wrapper to Fortran functions that compute the Weighted Alpha Shape of the molecule */
/* and get surface, volume and possibly derivatives of surface of volume, if needed  */

int alpha_geom(int switchval)
{
	int ival,i;
	int idx1,idx2,idx3,idx4;
	int iadd,ired,irem;
	int ntet_new, tet_new[MAX_NEW];
	double alpha_val=0;
	double coordat[3], radiusat;

	clock_t starttime, endtime;
	double timediff;

	ival = 1;

//	printf("Start regular triangulation :\n");
	starttime=clock();
	F77Name1(regular3)(&nredundant,list_redundant);	/* computes triangulation; stores tetrahedra*/
							/* Since it is a weighted Delaunay, we
							keep track of the redundant atoms */
	endtime=clock();
        timediff = ((double)(endtime-starttime))/CLOCKS_PER_SEC;
        printf("\nComputing time for Regular triangulation  : %lf \n",timediff);
	if(switchval == 3) return(ival);

	if(DelaunayFileName == YES)
	{
		write_delaunay();
	} 
//	printf("compute alpha complex ...\n");
	starttime=clock();
	F77Name1(alfcx)(&alpha_val, &nredalpha, list_redalpha); /* computes alpha shape, alpha=0 */
        endtime=clock();
        timediff = ((double)(endtime-starttime))/CLOCKS_PER_SEC;
        printf("Computing time for building Alpha Complex : %lf \n",timediff);
//	printf("end compute alpha complex ...\n");
	if(switchval == 4) return(ival);
	if(AlphaFileName == YES)
	{
		write_alphacx();
	} 
	F77Name2(readjust_nsphere)(&nat_orig, &nredalpha, list_redalpha);

//	printf("compute volume ...\n");
	starttime=clock();
	if(switchval == 0) {
		F77Name2(volume_only)(coef, &wsurf, &wvol, &surf, &vol, surface, volume); /* computes surface area and volume  */
	}
	else {
		if(switchval == 1){
			F77Name2(volume_deriv_coord)(coef, &wsurf, &wvol, &surf, &vol, surface, volume, 
			deriv_surf, deriv_vol); /* computes volume and surface */
		}
		else {
			F77Name2(volume_deriv_dist)(coef, &wsurf, &wvol, &surf, &vol, surface, volume, 
			deriv_surf, deriv_vol,deriv_surf_dist, deriv_vol_dist); /* computes volume and surface */
		}
	}
						 /* note that if you compute volume,
						    the surface comes out for "free"
						    (i.e. at no cost in computing time)
						*/
        endtime=clock();
        timediff = ((double)(endtime-starttime))/CLOCKS_PER_SEC;
        printf("Computing time for surface area           : %lf \n",timediff);


	return(ival);
}
/****************************************************************************************/
/* Print Delaunay triangulation                                                         */

void write_delaunay()
{
/*	This procedure outputs the Delaunay triangulation                               */

	int idx,i,nkeep;
	int idx1,idx2,idx3,idx4,iorient,ihull,ilink1,ilink2,ilink3,ilink4;
	int ilink5,ilink6,ilink7,ilink8;
	int ntetrahedron, ntriangle, nedges, nvertices;
	FILE *out;

	F77Name2(write_del)(delaunay_file);
	
}

/****************************************************************************************/
/* Print Alpha complex                                                         */

void write_alphacx()
{
/*	This procedure outputs the Alpha complex                               */

	int i,istat,nkeep;
	int idx1,idx2,idx3,idx4,iorient,ihull;
	int ntetrahedron, ntriangle, nedges, nvertices;
	FILE *out;

	
	F77Name2(write_alpha)(alpha_file);

}

/***************************************************************************************/
/* Print volume and surface area                                                       */

void print_volume(double vol, double surf)
{
	int iatom;
	FILE *out;

	if((out = fopen(output_file,"w"))== NULL) {
		printf(" Cannot open file %s to output results.\n",output_file);
		exit(1);
	}

	fprintf(out,"#   i        X          Y          Z        R         Vol         Surf\n");
	for(iatom = 0; iatom < npoint; iatom++) 
	{
		fprintf(out," %5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",iatom+1,F77Name2(xyz_vertex).coord[3*iatom+12],
		F77Name2(xyz_vertex).coord[3*iatom+13],F77Name2(xyz_vertex).coord[3*iatom+14],F77Name2(xyz_vertex).radius[iatom+4],volume[iatom+4],surface[iatom+4]);
	}
	fprintf(out,"\n");
	fprintf(out,"#                                         Sum      %10.4f %10.4f\n",vol,surf);

	fclose(out);
}
/***************************************************************************************/
/* Print volume, surface area and derivatives                                                       */

void print_deriv(double vol, double surf)
{
        int iatom;
        FILE *out;

        if((out = fopen(output_file,"w"))== NULL) {
                printf(" Cannot open file %s to output results.\n",output_file);
                exit(1);
        }

        fprintf(out,"#   i        X          Y          Z        R         Vol         Surf");
        fprintf(out,"      dVol/dx    dVol/dy    dVol/dz   dSurf/dx   dSurf/dy   dSurf/dz      \n");
        for(iatom = 0; iatom < npoint; iatom++)
        {
                fprintf(out," %5d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f",iatom+1,F77Name2(xyz_vertex).coord[3*iatom+12],
                F77Name2(xyz_vertex).coord[3*iatom+13],F77Name2(xyz_vertex).coord[3*iatom+14],F77Name2(xyz_vertex).radius[iatom+4],volume[iatom+4],surface[iatom+4]);
                fprintf(out," %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n",deriv_vol[3*iatom+12],deriv_vol[3*iatom+13],
                deriv_vol[3*iatom+14],deriv_surf[3*iatom+12],deriv_surf[3*iatom+13],deriv_surf[3*iatom+14]);
        }
        fprintf(out,"\n");
        fprintf(out,"#                                         Sum      %10.4f %10.4f\n",vol,surf);

        fclose(out);
}
/***************************************************************************************/
/* main program */

int main(int argc, char **argv)
{
	int i;
	int switchval;
	int ierr;
	int nvertices, ntetras;

/* Read command line flag, and interpret them using small procedure "command_line" */

	switchval= command_line(argv,argc); 

/* Read x,y,z and radius for each point considered */

	ierr=read_data(); 
	if(ierr == ERROR) {
		fprintf(stderr,"Problem reading input file\n");
				exit(1);
	}
	printf("Read the data ...\n");

/* Set up Alpha Shape calculation */

	ierr=alpha_setup();
	printf("Alpha set up\n");

/* The program computes a weighted sum of the contribution of each point to the 
   total surface area and volume of the molecule of interest. Here we set all weights to 1, 
   but it can be adjusted to specific needs */

	nvertices = F77Name2(vertex_zone).nvertex;
	for (i=0;i<nvertices;i++) coef[i]=i-3;
//	for (i=0;i<nvertices;i++) coef[i]=1.0;

/* Performs calculation: compute regular triangulation, filter to get alpha shape with 
   alpha = 0 then compute surface and volume. if switchval is set to 1, the procedure also 
   computes the derivatives of the surface and volume with respect to the atomic 
   coordinates                 */

	ierr=alpha_geom(switchval);

/* Write results on screen */

//	ntetras   = F77Name2(tetra_zone).ntetra;

	printf("\n");
	printf("Number of atoms                      : %d\n",nat_orig);
//	printf("Number of vertices                   : %d\n",nvertices);
	printf("Number of redundant atoms (Delaunay) : %d\n",nredundant);
//	printf("Number of tetrahedra in Delaunay     : %d\n",ntetras);
	if(nredundant > 0) 
	{
		printf("Redundant atoms in Regular triangulation          : ");
		for(i=0; i<nredundant; i++) printf(" %d ",list_redundant[i]);
		printf("\n");
	}
	if(switchval != 3) {
		printf("Number of redundant atoms (Alphacx)  : %d\n",nredalpha);
		if(nredalpha > 0) 
		{
			printf("Redundant atoms in Alpha Shape          : ");
			for(i=0; i<nredalpha; i++) printf(" %d ",list_redalpha[i]);
			printf("\n");
		}
	}
	if(switchval < 3) {
		printf("\n");
		printf("Total volume                : %12.5f\n",vol);
		printf("Total weighted volume       : %12.5f\n",wvol);
		printf("Total surface area          : %12.5f\n",surf);
		printf("Total weighted surface area : %12.5f\n",wsurf);
		printf("\n");

/* Write results in file. Only surface and volume info are written out (no derivatives info)
*/
		if(switchval != 0 ) {
			(void) print_deriv(vol,surf);
		}
		else {
			(void) print_volume(vol,surf);
		}
	}

/* Clear all gmp arrays */

	(void) clear_all_gmp_arrays(&nvertices);
	(void) clear_alf_gmp();

	return SAFE;
}
