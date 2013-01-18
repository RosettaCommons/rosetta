/**
 *
 *  @file    apbs_driver.h
 *  @author  Robert Konecny
 *  @brief   Header file for the main iAPBS driver
 *
 *  @version $Id: apbs_driver.h 389 2010-03-29 20:18:15Z rok $
 *
 */


/*! \def NATOMS 
    \brief Maximum number of atoms.
*/
#define NATOMS 150000

/**
 * @brief  Wrapper iAPBS function
 * @author Robert Konecny
 *
 * @param nat Number of atoms
 * @param x Atomic coordinate (x)
 * @param y Atomic coordinate (y)
 * @param z Atomic coordinate (z)
 * @param radius Atomic radii
 * @param charge Atomic charges
 * @param r_param Input APBS parameters (real values)
 *        0: pdie
 *        1: nsdie
 *        2: srad
 *        3: swin
 *        4: temp
 *        5: sdens
 *        6: gamma
 *        7: smpbe_vol
 *        8: smpbe_size
 * @param i_param Input APBS parameters (integer values)
 *        0: sim_type - 0=mg-manual, 1=mg-auto, 2=mg-parallel
 *        1: nlev - required if sim_type=mg-manual
 *        2: grid_centering_mode - 0=use center[0-2] values, 1=center on mol1
 *        3: coarse_centering_mode - 0=use ccenter[0-2] values, 1=center on mol 1
 *        4: fine_centering_mode - 0=use fcenter[0-2] values, 1= center on mol 1
 *        5: chgm - 0=sp10, 1=sp12, 2=sp14
 *        6: pbe_mode - 0=lpbe, 1=npbe, 2=lrpbe, 3=nrpbe, 4=smpbe
 *        7: bcfl - 0=zero, 1=sdh, 2=mdh, 4=focus
 *        8: srfm - 0=mol, 1=smol, 2=sp12, 3=sp14
 *        9: calcforce - 0=no, 1=yes
 *       10: calcenergy -  0=no, 1=yes
 *       11: write_pot - 0=no, 1=yes
 *       12: write_charge - 0=no, 1=yes
 *       13: write_smol - 0=no, 1=yes
 *       14: write_kappa - 0=no, 1=yes
 *       15: write_diel - 0=no, 1=yes
 *       16: write_atompot - 0=no, 1=yes
 *       17: use_pot - 0=no, 1=yes
 *       18: NOT USED
 *       19: apol_calcforce - 0=no, 1=yes
 *       20: apol_calcenergy - 0=no, 1=yes
 *       21: nIons - number of ions
 *       22: use_charge - 0=no, 1=yes
 *       23: use_kappa - 0=no, 1=yes
 *       24: use_diel - 0=no, 1=yes
 * @param grid Grid spacing
 * @param dime Grid dimensions
 * @param pdime Grid of processors to be used in calculation
 * @param glen Grid side lengths
 * @param center Grid center
 * @param cglen Coarse grid side lengths
 * @param fglen Fine grid side lengths
 * @param ccenter Coarse grid center
 * @param fcenter Fine grid center
 * @param ofrac Overlap fraction between procs
 * @param dbg Debug verbosity flag
 * @param ionq Mobile ion charge
 * @param ionc Mobile ion concentration
 * @param ionr Mobile ion radius
 * @param (out) esEnergy Electrostatic energy
 *     This must be allocated for the number of simulations by the caller
 * @param (out) npEnergy Non-polar energy
 *     This must be allocated for the number of simulations by the caller
 * @param (out) apbsdx Total electrostatic force per atom (x coordinate)
 * @param (out) apbsdy Total electrostatic force per atom (y coordinate)
 * @param (out) apbsdz Total electrostatic force per atom (z coordinate)
 * @param (out) apbsqfx Fixed charge force (x)
 * @param (out) apbsqfy Fixed charge force (y)
 * @param (out) apbsqfz Fixed charge force (z)
 * @param (out) apbsibx Ionic boundary force (x)
 * @param (out) apbsiby Ionic boundary force (y)
 * @param (out) apbsibz Ionic boundary force (z)
 * @param (out) apbsnpx Non-polar force (x)
 * @param (out) apbsnpy Non-polar force (y)
 * @param (out) apbsnpz Non-polar force (z)
 * @param (out) apbsdbx Dielectric boundary force (x)
 * @param (out) apbsdby Dielectric boundary force (y)
 * @param (out) apbsdbz Dielectric boundary force (z)
 *     The above apbsXXX[] must be allocated for at least the number of actual atoms and less than NATOMS
 *     if calcforce > 0.  Otherwise this can be null.
 * @param apbsgrid_meta grid meta data
 *        0: ndata (input)
 *             0 <= ndata <= number of grid-dimensioned data sets to be returned in apbsgrid**.
 *             If 0, "write" statements results will go to files.
 *        1: nx (out) actual x-dimension
 *        2: ny (out) actual y-dimension
 *        3: nz (out) actual z-dimension
 *        4: hx (out) actual x-grid resolution
 *        5: hy (out) actual y-grid resolution
 *        6: hz (out) actual z-grid resolution
 *        7: centx (out) actual x center coord
 *        8: centy (out) actual y center coord
 *        9: centz (out) actual z center coord
 *       10: minx (out) actual x min coord
 *       11: miny (out) actual y min coord
 *       12: minz (out) actual z min coord
 * @param (out) apbsgrid List of grid data sets
 *     The memory must be allocated by the user when apbsgrid_meta[0]>0,
 *     or the program shall terminate.
 *
 *     apbsgrid[0][0] points to an array of data corresponding to
 *     the first TRE value in i_param[11:16].
 *     e.g. If &i_param+11 = {0,1,0,1,0,0}, indicating "write_charge" and "write_kappa" = 1,
 *     apbsgrid should be allocated by the user to 2*x*ny*nz * sizeof(double).
 *
 *   The grid data is arranged in z, y, x in the order of increment frequency.
 *   e.g. apbsgrid[0][0] =: charge of grid location (0,0,0).
 *   abpsgrid[0][1] =: (0,0,1), apbsgrid[0][2] =: (0,0,2), ...,
 *   apbsgrid[0][nz+1] = (0,1,0), apbsgrid[0][nz+2] =: (0,1,1), ...
 *
 *      // This will print out the grid data in the same order as WRITE statement
 *      // would write to file.
 *      //
 *      int i, j, k, line=0;
 *	    double * data = apbsgrid[0];
 *	    for (i=0; i<nx; i++) {
 *	      for (j=0; j<ny; j++) {
 *	        for (k=0; k<nz; k++) {
 *	          u = k*(nx)*(ny)+j*(nx)+i;
 *	          printf("%12.6e ", data[u]);
 *	          icol++;
 *	          if (icol == 3) {
 *		        icol = 0;
 *		        lines++;
 *		         printf("\n");
 *            }
 *	        }
 *	      }
 *	    }
 * @return  1 if successful, 0 otherwise 
 *
 * All output data structures must be allocated memory by the caller if they are r
 */

VEXTERNC int apbsdrv_(int *nat,
		      double x[NATOMS], 
		      double y[NATOMS], 
		      double z[NATOMS], 
		      double radius[NATOMS], 
		      double charge[NATOMS], 
		      double r_param[9], 
		      int i_param[25],
		      double grid[3],
		      int dime[3], 
		      int pdime[3], 
		      double glen[3], 
		      double center[3], 
		      double cglen[3],
		      double fglen[3], 
		      double ccenter[3], 
		      double fcenter[3], 
		      double *ofrac, 
		      int *dbg, 
		      double ionq[MAXION], 
		      double ionc[MAXION], 
		      double ionr[MAXION], 
		      double esEnergy[NOSH_MAXCALC],
		      double npEnergy[NOSH_MAXCALC],
		      double apbsdx[NATOMS], 
		      double apbsdy[NATOMS], 
		      double apbsdz[NATOMS],
		      double apbsqfx[NATOMS], 
		      double apbsqfy[NATOMS],
		      double apbsqfz[NATOMS],
		      double apbsibx[NATOMS], 
		      double apbsiby[NATOMS],
		      double apbsibz[NATOMS],
		      double apbsnpx[NATOMS], 
		      double apbsnpy[NATOMS],
		      double apbsnpz[NATOMS],
		      double apbsdbx[NATOMS], 
		      double apbsdby[NATOMS],
		      double apbsdbz[NATOMS],
		      double apbsgrid_meta[13],
		      double * apbsgrid[]);

/**
 * @brief  Calculate forces from MG solution
 * @author Robert Konecny (based on forceMG)
 *
 * @param  mem  Memory management object
 * @param  nosh  Parameters from input file
 * @param  pbeparm  Generic PBE parameters
 * @param  mgparm  MG-specific parmaeters
 * @param  pmg  MG object
 * @param  nforce 0 => no forces, 1 => net forces, >1 => number of
 *                       forces (1 per atom)
 * @param  atomForce Pointer to array of force objects
 * @param  alist  List of atom lists
 * @param  debug verbosity flag
 * @return  1 if successful, 0 otherwise */
VEXTERNC int iforceMG(Vmem *mem, 
		      NOsh *nosh, 
		      PBEparm *pbeparm,  
		      MGparm *mgparm,
		      Vpmg *pmg, 
		      int *nforce, 
		      AtomForce **atomForce, 
		      Valist *alist[NOSH_MAXMOL], 
		      int debug);

/**
* @brief  Combine and pretty-print energy data
* @ingroup  Frontend
* @author  David Gohara
* @return  1 if successful, 0 otherwise */
//VEXTERNC double getElecEnergy(

double getElecEnergy(
	Vcom *com, /** Communications object */
	NOsh *nosh, /** Parameters from input file */
	double totEnergy[NOSH_MAXCALC], /** Array of energies
					  from different calculations */
	int iprint /** Index of energy statement to print */
	);


/**
* @brief Creates APBS input string
* @author Robert Konecny
* @param r_param Input APBS parameters (real values)
* @param i_param Input APBS parameters (integer values)
* @param grid Grid spacing
* @param dime Grid dimensions
* @param ionq Mobile ion charge
* @param ionc Mobile ion concentration
* @param ionr Mobile ion radius
* @param glen Grid side lengths
* @param center Grid center
* @param cglen Coarse grid side lengths
* @param fglen Fine grid side lengths
* @param ccenter Coarse grid center
* @param fcenter Fine grid center
* @param ofrac Overlap fraction between procs
* @param pdime Grid of processors to be used in calculation
* @param debug Debug verbosity level.
* @return the input string */
char *setupString(double r_param[9], int i_param[25], double grid[3], 
	int dime[3], double ionq[MAXION], double ionc[MAXION],
	double ionr[MAXION], double glen[3], double center[3],
	double cglen[3], double fglen[3],
	double ccenter[3], double fcenter[3], double *ofrac, 
	int pdime[3], int debug);
