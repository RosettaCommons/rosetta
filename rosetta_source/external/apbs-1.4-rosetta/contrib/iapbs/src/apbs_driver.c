/**
 *
 *  @file   apbs_driver.c
 *  @author Robert Konecny
 *  @brief  The main iAPBS driver code
 *  @note   Energy is returned in kJ/mol, forces in  kJ/(mol/A).
 *
 *  $Revision: 556 $
 *  $Id: apbs_driver.c 556 2012-01-10 03:03:33Z rok $
 *
 */
  
#include "apbs.h"  
#include "routines.h"
#include "generic/nosh.h"  
#include "generic/mgparm.h"  
#include "generic/pbeparm.h"  
#include "generic/femparm.h"
#include "generic/vhal.h"

#include "apbs_driver.h"


/*! \def MAX_BUF_SIZE
    \brief Buffer size for internal APBS string input.
*/
#define MAX_BUF_SIZE 4096

VEMBED(rcsid="$Id: apbs_driver.c 556 2012-01-10 03:03:33Z rok $")

/**
 * @brief  Wrapper iAPBS function
 * @author Robert Konecny
 *
 */
int apbsdrv_(
	     int *nat,
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
	     double esenergy[1],
	     double npenergy[1],
	     double apbsdx[NATOMS], double apbsdy[NATOMS], double apbsdz[NATOMS],
	     double apbsqfx[NATOMS], double apbsqfy[NATOMS], double apbsqfz[NATOMS],
	     double apbsibx[NATOMS], double apbsiby[NATOMS], double apbsibz[NATOMS],
	     double apbsnpx[NATOMS], double apbsnpy[NATOMS], double apbsnpz[NATOMS],
	     double apbsdbx[NATOMS], double apbsdby[NATOMS], double apbsdbz[NATOMS],
	     double apbsgrid_meta[13],
	     double *apbsgrid[]
	     )
{
    int i,k,j;

    // 1 if grid data to be written to files (traditional), 0 to return via apbsgrid**.
    int is_grid2file = 1;

    // Memory must be allocated for the outputs if they are to be returned in-memory.
    if( apbsgrid_meta[0] > 0 ) {
        is_grid2file = 0;
        VASSERT( apbsgrid != VNULL );
        for(i=0; i<apbsgrid_meta[0]; i++) {
            VASSERT( apbsgrid[i] != VNULL );
        }
    }

    NOsh *nosh = VNULL;

    MGparm *mgparm = VNULL;
    FEMparm *feparm = VNULL;
    PBEparm *pbeparm = VNULL;
    APOLparm *apolparm = VNULL;
    Vparam *param = VNULL;

    Vmem *mem = VNULL;
    Vcom *com = VNULL;
    Vio *sock = VNULL;
#ifdef HAVE_MC_H
    Vfetk *fetk[NOSH_MAXCALC];
    Gem *gm[NOSH_MAXMOL];
#else
    void *fetk[NOSH_MAXCALC];
    void *gm[NOSH_MAXMOL];
#endif
    Vpmg *pmg[NOSH_MAXCALC];
    Vpmgp *pmgp[NOSH_MAXCALC];
    Vpbe *pbe[NOSH_MAXCALC];
    Valist *alist[NOSH_MAXMOL];
    Vgrid *dielXMap[NOSH_MAXMOL],*dielYMap[NOSH_MAXMOL],*dielZMap[NOSH_MAXMOL];
    Vgrid *kappaMap[NOSH_MAXMOL];
    Vgrid *potMap[NOSH_MAXMOL];
    Vgrid *chargeMap[NOSH_MAXMOL];
    //    char *input_path = VNULL;
    char *output_path = VNULL;
    int rank, size, isolve, debug, natom;
    //    unsigned long int bytesTotal, highWater;
    size_t bytesTotal, highWater;
    Voutput_Format outputformat;

    int bufsize = MAX_BUF_SIZE;
    int rc = 0;
    double coord[3];
    char *inputString;

    /* These variables require some explaining... The energy double arrays
     * store energies from the various calculations.  The energy int array
     * stores either a flag (0,1) displaying whether energies were calculated
     * or if PCE_COMPS is used, the number of atom energies stored
     * for the given calculation.  Likewise, the
     * force double arrays store forces from the various calcualtions.  The
     * force int array stores an integer which either says no calculation was
     * performed (0) or gives the number of entries in the force array for each
     * calculation */
    double qfEnergy[NOSH_MAXCALC], qmEnergy[NOSH_MAXCALC];
    double dielEnergy[NOSH_MAXCALC], totEnergy[NOSH_MAXCALC];
    AtomForce *atomForce[NOSH_MAXCALC];
    double *atomEnergy[NOSH_MAXCALC];
    int nenergy[NOSH_MAXCALC], nforce[NOSH_MAXCALC];
    /* The real partition centers */
    double realCenter[3];



    /* ************** CHECK PARALLEL STATUS *************** */
    // init is done by the calling program
    //    VASSERT(Vcom_init(&argc, &argv)); 
    com = Vcom_ctor(1);
    rank = Vcom_rank(com);
    size = Vcom_size(com);
    startVio(); 
    Vnm_setIoTag(rank, size);
    Vnm_tprint( 0, "Hello world from PE %d\n", rank);

    /* A bit of array/pointer initialization */
    mem = Vmem_ctor("MAIN");
    for (i=0; i<NOSH_MAXCALC; i++) {
	pmg[i] = VNULL;
	pmgp[i] = VNULL;
	fetk[i] = VNULL;
	pbe[i] = VNULL;
	qfEnergy[i] = 0;
	qmEnergy[i] = 0;
	dielEnergy[i] = 0;
	totEnergy[i] = 0;
	atomForce[i] = VNULL;
	nenergy[i] = 0;
	nforce[i] = 0;
    }
    for (i=0; i<NOSH_MAXMOL; i++) {
	alist[i] = VNULL;
	dielXMap[i] = VNULL;
	dielYMap[i] = VNULL;
	dielZMap[i] = VNULL;
	kappaMap[i] = VNULL;
	potMap[i] = VNULL;
	chargeMap[i] = VNULL;
    }


    /* *************** PARSE INPUT FILE ******************* */
    nosh = NOsh_ctor(rank, size);

//    sock = Vio_ctor("FILE", "ASC", VNULL, input_path, "r");
//    Vnm_tprint( 1, "Parsing input file %s...\n", input_path);
    
    VASSERT( bufsize <= VMAX_BUFSIZE );
    sock = Vio_ctor("BUFF","ASC",VNULL,"0","r");

    debug = *dbg;

    /* generate input string */
    inputString = VNULL;
    inputString = setupString(r_param, i_param, grid, dime, ionq, ionc, 
		  ionr, glen, center, cglen, fglen, ccenter, fcenter, ofrac, 
		  pdime, debug);
    if(debug>2) Vnm_tprint(1, "debug: Input string:\n%s\n", inputString);
    Vio_bufTake(sock, inputString, bufsize);

    if (!NOsh_parseInput(nosh, sock)) {
	Vnm_tprint( 2, "Error while parsing input string.\n");
	VJMPERR1(0);
    } else if(debug>1) Vnm_tprint( 1, "Parsed input string.\n");
    
    sock->VIObuffer = VNULL;
    Vio_dtor(&sock);


    /* *************** LOAD PARAMETERS AND MOLECULES ******************* */
    //nosh->gotparm = 0; // not using param file for now
/*    
 * 
      param = loadParameter(nosh);
      if (loadMolecules(nosh, param, alist) != 1) {
	Vnm_tprint(2, "Error reading molecules!\n");
	VJMPERR1(0);
    }
*/

    /* alist fills nosh */
    alist[0] = Valist_ctor();

    alist[0]->center[0] = 0.;
    alist[0]->center[1] = 0.;
    alist[0]->center[2] = 0.;
    alist[0]->maxcrd[0] = -VLARGE;
    alist[0]->maxcrd[1] = -VLARGE;
    alist[0]->maxcrd[2] = -VLARGE;
    alist[0]->mincrd[0] = VLARGE;
    alist[0]->mincrd[1] = VLARGE;
    alist[0]->mincrd[2] = VLARGE;
    alist[0]->maxrad = 0.;
    alist[0]->charge = 0.;

    alist[0]->number = *nat;
    natom =  alist[0]->number;
    /* Allocate the necessary space for the atom array */
    alist[0]->atoms = Vmem_malloc(alist[0]->vmem, alist[0]->number,
	    (sizeof(Vatom)));
    VASSERT(alist[0]->atoms != VNULL);


    for (i=0; i<alist[0]->number; i++) {
	/* Fill atoms in the atom list */

	/* atom[i]->partID = 1; FIXME */

	if (x[i] < alist[0]->mincrd[0]) alist[0]->mincrd[0] = x[i];
	if (y[i] < alist[0]->mincrd[1]) alist[0]->mincrd[1] = y[i];
	if (z[i] < alist[0]->mincrd[2]) alist[0]->mincrd[2] = z[i];
	if (x[i] > alist[0]->maxcrd[0]) alist[0]->maxcrd[0] = x[i];
	if (y[i] > alist[0]->maxcrd[1]) alist[0]->maxcrd[1] = y[i];
	if (z[i] > alist[0]->maxcrd[2]) alist[0]->maxcrd[2] = z[i];
	if (radius[i] > alist[0]->maxrad) alist[0]->maxrad = radius[i];
	alist[0]->charge = alist[0]->charge + charge[i];

	coord[0] = x[i];
	coord[1] = y[i];
	coord[2] = z[i];
	Vatom_setPosition(&(alist[0]->atoms)[i], coord);
	Vatom_setCharge(&(alist[0]->atoms)[i], charge[i]);
	Vatom_setRadius(&(alist[0]->atoms)[i], radius[i]);
	Vatom_setAtomID(&(alist[0]->atoms)[i], i);
	// not necessary? Vatom_setPartID(&(alist[0]->atoms)[i, );
    }

    alist[0]->center[0] = 0.5*(alist[0]->maxcrd[0] + alist[0]->mincrd[0]);
    alist[0]->center[1] = 0.5*(alist[0]->maxcrd[1] + alist[0]->mincrd[1]);
    alist[0]->center[2] = 0.5*(alist[0]->maxcrd[2] + alist[0]->mincrd[2]);

    /* *************** SETUP CALCULATIONS *************** */
    if (NOsh_setupElecCalc(nosh, alist) != 1) {
	Vnm_tprint(2, "Error setting up ELEC calculations\n");
	VJMPERR1(0);
    }

    if ((rc = NOsh_setupApolCalc(nosh, alist)) == ACD_ERROR) {
	Vnm_tprint(2, "Error setting up APOL calculations\n");
	VJMPERR1(0);
    }

    /* ******************* CHECK APOL********************** */
    //if((nosh->gotparm == 0) && (rc == ACD_YES)){
    //	Vnm_print(1,"\nError you must provide a parameter file if you\n" \
    //		"     are performing an APOLAR calculation\n");
    //	VJMPERR1(0);
    //}

    /* *************** LOAD MAPS ******************* */
    if (loadDielMaps(nosh, dielXMap, dielYMap, dielZMap) != 1) {
	Vnm_tprint(2, "Error reading dielectric maps!\n");
	VJMPERR1(0);
    }
    if (loadKappaMaps(nosh, kappaMap) != 1) {
	Vnm_tprint(2, "Error reading kappa maps!\n");
	VJMPERR1(0);
    }
    if (loadPotMaps(nosh, potMap) != 1) {
      Vnm_tprint(2, "Error reading potential maps!\n");
      VJMPERR1(0);
    }
    if (loadChargeMaps(nosh, chargeMap) != 1) {
	Vnm_tprint(2, "Error reading charge maps!\n");
	VJMPERR1(0);
    }

    /* *************** Initialization ******************* */
    // Do this initialization only if calcforce is requested.
    for (i=0; i<nosh->ncalc; i++){
        pbeparm = nosh->calc[i]->pbeparm;
        if(pbeparm->calcforce >0) {
            for (i=0; i < natom; i++) {
                apbsdx[i] = 0.0;
                apbsdy[i] = 0.0;
                apbsdz[i] = 0.0;
                apbsqfx[i] = 0.0;
                apbsqfy[i] = 0.0;
                apbsqfz[i] = 0.0;
                apbsibx[i] = 0.0;
                apbsiby[i] = 0.0;
                apbsibz[i] = 0.0;
                apbsdbx[i] = 0.0;
                apbsdby[i] = 0.0;
                apbsdbz[i] = 0.0;
                apbsnpx[i] = 0.0;
                apbsnpy[i] = 0.0;
                apbsnpz[i] = 0.0;
            }
            // once is enough
            break;
        }
    }

    /* *************** DO THE CALCULATIONS ******************* */
    if(debug>1) Vnm_tprint( 1, "Preparing to run %d PBE calculations.\n",
	    nosh->ncalc);
    for (i=0; i<nosh->ncalc; i++) {
	if(debug>1) Vnm_tprint( 1, "----------------------------------------\n");

	switch (nosh->calc[i]->calctype) {  
	    case NCT_MG:
		/* What is this?  This seems like a very awkward way to find 
		   the right ELEC statement... */
		for (k=0; k<nosh->nelec; k++) {
		    if (nosh->elec2calc[k] >= i) {
			break;
		    }
		}
		if (Vstring_strcasecmp(nosh->elecname[k], "") == 0) {
		    if(debug>1) Vnm_tprint( 1, "CALCULATION #%d: MULTIGRID\n", i+1);
		} else {
		    if(debug>1) Vnm_tprint( 1, "CALCULATION #%d (%s): MULTIGRID\n", 
			    i+1, nosh->elecname[k]);
		}
		/* Useful local variables */
		mgparm = nosh->calc[i]->mgparm;
		pbeparm = nosh->calc[i]->pbeparm;
		natom =  alist[0]->number;

		/* Set up problem */
		if(debug>1) Vnm_tprint( 1, "  Setting up problem...\n");
		if (!initMG(i, nosh, mgparm, pbeparm, realCenter, pbe, 
			    alist, dielXMap, dielYMap, dielZMap, kappaMap, chargeMap, 
			    pmgp, pmg, potMap)) {
		    Vnm_tprint( 2, "Error setting up MG calculation!\n");
		    VJMPERR1(0);
		}

		/* Print problem parameters */
		if(debug>0) printMGPARM(mgparm, realCenter);
		if(debug>0) printPBEPARM(pbeparm);

		/* Solve PDE */
		if (solveMG(nosh, pmg[i], mgparm->type) != 1) {
		    Vnm_tprint(2, "Error solving PDE!\n");
		    VJMPERR1(0);
		}

		/* Set partition information for observables and I/O */
		if (setPartMG(nosh, mgparm, pmg[i]) != 1) {
		    Vnm_tprint(2, "Error setting partition info!\n");
		    VJMPERR1(0);
		}

		/* Write out energies */
		energyMG(nosh, i, pmg[i], 
			&(nenergy[i]), &(totEnergy[i]), &(qfEnergy[i]), 
			&(qmEnergy[i]), &(dielEnergy[i]));
		esenergy[0] = 0.0;
		esenergy[0] = getElecEnergy(com, nosh, totEnergy, i);
		
		//		if(debug>3) printElecEnergy(com, nosh, totEnergy, i);

		/* Write out forces */
		forceMG(mem, nosh, pbeparm, mgparm, pmg[i], &(nforce[i]), 
			&(atomForce[i]), alist);
		
		//if (pbeparm->calcforce > 0 && i == nosh->ncalc-1) {
		if (pbeparm->calcforce == PCF_TOTAL) {
		  //		  Vnm_tprint(2, "pcf_total\n");
		  apbsdx[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    (atomForce[i][0].qfForce[0] +
		     atomForce[i][0].ibForce[0] +
		     atomForce[i][0].dbForce[0]);
		  apbsdy[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    (atomForce[i][0].qfForce[1] +
		     atomForce[i][0].ibForce[1] +
		     atomForce[i][0].dbForce[1]);
		  apbsdz[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    (atomForce[i][0].qfForce[2] +
		     atomForce[i][0].ibForce[2] +
		     atomForce[i][0].dbForce[2]);

		  apbsqfx[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    atomForce[i][0].qfForce[0];
		  apbsqfy[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    atomForce[i][0].qfForce[1];
		  apbsqfz[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    atomForce[i][0].qfForce[2];

		  apbsibx[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    atomForce[i][0].ibForce[0];
		  apbsiby[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    atomForce[i][0].ibForce[1];
		  apbsibz[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    atomForce[i][0].ibForce[2];

		  apbsdbx[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    atomForce[i][0].dbForce[0];
		  apbsdby[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    atomForce[i][0].dbForce[1];
		  apbsdbz[0] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
		    atomForce[i][0].dbForce[2];
		}

		if (pbeparm->calcforce == PCF_COMPS) {
		  //		  Vnm_tprint(2, "pcf_comps\n");
		    for (j=0; j < natom; j++) {

			apbsdx[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    (atomForce[i][j].qfForce[0] +
			     atomForce[i][j].ibForce[0] +
			     atomForce[i][j].dbForce[0]);
			apbsdy[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    (atomForce[i][j].qfForce[1] +
			     atomForce[i][j].ibForce[1] +
			     atomForce[i][j].dbForce[1]);
			apbsdz[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    (atomForce[i][j].qfForce[2] +
			     atomForce[i][j].ibForce[2] +
			     atomForce[i][j].dbForce[2]);

			/* individual components */
			apbsqfx[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    atomForce[i][j].qfForce[0];
			apbsqfy[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    atomForce[i][j].qfForce[1];
			apbsqfz[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    atomForce[i][j].qfForce[2];

			apbsibx[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    atomForce[i][j].ibForce[0];
			apbsiby[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    atomForce[i][j].ibForce[1];
			apbsibz[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    atomForce[i][j].ibForce[2];

			apbsdbx[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    atomForce[i][j].dbForce[0];
			apbsdby[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    atomForce[i][j].dbForce[1];
			apbsdbz[j] = Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na *
			    atomForce[i][j].dbForce[2];
		    }
		    //		    if(debug>3) printElecForce(com, nosh, nforce, atomForce, i);
		}
		/* Write grid-dimensioned data to file if that's what user wants*/
		if( is_grid2file ) {
		  writedataMG(rank, nosh, pbeparm, pmg[i]);
		}
		/* Return in-memeory instead */
		else {		
		  for( k=0; k<pbeparm->numwrite; k++ ) {
		    int nx = pmg[i]->pmgp->nx;
		    int ny = pmg[i]->pmgp->ny;
		    int nz = pmg[i]->pmgp->nz;
		    double hx = pmg[i]->pmgp->hx;
		    double hy = pmg[i]->pmgp->hy;
		    double hz = pmg[i]->pmgp->hzed;
		    double centx = pmg[i]->pmgp->xcent;
		    double centy = pmg[i]->pmgp->ycent;
		    double centz = pmg[i]->pmgp->zcent;
		    double parm;
		    Vdata_Type data_type = pbeparm->writetype[k];
		    switch(data_type) {
		    case VDT_SMOL:
		      parm = pbeparm->srad;
		      break;
		    case VDT_SSPL:
		      parm = pbeparm->swin;
		      break;
		    case VDT_IVDW:
		      parm = pmg[i]->pbe->maxIonRadius;
		      break;
		    case VDT_DIELX:
		      centx +=0.5*hz;
		      parm = 0;
		      break;
		    case VDT_DIELY:
		      centy +=0.5*hz;
		      parm = 0;
		      break;
		    case VDT_DIELZ:
		      centz +=0.5*hz;
		      parm = 0;
		      break;
		    case VDT_CHARGE:
		    case VDT_POT:
		    case VDT_VDW:
		    case VDT_LAP:
		    case VDT_EDENS:
		    case VDT_NDENS:
		    case VDT_QDENS:
		    case VDT_KAPPA:
		    case VDT_ATOMPOT:
		      parm = 0;
		      break;
		    default:
		      Vnm_print(2, "Warning!  Skipping invalid data type for writing: %d\n", k);
		      continue;
		    }

		    apbsgrid_meta[1] = nx;
		    apbsgrid_meta[2] = ny;
		    apbsgrid_meta[3] = nz;
		    apbsgrid_meta[4] = hx;
		    apbsgrid_meta[5] = hy;
		    apbsgrid_meta[6] = hz;
		    apbsgrid_meta[7] = centx;
		    apbsgrid_meta[8] = centy;
		    apbsgrid_meta[9] = centz;
		    apbsgrid_meta[10] = centx - 0.5*(nx-1)*hx;
		    apbsgrid_meta[11] = centy - 0.5*(ny-1)*hy;
		    apbsgrid_meta[12] = centz - 0.5*(nz-1)*hz;

		    Vpmg_fillArray(pmg[i], apbsgrid[k], data_type, parm,
				   pbeparm->pbetype, pbeparm);
		  }
		}
		/* Write matrix */
		writematMG(rank, nosh, pbeparm, pmg[i]);

		/* If needed, cache atom energies */				
		nenergy[i] = 0;
		if ((pbeparm->calcenergy == PCE_COMPS) && (outputformat != OUTPUT_NULL)){
		    storeAtomEnergy(pmg[i], i, &(atomEnergy[i]), &(nenergy[i]));
		}

		/* clean up memory after final run - gets around APBS memory leak */
		if (i == nosh->elec2calc[k] ) {
		  if(debug>4) Vnm_tprint( 1, "Cleaning %d pmg memory segment\n", i);
		  Vpmg_dtor(&pmg[i]);
		  killMG(nosh, pbe, pmgp, pmg);
		}

		fflush(stdout);
		fflush(stderr);

		break;

		/* ***** Do FEM calculation ***** */
	    case NCT_FEM:
#ifdef HAVE_MC_H
		for (k=0; k<nosh->nelec; k++) {
		    if (nosh->elec2calc[k] >= i) break;
		}
		if (Vstring_strcasecmp(nosh->elecname[i+1], "") == 0) {
		    Vnm_tprint( 1, "CALCULATION #%d: FINITE ELEMENT\n", i+1);
		} else {
		    Vnm_tprint( 1, "CALCULATION #%d (%s): FINITE ELEMENT\n", i+1, nosh->elecname[k+1]);
		}

		/* Useful local variables */
		feparm = nosh->calc[i]->femparm;
		pbeparm = nosh->calc[i]->pbeparm;

		/* Warn the user about some things */
		Vnm_tprint(2, "#################### WARNING ###################\n");
		Vnm_tprint(2, "## FE support is currently very experimental! ##\n");
		Vnm_tprint(2, "#################### WARNING ###################\n");

		/* Set up problem */
		Vnm_tprint( 1, "  Setting up problem...\n");
		if (!initFE(i, nosh, feparm, pbeparm, pbe, alist, fetk)) {
		    Vnm_tprint( 2, "Error setting up FE calculation!\n");
		    VJMPERR1(0);
		}

		/* Print problem parameters */
		printFEPARM(i, nosh, feparm, fetk);
		printPBEPARM(pbeparm);

		/* Refine mesh */
		if (!preRefineFE(i, nosh, feparm, fetk)) {
		    Vnm_tprint( 2, "Error pre-refining mesh!\n");
		    VJMPERR1(0);
		}

		/* Solve-estimate-refine */
		Vnm_tprint(2, "\n\nWARNING!  DO NOT EXPECT PERFORMANCE OUT OF THE APBS/FEtk\n");
		Vnm_tprint(2, "INTERFACE AT THIS TIME.  THE FINITE ELEMENT SOLVER IS\n");
		Vnm_tprint(2, "CURRENTLY NOT OPTIMIZED FOR THE PB EQUATION.  IF YOU WANT\n");
		Vnm_tprint(2, "PERFORMANCE, PLEASE USE THE MULTIGRID-BASED METHODS, E.G.\n");
		Vnm_tprint(2, "MG-AUTO, MG-PARA, and MG-MANUAL (SEE DOCS.)\n\n");
		Vnm_tprint(1, "  Beginning solve-estimate-refine cycle:\n");
		for (isolve=0; isolve<feparm->maxsolve; isolve++) {
		    Vnm_tprint(1, "    Solve #%d...\n", isolve);
		    if (!solveFE(i, nosh, pbeparm, feparm, fetk)) {
			Vnm_tprint(2, "ERROR SOLVING EQUATION!\n");
			VJMPERR1(0);
		    }
		    if (!energyFE(nosh, i, fetk, &(nenergy[i]), 
				&(totEnergy[i]), &(qfEnergy[i]), 
				&(qmEnergy[i]), &(dielEnergy[i]))) {
			Vnm_tprint(2, "ERROR SOLVING EQUATION!\n");
			VJMPERR1(0);
		    }
		    /* We're not going to refine if we've hit the max number
		     * of solves */
		    if (isolve < (feparm->maxsolve)-1) {
			if (!postRefineFE(i, nosh, feparm, fetk)) break;
		    }
		    bytesTotal = Vmem_bytesTotal();
		    highWater = Vmem_highWaterTotal();
		    Vnm_tprint(1, "      Currently memory use:  %g MB\n", 
			    ((double)bytesTotal/(1024.)/(1024.)));
		    Vnm_tprint(1, "      High-water memory use:  %g MB\n", 
			    ((double)highWater/(1024.)/(1024.)));
		}

		Vnm_tprint(1, "  Writing FEM data to files.\n");
		if (!writedataFE(rank, nosh, pbeparm, fetk[i])) {
		    Vnm_tprint(2, "  Error while writing FEM data!\n");
		}
#else /* ifdef HAVE_MC_H */
		Vnm_print(2, "Error!  APBS not compiled with FEtk!\n");
		exit(2);
#endif /* ifdef HAVE_MC_H */
		break;

		/* Do an apolar calculation */
	    case NCT_APOL:
		/* Copied from NCT_MG. See the note above (top of loop) for
		   information about this loop.
		   */
		for (k=0; k<nosh->napol; k++) {
		    if (nosh->apol2calc[k] >= i) {
			break;
		    }
		}

		if (Vstring_strcasecmp(nosh->apolname[k], "") == 0) {
		    if(debug>1) Vnm_tprint( 1, "CALCULATION #%d: APOLAR\n", i+1);
		} else {
		    if(debug>1) Vnm_tprint( 1, "CALCULATION #%d (%s): APOLAR\n", 
			    i+1, nosh->apolname[k]);
		}

		apolparm = nosh->calc[i]->apolparm;
		rc = initAPOL(nosh, mem, param, apolparm, &(nforce[i]), &(atomForce[i]), 
			alist[(apolparm->molid)-1]);
		if(rc == 0) {
		    Vnm_tprint(2, "Error calculating apolar solvation quantities!\n");
		    VJMPERR1(0);
		}

		if(debug>3) printf("energyAPOL: %1.12E\n", apolparm->gamma*apolparm->sasa);
		npenergy[0] = 0.0;
		npenergy[0] = apolparm->gamma*apolparm->sasa;

		//if (pbeparm->calcforce > 0 && i == nosh->ncalc-1) {
		//		if (pbeparm->calcforce > 1) {
		if (apolparm->calcforce == ACF_COMPS) {
		    for (j=0; j < natom; j++) {
		      apbsnpx[j] = (atomForce[i][j]).sasaForce[0];
		      apbsnpy[j] = (atomForce[i][j]).sasaForce[1];
		      apbsnpz[j] = (atomForce[i][j]).sasaForce[2];
		    }
		}
		//if(debug>3) printApolEnergy(nosh, i);
		//if(debug>3) printApolForce(com, nosh, nforce, atomForce, i); 
		break;
	    default:
		Vnm_tprint(2, "  Unknown calculation type (%d)!\n", 
			   nosh->calc[i]->calctype);
		exit(2);
	}
    }

    //Clear out the parameter file memory
    if(param != VNULL) Vparam_dtor(&param);

    /* *************** HANDLE PRINT STATEMENTS ******************* */
    if (nosh->nprint > 0) {
	Vnm_tprint( 1, "----------------------------------------\n");
	Vnm_tprint( 1, "PRINT STATEMENTS\n");
    }
    for (i=0; i<nosh->nprint; i++) {
	/* Print energy */
	if (nosh->printwhat[i] == NPT_ENERGY) {
	    printEnergy(com, nosh, totEnergy, i);
	    /* Print force */
	} else if (nosh->printwhat[i] == NPT_FORCE) {
	    if(debug>6) printForce(com, nosh, nforce, atomForce, i);
	} else if (nosh->printwhat[i] == NPT_ELECENERGY) {
	    if(debug>3) printElecEnergy(com, nosh, totEnergy, i);
	    //esEnergy[0] = getElecEnergy(com, nosh, totEnergy, i);
	} else if (nosh->printwhat[i] == NPT_ELECFORCE) {
	    if(debug>6) printElecForce(com, nosh, nforce, atomForce, i);
	} else if (nosh->printwhat[i] == NPT_APOLENERGY) {
	    if(debug>3) printApolEnergy(nosh, i);
	} else if (nosh->printwhat[i] == NPT_APOLFORCE) {
	    if(debug>6) printApolForce(com, nosh, nforce, atomForce, i);
	} else {
	    Vnm_tprint( 2, "Undefined PRINT keyword!\n");
	    break;
	}
    } 
    if(debug>1) Vnm_tprint( 1, "----------------------------------------\n");

    /* *************** HANDLE LOGGING *********************** */

    if (outputformat == OUTPUT_FLAT) {
	Vnm_tprint(2," Writing data to flat file %s...\n\n", output_path);
	writedataFlat(nosh, com, output_path, totEnergy, qfEnergy, qmEnergy,
		dielEnergy, nenergy, atomEnergy, nforce, atomForce);
    }

    /* Destroy energy arrays if they still exist */

    for (i=0; i<nosh->ncalc; i++) {
	if (nenergy[i] > 0) Vmem_free(mem, nenergy[i], sizeof(double),
		(void **)&(atomEnergy[i]));    
    }

    /* *************** GARBAGE COLLECTION ******************* */

    //Vnm_tprint( 1, "CLEANING UP AND SHUTTING DOWN...\n");
    /* Clean up APBS structures */
    killForce(mem, nosh, nforce, atomForce);
    killEnergy();
    killMG(nosh, pbe, pmgp, pmg);
#ifdef HAVE_MC_H
    killFE(nosh, pbe, fetk, gm);
#endif
    killChargeMaps(nosh, chargeMap);
    killKappaMaps(nosh, kappaMap);
    killDielMaps(nosh, dielXMap, dielYMap, dielZMap);
    killMolecules(nosh, alist);
    NOsh_dtor(&nosh);

    /* Memory statistics */
    bytesTotal = Vmem_bytesTotal();
    highWater = Vmem_highWaterTotal();
    if(debug>1) Vnm_tprint( 1, "Final memory usage:  %4.3f MB total, %4.3f MB high water\n", 
	    (double)(bytesTotal)/(1024.*1024.), 
	    (double)(highWater)/(1024.*1024.));

    /* Clean up MALOC structures */
    Vcom_dtor(&com);
    Vmem_dtor(&mem);

    /* And now it's time to so "so long"... */
    //Vnm_tprint(1, "\n\n");
    //Vnm_tprint( 1, "Thanks for using APBS!\n\n");

    /* This should be last */
    Vnm_tstop(APBS_TIMER_WALL_CLOCK, "APBS WALL CLOCK");
    Vnm_flush(1);
    Vnm_flush(2);
//    Vcom_finalize();
    fflush(NULL);
    return 0;

VERROR1:
    Vcom_finalize();
    Vcom_dtor(&com);
    Vmem_dtor(&mem);
    return APBSRC;
}


/**
* @brief Creates APBS input string
* @author Robert Konecny
*/
char *setupString(double r_param[9], int i_param[25], double grid[3], 
	int dime[3], double ionq[MAXION], double ionc[MAXION],
	double ionr[MAXION], double glen[3], double center[3],
	double cglen[3], double fglen[3],
	double ccenter[3], double fcenter[3], double *ofrac, 
        int pdime[3], int debug)
{
    static char string[MAX_BUF_SIZE];
    static char mybuf[MAX_BUF_SIZE];
    int i;
    // read section
    strcpy(string, "read\n mol pqr ion.pqr\n");

    if(i_param[22] == 1){
      strcat(string, "charge dx iapbs-charge.dx\n");
    }
    if(i_param[23] == 1){
      strcat(string, "kappa dx iapbs-kappa.dx\n");
    }
    if(i_param[24] == 1){
      strcat(string, "diel dx iapbs-dielx.dx iapbs-diely.dx iapbs-dielz.dx\n");
    }
    if(i_param[17] == 1){
      strcat(string, "pot dx iapbs-pot.dx\n");
    }
    strcat(string, "end\n\n");

    // elec section
    strcat(string, "elec name elec\n");
    switch(i_param[0]){
	case 0: sprintf(mybuf, "mg-manual\n"); break;
	case 1: sprintf(mybuf, "mg-auto\n"); break;
	case 2: sprintf(mybuf, "mg-para\n"); break;
	default: printf("iAPBS: mg- keyword error\n"); exit(2);
    }
    strcat(string, mybuf);

    sprintf(mybuf, "dime %i %i %i\n", dime[0], dime[1], dime[2]);
    strcat(string,mybuf);
    if (i_param[0] == 0 && grid[0] > 0.) {
	sprintf(mybuf, "grid %.3f %.3f %.3f\n", grid[0], grid[1], grid[2]);
	strcat(string,mybuf);
	// nlev is deprecated now (2009/12)
	sprintf(mybuf, "# nlev %i\n", i_param[1]);
	strcat(string,mybuf);
    }
    if (i_param[0] == 0 && glen[0] > 0.) {
	sprintf(mybuf, "glen %.3f %.3f %.3f\n", glen[0], glen[1], glen[2]);
	strcat(string,mybuf);
    }

    if (i_param[0] == 0 && i_param[2] == 0) {
	sprintf(mybuf, "gcent %.3f %.3f %.3f\n", center[0], center[1], center[2]);
	strcat(string,mybuf);
    } else if (i_param[0] == 0 && i_param[2] == 1) {
      strcat(string, "gcent mol 1\n");
    }
    if (i_param[0] != 0 && cglen[0] > 0.) {
	sprintf(mybuf, "cglen %.3f %.3f %.3f\n", cglen[0], cglen[1], cglen[2]);
	strcat(string,mybuf);
    }
    if (i_param[0] != 0 && fglen[0] > 0.) {
	sprintf(mybuf, "fglen %.3f %.3f %.3f\n", fglen[0], fglen[1], fglen[2]);
	strcat(string,mybuf);
    }
    if (i_param[0] != 0 && i_param[3] == 0) {
	sprintf(mybuf, "cgcent %.3f %.3f %.3f\n", 
		ccenter[0], ccenter[1], ccenter[2]);
	strcat(string,mybuf);
    } else if (i_param[0] != 0 && i_param[3] == 1) {
      strcat(string, "cgcent mol 1\n");
    }
    if (i_param[0] != 0 && i_param[4] == 0) {
	sprintf(mybuf, "fgcent %.3f %.3f %.3f\n", 
		fcenter[0], fcenter[1], fcenter[2]);
	strcat(string,mybuf);
    } else if (i_param[0] != 0 && i_param[4] == 1) {
      strcat(string, "fgcent mol 1\n");
    }

    if (i_param[0] == 2) {
      sprintf(mybuf, "pdime %i %i %i\n", pdime[0], pdime[1], pdime[2]);
      strcat(string,mybuf);
      sprintf(mybuf, "ofrac %.3f\n", *ofrac);
      strcat(string,mybuf);
    }

    for (i=0; i < i_param[21]; i++) {
      sprintf(mybuf, "ion charge %.3f conc %.3f radius %.3f\n", 
	      ionq[i], ionc[i], ionr[i]);
      strcat(string,mybuf);
    }
    sprintf(mybuf, "pdie %.3f \nsdie %.3f\nsrad %.3f\nswin %.3f\n\
temp %.3f\nsdens %.3f\nmol 1\n", r_param[0], r_param[1], r_param[2],
	    r_param[3], r_param[4], r_param[5]);
    strcat(string, mybuf);

    switch(i_param[5]){
	case 0: sprintf(mybuf, "chgm spl0\n"); break;
	case 1: sprintf(mybuf, "chgm spl2\n"); break;
	case 2: sprintf(mybuf, "chgm spl4\n"); break;
	default: printf("iAPBS: chgm keyword error\n"); exit(2);
    }
    strcat(string, mybuf);
    switch(i_param[6]){
	case 0: sprintf(mybuf, "lpbe\n"); break;
	case 1: sprintf(mybuf, "npbe\n"); break;
	case 2: sprintf(mybuf, "lrpbe\n"); break;
	case 3: sprintf(mybuf, "nrpbe\n"); break;
        case 4: sprintf(mybuf, "smpbe vol %.3f size %.3f\n", 
			r_param[7], r_param[8]); break;
	default: printf("iAPBS: PBE keyword error\n"); exit(2);
    }
    strcat(string, mybuf);
    switch(i_param[7]){
	case 0: sprintf(mybuf, "bcfl zero\n"); break;
	case 1: sprintf(mybuf, "bcfl sdh\n"); break;
	case 2: sprintf(mybuf, "bcfl mdh\n"); break;
	case 4: sprintf(mybuf, "bcfl focus\n"); break;
	default: printf("iAPBS: bcfl keyword error\n"); exit(2);
    }
    strcat(string, mybuf);
    switch(i_param[8]){
	case 0: sprintf(mybuf, "srfm mol\n"); break;
	case 1: sprintf(mybuf, "srfm smol\n"); break;
	case 2: sprintf(mybuf, "srfm spl2\n"); break;
	case 3: sprintf(mybuf, "srfm spl4\n"); break;
	default: printf("iAPBS: srfm keyword error\n"); exit(2);
    }
    strcat(string, mybuf);
    switch(i_param[9]){
	case 0: sprintf(mybuf, "calcenergy no\n"); break;
	case 1: sprintf(mybuf, "calcenergy total\n"); break;
	case 2: sprintf(mybuf, "calcenergy comps\n"); break;
	default: printf("iAPBS: calcenergy keyword error\n"); exit(2);
    }
    strcat(string, mybuf);
    switch(i_param[10]){
        case 0: sprintf(mybuf, "calcforce no\n"); break;
	case 1: sprintf(mybuf, "calcforce total\n"); break;
	case 2: sprintf(mybuf, "calcforce comps\n"); break;
	default: printf("iAPBS: calcforce keyword error\n"); exit(2);
    }
    strcat(string, mybuf);

    if(i_param[11] == 1){
//      strcat(string, "write pot dx iapbs-pot\n");
      sprintf(string, "%swrite pot dx iapbs-pot-%d\n", string, getpid());
    }
    if(i_param[12] == 1){
      strcat(string, "write charge dx iapbs-charge\n");
    }
    if(i_param[13] == 1){
      strcat(string, "write smol dx iapbs-smol\n");
    }
    if(i_param[14] == 1){
      strcat(string, "write kappa dx iapbs-kappa\n");
    }
    if(i_param[15] == 1){
      strcat(string, "write dielx dx iapbs-dielx\n");
      strcat(string, "write diely dx iapbs-diely\n");
      strcat(string, "write dielz dx iapbs-dielz\n");
    }
    if(i_param[16] == 1){
      strcat(string, "write atompot dx iapbs-atompot\n");
    }


    if(i_param[22] == 1){
      strcat(string, "usemap charge 1\n");
    }
    if(i_param[23] == 1){
      strcat(string, "usemap kappa 1\n");
    }
    if(i_param[24] == 1){
      strcat(string, "usemap diel 1\n");
    }
    if(i_param[17] == 1){
      strcat(string, "usemap pot 1\n");
    }
    strcat(string, "end\n\n");

    // apolar section

    if(i_param[20] > 0){
      strcat(string, "apolar name npolar\n");
      strcat(string, "bconc 0.0\npress 0.0\ndpos 0.2\n");
      strcat(string, "mol 1\nsrfm sacc\n");
      strcat(string, "grid 0.2 0.2 0.2\n");
      sprintf(mybuf, "gamma %.3f\n", r_param[6]);
      strcat(string, mybuf);
      sprintf(mybuf, "srad %.3f\nswin %.3f\ntemp %.3f\nsdens %.3f\n", \
	      r_param[2], r_param[3], r_param[4], r_param[5]);
      strcat(string, mybuf);

      switch(i_param[20]){
      case 1: sprintf(mybuf, "calcenergy total\n"); break;
      case 2: sprintf(mybuf, "calcenergy comps\n"); break;
      default: printf("iAPBS: Unknown calcnpenergy option!\n"); exit(2);
      }
      strcat(string, mybuf);

      switch(i_param[19]){
      case 0: sprintf(mybuf, "calcforce no\n"); break;
      case 1: sprintf(mybuf, "calcforce total\n"); break;
      case 2: sprintf(mybuf, "calcforce comps\n"); break;
      default: printf("iAPBS: Unknown calcnpforce option!\n"); exit(2);
      }
      strcat(string, mybuf);

      strcat(string, "end\n\n");
    }

    // print section
    if (debug > 3) {
      strcat(string, "print elecEnergy elec end\n");
      if(i_param[20] > 0 ){
	strcat(string, "print apolEnergy npolar end\n");
      }
    }
    if (debug > 6) {
      if(i_param[10] > 0 ){
	strcat(string, "print elecForce elec end\n");
      }
      if(i_param[20] > 0  && i_param[19] > 0 ) {
	strcat(string, "print apolForce npolar end\n");
      }
    }

    strcat(string, "\nquit\n");
    //printf("%s", string);
    return string;
}

/**
 * @note Collect Elect Energy
 *
 *
 * @author Robert Konecny
 *
 */
//VPUBLIC double getElecEnergy(Vcom *com, NOsh *nosh, 
double getElecEnergy(Vcom *com, NOsh *nosh, 
	double totEnergy[NOSH_MAXCALC], int iprint)
{

    int iarg, calcid;
    double ltenergy, gtenergy, scalar;
    
    calcid = nosh->elec2calc[nosh->printcalc[iprint][0]];
    if (nosh->calc[calcid]->pbeparm->calcenergy != PCE_NO) {
	ltenergy = Vunit_kb * (1e-3) * Vunit_Na *
	    nosh->calc[calcid]->pbeparm->temp * totEnergy[calcid];
    } else {
	Vnm_tprint( 2, "  Didn't calculate energy in Calculation \
#%d\n", calcid+1);
	    return 0;
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
	calcid = nosh->elec2calc[nosh->printcalc[iprint][iarg]];
	/* Add or subtract? */
	if (nosh->printop[iprint][iarg-1] == 0) scalar = 1.0;
	else if (nosh->printop[iprint][iarg-1] == 1) scalar = -1.0;
	/* Accumulate */
	ltenergy += (scalar * Vunit_kb * (1e-3) * Vunit_Na *
		nosh->calc[calcid]->pbeparm->temp * totEnergy[calcid]);
    }

//    Vnm_tprint( 1, "  Local net energy (PE %d) = %1.12E kJ/mol\n",
//	    Vcom_rank(com), ltenergy);
//    Vnm_tprint( 0, "printEnergy:  Performing global reduction (sum)\n");
    Vcom_reduce(com, &ltenergy, &gtenergy, 1, 2, 0);
//    Vnm_tprint( 1, "  Global net ELEC energy = %1.12E kJ/mol\n", gtenergy);
    //printf("esenF %f\n", gtenergy);
    return gtenergy;
}

