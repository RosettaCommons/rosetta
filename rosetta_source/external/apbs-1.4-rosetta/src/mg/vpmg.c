/**
 *  @file    vpmg.c
 *  @author  Nathan Baker
 *  @brief   Class Vpmg methods
 *  @ingroup Vpmg
 *  @version $Id$
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 *  Nathan A. Baker (nathan.baker@pnnl.gov)
 *  Pacific Northwest National Laboratory
 *
 *  Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2012 Battelle Memorial Institute. Developed at the
 * Pacific Northwest National Laboratory, operated by Battelle Memorial
 * Institute, Pacific Northwest Division for the U.S. Department of Energy.
 *
 * Portions Copyright (c) 2002-2010, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2010, Nathan A. Baker.
 * Portions Copyright (c) 1999-2002, The Regents of the University of
 * California.
 * Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * Neither the name of the developer nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#include "vpmg.h"

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_VPMG)

VPUBLIC unsigned long int Vpmg_memChk(Vpmg *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

#endif /* if !defined(VINLINE_VPMG) */


VPUBLIC void Vpmg_printColComp(Vpmg *thee, char path[72], char title[72],
  char mxtype[3], int flag) {

    int nn, nxm2, nym2, nzm2, ncol, nrow, nonz;
    double *nzval;
    int *colptr, *rowind;

    /* Calculate the total number of unknowns */
    nxm2 = thee->pmgp->nx - 2;
    nym2 = thee->pmgp->ny - 2;
    nzm2 = thee->pmgp->nz - 2;
    nn = nxm2*nym2*nzm2;
    ncol = nn;
    nrow = nn;

    /* Calculate the number of non-zero matrix entries:
     *    nn       nonzeros on diagonal
     *    nn-1     nonzeros on first off-diagonal
     *    nn-nx    nonzeros on second off-diagonal
     *    nn-nx*ny nonzeros on third off-diagonal
     *
     *    7*nn-2*nx*ny-2*nx-2 TOTAL non-zeros
     */
    nonz = 7*nn - 2*nxm2*nym2 - 2*nxm2 - 2;
    nzval  = (double*)Vmem_malloc(thee->vmem, nonz, sizeof(double));
    rowind = (int*)Vmem_malloc(thee->vmem, nonz, sizeof(int));
    colptr = (int*)Vmem_malloc(thee->vmem, (ncol+1), sizeof(int));

#ifndef VAPBSQUIET
    Vnm_print(1, "Vpmg_printColComp:  Allocated space for %d nonzeros\n",
      nonz);
#endif

    bcolcomp(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
      nzval, rowind, colptr, &flag);


#if 0
    for (i=0; i<nn; i++) {
        Vnm_print(1, "nnz(%d) = %g\n", i, nzval[i]);
    }
#endif

    /* I do not understand why I need to pass nzval in this way, but it
     * works... */
    pcolcomp(&nrow, &ncol, &nonz, &(nzval[0]), rowind, colptr, path, title,
      mxtype);

    Vmem_free(thee->vmem, (ncol+1), sizeof(int), (void **)&colptr);
    Vmem_free(thee->vmem, nonz, sizeof(int), (void **)&rowind);
    Vmem_free(thee->vmem, nonz, sizeof(double), (void **)&nzval);

}

VPUBLIC Vpmg* Vpmg_ctor(Vpmgp *pmgp, Vpbe *pbe, int focusFlag,
        Vpmg *pmgOLD, MGparm *mgparm, PBEparm_calcEnergy energyFlag) {

    Vpmg *thee = VNULL;

    thee = (Vpmg*)Vmem_malloc(VNULL, 1, sizeof(Vpmg) );
    VASSERT(thee != VNULL);
    VASSERT( Vpmg_ctor2(thee, pmgp, pbe, focusFlag, pmgOLD, mgparm,
                energyFlag) );
    return thee;
}

VPUBLIC int Vpmg_ctor2(Vpmg *thee, Vpmgp *pmgp, Vpbe *pbe, int focusFlag,
                       Vpmg *pmgOLD, MGparm *mgparm, PBEparm_calcEnergy energyFlag) {

    int i, j, nion;
    double ionConc[MAXION], ionQ[MAXION], ionRadii[MAXION], zkappa2, zks2;
    double ionstr, partMin[3], partMax[3];

    /* Get the parameters */
    VASSERT(pmgp != VNULL);
    VASSERT(pbe != VNULL);
    thee->pmgp = pmgp;
    thee->pbe = pbe;

    /* Set up the memory */
    thee->vmem = Vmem_ctor("APBS:VPMG");



    /// @note  this is common to both replace/noreplace options
    /* Initialize ion concentrations and valencies in PMG routines */
    zkappa2 = Vpbe_getZkappa2(thee->pbe);
    ionstr = Vpbe_getBulkIonicStrength(thee->pbe);
    if (ionstr > 0.0) zks2 = 0.5/ionstr;
    else zks2 = 0.0;
    Vpbe_getIons(thee->pbe, &nion, ionConc, ionRadii, ionQ);

    /* TEMPORARY USEAQUA */
        /* Calculate storage requirements */
    if(mgparm->useAqua == 0){
            Vpmgp_size(thee->pmgp);
    }else{
        VABORT_MSG0("Aqua is currently disabled");
    }

    /* We need some additional storage if: nonlinear & newton OR cgmg */
    /* SMPBE Added - nonlin = 2 added since it mimics NPBE */
    if ( ( ((thee->pmgp->nonlin == NONLIN_NPBE) || (thee->pmgp->nonlin == NONLIN_SMPBE))
           && (thee->pmgp->meth == VSOL_Newton) ) || (thee->pmgp->meth == VSOL_CGMG) )
    {
        thee->pmgp->nrwk += (2*(thee->pmgp->nf));
    }


        if (thee->pmgp->iinfo > 1) {
            Vnm_print(2, "Vpmg_ctor2:  PMG chose nx = %d, ny = %d, nz = %d\n",
                            thee->pmgp->nx, thee->pmgp->ny, thee->pmgp->nz);
            Vnm_print(2, "Vpmg_ctor2:  PMG chose nlev = %d\n",
                            thee->pmgp->nlev);
            Vnm_print(2, "Vpmg_ctor2:  PMG chose nxc = %d, nyc = %d, nzc = %d\n",
                            thee->pmgp->nxc, thee->pmgp->nyc, thee->pmgp->nzc);
            Vnm_print(2, "Vpmg_ctor2:  PMG chose nf = %d, nc = %d\n",
                            thee->pmgp->nf, thee->pmgp->nc);
            Vnm_print(2, "Vpmg_ctor2:  PMG chose narr = %d, narrc = %d\n",
                            thee->pmgp->narr, thee->pmgp->narrc);
            Vnm_print(2, "Vpmg_ctor2:  PMG chose n_rpc = %d, n_iz = %d, n_ipc = %d\n",
                            thee->pmgp->n_rpc, thee->pmgp->n_iz, thee->pmgp->n_ipc);
            Vnm_print(2, "Vpmg_ctor2:  PMG chose nrwk = %d, niwk = %d\n",
                            thee->pmgp->nrwk, thee->pmgp->niwk);
        }



    /* Allocate boundary storage */
    thee->gxcf = (double *)Vmem_malloc(
        thee->vmem,
        10*(thee->pmgp->ny)*(thee->pmgp->nz),
        sizeof(double)
        );

    thee->gycf = (double *)Vmem_malloc(
        thee->vmem,
        10*(thee->pmgp->nx)*(thee->pmgp->nz),
        sizeof(double)
        );

    thee->gzcf = (double *)Vmem_malloc(
        thee->vmem,
        10*(thee->pmgp->nx)*(thee->pmgp->ny),
        sizeof(double)
        );



    /* Warn users if they are using BCFL_MAP that
       we do not include external energies */
    if (thee->pmgp->bcfl == BCFL_MAP)
        Vnm_print(2,"Vpmg_ctor2: \nWarning: External energies are not used in BCFL_MAP calculations!\n");

    if (focusFlag) {

        /* Overwrite any default or user-specified boundary condition
        * arguments; we are now committed to a calculation via focusing */
        if (thee->pmgp->bcfl != BCFL_FOCUS) {
            Vnm_print(2,
                      "Vpmg_ctor2: reset boundary condition flag to BCFL_FOCUS!\n");
            thee->pmgp->bcfl = BCFL_FOCUS;
        }

        /* Fill boundaries */
        Vnm_print(0, "Vpmg_ctor2:  Filling boundary with old solution!\n");
        focusFillBound(thee, pmgOLD);

        /* Calculate energetic contributions from region outside focusing
            * domain */
        if (energyFlag != PCE_NO) {

            if (mgparm->type == MCT_PARALLEL) {

                for (j=0; j<3; j++) {
                    partMin[j] = mgparm->partDisjCenter[j]
                    - 0.5*mgparm->partDisjLength[j];
                    partMax[j] = mgparm->partDisjCenter[j]
                        + 0.5*mgparm->partDisjLength[j];
                }

            } else {
                for (j=0; j<3; j++) {
                    partMin[j] = mgparm->center[j] - 0.5*mgparm->glen[j];
                    partMax[j] = mgparm->center[j] + 0.5*mgparm->glen[j];
                }
            }
            extEnergy(thee, pmgOLD, energyFlag, partMin, partMax,
                      mgparm->partDisjOwnSide);
        }

    } else {

        /* Ignore external energy contributions */
        thee->extQmEnergy = 0;
        thee->extDiEnergy = 0;
        thee->extQfEnergy = 0;
    }

    /* Allocate partition vector storage */
    thee->pvec = (double *)Vmem_malloc(
        thee->vmem,
        (thee->pmgp->nx)*(thee->pmgp->ny)*(thee->pmgp->nz),
        sizeof(double)
        );

    /* Allocate remaining storage */
    thee->iparm  = (   int *)Vmem_malloc(thee->vmem,                100, sizeof(   int));
    thee->rparm  = (double *)Vmem_malloc(thee->vmem,                100, sizeof(double));
    thee->iwork  = (   int *)Vmem_malloc(thee->vmem,   thee->pmgp->niwk, sizeof(   int));
    thee->rwork  = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->nrwk, sizeof(double));
    thee->charge = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->kappa  = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->pot    = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->epsx   = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->epsy   = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->epsz   = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->a1cf   = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->a2cf   = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->a3cf   = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->ccf    = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->fcf    = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->tcf    = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->u      = (double *)Vmem_malloc(thee->vmem,   thee->pmgp->narr, sizeof(double));
    thee->xf     = (double *)Vmem_malloc(thee->vmem, 5*(thee->pmgp->nx), sizeof(double));
    thee->yf     = (double *)Vmem_malloc(thee->vmem, 5*(thee->pmgp->ny), sizeof(double));
    thee->zf     = (double *)Vmem_malloc(thee->vmem, 5*(thee->pmgp->nz), sizeof(double));



    /* Packs parameters into the iparm and rparm arrays */
    Vpackmg(thee->iparm, thee->rparm, &(thee->pmgp->nrwk), &(thee->pmgp->niwk),
            &(thee->pmgp->nx), &(thee->pmgp->ny), &(thee->pmgp->nz),
            &(thee->pmgp->nlev), &(thee->pmgp->nu1), &(thee->pmgp->nu2),
            &(thee->pmgp->mgkey), &(thee->pmgp->itmax), &(thee->pmgp->istop),
            &(thee->pmgp->ipcon), &(thee->pmgp->nonlin), &(thee->pmgp->mgsmoo),
            &(thee->pmgp->mgprol), &(thee->pmgp->mgcoar), &(thee->pmgp->mgsolv),
            &(thee->pmgp->mgdisc), &(thee->pmgp->iinfo), &(thee->pmgp->errtol),
            &(thee->pmgp->ipkey), &(thee->pmgp->omegal), &(thee->pmgp->omegan),
            &(thee->pmgp->irite), &(thee->pmgp->iperf));



    /* Currently for SMPBE type calculations we do not want to apply a scale
        factor to the ionConc */
    /** @note  The fortran replacement functions are run along side the old
     *         fortran functions.  This is due to the use of common variables
     *         in the fortran sub-routines.  Once the fortran code has been
     *         successfully excised, these functions will no longer need to be
     *         called in tandem, and the fortran version may be dropped
     */
    switch(pmgp->ipkey){

        case IPKEY_SMPBE:

                        Vmypdefinitsmpbe(&nion, ionQ, ionConc, &pbe->smvolume, &pbe->smsize);
            break;



        case IPKEY_NPBE:

            /* Else adjust the inoConc by scaling factor zks2 */
            for (i=0; i<nion; i++)
                ionConc[i] = zks2 * ionConc[i];

                        Vmypdefinitnpbe(&nion, ionQ, ionConc);
            break;



        case IPKEY_LPBE:

            /* Else adjust the inoConc by scaling factor zks2 */
            for (i=0; i<nion; i++)
                            ionConc[i] = zks2 * ionConc[i];

            Vmypdefinitlpbe(&nion, ionQ, ionConc);
            break;



        default:
            Vnm_print(2, "PMG: Warning: PBE structure not initialized!\n");
            /* Else adjust the inoConc by scaling factor zks2 */
            for (i=0; i<nion; i++)
                ionConc[i] = zks2 * ionConc[i];
            break;
    }

    /* Set the default chargeSrc for 5th order splines */
    thee->chargeSrc = mgparm->chgs;

    /* Turn off restriction of observable calculations to a specific
    * partition */
    Vpmg_unsetPart(thee);

    /* The coefficient arrays have not been filled */
    thee->filled = 0;


    /*
     * TODO: Move the dtor out of here. The current ctor is done in routines.c,
     *       This was originally moved out to kill a memory leak. The dtor has
     *       has been removed from initMG and placed back here to keep memory
     *       usage low. killMG has been modified accordingly.
     */
    Vpmg_dtor(&pmgOLD);

    return 1;
}

VPUBLIC int Vpmg_solve(Vpmg *thee) {

    int i,
        nx,
        ny,
        nz,
        n;
    double zkappa2;

    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    n = nx*ny*nz;

    if (!(thee->filled)) {
        Vnm_print(2, "Vpmg_solve:  Need to call Vpmg_fillco()!\n");
        return 0;
    }

    /* Fill the "true solution" array */
    for (i=0; i<n; i++) {
        thee->tcf[i] = 0.0;
    }

    /* Fill the RHS array */
    for (i=0; i<n; i++) {
        thee->fcf[i] = thee->charge[i];
    }

    /* Fill the operator coefficient array. */
    for (i=0; i<n; i++) {
        thee->a1cf[i] = thee->epsx[i];
        thee->a2cf[i] = thee->epsy[i];
        thee->a3cf[i] = thee->epsz[i];
    }

    /* Fill the nonlinear coefficient array by multiplying the kappa
     * accessibility array (containing values between 0 and 1) by zkappa2. */
    zkappa2 = Vpbe_getZkappa2(thee->pbe);
    if (zkappa2 > VPMGSMALL) {
        for (i=0; i<n; i++) {
            thee->ccf[i] = zkappa2*thee->kappa[i];
        }
    } else {
        for (i=0; i<n; i++) {
            thee->ccf[i] = 0.0;
        }
    }

    switch(thee->pmgp->meth) {
        /* CGMG (linear) */
        case VSOL_CGMG:

            if (thee->pmgp->iinfo > 1)
                Vnm_print(2, "Driving with CGMGDRIV\n");

            VABORT_MSG0("CGMGDRIV is not currently supported");
            break;

        /* Newton (nonlinear) */
        case VSOL_Newton:

            if (thee->pmgp->iinfo > 1)
                Vnm_print(2, "Driving with NEWDRIV\n");

            Vnewdriv
                      (thee->iparm, thee->rparm, thee->iwork, thee->rwork,
                       thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
                       thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
                       thee->fcf, thee->tcf);
            break;

        /* MG (linear/nonlinear) */
        case VSOL_MG:

            if (thee->pmgp->iinfo > 1)
                Vnm_print(2, "Driving with MGDRIV\n");

            Vmgdriv(thee->iparm, thee->rparm, thee->iwork, thee->rwork,
                                        thee->u, thee->xf, thee->yf, thee->zf, thee->gxcf, thee->gycf,
                                        thee->gzcf, thee->a1cf, thee->a2cf, thee->a3cf, thee->ccf,
                                        thee->fcf, thee->tcf);
            break;

        /* CGHS (linear/nonlinear) */
        case VSOL_CG:

            if (thee->pmgp->iinfo > 1)
                Vnm_print(2, "Driving with NCGHSDRIV\n");

            VABORT_MSG0("NCGHSDRIV is not currently supported");
            break;

        /* SOR (linear/nonlinear) */
        case VSOL_SOR:

            if (thee->pmgp->iinfo > 1)
                Vnm_print(2, "Driving with NSORDRIV\n");

            VABORT_MSG0("NSORDRIV is not currently supported");
            break;

        /* GSRB (linear/nonlinear) */
        case VSOL_RBGS:

            if (thee->pmgp->iinfo > 1)
                Vnm_print(2, "Driving with NGSRBDRIV\n");

            VABORT_MSG0("NGSRBDRIV is not currently supported");
            break;

        /* WJAC (linear/nonlinear) */
        case VSOL_WJ:

            if (thee->pmgp->iinfo > 1)
                Vnm_print(2, "Driving with NWJACDRIV\n");

            VABORT_MSG0("NWJACDRIV is not currently supported");
            break;

        /* RICH (linear/nonlinear) */
        case VSOL_Richardson:

            if (thee->pmgp->iinfo > 1)
                Vnm_print(2, "Driving with NRICHDRIV\n");

            VABORT_MSG0("NRICHDRIV is not currently supported");
            break;

        /* CGMG (linear) TEMPORARY USEAQUA */
        case VSOL_CGMGAqua:

            if (thee->pmgp->iinfo > 1)
                Vnm_print(2, "Driving with CGMGDRIVAQUA\n");

            VABORT_MSG0("CGMGDRIVAQUA is not currently supported");
            break;

        /* Newton (nonlinear) TEMPORARY USEAQUA */
        case VSOL_NewtonAqua:

            if (thee->pmgp->iinfo > 1)
                Vnm_print(2, "Driving with NEWDRIVAQUA\n");

            VABORT_MSG0("NEWDRIVAQUA is not currently supported");
            break;

        /* Error handling */
        default:
            Vnm_print(2, "Vpmg_solve: invalid solver method key (%d)\n",
              thee->pmgp->key);
            return 0;
            break;
    }

    return 1;

}


VPUBLIC void Vpmg_dtor(Vpmg **thee) {

    if ((*thee) != VNULL) {
        Vpmg_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vpmg), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void Vpmg_dtor2(Vpmg *thee) {

    /* Clean up the storage */

    Vmem_free(thee->vmem,              100, sizeof(int),
      (void **)&(thee->iparm));
    Vmem_free(thee->vmem,              100, sizeof(double),
        (void **)&(thee->rparm));
    Vmem_free(thee->vmem, thee->pmgp->niwk, sizeof(int),
      (void **)&(thee->iwork));
    Vmem_free(thee->vmem, thee->pmgp->nrwk, sizeof(double),
      (void **)&(thee->rwork));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->charge));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->kappa));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
              (void **)&(thee->pot));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->epsx));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->epsy));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->epsz));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->a1cf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->a2cf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->a3cf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->ccf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->fcf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->tcf));
    Vmem_free(thee->vmem, thee->pmgp->narr, sizeof(double),
      (void **)&(thee->u));
    Vmem_free(thee->vmem, 5*(thee->pmgp->nx), sizeof(double),
      (void **)&(thee->xf));
    Vmem_free(thee->vmem, 5*(thee->pmgp->ny), sizeof(double),
      (void **)&(thee->yf));
    Vmem_free(thee->vmem, 5*(thee->pmgp->nz), sizeof(double),
      (void **)&(thee->zf));
    Vmem_free(thee->vmem, 10*(thee->pmgp->ny)*(thee->pmgp->nz), sizeof(double),
      (void **)&(thee->gxcf));
    Vmem_free(thee->vmem, 10*(thee->pmgp->nx)*(thee->pmgp->nz), sizeof(double),
      (void **)&(thee->gycf));
    Vmem_free(thee->vmem, 10*(thee->pmgp->nx)*(thee->pmgp->ny), sizeof(double),
      (void **)&(thee->gzcf));
    Vmem_free(thee->vmem, (thee->pmgp->nx)*(thee->pmgp->ny)*(thee->pmgp->nz),
      sizeof(double), (void **)&(thee->pvec));

    Vmem_dtor(&(thee->vmem));
}

VPUBLIC void Vpmg_setPart(Vpmg *thee, double lowerCorner[3],
        double upperCorner[3], int bflags[6]) {

    Valist *alist;
    Vatom *atom;
    int i, j, k, nx, ny, nz;
    double xmin, ymin, zmin, x, y, z, hx, hy, hzed, xok, yok, zok;
    double x0,x1,y0,y1,z0,z1;

    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmin = thee->pmgp->xcent - 0.5*hx*(nx-1);
    ymin = thee->pmgp->ycent - 0.5*hy*(ny-1);
    zmin = thee->pmgp->zcent - 0.5*hzed*(nz-1);

    xok = 0;
    yok = 0;
    zok = 0;

    /* We need have called Vpmg_fillco first */

    alist = thee->pbe->alist;

    Vnm_print(0, "Vpmg_setPart:  lower corner = (%g, %g, %g)\n",
      lowerCorner[0], lowerCorner[1], lowerCorner[2]);
    Vnm_print(0, "Vpmg_setPart:  upper corner = (%g, %g, %g)\n",
      upperCorner[0], upperCorner[1], upperCorner[2]);
    Vnm_print(0, "Vpmg_setPart:  actual minima = (%g, %g, %g)\n",
      xmin, ymin, zmin);
    Vnm_print(0, "Vpmg_setPart:  actual maxima = (%g, %g, %g)\n",
      xmin+hx*(nx-1), ymin+hy*(ny-1), zmin+hzed*(nz-1));
    Vnm_print(0, "Vpmg_setPart:  bflag[FRONT] = %d\n",
      bflags[VAPBS_FRONT]);
    Vnm_print(0, "Vpmg_setPart:  bflag[BACK] = %d\n",
      bflags[VAPBS_BACK]);
    Vnm_print(0, "Vpmg_setPart:  bflag[LEFT] = %d\n",
      bflags[VAPBS_LEFT]);
    Vnm_print(0, "Vpmg_setPart:  bflag[RIGHT] = %d\n",
      bflags[VAPBS_RIGHT]);
    Vnm_print(0, "Vpmg_setPart:  bflag[UP] = %d\n",
      bflags[VAPBS_UP]);
    Vnm_print(0, "Vpmg_setPart:  bflag[DOWN] = %d\n",
      bflags[VAPBS_DOWN]);

    /* Identify atoms as inside, outside, or on the border
       If on the border, use the bflags to determine if there
       is an adjacent processor - if so, this atom should be equally
       shared. */

    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);

        if ((atom->position[0] < upperCorner[0]) &&
            (atom->position[0] > lowerCorner[0])) xok = 1;
        else {
            if ((VABS(atom->position[0] - lowerCorner[0]) < VPMGSMALL) &&
                (bflags[VAPBS_LEFT] == 0)) xok = 1;
            else if ((VABS(atom->position[0] - lowerCorner[0]) < VPMGSMALL) &&
                (bflags[VAPBS_LEFT] == 1)) xok = 0.5;
            else if ((VABS(atom->position[0] - upperCorner[0]) < VPMGSMALL) &&
                (bflags[VAPBS_RIGHT] == 0)) xok = 1;
            else if ((VABS(atom->position[0] - upperCorner[0]) < VPMGSMALL) &&
                (bflags[VAPBS_RIGHT] == 1)) xok = 0.5;
            else xok = 0;
        }
        if ((atom->position[1] < upperCorner[1]) &&
            (atom->position[1] > lowerCorner[1])) yok = 1;
        else {
            if ((VABS(atom->position[1] - lowerCorner[1]) < VPMGSMALL) &&
                (bflags[VAPBS_BACK] == 0)) yok = 1;
            else if ((VABS(atom->position[1] - lowerCorner[1]) < VPMGSMALL) &&
                (bflags[VAPBS_BACK] == 1)) yok = 0.5;
            else if ((VABS(atom->position[1] - upperCorner[1]) < VPMGSMALL) &&
                (bflags[VAPBS_FRONT] == 0)) yok = 1;
            else if ((VABS(atom->position[1] - upperCorner[1]) < VPMGSMALL) &&
                (bflags[VAPBS_FRONT] == 1)) yok = 0.5;
            else yok = 0;
        }
        if ((atom->position[2] < upperCorner[2]) &&
            (atom->position[2] > lowerCorner[2])) zok = 1;
        else {
            if ((VABS(atom->position[2] - lowerCorner[2]) < VPMGSMALL) &&
                (bflags[VAPBS_DOWN] == 0)) zok = 1;
            else if ((VABS(atom->position[2] - lowerCorner[2]) < VPMGSMALL) &&
                (bflags[VAPBS_DOWN] == 1)) zok = 0.5;
            else if ((VABS(atom->position[2] - upperCorner[2]) < VPMGSMALL) &&
                (bflags[VAPBS_UP] == 0)) zok = 1;
            else if ((VABS(atom->position[2] - upperCorner[2]) < VPMGSMALL) &&
                (bflags[VAPBS_UP] == 1)) zok = 0.5;
            else zok = 0;
        }

        atom->partID = xok*yok*zok;
        /*
        Vnm_print(1, "DEBUG (%s, %d):  atom->position[0] - upperCorner[0] = %g\n",
                  __FILE__, __LINE__, atom->position[0] - upperCorner[0]);
        Vnm_print(1, "DEBUG (%s, %d):  atom->position[0] - lowerCorner[0] = %g\n",
                  __FILE__, __LINE__, atom->position[0] - lowerCorner[0]);
        Vnm_print(1, "DEBUG (%s, %d):  atom->position[1] - upperCorner[1] = %g\n",
                  __FILE__, __LINE__, atom->position[1] - upperCorner[1]);
        Vnm_print(1, "DEBUG (%s, %d):  atom->position[1] - lowerCorner[1] = %g\n",
                  __FILE__, __LINE__, atom->position[1] - lowerCorner[1]);
        Vnm_print(1, "DEBUG (%s, %d):  atom->position[2] - upperCorner[2] = %g\n",
                  __FILE__, __LINE__, atom->position[2] - upperCorner[2]);
        Vnm_print(1, "DEBUG (%s, %d):  atom->position[2] - lowerCorner[0] = %g\n",
                  __FILE__, __LINE__, atom->position[2] - lowerCorner[2]);
        Vnm_print(1, "DEBUG (%s, %d):  xok = %g, yok = %g, zok = %g\n",
                  __FILE__, __LINE__, xok, yok, zok);
         */

    }

    /* Load up pvec -
       For all points within h{axis}/2 of a border - use a gradient
       to determine the pvec weight.
       Points on the boundary depend on the presence of an adjacent
       processor. */

    for (i=0; i<(nx*ny*nz); i++) thee->pvec[i] = 0.0;

    for (i=0; i<nx; i++) {
        xok = 0.0;
        x = i*hx + xmin;
        if ( (x < (upperCorner[0]-hx/2)) &&
             (x > (lowerCorner[0]+hx/2))
           ) xok = 1.0;
        else if ( (VABS(x - lowerCorner[0]) < VPMGSMALL) &&
                  (bflags[VAPBS_LEFT] == 0)) xok = 1.0;
        else if ((VABS(x - lowerCorner[0]) < VPMGSMALL) &&
                 (bflags[VAPBS_LEFT] == 1)) xok = 0.5;
        else if ((VABS(x - upperCorner[0]) < VPMGSMALL) &&
                 (bflags[VAPBS_RIGHT] == 0)) xok = 1.0;
        else if ((VABS(x - upperCorner[0]) < VPMGSMALL) &&
                 (bflags[VAPBS_RIGHT] == 1)) xok = 0.5;
        else if ((x > (upperCorner[0] + hx/2)) || (x < (lowerCorner[0] - hx/2))) xok = 0.0;
        else if ((x < (upperCorner[0] + hx/2)) || (x > (lowerCorner[0] - hx/2))) {
            x0 = VMAX2(x - hx/2, lowerCorner[0]);
            x1 = VMIN2(x + hx/2, upperCorner[0]);
            xok = VABS(x1-x0)/hx;

            if (xok < 0.0) {
                if (VABS(xok) < VPMGSMALL) xok = 0.0;
                else {
                    Vnm_print(2, "Vpmg_setPart:  fell off x-interval (%1.12E)!\n",
                            xok);
                    VASSERT(0);
                }
            }
            if (xok > 1.0) {
                if (VABS(xok - 1.0) < VPMGSMALL) xok = 1.0;
                else {
                    Vnm_print(2, "Vpmg_setPart:  fell off x-interval (%1.12E)!\n",
                            xok);
                    VASSERT(0);
                }
            }

        } else xok = 0.0;

        for (j=0; j<ny; j++) {
            yok = 0.0;
            y = j*hy + ymin;
            if ((y < (upperCorner[1]-hy/2)) && (y > (lowerCorner[1]+hy/2))) yok = 1.0;
            else if ((VABS(y - lowerCorner[1]) < VPMGSMALL) &&
                     (bflags[VAPBS_BACK] == 0)) yok = 1.0;
            else if ((VABS(y - lowerCorner[1]) < VPMGSMALL) &&
                     (bflags[VAPBS_BACK] == 1)) yok = 0.5;
            else if ((VABS(y - upperCorner[1]) < VPMGSMALL) &&
                     (bflags[VAPBS_FRONT] == 0)) yok = 1.0;
            else if ((VABS(y - upperCorner[1]) < VPMGSMALL) &&
                     (bflags[VAPBS_FRONT] == 1)) yok = 0.5;
            else if ((y > (upperCorner[1] + hy/2)) || (y < (lowerCorner[1] - hy/2))) yok=0.0;
            else if ((y < (upperCorner[1] + hy/2)) || (y > (lowerCorner[1] - hy/2))){
                y0 = VMAX2(y - hy/2, lowerCorner[1]);
                y1 = VMIN2(y + hy/2, upperCorner[1]);
                yok = VABS(y1-y0)/hy;

                if (yok < 0.0) {
                    if (VABS(yok) < VPMGSMALL) yok = 0.0;
                    else {
                        Vnm_print(2, "Vpmg_setPart:  fell off y-interval (%1.12E)!\n",
                                yok);
                        VASSERT(0);
                    }
                }
                if (yok > 1.0) {
                    if (VABS(yok - 1.0) < VPMGSMALL) yok = 1.0;
                    else {
                        Vnm_print(2, "Vpmg_setPart:  fell off y-interval (%1.12E)!\n",
                                yok);
                        VASSERT(0);
                    }
                }
            }
            else yok=0.0;

            for (k=0; k<nz; k++) {
                zok = 0.0;
                z = k*hzed + zmin;
                if ((z < (upperCorner[2]-hzed/2)) && (z > (lowerCorner[2]+hzed/2))) zok = 1.0;
                else if ((VABS(z - lowerCorner[2]) < VPMGSMALL) &&
                         (bflags[VAPBS_DOWN] == 0)) zok = 1.0;
                else if ((VABS(z - lowerCorner[2]) < VPMGSMALL) &&
                         (bflags[VAPBS_DOWN] == 1)) zok = 0.5;
                else if ((VABS(z - upperCorner[2]) < VPMGSMALL) &&
                         (bflags[VAPBS_UP] == 0)) zok = 1.0;
                else if ((VABS(z - upperCorner[2]) < VPMGSMALL) &&
                         (bflags[VAPBS_UP] == 1)) zok = 0.5;
                else if ((z > (upperCorner[2] + hzed/2)) || (z < (lowerCorner[2] - hzed/2))) zok=0.0;
                else if ((z < (upperCorner[2] + hzed/2)) || (z > (lowerCorner[2] - hzed/2))){
                    z0 = VMAX2(z - hzed/2, lowerCorner[2]);
                    z1 = VMIN2(z + hzed/2, upperCorner[2]);
                    zok = VABS(z1-z0)/hzed;

                    if (zok < 0.0) {
                        if (VABS(zok) < VPMGSMALL) zok = 0.0;
                        else {
                            Vnm_print(2, "Vpmg_setPart:  fell off z-interval (%1.12E)!\n",
                                    zok);
                            VASSERT(0);
                        }
                    }
                    if (zok > 1.0) {
                        if (VABS(zok - 1.0) < VPMGSMALL) zok = 1.0;
                        else {
                            Vnm_print(2, "Vpmg_setPart:  fell off z-interval (%1.12E)!\n",
                                    zok);
                            VASSERT(0);
                        }
                    }
                }
                else zok = 0.0;

                if (VABS(xok*yok*zok) < VPMGSMALL) thee->pvec[IJK(i,j,k)] = 0.0;
                else thee->pvec[IJK(i,j,k)] = xok*yok*zok;

            }
        }
    }
}

VPUBLIC void Vpmg_unsetPart(Vpmg *thee) {

    int i, nx, ny, nz;
    Vatom *atom;
    Valist *alist;

    VASSERT(thee != VNULL);

    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    alist = thee->pbe->alist;

    for (i=0; i<(nx*ny*nz); i++) thee->pvec[i] = 1;
    for (i=0; i<Valist_getNumberAtoms(alist); i++) {
        atom = Valist_getAtom(alist, i);
        atom->partID = 1;
    }
}

VPUBLIC int Vpmg_fillArray(Vpmg *thee, double *vec, Vdata_Type type,
  double parm, Vhal_PBEType pbetype, PBEparm *pbeparm) {

    Vacc *acc = VNULL;
    Vpbe *pbe = VNULL;
    Vgrid *grid = VNULL;
    Vatom *atoms = VNULL;
    Valist *alist = VNULL;
    double position[3], hx, hy, hzed, xmin, ymin, zmin;
    double grad[3], eps, epsp, epss, zmagic;
    int i, j, k, l, nx, ny, nz, ichop;

    pbe = thee->pbe;
    acc = Vpbe_getVacc(pbe);
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    epsp = Vpbe_getSoluteDiel(pbe);
    epss = Vpbe_getSolventDiel(pbe);
    zmagic = Vpbe_getZmagic(pbe);

    if (!(thee->filled)) {
        Vnm_print(2, "Vpmg_fillArray:  need to call Vpmg_fillco first!\n");
        return 0;
    }

    switch (type) {

        case VDT_CHARGE:

            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->charge[i]/zmagic;
            break;

        case VDT_DIELX:

            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->epsx[i];
            break;

        case VDT_DIELY:

            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->epsy[i];
            break;

        case VDT_DIELZ:

            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->epsz[i];
            break;

        case VDT_KAPPA:

            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->kappa[i];
            break;

        case VDT_POT:

            for (i=0; i<nx*ny*nz; i++) vec[i] = thee->u[i];
            break;

        case VDT_ATOMPOT:
            alist = thee->pbe->alist;
            atoms = alist[pbeparm->molid-1].atoms;
            grid = Vgrid_ctor(nx, ny, nz, hx, hy,
                              hzed, xmin, ymin, zmin,thee->u);
            for (i=0; i<alist[pbeparm->molid-1].number;i++) {
                position[0] = atoms[i].position[0];
                position[1] = atoms[i].position[1];
                position[2] = atoms[i].position[2];

                Vgrid_value(grid, position, &vec[i]);
            }
            Vgrid_dtor(&grid);
            break;

        case VDT_SMOL:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;

                        vec[IJK(i,j,k)] = (Vacc_molAcc(acc,position,parm));
                    }
                }
            }
            break;

        case VDT_SSPL:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;

                        vec[IJK(i,j,k)] = Vacc_splineAcc(acc,position,parm,0);
                    }
                }
            }
            break;

        case VDT_VDW:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;

                        vec[IJK(i,j,k)] = Vacc_vdwAcc(acc,position);
                    }
                }
            }
            break;

        case VDT_IVDW:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;

                        vec[IJK(i,j,k)] = Vacc_ivdwAcc(acc,position,parm);
                    }
                }
            }
            break;

        case VDT_LAP:

            grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
              thee->u);
            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        if ((k==0) || (k==(nz-1)) ||
                            (j==0) || (j==(ny-1)) ||
                            (i==0) || (i==(nx-1))) {

                            vec[IJK(i,j,k)] = 0;

                        } else {
                                position[0] = i*hx + xmin;
                                position[1] = j*hy + ymin;
                                position[2] = k*hzed + zmin;
                                VASSERT(Vgrid_curvature(grid,position, 1,
                                  &(vec[IJK(i,j,k)])));
                        }
                    }
                }
            }
            Vgrid_dtor(&grid);
            break;

        case VDT_EDENS:

            grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
              thee->u);
            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;
                        VASSERT(Vgrid_gradient(grid, position, grad));
                        eps = epsp + (epss-epsp)*Vacc_molAcc(acc, position,
                          pbe->solventRadius);
                        vec[IJK(i,j,k)] = 0.0;
                        for (l=0; l<3; l++)
                          vec[IJK(i,j,k)] += eps*VSQR(grad[l]);
                    }
                }
            }
            Vgrid_dtor(&grid);
            break;

        case VDT_NDENS:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;
                        vec[IJK(i,j,k)] = 0.0;
                        if ( VABS(Vacc_ivdwAcc(acc,
                               position, pbe->maxIonRadius) - 1.0) < VSMALL) {
                            for (l=0; l<pbe->numIon; l++) {
                                if (pbetype == PBE_NPBE || pbetype == PBE_SMPBE /*  SMPBE Added */) {
                                    vec[IJK(i,j,k)] += (pbe->ionConc[l]
                                        * Vcap_exp(-pbe->ionQ[l]*thee->u[IJK(i,j,k)],
                                        &ichop));
                                } else if (pbetype == PBE_LPBE){
                                    vec[IJK(i,j,k)] += (pbe->ionConc[l]
                                        * (1 - pbe->ionQ[l]*thee->u[IJK(i,j,k)]));
                                }
                            }
                        }
                    }
                }
            }
            break;

        case VDT_QDENS:

            for (k=0; k<nz; k++) {
                for (j=0; j<ny; j++) {
                    for (i=0; i<nx; i++) {

                        position[0] = i*hx + xmin;
                        position[1] = j*hy + ymin;
                        position[2] = k*hzed + zmin;
                        vec[IJK(i,j,k)] = 0.0;
                        if ( VABS(Vacc_ivdwAcc(acc,
                               position, pbe->maxIonRadius) - 1.0) < VSMALL) {
                            for (l=0; l<pbe->numIon; l++) {
                                if (pbetype == PBE_NPBE || pbetype == PBE_SMPBE /*  SMPBE Added */) {
                                    vec[IJK(i,j,k)] += (pbe->ionConc[l]
                                        * pbe->ionQ[l]
                                        * Vcap_exp(-pbe->ionQ[l]*thee->u[IJK(i,j,k)],
                                        &ichop));
                                } else if (pbetype == PBE_LPBE) {
                                    vec[IJK(i,j,k)] += (pbe->ionConc[l]
                                        * pbe->ionQ[l]
                                        * (1 - pbe->ionQ[l]*thee->u[IJK(i,j,k)]));
                                }
                            }
                        }
                    }
                }
            }
            break;

        default:

            Vnm_print(2, "main:  Bogus data type (%d)!\n", type);
            return 0;
            break;

    }

    return 1;

}

VPRIVATE double Vpmg_polarizEnergy(Vpmg *thee,
                                   int extFlag
                                  ) {

    int i,
        j,
        k,
        ijk,
        nx,
        ny,
        nz,
        iatom;
    double xmin,
           ymin,
           zmin,
           //x, // gcc: not used
           //y,
           //z,
           hx,
           hy,
           hzed,
           epsp,
           lap,
           pt[3],
           T,
           pre,
           polq,
           dist2,
           dist,
           energy,
           q,
           *charge,
           *pos,
           eps_w;
    Vgrid *potgrid;
    Vpbe *pbe;
    Valist *alist;
    Vatom *atom;

    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->ymin;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    pbe = thee->pbe;
    epsp = Vpbe_getSoluteDiel(pbe);
    eps_w = Vpbe_getSolventDiel(pbe);
    alist = pbe->alist;
    charge = thee->charge;

    /* Calculate the prefactor for Coulombic calculations */
    T = Vpbe_getTemperature(pbe);
    pre = (Vunit_ec*Vunit_ec)/(4*VPI*Vunit_eps0*eps_w*Vunit_kb*T);
    pre = pre*(1.0e10);

    /* Set up Vgrid object with solution */
    potgrid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin, thee->u);

    /* Calculate polarization charge */
    energy = 0.0;
    for (i=1; i<(nx-1); i++) {
        pt[0] = xmin + hx*i;
        for (j=1; j<(ny-1); j++) {
            pt[1] = ymin + hy*j;
            for (k=1; k<(nz-1); k++) {
                pt[2] = zmin + hzed*k;

                /* Calculate polarization charge */
                VASSERT(Vgrid_curvature(potgrid, pt, 1, &lap));
                ijk = IJK(i,j,k);
                polq = charge[ijk] + epsp*lap*3.0;

                /* Calculate interaction energy with atoms */
                if (VABS(polq) > VSMALL) {
                    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
                        atom = Valist_getAtom(alist, iatom);
                        q = Vatom_getCharge(atom);
                        pos = Vatom_getPosition(atom);
                        dist2 = VSQR(pos[0]-pt[0]) + VSQR(pos[1]-pt[1]) \
                                + VSQR(pos[2]-pt[2]);
                        dist = VSQRT(dist2);

                        if (dist < VSMALL) {
                            Vnm_print(2, "Vpmg_polarizEnergy:  atom on grid point; ignoring!\n");
                        } else {
                            energy = energy + polq*q/dist;
                        }
                    }
                }
            }
        }
    }

    return pre*energy;
}

VPUBLIC double Vpmg_energy(Vpmg *thee,
                           int extFlag
                          ) {

    double totEnergy = 0.0,
           dielEnergy = 0.0,
           qmEnergy = 0.0,
           qfEnergy = 0.0;

    VASSERT(thee != VNULL);

    if ((thee->pmgp->nonlin) && (Vpbe_getBulkIonicStrength(thee->pbe) > 0.)) {
        Vnm_print(0, "Vpmg_energy:  calculating full PBE energy\n");
        qmEnergy = Vpmg_qmEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  qmEnergy = %1.12E kT\n", qmEnergy);
        qfEnergy = Vpmg_qfEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  qfEnergy = %1.12E kT\n", qfEnergy);
        dielEnergy = Vpmg_dielEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  dielEnergy = %1.12E kT\n", dielEnergy);
        totEnergy = qfEnergy - dielEnergy - qmEnergy;
    } else {
        Vnm_print(0, "Vpmg_energy:  calculating only q-phi energy\n");
        qfEnergy = Vpmg_qfEnergy(thee, extFlag);
        Vnm_print(0, "Vpmg_energy:  qfEnergy = %1.12E kT\n", qfEnergy);
        totEnergy = 0.5*qfEnergy;
    }

    return totEnergy;

}

VPUBLIC double Vpmg_dielEnergy(Vpmg *thee,
                               int extFlag
                              ) {

    double hx,
           hy,
           hzed,
           energy,
           nrgx,
           nrgy,
           nrgz,
           pvecx,
           pvecy,
           pvecz;
    int i,
        j,
        k,
        nx,
        ny,
        nz;

    VASSERT(thee != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    energy = 0.0;

    if (!thee->filled) {
        Vnm_print(2, "Vpmg_dielEnergy:  Need to call Vpmg_fillco!\n");
        VASSERT(0);
    }

    for (k=0; k<(nz-1); k++) {
        for (j=0; j<(ny-1); j++) {
            for (i=0; i<(nx-1); i++) {
                pvecx = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i+1,j,k)]);
                pvecy = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i,j+1,k)]);
                pvecz = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i,j,k+1)]);
                nrgx = thee->epsx[IJK(i,j,k)]*pvecx
                  * VSQR((thee->u[IJK(i,j,k)]-thee->u[IJK(i+1,j,k)])/hx);
                nrgy = thee->epsy[IJK(i,j,k)]*pvecy
                  * VSQR((thee->u[IJK(i,j,k)]-thee->u[IJK(i,j+1,k)])/hy);
                nrgz = thee->epsz[IJK(i,j,k)]*pvecz
                  * VSQR((thee->u[IJK(i,j,k)]-thee->u[IJK(i,j,k+1)])/hzed);
                energy += (nrgx + nrgy + nrgz);
            }
        }
    }

    energy = 0.5*energy*hx*hy*hzed;
    energy = energy/Vpbe_getZmagic(thee->pbe);

    if (extFlag == 1) energy += (thee->extDiEnergy);

    return energy;
}

VPUBLIC double Vpmg_dielGradNorm(Vpmg *thee) {

    double hx, hy, hzed, energy, nrgx, nrgy, nrgz, pvecx, pvecy, pvecz;
    int i, j, k, nx, ny, nz;

    VASSERT(thee != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    energy = 0.0;

    if (!thee->filled) {
        Vnm_print(2, "Vpmg_dielGradNorm:  Need to call Vpmg_fillco!\n");
        VASSERT(0);
    }

    for (k=1; k<nz; k++) {
        for (j=1; j<ny; j++) {
            for (i=1; i<nx; i++) {
                pvecx = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i-1,j,k)]);
                pvecy = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i,j-1,k)]);
                pvecz = 0.5*(thee->pvec[IJK(i,j,k)]+thee->pvec[IJK(i,j,k-1)]);
                nrgx = pvecx
                 * VSQR((thee->epsx[IJK(i,j,k)]-thee->epsx[IJK(i-1,j,k)])/hx);
                nrgy = pvecy
                 * VSQR((thee->epsy[IJK(i,j,k)]-thee->epsy[IJK(i,j-1,k)])/hy);
                nrgz = pvecz
                 * VSQR((thee->epsz[IJK(i,j,k)]-thee->epsz[IJK(i,j,k-1)])/hzed);
                energy += VSQRT(nrgx + nrgy + nrgz);
            }
        }
    }

    energy = energy*hx*hy*hzed;

    return energy;
}

VPUBLIC double Vpmg_qmEnergy(Vpmg *thee,
                             int extFlag
                            ) {

    double energy;

    if(thee->pbe->ipkey == IPKEY_SMPBE){
        energy = Vpmg_qmEnergySMPBE(thee,extFlag);
    }else{
        energy = Vpmg_qmEnergyNONLIN(thee,extFlag);
    }

    return energy;
}

VPRIVATE double Vpmg_qmEnergyNONLIN(Vpmg *thee,
                                    int extFlag
                                   ) {

    double hx,
           hy,
           hzed,
           energy,
           ionConc[MAXION],
           ionRadii[MAXION],
           ionQ[MAXION],
           zkappa2,
           ionstr,
           zks2;
    int i, /* Loop variable */
        j,
        nx,
        ny,
        nz,
        nion,
        ichop,
        nchop,
        len; /* Stores number of iterations for loops to avoid multiple recalculations */

    VASSERT(thee != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    zkappa2 = Vpbe_getZkappa2(thee->pbe);
    ionstr = Vpbe_getBulkIonicStrength(thee->pbe);

    /* Bail if we're at zero ionic strength */
    if (zkappa2 < VSMALL) {

#ifndef VAPBSQUIET
        Vnm_print(0, "Vpmg_qmEnergy:  Zero energy for zero ionic strength!\n");
#endif

        return 0.0;
    }
    zks2 = 0.5*zkappa2/ionstr;

    if (!thee->filled) {
        Vnm_print(2, "Vpmg_qmEnergy:  Need to call Vpmg_fillco()!\n");
        VASSERT(0);
    }

    energy = 0.0;
    nchop = 0;
    Vpbe_getIons(thee->pbe, &nion, ionConc, ionRadii, ionQ);
    if (thee->pmgp->nonlin) {
        Vnm_print(0, "Vpmg_qmEnergy:  Calculating nonlinear energy\n");
        for (i=0, len=nx*ny*nz; i<len; i++) {
            if (thee->pvec[i]*thee->kappa[i] > VSMALL) {
                for (j=0; j<nion; j++) {
                    energy += (thee->pvec[i]*thee->kappa[i]*zks2
                      * ionConc[j]
                      * (Vcap_exp(-ionQ[j]*thee->u[i], &ichop)-1.0));
                    nchop += ichop;
                }
            }
        }
        if (nchop > 0){
            Vnm_print(2, "Vpmg_qmEnergy:  Chopped EXP %d times!\n",nchop);
            Vnm_print(2, "\nERROR!  Detected large potential values in energy evaluation! \nERROR!  This calculation failed -- please report to the APBS developers!\n\n");
            VASSERT(0);
        }
    } else {
        /* Zkappa2 OK here b/c LPBE approx */
        Vnm_print(0, "Vpmg_qmEnergy:  Calculating linear energy\n");
        for (i=0, len=nx*ny*nz; i<len; i++) {
            if (thee->pvec[i]*thee->kappa[i] > VSMALL)
              energy += (thee->pvec[i]*zkappa2*thee->kappa[i]*VSQR(thee->u[i]));
        }
        energy = 0.5*energy;
    }
    energy = energy*hx*hy*hzed;
    energy = energy/Vpbe_getZmagic(thee->pbe);

    if (extFlag == 1) energy += thee->extQmEnergy;

    return energy;
}

VPUBLIC double Vpmg_qmEnergySMPBE(Vpmg *thee,
                                  int extFlag
                                 ) {

    double hx,
           hy,
           hzed,
           energy,
           ionConc[MAXION],
           ionRadii[MAXION],
           ionQ[MAXION],
           zkappa2,
           ionstr,
           zks2;
    int i,
        //j, // gcc: not used
        nx,
        ny,
        nz,
        nion,
        //ichop, // gcc: not used
        nchop,
        len; /* Loop variable */

    /* SMPB Modification (vchu, 09/21/06)*/
    /* variable declarations for SMPB energy terms */
    double a,
           k,
           z1,
           z2,
           z3,
           cb1,
           cb2,
           cb3,
           a1,
           a2,
           a3,
           c1,
           c2,
           c3,
           currEnergy,
           fracOccA,
           fracOccB,
           fracOccC,
           phi,
           gpark,
           denom;
           // Na; /**< @todo remove if no conflicts are caused - This constant is already defined in vpde.h.  no need to redefine. */
    int ichop1,
        ichop2,
        ichop3;

    VASSERT(thee != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    zkappa2 = Vpbe_getZkappa2(thee->pbe);
    ionstr = Vpbe_getBulkIonicStrength(thee->pbe);

    /* Bail if we're at zero ionic strength */
    if (zkappa2 < VSMALL) {

#ifndef VAPBSQUIET
        Vnm_print(0, "Vpmg_qmEnergySMPBE:  Zero energy for zero ionic strength!\n");
#endif

        return 0.0;
    }
    zks2 = 0.5*zkappa2/ionstr;

    if (!thee->filled) {
        Vnm_print(2, "Vpmg_qmEnergySMPBE:  Need to call Vpmg_fillco()!\n");
        VASSERT(0);
    }

    energy = 0.0;
    nchop = 0;
    Vpbe_getIons(thee->pbe, &nion, ionConc, ionRadii, ionQ);

    /* SMPB Modification (vchu, 09/21/06) */
    /* Extensive modification to the first part of the if statement
        where that handles the thee->pmgp->nonlin part. Basically, I've
        deleted all of the original code and written my own code that computes
        the electrostatic free energy in the SMPB framework. Definitely really hacky
        at this stage of the game, but gets the job done. The second part of the
        if statement (the part that handles linear poisson-boltzmann) has been deleted
        because there will be no linearized SMPB energy.. */

    z1 = ionQ[0];
    z2 = ionQ[1];
    z3 = ionQ[2];
    cb1 = ionConc[0];
    cb2 = ionConc[1];
    cb3 = ionConc[2];
    a  = thee->pbe->smvolume;
    k  = thee->pbe->smsize;

    /// @todo remove if no conflicts are caused
    // This constant is defined in vpde.h  Do not need to redefine
    //Na = 6.022045000e-04; /* Converts from Molar to N/A^3 */

    fracOccA = Na*cb1*VCUB(a);
    fracOccB = Na*cb2*VCUB(a);
    fracOccC = Na*cb3*VCUB(a);

    phi = (fracOccA/k) + fracOccB + fracOccC;

    if (thee->pmgp->nonlin) {
        Vnm_print(0, "Vpmg_qmEnergySMPBE:  Calculating nonlinear energy using SMPB functional!\n");
        for (i=0, len=nx*ny*nz; i<len; i++) {
            if (((k-1) > VSMALL) && (thee->pvec[i]*thee->kappa[i] > VSMALL)) {

                a1 = Vcap_exp(-1.0*z1*thee->u[i], &ichop1);
                a2 = Vcap_exp(-1.0*z2*thee->u[i], &ichop2);
                a3 = Vcap_exp(-1.0*z3*thee->u[i], &ichop3);

                nchop += ichop1 + ichop2 + ichop3;

                gpark = (1 - phi + (fracOccA/k)*a1);
                denom = VPOW(gpark, k) + VPOW(1-fracOccB-fracOccC, k-1)*(fracOccB*a2+fracOccC*a3);

                if (cb1 > VSMALL) {
                    c1 = Na*cb1*VPOW(gpark, k-1)*a1/denom;
                    if(c1 != c1) c1 = 0.;
                } else c1 = 0.;

                if (cb2 > VSMALL) {
                    c2 = Na*cb2*VPOW(1-fracOccB-fracOccC,k-1)*a2/denom;
                    if(c2 != c2) c2 = 0.;
                } else c2 = 0.;

                if (cb3 > VSMALL) {
                    c3 = Na*cb3*VPOW(1-fracOccB-fracOccC,k-1)*a3/denom;
                    if(c3 != c3) c3 = 0.;
                } else c3 = 0.;

                currEnergy = k*VLOG((1-(c1*VCUB(a)/k)-c2*VCUB(a)-c3*VCUB(a))/(1-phi))
                    -(k-1)*VLOG((1-c2*VCUB(a)-c3*VCUB(a))/(1-phi+(fracOccA/k)));

                energy += thee->pvec[i]*thee->kappa[i]*currEnergy;

            } else if (thee->pvec[i]*thee->kappa[i] > VSMALL){

                a1 = Vcap_exp(-1.0*z1*thee->u[i], &ichop1);
                a2 = Vcap_exp(-1.0*z2*thee->u[i], &ichop2);
                a3 = Vcap_exp(-1.0*z3*thee->u[i], &ichop3);

                nchop += ichop1 + ichop2 + ichop3;

                gpark = (1 - phi + (fracOccA)*a1);
                denom = gpark + (fracOccB*a2+fracOccC*a3);

                if (cb1 > VSMALL) {
                    c1 = Na*cb1*a1/denom;
                    if(c1 != c1) c1 = 0.;
                } else c1 = 0.;

                if (cb2 > VSMALL) {
                    c2 = Na*cb2*a2/denom;
                    if(c2 != c2) c2 = 0.;
                } else c2 = 0.;

                if (cb3 > VSMALL) {
                    c3 = Na*cb3*a3/denom;
                    if(c3 != c3) c3 = 0.;
                } else c3 = 0.;

                currEnergy = VLOG((1-c1*VCUB(a)-c2*VCUB(a)-c3*VCUB(a))/(1-fracOccA-fracOccB-fracOccC));

                energy += thee->pvec[i]*thee->kappa[i]*currEnergy;
            }
        }

        energy = -energy/VCUB(a);

        if (nchop > 0) Vnm_print(2, "Vpmg_qmEnergySMPBE:  Chopped EXP %d times!\n",
                                 nchop);

    } else {
        /* Zkappa2 OK here b/c LPBE approx */
        Vnm_print(0, "Vpmg_qmEnergySMPBE:  ERROR: NO LINEAR ENERGY!! Returning 0!\n");

        energy = 0.0;

    }
    energy = energy*hx*hy*hzed;

    if (extFlag == 1) energy += thee->extQmEnergy;

    return energy;
}

VPUBLIC double Vpmg_qfEnergy(Vpmg *thee,
                             int extFlag
                            ) {

    double energy = 0.0;

    VASSERT(thee != VNULL);

    if ((thee->useChargeMap) || (thee->chargeMeth == VCM_BSPL2)) {
        energy = Vpmg_qfEnergyVolume(thee, extFlag);
    } else {
        energy = Vpmg_qfEnergyPoint(thee, extFlag);
    }

    return energy;
}

VPRIVATE double Vpmg_qfEnergyPoint(Vpmg *thee,
                                   int extFlag
                                  ) {

    int iatom, nx, ny, nz, ihi, ilo, jhi, jlo, khi, klo;
    double xmax, ymax, zmax, xmin, ymin, zmin, hx, hy, hzed, ifloat, jfloat;
    double charge, kfloat, dx, dy, dz, energy, uval, *position;
    double *u;
    double *pvec;
    Valist *alist;
    Vatom *atom;
    Vpbe *pbe;

    pbe = thee->pbe;
    alist = pbe->alist;
    VASSERT(alist != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;

    u = thee->u;
    pvec = thee->pvec;

    energy = 0.0;

    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        /* Get atomic information */
        atom = Valist_getAtom(alist, iatom);

        position = Vatom_getPosition(atom);
        charge = Vatom_getCharge(atom);

        /* Figure out which vertices we're next to */
        ifloat = (position[0] - xmin)/hx;
        jfloat = (position[1] - ymin)/hy;
        kfloat = (position[2] - zmin)/hzed;
        ihi = (int)ceil(ifloat);
        ilo = (int)floor(ifloat);
        jhi = (int)ceil(jfloat);
        jlo = (int)floor(jfloat);
        khi = (int)ceil(kfloat);
        klo = (int)floor(kfloat);

        if (atom->partID > 0) {

            if ((ihi<nx) && (jhi<ny) && (khi<nz) &&
                (ilo>=0) && (jlo>=0) && (klo>=0)) {

                /* Now get trilinear interpolation constants */
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                uval =
                  dx*dy*dz*u[IJK(ihi,jhi,khi)]
                + dx*(1.0-dy)*dz*u[IJK(ihi,jlo,khi)]
                + dx*dy*(1.0-dz)*u[IJK(ihi,jhi,klo)]
                + dx*(1.0-dy)*(1.0-dz)*u[IJK(ihi,jlo,klo)]
                + (1.0-dx)*dy*dz*u[IJK(ilo,jhi,khi)]
                + (1.0-dx)*(1.0-dy)*dz*u[IJK(ilo,jlo,khi)]
                + (1.0-dx)*dy*(1.0-dz)*u[IJK(ilo,jhi,klo)]
                + (1.0-dx)*(1.0-dy)*(1.0-dz)*u[IJK(ilo,jlo,klo)];
                energy += (uval*charge*atom->partID);
            } else if (thee->pmgp->bcfl != BCFL_FOCUS) {
                Vnm_print(2, "Vpmg_qfEnergy:  Atom #%d at (%4.3f, %4.3f, \
%4.3f) is off the mesh (ignoring)!\n",
                iatom, position[0], position[1], position[2]);
            }
        }
    }

    if (extFlag) energy += thee->extQfEnergy;

    return energy;
}

VPUBLIC double Vpmg_qfAtomEnergy(Vpmg *thee, Vatom *atom) {

    int nx, ny, nz, ihi, ilo, jhi, jlo, khi, klo;
    double xmax, xmin, ymax, ymin, zmax, zmin, hx, hy, hzed, ifloat, jfloat;
    double charge, kfloat, dx, dy, dz, energy, uval, *position;
    double *u;


    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmax = thee->xf[nx-1];
    ymax = thee->yf[ny-1];
    zmax = thee->zf[nz-1];
    xmin = thee->xf[0];
    ymin = thee->yf[0];
    zmin = thee->zf[0];

    u = thee->u;

    energy = 0.0;


    position = Vatom_getPosition(atom);
    charge = Vatom_getCharge(atom);

    /* Figure out which vertices we're next to */
    ifloat = (position[0] - xmin)/hx;
    jfloat = (position[1] - ymin)/hy;
    kfloat = (position[2] - zmin)/hzed;
    ihi = (int)ceil(ifloat);
    ilo = (int)floor(ifloat);
    jhi = (int)ceil(jfloat);
    jlo = (int)floor(jfloat);
    khi = (int)ceil(kfloat);
    klo = (int)floor(kfloat);

    if (atom->partID > 0) {

        if ((ihi<nx) && (jhi<ny) && (khi<nz) &&
            (ilo>=0) && (jlo>=0) && (klo>=0)) {

            /* Now get trilinear interpolation constants */
            dx = ifloat - (double)(ilo);
            dy = jfloat - (double)(jlo);
            dz = kfloat - (double)(klo);
            uval =
              dx*dy*dz*u[IJK(ihi,jhi,khi)]
            + dx*(1.0-dy)*dz*u[IJK(ihi,jlo,khi)]
            + dx*dy*(1.0-dz)*u[IJK(ihi,jhi,klo)]
            + dx*(1.0-dy)*(1.0-dz)*u[IJK(ihi,jlo,klo)]
            + (1.0-dx)*dy*dz*u[IJK(ilo,jhi,khi)]
            + (1.0-dx)*(1.0-dy)*dz*u[IJK(ilo,jlo,khi)]
            + (1.0-dx)*dy*(1.0-dz)*u[IJK(ilo,jhi,klo)]
            + (1.0-dx)*(1.0-dy)*(1.0-dz)*u[IJK(ilo,jlo,klo)];
            energy += (uval*charge*atom->partID);
        } else if (thee->pmgp->bcfl != BCFL_FOCUS) {
            Vnm_print(2, "Vpmg_qfAtomEnergy:  Atom at (%4.3f, %4.3f, \
%4.3f) is off the mesh (ignoring)!\n",
            position[0], position[1], position[2]);
        }
    }

    return energy;
}

VPRIVATE double Vpmg_qfEnergyVolume(Vpmg *thee, int extFlag) {

    double hx, hy, hzed, energy;
    int i, nx, ny, nz;

    VASSERT(thee != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    if (!thee->filled) {
        Vnm_print(2, "Vpmg_qfEnergyVolume:  need to call Vpmg_fillco!\n");
        VASSERT(0);
    }

    energy = 0.0;
    Vnm_print(0, "Vpmg_qfEnergyVolume:  Calculating energy\n");
    for (i=0; i<(nx*ny*nz); i++) {
        energy += (thee->pvec[i]*thee->u[i]*thee->charge[i]);
    }
    energy = energy*hx*hy*hzed/Vpbe_getZmagic(thee->pbe);

    if (extFlag == 1) energy += thee->extQfEnergy;

    return energy;
}

VPRIVATE void Vpmg_splineSelect(int srfm,Vacc *acc,double *gpos,double win,
                                      double infrad,Vatom *atom,double *force){

    switch (srfm) {
        case VSM_SPLINE :
            Vacc_splineAccGradAtomNorm(acc, gpos, win, infrad, atom, force);
            break;
        case VSM_SPLINE3:
            Vacc_splineAccGradAtomNorm3(acc, gpos, win, infrad, atom, force);
            break;
        case VSM_SPLINE4 :
            Vacc_splineAccGradAtomNorm4(acc, gpos, win, infrad, atom, force);
            break;
        default:
            Vnm_print(2, "Vpmg_dbnbForce: Unknown surface method.\n");
            return;
    }

    return;
}

VPRIVATE void focusFillBound(Vpmg *thee,
                             Vpmg *pmgOLD
                            ) {

    Vpbe *pbe;
    double hxOLD,
           hyOLD,
           hzOLD,
           xminOLD,
           yminOLD,
           zminOLD,
           xmaxOLD,
           ymaxOLD,
           zmaxOLD,
           hxNEW,
           hyNEW,
           hzNEW,
           xminNEW,
           yminNEW,
           zminNEW,
           xmaxNEW,
           ymaxNEW,
           zmaxNEW,
           x,
           y,
           z,
           dx,
           dy,
           dz,
           ifloat,
           jfloat,
           kfloat,
           uval,
           eps_w,
           T,
           pre1,
           xkappa,
           size,
           *apos,
           charge,
           //pos[3], // gcc: not used
           uvalMin,
           uvalMax,
           *data;
    int nxOLD,
        nyOLD,
        nzOLD,
        nxNEW,
        nyNEW,
        nzNEW,
        i,
        j,
        k,
        ihi,
        ilo,
        jhi,
        jlo,
        khi,
        klo,
        nx,
        ny,
        nz;

    /* Calculate new problem dimensions */
    hxNEW = thee->pmgp->hx;
    hyNEW = thee->pmgp->hy;
    hzNEW = thee->pmgp->hzed;
    nx =  thee->pmgp->nx;
    ny =  thee->pmgp->ny;
    nz =  thee->pmgp->nz;
    nxNEW = thee->pmgp->nx;
    nyNEW = thee->pmgp->ny;
    nzNEW = thee->pmgp->nz;
    xminNEW = thee->pmgp->xcent - ((double)(nxNEW-1)*hxNEW)/2.0;
    xmaxNEW = thee->pmgp->xcent + ((double)(nxNEW-1)*hxNEW)/2.0;
    yminNEW = thee->pmgp->ycent - ((double)(nyNEW-1)*hyNEW)/2.0;
    ymaxNEW = thee->pmgp->ycent + ((double)(nyNEW-1)*hyNEW)/2.0;
    zminNEW = thee->pmgp->zcent - ((double)(nzNEW-1)*hzNEW)/2.0;
    zmaxNEW = thee->pmgp->zcent + ((double)(nzNEW-1)*hzNEW)/2.0;

    if(pmgOLD != VNULL){
        /* Relevant old problem parameters */
        hxOLD = pmgOLD->pmgp->hx;
        hyOLD = pmgOLD->pmgp->hy;
        hzOLD = pmgOLD->pmgp->hzed;
        nxOLD = pmgOLD->pmgp->nx;
        nyOLD = pmgOLD->pmgp->ny;
        nzOLD = pmgOLD->pmgp->nz;
        xminOLD = pmgOLD->pmgp->xcent - ((double)(nxOLD-1)*hxOLD)/2.0;
        xmaxOLD = pmgOLD->pmgp->xcent + ((double)(nxOLD-1)*hxOLD)/2.0;
        yminOLD = pmgOLD->pmgp->ycent - ((double)(nyOLD-1)*hyOLD)/2.0;
        ymaxOLD = pmgOLD->pmgp->ycent + ((double)(nyOLD-1)*hyOLD)/2.0;
        zminOLD = pmgOLD->pmgp->zcent - ((double)(nzOLD-1)*hzOLD)/2.0;
        zmaxOLD = pmgOLD->pmgp->zcent + ((double)(nzOLD-1)*hzOLD)/2.0;

        data = pmgOLD->u;
    }else{
        /* Relevant old problem parameters */
        hxOLD = thee->potMap->hx;
        hyOLD = thee->potMap->hy;
        hzOLD = thee->potMap->hzed;
        nxOLD = thee->potMap->nx;
        nyOLD = thee->potMap->ny;
        nzOLD = thee->potMap->nz;
        xminOLD = thee->potMap->xmin;
        xmaxOLD = thee->potMap->xmax;
        yminOLD = thee->potMap->ymin;
        ymaxOLD = thee->potMap->ymax;
        zminOLD = thee->potMap->zmin;
        zmaxOLD = thee->potMap->zmax;

        data = thee->potMap->data;
    }
    /* BOUNDARY CONDITION SETUP FOR POINTS OFF OLD MESH:
     * For each "atom" (only one for bcfl=1), we use the following formula to
     * calculate the boundary conditions:
     *    g(x) = \frac{q e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     *          * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     *          * 1/d
     * where d = ||x - x_0|| (in m) and a is the size of the atom (in m).
     * We only need to evaluate some of these prefactors once:
     *    pre1 = \frac{e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     * which gives the potential as
     *    g(x) = pre1 * q/d * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     */
    pbe = thee->pbe;
    eps_w = Vpbe_getSolventDiel(pbe);           /* Dimensionless */
    T = Vpbe_getTemperature(pbe);               /* K             */
    pre1 = (Vunit_ec)/(4*VPI*Vunit_eps0*eps_w*Vunit_kb*T);

    /* Finally, if we convert keep xkappa in A^{-1} and scale pre1 by
     * m/A, then we will only need to deal with distances and sizes in
     * Angstroms rather than meters.                                       */
    xkappa = Vpbe_getXkappa(pbe);              /* A^{-1}        */
    pre1 = pre1*(1.0e10);
    size = Vpbe_getSoluteRadius(pbe);
    apos = Vpbe_getSoluteCenter(pbe);
    charge = Vunit_ec*Vpbe_getSoluteCharge(pbe);

    /* Check for rounding error */
    if (VABS(xminOLD-xminNEW) < VSMALL) xminNEW = xminOLD;
    if (VABS(xmaxOLD-xmaxNEW) < VSMALL) xmaxNEW = xmaxOLD;
    if (VABS(yminOLD-yminNEW) < VSMALL) yminNEW = yminOLD;
    if (VABS(ymaxOLD-ymaxNEW) < VSMALL) ymaxNEW = ymaxOLD;
    if (VABS(zminOLD-zminNEW) < VSMALL) zminNEW = zminOLD;
    if (VABS(zmaxOLD-zmaxNEW) < VSMALL) zmaxNEW = zmaxOLD;


    /* Sanity check: make sure we're within the old mesh */
    Vnm_print(0, "VPMG::focusFillBound -- New mesh mins = %g, %g, %g\n",
              xminNEW, yminNEW, zminNEW);
    Vnm_print(0, "VPMG::focusFillBound -- New mesh maxs = %g, %g, %g\n",
              xmaxNEW, ymaxNEW, zmaxNEW);
    Vnm_print(0, "VPMG::focusFillBound -- Old mesh mins = %g, %g, %g\n",
              xminOLD, yminOLD, zminOLD);
    Vnm_print(0, "VPMG::focusFillBound -- Old mesh maxs = %g, %g, %g\n",
              xmaxOLD, ymaxOLD, zmaxOLD);

    /* The following is obsolete; we'll substitute analytical boundary
     * condition values when the new mesh falls outside the old */
    if ((xmaxNEW>xmaxOLD) || (ymaxNEW>ymaxOLD) || (zmaxNEW>zmaxOLD) ||
        (xminOLD>xminNEW) || (yminOLD>yminNEW) || (zminOLD>zminNEW)) {

        Vnm_print(2, "Vpmg::focusFillBound -- new mesh not contained in old!\n");
        Vnm_print(2, "Vpmg::focusFillBound -- old mesh min = (%g, %g, %g)\n",
                  xminOLD, yminOLD, zminOLD);
        Vnm_print(2, "Vpmg::focusFillBound -- old mesh max = (%g, %g, %g)\n",
                  xmaxOLD, ymaxOLD, zmaxOLD);
        Vnm_print(2, "Vpmg::focusFillBound -- new mesh min = (%g, %g, %g)\n",
                  xminNEW, yminNEW, zminNEW);
        Vnm_print(2, "Vpmg::focusFillBound -- new mesh max = (%g, %g, %g)\n",
                  xmaxNEW, ymaxNEW, zmaxNEW);
        fflush(stderr);
        VASSERT(0);
    }

    uvalMin	= VPMGSMALL;
    uvalMax = -VPMGSMALL;

    /* Fill the "i" boundaries (dirichlet) */
    for (k=0; k<nzNEW; k++) {
        for (j=0; j<nyNEW; j++) {
            /* Low X face */
            x = xminNEW;
            y = yminNEW + j*hyNEW;
            z = zminNEW + k*hzNEW;
            if ((x >= (xminOLD-VSMALL)) && (y >= (yminOLD-VSMALL)) && (z >= (zminOLD-VSMALL)) &&
                (x <= (xmaxOLD+VSMALL)) && (y <= (ymaxOLD+VSMALL)) && (z <= (zmaxOLD+VSMALL))) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(data[IJK(ihi,jhi,khi)])
                + dx*(1.0-dy)*dz*(data[IJK(ihi,jlo,khi)])
                + dx*dy*(1.0-dz)*(data[IJK(ihi,jhi,klo)])
                + dx*(1.0-dy)*(1.0-dz)*(data[IJK(ihi,jlo,klo)])
                + (1.0-dx)*dy*dz*(data[IJK(ilo,jhi,khi)])
                + (1.0-dx)*(1.0-dy)*dz*(data[IJK(ilo,jlo,khi)])
                + (1.0-dx)*dy*(1.0-dz)*(data[IJK(ilo,jhi,klo)])
                + (1.0-dx)*(1.0-dy)*(1.0-dz)*(data[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                Vnm_print(2, "focusFillBound (%s, %d):  Off old mesh at %g, %g \
                          %g!\n", __FILE__, __LINE__, x, y, z);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh lower corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xminOLD, yminOLD, zminOLD);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh upper corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xmaxOLD, ymaxOLD, zmaxOLD);
                VASSERT(0);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gxcf[IJKx(j,k,0)] = uval;
            if(uval < uvalMin) uvalMin = uval;
            if(uval > uvalMax) uvalMax = uval;

            /* High X face */
            x = xmaxNEW;
            if ((x >= (xminOLD-VSMALL)) && (y >= (yminOLD-VSMALL)) && (z >= (zminOLD-VSMALL)) &&
                (x <= (xmaxOLD+VSMALL)) && (y <= (ymaxOLD+VSMALL)) && (z <= (zmaxOLD+VSMALL))) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(data[IJK(ihi,jhi,khi)])
                + dx*(1.0-dy)*dz*(data[IJK(ihi,jlo,khi)])
                + dx*dy*(1.0-dz)*(data[IJK(ihi,jhi,klo)])
                + dx*(1.0-dy)*(1.0-dz)*(data[IJK(ihi,jlo,klo)])
                + (1.0-dx)*dy*dz*(data[IJK(ilo,jhi,khi)])
                + (1.0-dx)*(1.0-dy)*dz*(data[IJK(ilo,jlo,khi)])
                + (1.0-dx)*dy*(1.0-dz)*(data[IJK(ilo,jhi,klo)])
                + (1.0-dx)*(1.0-dy)*(1.0-dz)*(data[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                Vnm_print(2, "focusFillBound (%s, %d):  Off old mesh at %g, %g \
                          %g!\n", __FILE__, __LINE__, x, y, z);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh lower corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xminOLD, yminOLD, zminOLD);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh upper corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xmaxOLD, ymaxOLD, zmaxOLD);
                VASSERT(0);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gxcf[IJKx(j,k,1)] = uval;
            if(uval < uvalMin) uvalMin = uval;
            if(uval > uvalMax) uvalMax = uval;

            /* Zero Neumann conditions */
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gxcf[IJKx(j,k,2)] = 0.0;
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gxcf[IJKx(j,k,3)] = 0.0;
        }
    }

    /* Fill the "j" boundaries (dirichlet) */
    for (k=0; k<nzNEW; k++) {
        for (i=0; i<nxNEW; i++) {
            /* Low Y face */
            x = xminNEW + i*hxNEW;
            y = yminNEW;
            z = zminNEW + k*hzNEW;
            if ((x >= (xminOLD-VSMALL)) && (y >= (yminOLD-VSMALL)) && (z >= (zminOLD-VSMALL)) &&
                (x <= (xmaxOLD+VSMALL)) && (y <= (ymaxOLD+VSMALL)) && (z <= (zmaxOLD+VSMALL))) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(data[IJK(ihi,jhi,khi)])
                + dx*(1.0-dy)*dz*(data[IJK(ihi,jlo,khi)])
                + dx*dy*(1.0-dz)*(data[IJK(ihi,jhi,klo)])
                + dx*(1.0-dy)*(1.0-dz)*(data[IJK(ihi,jlo,klo)])
                + (1.0-dx)*dy*dz*(data[IJK(ilo,jhi,khi)])
                + (1.0-dx)*(1.0-dy)*dz*(data[IJK(ilo,jlo,khi)])
                + (1.0-dx)*dy*(1.0-dz)*(data[IJK(ilo,jhi,klo)])
                + (1.0-dx)*(1.0-dy)*(1.0-dz)*(data[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                Vnm_print(2, "focusFillBound (%s, %d):  Off old mesh at %g, %g \
                          %g!\n", __FILE__, __LINE__, x, y, z);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh lower corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xminOLD, yminOLD, zminOLD);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh upper corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xmaxOLD, ymaxOLD, zmaxOLD);
                VASSERT(0);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gycf[IJKy(i,k,0)] = uval;
            if(uval < uvalMin) uvalMin = uval;
            if(uval > uvalMax) uvalMax = uval;

            /* High Y face */
            y = ymaxNEW;
            if ((x >= (xminOLD-VSMALL)) && (y >= (yminOLD-VSMALL)) && (z >= (zminOLD-VSMALL)) &&
                (x <= (xmaxOLD+VSMALL)) && (y <= (ymaxOLD+VSMALL)) && (z <= (zmaxOLD+VSMALL))) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(data[IJK(ihi,jhi,khi)])
                + dx*(1.0-dy)*dz*(data[IJK(ihi,jlo,khi)])
                + dx*dy*(1.0-dz)*(data[IJK(ihi,jhi,klo)])
                + dx*(1.0-dy)*(1.0-dz)*(data[IJK(ihi,jlo,klo)])
                + (1.0-dx)*dy*dz*(data[IJK(ilo,jhi,khi)])
                + (1.0-dx)*(1.0-dy)*dz*(data[IJK(ilo,jlo,khi)])
                + (1.0-dx)*dy*(1.0-dz)*(data[IJK(ilo,jhi,klo)])
                + (1.0-dx)*(1.0-dy)*(1.0-dz)*(data[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                Vnm_print(2, "focusFillBound (%s, %d):  Off old mesh at %g, %g \
                          %g!\n", __FILE__, __LINE__, x, y, z);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh lower corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xminOLD, yminOLD, zminOLD);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh upper corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xmaxOLD, ymaxOLD, zmaxOLD);
                VASSERT(0);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gycf[IJKy(i,k,1)] = uval;
            if(uval < uvalMin) uvalMin = uval;
            if(uval > uvalMax) uvalMax = uval;

            /* Zero Neumann conditions */
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gycf[IJKy(i,k,2)] = 0.0;
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gycf[IJKy(i,k,3)] = 0.0;
        }
    }

    /* Fill the "k" boundaries (dirichlet) */
    for (j=0; j<nyNEW; j++) {
        for (i=0; i<nxNEW; i++) {
            /* Low Z face */
            x = xminNEW + i*hxNEW;
            y = yminNEW + j*hyNEW;
            z = zminNEW;
            if ((x >= (xminOLD-VSMALL)) && (y >= (yminOLD-VSMALL)) && (z >= (zminOLD-VSMALL)) &&
                (x <= (xmaxOLD+VSMALL)) && (y <= (ymaxOLD+VSMALL)) && (z <= (zmaxOLD+VSMALL))) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(data[IJK(ihi,jhi,khi)])
                + dx*(1.0-dy)*dz*(data[IJK(ihi,jlo,khi)])
                + dx*dy*(1.0-dz)*(data[IJK(ihi,jhi,klo)])
                + dx*(1.0-dy)*(1.0-dz)*(data[IJK(ihi,jlo,klo)])
                + (1.0-dx)*dy*dz*(data[IJK(ilo,jhi,khi)])
                + (1.0-dx)*(1.0-dy)*dz*(data[IJK(ilo,jlo,khi)])
                + (1.0-dx)*dy*(1.0-dz)*(data[IJK(ilo,jhi,klo)])
                + (1.0-dx)*(1.0-dy)*(1.0-dz)*(data[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                Vnm_print(2, "focusFillBound (%s, %d):  Off old mesh at %g, %g \
                          %g!\n", __FILE__, __LINE__, x, y, z);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh lower corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xminOLD, yminOLD, zminOLD);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh upper corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xmaxOLD, ymaxOLD, zmaxOLD);
                VASSERT(0);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gzcf[IJKz(i,j,0)] = uval;
            if(uval < uvalMin) uvalMin = uval;
            if(uval > uvalMax) uvalMax = uval;

            /* High Z face */
            z = zmaxNEW;
            if ((x >= (xminOLD-VSMALL)) && (y >= (yminOLD-VSMALL)) && (z >= (zminOLD-VSMALL)) &&
                (x <= (xmaxOLD+VSMALL)) && (y <= (ymaxOLD+VSMALL)) && (z <= (zmaxOLD+VSMALL))) {
                ifloat = (x - xminOLD)/hxOLD;
                jfloat = (y - yminOLD)/hyOLD;
                kfloat = (z - zminOLD)/hzOLD;
                ihi = (int)ceil(ifloat);
                if (ihi > (nxOLD-1)) ihi = nxOLD-1;
                ilo = (int)floor(ifloat);
                if (ilo < 0) ilo = 0;
                jhi = (int)ceil(jfloat);
                if (jhi > (nyOLD-1)) jhi = nyOLD-1;
                jlo = (int)floor(jfloat);
                if (jlo < 0) jlo = 0;
                khi = (int)ceil(kfloat);
                if (khi > (nzOLD-1)) khi = nzOLD-1;
                klo = (int)floor(kfloat);
                if (klo < 0) klo = 0;
                dx = ifloat - (double)(ilo);
                dy = jfloat - (double)(jlo);
                dz = kfloat - (double)(klo);
                nx = nxOLD; ny = nyOLD; nz = nzOLD;
                uval =  dx*dy*dz*(data[IJK(ihi,jhi,khi)])
                + dx*(1.0-dy)*dz*(data[IJK(ihi,jlo,khi)])
                + dx*dy*(1.0-dz)*(data[IJK(ihi,jhi,klo)])
                + dx*(1.0-dy)*(1.0-dz)*(data[IJK(ihi,jlo,klo)])
                + (1.0-dx)*dy*dz*(data[IJK(ilo,jhi,khi)])
                + (1.0-dx)*(1.0-dy)*dz*(data[IJK(ilo,jlo,khi)])
                + (1.0-dx)*dy*(1.0-dz)*(data[IJK(ilo,jhi,klo)])
                + (1.0-dx)*(1.0-dy)*(1.0-dz)*(data[IJK(ilo,jlo,klo)]);
                nx = nxNEW; ny = nyNEW; nz = nzNEW;
            } else {
                Vnm_print(2, "focusFillBound (%s, %d):  Off old mesh at %g, %g \
                          %g!\n", __FILE__, __LINE__, x, y, z);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh lower corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xminOLD, yminOLD, zminOLD);
                Vnm_print(2, "focusFillBound (%s, %d):  old mesh upper corner at \
                          %g %g %g.\n", __FILE__, __LINE__, xmaxOLD, ymaxOLD, zmaxOLD);
                VASSERT(0);
            }
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gzcf[IJKz(i,j,1)] = uval;
            if(uval < uvalMin) uvalMin = uval;
            if(uval > uvalMax) uvalMax = uval;

            /* Zero Neumann conditions */
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gzcf[IJKz(i,j,2)] = 0.0;
            nx = nxNEW; ny = nyNEW; nz = nzNEW;
            thee->gzcf[IJKz(i,j,3)] = 0.0;
        }
    }

    VWARN_MSG0(
        uvalMin >= SINH_MIN && uvalMax <= SINH_MAX,
        "Unusually large potential values\n"
        "    detected on the focusing boundary!\n"
        "    Convergence not guaranteed for NPBE/NRPBE calculations!"
        );
}

VPRIVATE void extEnergy(Vpmg *thee, Vpmg *pmgOLD, PBEparm_calcEnergy extFlag,
                        double partMin[3], double partMax[3], int bflags[6]) {

    Vatom *atom;
    double hxNEW, hyNEW, hzNEW;
    double lowerCorner[3], upperCorner[3];
    int nxNEW, nyNEW, nzNEW;
    int nxOLD, nyOLD, nzOLD;
    int i,j,k;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double hxOLD, hyOLD, hzOLD;
    double xval, yval, zval;
    double x,y,z;
    int nx, ny, nz;

    /* Set the new external energy contribution to zero.  Any external
     * contributions from higher levels will be included in the appropriate
     * energy function call. */
    thee->extQmEnergy = 0;
    thee->extQfEnergy = 0;
    thee->extDiEnergy = 0;

    /* New problem dimensions */
    hxNEW = thee->pmgp->hx;
    hyNEW = thee->pmgp->hy;
    hzNEW = thee->pmgp->hzed;
    nxNEW = thee->pmgp->nx;
    nyNEW = thee->pmgp->ny;
    nzNEW = thee->pmgp->nz;
    lowerCorner[0] = thee->pmgp->xcent - ((double)(nxNEW-1)*hxNEW)/2.0;
    upperCorner[0] = thee->pmgp->xcent + ((double)(nxNEW-1)*hxNEW)/2.0;
    lowerCorner[1] = thee->pmgp->ycent - ((double)(nyNEW-1)*hyNEW)/2.0;
    upperCorner[1] = thee->pmgp->ycent + ((double)(nyNEW-1)*hyNEW)/2.0;
    lowerCorner[2] = thee->pmgp->zcent - ((double)(nzNEW-1)*hzNEW)/2.0;
    upperCorner[2] = thee->pmgp->zcent + ((double)(nzNEW-1)*hzNEW)/2.0;

    Vnm_print(0, "VPMG::extEnergy:  energy flag = %d\n", extFlag);

    /* Old problem dimensions */
    nxOLD = pmgOLD->pmgp->nx;
    nyOLD = pmgOLD->pmgp->ny;
    nzOLD = pmgOLD->pmgp->nz;

    /* Create a partition based on the new problem dimensions */
    /* Vnm_print(1, "DEBUG (%s, %d):  extEnergy calling Vpmg_setPart for old PMG.\n",
     __FILE__, __LINE__); */
    Vpmg_setPart(pmgOLD, lowerCorner, upperCorner, bflags);


    Vnm_print(0,"VPMG::extEnergy:   Finding extEnergy dimensions...\n");
    Vnm_print(0,"VPMG::extEnergy    Disj part lower corner = (%g, %g, %g)\n",
              partMin[0], partMin[1], partMin[2]);
    Vnm_print(0,"VPMG::extEnergy    Disj part upper corner = (%g, %g, %g)\n",
              partMax[0], partMax[1], partMax[2]);

    /* Find the old dimensions */

    hxOLD = pmgOLD->pmgp->hx;
    hyOLD = pmgOLD->pmgp->hy;
    hzOLD = pmgOLD->pmgp->hzed;
    xmin =  pmgOLD->pmgp->xcent - 0.5*hxOLD*(nxOLD-1);
    ymin =  pmgOLD->pmgp->ycent - 0.5*hyOLD*(nyOLD-1);
    zmin =  pmgOLD->pmgp->zcent - 0.5*hzOLD*(nzOLD-1);
    xmax =  xmin+hxOLD*(nxOLD-1);
    ymax =  ymin+hyOLD*(nyOLD-1);
    zmax =  zmin+hzOLD*(nzOLD-1);

    Vnm_print(0,"VPMG::extEnergy    Old lower corner = (%g, %g, %g)\n",
              xmin, ymin, zmin);
    Vnm_print(0,"VPMG::extEnergy    Old upper corner = (%g, %g, %g)\n",
              xmax, ymax, zmax);

    /* Flip the partition, but do not include any points that will
     be included by another processor */

    nx = nxOLD;
    ny = nyOLD;
    nz = nzOLD;

    for(i=0; i<nx; i++) {
        xval = 1;
        x = i*hxOLD + xmin;
        if (x < partMin[0] && bflags[VAPBS_LEFT] == 1) xval = 0;
        else if (x > partMax[0] && bflags[VAPBS_RIGHT] == 1) xval = 0;

        for(j=0; j<ny; j++) {
            yval = 1;
            y = j*hyOLD + ymin;
            if (y < partMin[1] && bflags[VAPBS_BACK] == 1) yval = 0;
            else if (y > partMax[1] && bflags[VAPBS_FRONT] == 1) yval = 0;

            for(k=0; k<nz; k++) {
                zval = 1;
                z = k*hzOLD + zmin;
                if (z < partMin[2] && bflags[VAPBS_DOWN] == 1) zval = 0;
                else if (z > partMax[2] && bflags[VAPBS_UP] == 1) zval = 0;

                if (pmgOLD->pvec[IJK(i,j,k)] > VSMALL) pmgOLD->pvec[IJK(i,j,k)] = 1.0;
                pmgOLD->pvec[IJK(i,j,k)] = (1 - (pmgOLD->pvec[IJK(i,j,k)])) * (xval*yval*zval);
            }
        }
    }

    for (i=0; i<Valist_getNumberAtoms(thee->pbe->alist); i++) {
        xval=1;
        yval=1;
        zval=1;
        atom = Valist_getAtom(thee->pbe->alist, i);
        x = atom->position[0];
        y = atom->position[1];
        z = atom->position[2];
        if (x < partMin[0] && bflags[VAPBS_LEFT] == 1) xval = 0;
        else if (x > partMax[0] && bflags[VAPBS_RIGHT] == 1) xval = 0;
        if (y < partMin[1] && bflags[VAPBS_BACK] == 1) yval = 0;
        else if (y > partMax[1] && bflags[VAPBS_FRONT] == 1) yval = 0;
        if (z < partMin[2] && bflags[VAPBS_DOWN] == 1) zval = 0;
        else if (z > partMax[2] && bflags[VAPBS_UP] == 1) zval = 0;
        if (atom->partID > VSMALL) atom->partID = 1.0;
        atom->partID = (1 - atom->partID) * (xval*yval*zval);
    }

    /* Now calculate the energy on inverted subset of the domain */
    thee->extQmEnergy = Vpmg_qmEnergy(pmgOLD, 1);
    Vnm_print(0, "VPMG::extEnergy: extQmEnergy = %g kT\n", thee->extQmEnergy);
    thee->extQfEnergy = Vpmg_qfEnergy(pmgOLD, 1);
    Vnm_print(0, "VPMG::extEnergy: extQfEnergy = %g kT\n", thee->extQfEnergy);
    thee->extDiEnergy = Vpmg_dielEnergy(pmgOLD, 1);
    Vnm_print(0, "VPMG::extEnergy: extDiEnergy = %g kT\n", thee->extDiEnergy);
    Vpmg_unsetPart(pmgOLD);
}

VPRIVATE double bcfl1sp(double size, double *apos, double charge,
                        double xkappa, double pre1, double *pos) {

    double dist, val;

    dist = VSQRT(VSQR(pos[0]-apos[0]) + VSQR(pos[1]-apos[1])
                 + VSQR(pos[2]-apos[2]));
    if (xkappa > VSMALL) {
        val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
        / (1+xkappa*size);
    } else {
        val = pre1*(charge/dist);
    }

    return val;
}

VPRIVATE void bcfl1(double size, double *apos, double charge,
                    double xkappa, double pre1, double *gxcf, double *gycf, double *gzcf,
                    double *xf, double *yf, double *zf, int nx, int ny, int nz) {

    int i, j, k;
    double dist, val;
    double gpos[3];

    /* the "i" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        gpos[2] = zf[k];
        for (j=0; j<ny; j++) {
            gpos[1] = yf[j];
            gpos[0] = xf[0];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
                         + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gxcf[IJKx(j,k,0)] += val;
            gpos[0] = xf[nx-1];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
                         + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gxcf[IJKx(j,k,1)] += val;
        }
    }

    /* the "j" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        gpos[2] = zf[k];
        for (i=0; i<nx; i++) {
            gpos[0] = xf[i];
            gpos[1] = yf[0];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
                         + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gycf[IJKy(i,k,0)] += val;
            gpos[1] = yf[ny-1];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
                         + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gycf[IJKy(i,k,1)] += val;
        }
    }

    /* the "k" boundaries (dirichlet) */
    for (j=0; j<ny; j++) {
        gpos[1] = yf[j];
        for (i=0; i<nx; i++) {
            gpos[0] = xf[i];
            gpos[2] = zf[0];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
                         + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gzcf[IJKz(i,j,0)] += val;
            gpos[2] = zf[nz-1];
            dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
                         + VSQR(gpos[2]-apos[2]));
            if (xkappa > VSMALL) {
                val = pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                / (1+xkappa*size);
            } else {
                val = pre1*(charge/dist);
            }
            gzcf[IJKz(i,j,1)] += val;
        }
    }
}

VPRIVATE void bcfl2(double size, double *apos,
                    double charge, double *dipole, double *quad,
                    double xkappa, double eps_p, double eps_w, double T,
                    double *gxcf, double *gycf, double *gzcf,
                    double *xf, double *yf, double *zf,
                    int nx, int ny, int nz) {

    int i, j, k;
    double val;
    double gpos[3],tensor[3];
    double ux,uy,uz,xr,yr,zr;
    double qxx,qxy,qxz,qyx,qyy,qyz,qzx,qzy,qzz;
    double dist, pre;

    VASSERT(dipole != VNULL);
    ux = dipole[0];
    uy = dipole[1];
    uz = dipole[2];
    if (quad != VNULL) {
        /* The factor of 1/3 results from using a
         traceless quadrupole definition. See, for example,
         "The Theory of Intermolecular Forces" by A.J. Stone,
         Chapter 3. */
        qxx = quad[0] / 3.0;
        qxy = quad[1] / 3.0;
        qxz = quad[2] / 3.0;
        qyx = quad[3] / 3.0;
        qyy = quad[4] / 3.0;
        qyz = quad[5] / 3.0;
        qzx = quad[6] / 3.0;
        qzy = quad[7] / 3.0;
        qzz = quad[8] / 3.0;
    } else {
        qxx = 0.0;
        qxy = 0.0;
        qxz = 0.0;
        qyx = 0.0;
        qyy = 0.0;
        qyz = 0.0;
        qzx = 0.0;
        qzy = 0.0;
        qzz = 0.0;
    }

    pre = (Vunit_ec*Vunit_ec)/(4*VPI*Vunit_eps0*Vunit_kb*T);
    pre = pre*(1.0e10);

    /* the "i" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        gpos[2] = zf[k];
        for (j=0; j<ny; j++) {
            gpos[1] = yf[j];
            gpos[0] = xf[0];
            xr = gpos[0] - apos[0];
            yr = gpos[1] - apos[1];
            zr = gpos[2] - apos[2];
            dist = VSQRT(VSQR(xr) + VSQR(yr) + VSQR(zr));
            multipolebc(dist, xkappa, eps_p, eps_w, size, tensor);
            val = pre*charge*tensor[0];
            val -= pre*ux*xr*tensor[1];
            val -= pre*uy*yr*tensor[1];
            val -= pre*uz*zr*tensor[1];
            val += pre*qxx*xr*xr*tensor[2];
            val += pre*qyy*yr*yr*tensor[2];
            val += pre*qzz*zr*zr*tensor[2];
            val += pre*2.0*qxy*xr*yr*tensor[2];
            val += pre*2.0*qxz*xr*zr*tensor[2];
            val += pre*2.0*qyz*yr*zr*tensor[2];
            gxcf[IJKx(j,k,0)] += val;

            gpos[0] = xf[nx-1];
            xr = gpos[0] - apos[0];
            dist = VSQRT(VSQR(xr) + VSQR(yr) + VSQR(zr));
            multipolebc(dist, xkappa, eps_p, eps_w, size, tensor);
            val = pre*charge*tensor[0];
            val -= pre*ux*xr*tensor[1];
            val -= pre*uy*yr*tensor[1];
            val -= pre*uz*zr*tensor[1];
            val += pre*qxx*xr*xr*tensor[2];
            val += pre*qyy*yr*yr*tensor[2];
            val += pre*qzz*zr*zr*tensor[2];
            val += pre*2.0*qxy*xr*yr*tensor[2];
            val += pre*2.0*qxz*xr*zr*tensor[2];
            val += pre*2.0*qyz*yr*zr*tensor[2];
            gxcf[IJKx(j,k,1)] += val;
        }
    }

    /* the "j" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        gpos[2] = zf[k];
        for (i=0; i<nx; i++) {
            gpos[0] = xf[i];
            gpos[1] = yf[0];
            xr = gpos[0] - apos[0];
            yr = gpos[1] - apos[1];
            zr = gpos[2] - apos[2];
            dist = VSQRT(VSQR(xr) + VSQR(yr) + VSQR(zr));
            multipolebc(dist, xkappa, eps_p, eps_w, size, tensor);
            val = pre*charge*tensor[0];
            val -= pre*ux*xr*tensor[1];
            val -= pre*uy*yr*tensor[1];
            val -= pre*uz*zr*tensor[1];
            val += pre*qxx*xr*xr*tensor[2];
            val += pre*qyy*yr*yr*tensor[2];
            val += pre*qzz*zr*zr*tensor[2];
            val += pre*2.0*qxy*xr*yr*tensor[2];
            val += pre*2.0*qxz*xr*zr*tensor[2];
            val += pre*2.0*qyz*yr*zr*tensor[2];
            gycf[IJKy(i,k,0)] += val;

            gpos[1] = yf[ny-1];
            yr = gpos[1] - apos[1];
            dist = VSQRT(VSQR(xr) + VSQR(yr) + VSQR(zr));
            multipolebc(dist, xkappa, eps_p, eps_w, size, tensor);
            val = pre*charge*tensor[0];
            val -= pre*ux*xr*tensor[1];
            val -= pre*uy*yr*tensor[1];
            val -= pre*uz*zr*tensor[1];
            val += pre*qxx*xr*xr*tensor[2];
            val += pre*qyy*yr*yr*tensor[2];
            val += pre*qzz*zr*zr*tensor[2];
            val += pre*2.0*qxy*xr*yr*tensor[2];
            val += pre*2.0*qxz*xr*zr*tensor[2];
            val += pre*2.0*qyz*yr*zr*tensor[2];
            gycf[IJKy(i,k,1)] += val;
        }
    }

    /* the "k" boundaries (dirichlet) */
    for (j=0; j<ny; j++) {
        gpos[1] = yf[j];
        for (i=0; i<nx; i++) {
            gpos[0] = xf[i];
            gpos[2] = zf[0];
            xr = gpos[0] - apos[0];
            yr = gpos[1] - apos[1];
            zr = gpos[2] - apos[2];
            dist = VSQRT(VSQR(xr) + VSQR(yr) + VSQR(zr));
            multipolebc(dist, xkappa, eps_p, eps_w, size, tensor);
            val = pre*charge*tensor[0];
            val -= pre*ux*xr*tensor[1];
            val -= pre*uy*yr*tensor[1];
            val -= pre*uz*zr*tensor[1];
            val += pre*qxx*xr*xr*tensor[2];
            val += pre*qyy*yr*yr*tensor[2];
            val += pre*qzz*zr*zr*tensor[2];
            val += pre*2.0*qxy*xr*yr*tensor[2];
            val += pre*2.0*qxz*xr*zr*tensor[2];
            val += pre*2.0*qyz*yr*zr*tensor[2];
            gzcf[IJKz(i,j,0)] += val;

            gpos[2] = zf[nz-1];
            zr = gpos[2] - apos[2];
            dist = VSQRT(VSQR(xr) + VSQR(yr) + VSQR(zr));
            multipolebc(dist, xkappa, eps_p, eps_w, size, tensor);
            val = pre*charge*tensor[0];
            val -= pre*ux*xr*tensor[1];
            val -= pre*uy*yr*tensor[1];
            val -= pre*uz*zr*tensor[1];
            val += pre*qxx*xr*xr*tensor[2];
            val += pre*qyy*yr*yr*tensor[2];
            val += pre*qzz*zr*zr*tensor[2];
            val += pre*2.0*qxy*xr*yr*tensor[2];
            val += pre*2.0*qxz*xr*zr*tensor[2];
            val += pre*2.0*qyz*yr*zr*tensor[2];
            gzcf[IJKz(i,j,1)] += val;
        }
    }
}

VPRIVATE void bcCalcOrig(Vpmg *thee) {

    int nx, ny, nz;
    double size, *position, charge, xkappa, eps_w, T, pre1;
    double *dipole, *quadrupole, debye, eps_p;
    double xr,yr,zr,qave,*apos;
    double sdhcharge, sdhdipole[3], traced[9], sdhquadrupole[9];
    int i, j, k, iatom;
    Vpbe *pbe;
    Vatom *atom;
    Valist *alist;

    pbe = thee->pbe;
    alist = thee->pbe->alist;
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;

    /* Zero out the boundaries */
    /* the "i" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            thee->gxcf[IJKx(j,k,0)] = 0.0;
            thee->gxcf[IJKx(j,k,1)] = 0.0;
            thee->gxcf[IJKx(j,k,2)] = 0.0;
            thee->gxcf[IJKx(j,k,3)] = 0.0;
        }
    }

    /* the "j" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        for (i=0; i<nx; i++) {
            thee->gycf[IJKy(i,k,0)] = 0.0;
            thee->gycf[IJKy(i,k,1)] = 0.0;
            thee->gycf[IJKy(i,k,2)] = 0.0;
            thee->gycf[IJKy(i,k,3)] = 0.0;
        }
    }

    /* the "k" boundaries (dirichlet) */
    for (j=0; j<ny; j++) {
        for (i=0; i<nx; i++) {
            thee->gzcf[IJKz(i,j,0)] = 0.0;
            thee->gzcf[IJKz(i,j,1)] = 0.0;
            thee->gzcf[IJKz(i,j,2)] = 0.0;
            thee->gzcf[IJKz(i,j,3)] = 0.0;
        }
    }

    /* For each "atom" (only one for bcfl=1), we use the following formula to
    * calculate the boundary conditions:
    *    g(x) = \frac{q e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
    *          * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
    *          * 1/d
    * where d = ||x - x_0|| (in m) and a is the size of the atom (in m).
    * We only need to evaluate some of these prefactors once:
    *    pre1 = \frac{e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
    * which gives the potential as
    *    g(x) = pre1 * q/d * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
    */
    eps_w = Vpbe_getSolventDiel(pbe);           /* Dimensionless */
    eps_p = Vpbe_getSoluteDiel(pbe);           /* Dimensionless */
    T = Vpbe_getTemperature(pbe);               /* K             */
    pre1 = (Vunit_ec)/(4*VPI*Vunit_eps0*eps_w*Vunit_kb*T);

    /* Finally, if we convert keep xkappa in A^{-1} and scale pre1 by
        * m/A, then we will only need to deal with distances and sizes in
        * Angstroms rather than meters.                                       */
    xkappa = Vpbe_getXkappa(pbe);              /* A^{-1}        */
    pre1 = pre1*(1.0e10);

    switch (thee->pmgp->bcfl) {
        /*  If we have zero boundary conditions, we're done */
        case BCFL_ZERO:
            return;

            /*  For single DH sphere BC's, we only have one "atom" to deal with;
            *  get its information and */
        case BCFL_SDH:
            size = Vpbe_getSoluteRadius(pbe);
            position = Vpbe_getSoluteCenter(pbe);

            /*
             For AMOEBA SDH boundary conditions, we need to find the
             total monopole, dipole and traceless quadrupole moments
             of either the permanent multipoles, induced dipoles or
             non-local induced dipoles.
             */

            sdhcharge = 0.0;
            for (i=0; i<3; i++) sdhdipole[i] = 0.0;
            for (i=0; i<9; i++) sdhquadrupole[i] = 0.0;

            for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
                atom = Valist_getAtom(alist, iatom);
                apos = Vatom_getPosition(atom);
                xr = apos[0] - position[0];
                yr = apos[1] - position[1];
                zr = apos[2] - position[2];
                switch (thee->chargeSrc) {
                    case VCM_CHARGE:
                        charge = Vatom_getCharge(atom);
                        sdhcharge += charge;
                        sdhdipole[0] += xr * charge;
                        sdhdipole[1] += yr * charge;
                        sdhdipole[2] += zr * charge;
                        traced[0] = xr*xr*charge;
                        traced[1] = xr*yr*charge;
                        traced[2] = xr*zr*charge;
                        traced[3] = yr*xr*charge;
                        traced[4] = yr*yr*charge;
                        traced[5] = yr*zr*charge;
                        traced[6] = zr*xr*charge;
                        traced[7] = zr*yr*charge;
                        traced[8] = zr*zr*charge;
                        qave = (traced[0] + traced[4] + traced[8]) / 3.0;
                        sdhquadrupole[0] += 1.5*(traced[0] - qave);
                        sdhquadrupole[1] += 1.5*(traced[1]);
                        sdhquadrupole[2] += 1.5*(traced[2]);
                        sdhquadrupole[3] += 1.5*(traced[3]);
                        sdhquadrupole[4] += 1.5*(traced[4] - qave);
                        sdhquadrupole[5] += 1.5*(traced[5]);
                        sdhquadrupole[6] += 1.5*(traced[6]);
                        sdhquadrupole[7] += 1.5*(traced[7]);
                        sdhquadrupole[8] += 1.5*(traced[8] - qave);
#if defined(WITH_TINKER)
                    case VCM_PERMANENT:
                        charge = Vatom_getCharge(atom);
                        dipole = Vatom_getDipole(atom);
                        quadrupole = Vatom_getQuadrupole(atom);
                        sdhcharge += charge;
                        sdhdipole[0] += xr * charge;
                        sdhdipole[1] += yr * charge;
                        sdhdipole[2] += zr * charge;
                        traced[0] = xr*xr*charge;
                        traced[1] = xr*yr*charge;
                        traced[2] = xr*zr*charge;
                        traced[3] = yr*xr*charge;
                        traced[4] = yr*yr*charge;
                        traced[5] = yr*zr*charge;
                        traced[6] = zr*xr*charge;
                        traced[7] = zr*yr*charge;
                        traced[8] = zr*zr*charge;
                        sdhdipole[0] += dipole[0];
                        sdhdipole[1] += dipole[1];
                        sdhdipole[2] += dipole[2];
                        traced[0] += 2.0*xr*dipole[0];
                        traced[1] += xr*dipole[1] + yr*dipole[0];
                        traced[2] += xr*dipole[2] + zr*dipole[0];
                        traced[3] += yr*dipole[0] + xr*dipole[1];
                        traced[4] += 2.0*yr*dipole[1];
                        traced[5] += yr*dipole[2] + zr*dipole[1];
                        traced[6] += zr*dipole[0] + xr*dipole[2];
                        traced[7] += zr*dipole[1] + yr*dipole[2];
                        traced[8] += 2.0*zr*dipole[2];
                        qave = (traced[0] + traced[4] + traced[8]) / 3.0;
                        sdhquadrupole[0] += 1.5*(traced[0] - qave);
                        sdhquadrupole[1] += 1.5*(traced[1]);
                        sdhquadrupole[2] += 1.5*(traced[2]);
                        sdhquadrupole[3] += 1.5*(traced[3]);
                        sdhquadrupole[4] += 1.5*(traced[4] - qave);
                        sdhquadrupole[5] += 1.5*(traced[5]);
                        sdhquadrupole[6] += 1.5*(traced[6]);
                        sdhquadrupole[7] += 1.5*(traced[7]);
                        sdhquadrupole[8] += 1.5*(traced[8] - qave);
                        sdhquadrupole[0] += quadrupole[0];
                        sdhquadrupole[1] += quadrupole[1];
                        sdhquadrupole[2] += quadrupole[2];
                        sdhquadrupole[3] += quadrupole[3];
                        sdhquadrupole[4] += quadrupole[4];
                        sdhquadrupole[5] += quadrupole[5];
                        sdhquadrupole[6] += quadrupole[6];
                        sdhquadrupole[7] += quadrupole[7];
                        sdhquadrupole[8] += quadrupole[8];
                    case VCM_INDUCED:
                        dipole = Vatom_getInducedDipole(atom);
                        sdhdipole[0] += dipole[0];
                        sdhdipole[1] += dipole[1];
                        sdhdipole[2] += dipole[2];
                        traced[0] = 2.0*xr*dipole[0];
                        traced[1] = xr*dipole[1] + yr*dipole[0];
                        traced[2] = xr*dipole[2] + zr*dipole[0];
                        traced[3] = yr*dipole[0] + xr*dipole[1];
                        traced[4] = 2.0*yr*dipole[1];
                        traced[5] = yr*dipole[2] + zr*dipole[1];
                        traced[6] = zr*dipole[0] + xr*dipole[2];
                        traced[7] = zr*dipole[1] + yr*dipole[2];
                        traced[8] = 2.0*zr*dipole[2];
                        qave = (traced[0] + traced[4] + traced[8]) / 3.0;
                        sdhquadrupole[0] += 1.5*(traced[0] - qave);
                        sdhquadrupole[1] += 1.5*(traced[1]);
                        sdhquadrupole[2] += 1.5*(traced[2]);
                        sdhquadrupole[3] += 1.5*(traced[3]);
                        sdhquadrupole[4] += 1.5*(traced[4] - qave);
                        sdhquadrupole[5] += 1.5*(traced[5]);
                        sdhquadrupole[6] += 1.5*(traced[6]);
                        sdhquadrupole[7] += 1.5*(traced[7]);
                        sdhquadrupole[8] += 1.5*(traced[8] - qave);
                    case VCM_NLINDUCED:
                        dipole = Vatom_getNLInducedDipole(atom);
                        sdhdipole[0] += dipole[0];
                        sdhdipole[1] += dipole[1];
                        sdhdipole[2] += dipole[2];
                        traced[0] = 2.0*xr*dipole[0];
                        traced[1] = xr*dipole[1] + yr*dipole[0];
                        traced[2] = xr*dipole[2] + zr*dipole[0];
                        traced[3] = yr*dipole[0] + xr*dipole[1];
                        traced[4] = 2.0*yr*dipole[1];
                        traced[5] = yr*dipole[2] + zr*dipole[1];
                        traced[6] = zr*dipole[0] + xr*dipole[2];
                        traced[7] = zr*dipole[1] + yr*dipole[2];
                        traced[8] = 2.0*zr*dipole[2];
                        qave = (traced[0] + traced[4] + traced[8]) / 3.0;
                        sdhquadrupole[0] += 1.5*(traced[0] - qave);
                        sdhquadrupole[1] += 1.5*(traced[1]);
                        sdhquadrupole[2] += 1.5*(traced[2]);
                        sdhquadrupole[3] += 1.5*(traced[3]);
                        sdhquadrupole[4] += 1.5*(traced[4] - qave);
                        sdhquadrupole[5] += 1.5*(traced[5]);
                        sdhquadrupole[6] += 1.5*(traced[6]);
                        sdhquadrupole[7] += 1.5*(traced[7]);
                        sdhquadrupole[8] += 1.5*(traced[8] - qave);
#endif /* if defined(WITH_TINKER) */
                }
            }
            /* SDH dipole and traceless quadrupole values
             were checked against similar routines in TINKER
             for large proteins.

             debye=4.8033324;
             printf("%6.3f, %6.3f, %6.3f\n", sdhdipole[0]*debye,
             sdhdipole[1]*debye, sdhdipole[2]*debye);
             printf("%6.3f\n", sdhquadrupole[0]*debye);
             printf("%6.3f %6.3f\n", sdhquadrupole[3]*debye,
             sdhquadrupole[4]*debye);
             printf("%6.3f %6.3f %6.3f\n", sdhquadrupole[6]*debye,
             sdhquadrupole[7]*debye, sdhquadrupole[8]*debye);
             */

            bcfl2(size, position, sdhcharge, sdhdipole, sdhquadrupole,
                  xkappa, eps_p, eps_w, T, thee->gxcf, thee->gycf,
                  thee->gzcf, thee->xf, thee->yf, thee->zf, nx, ny, nz);
            break;

        case BCFL_MDH:
            for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
                atom = Valist_getAtom(alist, iatom);
                position = Vatom_getPosition(atom);
                charge = Vunit_ec*Vatom_getCharge(atom);
                dipole = VNULL;
                quadrupole = VNULL;
                size = Vatom_getRadius(atom);
                switch (thee->chargeSrc)
                {
                    case VCM_CHARGE:
                        ;
#if  defined(WITH_TINKER)
                    case VCM_PERMANENT:
                        dipole = Vatom_getDipole(atom);
                        quadrupole = Vatom_getQuadrupole(atom);

                    case VCM_INDUCED:
                        dipole = Vatom_getInducedDipole(atom);

                    case VCM_NLINDUCED:
                        dipole = Vatom_getNLInducedDipole(atom);
#endif
                }
                bcfl1(size, position, charge, xkappa, pre1,
                      thee->gxcf, thee->gycf, thee->gzcf,
                      thee->xf, thee->yf, thee->zf, nx, ny, nz);
            }
            break;

        case BCFL_UNUSED:
            Vnm_print(2, "bcCalc:  Invalid bcfl (%d)!\n", thee->pmgp->bcfl);
            VASSERT(0);

        case BCFL_FOCUS:
            Vnm_print(2, "VPMG::bcCalc -- not appropriate for focusing!\n");
            VASSERT(0);

        default:
            Vnm_print(2, "VPMG::bcCalc -- invalid boundary condition \
flag (%d)!\n", thee->pmgp->bcfl);
            VASSERT(0);
    }
}

/*
 Used by bcflnew
 */
VPRIVATE int gridPointIsValid(int i, int j, int k, int nx, int ny, int nz){

    int isValid = 0;

    if((k==0) || (k==nz-1)){
        isValid = 1;
    }else if((j==0) || (j==ny-1)){
        isValid = 1;
    }else if((i==0) || (i==nx-1)){
        isValid = 1;
    }

    return isValid;
}

/*
 Used by bcflnew
 */
#ifdef DEBUG_MAC_OSX_OCL
#include "mach_chud.h"
VPRIVATE void packAtomsOpenCL(float *ax, float *ay, float *az,
                        float *charge, float *size, Vpmg *thee){

    int i;
    int natoms;

    Vatom *atom = VNULL;
    Valist *alist = VNULL;

    alist = thee->pbe->alist;
    natoms = Valist_getNumberAtoms(alist);

    for(i=0;i<natoms;i++){
        atom = &(alist->atoms[i]);
        charge[i] = Vunit_ec*atom->charge;
        ax[i] = atom->position[0];
        ay[i] = atom->position[1];
        az[i] = atom->position[2];
        size[i] = atom->radius;
    }
}

/*
 Used by bcflnew
 */
VPRIVATE void packUnpackOpenCL(int nx, int ny, int nz, int ngrid,
                         float *gx, float *gy, float *gz, float *value,
                         Vpmg *thee, int pack){

    int i,j,k,igrid;
    int x0,x1,y0,y1,z0,z1;

    float gpos[3];
    double *xf, *yf, *zf;
    double *gxcf, *gycf, *gzcf;

    xf = thee->xf;
    yf = thee->yf;
    zf = thee->zf;

    gxcf = thee->gxcf;
    gycf = thee->gycf;
    gzcf = thee->gzcf;

    igrid = 0;
    for(k=0;k<nz;k++){
        gpos[2] = zf[k];
        for(j=0;j<ny;j++){
            gpos[1] = yf[j];
            for(i=0;i<nx;i++){
                gpos[0] = xf[i];
                if(gridPointIsValid(i, j, k, nx, ny, nz)){
                    if(pack != 0){
                        gx[igrid] = gpos[0];
                        gy[igrid] = gpos[1];
                        gz[igrid] = gpos[2];

                        value[igrid] = 0.0;
                    }else{
                        x0 = IJKx(j,k,0);
                        x1 = IJKx(j,k,1);
                        y0 = IJKy(i,k,0);
                        y1 = IJKy(i,k,1);
                        z0 = IJKz(i,j,0);
                        z1 = IJKz(i,j,1);

                        if(i==0){
                            gxcf[x0] += value[igrid];
                        }
                        if(i==nx-1){
                            gxcf[x1] += value[igrid];
                        }
                        if(j==0){
                            gycf[y0] += value[igrid];
                        }
                        if(j==ny-1){
                            gycf[y1] += value[igrid];
                        }
                        if(k==0){
                            gzcf[z0] += value[igrid];
                        }
                        if(k==nz-1){
                            gzcf[z1] += value[igrid];
                        }
                    }

                    igrid++;
                } //end is valid point
            } //end i
        } //end j
    } //end k

}

/*
 bcflnew is an optimized replacement for bcfl1. bcfl1 is still used when TINKER
 support is compiled in.
 bcflnew uses: packUnpack, packAtoms, gridPointIsValid
 */
VPRIVATE void bcflnewOpenCL(Vpmg *thee){

    int i,j,k, iatom, igrid;
    int x0, x1, y0, y1, z0, z1;

    int nx, ny, nz;
    int natoms, ngrid, ngadj;

    float dist, pre1, eps_w, eps_p, T, xkappa;

    float *ax, *ay, *az;
    float *charge, *size, *val;

    float *gx, *gy, *gz;

    Vpbe *pbe = thee->pbe;

    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;

    eps_w = Vpbe_getSolventDiel(pbe);           /* Dimensionless */
    eps_p = Vpbe_getSoluteDiel(pbe);           /* Dimensionless */
    T = Vpbe_getTemperature(pbe);               /* K             */
    pre1 = ((Vunit_ec)/(4*VPI*Vunit_eps0*eps_w*Vunit_kb*T))*(1.0e10);
    xkappa = Vpbe_getXkappa(pbe);

    natoms = Valist_getNumberAtoms(thee->pbe->alist);
    ngrid = 2*((nx*ny) + (ny*nz) + (nx*nz));
    ngadj = ngrid + (512 - (ngrid & 511));

    ax = (float*)malloc(natoms * sizeof(float));
    ay = (float*)malloc(natoms * sizeof(float));
    az = (float*)malloc(natoms * sizeof(float));

    charge = (float*)malloc(natoms * sizeof(float));
    size = (float*)malloc(natoms * sizeof(float));

    gx = (float*)malloc(ngrid * sizeof(float));
    gy = (float*)malloc(ngrid * sizeof(float));
    gz = (float*)malloc(ngrid * sizeof(float));

    val = (float*)malloc(ngrid * sizeof(float));

    packAtomsOpenCL(ax,ay,az,charge,size,thee);
    packUnpackOpenCL(nx,ny,nz,ngrid,gx,gy,gz,val,thee,1);

    runMDHCL(ngrid,natoms,ngadj,ax,ay,az,gx,gy,gz,charge,size,xkappa,pre1,val);

    packUnpackOpenCL(nx,ny,nz,ngrid,gx,gy,gz,val,thee,0);

    free(ax);
    free(ay);
    free(az);
    free(charge);
    free(size);

    free(gx);
    free(gy);
    free(gz);
    free(val);
}
#endif

VPRIVATE void packAtoms(double *ax, double *ay, double *az,
                        double *charge, double *size, Vpmg *thee){

    int i;
    int natoms;

    Vatom *atom = VNULL;
    Valist *alist = VNULL;

    alist = thee->pbe->alist;
    natoms = Valist_getNumberAtoms(alist);

    for(i=0;i<natoms;i++){
        atom = &(alist->atoms[i]);
        charge[i] = Vunit_ec*atom->charge;
        ax[i] = atom->position[0];
        ay[i] = atom->position[1];
        az[i] = atom->position[2];
        size[i] = atom->radius;
    }
}

/*
 Used by bcflnew
 */
VPRIVATE void packUnpack(int nx, int ny, int nz, int ngrid,
                         double *gx, double *gy, double *gz, double *value,
                         Vpmg *thee, int pack){

    int i,j,k,igrid;
    int x0,x1,y0,y1,z0,z1;

    double gpos[3];
    double *xf, *yf, *zf;
    double *gxcf, *gycf, *gzcf;

    xf = thee->xf;
    yf = thee->yf;
    zf = thee->zf;

    gxcf = thee->gxcf;
    gycf = thee->gycf;
    gzcf = thee->gzcf;

    igrid = 0;
    for(k=0;k<nz;k++){
        gpos[2] = zf[k];
        for(j=0;j<ny;j++){
            gpos[1] = yf[j];
            for(i=0;i<nx;i++){
                gpos[0] = xf[i];
                if(gridPointIsValid(i, j, k, nx, ny, nz)){
                    if(pack != 0){
                        gx[igrid] = gpos[0];
                        gy[igrid] = gpos[1];
                        gz[igrid] = gpos[2];

                        value[igrid] = 0.0;
                    }else{
                        x0 = IJKx(j,k,0);
                        x1 = IJKx(j,k,1);
                        y0 = IJKy(i,k,0);
                        y1 = IJKy(i,k,1);
                        z0 = IJKz(i,j,0);
                        z1 = IJKz(i,j,1);

                        if(i==0){
                            gxcf[x0] += value[igrid];
                        }
                        if(i==nx-1){
                            gxcf[x1] += value[igrid];
                        }
                        if(j==0){
                            gycf[y0] += value[igrid];
                        }
                        if(j==ny-1){
                            gycf[y1] += value[igrid];
                        }
                        if(k==0){
                            gzcf[z0] += value[igrid];
                        }
                        if(k==nz-1){
                            gzcf[z1] += value[igrid];
                        }
                    }

                    igrid++;
                } //end is valid point
            } //end i
        } //end j
    } //end k

}

VPRIVATE void bcflnew(Vpmg *thee){

    int i,j,k, iatom, igrid;
    int x0, x1, y0, y1, z0, z1;

    int nx, ny, nz;
    int natoms, ngrid;

    double dist, pre1, eps_w, eps_p, T, xkappa;

    double *ax, *ay, *az;
    double *charge, *size, *val;

    double *gx, *gy, *gz;

    Vpbe *pbe = thee->pbe;

    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;

    eps_w = Vpbe_getSolventDiel(pbe);           /* Dimensionless */
    eps_p = Vpbe_getSoluteDiel(pbe);           /* Dimensionless */
    T = Vpbe_getTemperature(pbe);               /* K             */
    pre1 = ((Vunit_ec)/(4*VPI*Vunit_eps0*eps_w*Vunit_kb*T))*(1.0e10);
    xkappa = Vpbe_getXkappa(pbe);

    natoms = Valist_getNumberAtoms(thee->pbe->alist);
    ngrid = 2*((nx*ny) + (ny*nz) + (nx*nz));

    ax = (double*)malloc(natoms * sizeof(double));
    ay = (double*)malloc(natoms * sizeof(double));
    az = (double*)malloc(natoms * sizeof(double));

    charge = (double*)malloc(natoms * sizeof(double));
    size = (double*)malloc(natoms * sizeof(double));

    gx = (double*)malloc(ngrid * sizeof(double));
    gy = (double*)malloc(ngrid * sizeof(double));
    gz = (double*)malloc(ngrid * sizeof(double));

    val = (double*)malloc(ngrid * sizeof(double));

    packAtoms(ax,ay,az,charge,size,thee);
    packUnpack(nx,ny,nz,ngrid,gx,gy,gz,val,thee,1);

    if(xkappa > VSMALL){
#pragma omp parallel for default(shared) private(igrid,iatom,dist)
        for(igrid=0;igrid<ngrid;igrid++){
            for(iatom=0; iatom<natoms; iatom++){
                dist = VSQRT(VSQR(gx[igrid]-ax[iatom]) + VSQR(gy[igrid]-ay[iatom])
                             + VSQR(gz[igrid]-az[iatom]));
                val[igrid] += pre1*(charge[iatom]/dist)*VEXP(-xkappa*(dist-size[iatom]))
                / (1+xkappa*size[iatom]);
            }
        }
    }else{
#pragma omp parallel for default(shared) private(igrid,iatom,dist)
        for(igrid=0;igrid<ngrid;igrid++){
            for(iatom=0; iatom<natoms; iatom++){
                dist = VSQRT(VSQR(gx[igrid]-ax[iatom]) + VSQR(gy[igrid]-ay[iatom])
                             + VSQR(gz[igrid]-az[iatom]));
                val[igrid] += pre1*(charge[iatom]/dist);
            }
        }
    }
    packUnpack(nx,ny,nz,ngrid,gx,gy,gz,val,thee,0);

    free(ax);
    free(ay);
    free(az);
    free(charge);
    free(size);

    free(gx);
    free(gy);
    free(gz);
    free(val);
}

VPRIVATE void multipolebc(double r, double kappa, double eps_p,
                          double eps_w, double rad, double tsr[3]) {
    double r2,r3,r5;
    double eps_r;
    double ka,ka2,ka3;
    double kr,kr2,kr3;

    /*
     Below an attempt is made to explain the potential outside of a
     multipole located at the center of spherical cavity of dieletric
     eps_p, with dielectric eps_w outside (and possibly kappa > 0).


     First, eps_p = 1.0
     eps_w = 1.0
     kappa = 0.0

     The general form for the potential of a traceless multipole tensor
     of rank n in vacuum is:

     V(r) = (-1)^n * u . n . Del^n (1/r)

     where
     u                     is a multipole of order n (3^n components)
     u . n. Del^n (1/r)    is the contraction of u with the nth
     derivative of 1/r

     for example, if n = 1, the dipole potential is
     V_vac(r) = (-1)*[ux*x + uy*y + uz*z]/r^3

     This function returns the parts of V(r) for multipoles of
     order 0, 1 and 2 that are independent of the contraction.

     For the vacuum example, this would be 1/r, -1/r^3 and 3/r^5
     respectively.

     *** Note that this requires that the quadrupole is
     traceless. If not, the diagonal quadrupole potential changes
     from
     qaa *  3*a^2/r^5
     to
     qaa * (3*a^2/r^5 - 1/r^3a )
     where we sum over the trace; a = x, y and z.

     (In other words, the -1/r^3 term cancels for a traceless quadrupole.
     qxx + qyy + qzz = 0
     such that
     -(qxx + qyy + qzz)/r^3 = 0

     For quadrupole with trace:
     qxx + qyy + qzz != 0
     such that
     -(qxx + qyy + qzz)/r^3 != 0
     )

     ========================================================================

     eps_p != 1 or eps_w != 1
     kappa = 0.0

     If the multipole is placed at the center of a sphere with
     dieletric eps_p in a solvent of dielectric eps_w, the potential
     outside the sphere is the solution to the Laplace equation:

     V(r) = 1/eps_w * (2*n+1)*eps_r/(n+(n+1)*eps_r)
     * (-1)^n * u . n . Del^n (1/r)
     where
     eps_r = eps_w / eps_p
     is the ratio of solvent to solute dielectric

     ========================================================================

     kappa > 0

     Finally, if the region outside the sphere is treated by the linearized
     PB equation with Debye-Huckel parameter kappa, the solution is:

     V(r) = kappa/eps_w * alpha_n(kappa*a) * K_n(kappa*r) * r^(n+1)/a^n
     * (-1)^n * u . n . Del^n (1/r)
     where
     alpha_n(x) is [(2n + 1) / x] / [(n*K_n(x)/eps_r) - x*K_n'(x)]
     K_n(x) are modified spherical Bessel functions of the third kind.
     K_n'(x) is the derivative of K_n(x)
     */

    eps_r = eps_w/eps_p;
    r2 = r*r;
    r3 = r2*r;
    r5 = r3*r2;
    tsr[0] = (1.0/eps_w)/r;
    tsr[1] = (1.0/eps_w)*(-1.0)/r3;
    tsr[2] = (1.0/eps_w)*(3.0)/r5;
    if (kappa < VSMALL) {
        tsr[1] = (3.0*eps_r)/(1.0 + 2.0*eps_r)*tsr[1];
        tsr[2] = (5.0*eps_r)/(2.0 + 3.0*eps_r)*tsr[2];
    } else {
        ka = kappa*rad;
        ka2 = ka*ka;
        ka3 = ka2*ka;
        kr = kappa*r;
        kr2 = kr*kr;
        kr3 = kr2*kr;
        tsr[0] = exp(ka-kr) / (1.0 + ka) * tsr[0];
        tsr[1] = 3.0*eps_r*exp(ka-kr)*(1.0 + kr) * tsr[1];
        tsr[1] = tsr[1] / (1.0 + ka + eps_r*(2.0 + 2.0*ka + ka2));
        tsr[2] = 5.0*eps_r*exp(ka-kr)*(3.0 + 3.0*kr + kr2) * tsr[2];
        tsr[2] = tsr[2]/(6.0+6.0*ka+2.0*ka2+eps_r*(9.0+9.0*ka+4.0*ka2+ka3));
    }
}

VPRIVATE void bcfl_sdh(Vpmg *thee){

    int i,j,k,iatom;
    int nx, ny, nz;

    double size, *position, charge, xkappa, eps_w, eps_p, T, pre, dist;
    double sdhcharge, sdhdipole[3], traced[9], sdhquadrupole[9];
    double *dipole, *quadrupole;

    double val, *apos, gpos[3], tensor[3], qave;
    double ux, uy, uz, xr, yr, zr;
    double qxx,qxy,qxz,qyx,qyy,qyz,qzx,qzy,qzz;

    double *xf, *yf, *zf;
    double *gxcf, *gycf, *gzcf;

    Vpbe *pbe;
    Vatom *atom;
    Valist *alist;

    pbe = thee->pbe;
    alist = thee->pbe->alist;
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;

    xf = thee->xf;
    yf = thee->yf;
    zf = thee->zf;

    gxcf = thee->gxcf;
    gycf = thee->gycf;
    gzcf = thee->gzcf;

    /* For each "atom" (only one for bcfl=1), we use the following formula to
     * calculate the boundary conditions:
     *    g(x) = \frac{q e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     *          * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     *          * 1/d
     * where d = ||x - x_0|| (in m) and a is the size of the atom (in m).
     * We only need to evaluate some of these prefactors once:
     *    pre1 = \frac{e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     * which gives the potential as
     *    g(x) = pre1 * q/d * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     */
    eps_w = Vpbe_getSolventDiel(pbe);           /* Dimensionless */
    eps_p = Vpbe_getSoluteDiel(pbe);           /* Dimensionless */
    T = Vpbe_getTemperature(pbe);               /* K             */

    pre = (Vunit_ec*Vunit_ec)/(4*VPI*Vunit_eps0*Vunit_kb*T);
    pre = pre*(1.0e10);

    /* Finally, if we convert keep xkappa in A^{-1} and scale pre1 by
     * m/A, then we will only need to deal with distances and sizes in
     * Angstroms rather than meters.                                       */
    xkappa = Vpbe_getXkappa(pbe);              /* A^{-1}        */

    /* Solute size and position */
    size = Vpbe_getSoluteRadius(pbe);
    position = Vpbe_getSoluteCenter(pbe);

    sdhcharge = 0.0;
    for (i=0; i<3; i++) sdhdipole[i] = 0.0;
    for (i=0; i<9; i++) sdhquadrupole[i] = 0.0;

    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        xr = apos[0] - position[0];
        yr = apos[1] - position[1];
        zr = apos[2] - position[2];
        switch (thee->chargeSrc) {
            case VCM_CHARGE:
                charge = Vatom_getCharge(atom);
                sdhcharge += charge;
                sdhdipole[0] += xr * charge;
                sdhdipole[1] += yr * charge;
                sdhdipole[2] += zr * charge;
                traced[0] = xr*xr*charge;
                traced[1] = xr*yr*charge;
                traced[2] = xr*zr*charge;
                traced[3] = yr*xr*charge;
                traced[4] = yr*yr*charge;
                traced[5] = yr*zr*charge;
                traced[6] = zr*xr*charge;
                traced[7] = zr*yr*charge;
                traced[8] = zr*zr*charge;
                qave = (traced[0] + traced[4] + traced[8]) / 3.0;
                sdhquadrupole[0] += 1.5*(traced[0] - qave);
                sdhquadrupole[1] += 1.5*(traced[1]);
                sdhquadrupole[2] += 1.5*(traced[2]);
                sdhquadrupole[3] += 1.5*(traced[3]);
                sdhquadrupole[4] += 1.5*(traced[4] - qave);
                sdhquadrupole[5] += 1.5*(traced[5]);
                sdhquadrupole[6] += 1.5*(traced[6]);
                sdhquadrupole[7] += 1.5*(traced[7]);
                sdhquadrupole[8] += 1.5*(traced[8] - qave);
#if defined(WITH_TINKER)
            case VCM_PERMANENT:
                charge = Vatom_getCharge(atom);
                dipole = Vatom_getDipole(atom);
                quadrupole = Vatom_getQuadrupole(atom);
                sdhcharge += charge;
                sdhdipole[0] += xr * charge;
                sdhdipole[1] += yr * charge;
                sdhdipole[2] += zr * charge;
                traced[0] = xr*xr*charge;
                traced[1] = xr*yr*charge;
                traced[2] = xr*zr*charge;
                traced[3] = yr*xr*charge;
                traced[4] = yr*yr*charge;
                traced[5] = yr*zr*charge;
                traced[6] = zr*xr*charge;
                traced[7] = zr*yr*charge;
                traced[8] = zr*zr*charge;
                sdhdipole[0] += dipole[0];
                sdhdipole[1] += dipole[1];
                sdhdipole[2] += dipole[2];
                traced[0] += 2.0*xr*dipole[0];
                traced[1] += xr*dipole[1] + yr*dipole[0];
                traced[2] += xr*dipole[2] + zr*dipole[0];
                traced[3] += yr*dipole[0] + xr*dipole[1];
                traced[4] += 2.0*yr*dipole[1];
                traced[5] += yr*dipole[2] + zr*dipole[1];
                traced[6] += zr*dipole[0] + xr*dipole[2];
                traced[7] += zr*dipole[1] + yr*dipole[2];
                traced[8] += 2.0*zr*dipole[2];
                qave = (traced[0] + traced[4] + traced[8]) / 3.0;
                sdhquadrupole[0] += 1.5*(traced[0] - qave);
                sdhquadrupole[1] += 1.5*(traced[1]);
                sdhquadrupole[2] += 1.5*(traced[2]);
                sdhquadrupole[3] += 1.5*(traced[3]);
                sdhquadrupole[4] += 1.5*(traced[4] - qave);
                sdhquadrupole[5] += 1.5*(traced[5]);
                sdhquadrupole[6] += 1.5*(traced[6]);
                sdhquadrupole[7] += 1.5*(traced[7]);
                sdhquadrupole[8] += 1.5*(traced[8] - qave);
                sdhquadrupole[0] += quadrupole[0];
                sdhquadrupole[1] += quadrupole[1];
                sdhquadrupole[2] += quadrupole[2];
                sdhquadrupole[3] += quadrupole[3];
                sdhquadrupole[4] += quadrupole[4];
                sdhquadrupole[5] += quadrupole[5];
                sdhquadrupole[6] += quadrupole[6];
                sdhquadrupole[7] += quadrupole[7];
                sdhquadrupole[8] += quadrupole[8];
            case VCM_INDUCED:
                dipole = Vatom_getInducedDipole(atom);
                sdhdipole[0] += dipole[0];
                sdhdipole[1] += dipole[1];
                sdhdipole[2] += dipole[2];
                traced[0] = 2.0*xr*dipole[0];
                traced[1] = xr*dipole[1] + yr*dipole[0];
                traced[2] = xr*dipole[2] + zr*dipole[0];
                traced[3] = yr*dipole[0] + xr*dipole[1];
                traced[4] = 2.0*yr*dipole[1];
                traced[5] = yr*dipole[2] + zr*dipole[1];
                traced[6] = zr*dipole[0] + xr*dipole[2];
                traced[7] = zr*dipole[1] + yr*dipole[2];
                traced[8] = 2.0*zr*dipole[2];
                qave = (traced[0] + traced[4] + traced[8]) / 3.0;
                sdhquadrupole[0] += 1.5*(traced[0] - qave);
                sdhquadrupole[1] += 1.5*(traced[1]);
                sdhquadrupole[2] += 1.5*(traced[2]);
                sdhquadrupole[3] += 1.5*(traced[3]);
                sdhquadrupole[4] += 1.5*(traced[4] - qave);
                sdhquadrupole[5] += 1.5*(traced[5]);
                sdhquadrupole[6] += 1.5*(traced[6]);
                sdhquadrupole[7] += 1.5*(traced[7]);
                sdhquadrupole[8] += 1.5*(traced[8] - qave);
            case VCM_NLINDUCED:
                dipole = Vatom_getNLInducedDipole(atom);
                sdhdipole[0] += dipole[0];
                sdhdipole[1] += dipole[1];
                sdhdipole[2] += dipole[2];
                traced[0] = 2.0*xr*dipole[0];
                traced[1] = xr*dipole[1] + yr*dipole[0];
                traced[2] = xr*dipole[2] + zr*dipole[0];
                traced[3] = yr*dipole[0] + xr*dipole[1];
                traced[4] = 2.0*yr*dipole[1];
                traced[5] = yr*dipole[2] + zr*dipole[1];
                traced[6] = zr*dipole[0] + xr*dipole[2];
                traced[7] = zr*dipole[1] + yr*dipole[2];
                traced[8] = 2.0*zr*dipole[2];
                qave = (traced[0] + traced[4] + traced[8]) / 3.0;
                sdhquadrupole[0] += 1.5*(traced[0] - qave);
                sdhquadrupole[1] += 1.5*(traced[1]);
                sdhquadrupole[2] += 1.5*(traced[2]);
                sdhquadrupole[3] += 1.5*(traced[3]);
                sdhquadrupole[4] += 1.5*(traced[4] - qave);
                sdhquadrupole[5] += 1.5*(traced[5]);
                sdhquadrupole[6] += 1.5*(traced[6]);
                sdhquadrupole[7] += 1.5*(traced[7]);
                sdhquadrupole[8] += 1.5*(traced[8] - qave);
#endif /* if defined(WITH_TINKER) */
        }
    }

    ux = sdhdipole[0];
    uy = sdhdipole[1];
    uz = sdhdipole[2];

    /* The factor of 1/3 results from using a
     traceless quadrupole definition. See, for example,
     "The Theory of Intermolecular Forces" by A.J. Stone,
     Chapter 3. */
    qxx = sdhquadrupole[0] / 3.0;
    qxy = sdhquadrupole[1] / 3.0;
    qxz = sdhquadrupole[2] / 3.0;
    qyx = sdhquadrupole[3] / 3.0;
    qyy = sdhquadrupole[4] / 3.0;
    qyz = sdhquadrupole[5] / 3.0;
    qzx = sdhquadrupole[6] / 3.0;
    qzy = sdhquadrupole[7] / 3.0;
    qzz = sdhquadrupole[8] / 3.0;

    for(k=0;k<nz;k++){
        gpos[2] = zf[k];
        for(j=0;j<ny;j++){
            gpos[1] = yf[j];
            for(i=0;i<nx;i++){
                gpos[0] = xf[i];
                if(gridPointIsValid(i, j, k, nx, ny, nz)){
                    xr = gpos[0] - position[0];
                    yr = gpos[1] - position[1];
                    zr = gpos[2] - position[2];

                    dist = VSQRT(VSQR(xr) + VSQR(yr) + VSQR(zr));
                    multipolebc(dist, xkappa, eps_p, eps_w, size, tensor);

                    val = pre*sdhcharge*tensor[0];
                    val -= pre*ux*xr*tensor[1];
                    val -= pre*uy*yr*tensor[1];
                    val -= pre*uz*zr*tensor[1];
                    val += pre*qxx*xr*xr*tensor[2];
                    val += pre*qyy*yr*yr*tensor[2];
                    val += pre*qzz*zr*zr*tensor[2];
                    val += pre*2.0*qxy*xr*yr*tensor[2];
                    val += pre*2.0*qxz*xr*zr*tensor[2];
                    val += pre*2.0*qyz*yr*zr*tensor[2];

                    if(i==0){
                        gxcf[IJKx(j,k,0)] = val;
                    }
                    if(i==nx-1){
                        gxcf[IJKx(j,k,1)] = val;
                    }
                    if(j==0){
                        gycf[IJKy(i,k,0)] = val;
                    }
                    if(j==ny-1){
                        gycf[IJKy(i,k,1)] = val;
                    }
                    if(k==0){
                        gzcf[IJKz(i,j,0)] = val;
                    }
                    if(k==nz-1){
                        gzcf[IJKz(i,j,1)] = val;
                    }
                } /* End grid point is valid */
            } /* End i loop */
        } /* End j loop */
    } /* End k loop */

}

VPRIVATE void bcfl_mdh(Vpmg *thee){

    int i,j,k,iatom;
    int nx, ny, nz;

    double val, *apos, gpos[3];
    double *dipole, *quadrupole;
    double size, charge, xkappa, eps_w, eps_p, T, pre1, dist;

    double *xf, *yf, *zf;
    double *gxcf, *gycf, *gzcf;

    Vpbe *pbe;
    Vatom *atom;
    Valist *alist;

    pbe = thee->pbe;
    alist = thee->pbe->alist;
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;

    xf = thee->xf;
    yf = thee->yf;
    zf = thee->zf;

    gxcf = thee->gxcf;
    gycf = thee->gycf;
    gzcf = thee->gzcf;

    /* For each "atom" (only one for bcfl=1), we use the following formula to
     * calculate the boundary conditions:
     *    g(x) = \frac{q e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     *          * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     *          * 1/d
     * where d = ||x - x_0|| (in m) and a is the size of the atom (in m).
     * We only need to evaluate some of these prefactors once:
     *    pre1 = \frac{e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     * which gives the potential as
     *    g(x) = pre1 * q/d * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     */
    eps_w = Vpbe_getSolventDiel(pbe);           /* Dimensionless */
    eps_p = Vpbe_getSoluteDiel(pbe);           /* Dimensionless */
    T = Vpbe_getTemperature(pbe);               /* K             */
    pre1 = (Vunit_ec)/(4*VPI*Vunit_eps0*eps_w*Vunit_kb*T);

    /* Finally, if we convert keep xkappa in A^{-1} and scale pre1 by
     * m/A, then we will only need to deal with distances and sizes in
     * Angstroms rather than meters.                                       */
    xkappa = Vpbe_getXkappa(pbe);              /* A^{-1}        */
    pre1 = pre1*(1.0e10);

    /* Finally, if we convert keep xkappa in A^{-1} and scale pre1 by
     * m/A, then we will only need to deal with distances and sizes in
     * Angstroms rather than meters.                                       */
    xkappa = Vpbe_getXkappa(pbe);              /* A^{-1}        */

    for(k=0;k<nz;k++){
        gpos[2] = zf[k];
        for(j=0;j<ny;j++){
            gpos[1] = yf[j];
            for(i=0;i<nx;i++){
                gpos[0] = xf[i];
                if(gridPointIsValid(i, j, k, nx, ny, nz)){

                    val = 0.0;

                    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
                        atom = Valist_getAtom(alist, iatom);
                        apos = Vatom_getPosition(atom);
                        charge = Vunit_ec*Vatom_getCharge(atom);
                        size = Vatom_getRadius(atom);

                        dist = VSQRT(VSQR(gpos[0]-apos[0]) + VSQR(gpos[1]-apos[1])
                                     + VSQR(gpos[2]-apos[2]));
                        if (xkappa > VSMALL) {
                            val += pre1*(charge/dist)*VEXP(-xkappa*(dist-size))
                            / (1+xkappa*size);
                        } else {
                            val += pre1*(charge/dist);
                        }

                    }

                    if(i==0){
                        gxcf[IJKx(j,k,0)] = val;
                    }
                    if(i==nx-1){
                        gxcf[IJKx(j,k,1)] = val;
                    }
                    if(j==0){
                        gycf[IJKy(i,k,0)] = val;
                    }
                    if(j==ny-1){
                        gycf[IJKy(i,k,1)] = val;
                    }
                    if(k==0){
                        gzcf[IJKz(i,j,0)] = val;
                    }
                    if(k==nz-1){
                        gzcf[IJKz(i,j,1)] = val;
                    }
                } /* End grid point is valid */
            } /* End i loop */
        } /* End j loop */
    } /* End k loop */

}

/* ///////////////////////////////////////////////////////////////////////////
 // Routine:  bcfl_mem
 //
 // Purpose:  Increment all the boundary points by the
 //           analytic expression for a membrane system in
 //           the presence of a membrane potential. This
 //           Boundary flag should only be used for systems
 //           that explicitly have membranes in the dielectric
 //           and solvent maps.
 //
 //           There should be several input variables add to this
 //           function such as membrane potential, membrane thickness
 //           and height.
 //
 // Args:     apos is a 3-vector
 //
 // Author: Michael Grabe
 /////////////////////////////////////////////////////////////////////////// */
VPRIVATE void bcfl_mem(double zmem, double L, double eps_m, double eps_w,
                    double V, double xkappa, double *gxcf, double *gycf, double *gzcf,
                    double *xf, double *yf, double *zf, int nx, int ny, int nz) {

    ///////////////////////////////////////////////////
    /* some definitions                              */
    /* L = total length of the membrane              */
    /* xkappa = inverse Debeye length                */
    /* zmem = z value of membrane bottom (Cytoplasm) */
    /* V = electrical potential inside the cell      */
    ///////////////////////////////////////////////////
    int i, j, k;
    double dist, val, z_low, z_high, z_shift;
    double A, B, C, D, edge_L, l;
    double G, z_0, z_rel;
    double gpos[3];

    Vnm_print(0, "Here is the value of kappa: %f\n",xkappa);
    Vnm_print(0, "Here is the value of L: %f\n",L);
    Vnm_print(0, "Here is the value of zmem: %f\n",zmem);
    Vnm_print(0, "Here is the value of mdie: %f\n",eps_m);
    Vnm_print(0, "Here is the value of memv: %f\n",V);

    /* no salt symmetric BC's at +/- infinity */
    // B=V/(edge_L - l*(1-eps_w/eps_m));
    // A=V + B*edge_L;
    // D=eps_w/eps_m*B;
    z_low = zmem;     /* this defines the bottom of the membrane */
    z_high = zmem + L;  /* this is the top of the membrane */

    /******************************************************/
    /* proper boundary conditions for V = 0 extracellular */
    /* and psi=-V cytoplasm.                              */
    /* Implicit in this formulation is that the membrane  */
    /* center be at z = 0                                 */
    /******************************************************/

    l=L/2;                     /* half of the membrane length */
    z_0 = z_low + l;           /* center of the membrane      */
    G=l*eps_w/eps_m*xkappa;
    A=-V/2*(1/(G+1))*exp(xkappa*l);
    B=V/2;
    C=-V/2*eps_w/eps_m*xkappa*(1/(G+1));
    D=-A;
    /* The analytic expression for the boundary conditions      */
    /* had the cytoplasmic surface of the membrane set to zero. */
    /* This requires an off-set of the BC equations.            */

    /* the "i" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        gpos[2] = zf[k];
        z_rel = gpos[2] - z_0;    /* relative position for BCs */

        for (j=0; j<ny; j++) {

            if (gpos[2] <= z_low) {                       /* cytoplasmic */

                val = A*exp(xkappa*z_rel) + V;
                gxcf[IJKx(j,k,0)] += val;    /* assign low side BC */
                gxcf[IJKx(j,k,1)] += val;    /* assign high side BC */

            }

            else if (gpos[2] > z_low && gpos[2] <= z_high) {  /* in membrane */

                val = B + C*z_rel;
                gxcf[IJKx(j,k,0)] += val;    /* assign low side BC */
                gxcf[IJKx(j,k,1)] += val;    /* assign high side BC */

            }

            else if (gpos[2] > z_high)  {                  /* extracellular */

                val = D*exp(-xkappa*z_rel);
                gxcf[IJKx(j,k,0)] += val;    /* assign low side BC */
                gxcf[IJKx(j,k,1)] += val;    /* assign high side BC */

            }

        }
    }

    /* the "j" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        gpos[2] = zf[k];
        z_rel = gpos[2] - z_0;
        for (i=0; i<nx; i++) {

            if (gpos[2] <= z_low) {                       /* cytoplasmic */

                val = A*exp(xkappa*z_rel) + V;
                gycf[IJKy(i,k,0)] += val;    /* assign low side BC */
                gycf[IJKy(i,k,1)] += val;    /* assign high side BC */
                //printf("%f \n",val);

            }

            else if (gpos[2] > z_low && gpos[2] <= z_high) {  /* in membrane */

                val = B + C*z_rel;
                gycf[IJKy(i,k,0)] += val;    /* assign low side BC */
                gycf[IJKy(i,k,1)] += val;    /* assign high side BC */
                //printf("%f \n",val);

            }
            else if (gpos[2] > z_high)  {                  /* extracellular */

                val = D*exp(-xkappa*z_rel);
                gycf[IJKy(i,k,0)] += val;    /* assign low side BC */
                gycf[IJKy(i,k,1)] += val;    /* assign high side BC */
                //printf("%f \n",val);

            }

        }
    }

    /* the "k" boundaries (dirichlet) */
    for (j=0; j<ny; j++) {
        for (i=0; i<nx; i++) {

            /* first assign the bottom boundary */

            gpos[2] = zf[0];
            z_rel = gpos[2] - z_0;

            if (gpos[2] <= z_low) {                       /* cytoplasmic */

                val = A*exp(xkappa*z_rel) + V;
                gzcf[IJKz(i,j,0)] += val;    /* assign low side BC */
                //printf("%f \n",val);

            }

            else if (gpos[2] > z_low && gpos[2] <= z_high) {  /* in membrane */

                val = B + C*z_rel;
                gzcf[IJKz(i,j,0)] += val;    /* assign low side BC */

            }

            else if (gpos[2] > z_high)  {                  /* extracellular */

                val = D*exp(-xkappa*z_rel);
                gzcf[IJKz(i,j,0)] += val;    /* assign low side BC */

            }

            /* now assign the top boundary */

            gpos[2] = zf[nz-1];
            z_rel = gpos[2] - z_0;

            if (gpos[2] <= z_low) {                       /* cytoplasmic */

                val = A*exp(xkappa*z_rel) + V;
                gzcf[IJKz(i,j,1)] += val;    /* assign high side BC */

            }

            else if (gpos[2] > z_low && gpos[2] <= z_high) {  /* in membrane */

                val = B + C*z_rel;
                gzcf[IJKz(i,j,1)] += val;    /* assign high side BC */

            }

            else if (gpos[2] > z_high)  {                  /* extracellular */

                val = D*exp(-xkappa*z_rel);
                gzcf[IJKz(i,j,1)] += val;    /* assign high side BC */
                //printf("%f \n",val);

            }

        }
    }
}

VPRIVATE void bcfl_map(Vpmg *thee){

    Vpbe *pbe;
    double position[3], pot, hx, hy, hzed;
    int i, j, k, nx, ny, nz, rc;


    VASSERT(thee != VNULL);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Reset the potential array */
    for (i=0; i<(nx*ny*nz); i++) thee->pot[i] = 0.0;

    /* Fill in the source term (atomic potentials) */
    Vnm_print(0, "Vpmg_fillco:  filling in source term.\n");
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                position[0] = thee->xf[i];
                position[1] = thee->yf[j];
                position[2] = thee->zf[k];
                rc = Vgrid_value(thee->potMap, position, &pot);
                if (!rc) {
                    Vnm_print(2, "fillcoChargeMap:  Error -- fell off of potential map at (%g, %g, %g)!\n",
                              position[0], position[1], position[2]);
                    VASSERT(0);
                }
                thee->pot[IJK(i,j,k)] = pot;
            }
        }
    }

}

#if  defined(WITH_TINKER)
VPRIVATE void bcfl_mdh_tinker(Vpmg *thee){

    int i,j,k,iatom;
    int nx, ny, nz;

    double val, *apos, gpos[3], tensor[9];
    double *dipole, *quadrupole;
    double size, charge, xkappa, eps_w, eps_p, T, pre1, dist;

    double ux,uy,uz,xr,yr,zr;
    double qxx,qxy,qxz,qyx,qyy,qyz,qzx,qzy,qzz;

    double *xf, *yf, *zf;
    double *gxcf, *gycf, *gzcf;

    Vpbe *pbe;
    Vatom *atom;
    Valist *alist;

    pbe = thee->pbe;
    alist = thee->pbe->alist;
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;

    xf = thee->xf;
    yf = thee->yf;
    zf = thee->zf;

    gxcf = thee->gxcf;
    gycf = thee->gycf;
    gzcf = thee->gzcf;

    /* For each "atom" (only one for bcfl=1), we use the following formula to
     * calculate the boundary conditions:
     *    g(x) = \frac{q e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     *          * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     *          * 1/d
     * where d = ||x - x_0|| (in m) and a is the size of the atom (in m).
     * We only need to evaluate some of these prefactors once:
     *    pre1 = \frac{e_c}{4*\pi*\eps_0*\eps_w*k_b*T}
     * which gives the potential as
     *    g(x) = pre1 * q/d * \frac{exp(-xkappa*(d - a))}{1+xkappa*a}
     */
    eps_w = Vpbe_getSolventDiel(pbe);           /* Dimensionless */
    eps_p = Vpbe_getSoluteDiel(pbe);           /* Dimensionless */
    T = Vpbe_getTemperature(pbe);               /* K             */
    pre1 = (Vunit_ec*Vunit_ec)/(4*VPI*Vunit_eps0*Vunit_kb*T);

    /* Finally, if we convert keep xkappa in A^{-1} and scale pre1 by
     * m/A, then we will only need to deal with distances and sizes in
     * Angstroms rather than meters.                                       */
    xkappa = Vpbe_getXkappa(pbe);              /* A^{-1}        */
    pre1 = pre1*(1.0e10);

    /* Finally, if we convert keep xkappa in A^{-1} and scale pre1 by
     * m/A, then we will only need to deal with distances and sizes in
     * Angstroms rather than meters.                                       */
    xkappa = Vpbe_getXkappa(pbe);              /* A^{-1}        */

    for(k=0;k<nz;k++){
        gpos[2] = zf[k];
        for(j=0;j<ny;j++){
            gpos[1] = yf[j];
            for(i=0;i<nx;i++){
                gpos[0] = xf[i];
                if(gridPointIsValid(i, j, k, nx, ny, nz)){

                    val = 0.0;

                    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
                        atom = Valist_getAtom(alist, iatom);
                        apos = Vatom_getPosition(atom);
                        size = Vatom_getRadius(atom);

                        charge = 0.0;

                        dipole = VNULL;
                        quadrupole = VNULL;

                        if (thee->chargeSrc == VCM_PERMANENT) {
                            charge = Vatom_getCharge(atom);
                            dipole = Vatom_getDipole(atom);
                            quadrupole = Vatom_getQuadrupole(atom);
                        } else if (thee->chargeSrc == VCM_INDUCED) {
                            dipole = Vatom_getInducedDipole(atom);
                        } else {
                            dipole = Vatom_getNLInducedDipole(atom);
                        }

                        ux = dipole[0];
                        uy = dipole[1];
                        uz = dipole[2];

                        if (quadrupole != VNULL) {
                            /* The factor of 1/3 results from using a
                             traceless quadrupole definition. See, for example,
                             "The Theory of Intermolecular Forces" by A.J. Stone,
                             Chapter 3. */
                            qxx = quadrupole[0] / 3.0;
                            qxy = quadrupole[1] / 3.0;
                            qxz = quadrupole[2] / 3.0;
                            qyx = quadrupole[3] / 3.0;
                            qyy = quadrupole[4] / 3.0;
                            qyz = quadrupole[5] / 3.0;
                            qzx = quadrupole[6] / 3.0;
                            qzy = quadrupole[7] / 3.0;
                            qzz = quadrupole[8] / 3.0;
                        } else {
                            qxx = 0.0;
                            qxy = 0.0;
                            qxz = 0.0;
                            qyx = 0.0;
                            qyy = 0.0;
                            qyz = 0.0;
                            qzx = 0.0;
                            qzy = 0.0;
                            qzz = 0.0;
                        }

                        xr = gpos[0] - apos[0];
                        yr = gpos[1] - apos[1];
                        zr = gpos[2] - apos[2];

                        dist = VSQRT(VSQR(xr) + VSQR(yr) + VSQR(zr));
                        multipolebc(dist, xkappa, eps_p, eps_w, size, tensor);

                        val += pre1*charge*tensor[0];
                        val -= pre1*ux*xr*tensor[1];
                        val -= pre1*uy*yr*tensor[1];
                        val -= pre1*uz*zr*tensor[1];
                        val += pre1*qxx*xr*xr*tensor[2];
                        val += pre1*qyy*yr*yr*tensor[2];
                        val += pre1*qzz*zr*zr*tensor[2];
                        val += pre1*2.0*qxy*xr*yr*tensor[2];
                        val += pre1*2.0*qxz*xr*zr*tensor[2];
                        val += pre1*2.0*qyz*yr*zr*tensor[2];

                    }

                    if(i==0){
                        gxcf[IJKx(j,k,0)] = val;
                    }
                    if(i==nx-1){
                        gxcf[IJKx(j,k,1)] = val;
                    }
                    if(j==0){
                        gycf[IJKy(i,k,0)] = val;
                    }
                    if(j==ny-1){
                        gycf[IJKy(i,k,1)] = val;
                    }
                    if(k==0){
                        gzcf[IJKz(i,j,0)] = val;
                    }
                    if(k==nz-1){
                        gzcf[IJKz(i,j,1)] = val;
                    }
                } /* End grid point is valid */
            } /* End i loop */
        } /* End j loop */
    } /* End k loop */

}
#endif

VPRIVATE void bcCalc(Vpmg *thee){

    int i, j, k;
    int nx, ny, nz;

    double zmem, eps_m, Lmem, memv, eps_w, xkappa;

    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;

    /* Zero out the boundaries */
    /* the "i" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            thee->gxcf[IJKx(j,k,0)] = 0.0;
            thee->gxcf[IJKx(j,k,1)] = 0.0;
            thee->gxcf[IJKx(j,k,2)] = 0.0;
            thee->gxcf[IJKx(j,k,3)] = 0.0;
        }
    }

    /* the "j" boundaries (dirichlet) */
    for (k=0; k<nz; k++) {
        for (i=0; i<nx; i++) {
            thee->gycf[IJKy(i,k,0)] = 0.0;
            thee->gycf[IJKy(i,k,1)] = 0.0;
            thee->gycf[IJKy(i,k,2)] = 0.0;
            thee->gycf[IJKy(i,k,3)] = 0.0;
        }
    }

    /* the "k" boundaries (dirichlet) */
    for (j=0; j<ny; j++) {
        for (i=0; i<nx; i++) {
            thee->gzcf[IJKz(i,j,0)] = 0.0;
            thee->gzcf[IJKz(i,j,1)] = 0.0;
            thee->gzcf[IJKz(i,j,2)] = 0.0;
            thee->gzcf[IJKz(i,j,3)] = 0.0;
        }
    }

    switch (thee->pmgp->bcfl) {
            /*  If we have zero boundary conditions, we're done */
        case BCFL_ZERO:
            return;
        case BCFL_SDH:
            bcfl_sdh(thee);
            break;
        case BCFL_MDH:
#if defined(WITH_TINKER)
            bcfl_mdh_tinker(thee);
#else

#ifdef DEBUG_MAC_OSX_OCL
#include "mach_chud.h"
            uint64_t mbeg = mach_absolute_time();

            /*
             * If OpenCL is available we use it, otherwise fall back to
             * normal route (CPU multithreaded w/ OpenMP)
             */
            if (kOpenCLAvailable == 1) bcflnewOpenCL(thee);
            else bcflnew(thee);

            mets_(&mbeg, "MDH");
#else
            /* bcfl_mdh(thee); */
            bcflnew(thee);
#endif	/* DEBUG_MAC_OSX_OCL */

#endif	/* WITH_TINKER */
            break;
        case BCFL_MEM:

            zmem  = Vpbe_getzmem(thee->pbe);
            Lmem  = Vpbe_getLmem(thee->pbe);
            eps_m = Vpbe_getmembraneDiel(thee->pbe);
            memv =  Vpbe_getmemv(thee->pbe);

            eps_w = Vpbe_getSolventDiel(thee->pbe);
            xkappa = Vpbe_getXkappa(thee->pbe);

            bcfl_mem(zmem, Lmem, eps_m, eps_w, memv, xkappa,
                  thee->gxcf, thee->gycf, thee->gzcf,
                  thee->xf, thee->yf, thee->zf, nx, ny, nz);
            break;
        case BCFL_UNUSED:
            Vnm_print(2, "bcCalc:  Invalid bcfl (%d)!\n", thee->pmgp->bcfl);
            VASSERT(0);
            break;
        case BCFL_FOCUS:
            Vnm_print(2, "VPMG::bcCalc -- not appropriate for focusing!\n");
            VASSERT(0);
            break;
        case BCFL_MAP:
            bcfl_map(thee);
            focusFillBound(thee,VNULL);
            break;
        default:
            Vnm_print(2, "VPMG::bcCalc -- invalid boundary condition \
                      flag (%d)!\n", thee->pmgp->bcfl);
            VASSERT(0);
            break;
    }
}

VPRIVATE void fillcoCoefMap(Vpmg *thee) {

    Vpbe *pbe;
    double ionstr, position[3], tkappa, eps, pot, hx, hy, hzed;
    int i, j, k, nx, ny, nz;
    double kappamax;
    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    ionstr = Vpbe_getBulkIonicStrength(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    if ((!thee->useDielXMap) || (!thee->useDielYMap)
        || (!thee->useDielZMap) || ((!thee->useKappaMap) && (ionstr>VPMGSMALL))) {

        Vnm_print(2, "fillcoCoefMap:  You need to use all coefficient maps!\n");
        VASSERT(0);

    }

    /* Scale the kappa map to values between 0 and 1
       Thus get the maximum value in the map - this
       is theoretically unnecessary, but a good check.*/
    kappamax = -1.00;
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                if (ionstr > VPMGSMALL) {
                     position[0] = thee->xf[i];
                     position[1] = thee->yf[j];
                     position[2] = thee->zf[k];
                     if (!Vgrid_value(thee->kappaMap, position, &tkappa)) {
                         Vnm_print(2, "Vpmg_fillco:  Off kappaMap at:\n");
                         Vnm_print(2, "Vpmg_fillco:  (x,y,z) = (%g,%g %g)\n",
                                   position[0], position[1], position[2]);
                         VASSERT(0);
                     }
                     if (tkappa > kappamax) {
                         kappamax = tkappa;
                     }
                     if (tkappa < 0.0){
                       Vnm_print(2, "Vpmg_fillcoCoefMap: Kappa map less than 0\n");
                       Vnm_print(2, "Vpmg_fillcoCoefMap: at (x,y,z) = (%g,%g %g)\n",
                                 position[0], position[1], position[2]);
                       VASSERT(0);
                     }
                }
            }
        }
    }

    if (kappamax > 1.0){
      Vnm_print(2, "Vpmg_fillcoCoefMap:  Maximum Kappa value\n");
      Vnm_print(2, "%g is greater than 1 - will scale appropriately!\n",
                kappamax);
    }
    else {
      kappamax = 1.0;
    }

    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {

                if (ionstr > VPMGSMALL) {
                     position[0] = thee->xf[i];
                     position[1] = thee->yf[j];
                     position[2] = thee->zf[k];
                     if (!Vgrid_value(thee->kappaMap, position, &tkappa)) {
                         Vnm_print(2, "Vpmg_fillco:  Off kappaMap at:\n");
                         Vnm_print(2, "Vpmg_fillco:  (x,y,z) = (%g,%g %g)\n",
                           position[0], position[1], position[2]);
                         VASSERT(0);
                     }
                     if (tkappa < VPMGSMALL) tkappa = 0.0;
                     thee->kappa[IJK(i,j,k)] = (tkappa / kappamax);
                }

                position[0] = thee->xf[i] + 0.5*hx;
                position[1] = thee->yf[j];
                position[2] = thee->zf[k];
                if (!Vgrid_value(thee->dielXMap, position, &eps)) {
                    Vnm_print(2, "Vpmg_fillco:  Off dielXMap at:\n");
                    Vnm_print(2, "Vpmg_fillco:  (x,y,z) = (%g,%g %g)\n",
                      position[0], position[1], position[2]);
                    VASSERT(0);
                 }
                 thee->epsx[IJK(i,j,k)] = eps;

                 position[0] = thee->xf[i];
                 position[1] = thee->yf[j] + 0.5*hy;
                 position[2] = thee->zf[k];
                 if (!Vgrid_value(thee->dielYMap, position, &eps)) {
                    Vnm_print(2, "Vpmg_fillco:  Off dielYMap at:\n");
                    Vnm_print(2, "Vpmg_fillco:  (x,y,z) = (%g,%g %g)\n",
                      position[0], position[1], position[2]);
                    VASSERT(0);
                 }
                 thee->epsy[IJK(i,j,k)] = eps;

                 position[0] = thee->xf[i];
                 position[1] = thee->yf[j];
                 position[2] = thee->zf[k] + 0.5*hzed;
                 if (!Vgrid_value(thee->dielZMap, position, &eps)) {
                    Vnm_print(2, "Vpmg_fillco:  Off dielZMap at:\n");
                    Vnm_print(2, "Vpmg_fillco:  (x,y,z) = (%g,%g %g)\n",
                      position[0], position[1], position[2]);
                    VASSERT(0);
                 }
                 thee->epsz[IJK(i,j,k)] = eps;
            }
        }
    }
}

VPRIVATE void fillcoCoefMol(Vpmg *thee) {

    if (thee->useDielXMap || thee->useDielYMap || thee->useDielZMap ||
      thee->useKappaMap)  {

        fillcoCoefMap(thee);

    } else {

        fillcoCoefMolDiel(thee);
        fillcoCoefMolIon(thee);

    }

}

VPRIVATE void fillcoCoefMolIon(Vpmg *thee) {

    Vacc *acc;
    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax, ionmask, ionstr;
    double xlen, ylen, zlen, irad;
    double hx, hy, hzed, *apos, arad;
    int i, nx, ny, nz, iatom;
    Vsurf_Meth surfMeth;

    VASSERT(thee != VNULL);
    surfMeth = thee->surfMeth;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    ionstr = Vpbe_getBulkIonicStrength(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* This is a floating point parameter related to the non-zero nature of the
     * bulk ionic strength.  If the ionic strength is greater than zero; this
     * parameter is set to 1.0 and later scaled by the appropriate pre-factors.
     * Otherwise, this parameter is set to 0.0 */
    if (ionstr > VPMGSMALL) ionmask = 1.0;
    else ionmask = 0.0;

    /* Reset the kappa array, marking everything accessible */
    for (i=0; i<(nx*ny*nz); i++) thee->kappa[i] = ionmask;

    if (ionstr < VPMGSMALL) return;

    /* Loop through the atoms and set kappa = 0.0 (inaccessible) if a point
     * is inside the ion-inflated van der Waals radii */
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        arad = Vatom_getRadius(atom);

        if (arad > VSMALL) {

            /* Make sure we're on the grid */
            if ((apos[0]<(xmin-irad-arad)) || (apos[0]>(xmax+irad+arad))  || \
                (apos[1]<(ymin-irad-arad)) || (apos[1]>(ymax+irad+arad))  || \
                (apos[2]<(zmin-irad-arad)) || (apos[2]>(zmax+irad+arad))) {
                if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
                    (thee->pmgp->bcfl != BCFL_MAP)) {
                    Vnm_print(2,
    "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n",
                      iatom, apos[0], apos[1], apos[2]);
                    Vnm_print(2, "Vpmg_fillco:  xmin = %g, xmax = %g\n",
                      xmin, xmax);
                    Vnm_print(2, "Vpmg_fillco:  ymin = %g, ymax = %g\n",
                      ymin, ymax);
                    Vnm_print(2, "Vpmg_fillco:  zmin = %g, zmax = %g\n",
                      zmin, zmax);
                }
                fflush(stderr);

            } else { /* if we're on the mesh */

                /* Mark ions */
                markSphere((irad+arad), apos,
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, ymin, zmin,
                        thee->kappa, 0.0);

            } /* endif (on the mesh) */
        }
    } /* endfor (over all atoms) */

}

VPRIVATE void fillcoCoefMolDiel(Vpmg *thee) {

    /* Always call NoSmooth to fill the epsilon arrays */
    fillcoCoefMolDielNoSmooth(thee);

    /* Call the smoothing algorithm as needed */
    if (thee->surfMeth == VSM_MOLSMOOTH) {
        fillcoCoefMolDielSmooth(thee);
    }
}

VPRIVATE void fillcoCoefMolDielNoSmooth(Vpmg *thee) {

    Vacc *acc;
    VaccSurf *asurf;
    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, position[3];
    double srad, epsw, epsp, deps, area;
    double hx, hy, hzed, *apos, arad;
    int i, nx, ny, nz, ntot, iatom, ipt;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    srad = Vpbe_getSolventRadius(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    epsp = Vpbe_getSoluteDiel(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Reset the arrays */
    ntot = nx*ny*nz;
    for (i=0; i<ntot; i++) {
        thee->epsx[i] = epsw;
        thee->epsy[i] = epsw;
        thee->epsz[i] = epsw;
    }

    /* Loop through the atoms and set a{123}cf = 0.0 (inaccessible)
     * if a point is inside the solvent-inflated van der Waals radii */
#pragma omp parallel for default(shared) private(iatom,atom,apos,arad)
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        arad = Vatom_getRadius(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
                (thee->pmgp->bcfl != BCFL_MAP)) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f,\
 %4.3f) is off the mesh (ignoring):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:  xmin = %g, xmax = %g\n",
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:  ymin = %g, ymax = %g\n",
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:  zmin = %g, zmax = %g\n",
                  zmin, zmax);
            }
            fflush(stderr);

        } else { /* if we're on the mesh */

            if (arad > VSMALL) {
                /* Mark x-shifted dielectric */
                markSphere((arad+srad), apos,
                        nx, ny, nz,
                        hx, hy, hzed,
                        (xmin+0.5*hx), ymin, zmin,
                        thee->epsx, epsp);

                /* Mark y-shifted dielectric */
                markSphere((arad+srad), apos,
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, (ymin+0.5*hy), zmin,
                        thee->epsy, epsp);

                /* Mark z-shifted dielectric */
                markSphere((arad+srad), apos,
                        nx, ny, nz,
                        hx, hy, hzed,
                        xmin, ymin, (zmin+0.5*hzed),
                        thee->epsz, epsp);
            }

        } /* endif (on the mesh) */
    } /* endfor (over all atoms) */

    area = Vacc_SASA(acc, srad);

    /* We only need to do the next step for non-zero solvent radii */
    if (srad > VSMALL) {

        /* Now loop over the solvent accessible surface points */

#pragma omp parallel for default(shared) private(iatom,atom,area,asurf,ipt,position)
        for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
            atom = Valist_getAtom(alist, iatom);
            area = Vacc_atomSASA(acc, srad, atom);
            if (area > 0.0 ) {
                asurf = Vacc_atomSASPoints(acc, srad, atom);

                /* Use each point on the SAS to reset the solvent accessibility */
                /* TODO:  Make sure we're not still wasting time here. */
                for (ipt=0; ipt<(asurf->npts); ipt++) {

                    position[0] = asurf->xpts[ipt];
                    position[1] = asurf->ypts[ipt];
                    position[2] = asurf->zpts[ipt];

                    /* Mark x-shifted dielectric */
                    markSphere(srad, position,
                               nx, ny, nz,
                               hx, hy, hzed,
                               (xmin+0.5*hx), ymin, zmin,
                               thee->epsx, epsw);

                    /* Mark y-shifted dielectric */
                    markSphere(srad, position,
                               nx, ny, nz,
                               hx, hy, hzed,
                               xmin, (ymin+0.5*hy), zmin,
                               thee->epsy, epsw);

                    /* Mark z-shifted dielectric */
                    markSphere(srad, position,
                               nx, ny, nz,
                               hx, hy, hzed,
                               xmin, ymin, (zmin+0.5*hzed),
                               thee->epsz, epsw);

                }
            }
        }
    }
}

VPRIVATE void fillcoCoefMolDielSmooth(Vpmg *thee) {

  /* This function smoothes using a 9 point method based on
     Bruccoleri, et al. J Comput Chem 18 268-276 (1997).  The nine points
     used are the shifted grid point and the 8 points that are 1/sqrt(2)
     grid spacings away.  The harmonic mean of the 9 points is then used to
     find the overall dielectric value for the point in question. The use of
     this function assumes that the non-smoothed values were placed in the
     dielectric arrays by the fillcoCoefMolDielNoSmooth function.*/

    Vpbe *pbe;
    double frac, epsw;
    int i, j, k, nx, ny, nz, numpts;

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;

    pbe = thee->pbe;
    epsw = Vpbe_getSolventDiel(pbe);

    /* Copy the existing diel arrays to work arrays */
    for (i=0; i<(nx*ny*nz); i++) {
        thee->a1cf[i] = thee->epsx[i];
        thee->a2cf[i] = thee->epsy[i];
        thee->a3cf[i] = thee->epsz[i];
        thee->epsx[i] = epsw;
        thee->epsy[i] = epsw;
        thee->epsz[i] = epsw;
    }

    /* Smooth the dielectric values */
    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++) {
            for (k=0; k<nz; k++) {

                /* Get the 8 points that are 1/sqrt(2) grid spacings away */

                /* Points for the X-shifted array */
                frac = 1.0/thee->a1cf[IJK(i,j,k)];
                frac += 1.0/thee->a2cf[IJK(i,j,k)];
                frac += 1.0/thee->a3cf[IJK(i,j,k)];
                numpts = 3;

                if (j > 0) {
                    frac += 1.0/thee->a2cf[IJK(i,j-1,k)];
                    numpts += 1;
                }
                if (k > 0) {
                    frac += 1.0/thee->a3cf[IJK(i,j,k-1)];
                    numpts += 1;
                }
                if (i < (nx-1)){
                    frac += 1.0/thee->a2cf[IJK(i+1,j,k)];
                    frac += 1.0/thee->a3cf[IJK(i+1,j,k)];
                    numpts += 2;
                    if (j > 0) {
                        frac += 1.0/thee->a2cf[IJK(i+1,j-1,k)];
                        numpts += 1;
                    }
                    if (k > 0) {
                        frac += 1.0/thee->a3cf[IJK(i+1,j,k-1)];
                        numpts += 1;
                    }
                }
                thee->epsx[IJK(i,j,k)] = numpts/frac;

                /* Points for the Y-shifted array */
                frac = 1.0/thee->a2cf[IJK(i,j,k)];
                frac += 1.0/thee->a1cf[IJK(i,j,k)];
                frac += 1.0/thee->a3cf[IJK(i,j,k)];
                numpts = 3;

                if (i > 0) {
                    frac += 1.0/thee->a1cf[IJK(i-1,j,k)];
                    numpts += 1;
                }
                if (k > 0) {
                    frac += 1.0/thee->a3cf[IJK(i,j,k-1)];
                    numpts += 1;
                }
                if (j < (ny-1)){
                    frac += 1.0/thee->a1cf[IJK(i,j+1,k)];
                    frac += 1.0/thee->a3cf[IJK(i,j+1,k)];
                    numpts += 2;
                    if (i > 0) {
                        frac += 1.0/thee->a1cf[IJK(i-1,j+1,k)];
                        numpts += 1;
                    }
                    if (k > 0) {
                        frac += 1.0/thee->a3cf[IJK(i,j+1,k-1)];
                        numpts += 1;
                    }
                }
                thee->epsy[IJK(i,j,k)] = numpts/frac;

                /* Points for the Z-shifted array */
                frac = 1.0/thee->a3cf[IJK(i,j,k)];
                frac += 1.0/thee->a1cf[IJK(i,j,k)];
                frac += 1.0/thee->a2cf[IJK(i,j,k)];
                numpts = 3;

                if (i > 0) {
                    frac += 1.0/thee->a1cf[IJK(i-1,j,k)];
                    numpts += 1;
                }
                if (j > 0) {
                    frac += 1.0/thee->a2cf[IJK(i,j-1,k)];
                    numpts += 1;
                }
                if (k < (nz-1)){
                    frac += 1.0/thee->a1cf[IJK(i,j,k+1)];
                    frac += 1.0/thee->a2cf[IJK(i,j,k+1)];
                    numpts += 2;
                    if (i > 0) {
                        frac += 1.0/thee->a1cf[IJK(i-1,j,k+1)];
                        numpts += 1;
                    }
                    if (j > 0) {
                        frac += 1.0/thee->a2cf[IJK(i,j-1,k+1)];
                        numpts += 1;
                    }
                }
                thee->epsz[IJK(i,j,k)] = numpts/frac;
            }
        }
    }
}


VPRIVATE void fillcoCoefSpline(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax, ionmask, ionstr, dist2;
    double xlen, ylen, zlen, position[3], itot, stot, ictot, ictot2, sctot;
    double irad, dx, dy, dz, epsw, epsp, w2i;
    double hx, hy, hzed, *apos, arad, sctot2;
    double dx2, dy2, dz2, stot2, itot2, rtot, rtot2, splineWin, w3i;
    double dist, value, sm, sm2;
    int i, j, k, nx, ny, nz, iatom;
    int imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);
    splineWin = thee->splineWin;
    w2i = 1.0/(splineWin*splineWin);
    w3i = 1.0/(splineWin*splineWin*splineWin);

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    ionstr = Vpbe_getBulkIonicStrength(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    epsp = Vpbe_getSoluteDiel(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* This is a floating point parameter related to the non-zero nature of the
     * bulk ionic strength.  If the ionic strength is greater than zero; this
     * parameter is set to 1.0 and later scaled by the appropriate pre-factors.
     * Otherwise, this parameter is set to 0.0 */
    if (ionstr > VPMGSMALL) ionmask = 1.0;
    else ionmask = 0.0;

    /* Reset the kappa, epsx, epsy, and epsz arrays */
    for (i=0; i<(nx*ny*nz); i++) {
        thee->kappa[i] = 1.0;
        thee->epsx[i] = 1.0;
        thee->epsy[i] = 1.0;
        thee->epsz[i] = 1.0;
    }

    /* Loop through the atoms and do assign the dielectric */
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        arad = Vatom_getRadius(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
                (thee->pmgp->bcfl != BCFL_MAP)) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f,\
 %4.3f) is off the mesh (ignoring):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:    xmin = %g, xmax = %g\n",
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:    ymin = %g, ymax = %g\n",
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:    zmin = %g, zmax = %g\n",
                  zmin, zmax);
            }
            fflush(stderr);

        } else if (arad > VPMGSMALL ) { /* if we're on the mesh */

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* MARK ION ACCESSIBILITY AND DIELECTRIC VALUES FOR LATER
             * ASSIGNMENT (Steps #1-3) */
            itot = irad + arad + splineWin;
            itot2 = VSQR(itot);
            ictot = VMAX2(0, (irad + arad - splineWin));
            ictot2 = VSQR(ictot);
            stot = arad + splineWin;
            stot2 = VSQR(stot);
            sctot = VMAX2(0, (arad - splineWin));
            sctot2 = VSQR(sctot);

           /* We'll search over grid points which are in the greater of
             * these two radii */
            rtot = VMAX2(itot, stot);
            rtot2 = VMAX2(itot2, stot2);
            dx = rtot + 0.5*hx;
            dy = rtot + 0.5*hy;
            dz = rtot + 0.5*hzed;
            imin = VMAX2(0,(int)floor((position[0] - dx)/hx));
            imax = VMIN2(nx-1,(int)ceil((position[0] + dx)/hx));
            jmin = VMAX2(0,(int)floor((position[1] - dy)/hy));
            jmax = VMIN2(ny-1,(int)ceil((position[1] + dy)/hy));
            kmin = VMAX2(0,(int)floor((position[2] - dz)/hzed));
            kmax = VMIN2(nz-1,(int)ceil((position[2] + dz)/hzed));
            for (i=imin; i<=imax; i++) {
                dx2 = VSQR(position[0] - hx*i);
                for (j=jmin; j<=jmax; j++) {
                    dy2 = VSQR(position[1] - hy*j);
                    for (k=kmin; k<=kmax; k++) {
                        dz2 = VSQR(position[2] - k*hzed);

                        /* ASSIGN CCF */
                        if (thee->kappa[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2 + dy2 + dx2;
                            if (dist2 >= itot2) {
                                ;
                            }
                            if (dist2 <= ictot2) {
                                thee->kappa[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 < itot2) && (dist2 > ictot2)) {
                                dist = VSQRT(dist2);
                                sm = dist - (arad + irad) + splineWin;
                                sm2 = VSQR(sm);
                                value = 0.75*sm2*w2i - 0.25*sm*sm2*w3i;
                                thee->kappa[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A1CF */
                        if (thee->epsx[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2+dy2+VSQR(position[0]-(i+0.5)*hx);
                            if (dist2 >= stot2) {
                                thee->epsx[IJK(i,j,k)] *= 1.0;
                            }
                            if (dist2 <= sctot2) {
                                thee->epsx[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist - arad + splineWin;
                                sm2 = VSQR(sm);
                                value = 0.75*sm2*w2i - 0.25*sm*sm2*w3i;
                                thee->epsx[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A2CF */
                        if (thee->epsy[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2+dx2+VSQR(position[1]-(j+0.5)*hy);
                            if (dist2 >= stot2) {
                                thee->epsy[IJK(i,j,k)] *= 1.0;
                            }
                            if (dist2 <= sctot2) {
                                thee->epsy[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist - arad + splineWin;
                                sm2 = VSQR(sm);
                                value = 0.75*sm2*w2i - 0.25*sm*sm2*w3i;
                                thee->epsy[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A3CF */
                        if (thee->epsz[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dy2+dx2+VSQR(position[2]-(k+0.5)*hzed);
                            if (dist2 >= stot2) {
                                thee->epsz[IJK(i,j,k)] *= 1.0;
                            }
                            if (dist2 <= sctot2) {
                                thee->epsz[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist - arad + splineWin;
                                sm2 = VSQR(sm);
                                value = 0.75*sm2*w2i - 0.25*sm*sm2*w3i;
                                thee->epsz[IJK(i,j,k)] *= value;
                            }
                        }


                    } /* k loop */
                } /* j loop */
            } /* i loop */
        } /* endif (on the mesh) */
    } /* endfor (over all atoms) */

    Vnm_print(0, "Vpmg_fillco:  filling coefficient arrays\n");
    /* Interpret markings and fill the coefficient arrays */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {

                thee->kappa[IJK(i,j,k)] = ionmask*thee->kappa[IJK(i,j,k)];
                thee->epsx[IJK(i,j,k)] = (epsw-epsp)*thee->epsx[IJK(i,j,k)]
                  + epsp;
                thee->epsy[IJK(i,j,k)] = (epsw-epsp)*thee->epsy[IJK(i,j,k)]
                  + epsp;
                thee->epsz[IJK(i,j,k)] = (epsw-epsp)*thee->epsz[IJK(i,j,k)]
                  + epsp;

            } /* i loop */
        } /* j loop */
    } /* k loop */

}

VPRIVATE void fillcoCoef(Vpmg *thee) {

    VASSERT(thee != VNULL);

    if (thee->useDielXMap || thee->useDielYMap ||
        thee->useDielZMap || thee->useKappaMap) {
        fillcoCoefMap(thee);
        return;
    }

    switch(thee->surfMeth) {
        case VSM_MOL:
            Vnm_print(0, "fillcoCoef:  Calling fillcoCoefMol...\n");
            fillcoCoefMol(thee);
            break;
        case VSM_MOLSMOOTH:
            Vnm_print(0, "fillcoCoef:  Calling fillcoCoefMol...\n");
            fillcoCoefMol(thee);
            break;
        case VSM_SPLINE:
            Vnm_print(0, "fillcoCoef:  Calling fillcoCoefSpline...\n");
            fillcoCoefSpline(thee);
            break;
        case VSM_SPLINE3:
            Vnm_print(0, "fillcoCoef:  Calling fillcoCoefSpline3...\n");
            fillcoCoefSpline3(thee);
            break;
        case VSM_SPLINE4:
            Vnm_print(0, "fillcoCoef:  Calling fillcoCoefSpline4...\n");
            fillcoCoefSpline4(thee);
            break;
        default:
            Vnm_print(2, "fillcoCoef:  Invalid surfMeth (%d)!\n",
              thee->surfMeth);
            VASSERT(0);
            break;
    }
}


VPRIVATE Vrc_Codes fillcoCharge(Vpmg *thee) {

    Vrc_Codes rc;

    VASSERT(thee != VNULL);

    if (thee->useChargeMap) {
        return fillcoChargeMap(thee);
    }

    switch(thee->chargeMeth) {
        case VCM_TRIL:
            Vnm_print(0, "fillcoCharge:  Calling fillcoChargeSpline1...\n");
            fillcoChargeSpline1(thee);
            break;
        case VCM_BSPL2:
            Vnm_print(0, "fillcoCharge:  Calling fillcoChargeSpline2...\n");
            fillcoChargeSpline2(thee);
            break;
        case VCM_BSPL4:
            switch (thee->chargeSrc) {
                case VCM_CHARGE:
                    Vnm_print(0, "fillcoCharge: Calling fillcoPermanentMultipole...\n");
                    fillcoPermanentMultipole(thee);
                    break;
#if defined(WITH_TINKER)
                case VCM_PERMANENT:
                    Vnm_print(0, "fillcoCharge: Calling fillcoPermanentMultipole...\n");
                    fillcoPermanentMultipole(thee);
                    break;
                case VCM_INDUCED:
                    Vnm_print(0, "fillcoCharge: Calling fillcoInducedDipole...\n");
                    fillcoInducedDipole(thee);
                    break;
                case VCM_NLINDUCED:
                     Vnm_print(0, "fillcoCharge: Calling fillcoNLInducedDipole...\n");
                     fillcoNLInducedDipole(thee);
                     break;
#endif /* if defined(WITH_TINKER) */
                default:
                    Vnm_print(2, "fillcoCharge:  Invalid chargeSource (%d)!\n",
                      thee->chargeSrc);
                    return VRC_FAILURE;
                    break;
            }
            break;
        default:
            Vnm_print(2, "fillcoCharge:  Invalid chargeMeth (%d)!\n",
              thee->chargeMeth);
            return VRC_FAILURE;
            break;
    }

    return VRC_SUCCESS;
}

VPRIVATE Vrc_Codes fillcoChargeMap(Vpmg *thee) {

    Vpbe *pbe;
    double position[3], charge, zmagic, hx, hy, hzed;
    int i, j, k, nx, ny, nz, rc;


    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    zmagic = Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Reset the charge array */
    for (i=0; i<(nx*ny*nz); i++) thee->charge[i] = 0.0;

    /* Fill in the source term (atomic charges) */
    Vnm_print(0, "Vpmg_fillco:  filling in source term.\n");
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {
                position[0] = thee->xf[i];
                position[1] = thee->yf[j];
                position[2] = thee->zf[k];
                rc = Vgrid_value(thee->chargeMap, position, &charge);
                if (!rc) {
                    Vnm_print(2, "fillcoChargeMap:  Error -- fell off of charge map at (%g, %g, %g)!\n",
                          position[0], position[1], position[2]);
                    return VRC_FAILURE;
                }
                /* Scale the charge to internal units */
                charge = charge*zmagic;
                thee->charge[IJK(i,j,k)] = charge;
            }
        }
    }

    return VRC_SUCCESS;
}

VPRIVATE void fillcoChargeSpline1(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, position[3], ifloat, jfloat, kfloat;
    double charge, dx, dy, dz, zmagic, hx, hy, hzed, *apos;
    int i, nx, ny, nz, iatom, ihi, ilo, jhi, jlo, khi, klo;


    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    zmagic = Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Reset the charge array */
    for (i=0; i<(nx*ny*nz); i++) thee->charge[i] = 0.0;

    /* Fill in the source term (atomic charges) */
    Vnm_print(0, "Vpmg_fillco:  filling in source term.\n");
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        charge = Vatom_getCharge(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
                (thee->pmgp->bcfl != BCFL_MAP)) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, \
%4.3f) is off the mesh (ignoring):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:    xmin = %g, xmax = %g\n",
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:    ymin = %g, ymax = %g\n",
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:    zmin = %g, zmax = %g\n",
                  zmin, zmax);
            }
            fflush(stderr);
        } else {

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* Scale the charge to be a delta function */
            charge = charge*zmagic/(hx*hy*hzed);

            /* Figure out which vertices we're next to */
            ifloat = position[0]/hx;
            jfloat = position[1]/hy;
            kfloat = position[2]/hzed;

            ihi = (int)ceil(ifloat);
            ilo = (int)floor(ifloat);
            jhi = (int)ceil(jfloat);
            jlo = (int)floor(jfloat);
            khi = (int)ceil(kfloat);
            klo = (int)floor(kfloat);

            /* Now assign fractions of the charge to the nearby verts */
            dx = ifloat - (double)(ilo);
            dy = jfloat - (double)(jlo);
            dz = kfloat - (double)(klo);
            thee->charge[IJK(ihi,jhi,khi)] += (dx*dy*dz*charge);
            thee->charge[IJK(ihi,jlo,khi)] += (dx*(1.0-dy)*dz*charge);
            thee->charge[IJK(ihi,jhi,klo)] += (dx*dy*(1.0-dz)*charge);
            thee->charge[IJK(ihi,jlo,klo)] += (dx*(1.0-dy)*(1.0-dz)*charge);
            thee->charge[IJK(ilo,jhi,khi)] += ((1.0-dx)*dy*dz *charge);
            thee->charge[IJK(ilo,jlo,khi)] += ((1.0-dx)*(1.0-dy)*dz *charge);
            thee->charge[IJK(ilo,jhi,klo)] += ((1.0-dx)*dy*(1.0-dz)*charge);
            thee->charge[IJK(ilo,jlo,klo)] += ((1.0-dx)*(1.0-dy)*(1.0-dz)*charge);
        } /* endif (on the mesh) */
    } /* endfor (each atom) */
}

VPRIVATE double bspline2(double x) {

    double m2m, m2, m3;

    if ((x >= 0.0) && (x <= 2.0)) m2m = 1.0 - VABS(x - 1.0);
    else m2m = 0.0;
    if ((x >= 1.0) && (x <= 3.0)) m2 = 1.0 - VABS(x - 2.0);
    else m2 = 0.0;

    if ((x >= 0.0) && (x <= 3.0)) m3 = 0.5*x*m2m + 0.5*(3.0-x)*m2;
    else m3 = 0.0;

    return m3;

}

VPRIVATE double dbspline2(double x) {

    double m2m, m2, dm3;

    if ((x >= 0.0) && (x <= 2.0)) m2m = 1.0 - VABS(x - 1.0);
    else m2m = 0.0;
    if ((x >= 1.0) && (x <= 3.0)) m2 = 1.0 - VABS(x - 2.0);
    else m2 = 0.0;

    dm3 = m2m - m2;

    return dm3;

}


VPRIVATE void fillcoChargeSpline2(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax, zmagic;
    double xlen, ylen, zlen, position[3], ifloat, jfloat, kfloat;
    double charge, hx, hy, hzed, *apos, mx, my, mz;
    int i, ii, jj, kk, nx, ny, nz, iatom;
    int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1, kp1, kp2;


    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    zmagic = Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Reset the charge array */
    for (i=0; i<(nx*ny*nz); i++) thee->charge[i] = 0.0;

    /* Fill in the source term (atomic charges) */
    Vnm_print(0, "Vpmg_fillco:  filling in source term.\n");
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        charge = Vatom_getCharge(atom);

        /* Make sure we're on the grid */
        if ((apos[0]<=(xmin-hx)) || (apos[0]>=(xmax+hx))  || \
            (apos[1]<=(ymin-hy)) || (apos[1]>=(ymax+hy))  || \
            (apos[2]<=(zmin-hzed)) || (apos[2]>=(zmax+hzed))) {
            if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
                (thee->pmgp->bcfl != BCFL_MAP)) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f, \
%4.3f) is off the mesh (for cubic splines!!) (ignoring this atom):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:    xmin = %g, xmax = %g\n",
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:    ymin = %g, ymax = %g\n",
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:    zmin = %g, zmax = %g\n",
                  zmin, zmax);
            }
            fflush(stderr);
        } else {

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* Scale the charge to be a delta function */
            charge = charge*zmagic/(hx*hy*hzed);

            /* Figure out which vertices we're next to */
            ifloat = position[0]/hx;
            jfloat = position[1]/hy;
            kfloat = position[2]/hzed;

            ip1   = (int)ceil(ifloat);
            ip2   = ip1 + 1;
            im1   = (int)floor(ifloat);
            im2   = im1 - 1;
            jp1   = (int)ceil(jfloat);
            jp2   = jp1 + 1;
            jm1   = (int)floor(jfloat);
            jm2   = jm1 - 1;
            kp1   = (int)ceil(kfloat);
            kp2   = kp1 + 1;
            km1   = (int)floor(kfloat);
            km2   = km1 - 1;

            /* This step shouldn't be necessary, but it saves nasty debugging
             * later on if something goes wrong */
            ip2 = VMIN2(ip2,nx-1);
            ip1 = VMIN2(ip1,nx-1);
            im1 = VMAX2(im1,0);
            im2 = VMAX2(im2,0);
            jp2 = VMIN2(jp2,ny-1);
            jp1 = VMIN2(jp1,ny-1);
            jm1 = VMAX2(jm1,0);
            jm2 = VMAX2(jm2,0);
            kp2 = VMIN2(kp2,nz-1);
            kp1 = VMIN2(kp1,nz-1);
            km1 = VMAX2(km1,0);
            km2 = VMAX2(km2,0);

            /* Now assign fractions of the charge to the nearby verts */
            for (ii=im2; ii<=ip2; ii++) {
                mx = bspline2(VFCHI(ii,ifloat));
                for (jj=jm2; jj<=jp2; jj++) {
                    my = bspline2(VFCHI(jj,jfloat));
                    for (kk=km2; kk<=kp2; kk++) {
                        mz = bspline2(VFCHI(kk,kfloat));
                        thee->charge[IJK(ii,jj,kk)] += (charge*mx*my*mz);
                    }
                }
            }

        } /* endif (on the mesh) */
    } /* endfor (each atom) */
}

VPUBLIC int Vpmg_fillco(Vpmg *thee,
                        Vsurf_Meth surfMeth,
                        double splineWin,
                        Vchrg_Meth chargeMeth,
                        int useDielXMap,
                        Vgrid *dielXMap,
                        int useDielYMap,
                        Vgrid *dielYMap,
                        int useDielZMap,
                        Vgrid *dielZMap,
                        int useKappaMap,
                        Vgrid *kappaMap,
                        int usePotMap,
                        Vgrid *potMap,
                        int useChargeMap,
                        Vgrid *chargeMap
                       ) {

    Vpbe *pbe;
    double xmin,
           xmax,
           ymin,
           ymax,
           zmin,
           zmax,
           xlen,
           ylen,
           zlen,
           hx,
           hy,
           hzed,
           epsw,
           epsp,
           ionstr;
    int i,
        nx,
        ny,
        nz,
        islap;
    Vrc_Codes rc;

    if (thee == VNULL) {
        Vnm_print(2, "Vpmg_fillco:  got NULL thee!\n");
        return 0;
    }

    thee->surfMeth = surfMeth;
    thee->splineWin = splineWin;
    thee->chargeMeth = chargeMeth;
    thee->useDielXMap = useDielXMap;
    if (thee->useDielXMap) thee->dielXMap = dielXMap;
    thee->useDielYMap = useDielYMap;
    if (thee->useDielYMap) thee->dielYMap = dielYMap;
    thee->useDielZMap = useDielZMap;
    if (thee->useDielZMap) thee->dielZMap = dielZMap;
    thee->useKappaMap = useKappaMap;
    if (thee->useKappaMap) thee->kappaMap = kappaMap;
    thee->usePotMap = usePotMap;
    if (thee->usePotMap) thee->potMap = potMap;
    thee->useChargeMap = useChargeMap;
    if (thee->useChargeMap) thee->chargeMap = chargeMap;

    /* Get PBE info */
    pbe = thee->pbe;
    ionstr = Vpbe_getBulkIonicStrength(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    epsp = Vpbe_getSoluteDiel(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    thee->pmgp->xmin = xmin;
    ymin = thee->pmgp->ycent - (ylen/2.0);
    thee->pmgp->ymin = ymin;
    zmin = thee->pmgp->zcent - (zlen/2.0);
    thee->pmgp->zmin = zmin;
    xmax = thee->pmgp->xcent + (xlen/2.0);
    thee->pmgp->xmax = xmax;
    ymax = thee->pmgp->ycent + (ylen/2.0);
    thee->pmgp->ymax = ymax;
    zmax = thee->pmgp->zcent + (zlen/2.0);
    thee->pmgp->zmax = zmax;
    thee->rparm[2] = xmin;
    thee->rparm[3] = xmax;
    thee->rparm[4] = ymin;
    thee->rparm[5] = ymax;
    thee->rparm[6] = zmin;
    thee->rparm[7] = zmax;

    /* This is a flag that gets set if the operator is a simple Laplacian;
     * i.e., in the case of a homogenous dielectric and zero ionic strength
     * The operator cannot be a simple Laplacian if maps are read in. */
    if(thee->useDielXMap || thee->useDielYMap || thee->useDielZMap ||
       thee->useKappaMap || thee->usePotMap){
        islap = 0;
    }else if ( (ionstr < VPMGSMALL) && (VABS(epsp-epsw) < VPMGSMALL) ){
        islap = 1;
    }else{
        islap = 0;
    }

    /* Fill the mesh point coordinate arrays */
    for (i=0; i<nx; i++) thee->xf[i] = xmin + i*hx;
    for (i=0; i<ny; i++) thee->yf[i] = ymin + i*hy;
    for (i=0; i<nz; i++) thee->zf[i] = zmin + i*hzed;

    /* Reset the tcf array */
    for (i=0; i<(nx*ny*nz); i++) thee->tcf[i] = 0.0;

    /* Fill in the source term (atomic charges) */
    Vnm_print(0, "Vpmg_fillco:  filling in source term.\n");
    rc = fillcoCharge(thee);
    switch(rc) {
        case VRC_SUCCESS:
            break;
        case VRC_WARNING:
            Vnm_print(2, "Vpmg_fillco:  non-fatal errors while filling charge map!\n");
            break;
        case VRC_FAILURE:
            Vnm_print(2, "Vpmg_fillco:  fatal errors while filling charge map!\n");
            return 0;
            break;
    }

    /* THE FOLLOWING NEEDS TO BE DONE IF WE'RE NOT USING A SIMPLE LAPLACIAN
     * OPERATOR */
    if (!islap) {
        Vnm_print(0, "Vpmg_fillco:  marking ion and solvent accessibility.\n");
        fillcoCoef(thee);
        Vnm_print(0, "Vpmg_fillco:  done filling coefficient arrays\n");

    } else { /* else (!islap) ==> It's a Laplacian operator! */

        for (i=0; i<(nx*ny*nz); i++) {
            thee->kappa[i] = 0.0;
            thee->epsx[i] = epsp;
            thee->epsy[i] = epsp;
            thee->epsz[i] = epsp;
        }

    } /* endif (!islap) */

    /* Fill the boundary arrays (except when focusing, bcfl = 4) */
    if (thee->pmgp->bcfl != BCFL_FOCUS) {
        Vnm_print(0, "Vpmg_fillco:  filling boundary arrays\n");
        bcCalc(thee);
        Vnm_print(0, "Vpmg_fillco:  done filling boundary arrays\n");
    }

    thee->filled = 1;

    return 1;
}


VPUBLIC int Vpmg_force(Vpmg *thee, double *force, int atomID,
  Vsurf_Meth srfm, Vchrg_Meth chgm) {

    int rc = 1;
    double qfF[3];                  /* Charge-field force */
    double dbF[3];                  /* Dielectric boundary force */
    double ibF[3];                  /* Ion boundary force */
    double npF[3];                  /* Non-polar boundary force */

    VASSERT(thee != VNULL);

    rc = rc && Vpmg_dbForce(thee, qfF, atomID, srfm);
    rc = rc && Vpmg_ibForce(thee, dbF, atomID, srfm);
    rc = rc && Vpmg_qfForce(thee, ibF, atomID, chgm);

    force[0] = qfF[0] + dbF[0] + ibF[0];
    force[1] = qfF[1] + dbF[1] + ibF[1];
    force[2] = qfF[2] + dbF[2] + ibF[2];

    return rc;

}

VPUBLIC int Vpmg_ibForce(Vpmg *thee, double *force, int atomID,
  Vsurf_Meth srfm) {

    Valist *alist;
    Vacc *acc;
    Vpbe *pbe;
    Vatom *atom;

    double *apos, position[3], arad, irad, zkappa2, hx, hy, hzed;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2;
    double rtot, dx, dx2, dy, dy2, dz, dz2, gpos[3], tgrad[3], fmag;
    double izmagic;
    int i, j, k, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;

    /* For nonlinear forces */
    int ichop, nchop, nion, m;
    double ionConc[MAXION], ionRadii[MAXION], ionQ[MAXION], ionstr;

    VASSERT(thee != VNULL);

    acc = thee->pbe->acc;
    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    arad = Vatom_getRadius(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* Check surface definition */
    if ((srfm != VSM_SPLINE) && (srfm!=VSM_SPLINE3) && (srfm!=VSM_SPLINE4)) {
        Vnm_print(2, "Vpmg_ibForce:  Forces *must* be calculated with \
spline-based surfaces!\n");
        Vnm_print(2, "Vpmg_ibForce:  Skipping ionic boundary force \
calculation!\n");
        return 0;
    }

    /* If we aren't in the current position, then we're done */
    if (atom->partID == 0) return 1;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    zkappa2 = Vpbe_getZkappa2(pbe);
    izmagic = 1.0/Vpbe_getZmagic(pbe);

    ionstr = Vpbe_getBulkIonicStrength(pbe);
    Vpbe_getIons(pbe, &nion, ionConc, ionRadii, ionQ);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;

    /* Sanity check: there is no force if there is zero ionic strength */
    if (zkappa2 < VPMGSMALL) {
#ifndef VAPBSQUIET
        Vnm_print(2, "Vpmg_ibForce:  No force for zero ionic strength!\n");
#endif
        return 1;
    }

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
        if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
            (thee->pmgp->bcfl != BCFL_MAP)) {
            Vnm_print(2, "Vpmg_ibForce:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n",
                  atom, apos[0], apos[1], apos[2]);
            Vnm_print(2, "Vpmg_ibForce:    xmin = %g, xmax = %g\n",
              xmin, xmax);
            Vnm_print(2, "Vpmg_ibForce:    ymin = %g, ymax = %g\n",
              ymin, ymax);
            Vnm_print(2, "Vpmg_ibForce:    zmin = %g, zmax = %g\n",
              zmin, zmax);
        }
        fflush(stderr);
    } else {

        /* Convert the atom position to grid reference frame */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;

        /* Integrate over points within this atom's (inflated) radius */
        rtot = (irad + arad + thee->splineWin);
        rtot2 = VSQR(rtot);
        dx = rtot + 0.5*hx;
        imin = VMAX2(0,(int)ceil((position[0] - dx)/hx));
        imax = VMIN2(nx-1,(int)floor((position[0] + dx)/hx));
        for (i=imin; i<=imax; i++) {
            dx2 = VSQR(position[0] - hx*i);
            if (rtot2 > dx2) dy = VSQRT(rtot2 - dx2) + 0.5*hy;
            else dy = 0.5*hy;
            jmin = VMAX2(0,(int)ceil((position[1] - dy)/hy));
            jmax = VMIN2(ny-1,(int)floor((position[1] + dy)/hy));
            for (j=jmin; j<=jmax; j++) {
                dy2 = VSQR(position[1] - hy*j);
                if (rtot2 > (dx2+dy2)) dz = VSQRT(rtot2-dx2-dy2)+0.5*hzed;
                else dz = 0.5*hzed;
                kmin = VMAX2(0,(int)ceil((position[2] - dz)/hzed));
                kmax = VMIN2(nz-1,(int)floor((position[2] + dz)/hzed));
                for (k=kmin; k<=kmax; k++) {
                    dz2 = VSQR(k*hzed - position[2]);
                    /* See if grid point is inside ivdw radius and set kappa
                     * accordingly (do spline assignment here) */
                    if ((dz2 + dy2 + dx2) <= rtot2) {
                        gpos[0] = i*hx + xmin;
                        gpos[1] = j*hy + ymin;
                        gpos[2] = k*hzed + zmin;

                        /* Select the correct function based on the surface definition
                         *	(now including the 7th order polynomial) */
                        Vpmg_splineSelect(srfm,acc, gpos,thee->splineWin, irad, atom, tgrad);

                        if (thee->pmgp->nonlin) {
                            /* Nonlinear forces */
                            fmag = 0.0;
                            nchop = 0;
                            for (m=0; m<nion; m++) {
                                fmag += (thee->kappa[IJK(i,j,k)])*ionConc[m]*(Vcap_exp(-ionQ[m]*thee->u[IJK(i,j,k)], &ichop)-1.0)/ionstr;
                                nchop += ichop;
                            }
                            /*          if (nchop > 0) Vnm_print(2, "Vpmg_ibForece:  Chopped EXP %d times!\n", nchop);*/
                            force[0] += (zkappa2*fmag*tgrad[0]);
                            force[1] += (zkappa2*fmag*tgrad[1]);
                            force[2] += (zkappa2*fmag*tgrad[2]);
                        } else {
                            /* Use of bulk factor (zkappa2) OK here becuase
                             * LPBE force approximation */
                            /* NAB -- did we forget a kappa factor here??? */
                            fmag = VSQR(thee->u[IJK(i,j,k)])*(thee->kappa[IJK(i,j,k)]);
                            force[0] += (zkappa2*fmag*tgrad[0]);
                            force[1] += (zkappa2*fmag*tgrad[1]);
                            force[2] += (zkappa2*fmag*tgrad[2]);
                        }
                    }
                } /* k loop */
            } /* j loop */
        } /* i loop */
    }
    force[0] = force[0] * 0.5 * hx * hy * hzed * izmagic;
    force[1] = force[1] * 0.5 * hx * hy * hzed * izmagic;
    force[2] = force[2] * 0.5 * hx * hy * hzed * izmagic;

    return 1;
}

VPUBLIC int Vpmg_dbForce(Vpmg *thee, double *dbForce, int atomID,
                         Vsurf_Meth srfm) {

    Vacc *acc;
    Vpbe *pbe;
    Vatom *atom;

    double *apos, position[3], arad, srad, hx, hy, hzed, izmagic, deps, depsi;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2, epsp;
    double rtot, dx, gpos[3], tgrad[3], dbFmag, epsw, kT;
    double *u, Hxijk, Hyijk, Hzijk, Hxim1jk, Hyijm1k, Hzijkm1;
    double dHxijk[3], dHyijk[3], dHzijk[3], dHxim1jk[3], dHyijm1k[3];
    double dHzijkm1[3];
    int i, j, k, l, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);
    if (!thee->filled) {
        Vnm_print(2, "Vpmg_dbForce:  Need to callVpmg_fillco!\n");
        return 0;
    }

    acc = thee->pbe->acc;
    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    arad = Vatom_getRadius(atom);
    srad = Vpbe_getSolventRadius(thee->pbe);

    /* Reset force */
    dbForce[0] = 0.0;
    dbForce[1] = 0.0;
    dbForce[2] = 0.0;

    /* Check surface definition */
    if ((srfm != VSM_SPLINE) && (srfm!=VSM_SPLINE3) && (srfm!=VSM_SPLINE4)) {
        Vnm_print(2, "Vpmg_dbForce:  Forces *must* be calculated with \
spline-based surfaces!\n");
        Vnm_print(2, "Vpmg_dbForce:  Skipping dielectric/apolar boundary \
force calculation!\n");
        return 0;
    }


    /* If we aren't in the current position, then we're done */
    if (atom->partID == 0) return 1;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    epsp = Vpbe_getSoluteDiel(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    kT = Vpbe_getTemperature(pbe)*(1e-3)*Vunit_Na*Vunit_kb;
    izmagic = 1.0/Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = thee->u;

    /* Sanity check: there is no force if there is zero ionic strength */
    if (VABS(epsp-epsw) < VPMGSMALL) {
        Vnm_print(0, "Vpmg_dbForce: No force for uniform dielectric!\n");
        return 1;
    }
    deps = (epsw - epsp);
    depsi = 1.0/deps;
    rtot = (arad + thee->splineWin + srad);

    /* Make sure we're on the grid */
    /* Grid checking modified by Matteo Rotter */
    if ((apos[0]<=xmin + rtot) || (apos[0]>=xmax - rtot)  || \
        (apos[1]<=ymin + rtot) || (apos[1]>=ymax - rtot)  || \
        (apos[2]<=zmin + rtot) || (apos[2]>=zmax - rtot)) {
        if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
            (thee->pmgp->bcfl != BCFL_MAP)) {
            Vnm_print(2, "Vpmg_dbForce:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n",
                      atomID, apos[0], apos[1], apos[2]);
            Vnm_print(2, "Vpmg_dbForce:    xmin = %g, xmax = %g\n",
                      xmin, xmax);
            Vnm_print(2, "Vpmg_dbForce:    ymin = %g, ymax = %g\n",
                      ymin, ymax);
            Vnm_print(2, "Vpmg_dbForce:    zmin = %g, zmax = %g\n",
                      zmin, zmax);
        }
        fflush(stderr);
    } else {

        /* Convert the atom position to grid reference frame */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;

        /* Integrate over points within this atom's (inflated) radius */
        rtot2 = VSQR(rtot);
        dx = rtot/hx;
        imin = (int)floor((position[0]-rtot)/hx);
        if (imin < 1) {
            Vnm_print(2, "Vpmg_dbForce:  Atom %d off grid!\n", atomID);
            return 0;
        }
        imax = (int)ceil((position[0]+rtot)/hx);
        if (imax > (nx-2)) {
            Vnm_print(2, "Vpmg_dbForce:  Atom %d off grid!\n", atomID);
            return 0;
        }
        jmin = (int)floor((position[1]-rtot)/hy);
        if (jmin < 1) {
            Vnm_print(2, "Vpmg_dbForce:  Atom %d off grid!\n", atomID);
            return 0;
        }
        jmax = (int)ceil((position[1]+rtot)/hy);
        if (jmax > (ny-2)) {
            Vnm_print(2, "Vpmg_dbForce:  Atom %d off grid!\n", atomID);
            return 0;
        }
        kmin = (int)floor((position[2]-rtot)/hzed);
        if (kmin < 1) {
            Vnm_print(2, "Vpmg_dbForce:  Atom %d off grid!\n", atomID);
            return 0;
        }
        kmax = (int)ceil((position[2]+rtot)/hzed);
        if (kmax > (nz-2)) {
            Vnm_print(2, "Vpmg_dbForce:  Atom %d off grid!\n", atomID);
            return 0;
        }
        for (i=imin; i<=imax; i++) {
            for (j=jmin; j<=jmax; j++) {
                for (k=kmin; k<=kmax; k++) {
                    /* i,j,k */
                    gpos[0] = (i+0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxijk = (thee->epsx[IJK(i,j,k)] - epsp)*depsi;

                    /* Select the correct function based on the surface definition
                        *	(now including the 7th order polynomial) */
                    Vpmg_splineSelect(srfm,acc, gpos, thee->splineWin, 0.,atom, dHxijk);
                    /*
                     switch (srfm) {
                         case VSM_SPLINE :
                             Vacc_splineAccGradAtomNorm(acc, gpos, thee->splineWin, 0.,
                                                        atom, dHxijk);
                             break;
                         case VSM_SPLINE4 :
                             Vacc_splineAccGradAtomNorm4(acc, gpos, thee->splineWin, 0.,
                                                         atom, dHxijk);
                             break;
                         default:
                             Vnm_print(2, "Vpmg_dbnbForce: Unknown surface method.\n");
                             return;
                     }
                     */
                    for (l=0; l<3; l++) dHxijk[l] *= Hxijk;
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j+0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijk = (thee->epsy[IJK(i,j,k)] - epsp)*depsi;

                    /* Select the correct function based on the surface definition
                        *	(now including the 7th order polynomial) */
                    Vpmg_splineSelect(srfm,acc, gpos, thee->splineWin, 0.,atom, dHyijk);

                    for (l=0; l<3; l++) dHyijk[l] *= Hyijk;
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k+0.5)*hzed + zmin;
                    Hzijk = (thee->epsz[IJK(i,j,k)] - epsp)*depsi;

                    /* Select the correct function based on the surface definition
                        *	(now including the 7th order polynomial) */
                    Vpmg_splineSelect(srfm,acc, gpos, thee->splineWin, 0.,atom, dHzijk);

                    for (l=0; l<3; l++) dHzijk[l] *= Hzijk;
                    /* i-1,j,k */
                    gpos[0] = (i-0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxim1jk = (thee->epsx[IJK(i-1,j,k)] - epsp)*depsi;

                    /* Select the correct function based on the surface definition
                        *	(now including the 7th order polynomial) */
                    Vpmg_splineSelect(srfm,acc, gpos, thee->splineWin, 0.,atom, dHxim1jk);

                    for (l=0; l<3; l++) dHxim1jk[l] *= Hxim1jk;
                    /* i,j-1,k */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j-0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijm1k = (thee->epsy[IJK(i,j-1,k)] - epsp)*depsi;

                    /* Select the correct function based on the surface definition
                        *	(now including the 7th order polynomial) */
                    Vpmg_splineSelect(srfm,acc, gpos, thee->splineWin, 0.,atom, dHyijm1k);

                    for (l=0; l<3; l++) dHyijm1k[l] *= Hyijm1k;
                    /* i,j,k-1 */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k-0.5)*hzed + zmin;
                    Hzijkm1 = (thee->epsz[IJK(i,j,k-1)] - epsp)*depsi;

                    /* Select the correct function based on the surface definition
                        *	(now including the 7th order polynomial) */
                    Vpmg_splineSelect(srfm,acc, gpos, thee->splineWin, 0.,atom, dHzijkm1);

                    for (l=0; l<3; l++) dHzijkm1[l] *= Hzijkm1;
                    /* *** CALCULATE DIELECTRIC BOUNDARY FORCES *** */
                    dbFmag = u[IJK(i,j,k)];
                    tgrad[0] =
                        (dHxijk[0]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                         +  dHxim1jk[0]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                        + (dHyijk[0]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                           +  dHyijm1k[0]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                        + (dHzijk[0]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                           + dHzijkm1[0]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    tgrad[1] =
                        (dHxijk[1]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                         +  dHxim1jk[1]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                        + (dHyijk[1]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                           +  dHyijm1k[1]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                        + (dHzijk[1]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                           + dHzijkm1[1]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    tgrad[2] =
                        (dHxijk[2]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                         +  dHxim1jk[2]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                        + (dHyijk[2]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                           +  dHyijm1k[2]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                        + (dHzijk[2]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                           + dHzijkm1[2]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    dbForce[0] += (dbFmag*tgrad[0]);
                    dbForce[1] += (dbFmag*tgrad[1]);
                    dbForce[2] += (dbFmag*tgrad[2]);

                } /* k loop */
            } /* j loop */
        } /* i loop */

        dbForce[0] = -dbForce[0]*hx*hy*hzed*deps*0.5*izmagic;
        dbForce[1] = -dbForce[1]*hx*hy*hzed*deps*0.5*izmagic;
        dbForce[2] = -dbForce[2]*hx*hy*hzed*deps*0.5*izmagic;
    }

    return 1;
}

VPUBLIC int Vpmg_qfForce(Vpmg *thee, double *force, int atomID,
  Vchrg_Meth chgm) {

    double tforce[3];

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* Check surface definition */
    if (chgm != VCM_BSPL2) {
        Vnm_print(2, "Vpmg_qfForce:  It is recommended that forces be \
calculated with the\n");
        Vnm_print(2, "Vpmg_qfForce:  cubic spline charge discretization \
scheme\n");
    }

    switch (chgm) {
        case VCM_TRIL:
            qfForceSpline1(thee, tforce, atomID);
            break;
        case VCM_BSPL2:
            qfForceSpline2(thee, tforce, atomID);
            break;
        case VCM_BSPL4:
            qfForceSpline4(thee, tforce, atomID);
            break;
        default:
            Vnm_print(2, "Vpmg_qfForce:  Undefined charge discretization \
method (%d)!\n", chgm);
            Vnm_print(2, "Vpmg_qfForce:  Forces not calculated!\n");
            return 0;
    }

    /* Assign forces */
    force[0] = tforce[0];
    force[1] = tforce[1];
    force[2] = tforce[2];

    return 1;
}


VPRIVATE void qfForceSpline1(Vpmg *thee, double *force, int atomID) {

    Vatom *atom;

    double *apos, position[3], hx, hy, hzed;
    double xmin, ymin, zmin, xmax, ymax, zmax;
    double dx, dy, dz;
    double *u, charge, ifloat, jfloat, kfloat;
    int nx, ny, nz, ihi, ilo, jhi, jlo, khi, klo;

    VASSERT(thee != VNULL);

    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    charge = Vatom_getCharge(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* If we aren't in the current position, then we're done */
    if (atom->partID == 0) return;

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = thee->u;

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax) || (apos[1]<=ymin) || \
        (apos[1]>=ymax) || (apos[2]<=zmin) || (apos[2]>=zmax)) {
        if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
            (thee->pmgp->bcfl != BCFL_MAP)) {
            Vnm_print(2, "Vpmg_qfForce:  Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n", atomID, apos[0], apos[1], apos[2]);
            Vnm_print(2, "Vpmg_qfForce:    xmin = %g, xmax = %g\n", xmin, xmax);
            Vnm_print(2, "Vpmg_qfForce:    ymin = %g, ymax = %g\n", ymin, ymax);
            Vnm_print(2, "Vpmg_qfForce:    zmin = %g, zmax = %g\n", zmin, zmax);
        }
        fflush(stderr);
    } else {

        /* Convert the atom position to grid coordinates */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;
        ifloat = position[0]/hx;
        jfloat = position[1]/hy;
        kfloat = position[2]/hzed;
        ihi = (int)ceil(ifloat);
        ilo = (int)floor(ifloat);
        jhi = (int)ceil(jfloat);
        jlo = (int)floor(jfloat);
        khi = (int)ceil(kfloat);
        klo = (int)floor(kfloat);
        VASSERT((ihi < nx) && (ihi >=0));
        VASSERT((ilo < nx) && (ilo >=0));
        VASSERT((jhi < ny) && (jhi >=0));
        VASSERT((jlo < ny) && (jlo >=0));
        VASSERT((khi < nz) && (khi >=0));
        VASSERT((klo < nz) && (klo >=0));
        dx = ifloat - (double)(ilo);
        dy = jfloat - (double)(jlo);
        dz = kfloat - (double)(klo);


#if 0
        Vnm_print(1, "Vpmg_qfForce: (DEBUG) u ~ %g\n",
          dx    *dy    *dz    *u[IJK(ihi,jhi,khi)]
         +dx    *dy    *(1-dz)*u[IJK(ihi,jhi,klo)]
         +dx    *(1-dy)*dz    *u[IJK(ihi,jlo,khi)]
         +dx    *(1-dy)*(1-dz)*u[IJK(ihi,jlo,klo)]
         +(1-dx)*dy    *dz    *u[IJK(ilo,jhi,khi)]
         +(1-dx)*dy    *(1-dz)*u[IJK(ilo,jhi,klo)]
         +(1-dx)*(1-dy)*dz    *u[IJK(ilo,jlo,khi)]
         +(1-dx)*(1-dy)*(1-dz)*u[IJK(ilo,jlo,klo)]);
#endif


        if ((dx > VPMGSMALL) && (VABS(1.0-dx) > VPMGSMALL)) {
            force[0] =
              -charge*(dy    *dz    *u[IJK(ihi,jhi,khi)]
                     + dy    *(1-dz)*u[IJK(ihi,jhi,klo)]
                     + (1-dy)*dz    *u[IJK(ihi,jlo,khi)]
                     + (1-dy)*(1-dz)*u[IJK(ihi,jlo,klo)]
                     - dy    *dz    *u[IJK(ilo,jhi,khi)]
                     - dy    *(1-dz)*u[IJK(ilo,jhi,klo)]
                     - (1-dy)*dz    *u[IJK(ilo,jlo,khi)]
                     - (1-dy)*(1-dz)*u[IJK(ilo,jlo,klo)])/hx;
        } else {
            force[0] = 0;
            Vnm_print(0,
              "Vpmg_qfForce:  Atom %d on x gridline; zero x-force\n", atomID);
        }
        if ((dy > VPMGSMALL) && (VABS(1.0-dy) > VPMGSMALL)) {
            force[1] =
              -charge*(dx    *dz    *u[IJK(ihi,jhi,khi)]
                     + dx    *(1-dz)*u[IJK(ihi,jhi,klo)]
                     - dx    *dz    *u[IJK(ihi,jlo,khi)]
                     - dx    *(1-dz)*u[IJK(ihi,jlo,klo)]
                     + (1-dx)*dz    *u[IJK(ilo,jhi,khi)]
                     + (1-dx)*(1-dz)*u[IJK(ilo,jhi,klo)]
                     - (1-dx)*dz    *u[IJK(ilo,jlo,khi)]
                     - (1-dx)*(1-dz)*u[IJK(ilo,jlo,klo)])/hy;
        } else {
            force[1] = 0;
            Vnm_print(0,
              "Vpmg_qfForce:  Atom %d on y gridline; zero y-force\n", atomID);
        }
        if ((dz > VPMGSMALL) && (VABS(1.0-dz) > VPMGSMALL)) {
            force[2] =
              -charge*(dy    *dx    *u[IJK(ihi,jhi,khi)]
                     - dy    *dx    *u[IJK(ihi,jhi,klo)]
                     + (1-dy)*dx    *u[IJK(ihi,jlo,khi)]
                     - (1-dy)*dx    *u[IJK(ihi,jlo,klo)]
                     + dy    *(1-dx)*u[IJK(ilo,jhi,khi)]
                     - dy    *(1-dx)*u[IJK(ilo,jhi,klo)]
                     + (1-dy)*(1-dx)*u[IJK(ilo,jlo,khi)]
                     - (1-dy)*(1-dx)*u[IJK(ilo,jlo,klo)])/hzed;
        } else {
            force[2] = 0;
            Vnm_print(0,
              "Vpmg_qfForce:  Atom %d on z gridline; zero z-force\n", atomID);
        }
    }
}

VPRIVATE void qfForceSpline2(Vpmg *thee, double *force, int atomID) {

    Vatom *atom;

    double *apos, position[3], hx, hy, hzed;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax;
    double mx, my, mz, dmx, dmy, dmz;
    double *u, charge, ifloat, jfloat, kfloat;
    int nx, ny, nz, im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1;
    int kp1, kp2, ii, jj, kk;

    VASSERT(thee != VNULL);

    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    charge = Vatom_getCharge(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* If we aren't in the current position, then we're done */
    if (atom->partID == 0) return;

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = thee->u;

    /* Make sure we're on the grid */
    if ((apos[0]<=(xmin+hx))   || (apos[0]>=(xmax-hx)) \
     || (apos[1]<=(ymin+hy))   || (apos[1]>=(ymax-hy)) \
     || (apos[2]<=(zmin+hzed)) || (apos[2]>=(zmax-hzed))) {
        if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
            (thee->pmgp->bcfl != BCFL_MAP)) {
            Vnm_print(2, "qfForceSpline2:  Atom #%d off the mesh \
                (ignoring)\n", atomID);
        }
        fflush(stderr);

    } else {

        /* Convert the atom position to grid coordinates */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;
        ifloat = position[0]/hx;
        jfloat = position[1]/hy;
        kfloat = position[2]/hzed;
        ip1 = (int)ceil(ifloat);
        ip2 = ip1 + 1;
        im1 = (int)floor(ifloat);
        im2 = im1 - 1;
        jp1 = (int)ceil(jfloat);
        jp2 = jp1 + 1;
        jm1 = (int)floor(jfloat);
        jm2 = jm1 - 1;
        kp1 = (int)ceil(kfloat);
        kp2 = kp1 + 1;
        km1 = (int)floor(kfloat);
        km2 = km1 - 1;

        /* This step shouldn't be necessary, but it saves nasty debugging
         * later on if something goes wrong */
        ip2 = VMIN2(ip2,nx-1);
        ip1 = VMIN2(ip1,nx-1);
        im1 = VMAX2(im1,0);
        im2 = VMAX2(im2,0);
        jp2 = VMIN2(jp2,ny-1);
        jp1 = VMIN2(jp1,ny-1);
        jm1 = VMAX2(jm1,0);
        jm2 = VMAX2(jm2,0);
        kp2 = VMIN2(kp2,nz-1);
        kp1 = VMIN2(kp1,nz-1);
        km1 = VMAX2(km1,0);
        km2 = VMAX2(km2,0);


        for (ii=im2; ii<=ip2; ii++) {
            mx = bspline2(VFCHI(ii,ifloat));
            dmx = dbspline2(VFCHI(ii,ifloat));
            for (jj=jm2; jj<=jp2; jj++) {
                my = bspline2(VFCHI(jj,jfloat));
                dmy = dbspline2(VFCHI(jj,jfloat));
                for (kk=km2; kk<=kp2; kk++) {
                    mz = bspline2(VFCHI(kk,kfloat));
                    dmz = dbspline2(VFCHI(kk,kfloat));

                    force[0] += (charge*dmx*my*mz*u[IJK(ii,jj,kk)])/hx;
                    force[1] += (charge*mx*dmy*mz*u[IJK(ii,jj,kk)])/hy;
                    force[2] += (charge*mx*my*dmz*u[IJK(ii,jj,kk)])/hzed;

                }
            }
        }

    }
}

VPRIVATE void qfForceSpline4(Vpmg *thee, double *force, int atomID) {

    Vatom *atom;
    double f, c, *u, *apos, position[3];

    /* Grid variables */
    int nx,ny,nz;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax;
    double hx, hy, hzed, ifloat, jfloat, kfloat;

    /* B-spline weights */
    double mx, my, mz, dmx, dmy, dmz;
    double mi, mj, mk;

    /* Loop indeces */
    int i, j, k, ii, jj, kk;
    int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1, kp1, kp2;

    /* field */
    double e[3];

    VASSERT(thee != VNULL);
    VASSERT(thee->filled);

    atom = Valist_getAtom(thee->pbe->alist, atomID);
    apos = Vatom_getPosition(atom);
    c = Vatom_getCharge(atom);

    for (i=0;i<3;i++){
        e[i] = 0.0;
    }

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = thee->u;

    /* Make sure we're on the grid */
    if ((apos[0]<=(xmin+2*hx))   || (apos[0]>=(xmax-2*hx)) \
        || (apos[1]<=(ymin+2*hy))   || (apos[1]>=(ymax-2*hy)) \
        || (apos[2]<=(zmin+2*hzed)) || (apos[2]>=(zmax-2*hzed))) {
        Vnm_print(2, "qfForceSpline4:  Atom off the mesh \
            (ignoring) %6.3f %6.3f %6.3f\n", apos[0], apos[1], apos[2]);
        fflush(stderr);
    } else {

        /* Convert the atom position to grid coordinates */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;
        ifloat = position[0]/hx;
        jfloat = position[1]/hy;
        kfloat = position[2]/hzed;
        ip1 = (int)ceil(ifloat);
        ip2 = ip1 + 2;
        im1 = (int)floor(ifloat);
        im2 = im1 - 2;
        jp1 = (int)ceil(jfloat);
        jp2 = jp1 + 2;
        jm1 = (int)floor(jfloat);
        jm2 = jm1 - 2;
        kp1 = (int)ceil(kfloat);
        kp2 = kp1 + 2;
        km1 = (int)floor(kfloat);
        km2 = km1 - 2;

        /* This step shouldn't be necessary, but it saves nasty debugging
            * later on if something goes wrong */
        ip2 = VMIN2(ip2,nx-1);
        ip1 = VMIN2(ip1,nx-1);
        im1 = VMAX2(im1,0);
        im2 = VMAX2(im2,0);
        jp2 = VMIN2(jp2,ny-1);
        jp1 = VMIN2(jp1,ny-1);
        jm1 = VMAX2(jm1,0);
        jm2 = VMAX2(jm2,0);
        kp2 = VMIN2(kp2,nz-1);
        kp1 = VMIN2(kp1,nz-1);
        km1 = VMAX2(km1,0);
        km2 = VMAX2(km2,0);

        for (ii=im2; ii<=ip2; ii++) {
            mi = VFCHI4(ii,ifloat);
            mx = bspline4(mi);
            dmx = dbspline4(mi);
            for (jj=jm2; jj<=jp2; jj++) {
                mj = VFCHI4(jj,jfloat);
                my = bspline4(mj);
                dmy = dbspline4(mj);
                for (kk=km2; kk<=kp2; kk++) {
                    mk = VFCHI4(kk,kfloat);
                    mz = bspline4(mk);
                    dmz = dbspline4(mk);
                    f = u[IJK(ii,jj,kk)];
                    /* Field */
                    e[0] += f*dmx*my*mz/hx;
                    e[1] += f*mx*dmy*mz/hy;
                    e[2] += f*mx*my*dmz/hzed;
                }
            }
        }
    }

    /* Monopole Force */
    force[0] = e[0]*c;
    force[1] = e[1]*c;
    force[2] = e[2]*c;

}

VPRIVATE void markFrac(
        double rtot, double *tpos,
        int nx, int ny, int nz,
        double hx, double hy, double hzed,
        double xmin, double ymin, double zmin,
        double *xarray, double *yarray, double *zarray) {

    int i, j, k, imin, imax, jmin, jmax, kmin, kmax;
    double dx, dx2, dy, dy2, dz, dz2, a000, a001, a010, a100, r2;
    double x, xp, xm, y, yp, ym, zp, z, zm, xspan, yspan, zspan;
    double rtot2, pos[3];

    /* Convert to grid reference frame */
    pos[0] = tpos[0] - xmin;
    pos[1] = tpos[1] - ymin;
    pos[2] = tpos[2] - zmin;

    rtot2 = VSQR(rtot);

    xspan = rtot + 2*hx;
    imin = VMAX2(0, (int)ceil((pos[0] - xspan)/hx));
    imax = VMIN2(nx-1, (int)floor((pos[0] + xspan)/hx));
    for (i=imin; i<=imax; i++) {
        x = hx*i;
        dx2 = VSQR(pos[0] - x);
        if (rtot2 > dx2) {
            yspan = VSQRT(rtot2 - dx2) + 2*hy;
        } else {
            yspan = 2*hy;
        }
        jmin = VMAX2(0,(int)ceil((pos[1] - yspan)/hy));
        jmax = VMIN2(ny-1,(int)floor((pos[1] + yspan)/hy));
        for (j=jmin; j<=jmax; j++) {
            y = hy*j;
            dy2 = VSQR(pos[1] - y);
            if (rtot2 > (dx2+dy2)) {
                zspan = VSQRT(rtot2-dx2-dy2) + 2*hzed;
            } else {
                zspan = 2*hzed;
            }
            kmin = VMAX2(0,(int)ceil((pos[2] - zspan)/hzed));
            kmax = VMIN2(nz-1,(int)floor((pos[2] + zspan)/hzed));
            for (k=kmin; k<=kmax; k++) {
                z = hzed*k;
                dz2 = VSQR(pos[2] - z);

                r2 = dx2 + dy2 + dz2;

                /* We need to determine the inclusion value a000 at (i,j,k) */
                if (r2 < rtot2) a000 = 1.0;
                else a000 = 0.0;

                /* We need to evaluate the values of x which intersect the
                 * sphere and determine if these are in the interval
                 * [(i,j,k), (i+1,j,k)] */
                if (r2 < (rtot2 - hx*hx)) a100 = 1.0;
                else if (r2 > (rtot2 + hx*hx)) a100 = 0.0;
                else if (rtot2 > (dy2 + dz2)) {
                    dx = VSQRT(rtot2 - dy2 - dz2);
                    xm = pos[0] - dx;
                    xp = pos[0] + dx;
                    if ((xm < x+hx) && (xm > x)) {
                        a100 = (xm - x)/hx;
                    } else if ((xp < x+hx) && (xp > x)) {
                        a100 = (xp - x)/hx;
                    }
                } else a100 = 0.0;

                /* We need to evaluate the values of y which intersect the
                 * sphere and determine if these are in the interval
                 * [(i,j,k), (i,j+1,k)] */
                if (r2 < (rtot2 - hy*hy)) a010 = 1.0;
                else if (r2 > (rtot2 + hy*hy)) a010 = 0.0;
                else if (rtot2 > (dx2 + dz2)) {
                    dy = VSQRT(rtot2 - dx2 - dz2);
                    ym = pos[1] - dy;
                    yp = pos[1] + dy;
                    if ((ym < y+hy) && (ym > y)) {
                        a010 = (ym - y)/hy;
                    } else if ((yp < y+hy) && (yp > y)) {
                        a010 = (yp - y)/hy;
                    }
                } else a010 = 0.0;

                /* We need to evaluate the values of y which intersect the
                 * sphere and determine if these are in the interval
                 * [(i,j,k), (i,j,k+1)] */
                if (r2 < (rtot2 - hzed*hzed)) a001 = 1.0;
                else if (r2 > (rtot2 + hzed*hzed)) a001 = 0.0;
                else if (rtot2 > (dx2 + dy2)) {
                    dz = VSQRT(rtot2 - dx2 - dy2);
                    zm = pos[2] - dz;
                    zp = pos[2] + dz;
                    if ((zm < z+hzed) && (zm > z)) {
                        a001 = (zm - z)/hzed;
                    } else if ((zp < z+hzed) && (zp > z)) {
                        a001 = (zp - z)/hzed;
                    }
                } else a001 = 0.0;

                if (a100 < xarray[IJK(i,j,k)]) xarray[IJK(i,j,k)] = a100;
                if (a010 < yarray[IJK(i,j,k)]) yarray[IJK(i,j,k)] = a010;
                if (a001 < zarray[IJK(i,j,k)]) zarray[IJK(i,j,k)] = a001;

            } /* k loop */
        } /* j loop */
    } /* i loop */
}

/*

 NOTE: This is the original version of the markSphere function. It's in here
 for reference and in case a reversion to the original code is needed.
 D. Gohara (2/14/08)
 */
/*
VPRIVATE void markSphere(
                         double rtot, double *tpos,
                         int nx, int ny, int nz,
                         double hx, double hy, double hzed,
                         double xmin, double ymin, double zmin,
                         double *array, double markVal) {

    int i, j, k, imin, imax, jmin, jmax, kmin, kmax;
    double dx, dx2, dy, dy2, dz, dz2;
    double rtot2, pos[3];

    // Convert to grid reference frame
    pos[0] = tpos[0] - xmin;
    pos[1] = tpos[1] - ymin;
    pos[2] = tpos[2] - zmin;

    rtot2 = VSQR(rtot);

    dx = rtot + 0.5*hx;
    imin = VMAX2(0,(int)ceil((pos[0] - dx)/hx));
    imax = VMIN2(nx-1,(int)floor((pos[0] + dx)/hx));
    for (i=imin; i<=imax; i++) {
        dx2 = VSQR(pos[0] - hx*i);
        if (rtot2 > dx2) {
            dy = VSQRT(rtot2 - dx2) + 0.5*hy;
        } else {
            dy = 0.5*hy;
        }
        jmin = VMAX2(0,(int)ceil((pos[1] - dy)/hy));
        jmax = VMIN2(ny-1,(int)floor((pos[1] + dy)/hy));
        for (j=jmin; j<=jmax; j++) {
            dy2 = VSQR(pos[1] - hy*j);
            if (rtot2 > (dx2+dy2)) {
                dz = VSQRT(rtot2-dx2-dy2)+0.5*hzed;
            } else {
                dz = 0.5*hzed;
            }
            kmin = VMAX2(0,(int)ceil((pos[2] - dz)/hzed));
            kmax = VMIN2(nz-1,(int)floor((pos[2] + dz)/hzed));
            for (k=kmin; k<=kmax; k++) {
                dz2 = VSQR(k*hzed - pos[2]);
                if ((dz2 + dy2 + dx2) <= rtot2) {
                    array[IJK(i,j,k)] = markVal;
                }
            } // k loop
        } // j loop
    } // i loop
}
*/
VPRIVATE void markSphere(double rtot, double *tpos,
                         int nx, int ny, int nz,
                         double hx, double hy, double hz,
                         double xmin, double ymin, double zmin,
                         double *array, double markVal) {

    int i, j, k;
    double fi,fj,fk;
    int imin, imax;
    int jmin, jmax;
    int kmin, kmax;
    double dx2, dy2, dz2;
    double xrange, yrange, zrange;
    double rtot2, posx, posy, posz;

    /* Convert to grid reference frame */
    posx = tpos[0] - xmin;
    posy = tpos[1] - ymin;
    posz = tpos[2] - zmin;

    rtot2 = VSQR(rtot);

    xrange = rtot + 0.5 * hx;
    yrange = rtot + 0.5 * hy;
    zrange = rtot + 0.5 * hz;

    imin = VMAX2(0, (int)ceil((posx - xrange)/hx));
    jmin = VMAX2(0, (int)ceil((posy - yrange)/hy));
    kmin = VMAX2(0, (int)ceil((posz - zrange)/hz));

    imax = VMIN2(nx-1, (int)floor((posx + xrange)/hx));
    jmax = VMIN2(ny-1, (int)floor((posy + yrange)/hy));
    kmax = VMIN2(nz-1, (int)floor((posz + zrange)/hz));

    for (i=imin,fi=imin; i<=imax; i++, fi+=1.) {
        dx2 = VSQR(posx - hx*fi);
        for (j=jmin,fj=jmin; j<=jmax; j++, fj+=1.) {
            dy2 = VSQR(posy - hy*fj);
            if((dx2 + dy2) > rtot2) continue;
            for (k=kmin, fk=kmin; k<=kmax; k++, fk+=1.) {
                dz2 = VSQR(posz - hz*fk);
                if ((dz2 + dy2 + dx2) <= rtot2) {
                    array[IJK(i,j,k)] = markVal;
                }
            }
        }
    }
}

VPRIVATE void zlapSolve(
        Vpmg *thee,
        double **solution,
        double **source,
        double **work1
        ) {

    /* NOTE:  this is an incredibly inefficient algorithm.  The next
     * improvement is to focus on only non-zero entries in the source term.
     * The best improvement is to use a fast sine transform */

    int n, nx, ny, nz, i, j, k, kx, ky, kz;
    double hx, hy, hzed, wx, wy, wz, xlen, ylen, zlen;
    double phix, phixp1, phixm1, phiy, phiym1, phiyp1, phiz, phizm1, phizp1;
    double norm, coef, proj, eigx, eigy, eigz;
    double ihx2, ihy2, ihzed2;
    double *u, *f, *phi;

    /* Snarf grid parameters */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    n = nx*ny*nz;
    hx = thee->pmgp->hx;
    ihx2 = 1.0/hx/hx;
    hy = thee->pmgp->hy;
    ihy2 = 1.0/hy/hy;
    hzed = thee->pmgp->hzed;
    ihzed2 = 1.0/hzed/hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Set solution and source array pointers */
    u = *solution;
    f = *source;
    phi = *work1;

    /* Zero out the solution vector */
    for (i=0; i<n; i++) thee->u[i] = 0.0;

    /* Iterate through the wavenumbers */
    for (kx=1; kx<(nx-1); kx++) {

        wx = (VPI*(double)kx)/((double)nx - 1.0);
        eigx = 2.0*ihx2*(1.0 - cos(wx));

        for (ky=1; ky<(ny-1); ky++) {

            wy = (VPI*(double)ky)/((double)ny - 1.0);
            eigy = 2.0*ihy2*(1.0 - cos(wy));

            for (kz=1; kz<(nz-1); kz++) {

                wz = (VPI*(double)kz)/((double)nz - 1.0);
                eigz = 2.0*ihzed2*(1.0 - cos(wz));

                /* Calculate the basis function.
                 * We could calculate each basis function as
                 *   phix(i) = sin(wx*i)
                 *   phiy(j) = sin(wy*j)
                 *   phiz(k) = sin(wz*k)
                 * However, this is likely to be very expensive.
                 * Therefore, we can use the fact that
                 *   phix(i+1) = (2-hx*hx*eigx)*phix(i) - phix(i-1)
                 * */
                for (i=1; i<(nx-1); i++) {
                    if (i == 1) {
                        phix = sin(wx*(double)i);
                        phixm1 = 0.0;
                    } else {
                        phixp1 = (2.0-hx*hx*eigx)*phix - phixm1;
                        phixm1 = phix;
                        phix = phixp1;
                    }
                    /* phix = sin(wx*(double)i); */
                    for (j=1; j<(ny-1); j++) {
                        if (j == 1) {
                            phiy = sin(wy*(double)j);
                            phiym1 = 0.0;
                        } else {
                            phiyp1 = (2.0-hy*hy*eigy)*phiy - phiym1;
                            phiym1 = phiy;
                            phiy = phiyp1;
                        }
                        /* phiy = sin(wy*(double)j); */
                        for (k=1; k<(nz-1); k++) {
                            if (k == 1) {
                                phiz = sin(wz*(double)k);
                                phizm1 = 0.0;
                            } else {
                                phizp1 = (2.0-hzed*hzed*eigz)*phiz - phizm1;
                                phizm1 = phiz;
                                phiz = phizp1;
                            }
                            /* phiz = sin(wz*(double)k); */

                            phi[IJK(i,j,k)] = phix*phiy*phiz;

                        }
                    }
                }

                /* Calculate the projection of the source function on this
                 * basis function */
                proj = 0.0;
                for (i=1; i<(nx-1); i++) {
                    for (j=1; j<(ny-1); j++) {
                        for (k=1; k<(nz-1); k++) {

                            proj += f[IJK(i,j,k)]*phi[IJK(i,j,k)];

                        } /* k loop */
                    } /* j loop */
                } /* i loop */

                /* Assemble the coefficient to weight the contribution of this
                 * basis function to the solution */
                /* The first contribution is the projection */
                coef = proj;
                /* The second contribution is the eigenvalue */
                coef = coef/(eigx + eigy + eigz);
                /* The third contribution is the normalization factor */
                coef = (8.0/xlen/ylen/zlen)*coef;
                /* The fourth contribution is from scaling the diagonal */
                /* coef = hx*hy*hzed*coef; */

                /* Evaluate the basis function at each grid point */
                for (i=1; i<(nx-1); i++) {
                    for (j=1; j<(ny-1); j++) {
                        for (k=1; k<(nz-1); k++) {

                            u[IJK(i,j,k)] += coef*phi[IJK(i,j,k)];

                        } /* k loop */
                    } /* j loop */
                } /* i loop */

            } /* kz loop */
        } /* ky loop */
    } /* kx loop */

}

VPUBLIC int Vpmg_solveLaplace(Vpmg *thee) {

    int i, j, k, ijk, nx, ny, nz, n, dilo, dihi, djlo, djhi, dklo, dkhi;
    double hx, hy, hzed, epsw, iepsw, scal, scalx, scaly, scalz;

    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    n = nx*ny*nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    epsw = Vpbe_getSolventDiel(thee->pbe);
    iepsw = 1.0/epsw;
    scal = hx*hy*hzed;
    scalx = hx*hy/hzed;
    scaly = hx*hzed/hy;
    scalz = hx*hy/hzed;

    if (!(thee->filled)) {
        Vnm_print(2, "Vpmg_solve:  Need to call Vpmg_fillco()!\n");
        return 0;
    }

    /* Load boundary conditions into the RHS array */
    for (i=1; i<(nx-1); i++) {

        if (i == 1) dilo = 1;
        else dilo = 0;
        if (i == nx-2) dihi = 1;
        else dihi = 0;

        for (j=1; j<(ny-1); j++) {

            if (j == 1) djlo = 1;
            else djlo = 0;
            if (j == ny-2) djhi = 1;
            else djhi = 0;

            for (k=1; k<(nz-1); k++) {

                if (k == 1) dklo = 1;
                else dklo = 0;
                if (k == nz-2) dkhi = 1;
                else dkhi = 0;

                /// @todo  Fix this to use Vmatrices
                thee->fcf[IJK(i,j,k)] = \
                      iepsw*scal*thee->charge[IJK(i,j,k)] \
                    + dilo*scalx*thee->gxcf[IJKx(j,k,0)] \
                    + dihi*scalx*thee->gxcf[IJKx(j,k,1)] \
                    + djlo*scaly*thee->gycf[IJKy(i,k,0)] \
                    + djhi*scaly*thee->gycf[IJKy(i,k,1)] \
                    + dklo*scalz*thee->gzcf[IJKz(i,j,0)] \
                    + dkhi*scalz*thee->gzcf[IJKz(i,j,1)] ;

            }
        }
    }

    /* Solve */
    zlapSolve( thee, &(thee->u), &(thee->fcf), &(thee->tcf) );

    /* Add boundary conditions to solution */
    /* i faces */
    for (j=0; j<ny; j++) {
        for (k=0; k<nz; k++) {
            thee->u[IJK(0,j,k)] = thee->gxcf[IJKx(j,k,0)];
            thee->u[IJK(nx-1,j,k)] = thee->gycf[IJKx(j,k,1)];
        }
    }
    /* j faces */
    for (i=0; i<nx; i++) {
        for (k=0; k<nz; k++) {
            thee->u[IJK(i,0,k)] = thee->gycf[IJKy(i,k,0)];
            thee->u[IJK(i,ny-1,k)] = thee->gycf[IJKy(i,k,1)];
        }
    }
    /* k faces */
    for (i=0; i<nx; i++) {
        for (j=0; j<ny; j++) {
            thee->u[IJK(i,j,0)] = thee->gzcf[IJKz(i,j,0)];
            thee->u[IJK(i,j,nz-1)] = thee->gzcf[IJKz(i,j,1)];
        }
    }

    return 1;

}

VPRIVATE double VFCHI4(int i, double f) {
  return (2.5+((double)(i)-(f)));
}

VPRIVATE double bspline4(double x) {

    double m, m2;
    static double one6 = 1.0/6.0;
    static double one8 = 1.0/8.0;
    static double one24 = 1.0/24.0;
    static double thirteen24 = 13.0/24.0;
    static double fourtyseven24 = 47.0/24.0;
    static double seventeen24 = 17.0/24.0;

    if ((x > 0.0) && (x <= 1.0)){
      m = x*x;
      return one24*m*m;
    } else if ((x > 1.0) && (x <= 2.0)){
      m = x - 1.0;
      m2 = m*m;
      return -one8 + one6*x + m2*(0.25 + one6*m - one6*m2);
    } else if ((x > 2.0) && (x <= 3.0)){
      m = x - 2.0;
      m2 = m*m;
      return -thirteen24 + 0.5*x + m2*(-0.25 - 0.5*m + 0.25*m2);
    } else if ((x > 3.0) && (x <= 4.0)){
      m = x - 3.0;
      m2 = m*m;
      return fourtyseven24 - 0.5*x + m2*(-0.25 + 0.5*m - one6*m2);
    } else if ((x > 4.0) && (x <= 5.0)){
      m = x - 4.0;
      m2 = m*m;
      return seventeen24 - one6*x + m2*(0.25 - one6*m + one24*m2);
    } else {
      return 0.0;
    }
}

VPUBLIC double dbspline4(double x) {

    double m, m2;
    static double one6 = 1.0/6.0;
    static double one3 = 1.0/3.0;
    static double two3 = 2.0/3.0;
    static double thirteen6 = 13.0/6.0;

    if ((x > 0.0) && (x <= 1.0)){
      m2 = x*x;
      return one6*x*m2;
    } else if ((x > 1.0) && (x <= 2.0)){
      m = x - 1.0;
      m2 = m*m;
      return -one3 + 0.5*x + m2*(0.5 - two3*m);
    } else if ((x > 2.0) && (x <= 3.0)){
      m = x - 2.0;
      m2 = m*m;
      return 1.5 - 0.5*x + m2*(-1.5 + m);
    } else if ((x > 3.0) && (x <= 4.0)){
      m = x - 3.0;
      m2 = m*m;
      return 1.0 - 0.5*x + m2*(1.5 - two3*m);
    } else if ((x > 4.0) && (x <= 5.0)){
      m = x - 4.0;
      m2 = m*m;
      return -thirteen6 + 0.5*x + m2*(-0.5 + one6*m);
    } else {
      return 0.0;
    }
}

VPUBLIC double d2bspline4(double x) {

    double m, m2;

    if ((x > 0.0) && (x <= 1.0)){
      return 0.5*x*x;
    } else if ((x > 1.0) && (x <= 2.0)){
      m = x - 1.0;
      m2 = m*m;
      return -0.5 + x - 2.0*m2;
    } else if ((x > 2.0) && (x <= 3.0)){
      m = x - 2.0;
      m2 = m*m;
      return 5.5 - 3.0*x + 3.0*m2;
    } else if ((x > 3.0) && (x <= 4.0)){
      m = x - 3.0;
      m2 = m*m;
      return -9.5 + 3.0*x - 2.0*m2;
    } else if ((x > 4.0) && (x <= 5.0)){
      m = x - 4.0;
      m2 = m*m;
      return 4.5 - x + 0.5*m2;
    } else {
      return 0.0;
    }
}

VPUBLIC double d3bspline4(double x) {

    if      ((x > 0.0) && (x <= 1.0)) return x;
    else if ((x > 1.0) && (x <= 2.0)) return 5.0 - 4.0 * x;
    else if ((x > 2.0) && (x <= 3.0)) return -15.0 + 6.0 * x;
    else if ((x > 3.0) && (x <= 4.0)) return 15.0 - 4.0 * x;
    else if ((x > 4.0) && (x <= 5.0)) return x - 5.0;
    else                              return 0.0;

}

VPUBLIC void fillcoPermanentMultipole(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    /* Coversions */
    double zmagic, f;
    /* Grid */
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, position[3], ifloat, jfloat, kfloat;
    double hx, hy, hzed, *apos;
    /* Multipole */
    double charge, *dipole,*quad;
    double c,ux,uy,uz,qxx,qyx,qyy,qzx,qzy,qzz,qave;
    /* B-spline weights */
    double mx,my,mz,dmx,dmy,dmz,d2mx,d2my,d2mz;
    double mi,mj,mk;
    /* Loop variables */
    int i, ii, jj, kk, nx, ny, nz, iatom;
    int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1, kp1, kp2;

    /* sanity check */
    double mir,mjr,mkr,mr2;
    double debye,mc,mux,muy,muz,mqxx,mqyx,mqyy,mqzx,mqzy,mqzz;

    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    zmagic = Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Conversion */
    f = zmagic/(hx*hy*hzed);

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Fill in the source term (permanent atomic multipoles) */
    Vnm_print(0, "fillcoPermanentMultipole:  filling in source term.\n");
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);

        c = Vatom_getCharge(atom)*f;

#if defined(WITH_TINKER)
        dipole = Vatom_getDipole(atom);
        ux = dipole[0]/hx*f;
        uy = dipole[1]/hy*f;
        uz = dipole[2]/hzed*f;
        quad = Vatom_getQuadrupole(atom);
        qxx = (1.0/3.0)*quad[0]/(hx*hx)*f;
        qyx = (2.0/3.0)*quad[3]/(hx*hy)*f;
        qyy = (1.0/3.0)*quad[4]/(hy*hy)*f;
        qzx = (2.0/3.0)*quad[6]/(hzed*hx)*f;
        qzy = (2.0/3.0)*quad[7]/(hzed*hy)*f;
        qzz = (1.0/3.0)*quad[8]/(hzed*hzed)*f;
#else
        ux = 0.0;
        uy = 0.0;
        uz = 0.0;
        qxx = 0.0;
        qyx = 0.0;
        qyy = 0.0;
        qzx = 0.0;
        qzy = 0.0;
        qzz = 0.0;
#endif /* if defined(WITH_TINKER) */

        /* check
        mc = 0.0;
        mux = 0.0;
        muy = 0.0;
        muz = 0.0;
        mqxx = 0.0;
        mqyx = 0.0;
        mqyy = 0.0;
        mqzx = 0.0;
        mqzy = 0.0;
        mqzz = 0.0; */

        /* Make sure we're on the grid */
        if ((apos[0]<=(xmin-2*hx)) || (apos[0]>=(xmax+2*hx))  || \
            (apos[1]<=(ymin-2*hy)) || (apos[1]>=(ymax+2*hy))  || \
            (apos[2]<=(zmin-2*hzed)) || (apos[2]>=(zmax+2*hzed))) {
            Vnm_print(2, "fillcoPermanentMultipole: Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring this atom):\n", iatom, apos[0], apos[1], apos[2]);
            Vnm_print(2, "fillcoPermanentMultipole: xmin = %g, xmax = %g\n", xmin, xmax);
            Vnm_print(2, "fillcoPermanentMultipole: ymin = %g, ymax = %g\n", ymin, ymax);
            Vnm_print(2, "fillcoPermanentMultipole: zmin = %g, zmax = %g\n", zmin, zmax);
            fflush(stderr);
        } else {

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* Figure out which vertices we're next to */
            ifloat = position[0]/hx;
            jfloat = position[1]/hy;
            kfloat = position[2]/hzed;

            ip1   = (int)ceil(ifloat);
            ip2   = ip1 + 2;
            im1   = (int)floor(ifloat);
            im2   = im1 - 2;
            jp1   = (int)ceil(jfloat);
            jp2   = jp1 + 2;
            jm1   = (int)floor(jfloat);
            jm2   = jm1 - 2;
            kp1   = (int)ceil(kfloat);
            kp2   = kp1 + 2;
            km1   = (int)floor(kfloat);
            km2   = km1 - 2;

            /* This step shouldn't be necessary, but it saves nasty debugging
             * later on if something goes wrong */
            ip2 = VMIN2(ip2,nx-1);
            ip1 = VMIN2(ip1,nx-1);
            im1 = VMAX2(im1,0);
            im2 = VMAX2(im2,0);
            jp2 = VMIN2(jp2,ny-1);
            jp1 = VMIN2(jp1,ny-1);
            jm1 = VMAX2(jm1,0);
            jm2 = VMAX2(jm2,0);
            kp2 = VMIN2(kp2,nz-1);
            kp1 = VMIN2(kp1,nz-1);
            km1 = VMAX2(km1,0);
            km2 = VMAX2(km2,0);

            /* Now assign fractions of the charge to the nearby verts */
            for (ii=im2; ii<=ip2; ii++) {
                mi = VFCHI4(ii,ifloat);
                mx = bspline4(mi);
                dmx = dbspline4(mi);
                d2mx = d2bspline4(mi);
                for (jj=jm2; jj<=jp2; jj++) {
                    mj = VFCHI4(jj,jfloat);
                    my = bspline4(mj);
                    dmy = dbspline4(mj);
                    d2my = d2bspline4(mj);
                    for (kk=km2; kk<=kp2; kk++) {
                        mk = VFCHI4(kk,kfloat);
                        mz = bspline4(mk);
                        dmz = dbspline4(mk);
                        d2mz = d2bspline4(mk);
                        charge = mx*my*mz*c -
                         dmx*my*mz*ux - mx*dmy*mz*uy - mx*my*dmz*uz +
                         d2mx*my*mz*qxx +
                         dmx*dmy*mz*qyx + mx*d2my*mz*qyy +
                         dmx*my*dmz*qzx + mx*dmy*dmz*qzy + mx*my*d2mz*qzz;
                        thee->charge[IJK(ii,jj,kk)] += charge;

                        /* sanity check - recalculate traceless multipoles
                           from the grid charge distribution for this
                           site.

                        mir = (mi - 2.5) * hx;
                        mjr = (mj - 2.5) * hy;
                        mkr = (mk - 2.5) * hzed;
                        mr2 = mir*mir+mjr*mjr+mkr*mkr;
                        mc += charge;
                        mux += mir * charge;
                        muy += mjr * charge;
                        muz += mkr * charge;
                        mqxx += (1.5*mir*mir - 0.5*mr2) * charge;
                        mqyx += 1.5*mjr*mir * charge;
                        mqyy += (1.5*mjr*mjr - 0.5*mr2) * charge;
                        mqzx += 1.5*mkr*mir * charge;
                        mqzy += 1.5*mkr*mjr * charge;
                        mqzz += (1.5*mkr*mkr - 0.5*mr2) * charge;
                         */
                    }
                }
            }
        } /* endif (on the mesh) */

        /* print out the Grid vs. Ideal Point Multipole. */

        /*
        debye = 4.8033324;
        mc = mc/f;
        mux = mux/f*debye;
        muy = muy/f*debye;
        muz = muz/f*debye;
        mqxx = mqxx/f*debye;
        mqyy = mqyy/f*debye;
        mqzz = mqzz/f*debye;
        mqyx = mqyx/f*debye;
        mqzx = mqzx/f*debye;
        mqzy = mqzy/f*debye;

        printf(" Grid v. Actual Permanent Multipole for Site %i\n",iatom);
        printf(" G: %10.6f\n",mc);
        printf(" A: %10.6f\n\n",c/f);
        printf(" G: %10.6f %10.6f %10.6f\n",mux,muy,muz);
        printf(" A: %10.6f %10.6f %10.6f\n\n",
                 (ux * hx / f) * debye,
                 (uy * hy / f) * debye,
                 (uz * hzed /f) * debye);
        printf(" G: %10.6f\n",mqxx);
        printf(" A: %10.6f\n",quad[0]*debye);
        printf(" G: %10.6f %10.6f\n",mqyx,mqyy);
        printf(" A: %10.6f %10.6f\n",quad[3]*debye,quad[4]*debye);
        printf(" G: %10.6f %10.6f %10.6f\n",mqzx,mqzy,mqzz);
        printf(" A: %10.6f %10.6f %10.6f\n\n",
                quad[6]*debye,quad[7]*debye,quad[8]*debye);  */

    } /* endfor (each atom) */
}

#if defined(WITH_TINKER)

VPUBLIC void fillcoInducedDipole(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    /* Conversions */
    double zmagic, f;
    /* Grid */
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, ifloat, jfloat, kfloat;
    double hx, hy, hzed, *apos, position[3];
    /* B-spline weights */
    double mx, my, mz, dmx, dmy, dmz;
    /* Dipole */
    double charge, *dipole, ux,uy,uz;
    double mi,mj,mk;
    /* Loop indeces */
    int i, ii, jj, kk, nx, ny, nz, iatom;
    int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1, kp1, kp2;

    double debye;
    double mux,muy,muz;
    double mir,mjr,mkr;

    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    zmagic = Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Conversion */
    f = zmagic/(hx*hy*hzed);

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Fill in the source term (induced dipoles) */
    Vnm_print(0, "fillcoInducedDipole:  filling in the source term.\n");
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);

        dipole = Vatom_getInducedDipole(atom);
        ux = dipole[0]/hx*f;
        uy = dipole[1]/hy*f;
        uz = dipole[2]/hzed*f;

        mux = 0.0;
        muy = 0.0;
        muz = 0.0;

        /* Make sure we're on the grid */
        if ((apos[0]<=(xmin-2*hx)) || (apos[0]>=(xmax+2*hx))  || \
            (apos[1]<=(ymin-2*hy)) || (apos[1]>=(ymax+2*hy))  || \
            (apos[2]<=(zmin-2*hzed)) || (apos[2]>=(zmax+2*hzed))) {
            Vnm_print(2, "fillcoInducedDipole: Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring this atom):\n", iatom, apos[0], apos[1], apos[2]);
            Vnm_print(2, "fillcoInducedDipole: xmin = %g, xmax = %g\n", xmin, xmax);
            Vnm_print(2, "fillcoInducedDipole: ymin = %g, ymax = %g\n", ymin, ymax);
            Vnm_print(2, "fillcoInducedDipole: zmin = %g, zmax = %g\n", zmin, zmax);
            fflush(stderr);
        } else {

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* Figure out which vertices we're next to */
            ifloat = position[0]/hx;
            jfloat = position[1]/hy;
            kfloat = position[2]/hzed;

            ip1   = (int)ceil(ifloat);
            ip2   = ip1 + 2;
            im1   = (int)floor(ifloat);
            im2   = im1 - 2;
            jp1   = (int)ceil(jfloat);
            jp2   = jp1 + 2;
            jm1   = (int)floor(jfloat);
            jm2   = jm1 - 2;
            kp1   = (int)ceil(kfloat);
            kp2   = kp1 + 2;
            km1   = (int)floor(kfloat);
            km2   = km1 - 2;

            /* This step shouldn't be necessary, but it saves nasty debugging
             * later on if something goes wrong */
            ip2 = VMIN2(ip2,nx-1);
            ip1 = VMIN2(ip1,nx-1);
            im1 = VMAX2(im1,0);
            im2 = VMAX2(im2,0);
            jp2 = VMIN2(jp2,ny-1);
            jp1 = VMIN2(jp1,ny-1);
            jm1 = VMAX2(jm1,0);
            jm2 = VMAX2(jm2,0);
            kp2 = VMIN2(kp2,nz-1);
            kp1 = VMIN2(kp1,nz-1);
            km1 = VMAX2(km1,0);
            km2 = VMAX2(km2,0);

            /* Now assign fractions of the dipole to the nearby verts */
            for (ii=im2; ii<=ip2; ii++) {
                mi = VFCHI4(ii,ifloat);
                mx = bspline4(mi);
                dmx = dbspline4(mi);
                for (jj=jm2; jj<=jp2; jj++) {
                    mj = VFCHI4(jj,jfloat);
                    my = bspline4(mj);
                    dmy = dbspline4(mj);
                    for (kk=km2; kk<=kp2; kk++) {
                        mk = VFCHI4(kk,kfloat);
                        mz = bspline4(mk);
                        dmz = dbspline4(mk);
                        charge = -dmx*my*mz*ux - mx*dmy*mz*uy - mx*my*dmz*uz;
                        thee->charge[IJK(ii,jj,kk)] += charge;

                        /*
                        mir = (mi - 2.5) * hx;
                        mjr = (mj - 2.5) * hy;
                        mkr = (mk - 2.5) * hzed;
                        mux += mir * charge;
                        muy += mjr * charge;
                        muz += mkr * charge;
                        */
                    }
                }
            }
        } /* endif (on the mesh) */

        /* check
        debye = 4.8033324;
        mux = mux/f*debye;
        muy = muy/f*debye;
        muz = muz/f*debye;

        printf(" Grid v. Actual Induced Dipole for Site %i\n",iatom);
        printf(" G: %10.6f %10.6f %10.6f\n",mux,muy,muz);
        printf(" A: %10.6f %10.6f %10.6f\n\n",
                 (ux * hx / f) * debye,
                 (uy * hy / f) * debye,
                 (uz * hzed /f) * debye);
         */

    } /* endfor (each atom) */
}

VPUBLIC void fillcoNLInducedDipole(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    /* Conversions */
    double zmagic, f;
    /* Grid */
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, ifloat, jfloat, kfloat;
    double hx, hy, hzed, *apos, position[3];
    /* B-spline weights */
    double mx, my, mz, dmx, dmy, dmz;
    /* Dipole */
    double charge, *dipole, ux,uy,uz;
    double mi,mj,mk;
    /* Loop indeces */
    int i, ii, jj, kk, nx, ny, nz, iatom;
    int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1, kp1, kp2;

    /* sanity check
    double debye;
    double mux,muy,muz;
    double mir,mjr,mkr;
     */

    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    zmagic = Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Conversion */
    f = zmagic/(hx*hy*hzed);

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Fill in the source term (non-local induced dipoles) */
    Vnm_print(0, "fillcoNLInducedDipole:  filling in source term.\n");
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);

        dipole = Vatom_getNLInducedDipole(atom);
        ux = dipole[0]/hx*f;
        uy = dipole[1]/hy*f;
        uz = dipole[2]/hzed*f;

        /*
        mux = 0.0;
        muy = 0.0;
        muz = 0.0;
         */

        /* Make sure we're on the grid */
        if ((apos[0]<=(xmin-2*hx)) || (apos[0]>=(xmax+2*hx))  || \
            (apos[1]<=(ymin-2*hy)) || (apos[1]>=(ymax+2*hy))  || \
            (apos[2]<=(zmin-2*hzed)) || (apos[2]>=(zmax+2*hzed))) {
            Vnm_print(2, "fillcoNLInducedDipole: Atom #%d at (%4.3f, %4.3f,%4.3f) is off the mesh (ignoring this atom):\n", iatom, apos[0], apos[1], apos[2]);
            Vnm_print(2, "fillcoNLInducedDipole: xmin = %g, xmax = %g\n", xmin, xmax);
            Vnm_print(2, "fillcoNLInducedDipole: ymin = %g, ymax = %g\n", ymin, ymax);
            Vnm_print(2, "fillcoNLInducedDipole: zmin = %g, zmax = %g\n", zmin, zmax);
            fflush(stderr);
        } else {

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* Figure out which vertices we're next to */
            ifloat = position[0]/hx;
            jfloat = position[1]/hy;
            kfloat = position[2]/hzed;

            ip1   = (int)ceil(ifloat);
            ip2   = ip1 + 2;
            im1   = (int)floor(ifloat);
            im2   = im1 - 2;
            jp1   = (int)ceil(jfloat);
            jp2   = jp1 + 2;
            jm1   = (int)floor(jfloat);
            jm2   = jm1 - 2;
            kp1   = (int)ceil(kfloat);
            kp2   = kp1 + 2;
            km1   = (int)floor(kfloat);
            km2   = km1 - 2;

            /* This step shouldn't be necessary, but it saves nasty debugging
             * later on if something goes wrong */
            ip2 = VMIN2(ip2,nx-1);
            ip1 = VMIN2(ip1,nx-1);
            im1 = VMAX2(im1,0);
            im2 = VMAX2(im2,0);
            jp2 = VMIN2(jp2,ny-1);
            jp1 = VMIN2(jp1,ny-1);
            jm1 = VMAX2(jm1,0);
            jm2 = VMAX2(jm2,0);
            kp2 = VMIN2(kp2,nz-1);
            kp1 = VMIN2(kp1,nz-1);
            km1 = VMAX2(km1,0);
            km2 = VMAX2(km2,0);

            /* Now assign fractions of the non local induced dipole
               to the nearby verts */
            for (ii=im2; ii<=ip2; ii++) {
                mi = VFCHI4(ii,ifloat);
                mx = bspline4(mi);
                dmx = dbspline4(mi);
                for (jj=jm2; jj<=jp2; jj++) {
                    mj = VFCHI4(jj,jfloat);
                    my = bspline4(mj);
                    dmy = dbspline4(mj);
                    for (kk=km2; kk<=kp2; kk++) {
                        mk = VFCHI4(kk,kfloat);
                        mz = bspline4(mk);
                        dmz = dbspline4(mk);
                        charge = -dmx*my*mz*ux - mx*dmy*mz*uy - mx*my*dmz*uz;
                        thee->charge[IJK(ii,jj,kk)] += charge;

                        /*
                        mir = (mi - 2.5) * hx;
                        mjr = (mj - 2.5) * hy;
                        mkr = (mk - 2.5) * hzed;
                        mux += mir * charge;
                        muy += mjr * charge;
                        muz += mkr * charge;
                         */
                    }
                }
            }
        } /* endif (on the mesh) */

        /*
        debye = 4.8033324;
        mux = mux/f*debye;
        muy = muy/f*debye;
        muz = muz/f*debye;

        printf(" Grid v. Actual Non-Local Induced Dipole for Site %i\n",iatom);
        printf(" G: %10.6f %10.6f %10.6f\n",mux,muy,muz);
        printf(" A: %10.6f %10.6f %10.6f\n\n",
                 (ux * hx / f) * debye,
                 (uy * hy / f) * debye,
                 (uz * hzed /f) * debye); */

    } /* endfor (each atom) */
}

VPUBLIC double Vpmg_qfPermanentMultipoleEnergy(Vpmg *thee, int atomID) {

    double *u;
    Vatom *atom;
    /* Grid variables */
    int nx, ny, nz;
    double xmax, xmin, ymax, ymin, zmax, zmin;
    double hx, hy, hzed, ifloat, jfloat, kfloat;
    double mi, mj, mk;
    double *position;
    /* B-spline weights */
    double mx, my, mz, dmx, dmy, dmz, d2mx, d2my, d2mz;
    /* Loop indeces */
    int ip1,ip2,im1,im2,jp1,jp2,jm1,jm2,kp1,kp2,km1,km2;
    int i,j,ii,jj,kk;
    /* Potential, field, field gradient and multipole components */
    double pot, rfe[3], rfde[3][3], energy;
    double f, charge, *dipole, *quad;
    double qxx, qyx, qyy, qzx, qzy, qzz;


    VASSERT(thee != VNULL);
    VASSERT(thee->filled);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmax = thee->xf[nx-1];
    ymax = thee->yf[ny-1];
    zmax = thee->zf[nz-1];
    xmin = thee->xf[0];
    ymin = thee->yf[0];
    zmin = thee->zf[0];

    u = thee->u;

    atom = Valist_getAtom(thee->pbe->alist, atomID);

    /* Currently all atoms must be in the same partition. */

    VASSERT(atom->partID != 0);

    /* Convert the atom position to grid coordinates */

    position = Vatom_getPosition(atom);
    ifloat = (position[0] - xmin)/hx;
    jfloat = (position[1] - ymin)/hy;
    kfloat = (position[2] - zmin)/hzed;
    ip1 = (int)ceil(ifloat);
    ip2 = ip1 + 2;
    im1 = (int)floor(ifloat);
    im2 = im1 - 2;
    jp1 = (int)ceil(jfloat);
    jp2 = jp1 + 2;
    jm1 = (int)floor(jfloat);
    jm2 = jm1 - 2;
    kp1 = (int)ceil(kfloat);
    kp2 = kp1 + 2;
    km1 = (int)floor(kfloat);
    km2 = km1 - 2;

    /* This step shouldn't be necessary, but it saves nasty debugging
     * later on if something goes wrong */
    ip2 = VMIN2(ip2,nx-1);
    ip1 = VMIN2(ip1,nx-1);
    im1 = VMAX2(im1,0);
    im2 = VMAX2(im2,0);
    jp2 = VMIN2(jp2,ny-1);
    jp1 = VMIN2(jp1,ny-1);
    jm1 = VMAX2(jm1,0);
    jm2 = VMAX2(jm2,0);
    kp2 = VMIN2(kp2,nz-1);
    kp1 = VMIN2(kp1,nz-1);
    km1 = VMAX2(km1,0);
    km2 = VMAX2(km2,0);

    /* Initialize observables to zero */
    energy = 0.0;
    pot = 0.0;
    for (i=0;i<3;i++){
       rfe[i] = 0.0;
       for (j=0;j<3;j++){
         rfde[i][j] = 0.0;
       }
    }

    for (ii=im2; ii<=ip2; ii++) {
      mi = VFCHI4(ii,ifloat);
      mx = bspline4(mi);
      dmx = dbspline4(mi);
      d2mx = d2bspline4(mi);
      for (jj=jm2; jj<=jp2; jj++) {
        mj = VFCHI4(jj,jfloat);
        my = bspline4(mj);
        dmy = dbspline4(mj);
        d2my = d2bspline4(mj);
        for (kk=km2; kk<=kp2; kk++) {
          mk = VFCHI4(kk,kfloat);
          mz = bspline4(mk);
          dmz = dbspline4(mk);
          d2mz = d2bspline4(mk);
          f = u[IJK(ii,jj,kk)];
          /* potential */
          pot  += f*mx*my*mz;
          /* field */
          rfe[0] += f*dmx*my*mz/hx;
          rfe[1] += f*mx*dmy*mz/hy;
          rfe[2] += f*mx*my*dmz/hzed;
          /* field gradient */
          rfde[0][0] += f*d2mx*my*mz/(hx*hx);
          rfde[1][0] += f*dmx*dmy*mz/(hy*hx);
          rfde[1][1] += f*mx*d2my*mz/(hy*hy);
          rfde[2][0] += f*dmx*my*dmz/(hx*hzed);
          rfde[2][1] += f*mx*dmy*dmz/(hy*hzed);
          rfde[2][2] += f*mx*my*d2mz/(hzed*hzed);
        }
      }
    }

    charge = Vatom_getCharge(atom);
    dipole = Vatom_getDipole(atom);
    quad = Vatom_getQuadrupole(atom);
    qxx = quad[0]/3.0;
    qyx = quad[3]/3.0;
    qyy = quad[4]/3.0;
    qzx = quad[6]/3.0;
    qzy = quad[7]/3.0;
    qzz = quad[8]/3.0;

    energy =   pot * charge
             - rfe[0] * dipole[0]
             - rfe[1] * dipole[1]
             - rfe[2] * dipole[2]
             +     rfde[0][0]*qxx
             + 2.0*rfde[1][0]*qyx +     rfde[1][1]*qyy
             + 2.0*rfde[2][0]*qzx + 2.0*rfde[2][1]*qzy + rfde[2][2]*qzz;

    return energy;
}

VPUBLIC void Vpmg_fieldSpline4(Vpmg *thee, int atomID, double field[3]) {

    Vatom *atom;
    double *u, f;
    /* Grid variables */
    int nx, ny, nz;
    double xmax, xmin, ymax, ymin, zmax, zmin;
    double hx, hy, hzed, ifloat, jfloat, kfloat;
    double *apos, position[3];
    /* B-Spline weights */
    double mx, my, mz, dmx, dmy, dmz;
    double mi, mj, mk;
    /* Loop indeces */
    int ip1,ip2,im1,im2,jp1,jp2,jm1,jm2,kp1,kp2,km1,km2;
    int i,j,ii,jj,kk;


    VASSERT (thee != VNULL);

    /* Get the mesh information */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xmax = thee->xf[nx-1];
    ymax = thee->yf[ny-1];
    zmax = thee->zf[nz-1];
    xmin = thee->xf[0];
    ymin = thee->yf[0];
    zmin = thee->zf[0];

    u = thee->u;

    atom = Valist_getAtom(thee->pbe->alist, atomID);

    /* Currently all atoms must be in the same partition. */

    VASSERT (atom->partID != 0);

    /* Convert the atom position to grid coordinates */

    apos = Vatom_getPosition(atom);
    position[0] = apos[0] - xmin;
    position[1] = apos[1] - ymin;
    position[2] = apos[2] - zmin;
    ifloat = position[0]/hx;
    jfloat = position[1]/hy;
    kfloat = position[2]/hzed;
    ip1 = (int)ceil(ifloat);
    ip2 = ip1 + 2;
    im1 = (int)floor(ifloat);
    im2 = im1 - 2;
    jp1 = (int)ceil(jfloat);
    jp2 = jp1 + 2;
    jm1 = (int)floor(jfloat);
    jm2 = jm1 - 2;
    kp1 = (int)ceil(kfloat);
    kp2 = kp1 + 2;
    km1 = (int)floor(kfloat);
    km2 = km1 - 2;

    /* This step shouldn't be necessary, but it saves nasty debugging
     * later on if something goes wrong */
    ip2 = VMIN2(ip2,nx-1);
    ip1 = VMIN2(ip1,nx-1);
    im1 = VMAX2(im1,0);
    im2 = VMAX2(im2,0);
    jp2 = VMIN2(jp2,ny-1);
    jp1 = VMIN2(jp1,ny-1);
    jm1 = VMAX2(jm1,0);
    jm2 = VMAX2(jm2,0);
    kp2 = VMIN2(kp2,nz-1);
    kp1 = VMIN2(kp1,nz-1);
    km1 = VMAX2(km1,0);
    km2 = VMAX2(km2,0);

    for (i=0;i<3;i++){
       field[i] = 0.0;
    }

    for (ii=im2; ii<=ip2; ii++) {
      mi = VFCHI4(ii,ifloat);
      mx = bspline4(mi);
      dmx = dbspline4(mi);
      for (jj=jm2; jj<=jp2; jj++) {
        mj = VFCHI4(jj,jfloat);
        my = bspline4(mj);
        dmy = dbspline4(mj);
        for (kk=km2; kk<=kp2; kk++) {
          mk = VFCHI4(kk,kfloat);
          mz = bspline4(mk);
          dmz = dbspline4(mk);
          f = u[IJK(ii,jj,kk)];

          field[0] += f*dmx*my*mz/hx;
          field[1] += f*mx*dmy*mz/hy;
          field[2] += f*mx*my*dmz/hzed;
        }
      }
    }
}

VPUBLIC void Vpmg_qfPermanentMultipoleForce(Vpmg *thee, int atomID,
                                           double force[3], double torque[3]) {

    Vatom *atom;
    double f, *u, *apos, position[3];

    /* Grid variables */
    int nx,ny,nz;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax;
    double hx, hy, hzed, ifloat, jfloat, kfloat;

    /* B-spline weights */
    double mx, my, mz, dmx, dmy, dmz, d2mx, d2my, d2mz, d3mx, d3my, d3mz;
    double mi, mj, mk;

    /* Loop indeces */
    int i, j, k, ii, jj, kk;
    int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1, kp1, kp2;

    /* Potential, field, field gradient and 2nd field gradient */
    double pot, e[3], de[3][3], d2e[3][3][3];

    /* Permanent multipole components */
    double *dipole, *quad;
    double c, ux, uy, uz, qxx, qxy, qxz, qyx, qyy, qyz, qzx, qzy, qzz;

    VASSERT(thee != VNULL);
    VASSERT(thee->filled);

    atom = Valist_getAtom(thee->pbe->alist, atomID);

    /* Currently all atoms must be in the same partition. */

    VASSERT(atom->partID != 0);

    apos = Vatom_getPosition(atom);

    c = Vatom_getCharge(atom);
    dipole = Vatom_getDipole(atom);
    ux = dipole[0];
    uy = dipole[1];
    uz = dipole[2];
    quad = Vatom_getQuadrupole(atom);
    qxx = quad[0]/3.0;
    qxy = quad[1]/3.0;
    qxz = quad[2]/3.0;
    qyx = quad[3]/3.0;
    qyy = quad[4]/3.0;
    qyz = quad[5]/3.0;
    qzx = quad[6]/3.0;
    qzy = quad[7]/3.0;
    qzz = quad[8]/3.0;

    /* Initialize observables */
    pot = 0.0;
    for (i=0;i<3;i++){
       e[i] = 0.0;
       for (j=0;j<3;j++){
          de[i][j] = 0.0;
          for (k=0;k<3;k++){
             d2e[i][j][k] = 0.0;
          }
       }
    }

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = thee->u;

    /* Make sure we're on the grid */
    if ((apos[0]<=(xmin+2*hx))   || (apos[0]>=(xmax-2*hx)) \
     || (apos[1]<=(ymin+2*hy))   || (apos[1]>=(ymax-2*hy)) \
     || (apos[2]<=(zmin+2*hzed)) || (apos[2]>=(zmax-2*hzed))) {
        Vnm_print(2, "qfPermanentMultipoleForce:  Atom off the mesh (ignoring) %6.3f %6.3f %6.3f\n", apos[0], apos[1], apos[2]);
        fflush(stderr);
    } else {

        /* Convert the atom position to grid coordinates */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;
        ifloat = position[0]/hx;
        jfloat = position[1]/hy;
        kfloat = position[2]/hzed;
        ip1 = (int)ceil(ifloat);
        ip2 = ip1 + 2;
        im1 = (int)floor(ifloat);
        im2 = im1 - 2;
        jp1 = (int)ceil(jfloat);
        jp2 = jp1 + 2;
        jm1 = (int)floor(jfloat);
        jm2 = jm1 - 2;
        kp1 = (int)ceil(kfloat);
        kp2 = kp1 + 2;
        km1 = (int)floor(kfloat);
        km2 = km1 - 2;

        /* This step shouldn't be necessary, but it saves nasty debugging
         * later on if something goes wrong */
        ip2 = VMIN2(ip2,nx-1);
        ip1 = VMIN2(ip1,nx-1);
        im1 = VMAX2(im1,0);
        im2 = VMAX2(im2,0);
        jp2 = VMIN2(jp2,ny-1);
        jp1 = VMIN2(jp1,ny-1);
        jm1 = VMAX2(jm1,0);
        jm2 = VMAX2(jm2,0);
        kp2 = VMIN2(kp2,nz-1);
        kp1 = VMIN2(kp1,nz-1);
        km1 = VMAX2(km1,0);
        km2 = VMAX2(km2,0);

        for (ii=im2; ii<=ip2; ii++) {
            mi = VFCHI4(ii,ifloat);
            mx = bspline4(mi);
            dmx = dbspline4(mi);
            d2mx = d2bspline4(mi);
            d3mx = d3bspline4(mi);
            for (jj=jm2; jj<=jp2; jj++) {
                mj = VFCHI4(jj,jfloat);
                my = bspline4(mj);
                dmy = dbspline4(mj);
                d2my = d2bspline4(mj);
                d3my = d3bspline4(mj);
                for (kk=km2; kk<=kp2; kk++) {
                    mk = VFCHI4(kk,kfloat);
                    mz = bspline4(mk);
                    dmz = dbspline4(mk);
                    d2mz = d2bspline4(mk);
                    d3mz = d3bspline4(mk);
                    f = u[IJK(ii,jj,kk)];
                    /* Potential */
                    pot  += f*mx*my*mz;
                    /* Field */
                    e[0] += f*dmx*my*mz/hx;
                    e[1] += f*mx*dmy*mz/hy;
                    e[2] += f*mx*my*dmz/hzed;
                    /* Field gradient */
                    de[0][0] += f*d2mx*my*mz/(hx*hx);
                    de[1][0] += f*dmx*dmy*mz/(hy*hx);
                    de[1][1] += f*mx*d2my*mz/(hy*hy);
                    de[2][0] += f*dmx*my*dmz/(hx*hzed);
                    de[2][1] += f*mx*dmy*dmz/(hy*hzed);
                    de[2][2] += f*mx*my*d2mz/(hzed*hzed);
                    /* 2nd Field Gradient
                       VxVxVa */
                    d2e[0][0][0] += f*d3mx*my*mz /(hx*hx*hx);
                    d2e[0][0][1] += f*d2mx*dmy*mz/(hx*hy*hx);
                    d2e[0][0][2] += f*d2mx*my*dmz/(hx*hx*hzed);
                    /* VyVxVa */
                    d2e[1][0][0] += f*d2mx*dmy*mz/(hx*hx*hy);
                    d2e[1][0][1] += f*dmx*d2my*mz/(hx*hy*hy);
                    d2e[1][0][2] += f*dmx*dmy*dmz/(hx*hy*hzed);
                    /* VyVyVa */
                    d2e[1][1][0] += f*dmx*d2my*mz/(hx*hy*hy);
                    d2e[1][1][1] += f*mx*d3my*mz /(hy*hy*hy);
                    d2e[1][1][2] += f*mx*d2my*dmz/(hy*hy*hzed);
                    /* VzVxVa */
                    d2e[2][0][0] += f*d2mx*my*dmz/(hx*hx*hzed);
                    d2e[2][0][1] += f*dmx*dmy*dmz/(hx*hy*hzed);
                    d2e[2][0][2] += f*dmx*my*d2mz/(hx*hzed*hzed);
                    /* VzVyVa */
                    d2e[2][1][0] += f*dmx*dmy*dmz/(hx*hy*hzed);
                    d2e[2][1][1] += f*mx*d2my*dmz/(hy*hy*hzed);
                    d2e[2][1][2] += f*mx*dmy*d2mz/(hy*hzed*hzed);
                    /* VzVzVa */
                    d2e[2][2][0] += f*dmx*my*d2mz/(hx*hzed*hzed);
                    d2e[2][2][1] += f*mx*dmy*d2mz/(hy*hzed*hzed);
                    d2e[2][2][2] += f*mx*my*d3mz /(hzed*hzed*hzed);
                }
            }
        }
    }

    /* Monopole Force */
    force[0] = e[0]*c;
    force[1] = e[1]*c;
    force[2] = e[2]*c;

    /* Dipole Force */
    force[0] -= de[0][0]*ux+de[1][0]*uy+de[2][0]*uz;
    force[1] -= de[1][0]*ux+de[1][1]*uy+de[2][1]*uz;
    force[2] -= de[2][0]*ux+de[2][1]*uy+de[2][2]*uz;

    /* Quadrupole Force */
    force[0] += d2e[0][0][0]*qxx
             +  d2e[1][0][0]*qyx*2.0+d2e[1][1][0]*qyy
             +  d2e[2][0][0]*qzx*2.0+d2e[2][1][0]*qzy*2.0+d2e[2][2][0]*qzz;
    force[1] += d2e[0][0][1]*qxx
             +  d2e[1][0][1]*qyx*2.0+d2e[1][1][1]*qyy
             +  d2e[2][0][1]*qzx*2.0+d2e[2][1][1]*qzy*2.0+d2e[2][2][1]*qzz;
    force[2] += d2e[0][0][2]*qxx
             +  d2e[1][0][2]*qyx*2.0+d2e[1][1][2]*qyy
             +  d2e[2][0][2]*qzx*2.0+d2e[2][1][2]*qzy*2.0+d2e[2][2][2]*qzz;

    /* Dipole Torque */
    torque[0] = uy * e[2] - uz * e[1];
    torque[1] = uz * e[0] - ux * e[2];
    torque[2] = ux * e[1] - uy * e[0];
    /* Quadrupole Torque */
    de[0][1] = de[1][0];
    de[0][2] = de[2][0];
    de[1][2] = de[2][1];
    torque[0] -= 2.0*(qyx*de[0][2] + qyy*de[1][2] + qyz*de[2][2]
                    - qzx*de[0][1] - qzy*de[1][1] - qzz*de[2][1]);
    torque[1] -= 2.0*(qzx*de[0][0] + qzy*de[1][0] + qzz*de[2][0]
                    - qxx*de[0][2] - qxy*de[1][2] - qxz*de[2][2]);
    torque[2] -= 2.0*(qxx*de[0][1] + qxy*de[1][1] + qxz*de[2][1]
                    - qyx*de[0][0] - qyy*de[1][0] - qyz*de[2][0]);


    /* printf(" qPhi Force %f %f %f\n", force[0], force[1], force[2]);
       printf(" qPhi Torque %f %f %f\n", torque[0], torque[1], torque[2]); */
}

VPUBLIC void Vpmg_ibPermanentMultipoleForce(Vpmg *thee, int atomID,
                                            double force[3]) {

    Valist *alist;
    Vacc *acc;
    Vpbe *pbe;
    Vatom *atom;
    Vsurf_Meth srfm;

    /* Grid variables */
    double *apos, position[3], arad, irad, zkappa2, hx, hy, hzed;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2;
    double rtot, dx, dx2, dy, dy2, dz, dz2, gpos[3], tgrad[3], fmag;
    double izmagic;
    int i, j, k, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);

    /* Nonlinear PBE is not implemented for AMOEBA */
    VASSERT(!thee->pmgp->nonlin);

    acc = thee->pbe->acc;
    srfm = thee->surfMeth;
    atom = Valist_getAtom(thee->pbe->alist, atomID);

    /* Currently all atoms must be in the same partition. */

    VASSERT(atom->partID != 0);
    apos = Vatom_getPosition(atom);
    arad = Vatom_getRadius(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    zkappa2 = Vpbe_getZkappa2(pbe);
    izmagic = 1.0/Vpbe_getZmagic(pbe);

    /* Should be a check for this further up. */
    VASSERT (zkappa2 > VPMGSMALL);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
        Vnm_print(2, "ibPermanentMultipoleForce:  Atom %d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n", atomID, apos[0], apos[1], apos[2]);
        Vnm_print(2, "ibPermanentMultipoleForce:  xmin = %g, xmax = %g\n", xmin, xmax);
        Vnm_print(2, "ibPermanentMultipoleForce:  ymin = %g, ymax = %g\n", ymin, ymax);
        Vnm_print(2, "ibPermanentMultipoleForce:  zmin = %g, zmax = %g\n", zmin, zmax);
        fflush(stderr);
    } else {

        /* Convert the atom position to grid reference frame */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;

        /* Integrate over points within this atom's (inflated) radius */
        rtot = (irad + arad + thee->splineWin);
        rtot2 = VSQR(rtot);
        dx = rtot + 0.5*hx;
        imin = VMAX2(0,(int)ceil((position[0] - dx)/hx));
        imax = VMIN2(nx-1,(int)floor((position[0] + dx)/hx));
        for (i=imin; i<=imax; i++) {
            dx2 = VSQR(position[0] - hx*i);
            if (rtot2 > dx2) dy = VSQRT(rtot2 - dx2) + 0.5*hy;
            else dy = 0.5*hy;
            jmin = VMAX2(0,(int)ceil((position[1] - dy)/hy));
            jmax = VMIN2(ny-1,(int)floor((position[1] + dy)/hy));
            for (j=jmin; j<=jmax; j++) {
                dy2 = VSQR(position[1] - hy*j);
                if (rtot2 > (dx2+dy2)) dz = VSQRT(rtot2-dx2-dy2)+0.5*hzed;
                else dz = 0.5*hzed;
                kmin = VMAX2(0,(int)ceil((position[2] - dz)/hzed));
                kmax = VMIN2(nz-1,(int)floor((position[2] + dz)/hzed));
                for (k=kmin; k<=kmax; k++) {
                    dz2 = VSQR(k*hzed - position[2]);
                    /* See if grid point is inside ivdw radius and set ccf
                     * accordingly (do spline assignment here) */
                    if ((dz2 + dy2 + dx2) <= rtot2) {
                        gpos[0] = i*hx + xmin;
                        gpos[1] = j*hy + ymin;
                        gpos[2] = k*hzed + zmin;
                        Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, irad, atom, tgrad);
                        fmag = VSQR(thee->u[IJK(i,j,k)])*thee->kappa[IJK(i,j,k)];
                        force[0] += (zkappa2*fmag*tgrad[0]);
                        force[1] += (zkappa2*fmag*tgrad[1]);
                        force[2] += (zkappa2*fmag*tgrad[2]);
                    }
                } /* k loop */
            } /* j loop */
        } /* i loop */
    }

    force[0] = force[0] * 0.5 * hx * hy * hzed * izmagic;
    force[1] = force[1] * 0.5 * hx * hy * hzed * izmagic;
    force[2] = force[2] * 0.5 * hx * hy * hzed * izmagic;

}

VPUBLIC void Vpmg_dbPermanentMultipoleForce(Vpmg *thee, int atomID,
                                            double force[3]) {

    Vacc *acc;
    Vpbe *pbe;
    Vatom *atom;
    Vsurf_Meth srfm;

    double *apos, position[3], arad, hx, hy, hzed, izmagic, deps, depsi;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2, epsp;
    double rtot, dx, gpos[3], tgrad[3], dbFmag, epsw, kT;
    double *u, Hxijk, Hyijk, Hzijk, Hxim1jk, Hyijm1k, Hzijkm1;
    double dHxijk[3], dHyijk[3], dHzijk[3], dHxim1jk[3], dHyijm1k[3];
    double dHzijkm1[3];
    int i, j, k, l, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);

    acc = thee->pbe->acc;
    srfm = thee->surfMeth;
    atom = Valist_getAtom(thee->pbe->alist, atomID);

    /* Currently all atoms must be in the same partition. */

    VASSERT(atom->partID != 0);
    arad = Vatom_getRadius(atom);
    apos = Vatom_getPosition(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    epsp = Vpbe_getSoluteDiel(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    kT = Vpbe_getTemperature(pbe)*(1e-3)*Vunit_Na*Vunit_kb;
    izmagic = 1.0/Vpbe_getZmagic(pbe);


    deps = (epsw - epsp);
    depsi = 1.0/deps;

    VASSERT(VABS(deps) > VPMGSMALL);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = thee->u;

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
        Vnm_print(2, "dbPermanentMultipoleForce:  Atom at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n", apos[0], apos[1], apos[2]);
        Vnm_print(2, "dbPermanentMultipoleForce:  xmin = %g, xmax = %g\n", xmin, xmax);
        Vnm_print(2, "dbPermanentMultipoleForce:  ymin = %g, ymax = %g\n", ymin, ymax);
        Vnm_print(2, "dbPermanentMultipoleForce:  zmin = %g, zmax = %g\n", zmin, zmax);
        fflush(stderr);
    } else {

        /* Convert the atom position to grid reference frame */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;

        /* Integrate over points within this atom's (inflated) radius */
        rtot = (arad + thee->splineWin);
        rtot2 = VSQR(rtot);
        dx = rtot/hx;
        imin = (int)floor((position[0]-rtot)/hx);
        if (imin < 1) {
            Vnm_print(2, "dbPermanentMultipoleForce:  Atom off grid!\n");
            return;
        }
        imax = (int)ceil((position[0]+rtot)/hx);
        if (imax > (nx-2)) {
            Vnm_print(2, "dbPermanentMultipoleForce:  Atom off grid!\n");
            return;
        }
        jmin = (int)floor((position[1]-rtot)/hy);
        if (jmin < 1) {
            Vnm_print(2, "dbPermanentMultipoleForce:  Atom off grid!\n");
            return;
        }
        jmax = (int)ceil((position[1]+rtot)/hy);
        if (jmax > (ny-2)) {
            Vnm_print(2, "dbPermanentMultipoleForce:  Atom off grid!\n");
            return;
        }
        kmin = (int)floor((position[2]-rtot)/hzed);
        if (kmin < 1) {
            Vnm_print(2, "dbPermanentMultipoleForce:  Atom off grid!\n");
            return;
        }
        kmax = (int)ceil((position[2]+rtot)/hzed);
        if (kmax > (nz-2)) {
            Vnm_print(2, "dbPermanentMultipoleForce:  Atom off grid!\n");
            return;
        }
        for (i=imin; i<=imax; i++) {
            for (j=jmin; j<=jmax; j++) {
                for (k=kmin; k<=kmax; k++) {
                    /* i,j,k */
                    gpos[0] = (i+0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxijk = (thee->epsx[IJK(i,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHxijk);
                    for (l=0; l<3; l++) dHxijk[l] *= Hxijk;
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j+0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijk = (thee->epsy[IJK(i,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHyijk);
                    for (l=0; l<3; l++) dHyijk[l] *= Hyijk;
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k+0.5)*hzed + zmin;
                    Hzijk = (thee->epsz[IJK(i,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHzijk);
                    for (l=0; l<3; l++) dHzijk[l] *= Hzijk;
                    /* i-1,j,k */
                    gpos[0] = (i-0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxim1jk = (thee->epsx[IJK(i-1,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHxim1jk);
                    for (l=0; l<3; l++) dHxim1jk[l] *= Hxim1jk;
                    /* i,j-1,k */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j-0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijm1k = (thee->epsy[IJK(i,j-1,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHyijm1k);
                    for (l=0; l<3; l++) dHyijm1k[l] *= Hyijm1k;
                    /* i,j,k-1 */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k-0.5)*hzed + zmin;
                    Hzijkm1 = (thee->epsz[IJK(i,j,k-1)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHzijkm1);
                    for (l=0; l<3; l++) dHzijkm1[l] *= Hzijkm1;
                    dbFmag = u[IJK(i,j,k)];
                    tgrad[0] =
                       (dHxijk[0]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[0]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[0]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[0]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[0]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[0]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    tgrad[1] =
                       (dHxijk[1]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[1]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[1]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[1]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[1]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[1]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    tgrad[2] =
                       (dHxijk[2]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[2]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[2]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[2]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[2]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[2]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                     force[0] += (dbFmag*tgrad[0]);
                     force[1] += (dbFmag*tgrad[1]);
                     force[2] += (dbFmag*tgrad[2]);
                } /* k loop */
            } /* j loop */
        } /* i loop */
        force[0] = -force[0]*hx*hy*hzed*deps*0.5*izmagic;
        force[1] = -force[1]*hx*hy*hzed*deps*0.5*izmagic;
        force[2] = -force[2]*hx*hy*hzed*deps*0.5*izmagic;
    }
}

VPUBLIC void Vpmg_qfDirectPolForce(Vpmg *thee, Vgrid* perm, Vgrid *induced,
                                   int atomID, double force[3], double torque[3]) {

    Vatom *atom;
    Vpbe *pbe;
    double f, fp, *u, *up, *apos, position[3];

    /* Grid variables */
    int nx,ny,nz;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax;
    double hx, hy, hzed, ifloat, jfloat, kfloat;

    /* B-spline weights */
    double mx, my, mz, dmx, dmy, dmz, d2mx, d2my, d2mz, d3mx, d3my, d3mz;
    double mi, mj, mk;

    /* Loop indeces */
    int i, j, k, ii, jj, kk;
    int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1, kp1, kp2;

    /* Permanent potential, field, field gradient and 2nd field gradient */
    double pot, e[3], de[3][3], d2e[3][3][3];
    /* Induced dipole field */
    double dep[3][3];

    /* Permanent multipole components */
    double *dipole, *quad;
    double c, ux, uy, uz, qxx, qxy, qxz, qyx, qyy, qyz, qzx, qzy, qzz;
    double uix, uiy, uiz;

    VASSERT(thee != VNULL);
    VASSERT(induced != VNULL); /* the potential due to permanent multipoles.*/
    VASSERT(induced != VNULL); /* the potential due to local induced dipoles.*/
    VASSERT(thee->pbe != VNULL);
    VASSERT(thee->pbe->alist != VNULL);

    atom = Valist_getAtom(thee->pbe->alist, atomID);
    VASSERT(atom->partID != 0); /* all atoms must be in the same partition.*/
    apos = Vatom_getPosition(atom);

    c = Vatom_getCharge(atom);
    dipole = Vatom_getDipole(atom);
    ux = dipole[0];
    uy = dipole[1];
    uz = dipole[2];
    quad = Vatom_getQuadrupole(atom);
    qxx = quad[0]/3.0;
    qxy = quad[1]/3.0;
    qxz = quad[2]/3.0;
    qyx = quad[3]/3.0;
    qyy = quad[4]/3.0;
    qyz = quad[5]/3.0;
    qzx = quad[6]/3.0;
    qzy = quad[7]/3.0;
    qzz = quad[8]/3.0;

    dipole = Vatom_getInducedDipole(atom);
    uix = dipole[0];
    uiy = dipole[1];
    uiz = dipole[2];

    /* Reset Field Gradients */
    pot = 0.0;
    for (i=0;i<3;i++){
       e[i] = 0.0;
       for (j=0;j<3;j++){
          de[i][j] = 0.0;
          dep[i][j] = 0.0;
          for (k=0;k<3;k++){
             d2e[i][j][k] = 0.0;
          }
       }
    }

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = induced->data;
    up = perm->data;

    /* Make sure we're on the grid */
    if ((apos[0]<=(xmin+2*hx))   || (apos[0]>=(xmax-2*hx)) \
     || (apos[1]<=(ymin+2*hy))   || (apos[1]>=(ymax-2*hy)) \
     || (apos[2]<=(zmin+2*hzed)) || (apos[2]>=(zmax-2*hzed))) {
        Vnm_print(2, "qfDirectPolForce:  Atom off the mesh (ignoring) %6.3f %6.3f %6.3f\n", apos[0], apos[1], apos[2]);
        fflush(stderr);

    } else {

        /* Convert the atom position to grid coordinates */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;
        ifloat = position[0]/hx;
        jfloat = position[1]/hy;
        kfloat = position[2]/hzed;
        ip1 = (int)ceil(ifloat);
        ip2 = ip1 + 2;
        im1 = (int)floor(ifloat);
        im2 = im1 - 2;
        jp1 = (int)ceil(jfloat);
        jp2 = jp1 + 2;
        jm1 = (int)floor(jfloat);
        jm2 = jm1 - 2;
        kp1 = (int)ceil(kfloat);
        kp2 = kp1 + 2;
        km1 = (int)floor(kfloat);
        km2 = km1 - 2;

        /* This step shouldn't be necessary, but it saves nasty debugging
         * later on if something goes wrong */
        ip2 = VMIN2(ip2,nx-1);
        ip1 = VMIN2(ip1,nx-1);
        im1 = VMAX2(im1,0);
        im2 = VMAX2(im2,0);
        jp2 = VMIN2(jp2,ny-1);
        jp1 = VMIN2(jp1,ny-1);
        jm1 = VMAX2(jm1,0);
        jm2 = VMAX2(jm2,0);
        kp2 = VMIN2(kp2,nz-1);
        kp1 = VMIN2(kp1,nz-1);
        km1 = VMAX2(km1,0);
        km2 = VMAX2(km2,0);

        for (ii=im2; ii<=ip2; ii++) {
            mi = VFCHI4(ii,ifloat);
            mx = bspline4(mi);
            dmx = dbspline4(mi);
            d2mx = d2bspline4(mi);
            d3mx = d3bspline4(mi);
            for (jj=jm2; jj<=jp2; jj++) {
                mj = VFCHI4(jj,jfloat);
                my = bspline4(mj);
                dmy = dbspline4(mj);
                d2my = d2bspline4(mj);
                d3my = d3bspline4(mj);
                for (kk=km2; kk<=kp2; kk++) {
                    mk = VFCHI4(kk,kfloat);
                    mz = bspline4(mk);
                    dmz = dbspline4(mk);
                    d2mz = d2bspline4(mk);
                    d3mz = d3bspline4(mk);
                    f = u[IJK(ii,jj,kk)];
                    fp = up[IJK(ii,jj,kk)];
                    /* The potential */
                    pot  += f*mx*my*mz;
                    /* The field */
                    e[0] += f*dmx*my*mz/hx;
                    e[1] += f*mx*dmy*mz/hy;
                    e[2] += f*mx*my*dmz/hzed;
                    /* The gradient of the field */
                    de[0][0] += f*d2mx*my*mz/(hx*hx);
                    de[1][0] += f*dmx*dmy*mz/(hy*hx);
                    de[1][1] += f*mx*d2my*mz/(hy*hy);
                    de[2][0] += f*dmx*my*dmz/(hx*hzed);
                    de[2][1] += f*mx*dmy*dmz/(hy*hzed);
                    de[2][2] += f*mx*my*d2mz/(hzed*hzed);
                    /* The gradient of the (permanent) field */
                    dep[0][0] += fp*d2mx*my*mz/(hx*hx);
                    dep[1][0] += fp*dmx*dmy*mz/(hy*hx);
                    dep[1][1] += fp*mx*d2my*mz/(hy*hy);
                    dep[2][0] += fp*dmx*my*dmz/(hx*hzed);
                    dep[2][1] += fp*mx*dmy*dmz/(hy*hzed);
                    dep[2][2] += fp*mx*my*d2mz/(hzed*hzed);
                    /* The 2nd gradient of the field
                       VxVxVa */
                    d2e[0][0][0] += f*d3mx*my*mz /(hx*hx*hx);
                    d2e[0][0][1] += f*d2mx*dmy*mz/(hx*hy*hx);
                    d2e[0][0][2] += f*d2mx*my*dmz/(hx*hx*hzed);
                    /* VyVxVa */
                    d2e[1][0][0] += f*d2mx*dmy*mz/(hx*hx*hy);
                    d2e[1][0][1] += f*dmx*d2my*mz/(hx*hy*hy);
                    d2e[1][0][2] += f*dmx*dmy*dmz/(hx*hy*hzed);
                    /* VyVyVa */
                    d2e[1][1][0] += f*dmx*d2my*mz/(hx*hy*hy);
                    d2e[1][1][1] += f*mx*d3my*mz /(hy*hy*hy);
                    d2e[1][1][2] += f*mx*d2my*dmz/(hy*hy*hzed);
                    /* VzVxVa */
                    d2e[2][0][0] += f*d2mx*my*dmz/(hx*hx*hzed);
                    d2e[2][0][1] += f*dmx*dmy*dmz/(hx*hy*hzed);
                    d2e[2][0][2] += f*dmx*my*d2mz/(hx*hzed*hzed);
                    /* VzVyVa */
                    d2e[2][1][0] += f*dmx*dmy*dmz/(hx*hy*hzed);
                    d2e[2][1][1] += f*mx*d2my*dmz/(hy*hy*hzed);
                    d2e[2][1][2] += f*mx*dmy*d2mz/(hy*hzed*hzed);
                    /* VzVzVa */
                    d2e[2][2][0] += f*dmx*my*d2mz/(hx*hzed*hzed);
                    d2e[2][2][1] += f*mx*dmy*d2mz/(hy*hzed*hzed);
                    d2e[2][2][2] += f*mx*my*d3mz /(hzed*hzed*hzed);
                }
            }
        }
    }

    /* force on permanent multipole due to induced reaction field */

    /* Monopole Force */
    force[0] = e[0]*c;
    force[1] = e[1]*c;
    force[2] = e[2]*c;

    /* Dipole Force */
    force[0] -= de[0][0]*ux+de[1][0]*uy+de[2][0]*uz;
    force[1] -= de[1][0]*ux+de[1][1]*uy+de[2][1]*uz;
    force[2] -= de[2][0]*ux+de[2][1]*uy+de[2][2]*uz;

    /* Quadrupole Force */
    force[0] += d2e[0][0][0]*qxx
             +  d2e[1][0][0]*qyx*2.0+d2e[1][1][0]*qyy
             +  d2e[2][0][0]*qzx*2.0+d2e[2][1][0]*qzy*2.0+d2e[2][2][0]*qzz;
    force[1] += d2e[0][0][1]*qxx
             +  d2e[1][0][1]*qyx*2.0+d2e[1][1][1]*qyy
             +  d2e[2][0][1]*qzx*2.0+d2e[2][1][1]*qzy*2.0+d2e[2][2][1]*qzz;
    force[2] += d2e[0][0][2]*qxx
             +  d2e[1][0][2]*qyx*2.0+d2e[1][1][2]*qyy
             +  d2e[2][0][2]*qzx*2.0+d2e[2][1][2]*qzy*2.0+d2e[2][2][2]*qzz;

    /* torque on permanent mulitpole due to induced reaction field */

    /* Dipole Torque */
    torque[0] = uy * e[2] - uz * e[1];
    torque[1] = uz * e[0] - ux * e[2];
    torque[2] = ux * e[1] - uy * e[0];

    /* Quadrupole Torque */
    /* Tx = -2.0*(Sum_a (Qya*dEaz) + Sum_b (Qzb*dEby))
       Ty = -2.0*(Sum_a (Qza*dEax) + Sum_b (Qxb*dEbz))
       Tz = -2.0*(Sum_a (Qxa*dEay) + Sum_b (Qyb*dEbx))  */
    de[0][1] = de[1][0];
    de[0][2] = de[2][0];
    de[1][2] = de[2][1];
    torque[0] -= 2.0*(qyx*de[0][2] + qyy*de[1][2] + qyz*de[2][2]
                    - qzx*de[0][1] - qzy*de[1][1] - qzz*de[2][1]);
    torque[1] -= 2.0*(qzx*de[0][0] + qzy*de[1][0] + qzz*de[2][0]
                    - qxx*de[0][2] - qxy*de[1][2] - qxz*de[2][2]);
    torque[2] -= 2.0*(qxx*de[0][1] + qxy*de[1][1] + qxz*de[2][1]
                    - qyx*de[0][0] - qyy*de[1][0] - qyz*de[2][0]);

    /* force on induced dipole due to permanent reaction field */

    force[0] -= dep[0][0]*uix+dep[1][0]*uiy+dep[2][0]*uiz;
    force[1] -= dep[1][0]*uix+dep[1][1]*uiy+dep[2][1]*uiz;
    force[2] -= dep[2][0]*uix+dep[2][1]*uiy+dep[2][2]*uiz;

    force[0] = 0.5 * force[0];
    force[1] = 0.5 * force[1];
    force[2] = 0.5 * force[2];
    torque[0] = 0.5 * torque[0];
    torque[1] = 0.5 * torque[1];
    torque[2] = 0.5 * torque[2];

    /* printf(" qPhi Force %f %f %f\n", force[0], force[1], force[2]);
      printf(" qPhi Torque %f %f %f\n", torque[0], torque[1], torque[2]); */
}

VPUBLIC void Vpmg_qfNLDirectPolForce(Vpmg *thee, Vgrid *perm, Vgrid *nlInduced,
                                     int atomID, double force[3], double torque[3]) {

    Vatom *atom;
    double *apos, *dipole, *quad, position[3], hx, hy, hzed;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax;
    double pot, e[3],de[3][3],dep[3][3],d2e[3][3][3];
    double mx, my, mz, dmx, dmy, dmz, mi, mj, mk;
    double d2mx, d2my, d2mz, d3mx, d3my, d3mz;
    double *u, *up, charge, ifloat, jfloat, kfloat;
    double f, fp, c, ux, uy, uz, qxx, qxy, qxz, qyx, qyy, qyz, qzx, qzy, qzz;
    double uix, uiy, uiz;
    int i,j,k,nx, ny, nz, im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1;
    int kp1, kp2, ii, jj, kk;

    VASSERT(thee != VNULL);
    VASSERT(perm != VNULL);      /* potential due to permanent multipoles. */
    VASSERT(nlInduced != VNULL); /* potential due to non-local induced dipoles */
    VASSERT(!thee->pmgp->nonlin); /* Nonlinear PBE is not implemented for AMOEBA */

    atom = Valist_getAtom(thee->pbe->alist, atomID);
    VASSERT(atom->partID != 0);   /* Currently all atoms must be in the same partition. */
    apos = Vatom_getPosition(atom);

    c = Vatom_getCharge(atom);
    dipole = Vatom_getDipole(atom);
    ux = dipole[0];
    uy = dipole[1];
    uz = dipole[2];
    quad = Vatom_getQuadrupole(atom);
    qxx = quad[0]/3.0;
    qxy = quad[1]/3.0;
    qxz = quad[2]/3.0;
    qyx = quad[3]/3.0;
    qyy = quad[4]/3.0;
    qyz = quad[5]/3.0;
    qzx = quad[6]/3.0;
    qzy = quad[7]/3.0;
    qzz = quad[8]/3.0;

    dipole = Vatom_getNLInducedDipole(atom);
    uix = dipole[0];
    uiy = dipole[1];
    uiz = dipole[2];

    /* Reset Field Gradients */
    pot = 0.0;
    for (i=0;i<3;i++){
       e[i] = 0.0;
       for (j=0;j<3;j++){
          de[i][j] = 0.0;
          dep[i][j] = 0.0;
          for (k=0;k<3;k++){
             d2e[i][j][k] = 0.0;
          }
       }
    }

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = nlInduced->data;
    up = perm->data;


    /* Make sure we're on the grid */
    if ((apos[0]<=(xmin+2*hx))   || (apos[0]>=(xmax-2*hx)) \
     || (apos[1]<=(ymin+2*hy))   || (apos[1]>=(ymax-2*hy)) \
     || (apos[2]<=(zmin+2*hzed)) || (apos[2]>=(zmax-2*hzed))) {
        Vnm_print(2, "qfNLDirectMultipoleForce:  Atom off the mesh (ignoring) %6.3f %6.3f %6.3f\n", apos[0], apos[1], apos[2]);
    } else {

        /* Convert the atom position to grid coordinates */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;
        ifloat = position[0]/hx;
        jfloat = position[1]/hy;
        kfloat = position[2]/hzed;
        ip1 = (int)ceil(ifloat);
        ip2 = ip1 + 2;
        im1 = (int)floor(ifloat);
        im2 = im1 - 2;
        jp1 = (int)ceil(jfloat);
        jp2 = jp1 + 2;
        jm1 = (int)floor(jfloat);
        jm2 = jm1 - 2;
        kp1 = (int)ceil(kfloat);
        kp2 = kp1 + 2;
        km1 = (int)floor(kfloat);
        km2 = km1 - 2;

        /* This step shouldn't be necessary, but it saves nasty debugging
         * later on if something goes wrong */
        ip2 = VMIN2(ip2,nx-1);
        ip1 = VMIN2(ip1,nx-1);
        im1 = VMAX2(im1,0);
        im2 = VMAX2(im2,0);
        jp2 = VMIN2(jp2,ny-1);
        jp1 = VMIN2(jp1,ny-1);
        jm1 = VMAX2(jm1,0);
        jm2 = VMAX2(jm2,0);
        kp2 = VMIN2(kp2,nz-1);
        kp1 = VMIN2(kp1,nz-1);
        km1 = VMAX2(km1,0);
        km2 = VMAX2(km2,0);

        for (ii=im2; ii<=ip2; ii++) {
            mi = VFCHI4(ii,ifloat);
            mx = bspline4(mi);
            dmx = dbspline4(mi);
            d2mx = d2bspline4(mi);
            d3mx = d3bspline4(mi);
            for (jj=jm2; jj<=jp2; jj++) {
                mj = VFCHI4(jj,jfloat);
                my = bspline4(mj);
                dmy = dbspline4(mj);
                d2my = d2bspline4(mj);
                d3my = d3bspline4(mj);
                for (kk=km2; kk<=kp2; kk++) {
                    mk = VFCHI4(kk,kfloat);
                    mz = bspline4(mk);
                    dmz = dbspline4(mk);
                    d2mz = d2bspline4(mk);
                    d3mz = d3bspline4(mk);
                    f = u[IJK(ii,jj,kk)];
                    fp = up[IJK(ii,jj,kk)];
                    /* The potential */
                    pot  += f*mx*my*mz;
                    /* The field */
                    e[0] += f*dmx*my*mz/hx;
                    e[1] += f*mx*dmy*mz/hy;
                    e[2] += f*mx*my*dmz/hzed;
                    /* The gradient of the field */
                    de[0][0] += f*d2mx*my*mz/(hx*hx);
                    de[1][0] += f*dmx*dmy*mz/(hy*hx);
                    de[1][1] += f*mx*d2my*mz/(hy*hy);
                    de[2][0] += f*dmx*my*dmz/(hx*hzed);
                    de[2][1] += f*mx*dmy*dmz/(hy*hzed);
                    de[2][2] += f*mx*my*d2mz/(hzed*hzed);
                    /* The gradient of the (permanent) field */
                    dep[0][0] += fp*d2mx*my*mz/(hx*hx);
                    dep[1][0] += fp*dmx*dmy*mz/(hy*hx);
                    dep[1][1] += fp*mx*d2my*mz/(hy*hy);
                    dep[2][0] += fp*dmx*my*dmz/(hx*hzed);
                    dep[2][1] += fp*mx*dmy*dmz/(hy*hzed);
                    dep[2][2] += fp*mx*my*d2mz/(hzed*hzed);
                    /* The 2nd gradient of the field */
                    /* VxVxVa */
                    d2e[0][0][0] += f*d3mx*my*mz /(hx*hx*hx);
                    d2e[0][0][1] += f*d2mx*dmy*mz/(hx*hy*hx);
                    d2e[0][0][2] += f*d2mx*my*dmz/(hx*hx*hzed);
                    /* VyVxVa */
                    d2e[1][0][0] += f*d2mx*dmy*mz/(hx*hx*hy);
                    d2e[1][0][1] += f*dmx*d2my*mz/(hx*hy*hy);
                    d2e[1][0][2] += f*dmx*dmy*dmz/(hx*hy*hzed);
                    /* VyVyVa */
                    d2e[1][1][0] += f*dmx*d2my*mz/(hx*hy*hy);
                    d2e[1][1][1] += f*mx*d3my*mz /(hy*hy*hy);
                    d2e[1][1][2] += f*mx*d2my*dmz/(hy*hy*hzed);
                    /* VzVxVa */
                    d2e[2][0][0] += f*d2mx*my*dmz/(hx*hx*hzed);
                    d2e[2][0][1] += f*dmx*dmy*dmz/(hx*hy*hzed);
                    d2e[2][0][2] += f*dmx*my*d2mz/(hx*hzed*hzed);
                    /* VzVyVa */
                    d2e[2][1][0] += f*dmx*dmy*dmz/(hx*hy*hzed);
                    d2e[2][1][1] += f*mx*d2my*dmz/(hy*hy*hzed);
                    d2e[2][1][2] += f*mx*dmy*d2mz/(hy*hzed*hzed);
                    /* VzVzVa */
                    d2e[2][2][0] += f*dmx*my*d2mz/(hx*hzed*hzed);
                    d2e[2][2][1] += f*mx*dmy*d2mz/(hy*hzed*hzed);
                    d2e[2][2][2] += f*mx*my*d3mz /(hzed*hzed*hzed);
                }
            }
        }
    }

    /* force on permanent multipole due to non-local induced reaction field */

    /* Monopole Force */
    force[0] = e[0]*c;
    force[1] = e[1]*c;
    force[2] = e[2]*c;

    /* Dipole Force */
    force[0] -= de[0][0]*ux+de[1][0]*uy+de[2][0]*uz;
    force[1] -= de[1][0]*ux+de[1][1]*uy+de[2][1]*uz;
    force[2] -= de[2][0]*ux+de[2][1]*uy+de[2][2]*uz;

    /* Quadrupole Force */
    force[0] += d2e[0][0][0]*qxx
             +  d2e[1][0][0]*qyx*2.0+d2e[1][1][0]*qyy
             +  d2e[2][0][0]*qzx*2.0+d2e[2][1][0]*qzy*2.0+d2e[2][2][0]*qzz;
    force[1] += d2e[0][0][1]*qxx
             +  d2e[1][0][1]*qyx*2.0+d2e[1][1][1]*qyy
             +  d2e[2][0][1]*qzx*2.0+d2e[2][1][1]*qzy*2.0+d2e[2][2][1]*qzz;
    force[2] += d2e[0][0][2]*qxx
             +  d2e[1][0][2]*qyx*2.0+d2e[1][1][2]*qyy
             +  d2e[2][0][2]*qzx*2.0+d2e[2][1][2]*qzy*2.0+d2e[2][2][2]*qzz;

    /* torque on permanent mulitpole due to non-local induced reaction field */

    /* Dipole Torque */
    torque[0] = uy * e[2] - uz * e[1];
    torque[1] = uz * e[0] - ux * e[2];
    torque[2] = ux * e[1] - uy * e[0];

    /* Quadrupole Torque */
    /* Tx = -2.0*(Sum_a (Qya*dEaz) + Sum_b (Qzb*dEby))
       Ty = -2.0*(Sum_a (Qza*dEax) + Sum_b (Qxb*dEbz))
       Tz = -2.0*(Sum_a (Qxa*dEay) + Sum_b (Qyb*dEbx))  */
    de[0][1] = de[1][0];
    de[0][2] = de[2][0];
    de[1][2] = de[2][1];
    torque[0] -= 2.0*(qyx*de[0][2] + qyy*de[1][2] + qyz*de[2][2]
                    - qzx*de[0][1] - qzy*de[1][1] - qzz*de[2][1]);
    torque[1] -= 2.0*(qzx*de[0][0] + qzy*de[1][0] + qzz*de[2][0]
                    - qxx*de[0][2] - qxy*de[1][2] - qxz*de[2][2]);
    torque[2] -= 2.0*(qxx*de[0][1] + qxy*de[1][1] + qxz*de[2][1]
                    - qyx*de[0][0] - qyy*de[1][0] - qyz*de[2][0]);

    /* force on non-local induced dipole due to permanent reaction field */

    force[0] -= dep[0][0]*uix+dep[1][0]*uiy+dep[2][0]*uiz;
    force[1] -= dep[1][0]*uix+dep[1][1]*uiy+dep[2][1]*uiz;
    force[2] -= dep[2][0]*uix+dep[2][1]*uiy+dep[2][2]*uiz;

    force[0] = 0.5 * force[0];
    force[1] = 0.5 * force[1];
    force[2] = 0.5 * force[2];
    torque[0] = 0.5 * torque[0];
    torque[1] = 0.5 * torque[1];
    torque[2] = 0.5 * torque[2];

    /* printf(" qPhi Force %f %f %f\n", force[0], force[1], force[2]);
       printf(" qPhi Torque %f %f %f\n", torque[0], torque[1], torque[2]); */
}

VPUBLIC void Vpmg_ibDirectPolForce(Vpmg *thee, Vgrid *perm, Vgrid *induced,
                                   int atomID, double force[3]) {

    Vatom *atom;
    Valist *alist;
    Vacc *acc;
    Vpbe *pbe;
    Vsurf_Meth srfm;

    double *apos, position[3], arad, irad, zkappa2, hx, hy, hzed;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2;
    double rtot, dx, dx2, dy, dy2, dz, dz2, gpos[3], tgrad[3], fmag;
    double izmagic;
    int i, j, k, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);
    VASSERT(perm != VNULL);        /* potential due to permanent multipoles.*/
    VASSERT(induced != VNULL);     /* potential due to induced dipoles. */
    VASSERT (!thee->pmgp->nonlin); /* Nonlinear PBE is not implemented for AMOEBA */

    acc = thee->pbe->acc;
    srfm = thee->surfMeth;
    atom = Valist_getAtom(thee->pbe->alist, atomID);
    VASSERT(atom->partID != 0);   /* Currently all atoms must be in the same partition. */
    apos = Vatom_getPosition(atom);
    arad = Vatom_getRadius(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    zkappa2 = Vpbe_getZkappa2(pbe);
    izmagic = 1.0/Vpbe_getZmagic(pbe);

    VASSERT (zkappa2 > VPMGSMALL); /* It is ok to run AMOEBA with no ions, but this is checked for higher up in the driver. */

    /* Mesh info */
    nx = induced->nx;
    ny = induced->ny;
    nz = induced->nz;
    hx = induced->hx;
    hy = induced->hy;
    hzed = induced->hzed;
    xmin = induced->xmin;
    ymin = induced->ymin;
    zmin = induced->zmin;
    xmax = induced->xmax;
    ymax = induced->ymax;
    zmax = induced->zmax;
    xlen = xmax-xmin;
    ylen = ymax-ymin;
    zlen = zmax-zmin;

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
        Vnm_print(2, "Vpmg_ibForce:  Atom at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n",
            apos[0], apos[1], apos[2]);
        Vnm_print(2, "Vpmg_ibForce:    xmin = %g, xmax = %g\n", xmin, xmax);
        Vnm_print(2, "Vpmg_ibForce:    ymin = %g, ymax = %g\n", ymin, ymax);
        Vnm_print(2, "Vpmg_ibForce:    zmin = %g, zmax = %g\n", zmin, zmax);
        fflush(stderr);
    } else {

        /* Convert the atom position to grid reference frame */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;

        /* Integrate over points within this atom's (inflated) radius */
        rtot = (irad + arad + thee->splineWin);
        rtot2 = VSQR(rtot);
        dx = rtot + 0.5*hx;
        imin = VMAX2(0,(int)ceil((position[0] - dx)/hx));
        imax = VMIN2(nx-1,(int)floor((position[0] + dx)/hx));
        for (i=imin; i<=imax; i++) {
            dx2 = VSQR(position[0] - hx*i);
            if (rtot2 > dx2) dy = VSQRT(rtot2 - dx2) + 0.5*hy;
            else dy = 0.5*hy;
            jmin = VMAX2(0,(int)ceil((position[1] - dy)/hy));
            jmax = VMIN2(ny-1,(int)floor((position[1] + dy)/hy));
            for (j=jmin; j<=jmax; j++) {
                dy2 = VSQR(position[1] - hy*j);
                if (rtot2 > (dx2+dy2)) dz = VSQRT(rtot2-dx2-dy2)+0.5*hzed;
                else dz = 0.5*hzed;
                kmin = VMAX2(0,(int)ceil((position[2] - dz)/hzed));
                kmax = VMIN2(nz-1,(int)floor((position[2] + dz)/hzed));
                for (k=kmin; k<=kmax; k++) {
                    dz2 = VSQR(k*hzed - position[2]);
                    /* See if grid point is inside ivdw radius and set ccf
                     * accordingly (do spline assignment here) */
                    if ((dz2 + dy2 + dx2) <= rtot2) {
                        gpos[0] = i*hx + xmin;
                        gpos[1] = j*hy + ymin;
                        gpos[2] = k*hzed + zmin;
                        Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, irad,
                          atom, tgrad);
                        fmag = induced->data[IJK(i,j,k)];
                        fmag *= perm->data[IJK(i,j,k)];
                        fmag *= thee->kappa[IJK(i,j,k)];
                        force[0] += (zkappa2*fmag*tgrad[0]);
                        force[1] += (zkappa2*fmag*tgrad[1]);
                        force[2] += (zkappa2*fmag*tgrad[2]);
                    }
                } /* k loop */
            } /* j loop */
        } /* i loop */
    }

    force[0] = force[0] * 0.5 * hx * hy * hzed * izmagic;
    force[1] = force[1] * 0.5 * hx * hy * hzed * izmagic;
    force[2] = force[2] * 0.5 * hx * hy * hzed * izmagic;

}

VPUBLIC void Vpmg_ibNLDirectPolForce(Vpmg *thee, Vgrid *perm, Vgrid *nlInduced,
                                     int atomID, double force[3]) {
     Vpmg_ibDirectPolForce(thee, perm, nlInduced, atomID, force);
}

VPUBLIC void Vpmg_dbDirectPolForce(Vpmg *thee, Vgrid *perm, Vgrid *induced,
                                   int atomID, double force[3]) {

    Vatom *atom;
    Vacc *acc;
    Vpbe *pbe;
    Vsurf_Meth srfm;

    double *apos, position[3], arad, hx, hy, hzed, izmagic, deps, depsi;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2, epsp;
    double rtot, dx, gpos[3], tgrad[3], dbFmag, epsw, kT;
    double *u, *up, Hxijk, Hyijk, Hzijk, Hxim1jk, Hyijm1k, Hzijkm1;
    double dHxijk[3], dHyijk[3], dHzijk[3], dHxim1jk[3], dHyijm1k[3];
    double dHzijkm1[3];
    int i, j, k, l, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);
    VASSERT(perm != VNULL);    /* permanent multipole PMG solution. */
    VASSERT(induced != VNULL); /* potential due to induced dipoles. */

    acc = thee->pbe->acc;
    atom = Valist_getAtom(thee->pbe->alist, atomID);
    VASSERT (atom->partID != 0);   /* Currently all atoms must be in the same partition. */
    apos = Vatom_getPosition(atom);
    arad = Vatom_getRadius(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    srfm = thee->surfMeth;
    epsp = Vpbe_getSoluteDiel(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    kT = Vpbe_getTemperature(pbe)*(1e-3)*Vunit_Na*Vunit_kb;
    izmagic = 1.0/Vpbe_getZmagic(pbe);

    deps = (epsw - epsp);
    depsi = 1.0/deps;
    VASSERT(VABS(deps) > VPMGSMALL);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    /* If the permanent and induced potentials are flipped the
       results are exactly the same. */
    u = induced->data;
    up = perm->data;

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
         Vnm_print(2, "Vpmg_dbDirectPolForce:  Atom at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n", apos[0], apos[1], apos[2]);
         Vnm_print(2, "Vpmg_dbDirectPolForce:    xmin = %g, xmax = %g\n", xmin, xmax);
         Vnm_print(2, "Vpmg_dbDirectPolForce:    ymin = %g, ymax = %g\n", ymin, ymax);
         Vnm_print(2, "Vpmg_dbDirectPolForce:    zmin = %g, zmax = %g\n", zmin, zmax);
         fflush(stderr);
    } else {

        /* Convert the atom position to grid reference frame */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;

        /* Integrate over points within this atom's (inflated) radius */
        rtot = (arad + thee->splineWin);
        rtot2 = VSQR(rtot);
        dx = rtot/hx;
        imin = (int)floor((position[0]-rtot)/hx);
        if (imin < 1) {
            Vnm_print(2, "Vpmg_dbDirectPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        imax = (int)ceil((position[0]+rtot)/hx);
        if (imax > (nx-2)) {
            Vnm_print(2, "Vpmg_dbDirectPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        jmin = (int)floor((position[1]-rtot)/hy);
        if (jmin < 1) {
            Vnm_print(2, "Vpmg_dbDirectPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        jmax = (int)ceil((position[1]+rtot)/hy);
        if (jmax > (ny-2)) {
            Vnm_print(2, "Vpmg_dbDirectPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        kmin = (int)floor((position[2]-rtot)/hzed);
        if (kmin < 1) {
            Vnm_print(2, "Vpmg_dbDirectPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        kmax = (int)ceil((position[2]+rtot)/hzed);
        if (kmax > (nz-2)) {
            Vnm_print(2, "Vpmg_dbDirectPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        for (i=imin; i<=imax; i++) {
            for (j=jmin; j<=jmax; j++) {
                for (k=kmin; k<=kmax; k++) {
                    /* i,j,k */
                    gpos[0] = (i+0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxijk = (thee->epsx[IJK(i,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHxijk);
                    for (l=0; l<3; l++) dHxijk[l] *= Hxijk;
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j+0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijk = (thee->epsy[IJK(i,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHyijk);
                    for (l=0; l<3; l++) dHyijk[l] *= Hyijk;
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k+0.5)*hzed + zmin;
                    Hzijk = (thee->epsz[IJK(i,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHzijk);
                    for (l=0; l<3; l++) dHzijk[l] *= Hzijk;
                    /* i-1,j,k */
                    gpos[0] = (i-0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxim1jk = (thee->epsx[IJK(i-1,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHxim1jk);
                    for (l=0; l<3; l++) dHxim1jk[l] *= Hxim1jk;
                    /* i,j-1,k */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j-0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijm1k = (thee->epsy[IJK(i,j-1,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHyijm1k);
                    for (l=0; l<3; l++) dHyijm1k[l] *= Hyijm1k;
                    /* i,j,k-1 */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k-0.5)*hzed + zmin;
                    Hzijkm1 = (thee->epsz[IJK(i,j,k-1)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHzijkm1);
                    for (l=0; l<3; l++) dHzijkm1[l] *= Hzijkm1;

                    dbFmag = up[IJK(i,j,k)];
                    tgrad[0] =
                       (dHxijk[0]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[0]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[0]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[0]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[0]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[0]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    tgrad[1] =
                       (dHxijk[1]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[1]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[1]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[1]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[1]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[1]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    tgrad[2] =
                       (dHxijk[2]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[2]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[2]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[2]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[2]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[2]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                     force[0] += (dbFmag*tgrad[0]);
                     force[1] += (dbFmag*tgrad[1]);
                     force[2] += (dbFmag*tgrad[2]);

                } /* k loop */
            } /* j loop */
        } /* i loop */

        force[0] = -force[0]*hx*hy*hzed*deps*0.5*izmagic;
        force[1] = -force[1]*hx*hy*hzed*deps*0.5*izmagic;
        force[2] = -force[2]*hx*hy*hzed*deps*0.5*izmagic;

    }
}

VPUBLIC void Vpmg_dbNLDirectPolForce(Vpmg *thee, Vgrid *perm, Vgrid *nlInduced,
                                     int atomID, double force[3]) {
     Vpmg_dbDirectPolForce(thee, perm, nlInduced, atomID, force);
}

VPUBLIC void Vpmg_qfMutualPolForce(Vpmg *thee, Vgrid *induced,
                              Vgrid *nlinduced, int atomID, double force[3]) {

    Vatom *atom;
    double *apos, *dipole, position[3], hx, hy, hzed;
    double *u, *unl;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax;
    double de[3][3], denl[3][3];
    double mx, my, mz, dmx, dmy, dmz, d2mx, d2my, d2mz, mi, mj, mk;
    double ifloat, jfloat, kfloat;
    double f, fnl, uix, uiy, uiz, uixnl, uiynl, uiznl;
    int i,j,k,nx, ny, nz, im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1;
    int kp1, kp2, ii, jj, kk;

    VASSERT(thee != VNULL);   /* PMG object with PBE info. */
    VASSERT(induced != VNULL); /* potential due to induced dipoles. */
    VASSERT(nlinduced != VNULL); /* potential due to non-local induced dipoles. */
    atom = Valist_getAtom(thee->pbe->alist, atomID);
    VASSERT(atom->partID != 0);    /* all atoms must be in the same partition. */
    apos = Vatom_getPosition(atom);
    dipole = Vatom_getInducedDipole(atom);
    uix = dipole[0];
    uiy = dipole[1];
    uiz = dipole[2];
    dipole = Vatom_getNLInducedDipole(atom);
    uixnl = dipole[0];
    uiynl = dipole[1];
    uiznl = dipole[2];
    u = induced->data;
    unl = nlinduced->data;

    for (i=0;i<3;i++){
       for (j=0;j<3;j++){
          de[i][j] = 0.0;
          denl[i][j] = 0.0;
       }
    }

    /* Mesh info */
    nx = induced->nx;
    ny = induced->ny;
    nz = induced->nz;
    hx = induced->hx;
    hy = induced->hy;
    hzed = induced->hzed;
    xmin = induced->xmin;
    ymin = induced->ymin;
    zmin = induced->zmin;
    xmax = induced->xmax;
    ymax = induced->ymax;
    zmax = induced->zmax;
    xlen = xmax-xmin;
    ylen = ymax-ymin;
    zlen = zmax-zmin;

    /* If we aren't in the current position, then we're done */
    if (atom->partID == 0) return;

    /* Make sure we're on the grid */
    if ((apos[0]<=(xmin+2*hx))   || (apos[0]>=(xmax-2*hx)) \
     || (apos[1]<=(ymin+2*hy))   || (apos[1]>=(ymax-2*hy)) \
     || (apos[2]<=(zmin+2*hzed)) || (apos[2]>=(zmax-2*hzed))) {
        Vnm_print(2, "qfMutualPolForce:  Atom off the mesh (ignoring) %6.3f %6.3f %6.3f\n", apos[0], apos[1], apos[2]);
        fflush(stderr);
    } else {

        /* Convert the atom position to grid coordinates */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;
        ifloat = position[0]/hx;
        jfloat = position[1]/hy;
        kfloat = position[2]/hzed;
        ip1 = (int)ceil(ifloat);
        ip2 = ip1 + 2;
        im1 = (int)floor(ifloat);
        im2 = im1 - 2;
        jp1 = (int)ceil(jfloat);
        jp2 = jp1 + 2;
        jm1 = (int)floor(jfloat);
        jm2 = jm1 - 2;
        kp1 = (int)ceil(kfloat);
        kp2 = kp1 + 2;
        km1 = (int)floor(kfloat);
        km2 = km1 - 2;

        /* This step shouldn't be necessary, but it saves nasty debugging
         * later on if something goes wrong */
        ip2 = VMIN2(ip2,nx-1);
        ip1 = VMIN2(ip1,nx-1);
        im1 = VMAX2(im1,0);
        im2 = VMAX2(im2,0);
        jp2 = VMIN2(jp2,ny-1);
        jp1 = VMIN2(jp1,ny-1);
        jm1 = VMAX2(jm1,0);
        jm2 = VMAX2(jm2,0);
        kp2 = VMIN2(kp2,nz-1);
        kp1 = VMIN2(kp1,nz-1);
        km1 = VMAX2(km1,0);
        km2 = VMAX2(km2,0);

        for (ii=im2; ii<=ip2; ii++) {
            mi = VFCHI4(ii,ifloat);
            mx = bspline4(mi);
            dmx = dbspline4(mi);
            d2mx = d2bspline4(mi);
            for (jj=jm2; jj<=jp2; jj++) {
                mj = VFCHI4(jj,jfloat);
                my = bspline4(mj);
                dmy = dbspline4(mj);
                d2my = d2bspline4(mj);
                for (kk=km2; kk<=kp2; kk++) {
                    mk = VFCHI4(kk,kfloat);
                    mz = bspline4(mk);
                    dmz = dbspline4(mk);
                    d2mz = d2bspline4(mk);
                    f = u[IJK(ii,jj,kk)];
                    fnl = unl[IJK(ii,jj,kk)];

                    /* The gradient of the reaction field
                       due to induced dipoles */
                    de[0][0] += f*d2mx*my*mz/(hx*hx);
                    de[1][0] += f*dmx*dmy*mz/(hy*hx);
                    de[1][1] += f*mx*d2my*mz/(hy*hy);
                    de[2][0] += f*dmx*my*dmz/(hx*hzed);
                    de[2][1] += f*mx*dmy*dmz/(hy*hzed);
                    de[2][2] += f*mx*my*d2mz/(hzed*hzed);

                    /* The gradient of the reaction field
                       due to non-local induced dipoles */
                    denl[0][0] += fnl*d2mx*my*mz/(hx*hx);
                    denl[1][0] += fnl*dmx*dmy*mz/(hy*hx);
                    denl[1][1] += fnl*mx*d2my*mz/(hy*hy);
                    denl[2][0] += fnl*dmx*my*dmz/(hx*hzed);
                    denl[2][1] += fnl*mx*dmy*dmz/(hy*hzed);
                    denl[2][2] += fnl*mx*my*d2mz/(hzed*hzed);
                }
            }
        }
    }

    /* mutual polarization force */
    force[0] = -(de[0][0]*uixnl + de[1][0]*uiynl + de[2][0]*uiznl);
    force[1] = -(de[1][0]*uixnl + de[1][1]*uiynl + de[2][1]*uiznl);
    force[2] = -(de[2][0]*uixnl + de[2][1]*uiynl + de[2][2]*uiznl);
    force[0] -=  denl[0][0]*uix + denl[1][0]*uiy + denl[2][0]*uiz;
    force[1] -=  denl[1][0]*uix + denl[1][1]*uiy + denl[2][1]*uiz;
    force[2] -=  denl[2][0]*uix + denl[2][1]*uiy + denl[2][2]*uiz;

    force[0] = 0.5 * force[0];
    force[1] = 0.5 * force[1];
    force[2] = 0.5 * force[2];

}

VPUBLIC void Vpmg_ibMutualPolForce(Vpmg *thee, Vgrid *induced, Vgrid *nlinduced,
                                   int atomID, double force[3]) {

    Vatom *atom;
    Valist *alist;
    Vacc *acc;
    Vpbe *pbe;
    Vsurf_Meth srfm;

    double *apos, position[3], arad, irad, zkappa2, hx, hy, hzed;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2;
    double rtot, dx, dx2, dy, dy2, dz, dz2, gpos[3], tgrad[3], fmag;
    double izmagic;
    int i, j, k, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);        /* We need a PMG object with PBE info. */
    VASSERT(induced != VNULL);     /* We need the potential due to induced dipoles. */
    VASSERT(nlinduced != VNULL);   /* We need the potential due to non-local induced dipoles. */
    VASSERT (!thee->pmgp->nonlin); /* Nonlinear PBE is not implemented for AMOEBA */

    atom = Valist_getAtom(thee->pbe->alist, atomID);
    VASSERT (atom->partID != 0);   /* Currently all atoms must be in the same partition. */

    acc = thee->pbe->acc;
    srfm = thee->surfMeth;
    apos = Vatom_getPosition(atom);
    arad = Vatom_getRadius(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* If we aren't in the current position, then we're done */
    if (atom->partID == 0) return;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    zkappa2 = Vpbe_getZkappa2(pbe);
    izmagic = 1.0/Vpbe_getZmagic(pbe);

    VASSERT (zkappa2 > VPMGSMALL); /* Should be a check for this further up.*/

    /* Mesh info */
    nx = induced->nx;
    ny = induced->ny;
    nz = induced->nz;
    hx = induced->hx;
    hy = induced->hy;
    hzed = induced->hzed;
    xmin = induced->xmin;
    ymin = induced->ymin;
    zmin = induced->zmin;
    xmax = induced->xmax;
    ymax = induced->ymax;
    zmax = induced->zmax;
    xlen = xmax-xmin;
    ylen = ymax-ymin;
    zlen = zmax-zmin;

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
        Vnm_print(2, "Vpmg_ibMutalPolForce:  Atom at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n", apos[0], apos[1], apos[2]);
        Vnm_print(2, "Vpmg_ibMutalPolForce:    xmin = %g, xmax = %g\n", xmin, xmax);
        Vnm_print(2, "Vpmg_ibMutalPolForce:    ymin = %g, ymax = %g\n", ymin, ymax);
        Vnm_print(2, "Vpmg_ibMutalPolForce:    zmin = %g, zmax = %g\n", zmin, zmax);
        fflush(stderr);
    } else {

        /* Convert the atom position to grid reference frame */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;

        /* Integrate over points within this atom's (inflated) radius */
        rtot = (irad + arad + thee->splineWin);
        rtot2 = VSQR(rtot);
        dx = rtot + 0.5*hx;
        imin = VMAX2(0,(int)ceil((position[0] - dx)/hx));
        imax = VMIN2(nx-1,(int)floor((position[0] + dx)/hx));
        for (i=imin; i<=imax; i++) {
            dx2 = VSQR(position[0] - hx*i);
            if (rtot2 > dx2) dy = VSQRT(rtot2 - dx2) + 0.5*hy;
            else dy = 0.5*hy;
            jmin = VMAX2(0,(int)ceil((position[1] - dy)/hy));
            jmax = VMIN2(ny-1,(int)floor((position[1] + dy)/hy));
            for (j=jmin; j<=jmax; j++) {
                dy2 = VSQR(position[1] - hy*j);
                if (rtot2 > (dx2+dy2)) dz = VSQRT(rtot2-dx2-dy2)+0.5*hzed;
                else dz = 0.5*hzed;
                kmin = VMAX2(0,(int)ceil((position[2] - dz)/hzed));
                kmax = VMIN2(nz-1,(int)floor((position[2] + dz)/hzed));
                for (k=kmin; k<=kmax; k++) {
                    dz2 = VSQR(k*hzed - position[2]);
                    /* See if grid point is inside ivdw radius and set ccf
                     * accordingly (do spline assignment here) */
                    if ((dz2 + dy2 + dx2) <= rtot2) {
                        gpos[0] = i*hx + xmin;
                        gpos[1] = j*hy + ymin;
                        gpos[2] = k*hzed + zmin;
                        Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, irad,
                          atom, tgrad);
                        fmag = induced->data[IJK(i,j,k)];
                        fmag *= nlinduced->data[IJK(i,j,k)];
                        fmag *= thee->kappa[IJK(i,j,k)];
                        force[0] += (zkappa2*fmag*tgrad[0]);
                        force[1] += (zkappa2*fmag*tgrad[1]);
                        force[2] += (zkappa2*fmag*tgrad[2]);
                    }
                } /* k loop */
            } /* j loop */
        } /* i loop */
    }

    force[0] = force[0] * 0.5 * hx * hy * hzed * izmagic;
    force[1] = force[1] * 0.5 * hx * hy * hzed * izmagic;
    force[2] = force[2] * 0.5 * hx * hy * hzed * izmagic;
}

VPUBLIC void Vpmg_dbMutualPolForce(Vpmg *thee, Vgrid *induced,
                                   Vgrid *nlinduced, int atomID,
                                   double force[3]) {

    Vatom *atom;
    Vacc *acc;
    Vpbe *pbe;
    Vsurf_Meth srfm;

    double *apos, position[3], arad, hx, hy, hzed, izmagic, deps, depsi;
    double xlen, ylen, zlen, xmin, ymin, zmin, xmax, ymax, zmax, rtot2, epsp;
    double rtot, dx, gpos[3], tgrad[3], dbFmag, epsw, kT;
    double *u, *unl, Hxijk, Hyijk, Hzijk, Hxim1jk, Hyijm1k, Hzijkm1;
    double dHxijk[3], dHyijk[3], dHzijk[3], dHxim1jk[3], dHyijm1k[3];
    double dHzijkm1[3];
    int i, j, k, l, nx, ny, nz, imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL); /* PMG object with PBE info. */
    VASSERT(induced != VNULL); /* potential due to induced dipoles.*/
    VASSERT(nlinduced != VNULL); /* potential due to non-local induced dipoles.*/

    acc = thee->pbe->acc;
    srfm = thee->surfMeth;
    atom = Valist_getAtom(thee->pbe->alist, atomID);
    VASSERT (atom->partID != 0); /* all atoms must be in the same partition.*/
    apos = Vatom_getPosition(atom);
    arad = Vatom_getRadius(atom);

    /* Reset force */
    force[0] = 0.0;
    force[1] = 0.0;
    force[2] = 0.0;

    /* Get PBE info */
    pbe = thee->pbe;
    acc = pbe->acc;
    epsp = Vpbe_getSoluteDiel(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    kT = Vpbe_getTemperature(pbe)*(1e-3)*Vunit_Na*Vunit_kb;
    izmagic = 1.0/Vpbe_getZmagic(pbe);

    deps = (epsw - epsp);
    depsi = 1.0/deps;
    VASSERT(VABS(deps) > VPMGSMALL);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;
    xmin = thee->pmgp->xmin;
    ymin = thee->pmgp->ymin;
    zmin = thee->pmgp->zmin;
    xmax = thee->pmgp->xmax;
    ymax = thee->pmgp->ymax;
    zmax = thee->pmgp->zmax;
    u = induced->data;
    unl = nlinduced->data;

    /* Make sure we're on the grid */
    if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
      (apos[1]<=ymin) || (apos[1]>=ymax)  || \
      (apos[2]<=zmin) || (apos[2]>=zmax)) {
        Vnm_print(2, "Vpmg_dbMutualPolForce:  Atom at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring):\n", apos[0], apos[1], apos[2]);
        Vnm_print(2, "Vpmg_dbMutualPolForce:    xmin = %g, xmax = %g\n", xmin, xmax);
        Vnm_print(2, "Vpmg_dbMutualPolForce:    ymin = %g, ymax = %g\n", ymin, ymax);
        Vnm_print(2, "Vpmg_dbMutualPolForce:    zmin = %g, zmax = %g\n", zmin, zmax);
        fflush(stderr);
    } else {

        /* Convert the atom position to grid reference frame */
        position[0] = apos[0] - xmin;
        position[1] = apos[1] - ymin;
        position[2] = apos[2] - zmin;

        /* Integrate over points within this atom's (inflated) radius */
        rtot = (arad + thee->splineWin);
        rtot2 = VSQR(rtot);
        dx = rtot/hx;
        imin = (int)floor((position[0]-rtot)/hx);
        if (imin < 1) {
            Vnm_print(2, "Vpmg_dbMutualPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        imax = (int)ceil((position[0]+rtot)/hx);
        if (imax > (nx-2)) {
            Vnm_print(2, "Vpmg_dbMutualPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        jmin = (int)floor((position[1]-rtot)/hy);
        if (jmin < 1) {
            Vnm_print(2, "Vpmg_dbMutualPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        jmax = (int)ceil((position[1]+rtot)/hy);
        if (jmax > (ny-2)) {
            Vnm_print(2, "Vpmg_dbMutualPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        kmin = (int)floor((position[2]-rtot)/hzed);
        if (kmin < 1) {
            Vnm_print(2, "Vpmg_dbMutualPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        kmax = (int)ceil((position[2]+rtot)/hzed);
        if (kmax > (nz-2)) {
            Vnm_print(2, "Vpmg_dbMutualPolForce:  Atom %d off grid!\n", atomID);
            return;
        }
        for (i=imin; i<=imax; i++) {
            for (j=jmin; j<=jmax; j++) {
                for (k=kmin; k<=kmax; k++) {
                    /* i,j,k */
                    gpos[0] = (i+0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxijk = (thee->epsx[IJK(i,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHxijk);
                    for (l=0; l<3; l++) dHxijk[l] *= Hxijk;
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j+0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijk = (thee->epsy[IJK(i,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHyijk);
                    for (l=0; l<3; l++) dHyijk[l] *= Hyijk;
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k+0.5)*hzed + zmin;
                    Hzijk = (thee->epsz[IJK(i,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHzijk);
                    for (l=0; l<3; l++) dHzijk[l] *= Hzijk;
                    /* i-1,j,k */
                    gpos[0] = (i-0.5)*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hxim1jk = (thee->epsx[IJK(i-1,j,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHxim1jk);
                    for (l=0; l<3; l++) dHxim1jk[l] *= Hxim1jk;
                    /* i,j-1,k */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = (j-0.5)*hy + ymin;
                    gpos[2] = k*hzed + zmin;
                    Hyijm1k = (thee->epsy[IJK(i,j-1,k)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHyijm1k);
                    for (l=0; l<3; l++) dHyijm1k[l] *= Hyijm1k;
                    /* i,j,k-1 */
                    gpos[0] = i*hx + xmin;
                    gpos[1] = j*hy + ymin;
                    gpos[2] = (k-0.5)*hzed + zmin;
                    Hzijkm1 = (thee->epsz[IJK(i,j,k-1)] - epsp)*depsi;
                    Vpmg_splineSelect(srfm, acc, gpos, thee->splineWin, 0.,
                            atom, dHzijkm1);
                    for (l=0; l<3; l++) dHzijkm1[l] *= Hzijkm1;
                    dbFmag = unl[IJK(i,j,k)];
                    tgrad[0] =
                       (dHxijk[0]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[0]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[0]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[0]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[0]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[0]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    tgrad[1] =
                       (dHxijk[1]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[1]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[1]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[1]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[1]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[1]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                    tgrad[2] =
                       (dHxijk[2]  *(u[IJK(i+1,j,k)]-u[IJK(i,j,k)])
                     +  dHxim1jk[2]*(u[IJK(i-1,j,k)]-u[IJK(i,j,k)]))/VSQR(hx)
                     + (dHyijk[2]  *(u[IJK(i,j+1,k)]-u[IJK(i,j,k)])
                     +  dHyijm1k[2]*(u[IJK(i,j-1,k)]-u[IJK(i,j,k)]))/VSQR(hy)
                     + (dHzijk[2]  *(u[IJK(i,j,k+1)]-u[IJK(i,j,k)])
                     + dHzijkm1[2]*(u[IJK(i,j,k-1)]-u[IJK(i,j,k)]))/VSQR(hzed);
                     force[0] += (dbFmag*tgrad[0]);
                     force[1] += (dbFmag*tgrad[1]);
                     force[2] += (dbFmag*tgrad[2]);
                } /* k loop */
            } /* j loop */
        } /* i loop */

        force[0] = -force[0]*hx*hy*hzed*deps*0.5*izmagic;
        force[1] = -force[1]*hx*hy*hzed*deps*0.5*izmagic;
        force[2] = -force[2]*hx*hy*hzed*deps*0.5*izmagic;
    }
}

#endif /* if defined(WITH_TINKER) */

VPRIVATE void fillcoCoefSpline4(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax, ionmask, ionstr, dist2;
    double xlen, ylen, zlen, position[3], itot, stot, ictot, ictot2, sctot;
    double irad, dx, dy, dz, epsw, epsp, w2i;
    double hx, hy, hzed, *apos, arad, sctot2;
    double dx2, dy2, dz2, stot2, itot2, rtot, rtot2, splineWin;
    double dist, value, denom, sm, sm2, sm3, sm4, sm5, sm6, sm7;
    double e, e2, e3, e4, e5, e6, e7;
    double b, b2, b3, b4, b5, b6, b7;
    double c0, c1, c2, c3, c4, c5, c6, c7;
    double ic0, ic1, ic2, ic3, ic4, ic5, ic6, ic7;
    int i, j, k, nx, ny, nz, iatom;
    int imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);
    splineWin = thee->splineWin;

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    ionstr = Vpbe_getBulkIonicStrength(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    epsp = Vpbe_getSoluteDiel(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* This is a floating point parameter related to the non-zero nature of the
     * bulk ionic strength.  If the ionic strength is greater than zero; this
     * parameter is set to 1.0 and later scaled by the appropriate pre-factors.
     * Otherwise, this parameter is set to 0.0 */
    if (ionstr > VPMGSMALL) ionmask = 1.0;
    else ionmask = 0.0;

    /* Reset the kappa, epsx, epsy, and epsz arrays */
    for (i=0; i<(nx*ny*nz); i++) {
        thee->kappa[i] = 1.0;
        thee->epsx[i] = 1.0;
        thee->epsy[i] = 1.0;
        thee->epsz[i] = 1.0;
    }

    /* Loop through the atoms and do assign the dielectric */
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        arad = Vatom_getRadius(atom);

        b = arad - splineWin;
        e = arad + splineWin;
        e2 = e * e;
        e3 = e2 * e;
        e4 = e3 * e;
        e5 = e4 * e;
        e6 = e5 * e;
        e7 = e6 * e;
        b2 = b * b;
        b3 = b2 * b;
        b4 = b3 * b;
        b5 = b4 * b;
        b6 = b5 * b;
        b7 = b6 * b;
        denom = e7  - 7.0*b*e6 + 21.0*b2*e5 - 35.0*e4*b3
              + 35.0*e3*b4 - 21.0*b5*e2  + 7.0*e*b6 - b7;
        c0 = b4*(35.0*e3 - 21.0*b*e2 + 7*e*b2 - b3)/denom;
        c1 = -140.0*b3*e3/denom;
        c2 = 210.0*e2*b2*(e + b)/denom;
        c3 = -140.0*e*b*(e2 + 3.0*b*e + b2)/denom;
        c4 =  35.0*(e3 + 9.0*b*e2 + + 9.0*e*b2 + b3)/denom;
        c5 = -84.0*(e2 + 3.0*b*e + b2)/denom;
        c6 =  70.0*(e + b)/denom;
        c7 = -20.0/denom;

        b = irad + arad - splineWin;
        e = irad + arad + splineWin;
        e2 = e * e;
        e3 = e2 * e;
        e4 = e3 * e;
        e5 = e4 * e;
        e6 = e5 * e;
        e7 = e6 * e;
        b2 = b * b;
        b3 = b2 * b;
        b4 = b3 * b;
        b5 = b4 * b;
        b6 = b5 * b;
        b7 = b6 * b;
        denom = e7  - 7.0*b*e6 + 21.0*b2*e5 - 35.0*e4*b3
              + 35.0*e3*b4 - 21.0*b5*e2  + 7.0*e*b6 - b7;
        ic0 = b4*(35.0*e3 - 21.0*b*e2 + 7*e*b2 - b3)/denom;
        ic1 = -140.0*b3*e3/denom;
        ic2 = 210.0*e2*b2*(e + b)/denom;
        ic3 = -140.0*e*b*(e2 + 3.0*b*e + b2)/denom;
        ic4 =  35.0*(e3 + 9.0*b*e2 + + 9.0*e*b2 + b3)/denom;
        ic5 = -84.0*(e2 + 3.0*b*e + b2)/denom;
        ic6 =  70.0*(e + b)/denom;
        ic7 = -20.0/denom;

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
                (thee->pmgp->bcfl != BCFL_MAP)) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f,\
 %4.3f) is off the mesh (ignoring):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:    xmin = %g, xmax = %g\n",
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:    ymin = %g, ymax = %g\n",
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:    zmin = %g, zmax = %g\n",
                  zmin, zmax);
            }
            fflush(stderr);

        } else if (arad > VPMGSMALL ) { /* if we're on the mesh */

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* MARK ION ACCESSIBILITY AND DIELECTRIC VALUES FOR LATER
             * ASSIGNMENT (Steps #1-3) */
            itot = irad + arad + splineWin;
            itot2 = VSQR(itot);
            ictot = VMAX2(0, (irad + arad - splineWin));
            ictot2 = VSQR(ictot);
            stot = arad + splineWin;
            stot2 = VSQR(stot);
            sctot = VMAX2(0, (arad - splineWin));
            sctot2 = VSQR(sctot);

           /* We'll search over grid points which are in the greater of
             * these two radii */
            rtot = VMAX2(itot, stot);
            rtot2 = VMAX2(itot2, stot2);
            dx = rtot + 0.5*hx;
            dy = rtot + 0.5*hy;
            dz = rtot + 0.5*hzed;
            imin = VMAX2(0,(int)floor((position[0] - dx)/hx));
            imax = VMIN2(nx-1,(int)ceil((position[0] + dx)/hx));
            jmin = VMAX2(0,(int)floor((position[1] - dy)/hy));
            jmax = VMIN2(ny-1,(int)ceil((position[1] + dy)/hy));
            kmin = VMAX2(0,(int)floor((position[2] - dz)/hzed));
            kmax = VMIN2(nz-1,(int)ceil((position[2] + dz)/hzed));
            for (i=imin; i<=imax; i++) {
                dx2 = VSQR(position[0] - hx*i);
                for (j=jmin; j<=jmax; j++) {
                    dy2 = VSQR(position[1] - hy*j);
                    for (k=kmin; k<=kmax; k++) {
                        dz2 = VSQR(position[2] - k*hzed);

                        /* ASSIGN CCF */
                        if (thee->kappa[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2 + dy2 + dx2;
                            if (dist2 >= itot2) {
                                ;
                            }
                            if (dist2 <= ictot2) {
                                thee->kappa[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 < itot2) && (dist2 > ictot2)) {
                                dist = VSQRT(dist2);
                                sm = dist;
                                sm2 = dist2;
                                sm3 = sm2 * sm;
                                sm4 = sm3 * sm;
                                sm5 = sm4 * sm;
                                sm6 = sm5 * sm;
                                sm7 = sm6 * sm;
                                value = ic0 + ic1*sm + ic2*sm2 + ic3*sm3
                                      + ic4*sm4 + ic5*sm5 + ic6*sm6 + ic7*sm7;
                                if (value > 1.0) {
                                   value = 1.0;
                                } else if (value < 0.0){
                                   value = 0.0;
                                }
                                thee->kappa[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A1CF */
                        if (thee->epsx[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2+dy2+VSQR(position[0]-(i+0.5)*hx);
                            if (dist2 >= stot2) {
                                thee->epsx[IJK(i,j,k)] *= 1.0;
                            }
                            if (dist2 <= sctot2) {
                                thee->epsx[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist;
                                sm2 = VSQR(sm);
                                sm3 = sm2 * sm;
                                sm4 = sm3 * sm;
                                sm5 = sm4 * sm;
                                sm6 = sm5 * sm;
                                sm7 = sm6 * sm;
                                value = c0 + c1*sm + c2*sm2 + c3*sm3
                                      + c4*sm4 + c5*sm5 + c6*sm6 + c7*sm7;
                                if (value > 1.0) {
                                   value = 1.0;
                                } else if (value < 0.0){
                                   value = 0.0;
                                }
                                thee->epsx[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A2CF */
                        if (thee->epsy[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2+dx2+VSQR(position[1]-(j+0.5)*hy);
                            if (dist2 >= stot2) {
                                thee->epsy[IJK(i,j,k)] *= 1.0;
                            }
                            if (dist2 <= sctot2) {
                                thee->epsy[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist;
                                sm2 = VSQR(sm);
                                sm3 = sm2 * sm;
                                sm4 = sm3 * sm;
                                sm5 = sm4 * sm;
                                sm6 = sm5 * sm;
                                sm7 = sm6 * sm;
                                value = c0 + c1*sm + c2*sm2 + c3*sm3
                                      + c4*sm4 + c5*sm5 + c6*sm6 + c7*sm7;
                                if (value > 1.0) {
                                   value = 1.0;
                                } else if (value < 0.0){
                                   value = 0.0;
                                }
                                thee->epsy[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A3CF */
                        if (thee->epsz[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dy2+dx2+VSQR(position[2]-(k+0.5)*hzed);
                            if (dist2 >= stot2) {
                                thee->epsz[IJK(i,j,k)] *= 1.0;
                            }
                            if (dist2 <= sctot2) {
                                thee->epsz[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist;
                                sm2 = dist2;
                                sm3 = sm2 * sm;
                                sm4 = sm3 * sm;
                                sm5 = sm4 * sm;
                                sm6 = sm5 * sm;
                                sm7 = sm6 * sm;
                                value = c0 + c1*sm + c2*sm2 + c3*sm3
                                      + c4*sm4 + c5*sm5 + c6*sm6 + c7*sm7;
                                if (value > 1.0) {
                                   value = 1.0;
                                } else if (value < 0.0){
                                   value = 0.0;
                                }
                                thee->epsz[IJK(i,j,k)] *= value;
                            }
                        }


                    } /* k loop */
                } /* j loop */
            } /* i loop */
        } /* endif (on the mesh) */
    } /* endfor (over all atoms) */

    Vnm_print(0, "Vpmg_fillco:  filling coefficient arrays\n");
    /* Interpret markings and fill the coefficient arrays */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {

                thee->kappa[IJK(i,j,k)] = ionmask*thee->kappa[IJK(i,j,k)];
                thee->epsx[IJK(i,j,k)] = (epsw-epsp)*thee->epsx[IJK(i,j,k)]
                  + epsp;
                thee->epsy[IJK(i,j,k)] = (epsw-epsp)*thee->epsy[IJK(i,j,k)]
                  + epsp;
                thee->epsz[IJK(i,j,k)] = (epsw-epsp)*thee->epsz[IJK(i,j,k)]
                  + epsp;

            } /* i loop */
        } /* j loop */
    } /* k loop */

}

VPUBLIC void fillcoPermanentInduced(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    /* Coversions */
    double zmagic, f;
    /* Grid */
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double xlen, ylen, zlen, position[3], ifloat, jfloat, kfloat;
    double hx, hy, hzed, *apos;
    /* Multipole */
    double charge, *dipole,*quad;
    double c,ux,uy,uz,qxx,qyx,qyy,qzx,qzy,qzz,qave;
    /* B-spline weights */
    double mx,my,mz,dmx,dmy,dmz,d2mx,d2my,d2mz;
    double mi,mj,mk;
    /* Loop variables */
    int i, ii, jj, kk, nx, ny, nz, iatom;
    int im2, im1, ip1, ip2, jm2, jm1, jp1, jp2, km2, km1, kp1, kp2;

    VASSERT(thee != VNULL);

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    zmagic = Vpbe_getZmagic(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Conversion */
    f = zmagic/(hx*hy*hzed);

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* Fill in the source term (permanent atomic multipoles
       and induced dipoles) */
    Vnm_print(0, "fillcoPermanentInduced:  filling in source term.\n");
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);

        c = Vatom_getCharge(atom)*f;

#if defined(WITH_TINKER)
        dipole = Vatom_getDipole(atom);
        ux = dipole[0]/hx*f;
        uy = dipole[1]/hy*f;
        uz = dipole[2]/hzed*f;
        dipole = Vatom_getInducedDipole(atom);
        ux = ux + dipole[0]/hx*f;
        uy = uy + dipole[1]/hy*f;
        uz = uz + dipole[2]/hzed*f;
        quad = Vatom_getQuadrupole(atom);
        qxx = (1.0/3.0)*quad[0]/(hx*hx)*f;
        qyx = (2.0/3.0)*quad[3]/(hx*hy)*f;
        qyy = (1.0/3.0)*quad[4]/(hy*hy)*f;
        qzx = (2.0/3.0)*quad[6]/(hzed*hx)*f;
        qzy = (2.0/3.0)*quad[7]/(hzed*hy)*f;
        qzz = (1.0/3.0)*quad[8]/(hzed*hzed)*f;
#else
        ux = 0.0;
        uy = 0.0;
        uz = 0.0;
        qxx = 0.0;
        qyx = 0.0;
        qyy = 0.0;
        qzx = 0.0;
        qzy = 0.0;
        qzz = 0.0;
#endif /* if defined(WITH_TINKER) */

        /* Make sure we're on the grid */
        if ((apos[0]<=(xmin-2*hx)) || (apos[0]>=(xmax+2*hx))  || \
            (apos[1]<=(ymin-2*hy)) || (apos[1]>=(ymax+2*hy))  || \
            (apos[2]<=(zmin-2*hzed)) || (apos[2]>=(zmax+2*hzed))) {
            Vnm_print(2, "fillcoPermanentMultipole: Atom #%d at (%4.3f, %4.3f, %4.3f) is off the mesh (ignoring this atom):\n", iatom, apos[0], apos[1], apos[2]);
            Vnm_print(2, "fillcoPermanentMultipole: xmin = %g, xmax = %g\n", xmin, xmax);
            Vnm_print(2, "fillcoPermanentMultipole: ymin = %g, ymax = %g\n", ymin, ymax);
            Vnm_print(2, "fillcoPermanentMultipole: zmin = %g, zmax = %g\n", zmin, zmax);
            fflush(stderr);
        } else {

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* Figure out which vertices we're next to */
            ifloat = position[0]/hx;
            jfloat = position[1]/hy;
            kfloat = position[2]/hzed;

            ip1   = (int)ceil(ifloat);
            ip2   = ip1 + 2;
            im1   = (int)floor(ifloat);
            im2   = im1 - 2;
            jp1   = (int)ceil(jfloat);
            jp2   = jp1 + 2;
            jm1   = (int)floor(jfloat);
            jm2   = jm1 - 2;
            kp1   = (int)ceil(kfloat);
            kp2   = kp1 + 2;
            km1   = (int)floor(kfloat);
            km2   = km1 - 2;

            /* This step shouldn't be necessary, but it saves nasty debugging
             * later on if something goes wrong */
            ip2 = VMIN2(ip2,nx-1);
            ip1 = VMIN2(ip1,nx-1);
            im1 = VMAX2(im1,0);
            im2 = VMAX2(im2,0);
            jp2 = VMIN2(jp2,ny-1);
            jp1 = VMIN2(jp1,ny-1);
            jm1 = VMAX2(jm1,0);
            jm2 = VMAX2(jm2,0);
            kp2 = VMIN2(kp2,nz-1);
            kp1 = VMIN2(kp1,nz-1);
            km1 = VMAX2(km1,0);
            km2 = VMAX2(km2,0);

            /* Now assign fractions of the charge to the nearby verts */
            for (ii=im2; ii<=ip2; ii++) {
                mi = VFCHI4(ii,ifloat);
                mx = bspline4(mi);
                dmx = dbspline4(mi);
                d2mx = d2bspline4(mi);
                for (jj=jm2; jj<=jp2; jj++) {
                    mj = VFCHI4(jj,jfloat);
                    my = bspline4(mj);
                    dmy = dbspline4(mj);
                    d2my = d2bspline4(mj);
                    for (kk=km2; kk<=kp2; kk++) {
                        mk = VFCHI4(kk,kfloat);
                        mz = bspline4(mk);
                        dmz = dbspline4(mk);
                        d2mz = d2bspline4(mk);
                        charge = mx*my*mz*c -
                         dmx*my*mz*ux - mx*dmy*mz*uy - mx*my*dmz*uz +
                         d2mx*my*mz*qxx +
                         dmx*dmy*mz*qyx + mx*d2my*mz*qyy +
                         dmx*my*dmz*qzx + mx*dmy*dmz*qzy + mx*my*d2mz*qzz;
                        thee->charge[IJK(ii,jj,kk)] += charge;

                    }
                }
            }
        } /* endif (on the mesh) */

    } /* endfor (each atom) */
}

VPRIVATE void fillcoCoefSpline3(Vpmg *thee) {

    Valist *alist;
    Vpbe *pbe;
    Vatom *atom;
    double xmin, xmax, ymin, ymax, zmin, zmax, ionmask, ionstr, dist2;
    double xlen, ylen, zlen, position[3], itot, stot, ictot, ictot2, sctot;
    double irad, dx, dy, dz, epsw, epsp, w2i;
    double hx, hy, hzed, *apos, arad, sctot2;
    double dx2, dy2, dz2, stot2, itot2, rtot, rtot2, splineWin;
    double dist, value, denom, sm, sm2, sm3, sm4, sm5;
    double e, e2, e3, e4, e5;
    double b, b2, b3, b4, b5;
    double c0, c1, c2, c3, c4, c5;
    double ic0, ic1, ic2, ic3, ic4, ic5;
    int i, j, k, nx, ny, nz, iatom;
    int imin, imax, jmin, jmax, kmin, kmax;

    VASSERT(thee != VNULL);
    splineWin = thee->splineWin;

    /* Get PBE info */
    pbe = thee->pbe;
    alist = pbe->alist;
    irad = Vpbe_getMaxIonRadius(pbe);
    ionstr = Vpbe_getBulkIonicStrength(pbe);
    epsw = Vpbe_getSolventDiel(pbe);
    epsp = Vpbe_getSoluteDiel(pbe);

    /* Mesh info */
    nx = thee->pmgp->nx;
    ny = thee->pmgp->ny;
    nz = thee->pmgp->nz;
    hx = thee->pmgp->hx;
    hy = thee->pmgp->hy;
    hzed = thee->pmgp->hzed;

    /* Define the total domain size */
    xlen = thee->pmgp->xlen;
    ylen = thee->pmgp->ylen;
    zlen = thee->pmgp->zlen;

    /* Define the min/max dimensions */
    xmin = thee->pmgp->xcent - (xlen/2.0);
    ymin = thee->pmgp->ycent - (ylen/2.0);
    zmin = thee->pmgp->zcent - (zlen/2.0);
    xmax = thee->pmgp->xcent + (xlen/2.0);
    ymax = thee->pmgp->ycent + (ylen/2.0);
    zmax = thee->pmgp->zcent + (zlen/2.0);

    /* This is a floating point parameter related to the non-zero nature of the
     * bulk ionic strength.  If the ionic strength is greater than zero; this
     * parameter is set to 1.0 and later scaled by the appropriate pre-factors.
     * Otherwise, this parameter is set to 0.0 */
    if (ionstr > VPMGSMALL) ionmask = 1.0;
    else ionmask = 0.0;

    /* Reset the kappa, epsx, epsy, and epsz arrays */
    for (i=0; i<(nx*ny*nz); i++) {
        thee->kappa[i] = 1.0;
        thee->epsx[i] = 1.0;
        thee->epsy[i] = 1.0;
        thee->epsz[i] = 1.0;
    }

    /* Loop through the atoms and do assign the dielectric */
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {

        atom = Valist_getAtom(alist, iatom);
        apos = Vatom_getPosition(atom);
        arad = Vatom_getRadius(atom);

        b = arad - splineWin;
        e = arad + splineWin;
        e2 = e * e;
        e3 = e2 * e;
        e4 = e3 * e;
        e5 = e4 * e;
        b2 = b * b;
        b3 = b2 * b;
        b4 = b3 * b;
        b5 = b4 * b;
        denom = pow((e - b), 5.0);
        c0 = -10.0*e2*b3 + 5.0*e*b4 - b5;
        c1 = 30.0*e2*b2;
        c2 = -30.0*(e2*b + e*b2);
        c3 = 10.0*(e2 + 4.0*e*b + b2);
        c4 = -15.0*(e + b);
        c5 = 6;
        c0 = c0/denom;
        c1 = c1/denom;
        c2 = c2/denom;
        c3 = c3/denom;
        c4 = c4/denom;
        c5 = c5/denom;

        b = irad + arad - splineWin;
        e = irad + arad + splineWin;
        e2 = e * e;
        e3 = e2 * e;
        e4 = e3 * e;
        e5 = e4 * e;
        b2 = b * b;
        b3 = b2 * b;
        b4 = b3 * b;
        b5 = b4 * b;
        denom = pow((e - b), 5.0);
        ic0 = -10.0*e2*b3 + 5.0*e*b4 - b5;
        ic1 = 30.0*e2*b2;
        ic2 = -30.0*(e2*b + e*b2);
        ic3 = 10.0*(e2 + 4.0*e*b + b2);
        ic4 = -15.0*(e + b);
        ic5 = 6;
        ic0 = c0/denom;
        ic1 = c1/denom;
        ic2 = c2/denom;
        ic3 = c3/denom;
        ic4 = c4/denom;
        ic5 = c5/denom;

        /* Make sure we're on the grid */
        if ((apos[0]<=xmin) || (apos[0]>=xmax)  || \
            (apos[1]<=ymin) || (apos[1]>=ymax)  || \
            (apos[2]<=zmin) || (apos[2]>=zmax)) {
            if ((thee->pmgp->bcfl != BCFL_FOCUS) &&
                (thee->pmgp->bcfl != BCFL_MAP)) {
                Vnm_print(2, "Vpmg_fillco:  Atom #%d at (%4.3f, %4.3f,\
 %4.3f) is off the mesh (ignoring):\n",
                  iatom, apos[0], apos[1], apos[2]);
                Vnm_print(2, "Vpmg_fillco:    xmin = %g, xmax = %g\n",
                  xmin, xmax);
                Vnm_print(2, "Vpmg_fillco:    ymin = %g, ymax = %g\n",
                  ymin, ymax);
                Vnm_print(2, "Vpmg_fillco:    zmin = %g, zmax = %g\n",
                  zmin, zmax);
            }
            fflush(stderr);

        } else if (arad > VPMGSMALL ) { /* if we're on the mesh */

            /* Convert the atom position to grid reference frame */
            position[0] = apos[0] - xmin;
            position[1] = apos[1] - ymin;
            position[2] = apos[2] - zmin;

            /* MARK ION ACCESSIBILITY AND DIELECTRIC VALUES FOR LATER
             * ASSIGNMENT (Steps #1-3) */
            itot = irad + arad + splineWin;
            itot2 = VSQR(itot);
            ictot = VMAX2(0, (irad + arad - splineWin));
            ictot2 = VSQR(ictot);
            stot = arad + splineWin;
            stot2 = VSQR(stot);
            sctot = VMAX2(0, (arad - splineWin));
            sctot2 = VSQR(sctot);

           /* We'll search over grid points which are in the greater of
             * these two radii */
            rtot = VMAX2(itot, stot);
            rtot2 = VMAX2(itot2, stot2);
            dx = rtot + 0.5*hx;
            dy = rtot + 0.5*hy;
            dz = rtot + 0.5*hzed;
            imin = VMAX2(0,(int)floor((position[0] - dx)/hx));
            imax = VMIN2(nx-1,(int)ceil((position[0] + dx)/hx));
            jmin = VMAX2(0,(int)floor((position[1] - dy)/hy));
            jmax = VMIN2(ny-1,(int)ceil((position[1] + dy)/hy));
            kmin = VMAX2(0,(int)floor((position[2] - dz)/hzed));
            kmax = VMIN2(nz-1,(int)ceil((position[2] + dz)/hzed));
            for (i=imin; i<=imax; i++) {
                dx2 = VSQR(position[0] - hx*i);
                for (j=jmin; j<=jmax; j++) {
                    dy2 = VSQR(position[1] - hy*j);
                    for (k=kmin; k<=kmax; k++) {
                        dz2 = VSQR(position[2] - k*hzed);

                        /* ASSIGN CCF */
                        if (thee->kappa[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2 + dy2 + dx2;
                            if (dist2 >= itot2) {
                                ;
                            }
                            if (dist2 <= ictot2) {
                                thee->kappa[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 < itot2) && (dist2 > ictot2)) {
                                dist = VSQRT(dist2);
                                sm = dist;
                                sm2 = dist2;
                                sm3 = sm2 * sm;
                                sm4 = sm3 * sm;
                                sm5 = sm4 * sm;
                                value = ic0 + ic1*sm + ic2*sm2 + ic3*sm3
                                      + ic4*sm4 + ic5*sm5;
                                if (value > 1.0) {
                                   value = 1.0;
                                } else if (value < 0.0){
                                   value = 0.0;
                                }
                                thee->kappa[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A1CF */
                        if (thee->epsx[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2+dy2+VSQR(position[0]-(i+0.5)*hx);
                            if (dist2 >= stot2) {
                                thee->epsx[IJK(i,j,k)] *= 1.0;
                            }
                            if (dist2 <= sctot2) {
                                thee->epsx[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist;
                                sm2 = VSQR(sm);
                                sm3 = sm2 * sm;
                                sm4 = sm3 * sm;
                                sm5 = sm4 * sm;
                                value = c0 + c1*sm + c2*sm2 + c3*sm3
                                      + c4*sm4 + c5*sm5;
                                if (value > 1.0) {
                                   value = 1.0;
                                } else if (value < 0.0){
                                   value = 0.0;
                                }
                                thee->epsx[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A2CF */
                        if (thee->epsy[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dz2+dx2+VSQR(position[1]-(j+0.5)*hy);
                            if (dist2 >= stot2) {
                                thee->epsy[IJK(i,j,k)] *= 1.0;
                            }
                            if (dist2 <= sctot2) {
                                thee->epsy[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist;
                                sm2 = VSQR(sm);
                                sm3 = sm2 * sm;
                                sm4 = sm3 * sm;
                                sm5 = sm4 * sm;
                                value = c0 + c1*sm + c2*sm2 + c3*sm3
                                      + c4*sm4 + c5*sm5;
                                if (value > 1.0) {
                                   value = 1.0;
                                } else if (value < 0.0){
                                   value = 0.0;
                                }
                                thee->epsy[IJK(i,j,k)] *= value;
                            }
                        }

                        /* ASSIGN A3CF */
                        if (thee->epsz[IJK(i,j,k)] > VPMGSMALL) {
                            dist2 = dy2+dx2+VSQR(position[2]-(k+0.5)*hzed);
                            if (dist2 >= stot2) {
                                thee->epsz[IJK(i,j,k)] *= 1.0;
                            }
                            if (dist2 <= sctot2) {
                                thee->epsz[IJK(i,j,k)] = 0.0;
                            }
                            if ((dist2 > sctot2) && (dist2 < stot2)) {
                                dist = VSQRT(dist2);
                                sm = dist;
                                sm2 = dist2;
                                sm3 = sm2 * sm;
                                sm4 = sm3 * sm;
                                sm5 = sm4 * sm;
                                value = c0 + c1*sm + c2*sm2 + c3*sm3
                                      + c4*sm4 + c5*sm5;
                                if (value > 1.0) {
                                   value = 1.0;
                                } else if (value < 0.0){
                                   value = 0.0;
                                }
                                thee->epsz[IJK(i,j,k)] *= value;
                            }
                        }


                    } /* k loop */
                } /* j loop */
            } /* i loop */
        } /* endif (on the mesh) */
    } /* endfor (over all atoms) */

    Vnm_print(0, "Vpmg_fillco:  filling coefficient arrays\n");
    /* Interpret markings and fill the coefficient arrays */
    for (k=0; k<nz; k++) {
        for (j=0; j<ny; j++) {
            for (i=0; i<nx; i++) {

                thee->kappa[IJK(i,j,k)] = ionmask*thee->kappa[IJK(i,j,k)];
                thee->epsx[IJK(i,j,k)] = (epsw-epsp)*thee->epsx[IJK(i,j,k)]
                  + epsp;
                thee->epsy[IJK(i,j,k)] = (epsw-epsp)*thee->epsy[IJK(i,j,k)]
                  + epsp;
                thee->epsz[IJK(i,j,k)] = (epsw-epsp)*thee->epsz[IJK(i,j,k)]
                  + epsp;

            } /* i loop */
        } /* j loop */
    } /* k loop */

}

VPRIVATE void bcolcomp(int *iparm, double *rparm, int *iwork, double *rwork,
        double *values, int *rowind, int *colptr, int *flag) {
    int nrow, ncol, nnzero, i;
    int nxc, nyc, nzc, nf, nc, narr, narrc, n_rpc;
    int n_iz, n_ipc, iretot, iintot;
    int nrwk, niwk, nx, ny, nz, nlev, ierror, maxlev, mxlv;
    int mgcoar, mgdisc, mgsolv;
    int k_iz;
    int k_ipc, k_rpc, k_ac, k_cc, k_fc, k_pc;

    WARN_UNTESTED;

    // Decode some parameters
    nrwk = VAT(iparm, 1);
    niwk = VAT(iparm, 2);
    nx   = VAT(iparm, 3);
    ny   = VAT(iparm, 4);
    nz   = VAT(iparm, 5);
    nlev = VAT(iparm, 6);

    // Some checks on input
    mxlv = Vmaxlev(nx, ny, nz);

    // Basic grid sizes, etc.
    mgcoar = VAT(iparm, 18);
    mgdisc = VAT(iparm, 19);
    mgsolv = VAT(iparm, 21);
    Vmgsz(&mgcoar, &mgdisc, &mgsolv,
            &nx, &ny, &nz,
            &nlev,
            &nxc, &nyc, &nzc,
            &nf, &nc,
            &narr, &narrc,
            &n_rpc, &n_iz, &n_ipc,
            &iretot, &iintot);

    // Split up the integer work array
    k_iz  = 1;
    k_ipc = k_iz + n_iz;

    // Split up the real work array
    k_rpc = 1;
    k_cc  = k_rpc + n_rpc;
    k_fc  = k_cc  + narr;
    k_pc  = k_fc  + narr;
    k_ac  = k_pc  + 27*narrc;

    bcolcomp2(iparm, rparm,
            &nx, &ny, &nz, RAT(iwork, k_iz),
            RAT(iwork, k_ipc), RAT(rwork, k_rpc),
            RAT(rwork, k_ac), RAT(rwork, k_cc),
            values, rowind, colptr, flag);
}

VPRIVATE void bcolcomp2(int *iparm, double *rparm,
        int *nx, int *ny, int *nz,
        int *iz, int *ipc, double *rpc,
        double *ac, double *cc, double *values,
        int *rowind, int *colptr, int *flag) {

    int nlev = 1;
    int lev = VAT(iparm, 6);

    MAT2(iz, 50, nlev);

    WARN_UNTESTED;

    /*
     * Build the multigrid data structure in iz
     *    THIS MAY HAVE BEEN DONE ALREADY, BUT IT'S OK TO DO IT AGAIN,
     *    RIGHT?
     *    call buildstr (nx,ny,nz,nlev,iz)
     *
     *    We're interested in the finest level
     */
    bcolcomp3(nx, ny, nz,
            RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
            RAT(ac, VAT2(iz, 7, lev)), RAT(cc, VAT2(iz, 1, lev)),
            values, rowind, colptr, flag);
}

/**************************************************************************
 * Routine:  bcolcomp3
 * Purpose:  Build a column-compressed matrix in Harwell-Boeing format
 * Args:     flag   0 ==> Use Poisson operator only
 *                  1 ==> Use linearization of full operator around current
 *                        solution
 * Author:   Nathan Baker (mostly ripped off from Harwell-Boeing format
 *           documentation)
 **************************************************************************/
VPRIVATE void bcolcomp3(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc,
        double *values, int *rowind, int *colptr, int *flag) {

    MAT2(ac, *nx * *ny * *nz, 1);

    WARN_UNTESTED;

    bcolcomp4(nx, ny, nz,
            ipc, rpc,
            RAT2(ac, 1, 1), cc,
            RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4),
            values, rowind, colptr, flag);
}



/**************************************************************************
 * Routine:  bcolcomp4
 * Purpose:  Build a column-compressed matrix in Harwell-Boeing format
 * Args:     flag   0 ==> Use Poisson operator only
 *                  1 ==> Use linearization of full operator around current
 *                        solution
 * Author:   Nathan Baker (mostly ripped off from Harwell-Boeing format
 *           documentation)
 **************************************************************************/
VPRIVATE void bcolcomp4(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *oC, double *cc, double *oE, double *oN, double *uC,
        double *values, int *rowind, int *colptr, int *flag) {

    int nxm2, nym2, nzm2;
    int ii, jj, kk, ll;
    int  i,  j,  k,  l;
    int inonz, iirow, nn, nrow, ncol, nonz, irow, n;

    int doit;

    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(uC, *nx, *ny, *nz);
    MAT3(cc, *nx, *ny, *nz);
    MAT3(oC, *nx, *ny, *nz);

    WARN_UNTESTED;

    // Get some column, row, and nonzero information
    n = *nx * *ny * *nz;
    nxm2 = *nx - 2;
    nym2 = *ny - 2;
    nzm2 = *nz - 2;
    nn   = nxm2 * nym2 * nzm2;
    ncol = nn;
    nrow = nn;
    nonz = 7 * nn - 2 * nxm2 * nym2 - 2 * nxm2 - 2;

    // Intialize some pointers
    inonz = 1;

    /*
     * Run over the dimensions of the matrix (non-zero only in the interior
     * of the mesh
     */
    for (k=2; k<=*nz-1; k++) {
        // Offset the index to the output grid index
        kk = k - 1;

        for (j=2; j<=*ny-1; j++) {
            // Offset the index to the output grid index
            jj = j - 1;

            for (i=2; i<=*nx-1; i++) {
                // Offset the index to the output grid index
                ii = i - 1;

                // Get the output (i,j,k) row number in natural ordering
                ll = (kk - 1) * nxm2 * nym2 + (jj - 1) * nxm2 + (ii - 1) + 1;
                l  = (k  - 1) * *nx  * *ny  + (j  - 1) * *nx  + (i  - 1) + 1;

                // Store where this column starts
                VAT(colptr,ll) = inonz;

                // SUB-DIAGONAL 3
                iirow = ll - nxm2 * nym2;
                irow  =  l - *nx  * *ny;

                doit = (iirow >= 1) && (iirow <= nn);
                doit = doit && (irow >= 1) && (irow <= n);

                if (doit) {
                    VAT(values, inonz) = -VAT3(uC, i, j, k-1);
                    VAT(rowind, inonz) = iirow;
                    inonz++;
                }



                // SUB-DIAGONAL 2
                iirow = ll - nxm2;
                irow  =  l - *nx;

                doit = (iirow >= 1) && (iirow <= nn);
                doit = doit && (irow >= 1) && (irow <= n);

                if (doit) {
                    VAT(values, inonz) = -VAT3(oN, i, j-1, k);
                    VAT(rowind, inonz) = iirow;
                    inonz++;
                }



                // SUB-DIAGONAL 1
                iirow = ll - 1;
                irow =   l - 1;

                doit = (iirow >= 1) && (iirow <= nn);
                doit = doit && (irow <= 1) && (irow <= n);
                if (doit) {
                    VAT(values, inonz) = -VAT3(oE, i-1, j, k);
                    VAT(rowind, inonz) = iirow;
                    inonz++;
                }



                // DIAGONAL
                iirow = ll;
                irow  =  l;

                if (*flag == 0) {
                    VAT(values, inonz) = VAT3(oC, i, j, k);
                } else if (*flag == 1) {
                    VAT(values, inonz) = VAT3(oC, i, j, k)
                                       + VAT3(cc, i, j, k);
                } else {
                    VABORT_MSG0("PMGF1");
                }

                VAT(rowind, inonz) = iirow;
                inonz++;

                // SUPER-DIAGONAL 1
                iirow = ll + 1;
                irow  =  l + 1;
                doit = (iirow >= 1) && (iirow <= nn);
                doit = doit && (irow >= 1) && (irow <= n);
                if (doit) {
                    VAT(values, inonz) = -VAT3(oE, i, j, k);
                    VAT(rowind, inonz) = iirow;
                    inonz++;
                }



                // SUPER-DIAGONAL 2
                iirow = ll + nxm2;
                irow  =  l + *nx;
                doit = (iirow >= 1) && (iirow <= nn);
                doit = doit && (irow >= 1) && (irow <= n);
                if (doit) {
                    VAT(values, inonz) = -VAT3(oN, i, j, k);
                    VAT(rowind, inonz) = iirow;
                    inonz++;
                }



                // SUPER-DIAGONAL 3
                iirow = ll + nxm2 * nym2;
                irow  =  l + *nx  * *ny;
                doit = (iirow >= 1) && (iirow <= nn);
                doit = doit && (irow >= 1) && (irow <= n);
                if (doit) {
                    VAT(values, inonz) = -VAT3(uC, i, j, k);
                    VAT(rowind, inonz) = iirow;
                    inonz++;
                }
            }
        }
    }

    VAT(colptr, ncol + 1) = inonz;

    if (inonz != (nonz + 1)) {
        VABORT_MSG2("BCOLCOMP4:  ERROR -- INONZ = %d, NONZ = %d", inonz, nonz);
    }
}



VPRIVATE void pcolcomp(int *nrow, int *ncol, int *nnzero,
        double *values, int *rowind, int *colptr,
        char *path, char *title, char *mxtype) {

    char key[] = "key";
    char ptrfmt[] = "(10I8)";
    char indfmt[] = "(10I8)";
    char valfmt[] = "(5E15.8)";
    char rhsfmt[] = "(5E15.8)";

    int i, totcrd, ptrcrd, indcrd, valcrd, neltvl, rhscrd;

    FILE *outFile;

    WARN_UNTESTED;

    // Open the file for reading
    outFile = fopen(path, "w");

    // Set some default values
    ptrcrd = (int)(*ncol   / 10 + 1) - 1;
    indcrd = (int)(*nnzero / 10 + 1) - 1;
    valcrd = (int)(*nnzero / 10 + 1) - 1;
    totcrd = ptrcrd + indcrd + valcrd;
    rhscrd = 0;
    neltvl = 0;

    // Print the header
    fprintf(outFile, "%72s%8s\n",
            title, key);
    fprintf(outFile, "%14d%14d%14d%14d%14d\n",
            totcrd, ptrcrd, indcrd, valcrd, rhscrd);
    fprintf(outFile, "%3s\n", mxtype);
    fprintf(outFile, "           %14d%14d%14d%14d\n",
            *nrow, *ncol, *nnzero, neltvl);
    fprintf(outFile, "%16s%16s%20s%20s\n",
            ptrfmt, indfmt, valfmt, rhsfmt);

    // Write the matrix structure
    for (i=1; i<=*ncol+1; i++)
        fprintf(outFile, "%8d", VAT(colptr, i));
    fprintf(outFile, "\n");

    for (i=1; i<=*nnzero; i++)
        fprintf(outFile, "%8d", VAT(rowind, i));
    fprintf(outFile, "\n");

    // Write out the values
    if (valcrd > 0) {
        for (i=1; i<=*nnzero; i++)
            fprintf(outFile, "%15.8e", VAT(values, i));
        fprintf(outFile, "\n");
    }

    // Close the file
    fclose (outFile);
}
