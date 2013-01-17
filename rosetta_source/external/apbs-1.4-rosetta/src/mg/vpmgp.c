/**
 *  @file    vpmgp.c
 *  @author  Nathan Baker
 *  @brief   Class Vpmgp methods
 *  @ingroup Vpmgp
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
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
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

#include "vpmgp.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VACC)
#endif /* if !defined(VINLINE_VACC) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_ctor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpmgp* Vpmgp_ctor(MGparm *mgparm) {

    Vpmgp *thee = VNULL;

    /* Set up the structure */
    thee = (Vpmgp*)Vmem_malloc(VNULL, 1, sizeof(Vpmgp) );
    VASSERT( thee != VNULL);
    VASSERT(Vpmgp_ctor2(thee,mgparm));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_ctor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpmgp_ctor2(Vpmgp *thee,MGparm *mgparm) {

    /* Specified parameters */
    thee->nx = mgparm->dime[0];
    thee->ny = mgparm->dime[1];
    thee->nz = mgparm->dime[2];
    thee->hx = mgparm->grid[0];
    thee->hy = mgparm->grid[1];
    thee->hzed = mgparm->grid[2];
    thee->xlen = ((double)(mgparm->dime[0]-1))*mgparm->grid[0];
    thee->ylen = ((double)(mgparm->dime[1]-1))*mgparm->grid[1];
    thee->zlen = ((double)(mgparm->dime[2]-1))*mgparm->grid[2];
    thee->nlev = mgparm->nlev;

    thee->nonlin = mgparm->nonlintype;
    thee->meth = mgparm->method;

#ifdef DEBUG_MAC_OSX_OCL
#include "mach_chud.h"
    if(kOpenCLAvailable)
        thee->meth = 4;
#endif

    if (thee->nonlin == NONLIN_LPBE) thee->ipkey = IPKEY_LPBE; /* LPBE case */
    else if(thee->nonlin == NONLIN_SMPBE) thee->ipkey = IPKEY_SMPBE; /* SMPBE case */
    else thee->ipkey = IPKEY_NPBE; /* NPBE standard case */

    /* Default parameters */
    if (mgparm->setetol) { /* If etol is set by the user in APBS input file, then use this custom-defined etol */
        thee->errtol = mgparm->etol;
        Vnm_print(1, "  Error tolerance (etol) is now set to user-defined \
value: %g \n", thee->errtol);
        Vnm_print(0, "Error tolerance (etol) is now set to user-defined \
value: %g \n", thee->errtol);
    } else thee->errtol = 1.0e-6;   /* Here are a few comments.  Mike had this set to
        * 1e-9; convential wisdom sets this at 1e-6 for
        * the PBE; Ray Luo sets this at 1e-3 for his
        * accelerated PBE (for dynamics, etc.) */
    thee->itmax = 200;
    thee->istop = 1;
    thee->iinfo = 1;         /* I'd recommend either 1 (for debugging LPBE) or 2 (for debugging NPBE), higher values give too much output */

    thee->bcfl = BCFL_SDH;
    thee->key = 0;
    thee->iperf = 0;
    thee->mgcoar = 2;
    thee->mgkey = 0;
    thee->nu1 = 2;
    thee->nu2 = 2;
    thee->mgprol = 0;
    thee->mgdisc = 0;
    thee->omegal = 19.4e-1;
    thee->omegan = 9.0e-1;
    thee->ipcon = 3;
    thee->irite = 8;
    thee->xcent = 0.0;
    thee->ycent = 0.0;
    thee->zcent = 0.0;

    /* Default value for all APBS runs */
    thee->mgsmoo = 1;
    if (thee->nonlin == NONLIN_NPBE || thee->nonlin == NONLIN_SMPBE) {
        /* SMPBE Added - SMPBE needs to mimic NPBE */
        Vnm_print(0, "Vpmp_ctor2:  Using meth = 1, mgsolv = 0\n");
        thee->mgsolv = 0;
    } else {
        /* Most rigorous (good for testing) */
        Vnm_print(0, "Vpmp_ctor2:  Using meth = 2, mgsolv = 1\n");
        thee->mgsolv = 1;
    }

    /* TEMPORARY USEAQUA */
    /* If we are using aqua, our solution method is either VSOL_CGMGAqua or VSOL_NewtonAqua
     * so we need to temporarily override the mgsolve value and set it to 0
     */
    if(mgparm->useAqua == 1) thee->mgsolv = 0;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmgp_dtor(Vpmgp **thee) {

    if ((*thee) != VNULL) {
        Vpmgp_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vpmgp), (void **)thee);
        (*thee) = VNULL;
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpmgp_dtor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpmgp_dtor2(Vpmgp *thee) { ; }


VPUBLIC void Vpmgp_size(
    Vpmgp *thee
    )
{

    int num_nf = 0;
    int num_narr = 2;
    int num_narrc = 27;
    int nxf, nyf, nzf, level, num_nf_oper, num_narrc_oper, n_band, nc_band, num_band, iretot;

    thee->nf = thee->nx * thee->ny * thee->nz;
    thee->narr = thee->nf;
    nxf = thee->nx;
    nyf = thee->ny;
    nzf = thee->nz;
    thee->nxc = thee->nx;
    thee->nyc = thee->ny;
    thee->nzc = thee->nz;

    for (level=2; level<=thee->nlev; level++) {
        Vpmgp_makeCoarse(1, nxf, nyf, nzf, &(thee->nxc), &(thee->nyc), &(thee->nzc)); /* NAB TO-DO -- implement this function and check which variables need to be passed by reference... */
        nxf = thee->nxc;
        nyf = thee->nyc;
        nzf = thee->nzc;
        thee->narr = thee->narr + (nxf * nyf * nzf);
    }

    thee->nc = thee->nxc * thee->nyc * thee->nzc;
    thee->narrc = thee->narr - thee->nf;

    /* Box or FEM discretization on fine grid? */
    switch (thee->mgdisc) { /* NAB TO-DO:  This needs to be changed into an enumeration */
    case 0:
        num_nf_oper = 4;
        break;
    case 1:
        num_nf_oper = 14;
        break;
    default:
        Vnm_print(2, "Vpmgp_size:  Invalid mgdisc value (%d)!\n", thee->mgdisc);
        VASSERT(0);
    }

    /* Galerkin or standard coarsening? */
    switch (thee->mgcoar) { /* NAB TO-DO:  This needs to be changed into an enumeration */
    case 0:
        if (thee->mgdisc != 0) {
            Vnm_print(2, "Vpmgp_size:  Invalid mgcoar value (%d); must be used with mgdisc 0!\n", thee->mgcoar);
            VASSERT(0);
        }
        num_narrc_oper = 4;
        break;
    case 1:
        if (thee->mgdisc != 0) {
            Vnm_print(2, "Vpmgp_size:  Invalid mgcoar value (%d); must be used with mgdisc 0!\n", thee->mgcoar);
            VASSERT(0);
        }
        num_narrc_oper = 14;
        break;
    case 2:
        num_narrc_oper = 14;
        break;
    default:
        Vnm_print(2, "Vpmgp_size:  Invalid mgcoar value (%d)!\n", thee->mgcoar);
        VASSERT(0);
    }

    /* LINPACK storage on coarse grid */
    switch (thee->mgsolv) { /* NAB TO-DO:  This needs to be changed into an enumeration */
    case 0:
        n_band = 0;
        break;
    case 1:
        if ( ( (thee->mgcoar == 0) || (thee->mgcoar == 1)) && (thee->mgdisc == 0) ) {
            num_band = 1 + (thee->nxc-2)*(thee->nyc-2);
        } else {
            num_band = 1 + (thee->nxc-2)*(thee->nyc-2) + (thee->nxc-2) + 1;
        }
        nc_band = (thee->nxc-2)*(thee->nyc-2)*(thee->nzc-2);
        n_band  = nc_band * num_band;
        break;
    default:
        Vnm_print(2, "Vpmgp_size:  Invalid mgsolv value (%d)!\n", thee->mgsolv);
        VASSERT(0);
    }

    /* Real storage parameters */
    thee->n_rpc = 100*(thee->nlev+1);

    /* Resulting total required for real storage */
    thee->nrwk = num_narr*thee->narr + (num_nf + num_nf_oper)*thee->nf + (num_narrc + num_narrc_oper)*thee->narrc + n_band + thee->n_rpc;

    /* Integer storage parameters */
    thee->n_iz = 50*(thee->nlev+1);
    thee->n_ipc = 100*(thee->nlev+1);
    thee->niwk = thee->n_iz + thee->n_ipc;
}

VPRIVATE int coarsenThis(int nOld) {

    int nOut;

    nOut = (nOld - 1) / 2 + 1;

    if (((nOut-1)*2) != (nOld-1)) {
        Vnm_print(2, "Vpmgp_makeCoarse:  Warning!  The grid dimensions you have chosen are not consistent with the nlev you have specified!\n");
        Vnm_print(2, "Vpmgp_makeCoarse:  This calculation will only work if you are running with mg-dummy type.\n");
    }
    if (nOut < 1) {
        Vnm_print(2, "D'oh!  You coarsened the grid below zero!  How did you do that?\n");
        VASSERT(0);
    }

    return nOut;
}

VPUBLIC void Vpmgp_makeCoarse(
    int numLevel,
    int nxOld,
    int nyOld,
    int nzOld,
    int *nxNew,
    int *nyNew,
    int *nzNew
    )
{
    int nxtmp, nytmp, nztmp, iLevel;

    for (iLevel=0; iLevel<numLevel; iLevel++) {
        nxtmp = *nxNew;
        nytmp = *nyNew;
        nztmp = *nzNew;
        *nxNew = coarsenThis(nxtmp);
        *nyNew = coarsenThis(nytmp);
        *nzNew = coarsenThis(nztmp);
    }


}
