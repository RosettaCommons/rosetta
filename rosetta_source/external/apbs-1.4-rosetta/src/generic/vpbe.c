/**
 *  @file    vpbe.c
 *  @ingroup Vpbe
 *  @author  Nathan Baker
 *  @brief   Class Vpbe methods
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

#include "vpbe.h"

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Private method declaration
/////////////////////////////////////////////////////////////////////////// */
#define MAX_SPLINE_WINDOW 0.5

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Inlineable methods
/////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VPBE)

VPUBLIC Valist* Vpbe_getValist(Vpbe *thee) {

   VASSERT(thee != VNULL);
   return thee->alist;

}

VPUBLIC Vacc* Vpbe_getVacc(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->acc;

}

VPUBLIC double Vpbe_getBulkIonicStrength(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->bulkIonicStrength;
}

VPUBLIC double Vpbe_getTemperature(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->T;

}

VPUBLIC double Vpbe_getSoluteDiel(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->soluteDiel;

}

VPUBLIC double* Vpbe_getSoluteCenter(Vpbe *thee) {

   VASSERT(thee != VNULL);
   return thee->soluteCenter;
}

VPUBLIC double Vpbe_getSolventDiel(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->solventDiel;
}

VPUBLIC double Vpbe_getSolventRadius(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->solventRadius;
}

VPUBLIC double Vpbe_getMaxIonRadius(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->maxIonRadius;
}

VPUBLIC double Vpbe_getXkappa(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->xkappa;
}

VPUBLIC double Vpbe_getDeblen(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->deblen;
}

VPUBLIC double Vpbe_getZkappa2(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->zkappa2;
}

VPUBLIC double Vpbe_getZmagic(Vpbe *thee) {

   VASSERT(thee != VNULL);
   VASSERT(thee->paramFlag);
   return thee->zmagic;
}

VPUBLIC double Vpbe_getSoluteRadius(Vpbe *thee) {

   VASSERT(thee != VNULL);
   return thee->soluteRadius;
}

VPUBLIC double Vpbe_getSoluteXlen(Vpbe *thee) {

   VASSERT(thee != VNULL);
   return thee->soluteXlen;
}

VPUBLIC double Vpbe_getSoluteYlen(Vpbe *thee) {

   VASSERT(thee != VNULL);
   return thee->soluteYlen;
}

VPUBLIC double Vpbe_getSoluteZlen(Vpbe *thee) {

   VASSERT(thee != VNULL);
   return thee->soluteZlen;
}

VPUBLIC double Vpbe_getSoluteCharge(Vpbe *thee) {

   VASSERT(thee != VNULL);
   return thee->soluteCharge;
}

/* ///////////////////////////////////////////////////////////////////////////
 // Routine:  Vpbe_getzmem
 // Purpose: This routine returns values stored in the structure thee.
 // Author:  Michael Grabe
 /////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getzmem(Vpbe *thee) {

    VASSERT(thee != VNULL);
    VASSERT(thee->param2Flag);
    return thee->z_mem;
}

/* ///////////////////////////////////////////////////////////////////////////
 // Routine:  Vpbe_getLmem
 // Purpose: This routine returns values stored in the structure thee.
 // Author:  Michael Grabe
 /////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getLmem(Vpbe *thee) {

    VASSERT(thee != VNULL);
    VASSERT(thee->param2Flag);
    return thee->L;
}

/* ///////////////////////////////////////////////////////////////////////////
 // Routine:  Vpbe_getmembraneDiel
 // Purpose: This routine returns values stored in the structure thee.
 // Author:  Michael Grabe
 /////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getmembraneDiel(Vpbe *thee) {

    VASSERT(thee != VNULL);
    VASSERT(thee->param2Flag);
    return thee->membraneDiel;
}

/* ///////////////////////////////////////////////////////////////////////////
 // Routine:  Vpbe_getmemv
 // Purpose: This routine returns values stored in the structure thee.
 // Author:  Michael Grabe
 /////////////////////////////////////////////////////////////////////////// */
VPUBLIC double Vpbe_getmemv(Vpbe *thee) {

    VASSERT(thee != VNULL);
    VASSERT(thee->param2Flag);
    return thee->V;
}

#endif /* if !defined(VINLINE_VPBE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpbe: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

VPUBLIC Vpbe* Vpbe_ctor(Valist *alist, int ionNum, double *ionConc,
                        double *ionRadii, double *ionQ, double T,
                        double soluteDiel, double solventDiel,
                        double solventRadius, int focusFlag, double sdens,
                        double z_mem, double L, double membraneDiel, double V ) {

    /* Set up the structure */
    Vpbe *thee = VNULL;
    thee = (Vpbe*)Vmem_malloc(VNULL, 1, sizeof(Vpbe) );
    VASSERT( thee != VNULL);
    VASSERT( Vpbe_ctor2(thee, alist, ionNum, ionConc, ionRadii, ionQ,
                        T, soluteDiel, solventDiel, solventRadius, focusFlag, sdens,
                        z_mem, L, membraneDiel, V) );

    return thee;
}


VPUBLIC int Vpbe_ctor2(Vpbe *thee, Valist *alist, int ionNum,
                       double *ionConc, double *ionRadii,
                       double *ionQ, double T, double soluteDiel,
                       double solventDiel, double solventRadius, int focusFlag,
                       double sdens, double z_mem, double L, double membraneDiel,
                       double V) {

    int i, iatom, inhash[3];
    double atomRadius;
    Vatom *atom;
    double center[3] = {0.0, 0.0, 0.0};
    double lower_corner[3] = {0.0, 0.0, 0.0};
    double upper_corner[3] = {0.0, 0.0, 0.0};
    double disp[3], dist, radius, charge, xmin, xmax, ymin, ymax, zmin, zmax;
    double x, y, z, netCharge;
    double nhash[3];
    const double N_A = 6.022045000e+23;
    const double e_c = 4.803242384e-10;
    const double k_B = 1.380662000e-16;
    const double pi  = 4. * VATAN(1.);

    /* Set up memory management object */
    thee->vmem = Vmem_ctor("APBS::VPBE");

    VASSERT(thee != VNULL);
    if (alist == VNULL) {
        Vnm_print(2, "Vpbe_ctor2: Got null pointer to Valist object!\n");
        return 0;
    }

    /* **** STUFF THAT GETS DONE FOR EVERYONE **** */
    /* Set pointers */
    thee->alist = alist;
    thee->paramFlag = 0;

    /* Determine solute center */
    center[0] = thee->alist->center[0];
    center[1] = thee->alist->center[1];
    center[2] = thee->alist->center[2];
    thee->soluteCenter[0] = center[0];
    thee->soluteCenter[1] = center[1];
    thee->soluteCenter[2] = center[2];

    /* Determine solute length and charge*/
    radius = 0;
    atom = Valist_getAtom(thee->alist, 0);
    xmin = Vatom_getPosition(atom)[0];
    xmax = Vatom_getPosition(atom)[0];
    ymin = Vatom_getPosition(atom)[1];
    ymax = Vatom_getPosition(atom)[1];
    zmin = Vatom_getPosition(atom)[2];
    zmax = Vatom_getPosition(atom)[2];
    charge = 0;
    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) {
        atom = Valist_getAtom(thee->alist, iatom);
        atomRadius = Vatom_getRadius(atom);
        x = Vatom_getPosition(atom)[0];
        y = Vatom_getPosition(atom)[1];
        z = Vatom_getPosition(atom)[2];
        if ((x+atomRadius) > xmax) xmax = x + atomRadius;
        if ((x-atomRadius) < xmin) xmin = x - atomRadius;
        if ((y+atomRadius) > ymax) ymax = y + atomRadius;
        if ((y-atomRadius) < ymin) ymin = y - atomRadius;
        if ((z+atomRadius) > zmax) zmax = z + atomRadius;
        if ((z-atomRadius) < zmin) zmin = z - atomRadius;
        disp[0] = (x - center[0]);
        disp[1] = (y - center[1]);
        disp[2] = (z - center[2]);
        dist = (disp[0]*disp[0]) + (disp[1]*disp[1]) + (disp[2]*disp[2]);
        dist = VSQRT(dist) + atomRadius;
        if (dist > radius) radius = dist;
        charge += Vatom_getCharge(Valist_getAtom(thee->alist, iatom));
    }
    thee->soluteRadius = radius;
    Vnm_print(0, "Vpbe_ctor2:  solute radius = %g\n", radius);
    thee->soluteXlen = xmax - xmin;
    thee->soluteYlen = ymax - ymin;
    thee->soluteZlen = zmax - zmin;
    Vnm_print(0, "Vpbe_ctor2:  solute dimensions = %g x %g x %g\n",
            thee->soluteXlen, thee->soluteYlen, thee->soluteZlen);
    thee->soluteCharge = charge;
    Vnm_print(0, "Vpbe_ctor2:  solute charge = %g\n", charge);

    /* Set parameters */
    thee->numIon = ionNum;
    if (thee->numIon >= MAXION) {
        Vnm_print(2, "Vpbe_ctor2:  Too many ion species (MAX = %d)!\n",
          MAXION);
        return 0;
    }
    thee->bulkIonicStrength = 0.0;
    thee->maxIonRadius = 0.0;
    netCharge = 0.0;
    for (i=0; i<thee->numIon; i++) {
        thee->ionConc[i] = ionConc[i];
        thee->ionRadii[i] = ionRadii[i];
        if (ionRadii[i] > thee->maxIonRadius) thee->maxIonRadius = ionRadii[i];
        thee->ionQ[i] = ionQ[i];
        thee->bulkIonicStrength += (0.5*ionConc[i]*VSQR(ionQ[i]));
        netCharge += (ionConc[i]*ionQ[i]);
    }
#ifndef VAPBSQUIET
    Vnm_print(1, "  Vpbe_ctor:  Using max ion radius (%g A) for exclusion \
function\n", thee->maxIonRadius);
#endif
    if (VABS(netCharge) > VSMALL) {
        Vnm_print(2, "Vpbe_ctor2:  You have a counterion charge imbalance!\n");
        Vnm_print(2, "Vpbe_ctor2:  Net charge conc. = %g M\n", netCharge);
        return 0;
    }
    thee->T = T;
    thee->soluteDiel = soluteDiel;
    thee->solventDiel = solventDiel;
    thee->solventRadius = solventRadius;

    /* Compute parameters:
     *
     * kappa^2 = (8 pi N_A e_c^2) I_s / (1000 eps_w k_B T)
     * kappa   = 0.325567 * I_s^{1/2}   angstroms^{-1}
     * deblen  = 1 / kappa
     *         = 3.071564378 * I_s^{1/2}   angstroms
     * \bar{kappa}^2 = eps_w * kappa^2
     * zmagic  = (4 * pi * e_c^2) / (k_B T)   (we scale the diagonal later)
     *         = 7046.528838
     */
    if (thee->T == 0.0) {
        Vnm_print(2, "Vpbe_ctor2:  You set the temperature to 0 K.\n");
        Vnm_print(2, "Vpbe_ctor2:  That violates the 3rd Law of Thermo.!");
        return 0;
    }
    if (thee->bulkIonicStrength == 0.) {
        thee->xkappa  = 0.;
        thee->deblen  = 0.;
        thee->zkappa2 = 0.;
    } else {
        thee->xkappa  = VSQRT( thee->bulkIonicStrength * 1.0e-16 *
            ((8.0 * pi * N_A * e_c*e_c) /
            (1000.0 * thee->solventDiel * k_B * T))
        );
        thee->deblen  = 1. / thee->xkappa;
        thee->zkappa2 = thee->solventDiel * VSQR(thee->xkappa);
    }
    Vnm_print(0, "Vpbe_ctor2:  bulk ionic strength = %g\n",
            thee->bulkIonicStrength);
    Vnm_print(0, "Vpbe_ctor2:  xkappa = %g\n", thee->xkappa);
    Vnm_print(0, "Vpbe_ctor2:  Debye length = %g\n", thee->deblen);
    Vnm_print(0, "Vpbe_ctor2:  zkappa2 = %g\n", thee->zkappa2);
    thee->zmagic  = ((4.0 * pi * e_c*e_c) / (k_B * thee->T)) * 1.0e+8;
    Vnm_print(0, "Vpbe_ctor2:  zmagic = %g\n", thee->zmagic);

    /* Compute accessibility objects:
     *   - Allow for extra room in the case of spline windowing
     *   - Place some limits on the size of the hash table in the case of very
     *     large molecules
     */
    if (thee->maxIonRadius > thee->solventRadius)
        radius = thee->maxIonRadius + MAX_SPLINE_WINDOW;
    else radius = thee->solventRadius + MAX_SPLINE_WINDOW;

    nhash[0] = (thee->soluteXlen)/0.5;
    nhash[1] = (thee->soluteYlen)/0.5;
    nhash[2] = (thee->soluteZlen)/0.5;
    for (i=0; i<3; i++) inhash[i] = (int)(nhash[i]);

    for (i=0;i<3;i++){
        if (inhash[i] < 3) inhash[i] = 3;
        if (inhash[i] > MAX_HASH_DIM) inhash[i] = MAX_HASH_DIM;
    }
    Vnm_print(0, "Vpbe_ctor2:  Constructing Vclist with %d x %d x %d table\n",
            inhash[0], inhash[1], inhash[2]);

    thee->clist = Vclist_ctor(thee->alist, radius, inhash,
            CLIST_AUTO_DOMAIN, lower_corner, upper_corner);

    VASSERT(thee->clist != VNULL);
    thee->acc = Vacc_ctor(thee->alist, thee->clist, sdens);

    VASSERT(thee->acc != VNULL);

    /* SMPBE Added */
    thee->smsize = 0.0;
    thee->smvolume = 0.0;
    thee->ipkey = 0;

    thee->paramFlag = 1;

    /*-----------------------------------------------------------*/
    /* added by Michael Grabe                                    */
    /*-----------------------------------------------------------*/

    thee->z_mem = z_mem;
    thee->L = L;
    thee->membraneDiel = membraneDiel;
    thee->V = V;

    //    if (V != VNULL) thee->param2Flag = 1;
    //    else thee->param2Flag = 0;

    /*-----------------------------------------------------------*/

    return 1;
}

VPUBLIC void Vpbe_dtor(Vpbe **thee) {
    if ((*thee) != VNULL) {
        Vpbe_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vpbe), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void Vpbe_dtor2(Vpbe *thee) {
    Vclist_dtor(&(thee->clist));
    Vacc_dtor(&(thee->acc));
    Vmem_dtor(&(thee->vmem));
}

VPUBLIC double Vpbe_getCoulombEnergy1(Vpbe *thee) {

    int i, j, k, natoms;

    double dist, *ipos, *jpos, icharge, jcharge;
    double energy = 0.0;
    double eps, T;
    Vatom *iatom, *jatom;
    Valist *alist;

    VASSERT(thee != VNULL);
    alist = Vpbe_getValist(thee);
    VASSERT(alist != VNULL);
    natoms = Valist_getNumberAtoms(alist);

    /* Do the sum */
    for (i=0; i<natoms; i++) {
        iatom = Valist_getAtom(alist,i);
        icharge = Vatom_getCharge(iatom);
        ipos = Vatom_getPosition(iatom);
        for (j=i+1; j<natoms; j++) {
            jatom = Valist_getAtom(alist,j);
            jcharge = Vatom_getCharge(jatom);
            jpos = Vatom_getPosition(jatom);
            dist = 0;
            for (k=0; k<3; k++) dist += ((ipos[k]-jpos[k])*(ipos[k]-jpos[k]));
            dist = VSQRT(dist);
            energy = energy + icharge*jcharge/dist;
        }
    }

    /* Convert the result to J */
    T = Vpbe_getTemperature(thee);
    eps = Vpbe_getSoluteDiel(thee);
    energy = energy*Vunit_ec*Vunit_ec/(4*Vunit_pi*Vunit_eps0*eps*(1.0e-10));

    /* Scale by Boltzmann energy */
    energy = energy/(Vunit_kb*T);

    return energy;
}

VPUBLIC unsigned long int Vpbe_memChk(Vpbe *thee) {

    unsigned long int memUse = 0;

    if (thee == VNULL) return 0;

    memUse = memUse + sizeof(Vpbe);
    memUse = memUse + (unsigned long int)Vacc_memChk(thee->acc);

    return memUse;
}

VPUBLIC int Vpbe_getIons(Vpbe *thee, int *nion, double ionConc[MAXION],
  double ionRadii[MAXION], double ionQ[MAXION]) {

    int i;

    VASSERT(thee != VNULL);

    *nion = thee->numIon;
    for (i=0; i<(*nion); i++) {
        ionConc[i] = thee->ionConc[i];
        ionRadii[i] = thee->ionRadii[i];
        ionQ[i] = thee->ionQ[i];
    }

    return *nion;
}
