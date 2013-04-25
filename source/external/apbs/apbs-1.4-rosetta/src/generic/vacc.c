/**
 *  @file    vacc.c
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @brief   Class Vacc methods
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

#include "vacc.h"

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_VACC)

VPUBLIC unsigned long int Vacc_memChk(Vacc *thee) {
    if (thee == VNULL)
        return 0;
    return Vmem_bytes(thee->mem);
}

#endif /* if !defined(VINLINE_VACC) */

/**
 * @brief  Determines if a point is within the union of the spheres centered
 *         at the atomic centers with radii equal to the sum of their van der
 *         Waals radii and the probe radius.  Does not include contributions
 *         from the specified atom.
 * @returns 1 if accessible (outside the inflated van der Waals radius), 0
 *          otherwise
 * @author  Nathan Baker
 */
VPRIVATE int ivdwAccExclus(
                           Vacc *thee,  /** Accessibility object */
                           double center[3],  /** Position to test */
                           double radius,  /** Radius of probe */
                           int atomID  /** ID of atom to ignore */
                           ) {

    int iatom;
    double dist2,
           *apos;
    Vatom *atom;
    VclistCell *cell;

    VASSERT(thee != VNULL);

    /* We can only test probes with radii less than the max specified */
    if (radius > Vclist_maxRadius(thee->clist)) {
        Vnm_print(2,
                  "Vacc_ivdwAcc: got radius (%g) bigger than max radius (%g)\n",
                  radius, Vclist_maxRadius(thee->clist));
        VASSERT(0);
    }

    /* Get the relevant cell from the cell list */
    cell = Vclist_getCell(thee->clist, center);

    /* If we have no cell, then no atoms are nearby and we're definitely
     * accessible */
    if (cell == VNULL) {
        return 1;
    }

    /* Otherwise, check for overlap with the atoms in the cell */
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom = cell->atoms[iatom];

        // We don't actually need to test this if the atom IDs do match; don't compute this if we're comparing atom against itself.
        if (atom->id == atomID) continue;

        apos = atom->position;
        dist2 = VSQR(center[0]-apos[0]) + VSQR(center[1]-apos[1])
                        + VSQR(center[2]-apos[2]);
        if (dist2 < VSQR(atom->radius+radius)){
            return 0;
        }
    }

    /* If we're still here, then the point is accessible */
    return 1;

}

VPUBLIC Vacc* Vacc_ctor(Valist *alist,
                        Vclist *clist,
                        double surf_density /* Surface density */
                        ) {


    Vacc *thee = VNULL;

    /* Set up the structure */
    thee = (Vacc*)Vmem_malloc(VNULL, 1, sizeof(Vacc) );
    VASSERT( thee != VNULL);
    VASSERT( Vacc_ctor2(thee, alist, clist, surf_density));
    return thee;
}

/** Check and store parameters passed to constructor */
VPRIVATE int Vacc_storeParms(Vacc *thee,
                             Valist *alist,
                             Vclist *clist,
                             double surf_density /* Surface density */
                             ) {

    int nsphere,
        iatom;
    double maxrad = 0.0,
           maxarea,
           rad;
    Vatom *atom;

    if (alist == VNULL) {
        Vnm_print(2, "Vacc_storeParms:  Got NULL Valist!\n");
        return 0;
    } else thee->alist = alist;
    if (clist == VNULL) {
        Vnm_print(2, "Vacc_storeParms:  Got NULL Vclist!\n");
        return 0;
    } else thee->clist = clist;
    thee->surf_density = surf_density;

    /* Loop through the atoms to determine the maximum radius */
    maxrad = 0.0;
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        atom = Valist_getAtom(alist, iatom);
        rad = Vatom_getRadius(atom);
        if (rad > maxrad) maxrad = rad;
    }
    maxrad = maxrad + Vclist_maxRadius(thee->clist);

    maxarea = 4.0*VPI*maxrad*maxrad;
    nsphere = (int)ceil(maxarea*surf_density);

    Vnm_print(0, "Vacc_storeParms:  Surf. density = %g\n", surf_density);
    Vnm_print(0, "Vacc_storeParms:  Max area = %g\n", maxarea);
    thee->refSphere = VaccSurf_refSphere(thee->mem, nsphere);
    Vnm_print(0, "Vacc_storeParms:  Using %d-point reference sphere\n",
            thee->refSphere->npts);

    return 1;
}

/** Allocate (and clear) space for storage */
VPRIVATE int Vacc_allocate(Vacc *thee) {

    int i,
        natoms;

    natoms = Valist_getNumberAtoms(thee->alist);

    thee->atomFlags = (int*)Vmem_malloc(thee->mem, natoms, sizeof(int));
    if (thee->atomFlags == VNULL) {
        Vnm_print(2,
               "Vacc_allocate:  Failed to allocate %d (int)s for atomFlags!\n",
                natoms);
        return 0;
    }
    for (i=0; i<natoms; i++) (thee->atomFlags)[i] = 0;

    return 1;
}


VPUBLIC int Vacc_ctor2(Vacc *thee,
                       Valist *alist,
                       Vclist *clist,
                       double surf_density
                       ) {

    /* Check and store parameters */
    if (!Vacc_storeParms(thee, alist, clist, surf_density)) {
        Vnm_print(2, "Vacc_ctor2:  parameter check failed!\n");
        return 0;
    }

    /* Set up memory management object */
    thee->mem = Vmem_ctor("APBS::VACC");
    if (thee->mem == VNULL) {
        Vnm_print(2, "Vacc_ctor2:  memory object setup failed!\n");
        return 0;
    }

    /* Setup and check probe */
    thee->surf = VNULL;

    /* Allocate space */
    if (!Vacc_allocate(thee)) {
        Vnm_print(2, "Vacc_ctor2:  memory allocation failed!\n");
        return 0;
    }

    return 1;
}


VPUBLIC void Vacc_dtor(Vacc **thee) {

    if ((*thee) != VNULL) {
        Vacc_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vacc), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void Vacc_dtor2(Vacc *thee) {

    int i,
        natoms;

    natoms = Valist_getNumberAtoms(thee->alist);
    Vmem_free(thee->mem, natoms, sizeof(int), (void **)&(thee->atomFlags));

    if (thee->refSphere != VNULL) {
        VaccSurf_dtor(&(thee->refSphere));
        thee->refSphere = VNULL;
    }
    if (thee->surf != VNULL) {
        for (i=0; i<natoms; i++) VaccSurf_dtor(&(thee->surf[i]));
        Vmem_free(thee->mem, natoms, sizeof(VaccSurf *),
                (void **)&(thee->surf));
        thee->surf = VNULL;
    }

    Vmem_dtor(&(thee->mem));
}

VPUBLIC double Vacc_vdwAcc(Vacc *thee,
                           double center[3]
                           ) {

    VclistCell *cell;
    Vatom *atom;
    int iatom;
    double *apos,
           dist2;

    /* Get the relevant cell from the cell list */
    cell = Vclist_getCell(thee->clist, center);

    /* If we have no cell, then no atoms are nearby and we're definitely
     * accessible */
    if (cell == VNULL) return 1.0;

    /* Otherwise, check for overlap with the atoms in the cell */
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom = cell->atoms[iatom];
        apos = Vatom_getPosition(atom);
        dist2 = VSQR(center[0]-apos[0]) + VSQR(center[1]-apos[1])
               + VSQR(center[2]-apos[2]);
        if (dist2 < VSQR(Vatom_getRadius(atom))) return 0.0;
    }

    /* If we're still here, then the point is accessible */
    return 1.0;
}

VPUBLIC double Vacc_ivdwAcc(Vacc *thee,
                            double center[3],
                            double radius
                            ) {

    return (double)ivdwAccExclus(thee, center, radius, -1);

}

VPUBLIC void Vacc_splineAccGradAtomNorm(Vacc *thee,
                                        double center[VAPBS_DIM],
                                        double win,
                                        double infrad,
                                        Vatom *atom,
                                        double *grad
                                        ) {

    int i;
    double dist,
           *apos,
           arad,
           sm,
           sm2,
           w2i, /* inverse of win squared */
           w3i, /* inverse of win cubed */
           mygrad,
           mychi = 1.0;           /* Char. func. value for given atom */

    VASSERT(thee != NULL);

    /* Inverse squared window parameter */
    w2i = 1.0/(win*win);
    w3i = 1.0/(win*win*win);

    /* The grad is zero by default */
    for (i=0; i<VAPBS_DIM; i++) grad[i] = 0.0;

    /* *** CALCULATE THE CHARACTERISTIC FUNCTION VALUE FOR THIS ATOM AND THE
     * *** MAGNITUDE OF THE FORCE *** */
    apos = Vatom_getPosition(atom);
    /* Zero-radius atoms don't contribute */
    if (Vatom_getRadius(atom) > 0.0) {
        arad = Vatom_getRadius(atom) + infrad;
        dist = VSQRT(VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
          + VSQR(apos[2]-center[2]));
        /* If we're inside an atom, the entire characteristic function
         * will be zero and the grad will be zero, so we can stop */
        if (dist < (arad - win)) return;
        /* Likewise, if we're outside the smoothing window, the characteristic
         * function is unity and the grad will be zero, so we can stop */
        else if (dist > (arad + win)) return;
        /* Account for floating point error at the border
         * NAB:  COULDN'T THESE TESTS BE COMBINED AS BELOW
         * (Vacc_splineAccAtom)? */
        else if ((VABS(dist - (arad - win)) < VSMALL) ||
                 (VABS(dist - (arad + win)) < VSMALL)) return;
        /* If we're inside the smoothing window */
        else {
            sm = dist - arad + win;
            sm2 = VSQR(sm);
            mychi = 0.75*sm2*w2i -0.25*sm*sm2*w3i;
            mygrad = 1.5*sm*w2i - 0.75*sm2*w3i;
        }
        /* Now assemble the grad vector */
        VASSERT(mychi > 0.0);
        for (i=0; i<VAPBS_DIM; i++)
            grad[i] = -(mygrad/mychi)*((center[i] - apos[i])/dist);
    }
}

VPUBLIC void Vacc_splineAccGradAtomUnnorm(Vacc *thee,
                                          double center[VAPBS_DIM],
                                          double win,
                                          double infrad,
                                          Vatom *atom,
                                          double *grad
                                          ) {

    int i;
    double dist,
           *apos,
           arad,
           sm,
           sm2,
           w2i, /* Inverse of win squared */
           w3i, /* Inverse of win cubed */
           mygrad,
           mychi = 1.0;           /* Char. func. value for given atom */

    VASSERT(thee != NULL);

    /* Inverse squared window parameter */
    w2i = 1.0/(win*win);
    w3i = 1.0/(win*win*win);

    /* The grad is zero by default */
    for (i=0; i<VAPBS_DIM; i++) grad[i] = 0.0;

    /* *** CALCULATE THE CHARACTERISTIC FUNCTION VALUE FOR THIS ATOM AND THE
     * *** MAGNITUDE OF THE FORCE *** */
    apos = Vatom_getPosition(atom);
    /* Zero-radius atoms don't contribute */
    if (Vatom_getRadius(atom) > 0.0) {
        arad = Vatom_getRadius(atom) + infrad;
        dist = VSQRT(VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
          + VSQR(apos[2]-center[2]));
        /* If we're inside an atom, the entire characteristic function
         * will be zero and the grad will be zero, so we can stop */
        if (dist < (arad - win)) return;
        /* Likewise, if we're outside the smoothing window, the characteristic
         * function is unity and the grad will be zero, so we can stop */
        else if (dist > (arad + win)) return;
        /* Account for floating point error at the border
         * NAB:  COULDN'T THESE TESTS BE COMBINED AS BELOW
         * (Vacc_splineAccAtom)? */
        else if ((VABS(dist - (arad - win)) < VSMALL) ||
                 (VABS(dist - (arad + win)) < VSMALL)) return;
        /* If we're inside the smoothing window */
        else {
            sm = dist - arad + win;
            sm2 = VSQR(sm);
            mychi = 0.75*sm2*w2i -0.25*sm*sm2*w3i;
            mygrad = 1.5*sm*w2i - 0.75*sm2*w3i;
        }
        /* Now assemble the grad vector */
        VASSERT(mychi > 0.0);
        for (i=0; i<VAPBS_DIM; i++)
            grad[i] = -(mygrad)*((center[i] - apos[i])/dist);
    }
}

VPUBLIC double Vacc_splineAccAtom(Vacc *thee,
                                  double center[VAPBS_DIM],
                                  double win,
                                  double infrad,
                                  Vatom *atom
                                  ) {

    double dist,
           *apos,
           arad,
           sm,
           sm2,
           w2i, /* Inverse of win squared */
           w3i, /* Inverse of win cubed */
           value,
           stot,
           sctot;

    VASSERT(thee != NULL);

    /* Inverse squared window parameter */
    w2i = 1.0/(win*win);
    w3i = 1.0/(win*win*win);

    apos = Vatom_getPosition(atom);
    /* Zero-radius atoms don't contribute */
    if (Vatom_getRadius(atom) > 0.0) {
        arad = Vatom_getRadius(atom) + infrad;
        stot = arad + win;
        sctot = VMAX2(0, (arad - win));
        dist = VSQRT(VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
          + VSQR(apos[2]-center[2]));
        /* If we're inside an atom, the entire characteristic function
         * will be zero */
        if ((dist < sctot) || (VABS(dist - sctot) < VSMALL)){
            value = 0.0;
        /* We're outside the smoothing window */
        } else if ((dist > stot) || (VABS(dist - stot) < VSMALL)) {
            value = 1.0;
        /* We're inside the smoothing window */
        } else {
            sm = dist - arad + win;
            sm2 = VSQR(sm);
            value = 0.75*sm2*w2i - 0.25*sm*sm2*w3i;
        }
    } else value = 1.0;

    return value;
}

/**
 * @brief  Fast spline-based surface computation subroutine
 * @returns  Spline value
 * @author  Todd Dolinsky and Nathan Baker
 */
VPRIVATE double splineAcc(
        Vacc *thee,  /** Accessibility object */
        double center[VAPBS_DIM],  /** Point at which the acc is to be
                                    * evaluated */
        double win,  /** Spline window */
        double infrad,  /** Radius to inflate atomic radius */
        VclistCell *cell  /** Cell of atom objects */
        ) {

    int atomID, iatom;
    Vatom *atom;
    double value = 1.0;

    VASSERT(thee != NULL);

    /* Now loop through the atoms assembling the characteristic function */
    for (iatom=0; iatom<cell->natoms; iatom++) {

        atom = cell->atoms[iatom];
        atomID = atom->id;

        /* Check to see if we've counted this atom already */
        if ( !(thee->atomFlags[atomID]) ) {

            thee->atomFlags[atomID] = 1;
            value *= Vacc_splineAccAtom(thee, center, win, infrad, atom);

            if (value < VSMALL) return value;
        }
    }

    return value;
}


VPUBLIC double Vacc_splineAcc(Vacc *thee, double center[VAPBS_DIM], double win,
  double infrad) {

    VclistCell *cell;
    Vatom *atom;
    int iatom, atomID;


    VASSERT(thee != NULL);

    if (Vclist_maxRadius(thee->clist) < (win + infrad)) {
        Vnm_print(2, "Vacc_splineAcc:  Vclist has max_radius=%g;\n",
                Vclist_maxRadius(thee->clist));
        Vnm_print(2, "Vacc_splineAcc:  Insufficient for win=%g, infrad=%g\n",
                win, infrad);
        VASSERT(0);
    }

    /* Get a cell or VNULL; in the latter case return 1.0 */
    cell = Vclist_getCell(thee->clist, center);
    if (cell == VNULL) return 1.0;

    /* First, reset the list of atom flags
     * NAB:  THIS SEEMS VERY INEFFICIENT */
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom = cell->atoms[iatom];
        atomID = atom->id;
        thee->atomFlags[atomID] = 0;
    }

    return splineAcc(thee, center, win, infrad, cell);
}

VPUBLIC void Vacc_splineAccGrad(Vacc *thee, double center[VAPBS_DIM],
        double win, double infrad, double *grad) {

    int iatom, i, atomID;
    double acc = 1.0;
    double tgrad[VAPBS_DIM];
    VclistCell *cell;
    Vatom *atom = VNULL;

    VASSERT(thee != NULL);

    if (Vclist_maxRadius(thee->clist) < (win + infrad)) {
        Vnm_print(2, "Vacc_splineAccGrad: Vclist max_radius=%g;\n",
                Vclist_maxRadius(thee->clist));
        Vnm_print(2, "Vacc_splineAccGrad: Insufficient for win=%g, infrad=%g\n",
                win, infrad);
        VASSERT(0);
    }

    /* Reset the gradient */
    for (i=0; i<VAPBS_DIM; i++) grad[i] = 0.0;

    /* Get the cell; check for nullity */
    cell = Vclist_getCell(thee->clist, center);
    if (cell == VNULL) return;

    /* Reset the list of atom flags */
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom = cell->atoms[iatom];
        atomID = atom->id;
        thee->atomFlags[atomID] = 0;
    }

    /* Get the local accessibility */
    acc = splineAcc(thee, center, win, infrad, cell);

    /* Accumulate the gradient of all local atoms */
    if (acc > VSMALL) {
        for (iatom=0; iatom<cell->natoms; iatom++) {
            atom = cell->atoms[iatom];
            Vacc_splineAccGradAtomNorm(thee, center, win, infrad, atom, tgrad);
        }
        for (i=0; i<VAPBS_DIM; i++) grad[i] += tgrad[i];
    }
    for (i=0; i<VAPBS_DIM; i++) grad[i] *= -acc;
}

VPUBLIC double Vacc_molAcc(Vacc *thee, double center[VAPBS_DIM],
        double radius) {

    double rc;

    /* ******* CHECK IF OUTSIDE ATOM+PROBE RADIUS SURFACE ***** */
    if (Vacc_ivdwAcc(thee, center, radius) == 1.0) {

        /* Vnm_print(2, "DEBUG:  ivdwAcc = 1.0\n"); */
        rc = 1.0;

    /* ******* CHECK IF INSIDE ATOM RADIUS SURFACE ***** */
    } else if (Vacc_vdwAcc(thee, center) == 0.0) {

        /* Vnm_print(2, "DEBUG:  vdwAcc = 0.0\n"); */
        rc = 0.0;

    /* ******* CHECK IF OUTSIDE MOLECULAR SURFACE ***** */
    } else {

        /* Vnm_print(2, "DEBUG:  calling fastMolAcc...\n"); */
        rc = Vacc_fastMolAcc(thee, center, radius);

    }

    return rc;

}

VPUBLIC double Vacc_fastMolAcc(Vacc *thee, double center[VAPBS_DIM],
        double radius) {

    Vatom *atom;
    VaccSurf *surf;
    VclistCell *cell;
    int ipt, iatom, atomID;
    double dist2, rad2;

    rad2 = radius*radius;

    /* Check to see if the SAS has been defined */
    if (thee->surf == VNULL) Vacc_SASA(thee, radius);

    /* Get the cell associated with this point */
    cell = Vclist_getCell(thee->clist, center);
    if (cell == VNULL) {
        Vnm_print(2, "Vacc_fastMolAcc:  unexpected VNULL VclistCell!\n");
        return 1.0;
    }

    /* Loop through all the atoms in the cell */
    for (iatom=0; iatom<cell->natoms; iatom++) {
        atom = cell->atoms[iatom];
        atomID = Vatom_getAtomID(atom);
        surf = thee->surf[atomID];
        /* Loop through all SAS points associated with this atom */
        for (ipt=0; ipt<surf->npts; ipt++) {
            /* See if we're within a probe radius of the point */
            dist2 = VSQR(center[0]-(surf->xpts[ipt]))
                + VSQR(center[1]-(surf->ypts[ipt]))
                + VSQR(center[2]-(surf->zpts[ipt]));
            if (dist2 < rad2) return 1.0;
        }
    }

    /* If all else failed, we are not inside the molecular surface */
    return 0.0;
}


#if defined(HAVE_MC_H)
VPUBLIC void Vacc_writeGMV(Vacc *thee, double radius, int meth, Gem *gm,
  char *iodev, char *iofmt, char *iohost, char *iofile) {

    double *accVals[MAXV], coord[3];
    Vio *sock;
    int ivert, icoord;

    for (ivert=0; ivert<MAXV; ivert++) accVals[ivert] = VNULL;
    accVals[0] = (void *)Vmem_malloc(thee->mem, Gem_numVV(gm), sizeof(double));
    accVals[1] = (void *)Vmem_malloc(thee->mem, Gem_numVV(gm), sizeof(double));
    for (ivert=0; ivert<Gem_numVV(gm); ivert++) {
        for (icoord=0;icoord<3;icoord++)
          coord[icoord] = VV_coord(Gem_VV(gm, ivert), icoord);
        if (meth == 0) {
            accVals[0][ivert] = Vacc_molAcc(thee, coord, radius);
            accVals[1][ivert] = Vacc_molAcc(thee, coord, radius);
        } else if (meth == 1) {
            accVals[0][ivert] = Vacc_ivdwAcc(thee, coord, radius);
            accVals[1][ivert] = Vacc_ivdwAcc(thee, coord, radius);
        } else if (meth == 2) {
            accVals[0][ivert] = Vacc_vdwAcc(thee, coord);
            accVals[1][ivert] = Vacc_vdwAcc(thee, coord);
        } else VASSERT(0);
    }
    sock = Vio_ctor(iodev, iofmt, iohost, iofile, "w");
    Gem_writeGMV(gm, sock, 1, accVals);
    Vio_dtor(&sock);
    Vmem_free(thee->mem, Gem_numVV(gm), sizeof(double),
      (void **)&(accVals[0]));
    Vmem_free(thee->mem, Gem_numVV(gm), sizeof(double),
      (void **)&(accVals[1]));
}
#endif /* defined(HAVE_MC_H) */

VPUBLIC double Vacc_SASA(Vacc *thee,
                         double radius
                         ) {

    int i,
        natom;
    double area;
           //*apos; // gcc says unused
    Vatom *atom;
    VaccSurf *asurf;

    time_t ts; // PCE: temp
    ts = clock();

    //unsigned long long mbeg; // gcc says unused

    natom = Valist_getNumberAtoms(thee->alist);

    /* Check to see if we need to build the surface */
    if (thee->surf == VNULL) {
        thee->surf = Vmem_malloc(thee->mem, natom, sizeof(VaccSurf *));

#if defined(DEBUG_MAC_OSX_OCL) || defined(DEBUG_MAC_OSX_STANDARD)
#include "mach_chud.h"
        machm_(&mbeg);
#pragma omp parallel for private(i,atom)
#endif
        for (i=0; i<natom; i++) {
            atom = Valist_getAtom(thee->alist, i);
            /* NOTE:  RIGHT NOW WE DO THIS FOR THE ENTIRE MOLECULE WHICH IS
             * INCREDIBLY INEFFICIENT, PARTICULARLY DURING FOCUSING!!! */
            thee->surf[i] = Vacc_atomSurf(thee, atom, thee->refSphere,
                                          radius);
        }
    }

    /* Calculate the area */
    area = 0.0;
    for (i=0; i<natom; i++) {
        atom = Valist_getAtom(thee->alist, i);
        asurf = thee->surf[i];
        /* See if this surface needs to be rebuilt */
        if (asurf->probe_radius != radius) {
            Vnm_print(2, "Vacc_SASA:  Warning -- probe radius changed from %g to %g!\n",
                      asurf->probe_radius, radius);
            VaccSurf_dtor2(asurf);
            thee->surf[i] = Vacc_atomSurf(thee, atom, thee->refSphere, radius);
            asurf = thee->surf[i];
        }
        area += (asurf->area);
    }

#if defined(DEBUG_MAC_OSX_OCL) || defined(DEBUG_MAC_OSX_STANDARD)
    mets_(&mbeg, "Vacc_SASA - Parallel");
#endif

    Vnm_print(0, "Vacc_SASA: Time elapsed: %f\n", ((double)clock() - ts) / CLOCKS_PER_SEC);
    return area;

}

VPUBLIC double Vacc_totalSASA(Vacc *thee, double radius) {

    return Vacc_SASA(thee, radius);

}

VPUBLIC double Vacc_atomSASA(Vacc *thee, double radius, Vatom *atom) {

    VaccSurf *asurf;
    int id;

    if (thee->surf == VNULL) Vacc_SASA(thee, radius);

    id = Vatom_getAtomID(atom);
    asurf = thee->surf[id];

    /* See if this surface needs to be rebuilt */
    if (asurf->probe_radius != radius) {
        Vnm_print(2, "Vacc_SASA:  Warning -- probe radius changed from %g to %g!\n",
                asurf->probe_radius, radius);
        VaccSurf_dtor2(asurf);
        thee->surf[id] = Vacc_atomSurf(thee, atom, thee->refSphere, radius);
        asurf = thee->surf[id];
    }

    return asurf->area;

}

VPUBLIC VaccSurf* VaccSurf_ctor(Vmem *mem, double probe_radius, int nsphere) {
    VaccSurf *thee;

    //thee = Vmem_malloc(mem, 1, sizeof(Vacc) );
    if (nsphere >= MAX_SPHERE_PTS) {
        Vnm_print(2, "VaccSurf_ctor:  Error!  The requested number of grid points (%d) exceeds the maximum (%d)!\n", nsphere, MAX_SPHERE_PTS);
        Vnm_print(2, "VaccSurf_ctor:  Please check the variable MAX_SPHERE_PTS to reset.\n");
        VASSERT(0);
    }
    thee = (VaccSurf*)calloc(1,sizeof(Vacc));
    VASSERT( VaccSurf_ctor2(thee, mem, probe_radius, nsphere) );

    return thee;
}

VPUBLIC int VaccSurf_ctor2(VaccSurf *thee, Vmem *mem, double probe_radius,
        int nsphere) {

    if (thee == VNULL)
        return 0;

    thee->mem = mem;
    thee->npts = nsphere;
    thee->probe_radius = probe_radius;
    thee->area = 0.0;

    if (thee->npts > 0) {
        thee->xpts = Vmem_malloc(thee->mem, thee->npts, sizeof(double));
        thee->ypts = Vmem_malloc(thee->mem, thee->npts, sizeof(double));
        thee->zpts = Vmem_malloc(thee->mem, thee->npts, sizeof(double));
        thee->bpts = Vmem_malloc(thee->mem, thee->npts, sizeof(char));
    } else {
        thee->xpts = VNULL;
        thee->ypts = VNULL;
        thee->zpts = VNULL;
        thee->bpts = VNULL;
    }

    return 1;
}

VPUBLIC void VaccSurf_dtor(VaccSurf **thee) {

    Vmem *mem;

    if ((*thee) != VNULL) {
        mem = (*thee)->mem;
        VaccSurf_dtor2(*thee);
        //Vmem_free(mem, 1, sizeof(VaccSurf), (void **)thee);
        free(*thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void VaccSurf_dtor2(VaccSurf *thee) {

    if (thee->npts > 0) {
        Vmem_free(thee->mem, thee->npts, sizeof(double),
                (void **)&(thee->xpts));
        Vmem_free(thee->mem, thee->npts, sizeof(double),
                (void **)&(thee->ypts));
        Vmem_free(thee->mem, thee->npts, sizeof(double),
                (void **)&(thee->zpts));
        Vmem_free(thee->mem, thee->npts, sizeof(char),
                (void **)&(thee->bpts));
    }

}

VPUBLIC VaccSurf* Vacc_atomSurf(Vacc *thee,
                                Vatom *atom,
                                VaccSurf *ref,
                                double prad
                               ) {

    VaccSurf *surf;
    int i,
        j,
        npts,
        atomID;
    double arad,
           rad,
           pos[3],
           *apos;
    char bpts[MAX_SPHERE_PTS];

    /* Get atom information */
    arad = Vatom_getRadius(atom);
    apos = Vatom_getPosition(atom);
    atomID = Vatom_getAtomID(atom);

    if (arad < VSMALL) {
        return VaccSurf_ctor(thee->mem, prad, 0);
    }

    rad = arad + prad;

    /* Determine which points will contribute */
    npts = 0;
    for (i=0; i<ref->npts; i++) {
        /* Reset point flag: zero-radius atoms do not contribute */
        pos[0] = rad*(ref->xpts[i]) + apos[0];
        pos[1] = rad*(ref->ypts[i]) + apos[1];
        pos[2] = rad*(ref->zpts[i]) + apos[2];
        if (ivdwAccExclus(thee, pos, prad, atomID)) {
            npts++;
            bpts[i] = 1;
        } else {
            bpts[i] = 0;
        }
    }

    /* Allocate space for the points */
    surf = VaccSurf_ctor(thee->mem, prad, npts);

    /* Assign the points */
    j = 0;
    for (i=0; i<ref->npts; i++) {
        if (bpts[i]) {
            surf->bpts[j] = 1;
            surf->xpts[j] = rad*(ref->xpts[i]) + apos[0];
            surf->ypts[j] = rad*(ref->ypts[i]) + apos[1];
            surf->zpts[j] = rad*(ref->zpts[i]) + apos[2];
            j++;
        }
    }

    /* Assign the area */
    surf->area = 4.0*VPI*rad*rad*((double)(surf->npts))/((double)(ref->npts));

    return surf;

}

VPUBLIC VaccSurf* VaccSurf_refSphere(Vmem *mem, int npts) {

    VaccSurf *surf;
    int nactual, i, itheta, ntheta, iphi, nphimax, nphi;
    double frac;
    double sintheta, costheta, theta, dtheta;
    double sinphi, cosphi, phi, dphi;

    /* Setup "constants" */
    frac = ((double)(npts))/4.0;
    ntheta = VRINT(VSQRT(Vunit_pi*frac));
    dtheta = Vunit_pi/((double)(ntheta));
    nphimax = 2*ntheta;

    /* Count the actual number of points to be used */
    nactual = 0;
    for (itheta=0; itheta<ntheta; itheta++) {
        theta = dtheta*((double)(itheta));
        sintheta = VSIN(theta);
        costheta = VCOS(theta);
        nphi = VRINT(sintheta*nphimax);
        nactual += nphi;
    }

    /* Allocate space for the points */
    surf = VaccSurf_ctor(mem, 1.0, nactual);

    /* Clear out the boolean array */
    for (i=0; i<nactual; i++) surf->bpts[i] = 1;

    /* Assign the points */
    nactual = 0;
    for (itheta=0; itheta<ntheta; itheta++) {
        theta = dtheta*((double)(itheta));
        sintheta = VSIN(theta);
        costheta = VCOS(theta);
        nphi = VRINT(sintheta*nphimax);
        if (nphi != 0) {
            dphi = 2*Vunit_pi/((double)(nphi));
            for (iphi=0; iphi<nphi; iphi++) {
                phi = dphi*((double)(iphi));
                sinphi = VSIN(phi);
                cosphi = VCOS(phi);
                surf->xpts[nactual] = cosphi * sintheta;
                surf->ypts[nactual] = sinphi * sintheta;
                surf->zpts[nactual] = costheta;
                nactual++;
            }
        }
    }

    surf->npts = nactual;

    return surf;
}

VPUBLIC VaccSurf* Vacc_atomSASPoints(Vacc *thee, double radius,
        Vatom *atom) {

    VaccSurf *asurf = VNULL;
    int id;

    if (thee->surf == VNULL) Vacc_SASA(thee, radius);
    id = Vatom_getAtomID(atom);

    asurf = thee->surf[id];

    /* See if this surface needs to be rebuilt */
    if (asurf->probe_radius != radius) {
        Vnm_print(2, "Vacc_SASA:  Warning -- probe radius changed from %g to %g!\n",
                asurf->probe_radius, radius);
        VaccSurf_dtor2(asurf);
        thee->surf[id] = Vacc_atomSurf(thee, atom, thee->refSphere, radius);
        asurf = thee->surf[id];
    }

    return asurf;

}

VPUBLIC void Vacc_splineAccGradAtomNorm4(Vacc *thee, double center[VAPBS_DIM],
                                         double win, double infrad, Vatom *atom, double *grad) {

    int i;
    double dist, *apos, arad, sm, sm2, sm3, sm4, sm5, sm6, sm7;
    double e, e2, e3, e4, e5, e6, e7;
    double b, b2, b3, b4, b5, b6, b7;
    double c0, c1, c2, c3, c4, c5, c6, c7;
    double denom, mygrad;
    double mychi = 1.0;           /* Char. func. value for given atom */

    VASSERT(thee != NULL);

    /* The grad is zero by default */
    for (i=0; i<VAPBS_DIM; i++) grad[i] = 0.0;

    /* *** CALCULATE THE CHARACTERISTIC FUNCTION VALUE FOR THIS ATOM AND THE
        * *** MAGNITUDE OF THE FORCE *** */
    apos = Vatom_getPosition(atom);
    /* Zero-radius atoms don't contribute */
    if (Vatom_getRadius(atom) > 0.0) {

        arad = Vatom_getRadius(atom);
        arad = arad + infrad;
        b = arad - win;
        e = arad + win;

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

        dist = VSQRT(VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
                     + VSQR(apos[2]-center[2]));

        /* If we're inside an atom, the entire characteristic function
            * will be zero and the grad will be zero, so we can stop */
        if (dist < (arad - win)) return;
        /* Likewise, if we're outside the smoothing window, the characteristic
            * function is unity and the grad will be zero, so we can stop */
        else if (dist > (arad + win)) return;
        /* Account for floating point error at the border
            * NAB:  COULDN'T THESE TESTS BE COMBINED AS BELOW
            * (Vacc_splineAccAtom)? */
        else if ((VABS(dist - (arad - win)) < VSMALL) ||
                 (VABS(dist - (arad + win)) < VSMALL)) return;
        /* If we're inside the smoothing window */
        else {
            sm = dist;
            sm2 = sm * sm;
            sm3 = sm2 * sm;
            sm4 = sm3 * sm;
            sm5 = sm4 * sm;
            sm6 = sm5 * sm;
            sm7 = sm6 * sm;
            mychi = c0 + c1*sm + c2*sm2 + c3*sm3
                + c4*sm4 + c5*sm5 + c6*sm6 + c7*sm7;
            mygrad = c1 + 2.0*c2*sm  + 3.0*c3*sm2 + 4.0*c4*sm3
                + 5.0*c5*sm4 + 6.0*c6*sm5 + 7.0*c7*sm6;
            if (mychi <= 0.0) {
                /* Avoid numerical round off errors */
                return;
            } else if (mychi > 1.0) {
                /* Avoid numerical round off errors */
                mychi = 1.0;
            }
        }
        /* Now assemble the grad vector */
        VASSERT(mychi > 0.0);
        for (i=0; i<VAPBS_DIM; i++)
            grad[i] = -(mygrad/mychi)*((center[i] - apos[i])/dist);
    }
}

VPUBLIC void Vacc_splineAccGradAtomNorm3(Vacc *thee, double center[VAPBS_DIM],
                                         double win, double infrad, Vatom *atom, double *grad) {

    int i;
    double dist, *apos, arad, sm, sm2, sm3, sm4, sm5;
    double e, e2, e3, e4, e5;
    double b, b2, b3, b4, b5;
    double c0, c1, c2, c3, c4, c5;
    double denom, mygrad;
    double mychi = 1.0;           /* Char. func. value for given atom */

    VASSERT(thee != NULL);

    /* The grad is zero by default */
    for (i=0; i<VAPBS_DIM; i++) grad[i] = 0.0;

    /* *** CALCULATE THE CHARACTERISTIC FUNCTION VALUE FOR THIS ATOM AND THE
        * *** MAGNITUDE OF THE FORCE *** */
    apos = Vatom_getPosition(atom);
    /* Zero-radius atoms don't contribute */
    if (Vatom_getRadius(atom) > 0.0) {

        arad = Vatom_getRadius(atom);
        arad = arad + infrad;
        b = arad - win;
        e = arad + win;

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

        dist = VSQRT(VSQR(apos[0]-center[0]) + VSQR(apos[1]-center[1])
                     + VSQR(apos[2]-center[2]));

        /* If we're inside an atom, the entire characteristic function
            * will be zero and the grad will be zero, so we can stop */
        if (dist < (arad - win)) return;
        /* Likewise, if we're outside the smoothing window, the characteristic
            * function is unity and the grad will be zero, so we can stop */
        else if (dist > (arad + win)) return;
        /* Account for floating point error at the border
            * NAB:  COULDN'T THESE TESTS BE COMBINED AS BELOW
            * (Vacc_splineAccAtom)? */
        else if ((VABS(dist - (arad - win)) < VSMALL) ||
                 (VABS(dist - (arad + win)) < VSMALL)) return;
        /* If we're inside the smoothing window */
        else {
            sm = dist;
            sm2 = sm * sm;
            sm3 = sm2 * sm;
            sm4 = sm3 * sm;
            sm5 = sm4 * sm;
            mychi = c0 + c1*sm + c2*sm2 + c3*sm3
                + c4*sm4 + c5*sm5;
            mygrad = c1 + 2.0*c2*sm  + 3.0*c3*sm2 + 4.0*c4*sm3
                + 5.0*c5*sm4;
            if (mychi <= 0.0) {
                /* Avoid numerical round off errors */
                return;
            } else if (mychi > 1.0) {
                /* Avoid numerical round off errors */
                mychi = 1.0;
            }
        }
        /* Now assemble the grad vector */
        VASSERT(mychi > 0.0);
        for (i=0; i<VAPBS_DIM; i++)
            grad[i] = -(mygrad/mychi)*((center[i] - apos[i])/dist);
    }
}

/* ///////////////////////////////////////////////////////////////////////////
   // Routine:  Vacc_atomdSAV
   //
   // Purpose:  Calculates the vector valued atomic derivative of volume
   //
   // Args:     radius  The radius of the solvent probe in Angstroms
   //           iatom   Index of the atom in thee->alist
   //
   // Author:   Jason Wagoner
   //           Nathan Baker (original FORTRAN routine from UHBD by Brock Luty)
   /////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vacc_atomdSAV(Vacc *thee,
                           double srad,
                           Vatom *atom,
                           double *dSA
                          ) {

    int ipt, iatom;

    double area;
    double *tPos, tRad, vec[3];
    double dx,dy,dz;
    VaccSurf *ref;
    dx = 0.0;
    dy = 0.0;
    dz = 0.0;
    /* Get the atom information */
    ref = thee->refSphere;
    iatom = Vatom_getAtomID(atom);

    dSA[0] = 0.0;
    dSA[1] = 0.0;
    dSA[2] = 0.0;

    tPos = Vatom_getPosition(atom);
    tRad = Vatom_getRadius(atom);

    if(tRad == 0.0) return;

    area = 4.0*VPI*(tRad+srad)*(tRad+srad)/((double)(ref->npts));
    for (ipt=0; ipt<ref->npts; ipt++) {
        vec[0] = (tRad+srad)*ref->xpts[ipt] + tPos[0];
        vec[1] = (tRad+srad)*ref->ypts[ipt] + tPos[1];
        vec[2] = (tRad+srad)*ref->zpts[ipt] + tPos[2];
        if (ivdwAccExclus(thee, vec, srad, iatom)) {
            dx = dx+vec[0]-tPos[0];
            dy = dy+vec[1]-tPos[1];
            dz = dz+vec[2]-tPos[2];
        }
    }

    if ((tRad+srad) != 0){
        dSA[0] = dx*area/(tRad+srad);
        dSA[1] = dy*area/(tRad+srad);
        dSA[2] = dz*area/(tRad+srad);
    }

}

/* Note: This is purely test code to make certain that the dSASA code is
         behaving properly. This function should NEVER be called by anyone
         other than an APBS developer at Wash U.
*/
VPRIVATE double Vacc_SASAPos(Vacc *thee, double radius) {

    int i, natom;
    double area;
    Vatom *atom;
    VaccSurf *asurf;

    natom = Valist_getNumberAtoms(thee->alist);

    /* Calculate the area */
    area = 0.0;
    for (i=0; i<natom; i++) {
        atom = Valist_getAtom(thee->alist, i);
        asurf = thee->surf[i];

        VaccSurf_dtor2(asurf);
        thee->surf[i] = Vacc_atomSurf(thee, atom, thee->refSphere, radius);
        asurf = thee->surf[i];
        area += (asurf->area);
    }

    return area;

}

VPRIVATE double Vacc_atomSASAPos(Vacc *thee,
                                 double radius,
                                 Vatom *atom, /* The atom being manipulated */
                                 int mode
                                ) {

    VaccSurf *asurf;
    int id;
    static int warned = 0;

    if ((thee->surf == VNULL) || (mode == 1)){
        if(!warned){
            Vnm_print(2, "WARNING: Recalculating entire surface!!!!\n");
            warned = 1;
        }
        Vacc_SASAPos(thee, radius); // reinitialize before we can do anything about doing a calculation on a repositioned atom
    }

    id = Vatom_getAtomID(atom);
    asurf = thee->surf[id];

    VaccSurf_dtor(&asurf);
    thee->surf[id] = Vacc_atomSurf(thee, atom, thee->refSphere, radius);
    asurf = thee->surf[id];

    //printf("%s: Time elapsed: %f\n", __func__, ((double)clock() - ts) / CLOCKS_PER_SEC);

    return asurf->area;
}

/* ///////////////////////////////////////////////////////////////////////////
   // Routine:  Vacc_atomdSASA
   //
   // Purpose:  Calculates the derivative of surface area with respect to atomic
   //           displacement using finite difference methods.
   //
   // Args:     radius  The radius of the solvent probe in Angstroms
   //           iatom   Index of the atom in thee->alist
   //
   // Author:   Jason Wagoner
   //			David Gohara
   //           Nathan Baker (original FORTRAN routine from UHBD by Brock Luty)
   /////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vacc_atomdSASA(Vacc *thee,
                            double dpos,
                            double srad,
                            Vatom *atom,
                            double *dSA
                           ) {

    double *temp_Pos,
           tPos[3],
           axb1,
           axt1,
           ayb1,
           ayt1,
           azb1,
           azt1;
    VaccSurf *ref;

    //printf("%s: entering\n", __func__);
    time_t ts;
    ts = clock();

    /* Get the atom information */
    ref = thee->refSphere;
    temp_Pos = Vatom_getPosition(atom); // Get a pointer to the position object.  You actually manipulate the atom doing this...

    tPos[0] = temp_Pos[0];
    tPos[1] = temp_Pos[1];
    tPos[2] = temp_Pos[2];

    /* Shift by pos -/+ on x */
    temp_Pos[0] -= dpos;
    axb1 = Vacc_atomSASAPos(thee, srad, atom,0);
    temp_Pos[0] = tPos[0];

    temp_Pos[0] += dpos;
    axt1 = Vacc_atomSASAPos(thee, srad, atom,0);
    temp_Pos[0] = tPos[0];

    /* Shift by pos -/+ on y */
    temp_Pos[1] -= dpos;
    ayb1 = Vacc_atomSASAPos(thee, srad, atom,0);
    temp_Pos[1] = tPos[1];

    temp_Pos[1] += dpos;
    ayt1 = Vacc_atomSASAPos(thee, srad, atom,0);
    temp_Pos[1] = tPos[1];

    /* Shift by pos -/+ on z */
    temp_Pos[2] -= dpos;
    azb1 = Vacc_atomSASAPos(thee, srad, atom,0);
    temp_Pos[2] = tPos[2];

    temp_Pos[2] += dpos;
    azt1 = Vacc_atomSASAPos(thee, srad, atom,0);
    temp_Pos[2] = tPos[2];

    /* Reset the atom SASA to zero displacement */
    Vacc_atomSASAPos(thee, srad, atom,0);

    /* Calculate the final value */
    dSA[0] = (axt1-axb1)/(2.0 * dpos);
    dSA[1] = (ayt1-ayb1)/(2.0 * dpos);
    dSA[2] = (azt1-azb1)/(2.0 * dpos);
}

/* Note: This is purely test code to make certain that the dSASA code is
         behaving properly. This function should NEVER be called by anyone
         other than an APBS developer at Wash U.
*/
VPUBLIC void Vacc_totalAtomdSASA(Vacc *thee, double dpos, double srad, Vatom *atom, double *dSA) {

    int iatom;
    double *temp_Pos, tRad;
    double tPos[3];
    double axb1,axt1,ayb1,ayt1,azb1,azt1;
    VaccSurf *ref;

    /* Get the atom information */
    ref = thee->refSphere;
    temp_Pos = Vatom_getPosition(atom);
    tRad = Vatom_getRadius(atom);
    iatom = Vatom_getAtomID(atom);

    dSA[0] = 0.0;
    dSA[1] = 0.0;
    dSA[2] = 0.0;

    tPos[0] = temp_Pos[0];
    tPos[1] = temp_Pos[1];
    tPos[2] = temp_Pos[2];

    /* Shift by pos -/+ on x */
    temp_Pos[0] -= dpos;
    axb1 = Vacc_atomSASAPos(thee, srad, atom, 1);
    temp_Pos[0] = tPos[0];

    temp_Pos[0] += dpos;
    axt1 = Vacc_atomSASAPos(thee, srad, atom, 1);
    temp_Pos[0] = tPos[0];

    /* Shift by pos -/+ on y */
    temp_Pos[1] -= dpos;
    ayb1 = Vacc_atomSASAPos(thee, srad, atom, 1);
    temp_Pos[1] = tPos[1];

    temp_Pos[1] += dpos;
    ayt1 = Vacc_atomSASAPos(thee, srad, atom, 1);
    temp_Pos[1] = tPos[1];

    /* Shift by pos -/+ on z */
    temp_Pos[2] -= dpos;
    azb1 = Vacc_atomSASAPos(thee, srad, atom, 1);
    temp_Pos[2] = tPos[2];

    temp_Pos[2] += dpos;
    azt1 = Vacc_atomSASAPos(thee, srad, atom, 1);
    temp_Pos[2] = tPos[2];

    /* Calculate the final value */
    dSA[0] = (axt1-axb1)/(2.0 * dpos);
    dSA[1] = (ayt1-ayb1)/(2.0 * dpos);
    dSA[2] = (azt1-azb1)/(2.0 * dpos);
}

/* Note: This is purely test code to make certain that the dSASA code is
         behaving properly. This function should NEVER be called by anyone
         other than an APBS developer at Wash U.
*/
VPUBLIC void Vacc_totalAtomdSAV(Vacc *thee, double dpos, double srad, Vatom *atom, double *dSA, Vclist *clist) {

    int iatom;
    double *temp_Pos, tRad;
    double tPos[3];
    double axb1,axt1,ayb1,ayt1,azb1,azt1;
    VaccSurf *ref;

    /* Get the atom information */
    ref = thee->refSphere;
    temp_Pos = Vatom_getPosition(atom);
    tRad = Vatom_getRadius(atom);
    iatom = Vatom_getAtomID(atom);

    dSA[0] = 0.0;
    dSA[1] = 0.0;
    dSA[2] = 0.0;

    tPos[0] = temp_Pos[0];
    tPos[1] = temp_Pos[1];
    tPos[2] = temp_Pos[2];

    /* Shift by pos -/+ on x */
    temp_Pos[0] -= dpos;
    axb1 = Vacc_totalSAV(thee,clist, VNULL, srad);
    temp_Pos[0] = tPos[0];

    temp_Pos[0] += dpos;
    axt1 = Vacc_totalSAV(thee,clist, VNULL, srad);
    temp_Pos[0] = tPos[0];

    /* Shift by pos -/+ on y */
    temp_Pos[1] -= dpos;
    ayb1 = Vacc_totalSAV(thee,clist, VNULL, srad);
    temp_Pos[1] = tPos[1];

    temp_Pos[1] += dpos;
    ayt1 = Vacc_totalSAV(thee,clist, VNULL, srad);
    temp_Pos[1] = tPos[1];

    /* Shift by pos -/+ on z */
    temp_Pos[2] -= dpos;
    azb1 = Vacc_totalSAV(thee,clist, VNULL, srad);
    temp_Pos[2] = tPos[2];

    temp_Pos[2] += dpos;
    azt1 = Vacc_totalSAV(thee,clist, VNULL, srad);
    temp_Pos[2] = tPos[2];

    /* Calculate the final value */
    dSA[0] = (axt1-axb1)/(2.0 * dpos);
    dSA[1] = (ayt1-ayb1)/(2.0 * dpos);
    dSA[2] = (azt1-azb1)/(2.0 * dpos);
}

VPUBLIC double Vacc_totalSAV(Vacc *thee, Vclist *clist, APOLparm *apolparm, double radius) {

    int i;
    int npts[3];

    double spacs[3], vec[3];
    double w, wx, wy, wz, len, fn, x, y, z, vol;
    double vol_density,sav;
    double *lower_corner, *upper_corner;

    sav = 0.0;
    vol = 1.0;
    vol_density = 2.0;

    lower_corner = clist->lower_corner;
    upper_corner = clist->upper_corner;

    for (i=0; i<3; i++) {
        len = upper_corner[i] - lower_corner[i];
        vol *= len;
        fn = len*vol_density + 1;
        npts[i] = (int)ceil(fn);
        spacs[i] = len/((double)(npts[i])-1.0);
        if (apolparm != VNULL) {
            if (apolparm->setgrid) {
                if (apolparm->grid[i] > spacs[i]) {
                    Vnm_print(2, "Vacc_totalSAV:  Warning, your GRID value (%g) is larger than the recommended value (%g)!\n",
                        apolparm->grid[i], spacs[i]);
                }
                spacs[i] = apolparm->grid[i];

            }
        }
    }

    for (x=lower_corner[0]; x<=upper_corner[0]; x=x+spacs[0]) {
        if ( VABS(x - lower_corner[0]) < VSMALL) {
            wx = 0.5;
        } else if ( VABS(x - upper_corner[0]) < VSMALL) {
            wx = 0.5;
        } else {
            wx = 1.0;
        }
        vec[0] = x;
        for (y=lower_corner[1]; y<=upper_corner[1]; y=y+spacs[1]) {
            if ( VABS(y - lower_corner[1]) < VSMALL) {
                wy = 0.5;
            } else if ( VABS(y - upper_corner[1]) < VSMALL) {
                wy = 0.5;
            } else {
                wy = 1.0;
            }
            vec[1] = y;
            for (z=lower_corner[2]; z<=upper_corner[2]; z=z+spacs[2]) {
                if ( VABS(z - lower_corner[2]) < VSMALL) {
                    wz = 0.5;
                } else if ( VABS(z - upper_corner[2]) < VSMALL) {
                    wz = 0.5;
                } else {
                    wz = 1.0;
                }
                vec[2] = z;

                w = wx*wy*wz;

                sav += (w*(1.0-Vacc_ivdwAcc(thee, vec, radius)));

            } /* z loop */
        } /* y loop */
    } /* x loop */

    w  = spacs[0]*spacs[1]*spacs[2];
    sav *= w;

    return sav;
}

int Vacc_wcaEnergyAtom(Vacc *thee, APOLparm *apolparm, Valist *alist,
                                 Vclist *clist, int iatom, double *value) {

    int i;
    int npts[3];
    int pad = 14;

        int xmin, ymin, zmin;
        int xmax, ymax, zmax;

        double sigma6, sigma12;

    double spacs[3], vec[3];
    double w, wx, wy, wz, len, fn, x, y, z, vol;
    double x2,y2,z2,r;
    double vol_density, energy, rho, srad;
    double psig, epsilon, watepsilon, sigma, watsigma, eni, chi;

    double *pos;
    double *lower_corner, *upper_corner;

    Vatom *atom = VNULL;
    VASSERT(apolparm != VNULL);

    energy = 0.0;
    vol = 1.0;
    vol_density = 2.0;

    lower_corner = clist->lower_corner;
    upper_corner = clist->upper_corner;

    atom = Valist_getAtom(alist, iatom);
    pos = Vatom_getPosition(atom);

    /* Note:  these are the original temporary water parameters... they have been
        replaced by entries in a parameter file:
    watsigma = 1.7683;
    watepsilon =  0.1521;
    watepsilon = watepsilon*4.184;
    */

    srad = apolparm->srad;
    rho = apolparm->bconc;
    watsigma = apolparm->watsigma;
    watepsilon = apolparm->watepsilon;
    psig = atom->radius;
    epsilon = atom->epsilon;
    sigma = psig + watsigma;
    epsilon = VSQRT((epsilon * watepsilon));

    /* parameters */
    sigma6 = VPOW(sigma,6);
    sigma12 = VPOW(sigma,12);
    /* OPLS-style radius:  double sigmar = sigma*VPOW(2, (1.0/6.0)); */

    xmin = pos[0] - pad;
    xmax = pos[0] + pad;
    ymin = pos[1] - pad;
    ymax = pos[1] + pad;
    zmin = pos[2] - pad;
    zmax = pos[2] + pad;

    for (i=0; i<3; i++) {
        len = (upper_corner[i] + pad) - (lower_corner[i] - pad);
        vol *= len;
        fn = len*vol_density + 1;
        npts[i] = (int)ceil(fn);
        spacs[i] = 0.5;
        if (apolparm->setgrid) {
            if (apolparm->grid[i] > spacs[i]) {
                Vnm_print(2, "Vacc_totalSAV:  Warning, your GRID value (%g) is larger than the recommended value (%g)!\n",
                    apolparm->grid[i], spacs[i]);
            }
            spacs[i] = apolparm->grid[i];
        }
    }

    for (x=xmin; x<=xmax; x=x+spacs[0]) {
        if ( VABS(x - xmin) < VSMALL) {
            wx = 0.5;
        } else if ( VABS(x - xmax) < VSMALL) {
            wx = 0.5;
        } else {
            wx = 1.0;
        }
        vec[0] = x;
        for (y=ymin; y<=ymax; y=y+spacs[1]) {
            if ( VABS(y - ymin) < VSMALL) {
                wy = 0.5;
            } else if ( VABS(y - ymax) < VSMALL) {
                wy = 0.5;
            } else {
                wy = 1.0;
            }
            vec[1] = y;
            for (z=zmin; z<=zmax; z=z+spacs[2]) {
                if ( VABS(z - zmin) < VSMALL) {
                    wz = 0.5;
                } else if ( VABS(z - zmax) < VSMALL) {
                    wz = 0.5;
                } else {
                    wz = 1.0;
                }
                vec[2] = z;

                w = wx*wy*wz;

                chi = Vacc_ivdwAcc(thee, vec, srad);

                if (VABS(chi) > VSMALL) {

                    x2 = VSQR(vec[0]-pos[0]);
                    y2 = VSQR(vec[1]-pos[1]);
                    z2 = VSQR(vec[2]-pos[2]);
                    r = VSQRT(x2+y2+z2);

                    if (r <= 14 && r >= sigma) {
                        eni = chi*rho*epsilon*(-2.0*sigma6/VPOW(r,6)+sigma12/VPOW(r,12));
                    }else if (r <= 14){
                        eni = -1.0*epsilon*chi*rho;
                    }else{
                        eni = 0.0;
                    }
                }else{
                    eni = 0.0;
                }

                energy += eni*w;

            } /* z loop */
        } /* y loop */
    } /* x loop */

    w  = spacs[0]*spacs[1]*spacs[2];
    energy *= w;

    *value = energy;

    return VRC_SUCCESS;
}

VPUBLIC int Vacc_wcaEnergy(Vacc *acc, APOLparm *apolparm, Valist *alist,
                             Vclist *clist){

    int iatom;
    int rc = 0;

    double energy = 0.0;
    double tenergy = 0.0;
    double rho = apolparm->bconc;

    /* Do a sanity check to make sure that watepsilon and watsigma are set
     * If not, return with an error. */
    if(apolparm->setwat == 0){
        Vnm_print(2,"Vacc_wcaEnergy: Error. No value was set for watsigma and watepsilon.\n");
        return VRC_FAILURE;
    }

    if (VABS(rho) < VSMALL) {
        apolparm->wcaEnergy = tenergy;
        return 1;
    }

    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++){
        rc = Vacc_wcaEnergyAtom(acc, apolparm, alist, clist, iatom, &energy);
        if(rc == 0) return 0;

        tenergy += energy;
    }

    apolparm->wcaEnergy = tenergy;

    return VRC_SUCCESS;

}

VPUBLIC int Vacc_wcaForceAtom(Vacc *thee,
                              APOLparm *apolparm,
                              Vclist *clist,
                              Vatom *atom,
                              double *force
                             ){
    int i,
        si,
        npts[3],
        pad = 14,
        xmin,
        ymin,
        zmin,
        xmax,
        ymax,
        zmax;

    double sigma6,
           sigma12,
           spacs[3],
           vec[3],
           fpt[3],
           w,
           wx,
           wy,
           wz,
           len,
           fn,
           x,
           y,
           z,
           vol,
           x2,
           y2,
           z2,
           r,
           vol_density,
           fo,
           rho,
           srad,
           psig,
           epsilon,
           watepsilon,
           sigma,
           watsigma,
           chi,
           *pos,
           *lower_corner,
           *upper_corner;

    /* Allocate needed variables now that we've asserted required conditions. */
    time_t ts;
    ts = clock();

    VASSERT(apolparm != VNULL);

    /* Do a sanity check to make sure that watepsilon and watsigma are set
     * If not, return with an error. */
    if(apolparm->setwat == 0){
        Vnm_print(2,"Vacc_wcaEnergy: Error. No value was set for watsigma and watepsilon.\n");
        return VRC_FAILURE;
    }

    vol = 1.0;
    vol_density = 2.0;

    lower_corner = clist->lower_corner;
    upper_corner = clist->upper_corner;

    pos = Vatom_getPosition(atom);

    srad = apolparm->srad;
    rho = apolparm->bconc;
    watsigma = apolparm->watsigma;
    watepsilon = apolparm->watepsilon;

    psig = atom->radius;
    epsilon = atom->epsilon;
    sigma = psig + watsigma;
    epsilon = VSQRT((epsilon * watepsilon));

    /* parameters */
    sigma6 = VPOW(sigma,6);
    sigma12 = VPOW(sigma,12);
    /* OPLS-style radius:  double sigmar = sigma*VPOW(2, (1.0/6.0));  */

    for (i=0; i<3; i++) {
        len = (upper_corner[i] + pad) - (lower_corner[i] - pad);
        vol *= len;
        fn = len*vol_density + 1;
        npts[i] = (int)ceil(fn);
        spacs[i] = 0.5;
        force[i] = 0.0;
        if (apolparm->setgrid) {
            if (apolparm->grid[i] > spacs[i]) {
                Vnm_print(2, "Vacc_totalSAV:  Warning, your GRID value (%g) is larger than the recommended value (%g)!\n",
                    apolparm->grid[i], spacs[i]);
            }
            spacs[i] = apolparm->grid[i];
        }
    }

    xmin = pos[0] - pad;
    xmax = pos[0] + pad;
    ymin = pos[1] - pad;
    ymax = pos[1] + pad;
    zmin = pos[2] - pad;
    zmax = pos[2] + pad;

    for (x=xmin; x<=xmax; x=x+spacs[0]) {
        if ( VABS(x - xmin) < VSMALL) {
            wx = 0.5;
        } else if ( VABS(x - xmax) < VSMALL) {
            wx = 0.5;
        } else {
            wx = 1.0;
        }
        vec[0] = x;
        for (y=ymin; y<=ymax; y=y+spacs[1]) {
            if ( VABS(y - ymin) < VSMALL) {
                wy = 0.5;
            } else if ( VABS(y - ymax) < VSMALL) {
                wy = 0.5;
            } else {
                wy = 1.0;
            }
            vec[1] = y;
            for (z=zmin; z<=zmax; z=z+spacs[2]) {
                if ( VABS(z - zmin) < VSMALL) {
                    wz = 0.5;
                } else if ( VABS(z - zmax) < VSMALL) {
                    wz = 0.5;
                } else {
                    wz = 1.0;
                }
                vec[2] = z;

                w = wx*wy*wz;

                chi = Vacc_ivdwAcc(thee, vec, srad);

                if (chi != 0.0) {

                    x2 = VSQR(vec[0]-pos[0]);
                    y2 = VSQR(vec[1]-pos[1]);
                    z2 = VSQR(vec[2]-pos[2]);
                    r = VSQRT(x2+y2+z2);

                    if (r <= 14 && r >= sigma){

                        fo = 12.0*chi*rho*epsilon*(sigma6/VPOW(r,7)-sigma12/VPOW(r,13));

                        fpt[0] = -1.0*(pos[0]-vec[0])*fo/r;
                        fpt[1] = -1.0*(pos[1]-vec[1])*fo/r;
                        fpt[2] = -1.0*(pos[2]-vec[2])*fo/r;

                    }else {
                        for (si=0; si < 3; si++) fpt[si] = 0.0;
                    }
                }else {
                    for (si=0; si < 3; si++) fpt[si] = 0.0;
                }

                for(i=0;i<3;i++){
                    force[i] += (w*fpt[i]);
                }

            } /* z loop */
        } /* y loop */
    } /* x loop */

    w  = spacs[0]*spacs[1]*spacs[2];
    for(i=0;i<3;i++) force[i] *= w;

    return VRC_SUCCESS;
}

