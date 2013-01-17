/**
 *  @file    vopot.c
 *  @author  Nathan Baker
 *  @brief   Class Vopot methods
 *  @ingroup Vopot
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

#include "vopot.h"

VEMBED(rcsid="$Id$")

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_ctor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vopot* Vopot_ctor(Vmgrid *mgrid, Vpbe *pbe, Vbcfl bcfl) {

    Vopot *thee = VNULL;

    thee = Vmem_malloc(VNULL, 1, sizeof(Vopot));
    VASSERT(thee != VNULL);
    VASSERT(Vopot_ctor2(thee, mgrid, pbe, bcfl));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_ctor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vopot_ctor2(Vopot *thee, Vmgrid *mgrid, Vpbe *pbe, Vbcfl bcfl) {

    if (thee == VNULL) return 0;
    thee->bcfl = bcfl;
    thee->mgrid = mgrid;
    thee->pbe = pbe;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_dtor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vopot_dtor(Vopot **thee) {

    if ((*thee) != VNULL) {
        Vopot_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vopot), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_dtor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vopot_dtor2(Vopot *thee) { return; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_pot
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
VPUBLIC int Vopot_pot(Vopot *thee, double pt[3], double *value) {

    Vatom *atom;
    int i, iatom;
    double u, T, charge, eps_w, xkappa, dist, size, val, *position;
    Valist *alist;

    VASSERT(thee != VNULL);

    eps_w = Vpbe_getSolventDiel(thee->pbe);
    xkappa = (1.0e10)*Vpbe_getXkappa(thee->pbe);
    T = Vpbe_getTemperature(thee->pbe);
    alist = Vpbe_getValist(thee->pbe);

    u = 0;

    /* See if we're on the mesh */
    if (Vmgrid_value(thee->mgrid, pt, &u)) {

        *value = u;

    } else {

        switch (thee->bcfl) {

            case BCFL_ZERO:
                u = 0;
                break;

            case BCFL_SDH:
                size = (1.0e-10)*Vpbe_getSoluteRadius(thee->pbe);
                position = Vpbe_getSoluteCenter(thee->pbe);
                charge = Vunit_ec*Vpbe_getSoluteCharge(thee->pbe);
                dist = 0;
                for (i=0; i<3; i++)
                  dist += VSQR(position[i] - pt[i]);
                dist = (1.0e-10)*VSQRT(dist);
                val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
                if (xkappa != 0.0)
                  val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
                val = val*Vunit_ec/(Vunit_kb*T);
                u = val;
                break;

            case BCFL_MDH:
                u = 0;
                for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
                    atom = Valist_getAtom(alist, iatom);
                    position = Vatom_getPosition(atom);
                    charge = Vunit_ec*Vatom_getCharge(atom);
                    size = (1e-10)*Vatom_getRadius(atom);
                    dist = 0;
                    for (i=0; i<3; i++)
                      dist += VSQR(position[i] - pt[i]);
                    dist = (1.0e-10)*VSQRT(dist);
                    val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
                    if (xkappa != 0.0)
                      val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
                    val = val*Vunit_ec/(Vunit_kb*T);
                    u = u + val;
                }
                break;

            case BCFL_UNUSED:
                Vnm_print(2, "Vopot_pot:  Invalid bcfl flag (%d)!\n",
                  thee->bcfl);
                return 0;

            case BCFL_FOCUS:
                Vnm_print(2, "Vopot_pot:  Invalid bcfl flag (%d)!\n",
                  thee->bcfl);
                return 0;

            default:
                Vnm_print(2, "Vopot_pot:  Bogus thee->bcfl flag (%d)!\n",
                  thee->bcfl);
                return 0;
                break;
        }

        *value = u;

    }

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_curvature
//
//   Notes:  cflag=0 ==> Reduced Maximal Curvature
//           cflag=1 ==> Mean Curvature (Laplace)
//           cflag=2 ==> Gauss Curvature
//           cflag=3 ==> True Maximal Curvature
//   If we are off the grid, we can still evaluate the Laplacian; assuming, we
//   are away from the molecular surface, it is simply equal to the DH factor.
//
// Authors:  Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vopot_curvature(Vopot *thee, double pt[3], int cflag,
  double *value) {

    Vatom *atom;
    int i, iatom;
    double u, T, charge, eps_w, xkappa, dist, size, val, *position, zkappa2;
    Valist *alist;

    VASSERT(thee != VNULL);

    eps_w = Vpbe_getSolventDiel(thee->pbe);
    xkappa = (1.0e10)*Vpbe_getXkappa(thee->pbe);
    zkappa2 = Vpbe_getZkappa2(thee->pbe);
    T = Vpbe_getTemperature(thee->pbe);
    alist = Vpbe_getValist(thee->pbe);

    u = 0;

    if (Vmgrid_curvature(thee->mgrid, pt, cflag, value)) return 1;
    else if (cflag != 1) {
        Vnm_print(2, "Vopot_curvature:  Off mesh!\n");
        return 1;
    } else {

        switch (thee->bcfl) {

            case BCFL_ZERO:
                u = 0;
                break;

            case BCFL_SDH:
                size = (1.0e-10)*Vpbe_getSoluteRadius(thee->pbe);
                position = Vpbe_getSoluteCenter(thee->pbe);
                charge = Vunit_ec*Vpbe_getSoluteCharge(thee->pbe);
                dist = 0;
                for (i=0; i<3; i++)
                  dist += VSQR(position[i] - pt[i]);
                dist = (1.0e-10)*VSQRT(dist);
                if (xkappa != 0.0)
                  u = zkappa2*(exp(-xkappa*(dist-size))/(1+xkappa*size));
                break;

            case BCFL_MDH:
                u = 0;
                for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
                    atom = Valist_getAtom(alist, iatom);
                    position = Vatom_getPosition(atom);
                    charge = Vunit_ec*Vatom_getCharge(atom);
                    size = (1e-10)*Vatom_getRadius(atom);
                    dist = 0;
                    for (i=0; i<3; i++)
                      dist += VSQR(position[i] - pt[i]);
                    dist = (1.0e-10)*VSQRT(dist);
                    if (xkappa != 0.0)
                      val = zkappa2*(exp(-xkappa*(dist-size))/(1+xkappa*size));
                    u = u + val;
                }
                break;

            case BCFL_UNUSED:
                Vnm_print(2, "Vopot_pot:  Invlid bcfl (%d)!\n", thee->bcfl);
                return 0;

            case BCFL_FOCUS:
                Vnm_print(2, "Vopot_pot:  Invlid bcfl (%d)!\n", thee->bcfl);
                return 0;

            default:
                Vnm_print(2, "Vopot_pot:  Bogus thee->bcfl flag (%d)!\n",
                  thee->bcfl);
                return 0;
                break;
        }

        *value = u;
    }

    return 1;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vopot_gradient
//
// Authors:  Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vopot_gradient(Vopot *thee, double pt[3], double grad[3]) {

    Vatom *atom;
    int iatom;
    double T, charge, eps_w, xkappa, size, val, *position;
    double dx, dy, dz, dist;
    Valist *alist;

    VASSERT(thee != VNULL);

    eps_w = Vpbe_getSolventDiel(thee->pbe);
    xkappa = (1.0e10)*Vpbe_getXkappa(thee->pbe);
    T = Vpbe_getTemperature(thee->pbe);
    alist = Vpbe_getValist(thee->pbe);


    if (!Vmgrid_gradient(thee->mgrid, pt, grad)) {

        switch (thee->bcfl) {

            case BCFL_ZERO:
                grad[0] = 0.0;
                grad[1] = 0.0;
                grad[2] = 0.0;
                break;

            case BCFL_SDH:
                grad[0] = 0.0;
                grad[1] = 0.0;
                grad[2] = 0.0;
                size = (1.0e-10)*Vpbe_getSoluteRadius(thee->pbe);
                position = Vpbe_getSoluteCenter(thee->pbe);
                charge = Vunit_ec*Vpbe_getSoluteCharge(thee->pbe);
                dx = position[0] - pt[0];
                dy = position[1] - pt[1];
                dz = position[2] - pt[2];
                dist = VSQR(dx) + VSQR(dy) + VSQR(dz);
                dist = (1.0e-10)*VSQRT(dist);
                val = (charge)/(4*VPI*Vunit_eps0*eps_w);
                if (xkappa != 0.0)
                  val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
                val = val*Vunit_ec/(Vunit_kb*T);
                grad[0] = val*dx/dist*(-1.0/dist/dist + xkappa/dist);
                grad[1] = val*dy/dist*(-1.0/dist/dist + xkappa/dist);
                grad[2] = val*dz/dist*(-1.0/dist/dist + xkappa/dist);
                break;

            case BCFL_MDH:
                grad[0] = 0.0;
                grad[1] = 0.0;
                grad[2] = 0.0;
                for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
                    atom = Valist_getAtom(alist, iatom);
                    position = Vatom_getPosition(atom);
                    charge = Vunit_ec*Vatom_getCharge(atom);
                    size = (1e-10)*Vatom_getRadius(atom);
                    dx = position[0] - pt[0];
                    dy = position[1] - pt[1];
                    dz = position[2] - pt[2];
                    dist = VSQR(dx) + VSQR(dy) + VSQR(dz);
                    dist = (1.0e-10)*VSQRT(dist);
                    val = (charge)/(4*VPI*Vunit_eps0*eps_w);
                    if (xkappa != 0.0)
                      val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
                    val = val*Vunit_ec/(Vunit_kb*T);
                    grad[0] += (val*dx/dist*(-1.0/dist/dist + xkappa/dist));
                    grad[1] += (val*dy/dist*(-1.0/dist/dist + xkappa/dist));
                    grad[2] += (val*dz/dist*(-1.0/dist/dist + xkappa/dist));
                }
                break;

            case BCFL_UNUSED:
                Vnm_print(2, "Vopot:  Invalid bcfl (%d)!\n", thee->bcfl);
                return 0;

            case BCFL_FOCUS:
                Vnm_print(2, "Vopot:  Invalid bcfl (%d)!\n", thee->bcfl);
                return 0;

            default:
                Vnm_print(2, "Vopot_pot:  Bogus thee->bcfl flag (%d)!\n",
                  thee->bcfl);
                return 0;
                break;
        }

        return 1;
    }

    return 1;

}

