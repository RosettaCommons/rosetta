/**
 *  @file    vgreen.c
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @brief   Class Vgreen methods
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

#include "vgreen.h"

/* Define wrappers for F77 treecode routines */
#ifdef HAVE_TREE
#  define F77TREEPEFORCE VF77_MANGLE(treepeforce, TREEPEFORCE)
#  define F77DIRECT_ENG_FORCE VF77_MANGLE(direct_eng_force, DIRECT_ENG_FORCE)
#  define F77CLEANUP VF77_MANGLE(mycleanup, MYCLEANUP)
#  define F77TREE_COMPP VF77_MANGLE(mytree_compp, MYTREE_COMP)
#  define F77TREE_COMPFP VF77_MANGLE(mytree_compfp, MYTREE_COMPFP)
#  define F77CREATE_TREE VF77_MANGLE(mycreate_tree, MYCREATE_TREE)
#  define F77INITLEVELS VF77_MANGLE(myinitlevels, MYINITLEVELS)
#  define F77SETUP VF77_MANGLE(mysetup, MYSETUP)
#endif  /* ifdef HAVE_TREE */

/* Some constants associated with the tree code */
#ifdef HAVE_TREE
    /**
     * @brief  Lower distance cutoff for electrostatic interactions
     * @ingroup  Vgreen */
#   define FMM_DIST_TOL VSMALL
    /**
     * @brief  Flag for energy and force evaluation:
     *         \li 1 =>  evaluate energy only
     *         \li 2 =>  evaluate energy and force
     * @ingroup  Vgreen */
#   define FMM_IFLAG 2
    /**
     * @brief  Order of multipole expansion
     * @ingroup  Vgreen */
#   define FMM_ORDER 4
    /**
     * @brief  Multipole acceptance criterion
     * @ingroup  Vgreen */
#   define FMM_THETA 0.5
    /**
     * @brief  Maximum number of particles per node
     * @ingroup  Vgreen */
#   define FMM_MAXPARNODE 150
    /**
     * @brief  Switch for oct-tree construction
     * @ingroup  Vgreen */
#   define FMM_SHRINK 1
    /**
     * @brief  Minimum treecode level
     * @ingroup  Vgreen */
#   define FMM_MINLEVEL 50000
    /**
     * @brief  Maximum treecode level
     * @ingroup  Vgreen */
#   define FMM_MAXLEVEL 0
#endif  /* ifdef HAVE_TREE */


/*
 * @brief  Setup treecode internal structures
 * @ingroup  Vgreen
 * @author  Nathan Baker
 * @param  thee  Vgreen object
 * @return  1 if successful, 0 otherwise
 */
VPRIVATE int treesetup(Vgreen *thee);

/*
 * @brief  Clean up treecode internal structures
 * @ingroup  Vgreen
 * @author  Nathan Baker
 * @param  thee  Vgreen object
 * @return  1 if successful, 0 otherwise
 */
VPRIVATE int treecleanup(Vgreen *thee);

/*
 * @brief  Calculate forces or potential
 * @ingroup  Vgreen
 * @author  Nathan Baker
 * @param  thee  Vgreen object
 * @return  1 if successful, 0 otherwise
 */
VPRIVATE int treecalc(Vgreen *thee, double *xtar, double *ytar, double *ztar,
        double *qtar, int numtars, double *tpengtar, double *x, double *y,
        double *z, double *q, int numpars, double *fx, double *fy, double *fz,
        int iflag, int farrdim, int arrdim);

#if !defined(VINLINE_VGREEN)

VPUBLIC Valist* Vgreen_getValist(Vgreen *thee) {

   VASSERT(thee != VNULL);
   return thee->alist;

}

VPUBLIC unsigned long int Vgreen_memChk(Vgreen *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

#endif /* if !defined(VINLINE_VGREEN) */

VPUBLIC Vgreen* Vgreen_ctor(Valist *alist) {

    /* Set up the structure */
    Vgreen *thee = VNULL;
    thee = (Vgreen *)Vmem_malloc(VNULL, 1, sizeof(Vgreen) );
    VASSERT( thee != VNULL);
    VASSERT( Vgreen_ctor2(thee, alist));

    return thee;
}

VPUBLIC int Vgreen_ctor2(Vgreen *thee, Valist *alist) {

    VASSERT( thee != VNULL );

    /* Memory management object */
    thee->vmem = Vmem_ctor("APBS:VGREEN");

    /* Set up the atom list and grid manager */
    if (alist == VNULL) {
        Vnm_print(2,"Vgreen_ctor2: got null pointer to Valist object!\n");
    }

    thee->alist = alist;

    /* Setup FMM tree (if applicable) */
#ifdef HAVE_TREE
    if (!treesetup(thee)) {
        Vnm_print(2, "Vgreen_ctor2:  Error setting up FMM tree!\n");
        return 0;
    }
#endif /* ifdef HAVE_TREE */

    return 1;
}

VPUBLIC void Vgreen_dtor(Vgreen **thee) {
    if ((*thee) != VNULL) {
        Vgreen_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vgreen), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void Vgreen_dtor2(Vgreen *thee) {

#ifdef HAVE_TREE
    treecleanup(thee);
#endif
    Vmem_dtor(&(thee->vmem));

}

VPUBLIC int Vgreen_helmholtz(Vgreen *thee, int npos, double *x, double *y,
  double *z, double *val, double kappa) {

    Vnm_print(2, "Error -- Vgreen_helmholtz not implemented yet!\n");
    return 0;
}

VPUBLIC int Vgreen_helmholtzD(Vgreen *thee, int npos, double *x, double *y,
  double *z, double *gradx, double *grady, double *gradz, double kappa) {

    Vnm_print(2, "Error -- Vgreen_helmholtzD not implemented yet!\n");
    return 0;

}

VPUBLIC int Vgreen_coulomb_direct(Vgreen *thee, int npos, double *x,
        double *y, double *z, double *val) {

    Vatom *atom;
    double *apos, charge, dist, dx, dy, dz, scale;
    double *q, qtemp, fx, fy, fz;
    int iatom, ipos;

    if (thee == VNULL) {
        Vnm_print(2, "Vgreen_coulomb:  Got NULL thee!\n");
        return 0;
    }

    for (ipos=0; ipos<npos; ipos++) val[ipos] = 0.0;

    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) {
        atom = Valist_getAtom(thee->alist, iatom);
        apos = Vatom_getPosition(atom);
        charge = Vatom_getCharge(atom);
        for (ipos=0; ipos<npos; ipos++) {
            dx = apos[0] - x[ipos];
            dy = apos[1] - y[ipos];
            dz = apos[2] - z[ipos];
            dist = VSQRT(VSQR(dx) + VSQR(dy) + VSQR(dz));
            if (dist > VSMALL) val[ipos] += (charge/dist);
        }
    }

    scale = Vunit_ec/(4*Vunit_pi*Vunit_eps0*1.0e-10);
    for (ipos=0; ipos<npos; ipos++) val[ipos] = val[ipos]*scale;

    return 1;
}

VPUBLIC int Vgreen_coulomb(Vgreen *thee, int npos, double *x, double *y,
  double *z, double *val) {

    Vatom *atom;
    double *apos, charge, dist, dx, dy, dz, scale;
    double *q, qtemp, fx, fy, fz;
    int iatom, ipos;

    if (thee == VNULL) {
        Vnm_print(2, "Vgreen_coulomb:  Got NULL thee!\n");
        return 0;
    }

    for (ipos=0; ipos<npos; ipos++) val[ipos] = 0.0;

#ifdef HAVE_TREE

    /* Allocate charge array (if necessary) */
    if (Valist_getNumberAtoms(thee->alist) > 1) {
        if (npos > 1) {
            q = VNULL;
            q = Vmem_malloc(thee->vmem, npos, sizeof(double));
            if (q == VNULL) {
                Vnm_print(2, "Vgreen_coulomb:  Error allocating charge array!\n");
                return 0;
            }
        } else {
            q = &(qtemp);
        }
        for (ipos=0; ipos<npos; ipos++) q[ipos] = 1.0;

        /* Calculate */
        treecalc(thee, x, y, z, q, npos, val, thee->xp, thee->yp, thee->zp,
          thee->qp, thee->np, &fx, &fy, &fz, 1, 1, thee->np);
    } else return Vgreen_coulomb_direct(thee, npos, x, y, z, val);

    /* De-allocate charge array (if necessary) */
    if (npos > 1) Vmem_free(thee->vmem, npos, sizeof(double), (void **)&q);

    scale = Vunit_ec/(4*Vunit_pi*Vunit_eps0*1.0e-10);
    for (ipos=0; ipos<npos; ipos++) val[ipos] = val[ipos]*scale;

    return 1;

#else /* ifdef HAVE_TREE */

    return Vgreen_coulomb_direct(thee, npos, x, y, z, val);

#endif

}

VPUBLIC int Vgreen_coulombD_direct(Vgreen *thee, int npos,
        double *x, double *y, double *z, double *pot, double *gradx,
        double *grady, double *gradz) {

    Vatom *atom;
    double *apos, charge, dist, dist2, idist3, dy, dz, dx, scale;
    double *q, qtemp;
    int iatom, ipos;

    if (thee == VNULL) {
        Vnm_print(2, "Vgreen_coulombD:  Got VNULL thee!\n");
        return 0;
    }

    for (ipos=0; ipos<npos; ipos++) {
        pot[ipos] = 0.0;
        gradx[ipos] = 0.0;
        grady[ipos] = 0.0;
        gradz[ipos] = 0.0;
    }

    for (iatom=0; iatom<Valist_getNumberAtoms(thee->alist); iatom++) {
        atom = Valist_getAtom(thee->alist, iatom);
        apos = Vatom_getPosition(atom);
        charge = Vatom_getCharge(atom);
        for (ipos=0; ipos<npos; ipos++) {
            dx = apos[0] - x[ipos];
            dy = apos[1] - y[ipos];
            dz = apos[2] - z[ipos];
            dist2 = VSQR(dx) + VSQR(dy) + VSQR(dz);
            dist = VSQRT(dist2);
            if (dist > VSMALL) {
                idist3 = 1.0/(dist*dist2);
                gradx[ipos] -= (charge*dx*idist3);
                grady[ipos] -= (charge*dy*idist3);
                gradz[ipos] -= (charge*dz*idist3);
                pot[ipos] += (charge/dist);
            }
        }
    }

    scale = Vunit_ec/(4*VPI*Vunit_eps0*(1.0e-10));
    for (ipos=0; ipos<npos; ipos++) {
        gradx[ipos] = gradx[ipos]*scale;
        grady[ipos] = grady[ipos]*scale;
        gradz[ipos] = gradz[ipos]*scale;
        pot[ipos] = pot[ipos]*scale;
    }

    return 1;
}

VPUBLIC int Vgreen_coulombD(Vgreen *thee, int npos, double *x, double *y,
        double *z, double *pot, double *gradx, double *grady, double *gradz) {

    Vatom *atom;
    double *apos, charge, dist, dist2, idist3, dy, dz, dx, scale;
    double *q, qtemp;
    int iatom, ipos;

    if (thee == VNULL) {
        Vnm_print(2, "Vgreen_coulombD:  Got VNULL thee!\n");
        return 0;
    }

    for (ipos=0; ipos<npos; ipos++) {
        pot[ipos] = 0.0;
        gradx[ipos] = 0.0;
        grady[ipos] = 0.0;
        gradz[ipos] = 0.0;
    }

#ifdef HAVE_TREE

    if (Valist_getNumberAtoms(thee->alist) > 1) {
        if (npos > 1) {
            q = VNULL;
            q = Vmem_malloc(thee->vmem, npos, sizeof(double));
            if (q == VNULL) {
                Vnm_print(2, "Vgreen_coulomb:  Error allocating charge array!\n");
                return 0;
            }
        } else {
            q = &(qtemp);
        }
        for (ipos=0; ipos<npos; ipos++) q[ipos] = 1.0;

        /* Calculate */
        treecalc(thee, x, y, z, q, npos, pot, thee->xp, thee->yp, thee->zp,
                thee->qp, thee->np, gradx, grady, gradz, 2, npos, thee->np);

        /* De-allocate charge array (if necessary) */
        if (npos > 1) Vmem_free(thee->vmem, npos, sizeof(double), (void **)&q);
    } else return Vgreen_coulombD_direct(thee, npos, x, y, z, pot,
            gradx, grady, gradz);

    scale = Vunit_ec/(4*VPI*Vunit_eps0*(1.0e-10));
    for (ipos=0; ipos<npos; ipos++) {
        gradx[ipos] = gradx[ipos]*scale;
        grady[ipos] = grady[ipos]*scale;
        gradz[ipos] = gradz[ipos]*scale;
        pot[ipos] = pot[ipos]*scale;
    }

    return 1;

#else /* ifdef HAVE_TREE */

    return Vgreen_coulombD_direct(thee, npos, x, y, z, pot,
            gradx, grady, gradz);

#endif

}

VPRIVATE int treesetup(Vgreen *thee) {

#ifdef HAVE_TREE

    double dist_tol = FMM_DIST_TOL;
    int iflag = FMM_IFLAG;
    double order = FMM_ORDER;
    int theta = FMM_THETA;
    int shrink = FMM_SHRINK;
    int maxparnode = FMM_MAXPARNODE;
    int minlevel = FMM_MINLEVEL;
    int maxlevel = FMM_MAXLEVEL;
    int level = 0;
    int one = 1;
    Vatom *atom;
    double xyzminmax[6], *pos;
    int i;

    /* Set up particle arrays with atomic coordinates and charges */
    Vnm_print(0, "treesetup:  Initializing FMM particle arrays...\n");
    thee->np = Valist_getNumberAtoms(thee->alist);
    thee->xp = VNULL;
    thee->xp = (double *)Vmem_malloc(thee->vmem, thee->np, sizeof(double));
    if (thee->xp == VNULL) {
        Vnm_print(2, "Vgreen_ctor2:  Failed to allocate %d*sizeof(double)!\n",
          thee->np);
        return 0;
    }
    thee->yp = VNULL;
    thee->yp = (double *)Vmem_malloc(thee->vmem, thee->np, sizeof(double));
    if (thee->yp == VNULL) {
        Vnm_print(2, "Vgreen_ctor2:  Failed to allocate %d*sizeof(double)!\n",
          thee->np);
        return 0;
    }
    thee->zp = VNULL;
    thee->zp = (double *)Vmem_malloc(thee->vmem, thee->np, sizeof(double));
    if (thee->zp == VNULL) {
        Vnm_print(2, "Vgreen_ctor2:  Failed to allocate %d*sizeof(double)!\n",
          thee->np);
        return 0;
    }
    thee->qp = VNULL;
    thee->qp = (double *)Vmem_malloc(thee->vmem, thee->np, sizeof(double));
    if (thee->qp == VNULL) {
        Vnm_print(2, "Vgreen_ctor2:  Failed to allocate %d*sizeof(double)!\n",
          thee->np);
        return 0;
    }
    for (i=0; i<thee->np; i++) {
        atom = Valist_getAtom(thee->alist, i);
        pos = Vatom_getPosition(atom);
        thee->xp[i] = pos[0];
        thee->yp[i] = pos[1];
        thee->zp[i] = pos[2];
        thee->qp[i] = Vatom_getCharge(atom);
    }

    Vnm_print(0, "treesetup:  Setting things up...\n");
    F77SETUP(thee->xp, thee->yp, thee->zp, &(thee->np), &order, &theta, &iflag,
            &dist_tol, xyzminmax, &(thee->np));


    Vnm_print(0, "treesetup:  Initializing levels...\n");
    F77INITLEVELS(&minlevel, &maxlevel);

    Vnm_print(0, "treesetup:  Creating tree...\n");
    F77CREATE_TREE(&one, &(thee->np), thee->xp, thee->yp, thee->zp, thee->qp,
      &shrink, &maxparnode, xyzminmax, &level, &(thee->np));

    return 1;

#else /* ifdef HAVE_TREE */

    Vnm_print(2, "treesetup:  Error!  APBS not linked with treecode!\n");
    return 0;

#endif /* ifdef HAVE_TREE */
}

VPRIVATE int treecleanup(Vgreen *thee) {

#ifdef HAVE_TREE

    Vmem_free(thee->vmem, thee->np, sizeof(double), (void **)&(thee->xp));
    Vmem_free(thee->vmem, thee->np, sizeof(double), (void **)&(thee->yp));
    Vmem_free(thee->vmem, thee->np, sizeof(double), (void **)&(thee->zp));
    Vmem_free(thee->vmem, thee->np, sizeof(double), (void **)&(thee->qp));
    F77CLEANUP();

    return 1;

#else /* ifdef HAVE_TREE */

    Vnm_print(2, "treecleanup:  Error!  APBS not linked with treecode!\n");
    return 0;

#endif /* ifdef HAVE_TREE */
}

VPRIVATE int treecalc(Vgreen *thee, double *xtar, double *ytar, double *ztar,
        double *qtar, int numtars, double *tpengtar, double *x, double *y,
        double *z, double *q, int numpars, double *fx, double *fy, double *fz,
        int iflag, int farrdim, int arrdim) {

#ifdef HAVE_TREE
    int i, level, err, maxlevel, minlevel, one;
    double xyzminmax[6];


    if (iflag != 1) {
        F77TREE_COMPFP(xtar, ytar, ztar, qtar, &numtars, tpengtar, x, y, z, q,
                fx, fy, fz, &numpars, &farrdim, &arrdim);
    } else {
        F77TREE_COMPP(xtar, ytar, ztar, qtar, &numtars, tpengtar, &farrdim, x,
                y, z, q, &numpars, &arrdim);
    }


    return 1;

#else /* ifdef HAVE_TREE */

    Vnm_print(2, "treecalc:  Error!  APBS not linked with treecode!\n");
    return 0;

#endif /* ifdef HAVE_TREE */
}

