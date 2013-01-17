/**
 *  @file    vfetk.c
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @brief   Class Vfetk methods
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

#include "vfetk.h"

/* Define the macro DONEUMANN to run with all-Neumann boundary conditions.
 * Set this macro at your own risk! */
/* #define DONEUMANN 1 */

/*
 * @brief  Calculuate the contribution to the charge-potential energy from one
 * atom
 * @ingroup Vfetk
 * @author  Nathan Baker
 * @param  thee  current Vfetk object
 * @param  iatom  current atom index
 * @param  color  simplex subset (partition) under consideration
 * @param  sol  current solution
 * @returns Per-atom energy
 */
VPRIVATE double Vfetk_qfEnergyAtom(
        Vfetk *thee,
        int iatom,
        int color,
        double *sol
        );

/*
 * @brief  Container for local variables
 * @ingroup  Vfetk
 * @bug  Not thread-safe
 */
VPRIVATE Vfetk_LocalVar var;

/*
 * @brief  MCSF-format cube mesh (all Dirichlet)
 * @ingroup Vfetk
 * @author  Based on mesh by Mike Holst
 */
VPRIVATE char *diriCubeString =
"mcsf_begin=1;\n\
\n\
dim=3;\n\
dimii=3;\n\
vertices=8;\n\
simplices=6;\n\
\n\
vert=[\n\
0 0 -0.5 -0.5 -0.5\n\
1 0  0.5 -0.5 -0.5\n\
2 0 -0.5  0.5 -0.5\n\
3 0  0.5  0.5 -0.5\n\
4 0 -0.5 -0.5  0.5\n\
5 0  0.5 -0.5  0.5\n\
6 0 -0.5  0.5  0.5\n\
7 0  0.5  0.5  0.5\n\
];\n\
\n\
simp=[\n\
0 0 0 0 1 0 1 0 5 1 2\n\
1 0 0 0 1 1 0 0 5 2 4\n\
2 0 0 0 1 0 1 1 5 3 2\n\
3 0 0 0 1 0 1 3 5 7 2\n\
4 0 0 1 1 0 0 2 5 7 6\n\
5 0 0 1 1 0 0 2 5 6 4\n\
];\n\
\n\
mcsf_end=1;\n\
\n\
";

/*
 * @brief  MCSF-format cube mesh (all Neumann)
 * @ingroup Vfetk
 * @author  Based on mesh by Mike Holst
 */
VPRIVATE char *neumCubeString =
"mcsf_begin=1;\n\
\n\
dim=3;\n\
dimii=3;\n\
vertices=8;\n\
simplices=6;\n\
\n\
vert=[\n\
0 0 -0.5 -0.5 -0.5\n\
1 0  0.5 -0.5 -0.5\n\
2 0 -0.5  0.5 -0.5\n\
3 0  0.5  0.5 -0.5\n\
4 0 -0.5 -0.5  0.5\n\
5 0  0.5 -0.5  0.5\n\
6 0 -0.5  0.5  0.5\n\
7 0  0.5  0.5  0.5\n\
];\n\
\n\
simp=[\n\
0 0 0 0 2 0 2 0 5 1 2\n\
1 0 0 0 2 2 0 0 5 2 4\n\
2 0 0 0 2 0 2 1 5 3 2\n\
3 0 0 0 2 0 2 3 5 7 2\n\
4 0 0 2 2 0 0 2 5 7 6\n\
5 0 0 2 2 0 0 2 5 6 4\n\
];\n\
\n\
mcsf_end=1;\n\
\n\
";

/*
 * @brief  Return the smoothed value of the dielectric coefficient at the
 * current point using a fast, chart-based method
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @returns  Value of dielectric coefficient
 * @bug  Not thread-safe
 */
VPRIVATE double diel();

/*
 * @brief  Return the smoothed value of the ion accessibility at the
 * current point using a fast, chart-based method
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @returns  Value of mobile ion coefficient
 * @bug  Not thread-safe
 */
VPRIVATE double ionacc();

/*
 * @brief  Smooths a mesh-based coefficient with a simple harmonic function
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  meth  Method for smoothing
 *   \li  0 ==> arithmetic mean (gives bad results)
 *   \li  1 ==> geometric mean
 * @param  nverts  Number of vertices
 * @param  dist  distance from point to each vertex
 * @param  coeff  coefficient value at each vertex
 * @note  Thread-safe
 * @return smoothed value of coefficieent at point of interest */
VPRIVATE double smooth(
        int nverts,
        double dist[VAPBS_NVS],
        double coeff[VAPBS_NVS],
        int meth
        );


/*
 * @brief  Return the analytical multi-sphere Debye-Huckel approximation (in
 * kT/e) at the specified point
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  pbe  Vpbe object
 * @param  d  Dimension of x
 * @param  x  Coordinates of point of interest (in &Aring;)
 * @note  Thread-safe
 * @returns  Multi-sphere Debye-Huckel potential in kT/e
 */
VPRIVATE double debye_U(
        Vpbe *pbe,
        int d,
        double x[]
        );

/*
 * @brief  Return the difference between the analytical multi-sphere
 * Debye-Huckel approximation and Coulomb's law (in kT/e) at the specified
 * point
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  pbe  Vpbe object
 * @param  d  Dimension of x
 * @param  x  Coordinates of point of interest (in &Aring;)
 * @note  Thread-safe
 * @returns  Multi-sphere Debye-Huckel potential in kT/e */
VPRIVATE double debye_Udiff(
        Vpbe *pbe,
        int d,
        double x[]
        );

/*
 * @brief  Calculate the Coulomb's
 * Debye-Huckel approximation and Coulomb's law (in kT/e) at the specified
 * point
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  pbe  Vpbe object
 * @param  d  Dimension of x
 * @param  x  Coordinates of point of interest (in &Aring;)
 * @param  eps  Dielectric constant
 * @param  U  Set to potential (in kT/e)
 * @param  dU  Set to potential gradient (in kT/e/&Aring;)
 * @param  d2U  Set to Laplacian of potential (in \f$kT e^{-1} \AA^{-2}\f$)
 * @returns  Multi-sphere Debye-Huckel potential in kT/e */
VPRIVATE void coulomb(
        Vpbe *pbe,
        int d,
        double x[],
        double eps,
        double *U,
        double dU[],
        double *d2U
        );

/*
 * @brief  2D linear master simplex information generator
 * @ingroup  Vfetk
 * @author  Mike Holst
 * @param dimIS  dunno
 * @param ndof  dunno
 * @param dof  dunno
 * @param c  dunno
 * @param cx  dunno
 * @note  Trust in Mike */
VPRIVATE void init_2DP1(
        int dimIS[],
        int *ndof,
        int dof[],
        double c[][VMAXP],
        double cx[][VMAXP],
        double cy[][VMAXP],
        double cz[][VMAXP]
        );

/*
 * @brief  3D linear master simplex information generator
 * @ingroup  Vfetk
 * @author  Mike Holst
 * @param  dimIS  dunno
 * @param ndof dunno
 * @param dof dunno
 * @param c dunno
 * @param cx dunno
 * @param cy dunno
 * @param cz dunno
 * @note  Trust in Mike */
VPRIVATE void init_3DP1(
        int dimIS[],
        int *ndof,
        int dof[],
        double c[][VMAXP],
        double cx[][VMAXP],
        double cy[][VMAXP],
        double cz[][VMAXP]
        );

/*
 * @brief  Setup coefficients of polynomials from integer table data
 * @ingroup  Vfetk
 * @author  Mike Holst
 * @param numP  dunno
 * @param c  dunno
 * @param cx  dunno
 * @param cy  dunno
 * @param cz  dunno
 * @param ic  dunno
 * @param icx  dunno
 * @param icy  dunno
 * @param icz  dunno
 * @note  Trust in Mike */
VPRIVATE void setCoef(
        int numP,
        double c[][VMAXP],
        double cx[][VMAXP],
        double cy[][VMAXP],
        double cz[][VMAXP],
        int ic[][VMAXP],
        int icx[][VMAXP],
        int icy[][VMAXP],
        int icz[][VMAXP]
        );

/*
 * @brief  Evaluate a collection of at most cubic polynomials at a
 * specified point in at most R^3.
 * @ingroup  Vfetk
 * @author  Mike Holst
 * @param numP  the number of polynomials to evaluate
 * @param p  the results of the evaluation
 * @param c   the coefficients of each polynomial
 * @param xv  the point (x,y,z) to evaluate the polynomials.
 * @note  Mike says:
 * <pre>
 *  Note that "VMAXP" must be >= 19 for cubic polynomials.
 *  The polynomials are build from the coefficients c[][] as
 *  follows.  To build polynomial "k", fix k and set:
 *
 *  c0=c[k][0], c1=c[k][1], .... , cp=c[k][p]
 *
 *  Then evaluate as:
 *
 *  p3(x,y,z) = c0 + c1*x + c2*y + c3*z
 *            + c4*x*x + c5*y*y + c6*z*z + c7*x*y + c8*x*z + c9*y*z
 *            + c10*x*x*x + c11*y*y*y + c12*z*z*z
 *            + c13*x*x*y + c14*x*x*z + c15*x*y*y
 *            + c16*y*y*z + c17*x*z*z + c18*y*z*z
 * </pre>
 */
VPRIVATE void polyEval(
        int numP,
        double p[],
        double c[][VMAXP],
        double xv[]
        );

/*
 * @brief  I have no clue what this variable does, but we need it to initialize
 * the simplices
 * @ingroup  Vfetk
 * @author  Mike Holst */
VPRIVATE int dim_2DP1 = 3;

/*
 * @brief  I have no clue what these variable do, but we need it to initialize
 * the simplices
 * @ingroup  Vfetk
 * @author  Mike Holst
 * @note  Mike says:
 * <pre>
 *  2D-P1 Basis:
 *
 *  p1(x,y) = c0 + c1*x + c2*y
 *
 *  Lagrange Point    Lagrange Basis Function Definition
 *  --------------    ----------------------------------
 *  (0, 0)            p[0](x,y) = 1 - x - y
 *  (1, 0)            p[1](x,y) = x
 *  (0, 1)            p[2](x,y) = y
 *  </pre>
 */
VPRIVATE int lgr_2DP1[3][VMAXP] = {
/*c0  c1  c2  c3
* ---------------------------------------------------------- */
/* 1   x   y   z
* ---------------------------------------------------------- */
{  2, -2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }
};
VPRIVATE int lgr_2DP1x[3][VMAXP] = {
/*c0 ---------------------------------------------------------------------- */
/* 1 ---------------------------------------------------------------------- */
{ -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }
};
VPRIVATE int lgr_2DP1y[3][VMAXP] = {
/*c0 ---------------------------------------------------------------------- */
/* 1 ---------------------------------------------------------------------- */
{ -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }
};
VPRIVATE int lgr_2DP1z[3][VMAXP] = {
/*c0 ---------------------------------------------------------------------- */
/* 1 ---------------------------------------------------------------------- */
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }
};


/*
 * @brief  I have no clue what these variable do, but we need it to initialize
 * the simplices
 * @ingroup  Vfetk
 * @author  Mike Holst
 * @note  Mike says:
 * <pre>
 * 3D-P1 Basis:
 *
 * p1(x,y,z) = c0 + c1*x + c2*y + c3*z
 *
 * Lagrange Point    Lagrange Basis Function Definition
 * --------------    ----------------------------------
 * (0, 0, 0)         p[0](x,y,z) = 1 - x - y - z
 * (1, 0, 0)         p[1](x,y,z) = x
 * (0, 1, 0)         p[2](x,y,z) = y
 * (0, 0, 1)         p[3](x,y,z) = z
 * </pre>
 */
VPRIVATE int dim_3DP1 = VAPBS_NVS;
VPRIVATE int lgr_3DP1[VAPBS_NVS][VMAXP] = {
/*c0  c1  c2  c3 ---------------------------------------------------------- */
/* 1   x   y   z ---------------------------------------------------------- */
{  2, -2, -2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }
};
VPRIVATE int lgr_3DP1x[VAPBS_NVS][VMAXP] = {
/*c0 ---------------------------------------------------------------------- */
/* 1 ---------------------------------------------------------------------- */
{ -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }
};
VPRIVATE int lgr_3DP1y[VAPBS_NVS][VMAXP] = {
/*c0 ---------------------------------------------------------------------- */
/* 1 ---------------------------------------------------------------------- */
{ -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 }
};
VPRIVATE int lgr_3DP1z[VAPBS_NVS][VMAXP] = {
/*c0 ---------------------------------------------------------------------- */
/* 1 ---------------------------------------------------------------------- */
{ -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
{  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0 },
};

/*
 * @brief  Another Holst variable
 * @ingroup  Vfetk
 * @author  Mike Holst
 * @note  Mike says: 1 = linear, 2 = quadratic */
VPRIVATE const int P_DEG=1;

/*
 * @brief  Another Holst variable
 * @ingroup  Vfetk
 * @author  Mike Holst */
VPRIVATE int numP;
VPRIVATE double c[VMAXP][VMAXP];
VPRIVATE double cx[VMAXP][VMAXP];
VPRIVATE double cy[VMAXP][VMAXP];
VPRIVATE double cz[VMAXP][VMAXP];

#if !defined(VINLINE_VFETK)

VPUBLIC Gem* Vfetk_getGem(Vfetk *thee) {

   VASSERT(thee != VNULL);
   return thee->gm;

}

VPUBLIC AM* Vfetk_getAM(Vfetk *thee) {

   VASSERT(thee != VNULL);
   return thee->am;
}

VPUBLIC Vpbe* Vfetk_getVpbe(Vfetk *thee) {

   VASSERT(thee != VNULL);
   return thee->pbe;

}

VPUBLIC Vcsm* Vfetk_getVcsm(Vfetk *thee) {

   VASSERT(thee != VNULL);
   return thee->csm;

}

VPUBLIC int Vfetk_getAtomColor(Vfetk *thee,
                               int iatom
                              ) {

    int natoms;

    VASSERT(thee != VNULL);

    natoms = Valist_getNumberAtoms(Vpbe_getValist(thee->pbe));
    VASSERT(iatom < natoms);

    return Vatom_getPartID(Valist_getAtom(Vpbe_getValist(thee->pbe), iatom));
}
#endif /* if !defined(VINLINE_VFETK) */

VPUBLIC Vfetk* Vfetk_ctor(Vpbe *pbe,
                          Vhal_PBEType type
                         ) {

    /* Set up the structure */
    Vfetk *thee = VNULL;
    thee = (Vfetk*)Vmem_malloc(VNULL, 1, sizeof(Vfetk) );
    VASSERT(thee != VNULL);
    VASSERT(Vfetk_ctor2(thee, pbe, type));

    return thee;
}

VPUBLIC int Vfetk_ctor2(Vfetk *thee,
                        Vpbe *pbe,
                        Vhal_PBEType type
                       ) {

    int i;
    double center[VAPBS_DIM];

    /* Make sure things have been properly initialized & store them */
    VASSERT(pbe != VNULL);
    thee->pbe = pbe;
    VASSERT(pbe->alist != VNULL);
    VASSERT(pbe->acc != VNULL);

    /* Store PBE type */
    thee->type = type;

    /* Set up memory management object */
    thee->vmem = Vmem_ctor("APBS::VFETK");

    /* Set up FEtk objects */
    Vnm_print(0, "Vfetk_ctor2:  Constructing PDE...\n");
    thee->pde = Vfetk_PDE_ctor(thee);
    Vnm_print(0, "Vfetk_ctor2:  Constructing Gem...\n");
    thee->gm = Gem_ctor(thee->vmem, thee->pde);
    Vnm_print(0, "Vfetk_ctor2:  Constructing Aprx...\n");
    thee->aprx = Aprx_ctor(thee->vmem, thee->gm, thee->pde);
    Vnm_print(0, "Vfetk_ctor2:  Constructing Aprx...\n");
    thee->am = AM_ctor(thee->vmem, thee->aprx);

    /* Reset refinement level */
    thee->level = 0;

    /* Set default solver variables */
    thee->lkey = VLT_MG;
    thee->lmax = 1000000;
    thee->ltol = 1e-5;
    thee->lprec = VPT_MG;
    thee->nkey = VNT_NEW;
    thee->nmax = 1000000;
    thee->ntol = 1e-5;
    thee->gues = VGT_ZERO;
    thee->pjac = -1;

    /* Store local copy of myself */
    var.fetk = thee;
    var.initGreen = 0;

    /* Set up the external Gem subdivision hook */
    Gem_setExternalUpdateFunction(thee->gm, Vfetk_externalUpdateFunction);

    /* Set up ion-related variables */
    var.zkappa2 = Vpbe_getZkappa2(var.fetk->pbe);
    var.ionstr = Vpbe_getBulkIonicStrength(var.fetk->pbe);
    if (var.ionstr > 0.0) var.zks2 = 0.5*var.zkappa2/var.ionstr;
    else var.zks2 = 0.0;
    Vpbe_getIons(var.fetk->pbe, &(var.nion), var.ionConc, var.ionRadii,
      var.ionQ);
    for (i=0; i<var.nion; i++) {
        var.ionConc[i] = var.zks2 * var.ionConc[i] * var.ionQ[i];
    }

    /* Set uninitialized objects to NULL */
    thee->pbeparm = VNULL;
    thee->feparm = VNULL;
    thee->csm = VNULL;

    return 1;
}

VPUBLIC void Vfetk_setParameters(Vfetk *thee,
                                 PBEparm *pbeparm,
                                 FEMparm *feparm
                                ) {

    VASSERT(thee != VNULL);
    thee->feparm = feparm;
    thee->pbeparm = pbeparm;
}

VPUBLIC void Vfetk_dtor(Vfetk **thee) {
    if ((*thee) != VNULL) {
        Vfetk_dtor2(*thee);
        //Vmem_free(VNULL, 1, sizeof(Vfetk), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void Vfetk_dtor2(Vfetk *thee) {
    Vcsm_dtor(&(thee->csm));
    AM_dtor(&(thee->am));
    Aprx_dtor(&(thee->aprx));
    Vfetk_PDE_dtor(&(thee->pde));
    Vmem_dtor(&(thee->vmem));
}

VPUBLIC double* Vfetk_getSolution(Vfetk *thee,
                                  int *length
                                 ) {

   int i;
   double *solution,
          *theAnswer;
   AM *am;

   VASSERT(thee != VNULL);

   /* Get the AM object */
   am = thee->am;
   /* Copy the solution into the w0 vector */
   Bvec_copy(am->w0, am->u);
   /* Add the Dirichlet conditions */
   Bvec_axpy(am->w0, am->ud, 1.);
   /* Get the data from the Bvec */
   solution = Bvec_addr(am->w0);
   /* Get the length of the data from the Bvec */
   *length = Bvec_numRT(am->w0);
   /* Make sure that we got scalar data (only one block) for the solution
    * to the FETK */
   VASSERT(1 == Bvec_numB(am->w0));
   /* Allocate space for the returned vector and copy the solution into it */
   theAnswer = VNULL;
   theAnswer = (double*)Vmem_malloc(VNULL, *length, sizeof(double));
   VASSERT(theAnswer != VNULL);
   for (i=0; i<(*length); i++) theAnswer[i] = solution[i];

   return theAnswer;
}


/**
 * @brief   Return the total electrostatic energy
 *
 *  Using the solution at the finest mesh level, get the electrostatic energy
 *  using the free energy functional for the Poisson-Boltzmann equation
 *  without removing any self-interaction terms (i.e., removing the reference
 *  state of isolated charges present in an infinite dielectric continuum with
 *  the same relative permittivity as the interior of the protein) and return
 *  the result in units of \f$k_B T\f$.  The argument color allows the user to
 *  control the partition on which this energy is calculated; if (color == -1)
 *  no restrictions are used.  The solution is obtained from the finest level
 *  of the passed AM object, but atomic data from the Vfetk object is used to
 *  calculate the energy.
 *
 * @ingroup Vfetk
 * @author  Nathan Baker
 * @return Total electrostatic energy in units of \f$k_B T\f$.
 */
VPUBLIC double Vfetk_energy(Vfetk *thee, /**< THe Vfetk object */
                            int color,  /**< Partition restriction for energy calculation; if
                                            non-negative, energy calculation is restricted to the
                                            specified partition (indexed by simplex and atom colors
                                            */
                            int nonlin  /**< If 1, the NPBE energy functional is used; otherwise,
                                            the LPBE energy functional is used. If -2, SMPBE is used. */
                           ) {

    double totEnergy = 0.0, /**< Total energy */
           qfEnergy = 0.0,  /**< @todo document */
           dqmEnergy = 0.0; /**< @todo docuemnt */

    VASSERT(thee != VNULL);

    if (nonlin && (Vpbe_getBulkIonicStrength(thee->pbe) > 0.)) {
        Vnm_print(0, "Vfetk_energy:  calculating full PBE energy\n");
        Vnm_print(0, "Vfetk_energy:  bulk ionic strength = %g M\n",
          Vpbe_getBulkIonicStrength(thee->pbe));
        dqmEnergy = Vfetk_dqmEnergy(thee, color);
        Vnm_print(0, "Vfetk_energy:  dqmEnergy = %g kT\n", dqmEnergy);
        qfEnergy = Vfetk_qfEnergy(thee, color);
        Vnm_print(0, "Vfetk_energy:  qfEnergy = %g kT\n", qfEnergy);

        totEnergy = qfEnergy - dqmEnergy;
    } else {
        Vnm_print(0, "Vfetk_energy:  calculating only q-phi energy\n");
        dqmEnergy = Vfetk_dqmEnergy(thee, color);
        Vnm_print(0, "Vfetk_energy:  dqmEnergy = %g kT (NOT USED)\n", dqmEnergy);
        qfEnergy = Vfetk_qfEnergy(thee, color);
        Vnm_print(0, "Vfetk_energy:  qfEnergy = %g kT\n", qfEnergy);
        totEnergy = 0.5*qfEnergy;
    }

    return totEnergy;

}


VPUBLIC double Vfetk_qfEnergy(Vfetk *thee,
                              int color
                             ) {

    double *sol,
           energy = 0.0;
    int nsol,
        iatom,
        natoms;
    AM *am;

    VASSERT(thee != VNULL);
    am = thee->am;

    /* Get the finest level solution */
    sol= VNULL;
    sol = Vfetk_getSolution(thee, &nsol);
    VASSERT(sol != VNULL);

    /* Make sure the number of entries in the solution array matches the
     * number of vertices currently in the mesh */
    if (nsol != Gem_numVV(thee->gm)) {
       Vnm_print(2, "Vfetk_qfEnergy: Number of unknowns in solution does not match\n");
       Vnm_print(2, "Vfetk_qfEnergy: number of vertices in mesh!!!  Bailing out!\n");
       VASSERT(0);
    }

    /* Now we do the sum over atoms... */
    natoms = Valist_getNumberAtoms(thee->pbe->alist);
    for (iatom=0; iatom<natoms; iatom++) {

        energy = energy + Vfetk_qfEnergyAtom(thee, iatom, color, sol);

    } /* end for iatom */

    /* Destroy the finest level solution */
    Vmem_free(VNULL, nsol, sizeof(double), (void **)&sol);

    /* Return the energy */
    return energy;
}

VPRIVATE double Vfetk_qfEnergyAtom(
        Vfetk *thee,
        int iatom,
        int color,
        double *sol) {

    Vatom *atom;
    double charge,
           phi[VAPBS_NVS],
           phix[VAPBS_NVS][3],
           *position,
           uval,
           energy = 0.0;
    int isimp,
        nsimps,
        icolor,
        ivert,
        usingColor;
    SS *simp;


    /* Get atom information */
    atom = Valist_getAtom(thee->pbe->alist, iatom);
    icolor = Vfetk_getAtomColor(thee, iatom);
    charge = Vatom_getCharge(atom);
    position = Vatom_getPosition(atom);

    /* Find out if we're using colors */
    usingColor = (color >= 0);

    if (usingColor && (icolor<0)) {
        Vnm_print(2, "Vfetk_qfEnergy: Atom colors not set!\n");
        VASSERT(0);
    }

    /* Check if this atom belongs to the specified partition */
    if ((icolor==color) || (!usingColor)) {
        /* Loop over the simps associated with this atom */
        nsimps =  Vcsm_getNumberSimplices(thee->csm, iatom);

        /* Get the first simp of the correct color; we can use just one
         * simplex for energy evaluations, but not for force
         * evaluations */
        for (isimp=0; isimp<nsimps; isimp++) {

            simp = Vcsm_getSimplex(thee->csm, isimp, iatom);

            /* If we've asked for a particular partition AND if the atom
             * is our partition, then compute the energy */
            if ((SS_chart(simp)==color)||(color<0)) {
                /* Get the value of each basis function evaluated at this
                 * point */
                Gem_pointInSimplexVal(thee->gm, simp, position, phi, phix);
                for (ivert=0; ivert<SS_dimVV(simp); ivert++) {
                    uval = sol[VV_id(SS_vertex(simp,ivert))];
                    energy += (charge*phi[ivert]*uval);
                } /* end for ivert */
                /* We only use one simplex of the appropriate color for
                 * energy calculations, so break here */
                break;
            } /* endif (color) */
        } /* end for isimp */
    }

   return energy;
}


VPUBLIC double Vfetk_dqmEnergy(Vfetk *thee,
                               int color) {

    return AM_evalJ(thee->am);

}

VPUBLIC void Vfetk_setAtomColors(Vfetk *thee) {

    SS *simp;
    Vatom *atom;
    int i,
        natoms;

    VASSERT(thee != VNULL);

    natoms = Valist_getNumberAtoms(thee->pbe->alist);
    for (i=0; i<natoms; i++) {
        atom = Valist_getAtom(thee->pbe->alist, i);
        simp = Vcsm_getSimplex(thee->csm, 0, i);
        Vatom_setPartID(atom, SS_chart(simp));
    }

}

VPUBLIC unsigned long int Vfetk_memChk(Vfetk *thee) {

    int memUse = 0;

    if (thee == VNULL) return 0;

    memUse = memUse + sizeof(Vfetk);
    memUse = memUse + Vcsm_memChk(thee->csm);

    return memUse;
}

/**
 * Generates a new cube mesh within the provided Vfetk object based on the
 * specified mesh type.  Creates a new copy of the mesh based on the global
 * variables at the top of the file and the mesh type, then recenters the
 * mesh based on the center and length variables provided to the function.
 */
VPUBLIC Vrc_Codes Vfetk_genCube(Vfetk *thee, /**< Vfetk object */
                                double center[3], /**< Center for mesh, which the new mesh will adjust to */
                                double length[3], /**< Mesh lengths, which the new mesh will adjust to */
                                Vfetk_MeshLoad meshType /**< Mesh boundary conditions */
                               ) {

    VASSERT(thee != VNULL);

    AM *am = VNULL; /* @todo - no idea what this is */
    Gem *gm = VNULL; /* Geometry manager */

    int skey = 0,  /* Simplex format */
        bufsize = 0, /* Buffer size */
        i, /* Loop counter */
        j; /* Loop counter */
    char *key = "r",  /* Read */
         *iodev = "BUFF",  /* Buffer */
         *iofmt = "ASC",  /* ASCII */
         *iohost = "localhost",  /* localhost (dummy) */
         *iofile = "0",  /*< socket 0 (dummy) */
         buf[VMAX_BUFSIZE]; /* Socket buffer */
    Vio *sock = VNULL; /* Socket object */
    VV *vx = VNULL; /* @todo - no idea what this is */
    double x;

    am = thee->am;
    VASSERT(am != VNULL);
    gm = thee->gm;
    VASSERT(gm != VNULL);

    /* @note This code is based on Gem_makeCube by Mike Holst */
    /* Write mesh string to buffer and read back */
    switch (meshType) {
        case VML_DIRICUBE:
            /* Create a new copy of the DIRICUBE mesh (see globals higher in this file) */
            bufsize = strlen(diriCubeString);
            VASSERT( bufsize <= VMAX_BUFSIZE );
            strncpy(buf, diriCubeString, VMAX_BUFSIZE);
            break;
        case VML_NEUMCUBE:
            /* Create a new copy of the NEUMCUBE mesh (see globals higher in this file) */
            bufsize = strlen(neumCubeString);
            Vnm_print(2, "Vfetk_genCube:  WARNING!  USING EXPERIMENTAL NEUMANN BOUNDARY CONDITIONS!\n");
            VASSERT( bufsize <= VMAX_BUFSIZE );
            strncpy(buf, neumCubeString, VMAX_BUFSIZE);
            break;
        case VML_EXTERNAL:
            Vnm_print(2, "Vfetk_genCube:  Got request for external mesh!\n");
            Vnm_print(2, "Vfetk_genCube:  How did we get here?\n");
            return VRC_FAILURE;
        default:
            Vnm_print(2, "Vfetk_genCube:  Unknown mesh type (%d)\n", meshType);
            return VRC_FAILURE;
    }

    VASSERT( VNULL != (sock=Vio_socketOpen(key,iodev,iofmt,iohost,iofile)) ); /* Open socket */
    Vio_bufTake(sock, buf, bufsize); /* Initialize internal buffer for socket */
    AM_read(am, skey, sock); /* Take the initial mesh from the socket and load
                                into internal AM data structure with simplex
                                format */
    Vio_connectFree(sock); /* Purge output buffers */
    Vio_bufGive(sock); /* Get pointer to output buffer?  No assignment of return value... */
    Vio_dtor(&sock); /* Destroy output buffer */

    /* @todo - could the following be done in a single pass? - PCE */
    /* Scale (unit) cube - for each vertex, set the new coordinates of that
       vertex based on the vertex length */
    for (i=0; i<Gem_numVV(gm); i++) {
        vx = Gem_VV(gm, i);
        for (j=0; j<3; j++) {
            x = VV_coord(vx, j);
            x *= length[j];
            VV_setCoord(vx, j, x);
        }
    }

    /* Add new center - for each vertex, set a new center for the vertex */
    for (i=0; i<Gem_numVV(gm); i++) {
        vx = Gem_VV(gm, i);
        for (j=0; j<3; j++) {
            x = VV_coord(vx, j);
            x += center[j];
            VV_setCoord(vx, j, x);
        }
    }

    return VRC_SUCCESS;
}

/**
 * If we have an external mesh, load that external mesh from the provided
 * socket.  If we specify a non-external mesh type, we generate a new mesh
 * cube based on templates.  We then create and store a new Vcsm object in
 * our Vfetk structure, which will carry the mesh data.
 */
VPUBLIC Vrc_Codes Vfetk_loadMesh(Vfetk *thee, /* Vfetk object to load into */
                                 double center[3], /* Center for mesh (if constructed) */
                                 double length[3], /* Mesh lengths (if constructed) */
                                 Vfetk_MeshLoad meshType, /* Type of mesh to load */
                                 Vio *sock /* Socket for external mesh data (NULL otherwise) */
                                ) {

    Vrc_Codes vrc; /* Function return codes - see vhal.h for enum */
    int skey = 0;  /* Simplex format */

    /* Load mesh from socket if external mesh, otherwise generate mesh */
    switch (meshType) {
        case VML_EXTERNAL:
            if (sock == VNULL) {
                Vnm_print(2, "Vfetk_loadMesh:  Got NULL socket!\n");
                return VRC_FAILURE;
            }
            AM_read(thee->am, skey, sock);
            Vio_connectFree(sock);
            Vio_bufGive(sock);
            Vio_dtor(&sock);
            break;
        case VML_DIRICUBE:
        case VML_NEUMCUBE:
            /* Create new mesh and store in thee */
            vrc = Vfetk_genCube(thee, center, length, meshType);
            if (vrc == VRC_FAILURE) return VRC_FAILURE;
            break;
        default:
            Vnm_print(2, "Vfetk_loadMesh:  unrecognized mesh type (%d)!\n",
                      meshType);
            return VRC_FAILURE;
    };

    /* Setup charge-simplex map */
    Vnm_print(0, "Vfetk_ctor2:  Constructing Vcsm...\n");
    thee->csm = VNULL;
    /* Construct a new Vcsm with the atom list and gem data */
    thee->csm = Vcsm_ctor(Vpbe_getValist(thee->pbe), thee->gm);
    VASSERT(thee->csm != VNULL);
    Vcsm_init(thee->csm);

    return VRC_SUCCESS;
}


VPUBLIC void Bmat_printHB(Bmat *thee,
                          char *fname
                         ) {

    Mat *Ablock;
    MATsym pqsym;
    int i, j, jj;
    int *IA, *JA;
    double *D, *L, *U;
    FILE *fp;

    char mmtitle[72];
    char mmkey[] = {"8charkey"};
    int totc = 0, ptrc = 0, indc = 0, valc = 0;
    char mxtyp[] = {"RUA"}; /* Real Unsymmetric Assembled */
    int nrow = 0, ncol = 0, numZ = 0;
    int numZdigits = 0, nrowdigits = 0;
    int nptrline = 8, nindline = 8, nvalline = 5;
    char ptrfmt[] = {"(8I10)          "}, ptrfmtstr[] = {"%10d"};
    char indfmt[] = {"(8I10)          "}, indfmtstr[] = {"%10d"};
    char valfmt[] = {"(5E16.8)            "}, valfmtstr[] = {"%16.8E"};

    VASSERT( thee->numB == 1 );             /* HARDWIRE FOR NOW */
    Ablock = thee->AD[0][0];

    VASSERT( Mat_format( Ablock ) == DRC_FORMAT );  /* HARDWIRE FOR NOW */

    pqsym = Mat_sym( Ablock );

    if ( pqsym == IS_SYM ) {
        mxtyp[1] = 'S';
    } else if ( pqsym == ISNOT_SYM ) {
        mxtyp[1] = 'U';
    } else {
        VASSERT( 0 ); /* NOT VALID */
    }

    nrow = Bmat_numRT( thee ); /* Number of rows */
    ncol = Bmat_numCT( thee ); /* Number of cols */
    numZ = Bmat_numZT( thee ); /* Number of entries */

    nrowdigits = (int) (log( nrow )/log( 10 )) + 1;
    numZdigits = (int) (log( numZ )/log( 10 )) + 1;

    nptrline = (int) ( 80 / (numZdigits + 1) );
    nindline = (int) ( 80 / (nrowdigits + 1) );

    sprintf(ptrfmt,"(%dI%d)",nptrline,numZdigits+1);
    sprintf(ptrfmtstr,"%%%dd",numZdigits+1);
    sprintf(indfmt,"(%dI%d)",nindline,nrowdigits+1);
    sprintf(indfmtstr,"%%%dd",nrowdigits+1);

    ptrc = (int) ( ( (ncol + 1) - 1 ) / nptrline ) + 1;
    indc = (int) ( (numZ - 1) / nindline ) + 1;
    valc = (int) ( (numZ - 1) / nvalline ) + 1;

    totc = ptrc + indc + valc;

    sprintf( mmtitle, "Sparse '%s' Matrix - Harwell-Boeing Format - '%s'",
             thee->name, fname );

   /* Step 0:  Open the file for writing */

    fp = fopen( fname, "w" );
    if (fp == VNULL) {
        Vnm_print(2,"Bmat_printHB:  Ouch couldn't open file <%s>\n",fname);
        return;
    }

    /* Step 1:  Print the header information */

    fprintf( fp, "%-72s%-8s\n", mmtitle, mmkey );
    fprintf( fp, "%14d%14d%14d%14d%14d\n", totc, ptrc, indc, valc, 0 );
    fprintf( fp, "%3s%11s%14d%14d%14d\n", mxtyp, " ", nrow, ncol, numZ );
    fprintf( fp, "%-16s%-16s%-20s%-20s\n", ptrfmt, indfmt, valfmt, "6E13.5" );

    IA = Ablock->IA;
    JA = Ablock->JA;
    D = Ablock->diag;
    L = Ablock->offL;
    U = Ablock->offU;

    if ( pqsym == IS_SYM ) {

        /* Step 2:  Print the pointer information */

        for (i=0; i<(ncol+1); i++) {
            fprintf( fp, ptrfmtstr, Ablock->IA[i] + (i+1) );
            if ( ( (i+1) % nptrline ) == 0 ) {
                fprintf( fp, "\n" );
            }
        }

        if ( ( (ncol+1) % nptrline ) != 0 ) {
            fprintf( fp, "\n" );
        }

        /* Step 3:  Print the index information */

        j = 0;
        for (i=0; i<ncol; i++) {
            fprintf( fp, indfmtstr, i+1); /* diagonal */
            if ( ( (j+1) % nindline ) == 0 ) {
                fprintf( fp, "\n" );
            }
            j++;
            for (jj=IA[i]; jj<IA[i+1]; jj++) {
                fprintf( fp, indfmtstr, JA[jj] + 1 ); /* lower triangle */
                if ( ( (j+1) % nindline ) == 0 ) {
                    fprintf( fp, "\n" );
                }
                j++;
            }
        }

        if ( ( j % nindline ) != 0 ) {
            fprintf( fp, "\n" );
        }

        /* Step 4:  Print the value information */

        j = 0;
        for (i=0; i<ncol; i++) {
            fprintf( fp, valfmtstr, D[i] );
            if ( ( (j+1) % nvalline ) == 0 ) {
                fprintf( fp, "\n" );
            }
            j++;
            for (jj=IA[i]; jj<IA[i+1]; jj++) {
                fprintf( fp, valfmtstr, L[jj] );
                if ( ( (j+1) % nvalline ) == 0 ) {
                    fprintf( fp, "\n" );
                }
                j++;
            }
        }

        if ( ( j % nvalline ) != 0 ) {
            fprintf( fp, "\n" );
        }

    } else { /* ISNOT_SYM */

        VASSERT( 0 ); /* NOT CODED YET */
    }

    /* Step 5:  Close the file */
    fclose( fp );
}

VPUBLIC PDE* Vfetk_PDE_ctor(Vfetk *fetk) {

    PDE *thee = VNULL;

    thee = (PDE*)Vmem_malloc(fetk->vmem, 1, sizeof(PDE));
    VASSERT(thee != VNULL);
    VASSERT(Vfetk_PDE_ctor2(thee, fetk));

    return thee;
}

VPUBLIC int Vfetk_PDE_ctor2(PDE *thee,
                            Vfetk *fetk
                           ) {

    int i;

    if (thee == VNULL) {
        Vnm_print(2, "Vfetk_PDE_ctor2:  Got NULL thee!\n");
        return 0;
    }

    /* Store a local copy of the Vfetk class */
    var.fetk = fetk;

    /* PDE-specific parameters and function pointers */
    thee->initAssemble = Vfetk_PDE_initAssemble;
    thee->initElement  = Vfetk_PDE_initElement;
    thee->initFace     = Vfetk_PDE_initFace;
    thee->initPoint    = Vfetk_PDE_initPoint;
    thee->Fu           = Vfetk_PDE_Fu;
    thee->Fu_v         = Vfetk_PDE_Fu_v;
    thee->DFu_wv       = Vfetk_PDE_DFu_wv;
    thee->delta        = Vfetk_PDE_delta;
    thee->u_D          = Vfetk_PDE_u_D;
    thee->u_T          = Vfetk_PDE_u_T;
    thee->Ju           = Vfetk_PDE_Ju;
    thee->vec          = 1; /* FIX! */
    thee->sym[0][0]    = 1;
    thee->est[0]       = 1.0;
    for (i=0; i<VMAX_BDTYPE; i++) thee->bmap[0][i] = i;

    /* Manifold-specific function pointers */
    thee->bisectEdge  = Vfetk_PDE_bisectEdge;
    thee->mapBoundary = Vfetk_PDE_mapBoundary;
    thee->markSimplex = Vfetk_PDE_markSimplex;
    thee->oneChart    = Vfetk_PDE_oneChart;

    /* Element-specific function pointers */
    thee->simplexBasisInit = Vfetk_PDE_simplexBasisInit;
    thee->simplexBasisForm = Vfetk_PDE_simplexBasisForm;

    return 1;
}

VPUBLIC void Vfetk_PDE_dtor(PDE **thee) {

    if ((*thee) != VNULL) {
        Vfetk_PDE_dtor2(*thee);
        /* TODO: The following line is commented out because at the moment,
            there is a seg fault when deallocating at the end of a run. Since
            this routine is called only once at the very end, we'll leave it
            commented out. However, this could be a memory leak.
         */
        /* Vmem_free(var.fetk->vmem, 1, sizeof(PDE), (void **)thee); */
        (*thee) = VNULL;
    }

}

VPUBLIC void Vfetk_PDE_dtor2(PDE *thee) {
    var.fetk = VNULL;
}

VPRIVATE double smooth(int nverts, double dist[VAPBS_NVS], double coeff[VAPBS_NVS], int meth) {

    int i;
    double weight;
    double num = 0.0;
    double den = 0.0;

    for (i=0; i<nverts; i++) {
        if (dist[i] < VSMALL) return coeff[i];
        weight = 1.0/dist[i];
        if (meth == 0) {
            num += (weight * coeff[i]);
            den += weight;
        } else if (meth == 1) {
            /* Small coefficients reset the average to 0; we need to break out
             * of the loop */
            if (coeff[i] < VSMALL) {
                num = 0.0;
                break;
            } else {
                num += weight; den += (weight/coeff[i]);
            }
        } else VASSERT(0);
    }

    return (num/den);

}

VPRIVATE double diel() {

    int i, j;
    double eps, epsp, epsw, dist[5], coeff[5], srad, swin, *vx;
    Vsurf_Meth srfm;
    Vacc *acc;
    PBEparm *pbeparm;

    epsp = Vpbe_getSoluteDiel(var.fetk->pbe);
    epsw = Vpbe_getSolventDiel(var.fetk->pbe);
    VASSERT(var.fetk->pbeparm != VNULL);
    pbeparm = var.fetk->pbeparm;
    srfm = pbeparm->srfm;
    srad = pbeparm->srad;
    swin = pbeparm->swin;
    acc = var.fetk->pbe->acc;

    eps = 0;

    if (VABS(epsp - epsw) < VSMALL) return epsp;
    switch (srfm) {
        case VSM_MOL:
            eps = ((epsw-epsp)*Vacc_molAcc(acc, var.xq, srad) + epsp);
            break;
        case VSM_MOLSMOOTH:
            for (i=0; i<var.nverts; i++) {
                dist[i] = 0;
                vx = var.vx[i];
                for (j=0; j<3; j++) {
                    dist[i] += VSQR(var.xq[j] - vx[j]);
                }
                dist[i] = VSQRT(dist[i]);
                coeff[i] = (epsw-epsp)*Vacc_molAcc(acc, var.xq, srad) + epsp;
            }
            eps = smooth(var.nverts, dist, coeff, 1);
            break;
        case VSM_SPLINE:
            eps = ((epsw-epsp)*Vacc_splineAcc(acc, var.xq, swin, 0.0) + epsp);
            break;
        default:
            Vnm_print(2, "Undefined surface method (%d)!\n", srfm);
            VASSERT(0);
    }

    return eps;
}

VPRIVATE double ionacc() {

    int i, j;
    double dist[5], coeff[5], irad, swin, *vx, accval;
    Vsurf_Meth srfm;
    Vacc *acc = VNULL;
    PBEparm *pbeparm = VNULL;

    VASSERT(var.fetk->pbeparm != VNULL);
    pbeparm = var.fetk->pbeparm;
    srfm = pbeparm->srfm;
    irad = Vpbe_getMaxIonRadius(var.fetk->pbe);
    swin = pbeparm->swin;
    acc = var.fetk->pbe->acc;

    if (var.zks2 < VSMALL) return 0.0;
    switch (srfm) {
        case VSM_MOL:
            accval = Vacc_ivdwAcc(acc, var.xq, irad);
            break;
        case VSM_MOLSMOOTH:
            for (i=0; i<var.nverts; i++) {
                dist[i] = 0;
                vx = var.vx[i];
                for (j=0; j<3; j++) {
                    dist[i] += VSQR(var.xq[j] - vx[j]);
                }
                dist[i] = VSQRT(dist[i]);
                coeff[i] = Vacc_ivdwAcc(acc, var.xq, irad);
            }
            accval = smooth(var.nverts, dist, coeff, 1);
            break;
        case VSM_SPLINE:
            accval = Vacc_splineAcc(acc, var.xq, swin, irad);
            break;
        default:
            Vnm_print(2, "Undefined surface method (%d)!\n", srfm);
            VASSERT(0);
    }

    return accval;
}

VPRIVATE double debye_U(Vpbe *pbe, int d, double x[]) {

    double size, *position, charge, xkappa, eps_w, dist, T, pot, val;
    int iatom, i;
    Valist *alist;
    Vatom *atom;

    eps_w = Vpbe_getSolventDiel(pbe);
    xkappa = (1.0e10)*Vpbe_getXkappa(pbe);
    T = Vpbe_getTemperature(pbe);
    alist = Vpbe_getValist(pbe);
    val = 0;
    pot = 0;

    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        atom = Valist_getAtom(alist, iatom);
        position = Vatom_getPosition(atom);
        charge = Vunit_ec*Vatom_getCharge(atom);
        size = (1e-10)*Vatom_getRadius(atom);
        dist = 0;
        for (i=0; i<d; i++) {
            dist += VSQR(position[i] - x[i]);
        }
        dist = (1.0e-10)*VSQRT(dist);
        val = (charge)/(4*VPI*Vunit_eps0*eps_w*dist);
        if (xkappa != 0.0) {
            val = val*(exp(-xkappa*(dist-size))/(1+xkappa*size));
        }
        pot = pot + val;
    }
    pot = pot*Vunit_ec/(Vunit_kb*T);

    return pot;
}

VPRIVATE double debye_Udiff(Vpbe *pbe, int d, double x[]) {

    double size, *position, charge, eps_p, dist, T, pot, val;
    double Ufull;
    int iatom, i;
    Valist *alist;
    Vatom *atom;

    Ufull = debye_U(pbe, d, x);

    eps_p = Vpbe_getSoluteDiel(pbe);
    T = Vpbe_getTemperature(pbe);
    alist = Vpbe_getValist(pbe);
    val = 0;
    pot = 0;

    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
        atom = Valist_getAtom(alist, iatom);
        position = Vatom_getPosition(atom);
        charge = Vunit_ec*Vatom_getCharge(atom);
        size = (1e-10)*Vatom_getRadius(atom);
        dist = 0;
        for (i=0; i<d; i++) {
            dist += VSQR(position[i] - x[i]);
        }
        dist = (1.0e-10)*VSQRT(dist);
        val = (charge)/(4*VPI*Vunit_eps0*eps_p*dist);
        pot = pot + val;
    }
    pot = pot*Vunit_ec/(Vunit_kb*T);

    pot = Ufull - pot;

    return pot;
}

VPRIVATE void coulomb(Vpbe *pbe, int d, double pt[], double eps, double *U,
  double dU[], double *d2U) {

    int iatom, i;
    double T, pot, fx, fy, fz, x, y, z, scale;
    double *position, charge, dist, dist2, val, vec[3], dUold[3], Uold;
    Valist *alist;
    Vatom *atom;

    /* Initialize variables */
    T = Vpbe_getTemperature(pbe);
    alist = Vpbe_getValist(pbe);
    pot = 0;  fx = 0; fy = 0; fz = 0;
    x = pt[0]; y = pt[1]; z = pt[2];

    /* Calculate */
    if (!Vgreen_coulombD(var.green, 1, &x, &y, &z, &pot, &fx, &fy, &fz)) {
        Vnm_print(2, "Error calculating Green's function!\n");
        VASSERT(0);
    }


    /* Scale the results */
    scale = Vunit_ec/(eps*Vunit_kb*T);
    *U = pot*scale;
    *d2U = 0.0;
    dU[0] = -fx*scale;
    dU[1] = -fy*scale;
    dU[2] = -fz*scale;

#if 0
    /* Compare with old results */
    val = 0.0;
    Uold = 0.0; dUold[0] = 0.0; dUold[1] = 0.0; dUold[2] = 0.0;
    for (iatom=0; iatom<Valist_getNumberAtoms(alist); iatom++) {
       atom = Valist_getAtom(alist, iatom);
       position = Vatom_getPosition(atom);
       charge = Vatom_getCharge(atom);
       dist2 = 0;
       for (i=0; i<d; i++) {
           vec[i] = (position[i] - pt[i]);
           dist2 += VSQR(vec[i]);
       }
       dist = VSQRT(dist2);

       /* POTENTIAL */
       Uold = Uold + charge/dist;

       /* GRADIENT */
       for (i=0; i<d; i++) dUold[i] = dUold[i] + vec[i]*charge/(dist2*dist);

    }
    Uold = Uold*VSQR(Vunit_ec)*(1.0e10)/(4*VPI*Vunit_eps0*eps*Vunit_kb*T);
    for (i=0; i<d; i++) {
        dUold[i] = dUold[i]*VSQR(Vunit_ec)*(1.0e10)/(4*VPI*Vunit_eps0*eps*Vunit_kb*T);
    }

    printf("Unew - Uold = %g - %g = %g\n", *U, Uold, (*U - Uold));
    printf("||dUnew - dUold||^2 = %g\n", (VSQR(dU[0] - dUold[0])
                + VSQR(dU[1] - dUold[1]) + VSQR(dU[2] - dUold[2])));
    printf("dUnew[0] = %g, dUold[0] = %g\n", dU[0], dUold[0]);
    printf("dUnew[1] = %g, dUold[1] = %g\n", dU[1], dUold[1]);
    printf("dUnew[2] = %g, dUold[2] = %g\n", dU[2], dUold[2]);

#endif

}

VPUBLIC void Vfetk_PDE_initAssemble(PDE *thee, int ip[], double rp[]) {

#if 1
    /* Re-initialize the Green's function oracle in case the atom list has
     * changed */
    if (var.initGreen) {
        Vgreen_dtor(&(var.green));
        var.initGreen = 0;
    }
    var.green = Vgreen_ctor(var.fetk->pbe->alist);
    var.initGreen = 1;
#else
    if (!var.initGreen) {
        var.green = Vgreen_ctor(var.fetk->pbe->alist);
        var.initGreen = 1;
    }
#endif

}

VPUBLIC void Vfetk_PDE_initElement(PDE *thee, int elementType, int chart,
  double tvx[][3], void *data) {

    int i, j;
    double epsp, epsw;

    /* We assume that the simplex has been passed in as the void *data * *
     * argument.  Store it */
    VASSERT(data != NULL);
    var.simp = (SS *)data;

    /* save the element type */
    var.sType = elementType;

    /* Grab the vertices from this simplex */
    var.nverts = thee->dim+1;
    for (i=0; i<thee->dim+1; i++) var.verts[i] = SS_vertex(var.simp, i);

    /* Vertex locations of this simplex */
    for (i=0; i<thee->dim+1; i++) {
        for (j=0; j<thee->dim; j++) {
            var.vx[i][j] = tvx[i][j];
        }
    }

    /* Set the dielectric constant for this element for use in the jump term *
     * of the residual-based error estimator.  The value is set to the average
     * * value of the vertices */
    var.jumpDiel = 0;  /* NOT IMPLEMENTED YET! */
}

VPUBLIC void Vfetk_PDE_initFace(PDE *thee, int faceType, int chart,
  double tnvec[]) {

    int i;

    /* unit normal vector of this face */
    for (i=0; i<thee->dim; i++) var.nvec[i] = tnvec[i];

    /* save the face type */
    var.fType = faceType;
}

VPUBLIC void Vfetk_PDE_initPoint(PDE *thee, int pointType, int chart,
  double txq[], double tU[], double tdU[][3]) {

    int i, j, ichop;
    double u2, coef2, eps_p;
    Vhal_PBEType pdetype;
    Vpbe *pbe = VNULL;

    eps_p = Vpbe_getSoluteDiel(var.fetk->pbe);
    pdetype = var.fetk->type;
    pbe = var.fetk->pbe;

    /* the point, the solution value and gradient, and the Coulomb value and *
     * gradient at the point */
    if ((pdetype == PBE_LRPBE) || (pdetype == PBE_NRPBE)) {
        coulomb(pbe, thee->dim, txq, eps_p, &(var.W), var.dW, &(var.d2W));
    }
    for (i=0; i<thee->vec; i++) {
        var.U[i] = tU[i];
        for (j=0; j<thee->dim; j++) {
            var.xq[j] = txq[j];
            var.dU[i][j] = tdU[i][j];
        }
    }

    /* interior form case */
    if (pointType == 0) {

        /* Get the dielectric values */
        var.diel  = diel();
        var.ionacc  = ionacc();
        var.A = var.diel;
        var.F = (var.diel - eps_p);

        switch (pdetype) {

            case PBE_LPBE:
                var.DB = var.ionacc*var.zkappa2*var.ionstr;
                var.B = var.DB*var.U[0];
                break;

            case PBE_NPBE:

                var.B  = 0;
                var.DB  = 0;
                if ((var.ionacc > VSMALL) && (var.zks2 > VSMALL)) {
                    for (i=0; i<var.nion; i++) {
                        u2 = -1.0 * var.U[0] * var.ionQ[i];

                        /* NONLINEAR TERM */
                        coef2 = -1.0 * var.ionacc * var.zks2 * var.ionConc[i];
                        var.B += (coef2 * Vcap_exp(u2, &ichop));
                        /* LINEARIZED TERM */
                        coef2 = -1.0 * var.ionQ[i] * coef2;
                        var.DB += (coef2 * Vcap_exp(u2, &ichop));
                    }
                }
                break;

            case PBE_LRPBE:
                var.DB = var.ionacc*var.zkappa2*var.ionstr;
                var.B = var.DB*(var.U[0]+var.W);
                break;

            case PBE_NRPBE:

                var.B  = 0;
                var.DB  = 0;
                if ((var.ionacc > VSMALL) && (var.zks2 > VSMALL)) {
                    for (i=0; i<var.nion; i++) {
                        u2 = -1.0 * (var.U[0] + var.W) * var.ionQ[i];

                        /* NONLINEAR TERM */
                        coef2 = -1.0 * var.ionacc * var.zks2 * var.ionConc[i];
                        var.B += (coef2 * Vcap_exp(u2, &ichop));

                        /* LINEARIZED TERM */
                        coef2 = -1.0 * var.ionQ[i] * coef2;
                        var.DB += (coef2 * Vcap_exp(u2, &ichop));
                    }
                }
                break;

            case PBE_SMPBE: /* SMPBE Temp */

                var.B  = 0;
                var.DB  = 0;
                if ((var.ionacc > VSMALL) && (var.zks2 > VSMALL)) {
                    for (i=0; i<var.nion; i++) {
                        u2 = -1.0 * var.U[0] * var.ionQ[i];

                        /* NONLINEAR TERM */
                        coef2 = -1.0 * var.ionacc * var.zks2 * var.ionConc[i];
                        var.B += (coef2 * Vcap_exp(u2, &ichop));
                        /* LINEARIZED TERM */
                        coef2 = -1.0 * var.ionQ[i] * coef2;
                        var.DB += (coef2 * Vcap_exp(u2, &ichop));
                    }
                }
                    break;
            default:
                Vnm_print(2, "Vfetk_PDE_initPoint:  Unknown PBE type (%d)!\n",
                  pdetype);
                VASSERT(0);
                break;
        }


    /* boundary form case */
    } else {
#ifdef DONEUMANN
        ;
#else
        Vnm_print(2, "Vfetk:  Whoa!  I just got a boundary point to evaluate (%d)!\n", pointType);
        Vnm_print(2, "Vfetk:  Did you do that on purpose?\n");
#endif
    }

#if 0 /* THIS IS VERY NOISY! */
    Vfetk_dumpLocalVar();
#endif

}

VPUBLIC void Vfetk_PDE_Fu(PDE *thee, int key, double F[]) {

    //Vnm_print(2, "Vfetk_PDE_Fu:  Setting error to zero!\n");

    F[0] = 0.;

}

VPUBLIC double Vfetk_PDE_Fu_v(
        PDE *thee,
        int key,
        double V[],
        double dV[][VAPBS_DIM]
        ) {

    Vhal_PBEType type;
    int i;
    double value = 0.;

    type = var.fetk->type;

    /* interior form case */
    if (key == 0) {

        for (i=0; i<thee->dim; i++) value += ( var.A * var.dU[0][i] * dV[0][i] );
        value += var.B * V[0];

        if ((type == PBE_LRPBE) || (type == PBE_NRPBE)) {
            for (i=0; i<thee->dim; i++) {
                if (var.F > VSMALL) value += (var.F * var.dW[i] * dV[0][i]);
            }
        }

    /* boundary form case */
    } else {
#ifdef DONEUMANN
        value = 0.0;
#else
        Vnm_print(2, "Vfetk:  Whoa! I was just asked to evaluate a boundary weak form for point type %d!\n", key);
#endif
    }

    var.Fu_v = value;
    return value;
}

VPUBLIC double Vfetk_PDE_DFu_wv(
        PDE *thee,
        int key,
        double W[],
        double dW[][VAPBS_DIM],
        double V[],
        double dV[][3]
        ) {

    Vhal_PBEType type;
    int i;
    double value = 0.;

    type = var.fetk->type;

    /* Interior form */
    if (key == 0) {
            value += var.DB * W[0] * V[0];
            for (i=0; i<thee->dim; i++) value += ( var.A * dW[0][i] * dV[0][i] );

    /* boundary form case */
    } else {
#ifdef DONEUMANN
        value = 0.0;
#else
        Vnm_print(2, "Vfetk:  Whoa! I was just asked to evaluate a boundary weak form for point type %d!\n", key);
#endif
    }

    var.DFu_wv = value;
    return value;
}

/** @brief  Maximum number of simplices in a simplex ring
 *  @ingroup  Vfetk */
#define VRINGMAX 1000
/** @brief  Maximum number of atoms associated with a vertex
 *  @ingroup  Vfetk */
#define VATOMMAX 1000000
VPUBLIC void Vfetk_PDE_delta(PDE *thee, int type, int chart, double txq[],
  void *user, double F[]) {

    int iatom, jatom, natoms, atomIndex, atomList[VATOMMAX], nAtomList;
    int gotAtom, numSring, isimp, ivert, sid;
    double *position, charge, phi[VAPBS_NVS], phix[VAPBS_NVS][3], value;
    Vatom *atom;
    Vhal_PBEType pdetype;
    SS *sring[VRINGMAX];
    VV *vertex = (VV *)user;

    pdetype = var.fetk->type;

    F[0] = 0.0;

    if ((pdetype == PBE_LPBE) || (pdetype == PBE_NPBE) || (pdetype == PBE_SMPBE) /* SMPBE Added */) {
        VASSERT( vertex != VNULL);
        numSring = 0;
        sring[numSring] = VV_firstSS(vertex);
        while (sring[numSring] != VNULL) {
            numSring++;
            sring[numSring] = SS_link(sring[numSring-1], vertex);
        }
        VASSERT( numSring > 0 );
        VASSERT( numSring <= VRINGMAX );

        /* Move around the simplex ring and determine the charge locations */
        F[0] = 0.;
        charge = 0.;
        nAtomList = 0;
        for (isimp=0; isimp<numSring; isimp++) {
            sid = SS_id(sring[isimp]);
            natoms = Vcsm_getNumberAtoms(Vfetk_getVcsm(var.fetk), sid);
            for (iatom=0; iatom<natoms; iatom++) {
                /* Get the delta function information * */
                atomIndex = Vcsm_getAtomIndex(Vfetk_getVcsm(var.fetk),
                  iatom, sid);
                gotAtom = 0;
                for (jatom=0; jatom<nAtomList; jatom++) {
                    if (atomList[jatom] == atomIndex) {
                        gotAtom = 1;
                        break;
                    }
                }
                if (!gotAtom) {
                    VASSERT(nAtomList < VATOMMAX);
                    atomList[nAtomList] = atomIndex;
                    nAtomList++;

                    atom = Vcsm_getAtom(Vfetk_getVcsm(var.fetk), iatom, sid);
                    charge = Vatom_getCharge(atom);
                    position = Vatom_getPosition(atom);

                    /* Get the test function value at the delta function I
                     * used to do a VASSERT to make sure the point was in the
                     * simplex (i.e., make sure round-off error isn't an
                     * issue), but round off errors became an issue */
                    if (!Gem_pointInSimplexVal(Vfetk_getGem(var.fetk),
                      sring[isimp], position, phi, phix)) {
                        if (!Gem_pointInSimplex(Vfetk_getGem(var.fetk),
                          sring[isimp], position)) {
                            Vnm_print(2, "delta: Both Gem_pointInSimplexVal \
and Gem_pointInSimplex detected misplaced point charge!\n");
                            Vnm_print(2, "delta: I think you have problems: \
phi = {");
                            for (ivert=0; ivert<Gem_dimVV(Vfetk_getGem(var.fetk)); ivert++) Vnm_print(2, "%e ", phi[ivert]);
                                                                                                            Vnm_print(2, "}\n");
                        }
                    }
                    value = 0;
                    for (ivert=0; ivert<Gem_dimVV(Vfetk_getGem(var.fetk)); ivert++) {
                        if (VV_id(SS_vertex(sring[isimp], ivert)) == VV_id(vertex)) value += phi[ivert];
                    }

                    F[0] += (value * Vpbe_getZmagic(var.fetk->pbe) * charge);
                } /* if !gotAtom */
            } /* for iatom */
        } /* for isimp */

    } else if ((pdetype == PBE_LRPBE) || (pdetype == PBE_NRPBE)) {
        F[0] = 0.0;
    } else { VASSERT(0); }

    var.delta = F[0];

}

VPUBLIC void Vfetk_PDE_u_D(PDE *thee, int type, int chart, double txq[],
  double F[]) {

    if ((var.fetk->type == PBE_LPBE) || (var.fetk->type == PBE_NPBE) || (var.fetk->type == PBE_SMPBE) /* SMPBE Added */) {
        F[0] = debye_U(var.fetk->pbe, thee->dim, txq);
    } else if ((var.fetk->type == PBE_LRPBE) || (var.fetk->type == PBE_NRPBE)) {
        F[0] = debye_Udiff(var.fetk->pbe, thee->dim, txq);
    } else VASSERT(0);

    var.u_D = F[0];

}

/**
 * The signature here doesn't match what's in mc's src/pde/mc/pde.h, which
 * g++ seems to dislike for GAMer integration.  Trying a change of function
 * signature to match to see if that makes g++ happy.  Also see vfetk.h for
 * similar signature change. - P. Ellis 11-8-2011
 */
VPUBLIC void Vfetk_PDE_u_T(PDE *thee, int type, int chart, double txq[],
  double F[]) {
/*VPUBLIC void Vfetk_PDE_u_T(sPDE *thee,
                           int type,
                           int chart,
                           double txq[],
                           double F[],
                           double dF[][3]
                          ) { */

    F[0] = 0.0;
    var.u_T = F[0];

}


VPUBLIC void Vfetk_PDE_bisectEdge(int dim, int dimII, int edgeType,
  int chart[], double vx[][3]) {

    int i;

    for (i=0; i<dimII; i++) vx[2][i] = .5 * (vx[0][i] + vx[1][i]);
    chart[2] = chart[0];

}

VPUBLIC void Vfetk_PDE_mapBoundary(int dim, int dimII, int vertexType,
  int chart, double vx[3]) {

}

VPUBLIC int Vfetk_PDE_markSimplex(int dim, int dimII, int simplexType,
  int faceType[VAPBS_NVS], int vertexType[VAPBS_NVS], int chart[], double vx[][3],
  void *simplex) {

    double targetRes, edgeLength, srad, swin, myAcc, refAcc;
    int i, natoms;
    Vsurf_Meth srfm;
    Vhal_PBEType type;
    FEMparm *feparm = VNULL;
    PBEparm *pbeparm = VNULL;
    Vpbe *pbe = VNULL;
    Vacc *acc = VNULL;
    Vcsm *csm = VNULL;
    SS *simp = VNULL;

    VASSERT(var.fetk->feparm != VNULL);
    feparm = var.fetk->feparm;
    VASSERT(var.fetk->pbeparm != VNULL);;
    pbeparm = var.fetk->pbeparm;
    pbe = var.fetk->pbe;
    csm = Vfetk_getVcsm(var.fetk);
    acc = pbe->acc;
    targetRes = feparm->targetRes;
    srfm = pbeparm->srfm;
    srad = pbeparm->srad;
    swin = pbeparm->swin;
    simp = (SS *)simplex;
    type = var.fetk->type;

    /* Check to see if this simplex is smaller than the target size */
    /* NAB WARNING:  I am providing face=-1 here to conform to the new MC API; however, I'm not sure if this is the correct behavior. */
    Gem_longestEdge(var.fetk->gm, simp, -1, &edgeLength);
    if (edgeLength < targetRes) return 0;

    /* For non-regularized PBE, check charge-simplex map */
    if ((type == PBE_LPBE) || (type == PBE_NPBE) || (type == PBE_SMPBE)  /* SMPBE Added */) {
        natoms = Vcsm_getNumberAtoms(csm, SS_id(simp));
        if (natoms > 0) {
            return 1;
        }
    }

    /* We would like to resolve the mesh between the van der Waals surface the
     * max distance from this surface where there could be coefficient
     * changes */
    switch(srfm) {
        case VSM_MOL:
            refAcc = Vacc_molAcc(acc, vx[0], srad);
            for (i=1; i<(dim+1); i++) {
                myAcc = Vacc_molAcc(acc, vx[i], srad);
                if (myAcc != refAcc) {
                    return 1;
                }
            }
            break;
        case VSM_MOLSMOOTH:
            refAcc = Vacc_molAcc(acc, vx[0], srad);
            for (i=1; i<(dim+1); i++) {
                myAcc = Vacc_molAcc(acc, vx[i], srad);
                if (myAcc != refAcc) {
                    return 1;
                }
            }
            break;
        case VSM_SPLINE:
            refAcc = Vacc_splineAcc(acc, vx[0], swin, 0.0);
            for (i=1; i<(dim+1); i++) {
                myAcc = Vacc_splineAcc(acc, vx[i], swin, 0.0);
                if (myAcc != refAcc) {
                    return 1;
                }
            }
            break;
        default:
            VASSERT(0);
            break;
    }

    return 0;
}

VPUBLIC void Vfetk_PDE_oneChart(int dim, int dimII, int objType, int chart[],
  double vx[][3], int dimV) {

}

VPUBLIC double Vfetk_PDE_Ju(PDE *thee, int key) {

    int i, ichop;
    double dielE, qmE, coef2, u2;
    double value = 0.;
    Vhal_PBEType type;

    type = var.fetk->type;

    /* interior form case */
    if (key == 0) {
        dielE = 0;
        for (i=0; i<3; i++) dielE += VSQR(var.dU[0][i]);
        dielE = dielE*var.diel;

        switch (type) {
            case PBE_LPBE:
                if ((var.ionacc > VSMALL) && (var.zkappa2 > VSMALL)) {
                    qmE = var.ionacc*var.zkappa2*VSQR(var.U[0]);
                } else qmE = 0;
                break;
            case PBE_NPBE:
                if ((var.ionacc > VSMALL) && (var.zks2 > VSMALL)) {
                    qmE = 0.;
                    for (i=0; i<var.nion; i++) {
                        coef2 = var.ionacc * var.zks2 * var.ionConc[i] * var.ionQ[i];
                        u2 = -1.0 * (var.U[0]) * var.ionQ[i];
                        qmE += (coef2 * (Vcap_exp(u2, &ichop) - 1.0));
                    }
                } else qmE = 0;
                break;
            case PBE_LRPBE:
                if ((var.ionacc > VSMALL) && (var.zkappa2 > VSMALL)) {
                    qmE = var.ionacc*var.zkappa2*VSQR((var.U[0] + var.W));
                } else qmE = 0;
                break;
            case PBE_NRPBE:
                if ((var.ionacc > VSMALL) && (var.zks2 > VSMALL)) {
                    qmE = 0.;
                    for (i=0; i<var.nion; i++) {
                        coef2 = var.ionacc * var.zks2 * var.ionConc[i] * var.ionQ[i];
                        u2 = -1.0 * (var.U[0] + var.W) * var.ionQ[i];
                        qmE += (coef2 * (Vcap_exp(u2, &ichop) - 1.0));
                    }
                } else qmE = 0;
                break;
            case PBE_SMPBE: /* SMPBE Temp */
                if ((var.ionacc > VSMALL) && (var.zks2 > VSMALL)) {
                    qmE = 0.;
                    for (i=0; i<var.nion; i++) {
                        coef2 = var.ionacc * var.zks2 * var.ionConc[i] * var.ionQ[i];
                        u2 = -1.0 * (var.U[0]) * var.ionQ[i];
                        qmE += (coef2 * (Vcap_exp(u2, &ichop) - 1.0));
                    }
                } else qmE = 0;
                break;
            default:
                Vnm_print(2, "Vfetk_PDE_Ju:  Invalid PBE type (%d)!\n", type);
                VASSERT(0);
                break;
        }

        value = 0.5*(dielE + qmE)/Vpbe_getZmagic(var.fetk->pbe);

    /* boundary form case */
    } else if (key == 1) {
        value = 0.0;

    /* how did we get here? */
    } else VASSERT(0);

    return value;

}

VPUBLIC void Vfetk_externalUpdateFunction(SS **simps, int num) {

    Vcsm *csm = VNULL;
    int rc;

    VASSERT(var.fetk != VNULL);
    csm = Vfetk_getVcsm(var.fetk);
    VASSERT(csm != VNULL);

    rc = Vcsm_update(csm, simps, num);

    if (!rc) {
        Vnm_print(2, "Error while updating charge-simplex map!\n");
        VASSERT(0);
    }
}

VPRIVATE void polyEval(int numP, double p[], double c[][VMAXP], double xv[]) {
    int i;
    double x, y, z;

    x = xv[0];
    y = xv[1];
    z = xv[2];
    for (i=0; i<numP; i++) {
        p[i] = c[i][0]
             + c[i][1]  * x
             + c[i][2]  * y
             + c[i][3]  * z
             + c[i][4]  * x*x
             + c[i][5]  * y*y
             + c[i][6]  * z*z
             + c[i][7]  * x*y
             + c[i][8]  * x*z
             + c[i][9]  * y*z
             + c[i][10] * x*x*x
             + c[i][11] * y*y*y
             + c[i][12] * z*z*z
             + c[i][13] * x*x*y
             + c[i][14] * x*x*z
             + c[i][15] * x*y*y
             + c[i][16] * y*y*z
             + c[i][17] * x*z*z
             + c[i][18] * y*z*z;
    }
}

VPRIVATE void setCoef(int numP, double c[][VMAXP], double cx[][VMAXP],
  double cy[][VMAXP], double cz[][VMAXP], int ic[][VMAXP], int icx[][VMAXP],
  int icy[][VMAXP], int icz[][VMAXP]) {

    int i, j;
    for (i=0; i<numP; i++) {
        for (j=0; j<VMAXP; j++) {
            c[i][j]  = 0.5 * (double)ic[i][j];
            cx[i][j] = 0.5 * (double)icx[i][j];
            cy[i][j] = 0.5 * (double)icy[i][j];
            cz[i][j] = 0.5 * (double)icz[i][j];
        }
    }
}

VPUBLIC int Vfetk_PDE_simplexBasisInit(int key, int dim, int comp, int *ndof,
  int dof[]) {

    int qorder, bump, dimIS[VAPBS_NVS];

    /* necessary quadrature order to return at the end */
    qorder = P_DEG;

    /* deal with bump function requests */
    if ((key == 0) || (key == 1)) {
        bump = 0;
    } else if ((key == 2) || (key == 3)) {
        bump = 1;
    } else { VASSERT(0); }

    /* for now use same element for all components, both trial and test */
    if (dim==2) {
        /* 2D simplex dimensions */
        dimIS[0] = 3;  /* number of vertices             */
        dimIS[1] = 3;  /* number of edges                */
        dimIS[2] = 0;  /* number of faces (3D only)      */
        dimIS[3] = 1;  /* number of simplices (always=1) */
        if (bump==0) {
            if (P_DEG==1) {
                init_2DP1(dimIS, ndof, dof, c, cx, cy, cz);
            } else if (P_DEG==2) {
                init_2DP1(dimIS, ndof, dof, c, cx, cy, cz);
            } else if (P_DEG==3) {
                init_2DP1(dimIS, ndof, dof, c, cx, cy, cz);
            } else Vnm_print(2, "..bad order..");
        } else if (bump==1) {
            if (P_DEG==1) {
                init_2DP1(dimIS, ndof, dof, c, cx, cy, cz);
            } else Vnm_print(2, "..bad order..");
        } else Vnm_print(2, "..bad bump..");
    } else if (dim==3) {
        /* 3D simplex dimensions */
        dimIS[0] = 4;  /* number of vertices             */
        dimIS[1] = 6;  /* number of edges                */
        dimIS[2] = 4;  /* number of faces (3D only)      */
        dimIS[3] = 1;  /* number of simplices (always=1) */
        if (bump==0) {
            if (P_DEG==1) {
                init_3DP1(dimIS, ndof, dof, c, cx, cy, cz);
            } else if (P_DEG==2) {
                init_3DP1(dimIS, ndof, dof, c, cx, cy, cz);
            } else if (P_DEG==3) {
                init_3DP1(dimIS, ndof, dof, c, cx, cy, cz);
            } else Vnm_print(2, "..bad order..");
        } else if (bump==1) {
            if (P_DEG==1) {
                init_3DP1(dimIS, ndof, dof, c, cx, cy, cz);
            } else Vnm_print(2, "..bad order..");
        } else Vnm_print(2, "..bad bump..");
    } else Vnm_print(2, "..bad dimension..");

    /* save number of DF */
    numP = *ndof;

    /* return the required quarature order */
    return qorder;
}

VPUBLIC void Vfetk_PDE_simplexBasisForm(int key, int dim, int comp, int pdkey,
  double xq[], double basis[]) {

    if (pdkey == 0) {
        polyEval(numP, basis, c, xq);
    } else if (pdkey == 1) {
        polyEval(numP, basis, cx, xq);
    } else if (pdkey == 2) {
        polyEval(numP, basis, cy, xq);
    } else if (pdkey == 3) {
        polyEval(numP, basis, cz, xq);
    } else { VASSERT(0); }
}

VPRIVATE void init_2DP1(int dimIS[], int *ndof, int dof[], double c[][VMAXP],
  double cx[][VMAXP], double cy[][VMAXP], double cz[][VMAXP]) {

    int i;

    /* dof number and locations */
    dof[0] = 1;
    dof[1] = 0;
    dof[2] = 0;
    dof[3] = 0;
    *ndof  = 0;
    for (i=0; i<VAPBS_NVS; i++) *ndof += dimIS[i] * dof[i];
    VASSERT( *ndof == dim_2DP1 );
    VASSERT( *ndof <= VMAXP );

    /* coefficients of the polynomials */
    setCoef( *ndof, c, cx, cy, cz, lgr_2DP1, lgr_2DP1x, lgr_2DP1y, lgr_2DP1z );
}

VPRIVATE void init_3DP1(int dimIS[], int *ndof, int dof[], double c[][VMAXP],
  double cx[][VMAXP], double cy[][VMAXP], double cz[][VMAXP]) {

    int i;

    /* dof number and locations */
    dof[0] = 1;
    dof[1] = 0;
    dof[2] = 0;
    dof[3] = 0;
    *ndof  = 0;
    for (i=0; i<VAPBS_NVS; i++) *ndof += dimIS[i] * dof[i];
    VASSERT( *ndof == dim_3DP1 );
    VASSERT( *ndof <= VMAXP );

    /* coefficients of the polynomials */
    setCoef( *ndof, c, cx, cy, cz, lgr_3DP1, lgr_3DP1x, lgr_3DP1y, lgr_3DP1z );
}

VPUBLIC void Vfetk_dumpLocalVar() {

    int i;

    Vnm_print(1, "DEBUG: nvec = (%g, %g, %g)\n", var.nvec[0], var.nvec[1],
      var.nvec[2]);
    Vnm_print(1, "DEBUG: nverts = %d\n", var.nverts);
    for (i=0; i<var.nverts; i++) {
        Vnm_print(1, "DEBUG: verts[%d] ID = %d\n", i, VV_id(var.verts[i]));
        Vnm_print(1, "DEBUG: vx[%d] = (%g, %g, %g)\n", i, var.vx[i][0],
          var.vx[i][1], var.vx[i][2]);
    }
    Vnm_print(1, "DEBUG: simp ID = %d\n", SS_id(var.simp));
    Vnm_print(1, "DEBUG: sType = %d\n", var.sType);
    Vnm_print(1, "DEBUG: fType = %d\n", var.fType);
    Vnm_print(1, "DEBUG: xq = (%g, %g, %g)\n", var.xq[0], var.xq[1], var.xq[2]);
    Vnm_print(1, "DEBUG: U[0] = %g\n", var.U[0]);
    Vnm_print(1, "DEBUG: dU[0] = (%g, %g, %g)\n", var.dU[0][0], var.dU[0][1],
      var.dU[0][2]);
    Vnm_print(1, "DEBUG: W = %g\n", var.W);
    Vnm_print(1, "DEBUG: d2W = %g\n", var.d2W);
    Vnm_print(1, "DEBUG: dW = (%g, %g, %g)\n", var.dW[0], var.dW[1], var.dW[2]);
    Vnm_print(1, "DEBUG: diel = %g\n", var.diel);
    Vnm_print(1, "DEBUG: ionacc = %g\n", var.ionacc);
    Vnm_print(1, "DEBUG: A = %g\n", var.A);
    Vnm_print(1, "DEBUG: F = %g\n", var.F);
    Vnm_print(1, "DEBUG: B = %g\n", var.B);
    Vnm_print(1, "DEBUG: DB = %g\n", var.DB);
    Vnm_print(1, "DEBUG: nion = %d\n", var.nion);
    for (i=0; i<var.nion; i++) {
        Vnm_print(1, "DEBUG: ionConc[%d] = %g\n", i, var.ionConc[i]);
        Vnm_print(1, "DEBUG: ionQ[%d] = %g\n", i, var.ionQ[i]);
        Vnm_print(1, "DEBUG: ionRadii[%d] = %g\n", i, var.ionRadii[i]);
    }
    Vnm_print(1, "DEBUG: zkappa2 = %g\n", var.zkappa2);
    Vnm_print(1, "DEBUG: zks2 = %g\n", var.zks2);
    Vnm_print(1, "DEBUG: Fu_v = %g\n", var.Fu_v);
    Vnm_print(1, "DEBUG: DFu_wv = %g\n", var.DFu_wv);
    Vnm_print(1, "DEBUG: delta = %g\n", var.delta);
    Vnm_print(1, "DEBUG: u_D = %g\n", var.u_D);
    Vnm_print(1, "DEBUG: u_T = %g\n", var.u_T);

};

VPUBLIC int Vfetk_fillArray(Vfetk *thee, Bvec *vec, Vdata_Type type) {

    int i, j, ichop;
    double coord[3], chi, q, conc, val;
    VV *vert;
    Bvec *u, *u_d;
    AM *am;
    Gem *gm;
    PBEparm *pbeparm;
    Vacc *acc;
    Vpbe *pbe;

    gm = thee->gm;
    am = thee->am;
    pbe = thee->pbe;
    pbeparm = thee->pbeparm;
    acc = pbe->acc;

    /* Make sure vec has enough rows to accomodate the vertex data */
    if (Bvec_numRB(vec, 0) != Gem_numVV(gm)) {
        Vnm_print(2, "Vfetk_fillArray:  insufficient space in Bvec!\n");
        Vnm_print(2, "Vfetk_fillArray:  Have %d, need %d!\n", Bvec_numRB(vec, 0),
          Gem_numVV(gm));
        return 0;
    }

    switch (type) {

        case VDT_CHARGE:
            Vnm_print(2, "Vfetk_fillArray:  can't write out charge distribution!\n");
            return 0;
            break;

        case VDT_POT:
            u = am->u;
            u_d = am->ud;
            /* Copy in solution */
            Bvec_copy(vec, u);
            /* Add dirichlet condition */
            Bvec_axpy(vec, u_d, 1.0);
            break;

        case VDT_SMOL:
            for (i=0; i<Gem_numVV(gm); i++) {
                vert = Gem_VV(gm, i);
                for (j=0; j<3; j++) coord[j] = VV_coord(vert, j);
                chi = Vacc_molAcc(acc, coord, pbe->solventRadius);
                Bvec_set(vec, i, chi);
            }
            break;

        case VDT_SSPL:
            for (i=0; i<Gem_numVV(gm); i++) {
                vert = Gem_VV(gm, i);
                for (j=0; j<3; j++) coord[j] = VV_coord(vert, j);
                chi = Vacc_splineAcc(acc, coord, pbeparm->swin, 0.0);
                Bvec_set(vec, i, chi);
            }
            break;

        case VDT_VDW:
            for (i=0; i<Gem_numVV(gm); i++) {
                vert = Gem_VV(gm, i);
                for (j=0; j<3; j++) coord[j] = VV_coord(vert, j);
                chi = Vacc_vdwAcc(acc, coord);
                Bvec_set(vec, i, chi);
            }
            break;

        case VDT_IVDW:
            for (i=0; i<Gem_numVV(gm); i++) {
                vert = Gem_VV(gm, i);
                for (j=0; j<3; j++) coord[j] = VV_coord(vert, j);
                chi = Vacc_ivdwAcc(acc, coord, pbe->maxIonRadius);
                Bvec_set(vec, i, chi);
            }
            break;

        case VDT_LAP:
            Vnm_print(2, "Vfetk_fillArray:  can't write out Laplacian!\n");
            return 0;
            break;

        case VDT_EDENS:
            Vnm_print(2, "Vfetk_fillArray:  can't write out energy density!\n");
            return 0;
            break;

        case VDT_NDENS:
            u = am->u;
            u_d = am->ud;
            /* Copy in solution */
            Bvec_copy(vec, u);
            /* Add dirichlet condition */
            Bvec_axpy(vec, u_d, 1.0);
            /* Load up ions */
            ichop = 0;
            for (i=0; i<Gem_numVV(gm); i++) {
                val = 0;
                for (j=0; j<pbe->numIon; j++) {
                    q = pbe->ionQ[j];
                    conc = pbe->ionConc[j];
                    if (thee->type == PBE_NPBE || thee->type == PBE_SMPBE /* SMPBE Added */) {
                        val += (conc*Vcap_exp(-q*Bvec_val(vec, i), &ichop));
                    } else if (thee->type == PBE_LPBE) {
                        val += (conc * ( 1 - q*Bvec_val(vec, i)));
                    }
                }
                Bvec_set(vec, i, val);
            }
            break;

        case VDT_QDENS:
            u = am->u;
            u_d = am->ud;
            /* Copy in solution */
            Bvec_copy(vec, u);
            /* Add dirichlet condition */
            Bvec_axpy(vec, u_d, 1.0);
            /* Load up ions */
            ichop = 0;
            for (i=0; i<Gem_numVV(gm); i++) {
                val = 0;
                for (j=0; j<pbe->numIon; j++) {
                    q = pbe->ionQ[j];
                    conc = pbe->ionConc[j];
                    if (thee->type == PBE_NPBE || thee->type == PBE_SMPBE /* SMPBE Added */) {
                        val += (q*conc*Vcap_exp(-q*Bvec_val(vec, i), &ichop));
                    } else if (thee->type == PBE_LPBE) {
                        val += (q*conc*(1 - q*Bvec_val(vec, i)));
                    }
                }
                Bvec_set(vec, i, val);
            }
            break;

        case VDT_DIELX:
            Vnm_print(2, "Vfetk_fillArray:  can't write out x-shifted diel!\n");
            return 0;
            break;

        case VDT_DIELY:
            Vnm_print(2, "Vfetk_fillArray:  can't write out y-shifted diel!\n");
            return 0;
            break;

        case VDT_DIELZ:
            Vnm_print(2, "Vfetk_fillArray:  can't write out z-shifted diel!\n");
            return 0;
            break;

        case VDT_KAPPA:
            Vnm_print(2, "Vfetk_fillArray:  can't write out kappa!\n");
            return 0;
            break;

        default:
            Vnm_print(2, "Vfetk_fillArray:  invalid data type (%d)!\n", type);
            return 0;
            break;
    }

    return 1;
}

VPUBLIC int Vfetk_write(Vfetk *thee,  const char *iodev, const char *iofmt,
  const char *thost, const char *fname, Bvec *vec, Vdata_Format format) {

    int i, j, ichop;
    Aprx *aprx;
    Gem *gm;
    Vio *sock;

    VASSERT(thee != VNULL);
    aprx = thee->aprx;
    gm = thee->gm;

    sock = Vio_ctor(iodev,iofmt,thost,fname,"w");
    if (sock == VNULL) {
        Vnm_print(2, "Vfetk_write: Problem opening virtual socket %s\n",
          fname);
        return 0;
    }
    if (Vio_connect(sock, 0) < 0) {
        Vnm_print(2, "Vfetk_write: Problem connecting to virtual socket %s\n",
          fname);
        return 0;
    }

    /* Make sure vec has enough rows to accomodate the vertex data */
    if (Bvec_numRB(vec, 0) != Gem_numVV(gm)) {
        Vnm_print(2, "Vfetk_fillArray:  insufficient space in Bvec!\n");
        Vnm_print(2, "Vfetk_fillArray:  Have %d, need %d!\n", Bvec_numRB(vec, 0),
          Gem_numVV(gm));
        return 0;
    }

    switch (format) {

        case VDF_DX:
            Aprx_writeSOL(aprx, sock, vec, "DX");
            break;
        case VDF_AVS:
            Aprx_writeSOL(aprx, sock, vec, "UCD");
            break;
        case VDF_UHBD:
            Vnm_print(2, "Vfetk_write:  UHBD format not supported!\n");
            return 0;
        default:
            Vnm_print(2, "Vfetk_write:  Invalid data format (%d)!\n", format);
            return 0;
    }


    Vio_connectFree(sock);
    Vio_dtor(&sock);

    return 1;
}

#endif
