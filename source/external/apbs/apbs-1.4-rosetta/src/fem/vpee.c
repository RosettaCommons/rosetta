/**
 *  @file    vpee.c
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @brief   Class Vpee methods
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

#include "vpee.h"

VPRIVATE int Vpee_userDefined(Vpee *thee,
                              SS *sm
                              );
VPRIVATE int Vpee_ourSimp(Vpee *thee,
                          SS *sm,
                          int rcol
                          );
VEXTERNC double Aprx_estNonlinResid(Aprx *thee,
                                    SS *sm,
                                    Bvec *u,
                                    Bvec *ud,
                                    Bvec *f
                                    );
VEXTERNC double Aprx_estLocalProblem(Aprx *thee,
                                     SS *sm,
                                     Bvec *u,
                                     Bvec *ud,
                                     Bvec *f);
VEXTERNC double Aprx_estDualProblem(Aprx *thee,
                                    SS *sm,
                                    Bvec *u,
                                    Bvec *ud,
                                    Bvec *f
                                    );

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee: Non-inlineable methods
/////////////////////////////////////////////////////////////////////////// */

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpee_ctor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vpee* Vpee_ctor(Gem *gm,
                        int localPartID,
                        int killFlag,
                        double killParam
                        ) {

    Vpee *thee = VNULL;

    /* Set up the structure */
    thee = Vmem_malloc(VNULL, 1, sizeof(Vpee) );
    VASSERT( thee != VNULL);
    VASSERT( Vpee_ctor2(thee, gm, localPartID, killFlag, killParam));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpee_ctor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpee_ctor2(Vpee *thee,
                       Gem *gm,
                       int localPartID,
                       int killFlag,
                       double killParam
                       ) {

    int ivert,
        nLocalVerts;
    SS *simp;
    VV *vert;
    double radius,
           dx,
           dy,
           dz;

    VASSERT(thee != VNULL);

    /* Sanity check on input values */
    if (killFlag == 0) {
        Vnm_print(0, "Vpee_ctor2: No error attenuation outside partition.\n");
    } else if (killFlag == 1) {
        Vnm_print(0, "Vpee_ctor2: Error outside local partition ignored.\n");
    } else if (killFlag == 2) {
        Vnm_print(0, "Vpee_ctor2: Error ignored outside sphere with radius %4.3f times the radius of the circumscribing sphere\n", killParam);
        if (killParam < 1.0) {
          Vnm_print(2, "Vpee_ctor2: Warning! Parameter killParam = %4.3 < 1.0!\n",
            killParam);
          Vnm_print(2, "Vpee_ctor2: This may result in non-optimal marking and refinement!\n");
        }
    } else if (killFlag == 3) {
        Vnm_print(0, "Vpee_ctor2: Error outside local partition and immediate neighbors ignored [NOT IMPLEMENTED].\n");
    } else {
        Vnm_print(2, "Vpee_ctor2: UNRECOGNIZED killFlag PARAMETER! BAILING!.\n");
        VASSERT(0);
    }

    thee->gm = gm;
    thee->localPartID = localPartID;
    thee->killFlag = killFlag;
    thee->killParam = killParam;
    thee->mem = Vmem_ctor("APBS::VPEE");

    /* Now, figure out the center of geometry for the local partition.  The
     * general plan is to loop through the vertices, loop through the
     * vertices' simplex lists and find the vertices with simplices containing
     * chart values we're interested in. */
    thee->localPartCenter[0] = 0.0;
    thee->localPartCenter[1] = 0.0;
    thee->localPartCenter[2] = 0.0;
    nLocalVerts = 0;
    for (ivert=0; ivert<Gem_numVV(thee->gm); ivert++) {
        vert = Gem_VV(thee->gm, ivert);
        simp = VV_firstSS(vert);
        VASSERT(simp != VNULL);
        while (simp != VNULL) {
            if (SS_chart(simp) == thee->localPartID) {
                thee->localPartCenter[0] += VV_coord(vert, 0);
                thee->localPartCenter[1] += VV_coord(vert, 1);
                thee->localPartCenter[2] += VV_coord(vert, 2);
                nLocalVerts++;
                break;
            }
            simp = SS_link(simp, vert);
        }
    }
    VASSERT(nLocalVerts > 0);
    thee->localPartCenter[0] =
      thee->localPartCenter[0]/((double)(nLocalVerts));
    thee->localPartCenter[1] =
      thee->localPartCenter[1]/((double)(nLocalVerts));
    thee->localPartCenter[2] =
      thee->localPartCenter[2]/((double)(nLocalVerts));
    Vnm_print(0, "Vpee_ctor2: Part %d centered at (%4.3f, %4.3f, %4.3f)\n",
      thee->localPartID, thee->localPartCenter[0], thee->localPartCenter[1],
      thee->localPartCenter[2]);


    /* Now, figure out the radius of the sphere circumscribing the local
     * partition.  We need to keep track of vertices so we don't double count
     * them. */
    thee->localPartRadius = 0.0;
    for (ivert=0; ivert<Gem_numVV(thee->gm); ivert++) {
        vert = Gem_VV(thee->gm, ivert);
        simp = VV_firstSS(vert);
        VASSERT(simp != VNULL);
        while (simp != VNULL) {
            if (SS_chart(simp) == thee->localPartID) {
                dx = thee->localPartCenter[0] - VV_coord(vert, 0);
                dy = thee->localPartCenter[1] - VV_coord(vert, 1);
                dz = thee->localPartCenter[2] - VV_coord(vert, 2);
                radius = dx*dx + dy*dy + dz*dz;
                if (radius > thee->localPartRadius) thee->localPartRadius =
                  radius;
                break;
            }
            simp = SS_link(simp, vert);
        }
    }
    thee->localPartRadius = VSQRT(thee->localPartRadius);
    Vnm_print(0, "Vpee_ctor2: Part %d has circumscribing sphere of radius %4.3f\n",
      thee->localPartID, thee->localPartRadius);

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpee_dtor
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpee_dtor(Vpee **thee) {

    if ((*thee) != VNULL) {
        Vpee_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vpee), (void **)thee);
        (*thee) = VNULL;
    }

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpee_dtor2
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vpee_dtor2(Vpee *thee) {
    Vmem_dtor(&(thee->mem));
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpee_markRefine
//
// Author:   Nathan Baker (and Michael Holst: the author of AM_markRefine, on
//           which this is based)
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpee_markRefine(Vpee *thee,
                            AM *am,
                            int level,
                            int akey,
                            int rcol,
                            double etol,
                            int bkey
                            ) {

    Aprx *aprx;
    int marked = 0,
        markMe,
        i,
        smid,
        count,
        currentQ;
    double minError = 0.0,
           maxError = 0.0,
           errEst = 0.0,
           mlevel,
           barrier;
    SS *sm;


    VASSERT(thee != VNULL);

    /* Get the Aprx object from AM */
    aprx = am->aprx;

    /* input check and some i/o */
    if ( ! ((-1 <= akey) && (akey <= 4)) ) {
        Vnm_print(0,"Vpee_markRefine: bad refine key; simplices marked = %d\n",
            marked);
        return marked;
    }

    /* For uniform markings, we have no effect */
    if ((-1 <= akey) && (akey <= 0)) {
        marked = Gem_markRefine(thee->gm, akey, rcol);
        return marked;
    }

    /* Informative I/O */
    if (akey == 2) {
        Vnm_print(0,"Vpee_estRefine: using Aprx_estNonlinResid().\n");
    } else if (akey == 3) {
        Vnm_print(0,"Vpee_estRefine: using Aprx_estLocalProblem().\n");
    } else if (akey == 4) {
        Vnm_print(0,"Vpee_estRefine: using Aprx_estDualProblem().\n");
    } else {
        Vnm_print(0,"Vpee_estRefine: bad key given; simplices marked = %d\n",
            marked);
        return marked;
    }
    if (thee->killFlag == 0) {
        Vnm_print(0, "Vpee_markRefine: No error attenuation -- simplices in all partitions will be marked.\n");
    } else if (thee->killFlag == 1) {
        Vnm_print(0, "Vpee_markRefine: Maximum error attenuation -- only simplices in local partition will be marked.\n");
    } else if (thee->killFlag == 2) {
        Vnm_print(0, "Vpee_markRefine: Spherical error attenutation -- simplices within a sphere of %4.3f times the size of the partition will be marked\n",
          thee->killParam);
    } else if (thee->killFlag == 2) {
        Vnm_print(0, "Vpee_markRefine: Neighbor-based error attenuation -- simplices in the local and neighboring partitions will be marked [NOT IMPLEMENTED]!\n");
        VASSERT(0);
    } else {
        Vnm_print(2,"Vpee_markRefine: bogus killFlag given; simplices marked = %d\n",
            marked);
        return marked;
    }

    /* set the barrier type */
    mlevel = (etol*etol) / Gem_numSS(thee->gm);
    if (bkey == 0) {
        barrier = (etol*etol);
        Vnm_print(0,"Vpee_estRefine: forcing [err per S] < [TOL] = %g\n",
            barrier);
    } else if (bkey == 1) {
        barrier = mlevel;
        Vnm_print(0,"Vpee_estRefine: forcing [err per S] < [(TOL^2/numS)^{1/2}] = %g\n",
            VSQRT(barrier));
    } else {
        Vnm_print(0,"Vpee_estRefine: bad bkey given; simplices marked = %d\n",
            marked);
        return marked;
    }

    /* timer */
    Vnm_tstart(30, "error estimation");

    /* count = num generations to produce from marked simplices (minimally) */
    count = 1; /* must be >= 1 */

    /* check the refinement Q for emptyness */
    currentQ = 0;
    if (Gem_numSQ(thee->gm,currentQ) > 0) {
        Vnm_print(0,"Vpee_markRefine: non-empty refinement Q%d....clearing..",
            currentQ);
        Gem_resetSQ(thee->gm,currentQ);
        Vnm_print(0,"..done.\n");
    }
    if (Gem_numSQ(thee->gm,!currentQ) > 0) {
        Vnm_print(0,"Vpee_markRefine: non-empty refinement Q%d....clearing..",
            !currentQ);
        Gem_resetSQ(thee->gm,!currentQ);
        Vnm_print(0,"..done.\n");
    }
    VASSERT( Gem_numSQ(thee->gm,currentQ)  == 0 );
    VASSERT( Gem_numSQ(thee->gm,!currentQ) == 0 );

    /* clear everyone's refinement flags */
    Vnm_print(0,"Vpee_markRefine: clearing all simplex refinement flags..");
    for (i=0; i<Gem_numSS(thee->gm); i++) {
        if ( (i>0) && (i % VPRTKEY) == 0 ) Vnm_print(0,"[MS:%d]",i);
        sm = Gem_SS(thee->gm,i);
        SS_setRefineKey(sm,currentQ,0);
        SS_setRefineKey(sm,!currentQ,0);
        SS_setRefinementCount(sm,0);
    }
    Vnm_print(0,"..done.\n");

    /* NON-ERROR-BASED METHODS */
    /* Simplex flag clearing */
    if (akey == -1) return marked;
    /* Uniform & user-defined refinement*/
    if ((akey == 0) || (akey == 1)) {
        smid = 0;
        while ( smid < Gem_numSS(thee->gm)) {
            /* Get the simplex and find out if it's markable */
            sm = Gem_SS(thee->gm,smid);
            markMe = Vpee_ourSimp(thee, sm, rcol);
            if (markMe) {
                if (akey == 0) {
                    marked++;
                    Gem_appendSQ(thee->gm,currentQ, sm);
                    SS_setRefineKey(sm,currentQ,1);
                    SS_setRefinementCount(sm,count);
                } else if (Vpee_userDefined(thee, sm)) {
                    marked++;
                    Gem_appendSQ(thee->gm,currentQ, sm);
                    SS_setRefineKey(sm,currentQ,1);
                    SS_setRefinementCount(sm,count);
                }
            }
            smid++;
        }
    }

    /* ERROR-BASED METHODS */
    /* gerror = global error accumulation */
    aprx->gerror = 0.;

    /* traverse the simplices and process the error estimates */
    Vnm_print(0,"Vpee_markRefine: estimating error..");
    smid = 0;
    while ( smid < Gem_numSS(thee->gm)) {

        /* Get the simplex and find out if it's markable */
        sm = Gem_SS(thee->gm,smid);
        markMe = Vpee_ourSimp(thee, sm, rcol);

        if ( (smid>0) && (smid % VPRTKEY) == 0 ) Vnm_print(0,"[MS:%d]",smid);

        /* Produce an error estimate for this element if it is in the set */
        if (markMe) {
            if (akey == 2) {
                errEst = Aprx_estNonlinResid(aprx, sm, am->u,am->ud,am->f);
            } else if (akey == 3) {
                errEst = Aprx_estLocalProblem(aprx, sm, am->u,am->ud,am->f);
            } else if (akey == 4) {
                errEst = Aprx_estDualProblem(aprx, sm, am->u,am->ud,am->f);
            }
            VASSERT( errEst >= 0. );

            /* if error estimate above tol, mark element for refinement */
            if ( errEst > barrier ) {
                marked++;
                Gem_appendSQ(thee->gm,currentQ, sm); /*add to refinement Q*/
                SS_setRefineKey(sm,currentQ,1);      /* note now on refine Q */
                SS_setRefinementCount(sm,count);     /* refine X many times? */
            }

            /* keep track of min/max errors over the mesh */
            minError = VMIN2( VSQRT(VABS(errEst)), minError );
            maxError = VMAX2( VSQRT(VABS(errEst)), maxError );

            /* store the estimate */
            Bvec_set( aprx->wev, smid, errEst );

            /* accumlate into global error (errEst is SQUAREd already) */
            aprx->gerror += errEst;

        /* otherwise store a zero for the estimate */
        } else {
            Bvec_set( aprx->wev, smid, 0. );
        }

        smid++;
    }

    /* do some i/o */
    Vnm_print(0,"..done.  [marked=<%d/%d>]\n",marked,Gem_numSS(thee->gm));
    Vnm_print(0,"Vpee_estRefine: TOL=<%g>  Global_Error=<%g>\n",
        etol, aprx->gerror);
    Vnm_print(0,"Vpee_estRefine: (TOL^2/numS)^{1/2}=<%g>  Max_Ele_Error=<%g>\n",
        VSQRT(mlevel),maxError);
    Vnm_tstop(30, "error estimation");

    /* check for making the error tolerance */
    if ((bkey == 1) && (aprx->gerror <= etol)) {
        Vnm_print(0,
            "Vpee_estRefine: *********************************************\n");
        Vnm_print(0,
            "Vpee_estRefine: Global Error criterion met; setting marked=0.\n");
        Vnm_print(0,
            "Vpee_estRefine: *********************************************\n");
        marked = 0;
    }


    /* return */
    return marked;

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpee_numSS
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vpee_numSS(Vpee *thee) {
    int num = 0;
    int isimp;

    for (isimp=0; isimp<Gem_numSS(thee->gm); isimp++) {
        if (SS_chart(Gem_SS(thee->gm, isimp)) == thee->localPartID) num++;
    }

    return num;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpee_userDefined
//
// Purpose:  Reduce code bloat by wrapping up the common steps for getting the
//           user-defined error estimate
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int Vpee_userDefined(Vpee *thee,
                              SS *sm
                              ) {

    int ivert,
        icoord,
        chart[4],
        fType[4],
        vType[4];
    double vx[4][3];

    for (ivert=0; ivert<Gem_dimVV(thee->gm); ivert++) {
        fType[ivert] = SS_faceType(sm,ivert);
        vType[ivert] = VV_type(SS_vertex(sm,ivert) );
        chart[ivert] = VV_chart(SS_vertex(sm,ivert) );
        for (icoord=0; icoord<Gem_dimII(thee->gm); icoord++) {
            vx[ivert][icoord] = VV_coord(SS_vertex(sm,ivert), icoord );
        }
    }
    return thee->gm->pde->markSimplex(Gem_dim(thee->gm), Gem_dimII(thee->gm),
             SS_type(sm), fType, vType, chart, vx, sm);
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vpee_ourSimp
//
// Purpose:  Reduce code bloat by wrapping up the common steps for determining
//           whether the given simplex can be marked (i.e., belongs to our
//           partition or overlap region)
//
// Returns:  1 if could be marked, 0 otherwise
//
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPRIVATE int Vpee_ourSimp(Vpee *thee,
                          SS *sm,
                          int rcol
                          ) {

    int ivert;
    double dist,
           dx,
           dy,
           dz;

    if (thee->killFlag == 0) return 1;
    else if (thee->killFlag == 1) {
        if ((SS_chart(sm) == rcol) || (rcol < 0)) return 1;
    } else if (thee->killFlag == 2) {
        if (rcol < 0) return 1;
        else {
            /* We can only do distance-based searches on the local partition */
            VASSERT(rcol == thee->localPartID);
            /* Find the closest distance between this simplex and the
             * center of the local partition and check it against
             * (thee->localPartRadius*thee->killParam) */
            dist = 0;
            for (ivert=0; ivert<SS_dimVV(sm); ivert++) {
                dx = VV_coord(SS_vertex(sm, ivert), 0) -
                  thee->localPartCenter[0];
                dy = VV_coord(SS_vertex(sm, ivert), 1) -
                  thee->localPartCenter[1];
                dz = VV_coord(SS_vertex(sm, ivert), 2) -
                  thee->localPartCenter[2];
                dist = VSQRT((dx*dx + dy*dy + dz*dz));
            }
            if (dist < thee->localPartRadius*thee->killParam) return 1;
        }
    } else if (thee->killFlag == 3) VASSERT(0);
    else VASSERT(0);

    return 0;

}
