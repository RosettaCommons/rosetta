/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief
 *  @version $Id:
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (nathan.baker@pnl.gov)
 * Pacific Northwest National Laboratory
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2012 Battelle Memorial Institute. Developed at the Pacific Northwest National Laboratory, operated by Battelle Memorial Institute, Pacific Northwest Division for the U.S. Department Energy.  Portions Copyright (c) 2002-2010, Washington University in St. Louis.  Portions Copyright (c) 2002-2010, Nathan A. Baker.  Portions Copyright (c) 1999-2002, The Regents of the University of California. Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
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
 * @endverbatim
 */

#include "mgsubd.h"

VPUBLIC void Vbuildops(
        int *nx, int *ny, int *nz,
        int *nlev, int *ipkey, int *iinfo,
        int *ido, int *iz,
        int *mgprol, int *mgcoar, int *mgsolv, int *mgdisc, int *ipc,
        double *rpc, double *pc, double *ac, double *cc, double *fc,
        double *xf, double *yf, double *zf,
        double *gxcf, double *gycf, double *gzcf,
        double *a1cf, double *a2cf, double *a3cf,
        double *ccf, double *fcf, double *tcf
        ) {

    // @todo Document this function
    int lev = 0;
    int nxx = 0;
    int nyy = 0;
    int nzz = 0;
    int nxold = 0;
    int nyold = 0;
    int nzold = 0;
    int numdia = 0;
    int key = 0;

    // Utility variables
    int i;

    MAT2(iz, 50, *nlev);

    // Setup
    nxx = *nx;
    nyy = *ny;
    nzz = *nz;

    // Build the operator a on the finest level
    if (*ido == 0 || *ido == 2) {

        lev = 1;

        // Some i/o
        if (*iinfo > 0)
            VMESSAGE3("Fine: (%03d, %03d, %03d)", nxx, nyy, nzz);

        // Finest level discretization
        VbuildA(&nxx, &nyy, &nzz,
                ipkey, mgdisc, &numdia,
                 RAT(ipc, VAT2(iz, 5,lev)),  RAT(rpc, VAT2(iz, 6,lev)),
                  RAT(ac, VAT2(iz, 7,lev)),   RAT(cc, VAT2(iz, 1,lev)),   RAT(fc, VAT2(iz,  1,lev)),
                  RAT(xf, VAT2(iz, 8,lev)),   RAT(yf, VAT2(iz, 9,lev)),   RAT(zf, VAT2(iz, 10,lev)),
                RAT(gxcf, VAT2(iz, 2,lev)), RAT(gycf, VAT2(iz, 3,lev)), RAT(gzcf, VAT2(iz,  4,lev)),
                RAT(a1cf, VAT2(iz, 1,lev)), RAT(a2cf, VAT2(iz, 1,lev)), RAT(a3cf, VAT2(iz,  1,lev)),
                 RAT(ccf, VAT2(iz, 1,lev)),  RAT(fcf, VAT2(iz, 1,lev)));

        VMESSAGE2("Operator stencil (lev, numdia) = (%d, %d)", lev, numdia);

        // Now initialize the differential operator offset
        VAT2(iz, 7, lev+1) = VAT2(iz, 7, lev) + numdia * nxx * nyy * nzz;

        // Debug
        if (*iinfo > 7) {
            Vprtmatd(&nxx, &nyy, &nzz,
                    RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)), RAT(ac, VAT2(iz, 7,lev)));
        }
    }

    // Build the (nlev-1) level operators
    if (*ido == 1 || *ido == 2 || *ido == 3) {

        for (lev=2; lev<=*nlev; lev++) {
            nxold = nxx;
            nyold = nyy;
            nzold = nzz;
            i = 1;


            Vmkcors(&i, &nxold, &nyold, &nzold, &nxx, &nyy, &nzz);
            if (*ido != 3) {

                // Build the interpolation operator on this level
                VbuildP(&nxold, &nyold, &nzold,
                        &nxx, &nyy, &nzz,
                        mgprol,
                        RAT(ipc, VAT2(iz,  5,lev-1)), RAT(rpc, VAT2(iz, 6,lev-1)),
                         RAT(pc, VAT2(iz, 11,lev-1)),  RAT(ac, VAT2(iz, 7,lev-1)),
                         RAT(xf, VAT2(iz,  8,lev-1)),  RAT(yf, VAT2(iz, 9,lev-1)), RAT(zf, VAT2(iz, 10,lev-1)));

                // Differential operator this level with standard disc.
                if (*mgcoar == 0) {

                    // Some i/o
                    if (*iinfo > 0)
                        VMESSAGE3("Stand: (%03d, %03d, %03d)", nxx, nyy, nzz);



                    Vbuildcopy0(&nxx, &nyy, &nzz,
                            &nxold, &nyold, &nzold,
                              RAT(xf, VAT2(iz, 8,lev  )),   RAT(yf, VAT2(iz, 9,lev  )),   RAT(zf, VAT2(iz, 10,lev  )),
                            RAT(gxcf, VAT2(iz, 2,lev  )), RAT(gycf, VAT2(iz, 3,lev  )), RAT(gzcf, VAT2(iz,  4,lev  )),
                            RAT(a1cf, VAT2(iz, 1,lev  )), RAT(a2cf, VAT2(iz, 1,lev  )), RAT(a3cf, VAT2(iz,  1,lev  )),
                             RAT(ccf, VAT2(iz, 1,lev  )),  RAT(fcf, VAT2(iz, 1,lev  )),  RAT(tcf, VAT2(iz,  1,lev  )),
                              RAT(xf, VAT2(iz, 8,lev-1)),   RAT(yf, VAT2(iz, 9,lev-1)),   RAT(zf, VAT2(iz, 10,lev-1)),
                            RAT(gxcf, VAT2(iz, 2,lev-1)), RAT(gycf, VAT2(iz, 3,lev-1)), RAT(gzcf, VAT2(iz,  4,lev-1)),
                            RAT(a1cf, VAT2(iz, 1,lev-1)), RAT(a2cf, VAT2(iz, 1,lev-1)), RAT(a3cf, VAT2(iz,  1,lev-1)),
                             RAT(ccf, VAT2(iz, 1,lev-1)),  RAT(fcf, VAT2(iz, 1,lev-1)),  RAT(tcf, VAT2(iz,  1,lev-1)));

                    VbuildA(&nxx, &nyy, &nzz,
                            ipkey, mgdisc, &numdia,
                             RAT(ipc, VAT2(iz, 5,lev)),  RAT(rpc, VAT2(iz, 6,lev)),
                              RAT(ac, VAT2(iz, 7,lev)),   RAT(cc, VAT2(iz, 1,lev)),   RAT(fc, VAT2(iz,  1,lev)),
                              RAT(xf, VAT2(iz, 8,lev)),   RAT(yf, VAT2(iz, 9,lev)),   RAT(zf, VAT2(iz, 10,lev)),
                            RAT(gxcf, VAT2(iz, 2,lev)), RAT(gycf, VAT2(iz, 3,lev)), RAT(gzcf, VAT2(iz,  4,lev)),
                            RAT(a1cf, VAT2(iz, 1,lev)), RAT(a2cf, VAT2(iz, 1,lev)), RAT(a3cf, VAT2(iz,  1,lev)),
                             RAT(ccf, VAT2(iz, 1,lev)),  RAT(fcf, VAT2(iz, 1,lev)));
                }

                // Differential operator this level with harmonic disc.
                else if (*mgcoar == 1) {

                    // Some i/o
                    if (*iinfo > 0)
                        VMESSAGE3("Harm0: (%03d, %03d, %03d)", nxx, nyy, nzz);

                    Vbuildharm0(&nxx, &nyy, &nzz, &nxold, &nyold, &nzold,
                              RAT(xf, VAT2(iz, 8, lev  )),   RAT(yf, VAT2(iz, 9, lev  )),   RAT(zf, VAT2(iz, 10, lev  )),
                            RAT(gxcf, VAT2(iz, 2, lev  )), RAT(gycf, VAT2(iz, 3, lev  )), RAT(gzcf, VAT2(iz,  4, lev  )),
                            RAT(a1cf, VAT2(iz, 1, lev  )), RAT(a2cf, VAT2(iz, 1, lev  )), RAT(a3cf, VAT2(iz,  1, lev  )),
                             RAT(ccf, VAT2(iz, 1, lev  )),  RAT(fcf, VAT2(iz, 1, lev  )),  RAT(tcf, VAT2(iz,  1, lev  )),
                              RAT(xf, VAT2(iz, 8, lev-1)),   RAT(yf, VAT2(iz, 9, lev-1)),   RAT(zf, VAT2(iz, 10, lev-1)),
                            RAT(gxcf, VAT2(iz, 2, lev-1)), RAT(gycf, VAT2(iz, 3, lev-1)), RAT(gzcf, VAT2(iz,  4, lev-1)),
                            RAT(a1cf, VAT2(iz, 1, lev-1)), RAT(a2cf, VAT2(iz, 1, lev-1)), RAT(a3cf, VAT2(iz,  1, lev-1)),
                             RAT(ccf, VAT2(iz, 1, lev-1)),  RAT(fcf, VAT2(iz, 1, lev-1)),  RAT(tcf, VAT2(iz,  1, lev-1)));

                    VbuildA(&nxx, &nyy, &nzz,
                            ipkey, mgdisc, &numdia,
                             RAT(ipc, VAT2(iz, 5,lev)),  RAT(rpc, VAT2(iz, 6,lev)),
                              RAT(ac, VAT2(iz, 7,lev)),   RAT(cc, VAT2(iz, 1,lev)),   RAT(fc, VAT2(iz,  1,lev)),
                              RAT(xf, VAT2(iz, 8,lev)),   RAT(yf, VAT2(iz, 9,lev)),   RAT(zf, VAT2(iz, 10,lev)),
                            RAT(gxcf, VAT2(iz, 2,lev)), RAT(gycf, VAT2(iz, 3,lev)), RAT(gzcf, VAT2(iz,  4,lev)),
                            RAT(a1cf, VAT2(iz, 1,lev)), RAT(a2cf, VAT2(iz, 1,lev)), RAT(a3cf, VAT2(iz,  1,lev)),
                             RAT(ccf, VAT2(iz, 1,lev)),  RAT(fcf, VAT2(iz, 1,lev)));
                }

                // Differential operator with galerkin formulation ***
                else if (*mgcoar == 2) {

                    // Some i/o
                    if (*iinfo > 0)
                        VMESSAGE3("Galer: (%03d, %03d, %03d)", nxx, nyy, nzz);

                    Vbuildgaler0(&nxold, &nyold, &nzold,
                            &nxx, &nyy, &nzz,
                            ipkey, &numdia,
                             RAT(pc, VAT2(iz, 11,lev-1)),
                            RAT(ipc, VAT2(iz,  5,lev-1)), RAT(rpc, VAT2(iz, 6,lev-1)),
                             RAT(ac, VAT2(iz,  7,lev-1)),  RAT(cc, VAT2(iz, 1,lev-1)), RAT(fc, VAT2(iz, 1,lev-1)),
                            RAT(ipc, VAT2(iz,  5,lev  )), RAT(rpc, VAT2(iz, 6,lev  )),
                             RAT(ac, VAT2(iz,  7,lev  )),  RAT(cc, VAT2(iz, 1,lev  )), RAT(fc, VAT2(iz, 1,lev  )));



                    Vextrac(&nxold, &nyold, &nzold,
                            &nxx, &nyy, &nzz,
                            RAT(tcf, VAT2(iz, 1,lev-1)), RAT(tcf, VAT2(iz, 1,lev)));
                }
                else {
                    VABORT_MSG1("Bad mgcoar value given: %d", *mgcoar);
                }

                // Now initialize the differential operator offset
                // Vnm_print(0, "BUILDOPS: operator stencil (lev,numdia) = (%d, %d)\n",
                //		lev,numdia);
                VAT2(iz, 7, lev+1) = VAT2(iz, 7,lev) + numdia * nxx * nyy * nzz;

                // Debug
                if (*iinfo > 8) {

                    Vprtmatd(&nxx, &nyy, &nzz,
                            RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)), RAT(ac, VAT2(iz, 7,lev)));
                }
            }
        }

        // Build a sparse format coarse grid operator
        if (*mgsolv == 1) {
            lev = *nlev;

            Vbuildband(&key, &nxx, &nyy, &nzz,
                    RAT(ipc, VAT2(iz, 5,lev  )), RAT(rpc, VAT2(iz, 6,lev  )), RAT(ac, VAT2(iz, 7,lev  )),
                    RAT(ipc, VAT2(iz, 5,lev+1)), RAT(rpc, VAT2(iz, 6,lev+1)), RAT(ac, VAT2(iz, 7,lev+1)));

            if (key == 1) {
                VERRMSG0("Changing your mgsolv to iterative");
                *mgsolv = 0;
            }
        }
    }
}



VPUBLIC void Vbuildstr(int *nx, int *ny, int *nz, int *nlev, int *iz) {

    int nxold, nyold, nzold;
    int nxnew, nynew, nznew;
    int n, lev;

    // Utility variable
    int numlev;

    MAT2(iz, 50, *nlev);

    // Setup
    nxnew  = *nx;
    nynew  = *ny;
    nznew  = *nz;
    n      = nxnew * nynew * nznew;

    // Start with level 1
    lev = 1;

    // Mark beginning of everything at level 1 ***
    VAT2(iz,  1, lev) = 1;
    VAT2(iz,  2, lev) = 1;
    VAT2(iz,  3, lev) = 1;
    VAT2(iz,  4, lev) = 1;
    VAT2(iz,  5, lev) = 1;
    VAT2(iz,  6, lev) = 1;
    VAT2(iz,  7, lev) = 1;
    VAT2(iz,  8, lev) = 1;
    VAT2(iz,  9, lev) = 1;
    VAT2(iz, 10, lev) = 1;
    VAT2(iz, 11, lev) = 1;

    // Mark beginning of everything at level 2
    VAT2(iz,  1, lev + 1) = VAT2(iz,  1, lev) + n;
    VAT2(iz,  2, lev + 1) = VAT2(iz,  2, lev) + 4 * nynew * nznew;
    VAT2(iz,  3, lev + 1) = VAT2(iz,  3, lev) + 4 * nxnew * nznew;
    VAT2(iz,  4, lev + 1) = VAT2(iz,  4, lev) + 4 * nxnew * nynew;
    VAT2(iz,  5, lev + 1) = VAT2(iz,  5, lev) + 100;
    VAT2(iz,  6, lev + 1) = VAT2(iz,  6, lev) + 100;
    VAT2(iz,  8, lev + 1) = VAT2(iz,  8, lev) + nxnew;
    VAT2(iz,  9, lev + 1) = VAT2(iz,  9, lev) + nynew;
    VAT2(iz, 10, lev + 1) = VAT2(iz, 10, lev) + nznew;

    /* ******************************************************************
     * ***NOTE: we mark operator offsets as we build the operators    ***
     * ***VAT2(iz, 7,lev+1)  = VAT2(iz, 7, lev)  + 4*n              ***
     * ******************************************************************
     * ***NOTE: we mark prolongation operator offsets lagging a level ***
     * ***VAT2(iz, 11, lev)   = VAT2(iz, 11,lev-1) + 27*nsmall      ***
     * ******************************************************************/

    // Mark the beginning of everything at (nlev-1) more
    for (lev=2; lev<=*nlev; lev++) {
       nxold = nxnew;
       nyold = nynew;
       nzold = nznew;
       numlev = 1;
       Vmkcors(&numlev, &nxold, &nyold, &nzold, &nxnew, &nynew, &nznew);
       n = nxnew * nynew * nznew;

       // Mark the beginning of everything at level (lev+1)
       VAT2(iz,  1, lev + 1) = VAT2(iz,  1, lev) + n;
       VAT2(iz,  2, lev + 1) = VAT2(iz,  2, lev) + 4 * nynew * nznew;
       VAT2(iz,  3, lev + 1) = VAT2(iz,  3, lev) + 4 * nxnew * nznew;
       VAT2(iz,  4, lev + 1) = VAT2(iz,  4, lev) + 4 * nxnew * nynew;
       VAT2(iz,  5, lev + 1) = VAT2(iz,  5, lev) + 100;
       VAT2(iz,  6, lev + 1) = VAT2(iz,  6, lev) + 100;
       VAT2(iz,  7, lev + 1) = VAT2(iz,  7, lev) + 4 * n;
       VAT2(iz,  8, lev + 1) = VAT2(iz,  8, lev) + nxnew;
       VAT2(iz,  9, lev + 1) = VAT2(iz,  9, lev) + nynew;
       VAT2(iz, 10, lev + 1) = VAT2(iz, 10, lev) + nznew;

       // Mark prolongation operator storage for previous level
       VAT2(iz, 11, lev) = VAT2(iz, 11, lev - 1) + 27 * n;

       /* ****************************************************************
        * *** NOTE: we mark operator offsets as we build the operators ***
        * ***       VAT2(iz, 7, lev + 1) = VAT2(iz, 7, lev) + 4 * n  ***
        ******************************************************************/
    }
}


VPUBLIC void Vbuildgaler0(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        int *ipkey, int *numdia,
        double *pcFF, int   *ipcFF, double *rpcFF,
        double *acFF, double *ccFF, double *fcFF,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc) {

    int numdia_loc;

    // Call the algebraic galerkin routine
    numdia_loc = VAT(ipcFF, 11);
    VbuildG(nxf, nyf, nzf,
            nxc, nyc, nzc,
            &numdia_loc,
            pcFF, acFF, ac);

    // Note how many nonzeros in this new discretization stencil
    VAT(ipc, 11) = 27;
    *numdia = 14;

    // Save the problem key with this new operator
    VAT(ipc, 10) = *ipkey;

    // Restrict the helmholtz term and source function
    Vrestrc(nxf, nyf, nzf,
            nxc, nyc, nzc,
            ccFF, cc, pcFF);

    Vrestrc(nxf, nyf, nzf,
            nxc, nyc, nzc,
            fcFF, fc, pcFF);
}



VPUBLIC void Vmkcors(int *numlev,
        int *nxold, int *nyold, int *nzold,
        int *nxnew, int *nynew, int *nznew) {
    int nxtmp, nytmp, nztmp; // Temporary variables to hold current x,y,z values
    int i;                   // Index used in for loops

    // Determine the coarser grid
    *nxnew = *nxold;
    *nynew = *nyold;
    *nznew = *nzold;

    for (i=1; i<=*numlev; i++) {
       nxtmp = *nxnew;
       nytmp = *nynew;
       nztmp = *nznew;
       Vcorsr(&nxtmp,nxnew);
       Vcorsr(&nytmp,nynew);
       Vcorsr(&nztmp,nznew);
    }
}



VPUBLIC void Vcorsr(int *nold, int *nnew) {

    // Find the coarser grid size ***
    *nnew = (*nold - 1) / 2 + 1;

    // Check a few things
    if ((*nnew - 1) * 2 != *nold - 1) {
        Vnm_print(2, "Vcorsr:  WARNING!  The grid dimensions are not\n");
        Vnm_print(2, "Vcorsr:  consistent with nlev!  Your\n");
        Vnm_print(2, "Vcorsr:  calculation will only work if you\n");
        Vnm_print(2, "Vcorsr:  are performing a mg-dummy run.\n");
    }
    if (*nnew < 1) {
        /// @todo:  Determine if this should abort at this point
        Vnm_print(2, "Vcorsr:  ERROR!  The grid dimensions are not\n");
        Vnm_print(2, "Vcorsr:  consistent with nlev!\n");
        Vnm_print(2, "Vcorsr:  Grid coarsened below zero.\n");
    }
}



VPUBLIC void Vmkfine(int *numlev,
        int *nxold, int *nyold, int *nzold,
        int *nxnew, int *nynew, int *nznew) {

    int nxtmp, nytmp, nztmp, i;

    // Determine the finer grid
    *nxnew = *nxold;
    *nynew = *nyold;
    *nznew = *nzold;

    for (i=1; i<=*numlev; i++) {

        nxtmp = *nxnew;
        nytmp = *nynew;
        nztmp = *nznew;

        Vfiner(&nxtmp, nxnew);
        Vfiner(&nytmp, nynew);
        Vfiner(&nztmp, nznew);
    }

}



VPUBLIC void Vfiner(int *nold, int *nnew) {
    // Find the coarser grid size ***
    *nnew = (*nold - 1) * 2 + 1;
}



VPUBLIC int Vivariv(int *nu, int *level) {

    /// @todo:  Determine if the Variable V-cycle will ever be used
    int ivariv;

    // Variable V-cycle
    // ivariv = *nu * VPOW(2, level - 1);

    // Standard V-cycle
    ivariv = *nu;

    return ivariv;
}



VPUBLIC int Vmaxlev(int n1, int n2, int n3) {

    int n1c;
    int n2c;
    int n3c;
    int lev;
    int iden;
    int idone;

    // Fine the common level
    idone = 0;
    lev = 0;
    do {
        lev += 1;
        iden = (int)VPOW(2, lev - 1);
        n1c = (n1 - 1) / iden + 1;
        n2c = (n2 - 1) / iden + 1;
        n3c = (n3 - 1) / iden + 1;
        if (((n1c - 1) * iden != (n1 - 1)) || (n1c <= 2) )
            idone = 1;
        if (((n2c - 1) * iden != (n2 - 1)) || (n2c <= 2) )
            idone = 1;
        if (((n3c - 1) * iden != (n3 - 1)) || (n3c <= 2) )
            idone = 1;

    } while (idone != 1);

    return lev - 1;
}



VPUBLIC void Vprtstp(int iok, int iters,
        double rsnrm, double rsden, double orsnrm) {

    double relres = 0.0;
    double contrac = 0.0;

    // Initializing timer
    if (iters == -99) {
        // Vnm_tstart(40, "MG iteration");
        cputme = 0.0;
        return;
    }

    // Setup for the iteration
    else if (iters == -1) {
        Vnm_tstop(40, "MG iteration");
        return;
    }

    // During the iteration
    else {

        // Stop the timer
        // Vnm_tstop(40, "MG iteration");

        // Relative residual
        if (rsden == 0.0) {
            relres = 1.0e6;
            VERRMSG0("Vprtstp: avoided division by zero\n");
        } else {
            relres = rsnrm / rsden;
        }

        // Contraction number
        if (orsnrm == 0.0) {
            contrac = 1.0e6;
            VERRMSG0("avoided division by zero\n");
        } else {
            contrac = rsnrm / orsnrm;
        }

        // The i/o
        if (iok == 1 || iok == 2) {
            VMESSAGE1("iteration = %d", iters);
            VMESSAGE1("relative residual = %e", relres);
            VMESSAGE1("contraction number = %e", contrac);
        }
    }
}



VPUBLIC void Vpackmg(int *iparm, double *rparm, int *nrwk, int *niwk,
        int *nx, int *ny, int *nz, int *nlev, int *nu1, int *nu2, int *mgkey,
        int *itmax, int *istop, int *ipcon, int *nonlin, int *mgsmoo, int *mgprol,
        int *mgcoar, int *mgsolv, int *mgdisc, int *iinfo, double *errtol,
        int *ipkey, double *omegal, double *omegan, int *irite, int *iperf) {

    /// @todo  Convert this into a struct

    // Encode iparm parameters ***
    VAT(iparm, 1)  = *nrwk;
    VAT(iparm, 2)  = *niwk;
    VAT(iparm, 3)  = *nx;
    VAT(iparm, 4)  = *ny;
    VAT(iparm, 5)  = *nz;
    VAT(iparm, 6)  = *nlev;
    VAT(iparm, 7)  = *nu1;
    VAT(iparm, 8)  = *nu2;
    VAT(iparm, 9)  = *mgkey;
    VAT(iparm, 10) = *itmax;
    VAT(iparm, 11) = *istop;
    VAT(iparm, 12) = *iinfo;
    VAT(iparm, 13) = *irite;
    VAT(iparm, 14) = *ipkey;
    VAT(iparm, 15) = *ipcon;
    VAT(iparm, 16) = *nonlin;
    VAT(iparm, 17) = *mgprol;
    VAT(iparm, 18) = *mgcoar;
    VAT(iparm, 19) = *mgdisc;
    VAT(iparm, 20) = *mgsmoo;
    VAT(iparm, 21) = *mgsolv;
    VAT(iparm, 22) = *iperf;

    // Encode rparm parameters
    VAT(rparm, 1)  = *errtol;
    VAT(rparm, 9)  = *omegal;
    VAT(rparm, 10) = *omegan;
}



VEXTERNC void Vbuildcopy0(int *nx, int *ny, int *nz,
        int *nxf, int *nyf, int *nzf,
        double *xc, double *yc, double *zc,
        double *gxc, double *gyc, double *gzc,
        double *a1c, double *a2c, double *a3c,
        double *cc, double *fc, double *tc,
        double *xf, double *yf, double *zf,
        double *gxcf, double *gycf, double *gzcf,
        double *a1cf, double *a2cf, double *a3cf,
        double *ccf, double *fcf, double *tcf) {


    int    i,    j,    k;
    int   ii,   jj,   kk;
    int iadd, jadd, kadd;

    MAT3( gxc,  *ny,  *nz,    4);
    MAT3( gyc,  *nx,  *nz,    4);
    MAT3( gzc,  *nx,  *ny,    4);
    MAT3( a1c,  *nx,  *ny,  *nz);
    MAT3( a2c,  *nx,  *ny,  *nz);
    MAT3( a3c,  *nx,  *ny,  *nz);
    MAT3(  cc,  *nx,  *ny,  *nz);
    MAT3(  fc,  *nx,  *ny,  *nz);
    MAT3(  tc,  *nx,  *ny,  *nz);
    MAT3(gxcf, *nyf, *nzf,    4);
    MAT3(gycf, *nxf, *nzf,    4);
    MAT3(gzcf, *nxf, *nyf,    4);
    MAT3(a1cf, *nxf, *nyf, *nzf);
    MAT3(a2cf, *nxf, *nyf, *nzf);
    MAT3(a3cf, *nxf, *nyf, *nzf);
    MAT3( tcf, *nxf, *nyf, *nzf);
    MAT3( ccf, *nxf, *nyf, *nzf);
    MAT3( fcf, *nxf, *nyf, *nzf);

    WARN_UNTESTED;

    // How far to step into the coefficient arrays
    iadd   = (*nxf - 1) / (*nx - 1);
    jadd   = (*nyf - 1) / (*ny - 1);
    kadd   = (*nzf - 1) / (*nz - 1);

    if (iadd != 2 || jadd != 2 || kadd != 2) {
        Vnm_print(2, "Vbuildcopy0:  Problem with grid dimensions...\n");
    }

    // Compute the coarse grid pde coefficients
    for (k=1; k<=*nz; k++) {
        kk = 2 * k - 1;
        VAT(zc, k) = VAT(zf, kk);

        for (j=1; j<=*ny; j++) {
            jj = 2 * j - 1;
            VAT(yc, j) = VAT(yf, jj);

            for (i=1; i<=*nx; i++){
                   ii = 2 * i - 1;
                   VAT(xc, i) = VAT(xf, ii);

                   // True solution
                   VAT3( tc, i, j, k) = VAT3( tcf, ii, jj, kk);

                   // Helmholtz coefficient
                   VAT3( cc, i, j, k) = VAT3( ccf, ii, jj, kk);

                   // Source function
                   VAT3( fc, i, j, k) = VAT3( fcf, ii, jj, kk);

                   // East/West neighbor
                   VAT3(a1c, i, j, k) = VAT3(a1cf, ii, jj, kk);

                   // North/South neighbor
                   VAT3(a2c, i, j, k) = VAT3(a2cf, ii, jj, kk);

                   // Up/Down neighbor
                   VAT3(a3c, i, j, k) = VAT3(a3cf, ii, jj, kk);
            }
        }
    }

    // The (i=1) and (i=nx) boundaries
    for (k=1; k<=*nz; k++) {
        kk = 2 * k - 1;

        for (j=1; j<=*ny; j++) {
            jj = 2 * j - 1;

                VAT3(gxc, j, k, 1) = VAT3(gxcf, jj, kk, 1);
                VAT3(gxc, j, k, 2) = VAT3(gxcf, jj, kk, 2);
                VAT3(gxc, j, k, 3) = VAT3(gxcf, jj, kk, 3);
                VAT3(gxc, j, k, 4) = VAT3(gxcf, jj, kk, 4);
        }
    }

    // The (j=1) and (j=ny) boundaries
    for (k=1; k<=*nz; k++) {
        kk = 2 * k - 1;

        for (i=1; i<=*nx; i++) {
            ii = 2 * i - 1;

                VAT3(gyc, i, k, 1) = VAT3(gycf, ii, kk, 1);
                VAT3(gyc, i, k, 2) = VAT3(gycf, ii, kk, 2);
                VAT3(gyc, i, k, 3) = VAT3(gycf, ii, kk, 3);
                VAT3(gyc, i, k, 4) = VAT3(gycf, ii, kk, 4);
        }
    }

    // The (k=1) and (k=nz) boundaries
    for (j=1; j<=*ny; j++) {
        jj = 2 * j - 1;

        for (i=1; i<=*nx; i++) {
            ii = 2 * i - 1;

                VAT3(gzc, i, j, 1) = VAT3(gzcf, ii, jj, 1);
                VAT3(gzc, i, j, 2) = VAT3(gzcf, ii, jj, 2);
                VAT3(gzc, i, j, 3) = VAT3(gzcf, ii, jj, 3);
                VAT3(gzc, i, j, 4) = VAT3(gzcf, ii, jj, 4);
        }
    }
}

VPUBLIC void Vbuildharm0(int *nx, int *ny, int *nz,
        int *nxf, int *nyf, int *nzf,
        double *xc, double *yc, double *zc,
        double *gxc, double *gyc, double *gzc,
        double *a1c, double *a2c, double *a3c,
        double *cc, double *fc, double *tc,
        double *xf, double *yf, double *zf,
        double *gxcf, double *gycf, double *gzcf,
        double *a1cf, double *a2cf, double *a3cf,
        double *ccf, double *fcf, double *tcf) {
#if 1
    Vnm_print(2, "WARNING:  FUNCTION IS NOT FULLY IMPLEMENTED YET!!!");
#else
    int    i,    j,    k;
    int   ii,   jj,   kk;
    int iadd, jadd, kadd;

    MAT3( gxc, *ny, *nz, 4);
    MAT3( gyc, *nx, *nz, 4);
    MAT3( gzc, *nx, *ny, 4);

    MAT3( a1c, *nx, *ny, *nz);
    MAT3( a2c, *nx, *ny, *nz);
    MAT3( a3c, *nx, *ny, *nz);

    MAT3(  cc, *nx, *ny, *nz);
    MAT3(  fc, *nx, *ny, *nz);
    MAT3(  tc, *nx, *ny, *nz);

    MAT3(gxcf, *nyf, *nzf, 4);
    MAT3(gycf, *nxf, *nzf, 4);
    MAT3(gzcf, *nxf, *nyf, 4);
    MAT3(a1cf, *nxf, *nyf, *nzf);
    MAT3(a2cf, *nxf, *nyf, *nzf);
    MAT3(a3cf, *nxf, *nyf, *nzf);
    MAT3( tcf, *nxf, *nyf, *nzf);
    MAT3( ccf, *nxf, *nyf, *nzf);
    MAT3( fcf, *nxf, *nyf, *nzf);

    // Statement functions
    /// @todo  Figure out where the harmo ad arith functions come from
    double a, b, c, d, e, f, g, h;

    // How far to step into the coefficient arrays
    iadd = (*nxf - 1) / (*nx - 1);
    jadd = (*nyf - 1) / (*ny - 1);
    kadd = (*nzf - 1) / (*nz - 1);
    if (iadd !=2 || jadd != 2 || kadd != 2) {
       Vnm_print(2, "BUILDHARM0: problem with grid dimensions...\n");
    }

    // Compute the coarse grid pde coefficients
    for (k=1; k<=*nz; k++) {
        kk = 2 * k - 1;
        VAT(zc, k) = VAT(zf, kk);

        for (j=1; j<=*ny; j++) {
            jj = 2 * j - 1;
            VAT(yc, j) = VAT(yf,jj);

            for (i=1; i<=*nx; i++) {
                ii = 2 * i - 1;
                VAT(xc, i) = VAT(xf, ii);

                // True solution
                VAT3(tc, i, j, k) = VAT3(tcf, ii, jj, kk);

                // Helmholtz coefficient
                VAT3(cc, i, j, k) = VAT3(ccf, ii, jj, kk);

                /*  Commented out in original fortran code
                cc(i,j,k) = (
                    +0.5e0 * ccf(ii,jj,kk)
                    +0.5e0 * ARITH6( ccf(max0(1,ii-1),jj,kk),
                        ccf(min0(nxf,ii+1),jj,kk),
                        ccf(ii,max0(1,jj-1),kk),
                        ccf(ii,min0(nyf,jj+1),kk),
                        ccf(ii,jj,max0(nzf,kk-1)),
                        ccf(ii,jj,min0(nzf,kk+1)) )
                )
                */

                // Source function
                VAT3(fc, i, j, k) = VAT3(fcf, ii, jj, kk);
                /*
                fc(i,j,k) = (
                   +0.5e0 * fcf(ii,jj,kk)
                   +0.5e0 * ARITH6( fcf(max0(1,ii-1),jj,kk),
                                    fcf(min0(nxf,ii+1),jj,kk),
                                    fcf(ii,max0(1,jj-1),kk),
                                    fcf(ii,min0(nyf,jj+1),kk),
                                    fcf(ii,jj,max0(nzf,kk-1)),
                                    fcf(ii,jj,min0(nzf,kk+1)) )
                 )
                 */

                // East/West neighbor
                VAT3(a1c, i, j, k) = (
                    +0.500 * HARMO2(VAT3(a1cf, ii, jj, kk),
                                    VAT3(a1cf, VMIN2(*nxf, ii+1), jj, kk))
                    +0.125 * HARMO2(VAT3(a1cf, ii, jj, VMAX2(1, kk-1)),
                                    VAT3(a1cf, VMIN2(*nxf, ii+1), jj, VMAX2(1, kk-1)))
                    +0.125 * HARMO2(VAT3(a1cf, ii, jj, VMIN2(*nzf, kk+1)),
                                    VAT3(a1cf, VMIN2(*nxf, ii+1), jj, VMIN2(*nzf, kk+1)))
                    +0.125 * HARMO2(VAT3(a1cf, ii, VMAX2(1, jj-1), kk),
                                    VAT3(a1cf, VMIN2(*nxf, ii+1), VMAX2(1, jj-1), kk))
                    +0.125 * HARMO2(VAT3(a1cf, ii, VMIN2(*nyf, jj+1), kk),
                                    VAT3(a1cf, VMIN2(*nxf, ii+1), VMIN2(*nyf, jj+1), kk))
               );

                // North/South neighbor
                VAT3(a2c, i, j, k) = (
                    +0.500 * HARMO2(VAT3(a2cf, ii, jj, kk),
                                    VAT3(a2cf, ii, VMIN2(*nyf, jj+1), kk))
                    +0.125 * HARMO2(VAT3(a2cf, ii, jj, VMAX2(1, kk-1)),
                                    VAT3(a2cf, ii, VMIN2(*nyf, jj+1), VMAX2(1, kk-1)))
                    +0.125 * HARMO2(VAT3(a2cf, ii, jj, VMIN2(*nzf, kk+1)),
                                    VAT3(a2cf, ii, VMIN2(*nyf, jj+1), VMIN2(*nzf, kk+1)))
                    +0.125 * HARMO2(VAT3(a2cf, VMAX2(1, ii-1), jj, kk),
                                    VAT3(a2cf, VMAX2(1, ii-1), VMIN2(nyf, jj+1), kk))
                    +0.125 * HARMO2(VAT3(a2cf, VMIN2(*nxf, ii+1), jj, kk),
                                    VAT3(a2cf, VMIN2(*nxf, ii+1), VMIN2(*nyf, jj+1), kk))
               );

                // Up/Down neighbor
                VAT3(a3c, i, j, k) = (
                    +0.500 * HARMO2(VAT3(a3cf, ii, jj, kk),
                                    VAT3(a3cf, ii, jj, VMIN2(*nzf, kk+1)))
                    +0.125 * HARMO2(VAT3(a3cf, ii, VMAX2(1, jj-1), kk),
                                    VAT3(a3cf, ii, VMAX2(1, jj-1), VMIN2(*nzf, kk+1)))
                    +0.125 * HARMO2(VAT3(a3cf, ii, VMIN2(*nyf, jj+1), kk),
                                    VAT3(a3cf, ii, VMIN2(*nyf, jj+1), VMIN2(*nzf, kk+1)))
                    +0.125 * HARMO2(VAT3(a3cf, VMAX2(1, ii-1), jj, kk),
                                    VAT3(a3cf, VMAX2(1, ii-1), jj, VMIN2(*nzf, kk+1)))
                    +0.125 * HARMO2(VAT3(a3cf, VMIN2(*nxf, ii+1), jj, kk),
                                    VAT3(a3cf, VMIN2(*nxf, ii+1), jj, VMIN2(*nzf, kk+1)))
               );
            }
        }
    }

    // The (i=1) and (i=nx) boundaries ***
    for (k=1; k<=*nz; k++) {
       kk = 2 * k - 1;

       for (j=1; j<=*ny; j++) {
          jj = 2 * j - 1;

          VAT3(gxc, j, k, 1) = VAT3(gxcf, jj, kk, 1);
          VAT3(gxc, j, k, 2) = VAT3(gxcf, jj, kk, 2);
          VAT3(gxc, j, k, 3) = VAT3(gxcf, jj, kk, 3);
          VAT3(gxc, j, k, 4) = VAT3(gxcf, jj, kk, 4);
       }
    }

    // The (j=1) and (j=ny) boundaries
    for (k=1; k<=*nz; k++) {
       kk = 2 * k - 1;

       for (i=1; i<=*nx; i++) {
          ii = 2 * i - 1;
          VAT3(gyc, i, k, 1) = VAT3(gycf, ii, kk, 1);
          VAT3(gyc, i, k, 2) = VAT3(gycf, ii, kk, 2);
          VAT3(gyc, i, k, 3) = VAT3(gycf, ii, kk, 3);
          VAT3(gyc, i, k, 4) = VAT3(gycf, ii, kk, 4);
       }
    }

    // The (k=1) and (k=nz) boundaries
    for (j=1; j<=*ny; j++) {
       jj = 2 * j - 1;

       for (i=1; i<=*nx; i++) {
          ii = 2 * i - 1;

          VAT3(gzc, i, j, 1) = VAT3(gzcf, ii, jj, 1);
          VAT3(gzc, i, j, 2) = VAT3(gzcf, ii, jj, 2);
          VAT3(gzc, i, j, 3) = VAT3(gzcf, ii, jj, 3);
          VAT3(gzc, i, j, 4) = VAT3(gzcf, ii, jj, 4);
       }
    }
#endif
}



VPUBLIC void Vbuildalg(int *nx, int *ny, int *nz,
        int *mode, int *nlev, int *iz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *y, double *tmp) {

    int   nxx,   nyy,   nzz;
    int nxold, nyold, nzold;
    int lev, numlev;

    MAT2(iz, 50, *nlev);

    // Setup
    nxx = *nx;
    nyy = *ny;
    nzz = *nz;

    // Build the rhs the finest level
    lev = 1;
    if ((*mode == 1) || (*mode == 2)) {
        Vnmatvec(&nxx, &nyy, &nzz,
                 RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                 RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                 RAT(  x, VAT2(iz, 1, lev)), RAT(  y, VAT2(iz, 1, lev)),
                 tmp);
    } else {
        Vmatvec(&nxx, &nyy, &nzz,
                RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                RAT(  x, VAT2(iz, 1, lev)), RAT(  y, VAT2(iz, 1, lev)));
    }

    // Build the (nlev-1) level rhs function
    for (lev=2; lev <= *nlev; lev++) {
        nxold = nxx;
        nyold = nyy;
        nzold = nzz;

        numlev = 1;
        Vmkcors(&numlev, &nxold, &nyold, &nzold, &nxx, &nyy, &nzz);

        // Build the rhs on this level
        if ((*mode == 1) || (*mode == 2)) {
            Vnmatvec(&nxx, &nyy, &nzz,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                    RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                    RAT(  x, VAT2(iz, 1, lev)), RAT(  y, VAT2(iz, 1, lev)),
                    tmp);
        } else {
            Vmatvec(&nxx, &nyy, &nzz,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                    RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                    RAT(  x, VAT2(iz, 1, lev)), RAT(  y, VAT2(iz, 1, lev)));
        }
    }
}
