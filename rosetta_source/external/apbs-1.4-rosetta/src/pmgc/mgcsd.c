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

#include "mgcsd.h"

VEXTERNC void Vmvcs(int *nx, int *ny, int *nz,
        double *x,
        int *iz,
        double *w0, double *w1, double *w2, double *w3,
        int *istop, int *itmax, int *iters, int *ierror,
        int *nlev, int *ilev, int *nlev_real,
        int *mgsolv, int *iok, int *iinfo,
        double *epsiln, double *errtol, double *omega,
        int *nu1, int *nu2,
        int *mgsmoo,
        int *ipc, double *rpc,
        double *pc, double *ac, double *cc, double *fc, double *tru) {

    int level;       // @todo: doc
    int lev;         // @todo: doc
    int itmax_s;     // @todo: doc
    int iters_s;     // @todo: doc
    int nuuu;        // @todo: doc
    int mgsmoo_s;    // @todo: doc
    int iresid;      // @todo: doc
    int nxf;         // @todo: doc
    int nyf;         // @todo: doc
    int nzf;         // @todo: doc
    int nxc;         // @todo: doc
    int nyc;         // @todo: doc
    int nzc;         // @todo: doc
    int lpv;         // @todo: doc
    int n;           // @todo: doc
    int m;           // @todo: doc
    int iadjoint;    // @todo: doc
    double errtol_s; // @todo: doc
    double rsden;    // @todo: doc
    double rsnrm;    // @todo: doc
    double orsnrm;   // @todo: doc
    double xnum;     // @todo: doc
    double xden;     // @todo: doc
    double xdamp;    // @todo: doc
    int lda;         // @todo: doc

    double alpha;     // A utility variable used to pass a parameter to xaxpy
    int numlev;       // A utility variable used to pass a parameter to mkcors

    MAT2(iz, 50, 1);

    // Recover level information
    level = 1;
    lev = (*ilev - 1) + level;

    // Recover grid sizes
    nxf = *nx;
    nyf = *ny;
    nzf = *nz;
    numlev = *nlev - 1;
    Vmkcors(&numlev, &nxf, &nyf, &nzf, &nxc, &nyc, &nzc);

    // Do some i/o if requested
    if (*iinfo > 1) {
        VMESSAGE0("Starting mvcs operation");
        VMESSAGE3("Fine Grid Size:   (%d, %d, %d)", nxf, nyf, nzf);
        VMESSAGE3("Coarse Grid Size: (%d, %d, %d)", nxc, nyc, nzc);
    }

    if (*iok != 0) {
        Vprtstp(*iok, -1, 0.0, 0.0, 0.0);

    }

    /*    **************************************************************
     *    *** Note: if (iok != 0) then:  use a stopping test.        ***
     *    ***       else:  use just the itmax to stop iteration.     ***
     *    **************************************************************
     *    *** istop=0 most efficient (whatever it is)                ***
     *    *** istop=1 relative residual                              ***
     *    *** istop=2 rms difference of successive iterates          ***
     *    *** istop=3 relative true error (provided for testing)     ***
     *    **************************************************************/

    // Compute denominator for stopping criterion
    if (*iok != 0) {
        if (*istop == 0) {
            rsden = 1.0;
        }
        else if (*istop == 1) {
            rsden = Vxnrm1(&nxf, &nyf, &nzf, RAT(fc, VAT2(iz, 1,lev)));
        }
        else if (*istop == 2) {
            rsden = VSQRT(nxf * nyf * nzf);
        }
        else if (*istop == 3) {
            rsden = Vxnrm2(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)));
        }
        else if (*istop == 4) {
            rsden = Vxnrm2(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)));
        }
        else if (*istop == 5) {
            Vmatvec(&nxf, &nyf, &nzf,
                RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                 RAT(ac, VAT2(iz, 7,lev)),  RAT(cc, VAT2(iz, 1,lev)),
                RAT(tru, VAT2(iz, 1,lev)),  w1);
            rsden = VSQRT(Vxdot(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)), w1));
        }
        else {
            VABORT_MSG1("Bad istop value: %d", *istop);
        }

        if (rsden == 0.0) {
            rsden = 1.0;
            VERRMSG0("rhs is zero on finest level");
        }
        rsnrm = rsden;
        orsnrm = rsnrm;
        iters_s = 0;

        Vprtstp(*iok, 0, rsnrm, rsden, orsnrm);
    }



    /* *********************************************************************
     * *** solve directly if nlev = 1
     * *********************************************************************/

    // Solve directly if on the coarse grid
    if (*nlev == 1) {

        // Use iterative method?
        if (*mgsolv == 0) {

            // solve on coarsest grid with cghs, mgsmoo_s=4 (no residual)
            iresid = 0;
            iadjoint = 0;
            itmax_s  = 100;
            iters_s  = 0;
            errtol_s = *epsiln;
            mgsmoo_s = 4;

            Vazeros(&nxf, &nyf, &nzf, RAT(x, VAT2(iz, 1,lev)));

            Vsmooth(&nxf, &nyf, &nzf,
                    RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                     RAT(ac, VAT2(iz, 7,lev)), RAT(cc, VAT2(iz, 1,lev)), RAT(fc, VAT2(iz, 1,lev)),
                      RAT(x, VAT2(iz, 1,lev)), w1, w2, w3,
                      &itmax_s, &iters_s, &errtol_s, omega,
                      &iresid, &iadjoint, &mgsmoo_s);

            // Check for trouble on the coarse grid
            VWARN_MSG2(iters_s <= itmax_s,
                "Exceeded maximum iterations: iters_s=%d, itmax_s=%d",
                iters_s, itmax_s);

        } else if (*mgsolv == 1) {

            // Use direct method?

            // Setup lpv to access the factored/banded operator
            lpv = lev + 1;

            // setup for banded format
            n   = *RAT(ipc, (VAT2(iz, 5,lpv) - 1) + 1);
            m   = *RAT(ipc, (VAT2(iz, 5,lpv) - 1) + 2);
            lda = *RAT(ipc, (VAT2(iz, 5,lpv) - 1) + 3);

            // Call dpbsl to solve
            Vxcopy_small(&nxf, &nyf, &nzf, RAT(fc, VAT2(iz, 1,lev)), w1);
            Vdpbsl(RAT(ac, VAT2(iz, 7,lpv)), &lda, &n, &m, w1);
            Vxcopy_large(&nxf, &nyf, &nzf, w1, RAT(x, VAT2(iz, 1,lev)));
            VfboundPMG00(&nxf, &nyf, &nzf, RAT(x, VAT2(iz, 1,lev)));

        } else {
            VABORT_MSG1("Invalid coarse solver requested: %d", *mgsolv);
        }


        // Compute the stopping test
        *iters = 1;
        if (*iok != 0) {

            orsnrm = rsnrm;

            if (*istop == 0) {

                Vmresid(&nxf, &nyf, &nzf,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                    RAT( ac, VAT2(iz, 7, lev)), RAT(cc , VAT2(iz, 1, lev)),
                    RAT( fc, VAT2(iz, 1,lev)),
                    RAT(  x, VAT2(iz, 1, lev)), w1);

                rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
            }

            else if (*istop == 1) {

                Vmresid(&nxf, &nyf, &nzf,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                    RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                    RAT( fc, VAT2(iz, 1, lev)), RAT(  x, VAT2(iz, 1, lev)),
                    w1);
                rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
            }

            else if (*istop == 2) {

                alpha = -1.0;

                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)), w1);
                Vxaxpy(&nxf, &nyf, &nzf, &alpha,
                        RAT(x, VAT2(iz, 1,lev)), w1);
                rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
                Vxcopy(&nxf, &nyf, &nzf, RAT(x, VAT2(iz, 1,lev)), RAT(tru, VAT2(iz, 1,lev)));
            }

            else if (*istop == 3) {

                alpha = -1.0;

                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)), w1);
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1,lev)), w1);
                rsnrm = Vxnrm2(&nxf, &nyf, &nzf, w1);
            }

            else if (*istop == 4) {

                alpha = -1.0;

                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)), w1);
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1,lev)), w1);
                rsnrm = Vxnrm2(&nxf, &nyf, &nzf, w1);
            }

            else if (*istop == 5) {

                alpha = -1.0;

                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)), w1);
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1,lev)), w1);

                Vmatvec(&nxf, &nyf, &nzf,
                        RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                         RAT(ac, VAT2(iz, 7,lev)),  RAT(cc, VAT2(iz, 1,lev)),
                         w1, w2);
                rsnrm = VSQRT(Vxdot(&nxf, &nyf, &nzf, w1, w2));
            }

            else {
                VABORT_MSG1("Bad istop value: %d\n", *istop);
            }
            Vprtstp(*iok, *iters, rsnrm, rsden, orsnrm);
        }
        return;
    }


    /* *********************************************************************
     * *** begin mg iteration (note nxf,nyf,nzf changes during loop)
     * *********************************************************************/

    // Setup for the v-cycle looping
    *iters = 0;
    do {

        // Finest level initialization
        level = 1;
        lev   = (*ilev - 1) + level;

        // nu1 pre-smoothings on fine grid (with residual)
        iresid = 1;
        iadjoint = 0;
        iters_s  = 0;
        errtol_s = 0.0;
        nuuu = Vivariv(nu1, &lev);

        Vsmooth(&nxf, &nyf, &nzf,
                RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                 RAT(ac, VAT2(iz, 7,lev)),  RAT(cc, VAT2(iz, 1,lev)),  RAT(fc, VAT2(iz, 1,lev)),
                  RAT(x, VAT2(iz, 1,lev)), w2, w3, w1,
                &nuuu, &iters_s,
                &errtol_s, omega,
                &iresid, &iadjoint, mgsmoo);

        Vxcopy(&nxf, &nyf, &nzf, w1, RAT(w0, VAT2(iz, 1,lev)));



        /* *********************************************************************
         * begin cycling down to coarse grid
         * *********************************************************************/

        // Go down grids: restrict resid to coarser and smooth
        for (level=2; level<=*nlev; level++) {

            lev = (*ilev - 1) + level;

            // Find new grid size
            numlev = 1;
            Vmkcors(&numlev, &nxf, &nyf, &nzf, &nxc, &nyc, &nzc);

            // Restrict residual to coarser grid ***
            Vrestrc(&nxf, &nyf, &nzf,
                    &nxc, &nyc, &nzc,
                    w1, RAT(w0, VAT2(iz, 1,lev)), RAT(pc, VAT2(iz, 11,lev-1)));

            /// New grid size
            nxf = nxc;
            nyf = nyc;
            nzf = nzc;

            // if not on coarsest level yet...
            if (level != *nlev) {

                // nu1 pre-smoothings on this level (with residual)
                // (w1 has residual...)
                Vazeros(&nxf, &nyf, &nzf, RAT(x, VAT2(iz, 1,lev)));
                iresid = 1;
                iadjoint = 0;
                iters_s  = 0;
                errtol_s = 0.0;
                nuuu = Vivariv(nu1, &lev);
                Vsmooth(&nxf, &nyf, &nzf,
                       RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                        RAT(ac, VAT2(iz, 7,lev)),  RAT(cc, VAT2(iz, 1,lev)), RAT(w0, VAT2(iz, 1,lev)),
                        RAT(x, VAT2(iz, 1,lev)), w2, w3, w1,
                        &nuuu, &iters_s,
                        &errtol_s, omega ,
                        &iresid, &iadjoint, mgsmoo);
            }
            // End of cycling down to coarse grid loop
        }



        /* *********************************************************************
         * begin coarse grid
         * *********************************************************************/

        // Coarsest level
        level = *nlev;
        lev = (*ilev - 1) + level;

        // Use iterative method?
        if (*mgsolv == 0) {

            // solve on coarsest grid with cghs, mgsmoo_s=4 (no residual)
            iresid = 0;
            iadjoint = 0;
            itmax_s  = 100;
            iters_s  = 0;
            errtol_s = *epsiln;
            mgsmoo_s = 4;
            Vazeros(&nxf, &nyf, &nzf, RAT(x, VAT2(iz, 1,lev)));
            Vsmooth(&nxf, &nyf, &nzf,
                    RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                     RAT(ac, VAT2(iz, 7,lev)), RAT(cc, VAT2(iz, 1,lev)), RAT(w0, VAT2(iz, 1,lev)),
                       RAT(x, VAT2(iz, 1,lev)), w1, w2, w3,
                    &itmax_s, &iters_s,
                    &errtol_s, omega,
                    &iresid, &iadjoint, &mgsmoo_s);

            // Check for trouble on the coarse grid
            VWARN_MSG2(iters_s <= itmax_s,
                "Exceeded maximum iterations: iters_s=%d, itmax_s=%d",
                iters_s, itmax_s);
        } else if (*mgsolv == 1) {

            // use direct method?

            // Setup lpv to access the factored/banded operator
            lpv = lev + 1;

            // Setup for banded format
            n   = VAT(ipc, (VAT2(iz, 5, lpv) - 1) + 1);
            m   = VAT(ipc, (VAT2(iz, 5, lpv) - 1) + 2);
            lda = VAT(ipc, (VAT2(iz, 5, lpv) - 1) + 3);

            // Call dpbsl to solve
            Vxcopy_small(&nxf, &nyf, &nzf, RAT(w0, VAT2(iz, 1,lev)), w1);
            Vdpbsl(RAT(ac, VAT2(iz, 7,lpv)), &lda, &n, &m, w1);
            Vxcopy_large(&nxf, &nyf, &nzf, w1, RAT(x, VAT2(iz, 1,lev)));
            VfboundPMG00(&nxf, &nyf, &nzf, RAT(x, VAT2(iz, 1,lev)));

        } else {
            VABORT_MSG1("Invalid coarse solver requested: %d", *mgsolv);
        }


        /* *********************************************************************
         * begin cycling back to fine grid
         * *********************************************************************/


        // Move up grids: interpolate resid to finer and smooth
        for (level=*nlev-1; level>=1; level--) {

            lev = (*ilev - 1) + level;

            // Find new grid size
            numlev = 1;
            Vmkfine(&numlev,
                    &nxf, &nyf, &nzf,
                    &nxc, &nyc, &nzc);

            // Interpolate to next finer grid
            VinterpPMG(&nxf, &nyf, &nzf,
                    &nxc, &nyc, &nzc,
                    RAT(x, VAT2(iz, 1,lev+1)), w1, RAT(pc, VAT2(iz, 11,lev)));

            /* Compute the hackbusch/reusken damping parameter
             * which is equivalent to the standard linear cg steplength
             */
            Vmatvec(&nxf, &nyf, &nzf,
                    RAT(ipc, VAT2(iz, 5,lev+1)), RAT(rpc, VAT2(iz, 6,lev+1)),
                     RAT(ac, VAT2(iz, 7,lev+1)),  RAT(cc, VAT2(iz, 1,lev+1)),
                      RAT(x, VAT2(iz, 1,lev+1)),  w2);

            xnum = Vxdot(&nxf, &nyf, &nzf,
                    RAT(x, VAT2(iz, 1,lev+1)), RAT(w0, VAT2(iz, 1,lev+1)));

            xden = Vxdot(&nxf, &nyf, &nzf,
                    RAT(x, VAT2(iz, 1,lev+1)), w2);
            xdamp = xnum / xden;

            // New grid size
            nxf = nxc;
            nyf = nyc;
            nzf = nzc;

            // perform the coarse grid correction
            // xdamp = 1.0d0
            Vxaxpy(&nxf, &nyf, &nzf,
                    &xdamp, w1, RAT(x, VAT2(iz, 1,lev)));

            // nu2 post-smoothings for correction (no residual)
            iresid = 0;
            iadjoint = 1;
            iters_s  = 0;
            errtol_s = 0.0;
            nuuu = Vivariv(nu2, &lev);
            if (level == 1) {
                Vsmooth(&nxf, &nyf, &nzf,
                        RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                         RAT(ac, VAT2(iz, 7,lev)), RAT(cc, VAT2(iz, 1,lev)), RAT(fc, VAT2(iz, 1,lev)),
                          RAT(x, VAT2(iz, 1,lev)),w1,w2,w3,
                         &nuuu, &iters_s, &errtol_s, omega,
                         &iresid, &iadjoint, mgsmoo);
            } else {
                Vsmooth(&nxf, &nyf, &nzf,
                        RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                         RAT(ac, VAT2(iz, 7,lev)),  RAT(cc, VAT2(iz, 1,lev)), RAT(w0, VAT2(iz, 1,lev)),
                          RAT(x, VAT2(iz, 1,lev)), w1, w2, w3,
                         &nuuu, &iters_s, &errtol_s, omega,
                         &iresid, &iadjoint, mgsmoo);
            }
        }

        /* *********************************************************************
         * iteration complete: do some i/o
         * *********************************************************************/

        // Increment the iteration counter
        (*iters)++;

        // Compute/check the current stopping test
        if (iok != 0) {
            orsnrm = rsnrm;
            if (*istop == 0) {
                Vmresid(&nxf, &nyf, &nzf,
                        RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                         RAT(ac, VAT2(iz, 7,lev)), RAT(cc, VAT2(iz, 1,lev)), RAT(fc, VAT2(iz, 1,lev)),
                          RAT(x, VAT2(iz, 1,lev)), w1);
                rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
            } else if(*istop == 1) {
                Vmresid(&nxf, &nyf, &nzf,
                        RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                         RAT(ac, VAT2(iz, 7,lev)),  RAT(cc, VAT2(iz, 1,lev)), RAT(fc, VAT2(iz, 1,lev)),
                          RAT(x, VAT2(iz, 1,lev)), w1);
                rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
            } else if (*istop == 2) {
                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)), w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1,lev)), w1);
                rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
                Vxcopy(&nxf, &nyf, &nzf, RAT(x, VAT2(iz, 1,lev)), RAT(tru, VAT2(iz, 1,lev)));
            } else if (*istop == 3) {
                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)), w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1,lev)), w1);
                rsnrm = Vxnrm2(&nxf, &nyf, &nzf, w1);
            } else if (*istop == 4) {
                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)), w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1,lev)), w1);
                rsnrm = Vxnrm2(&nxf, &nyf, &nzf, w1);
            } else if (*istop == 5) {
                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1,lev)), w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1,lev)), w1);
                Vmatvec(&nxf, &nyf, &nzf,
                        RAT(ipc, VAT2(iz, 5,lev)), RAT(rpc, VAT2(iz, 6,lev)),
                         RAT(ac, VAT2(iz, 7,lev)),  RAT(cc, VAT2(iz, 1,lev)),
                         w1, w2);
                rsnrm = VSQRT(Vxdot(&nxf, &nyf, &nzf, w1, w2));
            } else {
                VABORT_MSG1("Bad istop value: %d", *istop);
            }
            Vprtstp(*iok, *iters, rsnrm, rsden, orsnrm);
        }
    } while (*iters<*itmax && (rsnrm/rsden) > *errtol);

    *ierror = *iters < *itmax ? 0 : 1;
}
