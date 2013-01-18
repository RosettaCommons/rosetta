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

#include "mgfasd.h"

VPUBLIC void Vfmvfas(int *nx, int *ny, int *nz,
        double *x,
        int *iz,
        double *w0, double *w1, double *w2, double *w3, double *w4,
        int *istop, int *itmax, int *iters, int *ierror,
        int *nlev, int *ilev, int *nlev_real,
        int *mgsolv, int *iok, int *iinfo,
        double *epsiln, double *errtol, double *omega,
        int *nu1, int *nu2, int *mgsmoo,
        int *ipc, double *rpc,
        double *pc, double *ac, double *cc, double *fc, double *tru) {

    // Other Declarations
    int level, itmxd, nlevd, iterd, iokd;
    int nxf, nyf, nzf;
    int nxc, nyc, nzc;
    int istpd, iinfod;
    double errd;

    // Utility variables
    int numlev;

    MAT2(iz, 50, *nlev);

    // Recover gridsizes
    nxf = *nx;
    nyf = *ny;
    nzf = *nz;

    numlev = *nlev - 1;
    Vmkcors(&numlev, &nxf, &nyf, &nzf, &nxc, &nyc, &nzc);

    // Move up grids: interpolate solution to finer, do v cycle
    if (*iinfo != 0) {

        Vnm_print(2, "Vfmvfas: starting: (%d,%d,%d) (%d,%d,%d)\n",
                nxf, nyf, nzf, nxc, nyc, nzc);
    }

    for (level = *nlev_real; level >= *ilev + 1; level--) {

        // Call mv cycle
        errd   = 1.0e-5;
        itmxd  = 1;
        nlevd  = *nlev_real - level + 1;
        iterd  = 0;
        iokd   = 2;
        iinfod = *iinfo;
        istpd  = *istop;
        if (*iinfo >= 2)
            iokd = 2;

        Vmvfas(&nxc, &nyc, &nzc,
                x,
                iz,
                w0, w1, w2, w3, w4,
                &istpd, &itmxd, &iterd, ierror,
                &nlevd, &level, nlev_real,
                mgsolv, &iokd, &iinfod,
                epsiln, errtol, omega,
                nu1, nu2, mgsmoo,
                ipc, rpc,
                pc, ac, cc, fc, tru);

        // Find new grid size
        numlev = 1;
        Vmkfine(&numlev, &nxc, &nyc, &nzc, &nxf, &nyf, &nzf);

        // Interpolate to next finer grid
        VinterpPMG(&nxc, &nyc, &nzc,
                &nxf, &nyf, &nzf,
                RAT( x, VAT2(iz,  1,   level)),
                RAT( x, VAT2(iz,  1, level-1)),
                RAT(pc, VAT2(iz, 11, level-1)));

        // New grid size
        nxc = nxf;
        nyc = nyf;
        nzc = nzf;
    }



    // Call mv cycle
    level = *ilev;

    Vmvfas(&nxf, &nyf, &nzf,
            x, iz,
            w0, w1, w2, w3, w4,
            istop, itmax, iters,
            ierror, nlev, &level, nlev_real,
            mgsolv, iok, iinfo,
            epsiln, errtol, omega,
            nu1, nu2, mgsmoo,
            ipc, rpc,
            pc, ac, cc, fc, tru);
}



VPUBLIC void Vmvfas(int *nx, int *ny, int *nz,
        double *x,
        int *iz,
        double *w0, double *w1, double *w2, double *w3, double *w4,
        int *istop, int *itmax, int *iters, int *ierror,
        int *nlev, int *ilev, int *nlev_real,
        int *mgsolv, int *iok, int *iinfo,
        double *epsiln, double *errtol, double *omega,
        int *nu1, int *nu2, int *mgsmoo,
        int *ipc, double *rpc,
        double *pc, double *ac, double *cc, double *fc, double *tru) {

    // Other declarations
    int level, lev;
    int itmax_s, iters_s, nuuu, ivariv, mgsmoo_s, iresid;
    int nxf, nyf, nzf;
    int nxc, nyc, nzc;
    int iadjoint;
    double errtol_s;
    double rsden, rsnrm, orsnrm;
    double xdamp;

    int numlev;
    double alpha;

    MAT2(iz, 50, *nlev);

    WARN_UNTESTED;

    // Recover level information
    level = 1;
    lev   = (*ilev - 1) + level;

    // Recover gridsizes
    nxf = *nx;
    nyf = *ny;
    nzf = *nz;

    numlev = *nlev - 1;
    Vmkcors(&numlev, &nxf, &nyf, &nzf, &nxc, &nyc, &nzc);

    // Do some i/o if requested
    if (*iinfo != 0) {
        Vnm_print(2, "Vmvfas: starting:  (%d, %d, %d) (%d, %d, %d)\n",
                nxf, nyf, nzf, nxc, nyc, nzc);
    }

    // Initial wall clock
    if (*iok != 0) {
        Vprtstp(*iok, -1, 0.0, 0.0, 0.0);
    }

    /**************************************************************
     *** note: if (iok != 0) then:  use a stopping test.        ***
     ***       else:  use just the itmax to stop iteration.     ***
     **************************************************************
     *** istop=0 most efficient (whatever it is)                ***
     *** istop=1 relative residual                              ***
     *** istop=2 rms difference of successive iterates          ***
     *** istop=3 relative true error (provided for testing)     ***
     **************************************************************/

    // Compute denominator for stopping criterion
    if (*iok != 0) {
        if (*istop == 0) {
            rsden = 1.0;
        } else if (*istop == 1) {

            // Compute initial residual with zero initial guess
            // this is analogous to the linear case where one can
            // simply take norm of rhs for a zero initial guess
            Vazeros(&nxf, &nyf, &nzf, w1);

            Vnmresid(&nxf, &nyf, &nzf,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                    RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                    RAT( fc, VAT2(iz, 1, lev)),
                    w1, w2, w3);

            rsden = Vxnrm1(&nxf, &nyf, &nzf, w2);

         } else if (*istop == 2) {
            rsden = VSQRT(nxf * nyf * nzf);
         } else if (*istop == 3) {
            rsden = Vxnrm2(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1, lev)));
         } else if (*istop == 4) {
            rsden = Vxnrm2(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1, lev)));
         } else if (*istop == 5) {
            Vnmatvec(&nxf, &nyf, &nzf,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                    RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                    RAT(tru, VAT2(iz, 1, lev)),
                    w1, w2);
            rsden = VSQRT(Vxdot(&nxf, &nyf, &nzf,RAT(tru, VAT2(iz, 1, lev)), w1));
         } else {
            Vnm_print(2, "Vmvfas: bad istop value: %d\n", *istop);
         }
         if (rsden == 0.0) {
            rsden = 1.0;
            Vnm_print(2, "Vmfas: rhs is zero on finest level\n");
         }
         rsnrm = rsden;
         orsnrm = rsnrm;

         Vprtstp(*iok, 0, rsnrm, rsden, orsnrm);
    }



    /********************************************************************
     *** solve directly if nlev = 1
     ********************************************************************/

    // Solve directly if on the coarse grid
    if (*nlev == 1) {

        //solve with ncghs, mgsmoo_s=4 (no residual)
        iresid = 0;
        iadjoint = 0;
        itmax_s  = 100;
        iters_s  = 0;
        errtol_s = *epsiln;
        mgsmoo_s = *mgsmoo;
        Vazeros(&nxf, &nyf, &nzf,RAT(x, VAT2(iz, 1, lev)));
        Vnsmooth(&nxf, &nyf, &nzf,
                RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                RAT( fc, VAT2(iz, 1, lev)),
                RAT(  x, VAT2(iz, 1, lev)),
                w1, w2, w3,
                &itmax_s, &iters_s, &errtol_s, omega,
                &iresid, &iadjoint, &mgsmoo_s);

        // Compute the stopping test
        *iters = 1;
        if (*iok != 0) {
            orsnrm = rsnrm;
            if (*istop == 0) {
                    Vnmresid(&nxf, &nyf, &nzf,
                            RAT(ipc, VAT2(iz, 5, lev)),
                            RAT(rpc, VAT2(iz, 6, lev)),
                            RAT( ac, VAT2(iz, 7, lev)),
                            RAT( cc, VAT2(iz, 1, lev)),
                            RAT( fc, VAT2(iz, 1, lev)),
                            RAT(  x, VAT2(iz, 1, lev)),
                            w1, w2);
                    rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
            } else if (*istop == 1) {
                    Vnmresid(&nxf, &nyf, &nzf,
                            RAT(ipc, VAT2(iz, 5, lev)),
                            RAT(rpc, VAT2(iz, 6, lev)),
                            RAT( ac, VAT2(iz, 7, lev)),
                            RAT( cc, VAT2(iz, 1, lev)),
                            RAT( fc, VAT2(iz, 1, lev)),
                            RAT(  x, VAT2(iz, 1, lev)),
                            w1, w2);
                    rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
            } else if (*istop == 2) {
                Vxcopy(&nxf, &nyf, &nzf,
                        RAT(tru, VAT2(iz, 1, lev)),
                        w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf,
                        &alpha, RAT(x, VAT2(iz, 1, lev)), w1);
                rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
                Vxcopy(&nxf, &nyf, &nzf,
                        RAT(x, VAT2(iz, 1, lev)),
                        RAT(tru, VAT2(iz, 1, lev)));
            } else if (*istop == 3) {
                Vxcopy(&nxf, &nyf, &nzf,
                        RAT(tru, VAT2(iz, 1, lev)),
                        w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf,
                        &alpha, RAT(x, VAT2(iz, 1, lev)), w1);
                rsnrm = Vxnrm2(&nxf, &nyf, &nzf, w1);
            } else if (*istop == 4) {
                Vxcopy(&nxf, &nyf, &nzf,
                            RAT(tru, VAT2(iz, 1, lev)),
                            w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf,
                            &alpha, RAT(x, VAT2(iz, 1, lev)), w1);
                rsnrm = Vxnrm2(&nxf, &nyf, &nzf, w1);
            } else if (*istop == 5) {
                Vxcopy(&nxf, &nyf, &nzf,
                        RAT(tru, VAT2(iz, 1, lev)),
                        w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf,
                        &alpha, RAT(x, VAT2(iz, 1, lev)), w1);
                Vnmatvec(&nxf, &nyf, &nzf,
                        RAT(ipc, VAT2(iz, 5, lev)),
                        RAT(rpc, VAT2(iz, 6, lev)),
                        RAT( ac, VAT2(iz, 7, lev)),
                        RAT( cc, VAT2(iz, 1, lev)),
                        w1, w2, w3);
                rsnrm = VSQRT(Vxdot(&nxf, &nyf, &nzf, w1, w2));
            } else {
                Vnm_print(2, "Vmvcs: bad istop value: %d\n", *istop);
            }

            Vprtstp(*iok, *iters, rsnrm, rsden, orsnrm);
        }
        return;
    }

    /*********************************************************************
     *** begin mg iteration (note nxf,nyf,nzf changes during loop)
     *********************************************************************/

    // setup for the v-cycle looping ***
    *iters = 0;
    while(1) {

        // Finest level initialization
        level = 1;
        lev   = (*ilev - 1) + level;

        // Nu1 pre-smoothings on fine grid (with residual)
        iresid = 1;
        iadjoint = 0;
        iters_s  = 0;
        errtol_s = 0.0;
        nuuu = Vivariv(nu1, &lev);
        Vnsmooth(&nxf, &nyf, &nzf,
                RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                RAT( fc, VAT2(iz, 1, lev)), RAT(  x, VAT2(iz, 1, lev)),
                w2, w3, w1,
                &nuuu, &iters_s, &errtol_s, omega,
                &iresid, &iadjoint, mgsmoo);
        Vxcopy(&nxf, &nyf, &nzf, w1, RAT(w0, VAT2(iz, 1, lev)));

        /* *********************************************************************
         * begin cycling down to coarse grid
         * ********************************************************************/

        // Go down grids: restrict resid to coarser and smooth
        for (level = 2; level <= *nlev; level++) {
            lev = (*ilev - 1) + level;

            // Find new grid size
            numlev = 1;
            Vmkcors(&numlev, &nxf, &nyf, &nzf, &nxc, &nyc, &nzc);

            // Restrict residual to coarser grid
            Vrestrc(&nxf, &nyf, &nzf, &nxc, &nyc, &nzc,
                    w1, RAT(w0, VAT2(iz, 1, lev)),
                    RAT(pc, VAT2(iz, 11, lev-1)));

            // Restrict (extract) solution to coarser grid
            Vextrac(&nxf, &nyf, &nzf, &nxc, &nyc, &nzc,
                    RAT(x, VAT2(iz, 1, lev-1)), RAT(w4, VAT2(iz, 1, lev)));

            // New grid size
            nxf = nxc;
            nyf = nyc;
            nzf = nzc;

            // Apply coarse grid operator to coarse grid soln
            Vnmatvec(&nxf, &nyf, &nzf,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                    RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                    RAT( w4, VAT2(iz, 1, lev)), RAT( fc, VAT2(iz, 1, lev)),
                    w3);

            // Build coarse grid right hand side
            alpha = 1.0;
            Vxaxpy(&nxf, &nyf, &nzf, &alpha,
                    RAT(w0, VAT2(iz, 1, lev)), RAT(fc, VAT2(iz, 1, lev)));

            // If not on coarsest level yet...
            if (level != *nlev) {

                // nu1 pre-smoothings on this level (with residual)
               Vxcopy(&nxf, &nyf, &nzf,
                       RAT(w4, VAT2(iz, 1, lev)), RAT(x, VAT2(iz, 1, lev)));
               iresid = 1;
               iadjoint = 0;
               iters_s  = 0;
               errtol_s = 0.0;
               nuuu = Vivariv (nu1,&lev);
               Vnsmooth(&nxf, &nyf, &nzf,
                       RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                       RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                       RAT( fc, VAT2(iz, 1, lev)),
                       RAT(  x, VAT2(iz, 1, lev)),
                       w2, w3, w1,
                       &nuuu, &iters_s, &errtol_s, omega,
                       &iresid, &iadjoint, mgsmoo);
            }

        // End of cycling down to coarse grid loop
        }



        /* *********************************************************************
        * begin coarse grid
        * *********************************************************************/

        // Coarsest level
        level = *nlev;
        lev   = (*ilev - 1) + level;

        // Solve on coarsest grid with ncghs, mgsmoo_s=4 (no residual)
        iresid   = 0;
        iadjoint = 0;
        itmax_s  = 100;
        iters_s  = 0;
        errtol_s = *epsiln;
        mgsmoo_s = *mgsmoo;
        Vazeros(&nxf, &nyf, &nzf, RAT(x, VAT2(iz, 1, lev)));
        Vnsmooth(&nxf, &nyf, &nzf,
                RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                RAT( fc, VAT2(iz, 1, lev)),
                RAT(  x, VAT2(iz, 1, lev)),
                w1, w2, w3,
                &itmax_s, &iters_s, &errtol_s, omega,
                &iresid, &iadjoint, &mgsmoo_s);

        /* *********************************************************************
         * begin cycling back to fine grid
         * ********************************************************************/

        // Move up grids: interpolate resid to finer and smooth
        for (level = *nlev - 1; level >= 1; level--)  {
            lev = (*ilev - 1) + level;

            // Find new grid size
            numlev = 1;
            Vmkfine(&numlev, &nxf, &nyf, &nzf, &nxc, &nyc, &nzc);

            // Form difference of new approx at the coarse grid
            alpha = -1.0;
            Vxaxpy(&nxf, &nyf, &nzf, &alpha,
                RAT(w4, VAT2(iz, 1, lev + 1)), RAT(x, VAT2(iz, 1, lev + 1)));

            // Call the line search (on the coarser level)
            Vlinesearch(&nxf, &nyf, &nzf, &xdamp,
                    RAT(ipc, VAT2(iz, 5, lev + 1)),
                    RAT(rpc, VAT2(iz, 6, lev + 1)),
                    RAT( ac, VAT2(iz, 7, lev + 1)),
                    RAT( cc, VAT2(iz, 1, lev + 1)),
                    RAT( fc, VAT2(iz, 1, lev + 1)),
                    RAT(  x, VAT2(iz, 1, lev + 1)),
                    RAT( w4, VAT2(iz, 1, lev + 1)),
                    RAT( w0, VAT2(iz, 1, lev + 1)),
                    w1, w2, w3);

            // Interpolate to next finer grid
            VinterpPMG(&nxf, &nyf, &nzf, &nxc, &nyc, &nzc,
                    RAT(x, VAT2(iz, 1, lev + 1)),
                    w1,
                    RAT(pc, VAT2(iz, 11, lev)));

            // New grid size
            nxf = nxc;
            nyf = nyc;
            nzf = nzc;

            // Perform the coarse grid correction
            Vxaxpy(&nxf, &nyf, &nzf, &xdamp, w1, RAT(x, VAT2(iz, 1, lev)));

            // nu2 post-smoothings for correction (no residual) ***
            iresid = 0;
            iadjoint = 1;
            iters_s  = 0;
            errtol_s = 0.0;
            nuuu = Vivariv(nu2, &lev);
            Vnsmooth(&nxf, &nyf, &nzf,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                    RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                    RAT( fc, VAT2(iz, 1, lev)),
                    RAT(  x, VAT2(iz, 1, lev)),
                    w1, w2, w3,
                    &nuuu, &iters_s, &errtol_s, omega,
                    &iresid, &iadjoint, mgsmoo);
        }


        /* *********************************************************************
         * iteration complete: do some i/o
         * ********************************************************************/

        // Increment the iteration counter
         iters = iters + 1;

        // Compute/check the current stopping test
        if (*iok != 0) {
            orsnrm = rsnrm;
            if (*istop == 0) {
                Vnmresid(&nxf, &nyf, &nzf,
                       RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                       RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                       RAT( fc, VAT2(iz, 1, lev)), RAT(  x, VAT2(iz, 1, lev)),
                       w1, w2);
                       rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
            } else if (*istop == 1) {
                Vnmresid(&nxf, &nyf, &nzf,
                       RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                       RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                       RAT( fc, VAT2(iz, 1, lev)), RAT(  x, VAT2(iz, 1, lev)),
                       w1, w2);
                       rsnrm = Vxnrm1(&nxf, &nyf, &nzf,w1);
            } else if (*istop == 2) {
                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1, lev)), w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1, lev)), w1);
                rsnrm = Vxnrm1(&nxf, &nyf, &nzf, w1);
                Vxcopy(&nxf, &nyf, &nzf,
                       RAT(  x, VAT2(iz, 1, lev)),
                       RAT(tru, VAT2(iz, 1, lev)));
            } else if (*istop == 3) {
                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1, lev)), w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1, lev)), w1);
                       rsnrm = Vxnrm2(&nxf, &nyf, &nzf,w1);
            } else if (*istop == 4) {
                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1, lev)), w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1, lev)), w1);
                       rsnrm = Vxnrm2(&nxf, &nyf, &nzf,w1);
            } else if (*istop == 5) {
                Vxcopy(&nxf, &nyf, &nzf, RAT(tru, VAT2(iz, 1, lev)), w1);
                alpha = -1.0;
                Vxaxpy(&nxf, &nyf, &nzf, &alpha, RAT(x, VAT2(iz, 1, lev)), w1);
                Vnmatvec(&nxf, &nyf, &nzf,
                       RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                       RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)),
                       w1, w2, w3);
                       rsnrm = VSQRT(Vxdot(&nxf, &nyf, &nzf,w1,w2));
            } else {
                VABORT_MSG1("Bad istop value: %d", *istop);
            }

            Vprtstp(*iok, *iters, rsnrm, rsden, orsnrm);

            if ((rsnrm / rsden) <= *errtol)
                break;
         }

         if (*iters >= *itmax) {
             *ierror = 1;
             break;
         }
    }
}
