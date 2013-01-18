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

#include "powerd.h"

VPUBLIC void Vpower(int *nx, int *ny, int *nz,
        int *iz, int *ilev,
        int *ipc, double *rpc, double *ac, double *cc,
        double *w1, double *w2, double *w3, double *w4,
        double *eigmax, double *eigmax_model, double *tol,
        int *itmax, int *iters, int *iinfo) {

    int lev, level;
    double denom, fac, rho, oldrho, error, relerr;

    /// @todo  Just use a constant definition of PI here
    double pi = 4.0 * atan( 1.0 );

    // Utility variables
    int skipIters = 0;
    double alpha;

    MAT2(iz, 50, 1);

    WARN_UNTESTED;

    // Recover level information
    level = 1;
    lev   = (*ilev - 1) + level;

    // Seed vector: random to contain all components

    Vaxrand(nx, ny, nz, w1);

    Vazeros(nx, ny, nz, w2);
    Vazeros(nx, ny, nz, w3);
    Vazeros(nx, ny, nz, w4);

    // Compute raleigh quotient with the seed vector
    denom = Vxnrm2(nx, ny, nz, w1);
    fac = 1.0 / denom;
    Vxscal(nx, ny, nz, &fac, w1);
    Vmatvec(nx, ny, nz,
            RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6,lev)),
             RAT(ac, VAT2(iz, 7, lev)),  RAT(cc, VAT2(iz, 1,lev)), w1, w2);
    oldrho = Vxdot(nx, ny, nz, w1, w2);

    // I/O
    if (oldrho == 0.0) {
        if (*iinfo > 3)  {
            Vnm_print(2, "POWER: iter: estimate = %d %g\n", *iters, oldrho);
        }
        rho = oldrho;
    } else {

        // Main iteration
        *iters = 0;
        while(1) {
            (*iters)++;

            // Apply the matrix A
            Vmatvec(nx, ny, nz,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                     RAT(ac, VAT2(iz, 7, lev)),  RAT(cc, VAT2(iz, 1, lev)), w1, w2);

            Vxcopy(nx, ny, nz, w2, w1);

            // Normalize the new vector
            denom = Vxnrm2(nx, ny, nz, w1);
            fac = 1.0 / denom;
            Vxscal(nx, ny, nz, &fac, w1);

            // Compute the new raleigh quotient
            Vmatvec(nx, ny, nz,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                    RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)), w1, w2);
            rho = Vxdot(nx, ny, nz, w1, w2);

            // Stopping test ***
            // w2=A*x, w1=x, stop = 2-norm(A*x-lamda*x)

            Vxcopy(nx, ny, nz, w1, w3);
            Vxcopy(nx, ny, nz, w2, w4);
            Vxscal(nx, ny, nz, &rho, w3);
            alpha = -1.0;
            Vxaxpy(nx, ny, nz, &alpha, w3, w4);
            error = Vxnrm2(nx, ny, nz, w4);
            relerr = VABS(rho - oldrho ) / VABS( rho );

            // I/O
            if (*iinfo > 3) {

                Vnm_print(2, "POWER:  iters  =%d\n", *iters);
                Vnm_print(2, "        error  =%g\n", error);
                Vnm_print(2, "        relerr =%g\n", relerr);
                Vnm_print(2, "        rho    =%g\n", rho);
            }

            if( relerr < *tol || *iters == *itmax)
                break;

            oldrho = rho;
        }
    }

    // Return some stuff ***
    *eigmax = rho;
    fac = VPOW(2.0, *ilev - 1);
    *eigmax_model = fac * (6.0 - 2.0 * VCOS((*nx - 2) * pi / (*nx - 1))
                                     - 2.0 * VCOS((*ny - 2) * pi / (*ny - 1)));
}


VPUBLIC void Vipower(int *nx,int *ny,int *nz,
        double *u, int *iz,
        double *w0, double *w1, double *w2, double *w3, double *w4,
        double *eigmin, double *eigmin_model, double *tol,
        int *itmax, int *iters,
        int *nlev, int *ilev, int *nlev_real, int *mgsolv,
        int *iok, int *iinfo, double *epsiln, double *errtol, double *omega,
        int *nu1, int *nu2, int *mgsmoo,
        int *ipc, double *rpc,
        double *pc, double *ac, double *cc, double *tru) {

    int level, lev;
    double denom, fac, rho, oldrho;
    double error, relerr, errtol_s;
    int itmax_s, iters_s, ierror_s, iok_s, iinfo_s, istop_s;
    int nu1_s, nu2_s, mgsmoo_s;

    /// @todo  Just use a constant definition of PI here
    double pi = 4.0 * atan( 1.0 );

    // Utility variables
    double alpha;

    MAT2(iz, 50, 1);

    WARN_UNTESTED;

    // Recover level information
    level = 1;
    lev   = (*ilev - 1) + level;

    // Seed vector: random to contain all components
    Vaxrand(nx, ny, nz, w1);
    Vazeros(nx, ny, nz, w2);
    Vazeros(nx, ny, nz, w3);
    Vazeros(nx, ny, nz, w4);
    Vazeros(nx, ny, nz, RAT(w0, VAT2(iz, 1, lev)));
    Vazeros(nx, ny, nz, RAT( u, VAT2(iz, 1, lev)));

    // Compute raleigh quotient with the seed vector ***
    denom = Vxnrm2(nx, ny, nz, w1);
    fac = 1.0 / denom;
    Vxscal(nx, ny, nz, &fac, w1);
    Vmatvec(nx, ny, nz,
            RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
            RAT( ac, VAT2(iz, 7, lev)), RAT( cc, VAT2(iz, 1, lev)), w1, w2);
    oldrho = Vxdot(nx, ny, nz, w1, w2);

    // I/O
    if (oldrho == 0.0) {
           if (*iinfo > 3) {
               Vnm_print(2, "Vipower: iters=%d\n",    *iters);
               Vnm_print(2, "         estimate=%f\n", oldrho);
        }
       rho = oldrho;
    } else {

        //main iteration
        *iters = 0;
        while (1) {
               (*iters)++;

            // Apply the matrix A^{-1} (using MG solver)
               itmax_s = 100;
               iters_s = 0;
               ierror_s = 0;
               iok_s = 0;
               iinfo_s = 0;
               istop_s = 0;
            mgsmoo_s = 1;
            nu1_s = 1;
            nu2_s = 1;
            errtol_s = *epsiln;

            Vxcopy(nx, ny, nz, w1, RAT(w0, VAT2(iz, 1,lev)));
            Vmvcs(nx, ny, nz, u, iz,
                    w1, w2, w3, w4,
                    &istop_s, &itmax_s, &iters_s, &ierror_s,
                    nlev, ilev, nlev_real, mgsolv,
                    &iok_s, &iinfo_s, epsiln,
                    &errtol_s, omega, &nu1_s, &nu2_s, &mgsmoo_s,
                    ipc, rpc, pc, ac, cc, w0, tru);
            Vxcopy(nx, ny, nz, RAT(u, VAT2(iz, 1, lev)), w1);

            // Normalize the new vector
            denom = Vxnrm2(nx, ny, nz, w1);
            fac = 1.0 / denom;
            Vxscal(nx, ny, nz, &fac, w1);

            // Compute the new raleigh quotient
            Vmatvec(nx, ny, nz,
                    RAT(ipc, VAT2(iz, 5, lev)), RAT(rpc, VAT2(iz, 6, lev)),
                     RAT(ac, VAT2(iz, 7,lev)),   RAT(cc, VAT2(iz, 1, lev)), w1, w2);
            rho = Vxdot(nx, ny, nz, w1, w2);

            // Stopping test
            // w2=A*x, w1=x, stop = 2-norm(A*x-lamda*x) ***
            Vxcopy(nx, ny, nz, w1, w3);
            Vxcopy(nx, ny, nz, w2, w4);
            Vxscal(nx, ny, nz, &rho, w3);
            alpha = -1.0;
            Vxaxpy(nx, ny, nz, &alpha, w3, w4);
            error = Vxnrm2(nx, ny, nz, w4);
            relerr = VABS(rho - oldrho ) / VABS( rho );

            // I/O
            if (*iinfo > 3) {

                Vnm_print(2, "POWER:  iters  =%d\n", *iters);
                Vnm_print(2, "        error  =%g\n", error);
                Vnm_print(2, "        relerr =%g\n", relerr);
                Vnm_print(2, "        rho    =%g\n", rho);
            }

            if (relerr < *tol || *iters == *itmax)
                break;

            oldrho = rho;
        }
    }

    // Return some stuff
    *eigmin = rho;
    fac = VPOW(2.0, *ilev - 1);
    *eigmin_model = fac * (6.0 - 2.0 * VCOS(pi / (*nx - 1))
                               - 2.0 * VCOS(pi / (*ny - 1))
                               - 2.0 * VCOS(pi / (*nz - 1)));
}

VEXTERNC void Vmpower(int *nx, int *ny, int *nz,
        double *u, int *iz,
        double *w0, double *w1, double *w2, double *w3, double *w4,
        double *eigmax, double *tol,
        int *itmax, int *iters, int *nlev, int *ilev, int *nlev_real,
        int *mgsolv, int *iok, int *iinfo,
        double *epsiln, double *errtol, double *omega,
        int *nu1, int *nu2, int *mgsmoo, int *ipc, double *rpc,
        double *pc, double *ac, double *cc, double *fc, double *tru) {

    // Local variables
    int lev, level;
    double denom, fac, rho, oldrho, error;
    double relerr;
    int itmax_s, iters_s, ierror_s, iok_s, iinfo_s, istop_s;
    double alpha;

    MAT2(iz, 50, 1);

    // Recover level information
    level = 1;
    lev   = (*ilev - 1) + level;

    // Seed vector: random to contain all components
    Vaxrand(nx, ny, nz, w1);
    Vazeros(nx, ny, nz, w2);
    Vazeros(nx, ny, nz, w3);
    Vazeros(nx, ny, nz, w4);
    Vazeros(nx, ny, nz, RAT(u, VAT2(iz, 1, lev)));

    // NOTE: we destroy "fc" on this level due to lack of vectors... ***
    Vazeros(nx,ny,nz,RAT(fc, VAT2(iz, 1, lev)));

    // Normalize the seed vector
    denom = Vxnrm2(nx, ny, nz, w1);
    fac = 1.0 / denom;
    Vxscal(nx, ny, nz, &fac, w1);

    // Compute raleigh quotient with the seed vector
    Vxcopy(nx, ny, nz, w1, RAT(u, VAT2(iz, 1, lev)));
    itmax_s = 1;
    iters_s = 0;
    ierror_s = 0;
    iok_s = 0;
    iinfo_s = 0;
    istop_s = 1;
    Vmvcs(nx, ny, nz,
            u, iz, w0, w2, w3, w4,
            &istop_s, &itmax_s, &iters_s, &ierror_s,
            nlev, ilev, nlev_real, mgsolv,
            &iok_s, &iinfo_s,
            epsiln, errtol, omega, nu1, nu2, mgsmoo,
            ipc, rpc,
            pc, ac, cc, fc, tru);
    oldrho = Vxdot(nx, ny, nz, w1, RAT(u, VAT2(iz, 1, lev)));

    // I/O
    if (oldrho == 0.0) {
       if (*iinfo > 3) {
           Vnm_print(2, "Vmp0ower: iter=%d, estimate=%f", *iters, oldrho);
       }
       rho = oldrho;

    } else {

        // Main iteration
        *iters = 0;
        while (1) {
            (*iters)++;

            // Apply the matrix M
           Vxcopy(nx, ny, nz, w1, RAT(u, VAT2(iz, 1, lev)));
           itmax_s = 1;
           iters_s = 0;
           ierror_s = 0;
           iok_s = 0;
           iinfo_s = 0;
           istop_s = 1;
           Vmvcs(nx, ny, nz,
                   u, iz, w1, w2, w3, w4,
                   &istop_s, &itmax_s, &iters_s, &ierror_s,
                   nlev, ilev, nlev_real, mgsolv,
                   &iok_s, &iinfo_s,
                   epsiln, errtol, omega, nu1, nu2, mgsmoo,
                   ipc, rpc,
                   pc, ac, cc, fc, tru);
           Vxcopy(nx, ny, nz, RAT(u, VAT2(iz, 1, lev)), w1);

           // Normalize the new vector
           denom = Vxnrm2(nx, ny, nz, w1);
           fac = 1.0 / denom;
           Vxscal(nx, ny, nz, &fac, w1);

           // Compute the new raleigh quotient
           Vxcopy(nx, ny, nz, w1, RAT(u, VAT2(iz, 1, lev)));
           itmax_s = 1;
           iters_s = 0;
           ierror_s = 0;
           iok_s = 0;
           iinfo_s = 0;
           istop_s = 1;
           Vmvcs(nx, ny, nz,
                   u, iz, w0, w2, w3, w4,
                   &istop_s, &itmax_s, &iters_s, &ierror_s,
                   nlev, ilev, nlev_real, mgsolv,
                   &iok_s, &iinfo_s,
                   epsiln, errtol, omega, nu1, nu2, mgsmoo,
                   ipc, rpc,
                   pc, ac, cc, fc, tru);
           Vxcopy(nx, ny, nz, RAT(u, VAT2(iz, 1, lev)), w2);
           rho = Vxdot(nx, ny, nz, w1, w2);

           // Stopping test
           // w2=A*x, w1=x, stop = 2-norm(A*x-lamda*x)
           alpha = -1.0;
           Vxcopy(nx, ny, nz, w1, w3);
           Vxcopy(nx, ny, nz, w2, w4);
           Vxscal(nx, ny, nz, &rho, w3);
           Vxaxpy(nx, ny, nz, &alpha, w3, w4);
           error = Vxnrm2(nx, ny, nz, w4);
           relerr = VABS( rho - oldrho ) / VABS( rho );

           // I/O
           if (*iinfo > 3) {
               Vnm_print(2, "Vmpower: iter=%d; error=%f; relerr=%f; estimate=%f",
                       *iters, error, relerr, rho);
           }

           if ((relerr < *tol) || (*iters == *itmax)) {
               break;
           }

           oldrho = rho;
        }
    }

    *eigmax = rho;
}
