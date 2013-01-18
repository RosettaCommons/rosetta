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

#include "mikpckd.h"

VPUBLIC void Vxcopy(int *nx, int *ny, int *nz, double *x, double *y) {

    MAT3(x, *nx, *ny, *nz);
    MAT3(y, *nx, *ny, *nz);

    // The indices used to traverse the matrices
    int i, j, k;

    /// @todo:  Once the refactor begins, this will need to be unrolled
    #pragma omp parallel for private(i,j,k)
    for(k=2; k<=*nz-1; k++)
        for(j=2; j<=*ny-1; j++)
            for(i=2; i<=*nx-1; i++)
                VAT3(y, i, j, k) = VAT3(x, i, j, k);
}



VPUBLIC void Vxcopy_small(int *nx, int *ny, int *nz, double *x, double *y) {

    MAT3(x, *nx, *ny, *nz);
    MAT3(y, *nx - 2, *ny - 2, *nz - 2);

    // The indices used to traverse the matrices
    int i, j, k;

    for(k=2; k<=*nz-1; k++)
        for(j=2; j<=*ny-1; j++)
            for(i=2; i<=*nx-1; i++)
                VAT3(y, i - 1, j - 1, k - 1) = VAT3(x, i, j, k);
}



VPUBLIC void Vxcopy_large(int *nx, int *ny, int *nz, double *x, double *y) {

    /** @note This function is exactly equivalent to calling xcopy_small with
     *        the matrix arguments reversed.
     *  @todo Replace this.  it's worthless
     */

    MAT3(x, *nx - 2, *ny - 2, *nz - 2);
    MAT3(y, *nx, *ny, *nz);

    // The indices used to traverse the matrices
    int i, j, k;

    for(k=2; k<=*nz-1; k++)
        for(j=2; j<=*ny-1; j++)
            for(i=2; i<=*nx-1; i++)
                VAT3(y, i, j, k) = VAT3(x, i - 1, j - 1, k - 1);
}



VPUBLIC void Vxaxpy(int *nx, int *ny, int *nz,
        double *alpha, double *x, double *y) {

    // Create the wrappers
    MAT3(x, *nx, *ny, *nz);
    MAT3(y, *nx, *ny, *nz);

    // The indices used to traverse the matrices
    int i, j, k;

    /// @todo parallel optimization
    for(k=2; k<=*nz-1; k++)
        for(j=2; j<=*ny-1; j++)
            for(i=2; i<=*nx-1; i++)
                VAT3(y, i, j, k) += *alpha * VAT3(x, i, j, k);
}



VPUBLIC double Vxnrm1(int *nx, int *ny, int *nz,
        double *x) {

    double xnrm1 = 0.0;  ///< Accumulates the calculated normal value

    MAT3(x, *nx, *ny, *nz);

    // The indices used to traverse the matrices
    int i, j, k;

    /// @todo parallel optimization
    for(k=2; k<=*nz-1; k++)
        for(j=2; j<=*ny-1; j++)
            for(i=2; i<=*nx-1; i++)
                xnrm1 += VABS(VAT3(x, i, j, k));

    return xnrm1;
}



VPUBLIC double Vxnrm2(int *nx, int *ny, int *nz,
        double *x) {

    double xnrm2 = 0.0;  ///< Accumulates the calculated normal value

    MAT3(x, *nx, *ny, *nz);

    // The indices used to traverse the matrices
    int i, j, k;

    /// @todo parallel optimization
    for(k=2; k<=*nz-1; k++)
        for(j=2; j<=*ny-1; j++)
            for(i=2; i<=*nx-1; i++)
                xnrm2 += VAT3(x, i, j, k) * VAT3(x, i, j, k);

    return VSQRT(xnrm2);
}



VPUBLIC double Vxdot(int *nx, int *ny, int *nz,
        double *x, double *y) {

    int i, j, k;

    // Initialize
    double xdot = 0.0;

    MAT3(x, *nx, *ny, *nz);
    MAT3(y, *nx, *ny, *nz);

    // Do it
    for(k=2; k<=*nz-1; k++)
        for(j=2; j<=*ny-1; j++)
            for(i=2; i<=*nx-1; i++)
                xdot += VAT3(x, i, j, k) * VAT3(y, i, j, k);

    return xdot;
}



VPUBLIC void Vazeros(int *nx, int *ny, int *nz, double *x) {

    int i, n;
    int nproc = 1;

    n = *nx * *ny * *nz;

    #pragma omp parallel for private(i)
    for (i=1; i<=n; i++)
        VAT(x, i) = 0.0;
}



VPUBLIC void VfboundPMG(int *ibound, int *nx, int *ny, int *nz,
        double *x, double *gxc, double *gyc, double *gzc) {

    int i, j, k;

    // Create and bind the wrappers for the source data
    MAT3(  x, *nx, *ny, *nz);
    MAT3(gxc, *ny, *nz,   2);
    MAT3(gyc, *nx, *nz,   2);
    MAT3(gzc, *nx, *ny,   2);

    // Dirichlet test
    if (ibound == 0) {

        // Dero dirichlet
        VfboundPMG00(nx, ny, nz, x);

    } else {

        // Nonzero dirichlet

        // The (i=1) and (i=nx) boundaries
        for (k=1; k<=*nz; k++) {
            for (j=1; j<=*ny; j++) {
                   VAT3(x,   1, j, k) = VAT3(gxc, j, k, 1);
                   VAT3(x, *nx, j, k) = VAT3(gxc, j, k, 2);
            }
        }

        // The (j=1) and (j=ny) boundaries
        for (k=1; k<=*nz; k++) {
            for (i=1; i<=*nx; i++) {
                   VAT3(x, i,   1 ,k) = VAT3(gyc, i, k, 1);
                   VAT3(x, i, *ny, k) = VAT3(gyc, i, k, 2);
            }
        }

        // The (k=1) and (k=nz) boundaries
        for (j=1; j<=*ny; j++) {
            for (i=1; i<=*nx; i++) {
                   VAT3(x, i, j,   1) = VAT3(gzc, i, j, 1);
                   VAT3(x, i, j, *nz) = VAT3(gzc, i, j, 2);
            }
        }
    }
}



VPUBLIC void VfboundPMG00(int *nx, int *ny, int *nz, double *x) {

    int i, j, k;

    MAT3(  x, *nx, *ny, *nz);

    // The (i=1) and (i=nx) boundaries
    for (k=1; k<=*nz; k++) {
        for (j=1; j<=*ny; j++) {
            VAT3(x,   1, j, k) = 0.0;
            VAT3(x, *nx, j, k) = 0.0;
        }
    }

    // The (j=1) and (j=ny) boundaries
    for (k=1; k<=*nz; k++) {
        for(i=1; i<=*nx; i++) {
            VAT3(x, i,   1, k) = 0.0;
            VAT3(x, i, *ny, k) = 0.0;
        }
    }

    // The (k=1) and (k=nz) boundaries
    for (j=1; j<=*ny; j++) {
        for (i=1; i<=*nx; i++) {
            VAT3(x, i, j,   1) = 0.0;
            VAT3(x, i, j, *nz) = 0.0;
        }
    }
}



VPUBLIC void Vaxrand(int *nx, int *ny, int *nz, double *x) {

    int n, i, ii, ipara, ivect, iflag;
    int nproc = 1;
    double xdum;

    WARN_UNTESTED;

    // Find parallel loops (ipara), remainder (ivect)
    n = *nx * *ny * *nz;
    ipara = n / nproc;
    ivect = n % nproc;
    iflag = 1;
    xdum  = (double)(VRAND);

    // Do parallel loops
    for (ii=1; ii<=nproc; ii++)
        for (i=1+(ipara*(ii-1)); i<=ipara*ii; i++)
            VAT(x, i) = (double)(VRAND);

    // Do vector loops
    for (i=ipara*nproc+1; i<=n; i++)
        VAT(x, i) = (double)(VRAND);
}



VPUBLIC void Vxscal(int *nx, int *ny, int *nz, double *fac, double *x) {

    int i, j, k;

    MAT3(x, *nx, *ny, *nz);

    for (k=2; k<=*nz-1; k++)
        for (j=2; j<=*ny-1; j++)
            for (i=2; i<=*nx-1; i++)
                VAT3(x, i, j, k) *= *fac;
}



VPUBLIC void Vprtmatd(int *nx, int *ny, int *nz,
        int *ipc, double *rpc, double *ac) {

    int numdia;

    MAT2(ac, *nx * *ny * *nz, 1);

    WARN_UNTESTED;

    // Do the printing
    numdia = VAT(ipc, 11);
    if (numdia == 7) {
       Vprtmatd7(nx, ny, nz,
               ipc, rpc,
               RAT2(ac, 1, 1), RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4));
    } else if (numdia == 27) {
       Vprtmatd27(nx, ny, nz,
               ipc, rpc,
               RAT2(ac, 1,  1), RAT2(ac, 1,  2), RAT2(ac, 1,  3), RAT2(ac, 1,  4),
               RAT2(ac, 1,  5), RAT2(ac, 1,  6),
               RAT2(ac, 1,  7), RAT2(ac, 1,  8), RAT2(ac, 1,  9), RAT2(ac, 1, 10),
               RAT2(ac, 1, 11), RAT2(ac, 1, 12), RAT2(ac, 1, 13), RAT2(ac, 1, 14));
    } else {
       Vnm_print(2, "Vprtmatd: invalid stencil type given: %d\n", numdia);
    }
}



VPUBLIC void Vprtmatd7(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *oC, double *oE, double *oN, double *uC) {

    int n, i, j, k;

    MAT3(oC, *nx, *ny, *nz);
    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(uC, *nx, *ny, *nz);

    WARN_UNTESTED;

    // Recover matrix dimension
    n = (*nx - 2) * (*ny - 2) * (*nz - 2);

    Vnm_print(2, "Vprtmatd7: Dimension of matrix = %d\n", n);
    Vnm_print(2, "Vprtmatd7: Begin diagonal matrix\n");

    // Handle first block
    for (k=2; k<=*nz-1; k++)
        for (j=2; j<=*ny-1; j++)
            for (i=2; i<=*nx-1; i++)
                Vnm_print(2, "Vprtmatd7: (%d,%d,%d) - oC=%g, oE=%g, oN=%g, uC=%g\n",
                        VAT3(oC,i,j,k), VAT3(oE,i,j,k), VAT3(oN,i,j,k), VAT3(uC,i,j,k));

    // Finish up
    Vnm_print(2, "Vprtmatd7: End diagonal matrix\n");
}



VEXTERNC void Vprtmatd27(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *oC, double *oE, double *oN, double *uC,
        double *oNE, double *oNW,
        double *uE, double *uW, double *uN, double *uS,
        double *uNE, double *uNW, double *uSE, double *uSW) {

    int n, i, j, k;

    MAT3( oC, *nx, *ny, *nz);
    MAT3( oE, *nx, *ny, *nz);
    MAT3( oN, *nx, *ny, *nz);
    MAT3( uC, *nx, *ny, *nz);
    MAT3(oNE, *nx, *ny, *nz);
    MAT3(oNW, *nx, *ny, *nz);
    MAT3( uE, *nx, *ny, *nz);
    MAT3( uW, *nx, *ny, *nz);
    MAT3( uN, *nx, *ny, *nz);
    MAT3( uS, *nx, *ny, *nz);
    MAT3(uNE, *nx, *ny, *nz);
    MAT3(uNW, *nx, *ny, *nz);
    MAT3(uSE, *nx, *ny, *nz);
    MAT3(uSW, *nx, *ny, *nz);

    WARN_UNTESTED;

    // Recover matrix dimension
    n = (*nx - 2) * (*ny - 2) * (*nz - 2);

    Vnm_print(2, "Vprtmatd27: Dimension of matrix = %d\n", n);
    Vnm_print(2, "Vprtmatd27: Begin diagonal matrix\n");

    // Handle first block
    for (k=2; k<=*nz-1; k++)
        for (j=2; j<=*ny-1; j++)
            for (i=2; i<=*nx-1; i++)
                Vnm_print(2, "Vprtmatd7: (%d,%d,%d) -  oC = %g, oE = %g, "
                        "oNW = %g, oN = %g, oNE = %g, uSW = %g, uS = %g, "
                        "uSE = %g, uW = %g, uC = %g, uE = %g, uNW = %g, "
                        "uN = %g, uNE = %g\n",
                        VAT3( oC, i, j, k), VAT3( oE, i, j, k),
                        VAT3(oNW, i, j, k), VAT3( oN, i, j, k),
                        VAT3(oNE, i, j, k), VAT3(uSW, i, j, k),
                        VAT3( uS, i, j, k), VAT3(uSE, i, j, k),
                        VAT3( uW, i, j, k), VAT3( uC, i, j, k),
                        VAT3( uE, i, j, k), VAT3(uNW, i, j, k),
                        VAT3( uN, i, j, k), VAT3(uNE, i, j, k) );

    // Finish up
    Vnm_print(2, "Vprtmatd27: End diagonal matrix\n");
}

VPUBLIC void Vlinesearch(int *nx, int *ny, int *nz,
        double *alpha,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *p, double *x, double *r,
        double *ap, double *zk, double *zkp1) {
    VABORT_MSG0("Not translated yet");
}
