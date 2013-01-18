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

#include "matvecd.h"

VPUBLIC void Vmatvec(int *nx, int *ny, int *nz,
        int    *ipc, double *rpc,
        double  *ac, double  *cc,
        double   *x, double   *y) {

    int numdia;

    // Do in one step
    numdia = VAT(ipc, 11);

    if (numdia == 7) {
        Vmatvec7(nx, ny, nz,
                ipc, rpc,
                ac, cc,
                x, y);
    } else  if (numdia == 27) {
        Vmatvec27(nx, ny, nz,
                ipc, rpc,
                ac, cc,
                x, y);
    } else {
        Vnm_print(2, "MATVEC: invalid stencil type given...");
    }
}

VPUBLIC void Vmatvec7(int *nx, int *ny, int *nz,
        int    *ipc, double *rpc,
        double  *ac, double  *cc,
        double   *x, double   *y) {

    MAT2(ac, *nx * *ny * *nz, 1);

    Vmatvec7_1s(nx, ny, nz,
                ipc,     rpc,
            RAT2(ac, 1, 1),      cc,
            RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4),
                  x,      y);
}



VEXTERNC void Vmatvec7_1s(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *oC, double *cc,
        double *oE, double *oN, double *uC,
        double  *x, double  *y) {

    int i, j, k;

    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(uC, *nx, *ny, *nz);
    MAT3(cc, *nx, *ny, *nz);
    MAT3(oC, *nx, *ny, *nz);
    MAT3(x, *nx, *ny, *nz);
    MAT3(y, *nx, *ny, *nz);

    // Do it
    #pragma omp parallel for private(i, j, k)
    for (k=2; k<=*nz-1; k++) {
        for (j=2; j<=*ny-1; j++) {
            for(i=2; i<=*nx-1; i++) {
                VAT3(y, i, j, k) =
                           - VAT3( oN,   i,   j,   k)                * VAT3(x,   i, j+1,  k)
                           - VAT3( oN,   i, j-1,   k)                * VAT3(x,   i, j-1,  k)
                           - VAT3( oE,   i,   j,   k)                * VAT3(x, i+1,   j,  k)
                           - VAT3( oE, i-1,   j,   k)                * VAT3(x, i-1,   j,  k)
                           - VAT3( uC,   i,   j, k-1)                * VAT3(x,   i,   j,k-1)
                           - VAT3( uC,   i,   j,   k)                * VAT3(x,   i,   j,k+1)
                           + (VAT3(oC,   i,   j,   k) + VAT3(cc, i, j, k)) * VAT3(x,   i,   j,  k);
            }
        }
    }
}



VPUBLIC void Vmatvec27(int *nx, int *ny, int *nz,
        int    *ipc, double *rpc,
        double  *ac, double  *cc,
        double   *x, double   *y) {

    MAT2(ac, *nx * *ny * *nz, 1);

    Vmatvec27_1s(nx, ny, nz,
            ipc, rpc,
            RAT2(ac, 1, 1), cc,
            RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4),
            RAT2(ac, 1, 5), RAT2(ac, 1, 6),
            RAT2(ac, 1, 7), RAT2(ac, 1, 8), RAT2(ac, 1, 9), RAT2(ac, 1,10),
            RAT2(ac, 1,11), RAT2(ac, 1,12), RAT2(ac, 1,13), RAT2(ac, 1,14),
             x,        y);
}



VPUBLIC void Vmatvec27_1s(int *nx, int *ny, int *nz,
        int    *ipc, double *rpc,
        double  *oC, double  *cc,
        double  *oE, double  *oN, double  *uC,
        double *oNE, double *oNW,
        double  *uE, double  *uW, double  *uN, double  *uS,
        double *uNE, double *uNW, double *uSE, double *uSW,
        double   *x, double   *y) {

    int i, j, k;

    double tmpO, tmpU, tmpD;

    MAT3(cc, *nx, *ny, *nz);
    MAT3(x, *nx, *ny, *nz);
    MAT3(y, *nx, *ny, *nz);

    MAT3(oC, *nx, *ny, *nz);
    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(oNE, *nx, *ny, *nz);
    MAT3(oNW, *nx, *ny, *nz);

    MAT3(uC, *nx, *ny, *nz);
    MAT3(uE, *nx, *ny, *nz);
    MAT3(uW, *nx, *ny, *nz);
    MAT3(uN, *nx, *ny, *nz);
    MAT3(uS, *nx, *ny, *nz);
    MAT3(uNE, *nx, *ny, *nz);
    MAT3(uNW, *nx, *ny, *nz);
    MAT3(uSE, *nx, *ny, *nz);
    MAT3(uSW, *nx, *ny, *nz);

    // Do it
    #pragma omp parallel for private(i, j, k, tmpO, tmpU, tmpD)
    for (k=2; k<=*nz-1; k++) {
        for (j=2; j<=*ny-1; j++) {
            for(i=2; i<=*nx-1; i++) {
                tmpO =
                     - VAT3(  oN,   i,   j,   k) * VAT3(x,   i, j+1,   k)
                     - VAT3(  oN,   i, j-1,   k) * VAT3(x,   i, j-1,   k)
                     - VAT3(  oE,   i,   j,   k) * VAT3(x, i+1,   j,   k)
                     - VAT3(  oE, i-1,   j,   k) * VAT3(x, i-1,   j,   k)
                     - VAT3( oNE,   i,   j,   k) * VAT3(x, i+1, j+1,   k)
                     - VAT3( oNW,   i,   j,   k) * VAT3(x, i-1, j+1,   k)
                     - VAT3( oNW, i+1, j-1,   k) * VAT3(x, i+1, j-1,   k)
                     - VAT3( oNE, i-1, j-1,   k) * VAT3(x, i-1, j-1,   k);

                tmpU =
                     - VAT3(  uC,   i,   j,   k) * VAT3(x,   i,   j, k+1)
                     - VAT3(  uN,   i,   j,   k) * VAT3(x,   i, j+1, k+1)
                     - VAT3(  uS,   i,   j,   k) * VAT3(x,   i, j-1, k+1)
                     - VAT3(  uE,   i,   j,   k) * VAT3(x, i+1,   j, k+1)
                     - VAT3(  uW,   i,   j,   k) * VAT3(x, i-1,   j, k+1)
                     - VAT3( uNE,   i,   j,   k) * VAT3(x, i+1, j+1, k+1)
                     - VAT3( uNW,   i,   j,   k) * VAT3(x, i-1, j+1, k+1)
                     - VAT3( uSE,   i,   j,   k) * VAT3(x, i+1, j-1, k+1)
                     - VAT3( uSW,   i,   j,   k) * VAT3(x, i-1, j-1, k+1);

                tmpD =
                     - VAT3(  uC,   i,   j, k-1) * VAT3(x,   i,   j, k-1)
                     - VAT3(  uS,   i, j+1, k-1) * VAT3(x,   i, j+1, k-1)
                     - VAT3(  uN,   i, j-1, k-1) * VAT3(x,   i, j-1, k-1)
                     - VAT3(  uW, i+1,   j, k-1) * VAT3(x, i+1,   j, k-1)
                     - VAT3(  uE, i-1,   j, k-1) * VAT3(x, i-1,   j, k-1)
                     - VAT3( uSW, i+1, j+1, k-1) * VAT3(x, i+1, j+1, k-1)
                     - VAT3( uSE, i-1, j+1, k-1) * VAT3(x, i-1, j+1, k-1)
                     - VAT3( uNW, i+1, j-1, k-1) * VAT3(x, i+1, j-1, k-1)
                     - VAT3( uNE, i-1, j-1, k-1) * VAT3(x, i-1, j-1, k-1);

                VAT3(y, i, j, k) = tmpO + tmpU + tmpD
                           + (VAT3(oC, i, j, k) + VAT3(cc, i, j, k)) * VAT3(x, i, j, k);
            }
        }
    }
}



VEXTERNC void Vnmatvec(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *x, double *y, double *w1) {

    int numdia;

    // Do in one step
    numdia = VAT(ipc, 11);

    if (numdia == 7) {
        Vnmatvec7(nx, ny, nz,
                ipc, rpc,
                ac, cc,
                x, y, w1);
    } else  if (numdia == 27) {
        Vnmatvec27(nx, ny, nz,
                ipc, rpc,
                ac, cc,
                x, y, w1);
    } else {
        Vnm_print(2, "MATVEC: invalid stencil type given...");
    }
}



VPUBLIC void Vnmatvec7(int *nx, int *ny, int *nz,
        int    *ipc, double *rpc,
        double  *ac, double  *cc,
        double   *x, double   *y, double *w1) {

    MAT2(ac, *nx * *ny * *nz, 1);

    WARN_UNTESTED;

    Vnmatvecd7_1s(nx, ny, nz,
            ipc, rpc,
            RAT2(ac, 1, 1), cc,
            RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4),
            x, y, w1);
}



VPUBLIC void Vnmatvecd7_1s(int *nx, int *ny, int *nz,
        int    *ipc, double *rpc,
        double  *oC, double  *cc,
        double  *oE, double  *oN, double *uC,
        double   *x, double   *y, double *w1) {

    int i, j, k;
    int ipkey;

    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(uC, *nx, *ny, *nz);
    MAT3(cc, *nx, *ny, *nz);
    MAT3(oC, *nx, *ny, *nz);
    MAT3( x, *nx, *ny, *nz);
    MAT3( y, *nx, *ny, *nz);
    MAT3(w1, *nx, *ny, *nz);

    WARN_UNTESTED;

    // first get vector nonlinear term to avoid subroutine calls
    ipkey = VAT(ipc, 10);
    Vc_vec(cc, x, w1, nx, ny, nz, &ipkey);

    // The operator
    #pragma omp parallel for private(i, j, k)
    for (k=2; k<=*nz-1; k++)
        for (j=2; j<=*ny-1; j++)
            for(i=2; i<=*nx-1; i++)
                VAT3(y, i, j, k) =
                        -  VAT3(oN,   i,   j,   k) * VAT3(x,   i, j+1,   k)
                        -  VAT3(oN,   i, j-1,   k) * VAT3(x,   i, j-1,   k)
                        -  VAT3(oE,   i,   j,   k) * VAT3(x, i+1,   j,   k)
                        -  VAT3(oE, i-1,   j,   k) * VAT3(x, i-1,   j,   k)
                        -  VAT3(uC,   i,   j, k-1) * VAT3(x,   i,   j, k-1)
                        -  VAT3(uC,   i,   j,   k) * VAT3(x,   i,   j, k+1)
                        +  VAT3(oC,   i,   j,   k) * VAT3(x,   i,   j,   k)
                        +  VAT3(w1,   i,   j,   k);
}


VPUBLIC void Vnmatvec27(int *nx, int *ny, int *nz,
        int    *ipc, double *rpc,
        double  *ac, double  *cc,
        double   *x, double   *y, double *w1) {

    MAT2(ac, *nx * *ny * *nz, 1);

    WARN_UNTESTED;

    // Do in one step
    Vnmatvecd27_1s(nx, ny, nz,
            ipc, rpc,
            RAT2(ac, 1, 1), cc,
            RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4),
            RAT2(ac, 1, 5), RAT2(ac, 1, 6),
            RAT2(ac, 1, 7), RAT2(ac, 1, 8), RAT2(ac, 1, 9), RAT2(ac, 1,10),
            RAT2(ac, 1,11), RAT2(ac, 1,12), RAT2(ac, 1,13), RAT2(ac, 1,14),
            x, y, w1);
}



VPUBLIC void Vnmatvecd27_1s(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *oC, double *cc,
        double *oE, double *oN, double *uC,
        double *oNE, double *oNW,
        double *uE, double *uW, double *uN, double *uS,
        double *uNE, double *uNW, double *uSE, double *uSW,
        double *x, double *y, double *w1) {

    int i, j, k;
    int ipkey;

    double tmpO, tmpU, tmpD;

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
    MAT3( cc, *nx, *ny, *nz);
    MAT3( oC, *nx, *ny, *nz);
    MAT3(  x, *nx, *ny, *nz);
    MAT3(  y, *nx, *ny, *nz);
    MAT3( w1, *nx, *ny, *nz);

    WARN_UNTESTED;

    // First get vector noNlinear term to avoid subroutine calls
    ipkey = VAT(ipc, 10);
    Vc_vec(cc, x, w1, nx, ny, nz, &ipkey);

    // The operator
    #pragma omp parallel for private(i, j, k, tmpO, tmpU, tmpD)
    for (k=2; k<=*nz-1; k++) {
        for (j=2; j<=*ny-1; j++) {
            for(i=2; i<=*nx-1; i++) {

                tmpO =
                  -  VAT3( oN,   i,   j,   k) * VAT3(x,   i, j+1,   k)
                  -  VAT3( oN,   i, j-1,   k) * VAT3(x,   i, j-1,   k)
                  -  VAT3( oE,   i,   j,   k) * VAT3(x, i+1,   j,   k)
                  -  VAT3( oE, i-1,   j,   k) * VAT3(x, i-1,   j,   k)
                  -  VAT3(oNE,   i,   j,   k) * VAT3(x, i+1, j+1,   k)
                  -  VAT3(oNW,   i,   j,   k) * VAT3(x, i-1, j+1,   k)
                  -  VAT3(oNW, i+1, j-1,   k) * VAT3(x, i+1, j-1,   k)
                  -  VAT3(oNE, i-1, j-1,   k) * VAT3(x, i-1, j-1,   k);

                tmpU =
                  -  VAT3( uC,   i,   j,   k) * VAT3(x,   i,   j, k+1)
                  -  VAT3( uN,   i,   j,   k) * VAT3(x,   i, j+1, k+1)
                  -  VAT3( uS,   i,   j,   k) * VAT3(x,   i, j-1, k+1)
                  -  VAT3( uE,   i,   j,   k) * VAT3(x, i+1,   j, k+1)
                  -  VAT3( uW,   i,   j,   k) * VAT3(x, i-1,   j, k+1)
                  -  VAT3(uNE,   i,   j,   k) * VAT3(x, i+1, j+1, k+1)
                  -  VAT3(uNW,   i,   j,   k) * VAT3(x, i-1, j+1, k+1)
                  -  VAT3(uSE,   i,   j,   k) * VAT3(x, i+1, j-1, k+1)
                  -  VAT3(uSW,   i,   j,   k) * VAT3(x, i-1, j-1, k+1);

                tmpD =
                  -  VAT3( uC,   i,   j, k-1) * VAT3(x,   i,   j, k-1)
                  -  VAT3( uS,   i, j+1, k-1) * VAT3(x,   i, j+1, k-1)
                  -  VAT3( uN,   i, j-1, k-1) * VAT3(x,   i, j-1, k-1)
                  -  VAT3( uW, i+1,   j, k-1) * VAT3(x, i+1,   j, k-1)
                  -  VAT3( uE, i-1,   j, k-1) * VAT3(x, i-1,   j, k-1)
                  -  VAT3(uSW, i+1, j+1, k-1) * VAT3(x, i+1, j+1, k-1)
                  -  VAT3(uSE, i-1, j+1, k-1) * VAT3(x, i-1, j+1, k-1)
                  -  VAT3(uNW, i+1, j-1, k-1) * VAT3(x, i+1, j-1, k-1)
                  -  VAT3(uNE, i-1, j-1, k-1) * VAT3(x, i-1, j-1, k-1);

                VAT3(y, i, j, k) = tmpO + tmpU + tmpD
                        + VAT3(oC, i, j, k) * VAT3(x, i, j, k)
                        + VAT3(w1, i, j, k);
            }
        }
    }
}



VPUBLIC void Vmresid(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *r) {

    int numdia;

    // Do in one step
    numdia = VAT(ipc, 11);
    if (numdia == 7) {
         Vmresid7(nx, ny, nz, ipc, rpc, ac, cc, fc, x, r);
   } else if (numdia == 27) {
         Vmresid27(nx, ny, nz, ipc, rpc, ac, cc, fc, x, r);
   } else {
       Vnm_print(2, "Vmresid: invalid stencil type given...\n");
   }
}



VPUBLIC void Vmresid7(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *r) {

    MAT2(ac, *nx * *ny * *nz, 1);

    // Do in one step
    Vmresid7_1s(nx, ny, nz,
            ipc, rpc,
            RAT2(ac, 1,1), cc, fc,
            RAT2(ac, 1,2), RAT2(ac, 1,3), RAT2(ac, 1,4),
            x,r);
}

VPUBLIC void Vmresid7_1s(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *oC, double *cc, double *fc,
        double *oE, double *oN, double *uC,
        double *x, double *r) {

    int i, j, k;

    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(uC, *nx, *ny, *nz);
    MAT3(cc, *nx, *ny, *nz);
    MAT3(fc, *nx, *ny, *nz);
    MAT3(oC, *nx, *ny, *nz);
    MAT3(x, *nx, *ny, *nz);
    MAT3(r, *nx, *ny, *nz);

    // Do it
    #pragma omp parallel for private(i, j, k)
    for (k=2; k<=*nz-1; k++) {
        for (j=2; j<=*ny-1; j++) {
            for(i=2; i<=*nx-1; i++) {
                VAT3(r, i,j,k) =  VAT3(fc,   i,   j,   k)
                         + VAT3( oN,   i,   j,   k)                * VAT3(x,   i, j+1,   k)
                         + VAT3( oN,   i, j-1,   k)                * VAT3(x,   i, j-1,   k)
                         + VAT3( oE,   i,   j,   k)                * VAT3(x, i+1,   j,   k)
                         + VAT3( oE, i-1,   j,   k)                * VAT3(x, i-1,   j,   k)
                         + VAT3( uC,   i,   j, k-1)                * VAT3(x,   i,   j, k-1)
                         + VAT3( uC,   i,   j,   k)                * VAT3(x,   i,   j, k+1)
                         - (VAT3(oC,   i,   j,   k) + VAT3(cc, i, j, k)) * VAT3(x,   i,   j,   k);
            }
        }
    }
}



VPUBLIC void Vmresid27(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *r) {

    MAT2(ac, *nx * *ny * *nz, 1);

    // Do in one step
    Vmresid27_1s(nx,ny,nz,
            ipc, rpc,
            RAT2(ac, 1, 1),       cc,       fc,
            RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4),
            RAT2(ac, 1, 5), RAT2(ac, 1, 6),
            RAT2(ac, 1, 7), RAT2(ac, 1, 8), RAT2(ac, 1, 9), RAT2(ac, 1,10),
            RAT2(ac, 1,11), RAT2(ac, 1,12), RAT2(ac, 1,13), RAT2(ac, 1,14),
            x,r);
}



VPUBLIC void Vmresid27_1s(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double  *oC, double  *cc, double  *fc,
        double  *oE, double  *oN, double  *uC,
        double *oNE, double *oNW,
        double  *uE, double  *uW, double  *uN, double  *uS,
        double *uNE, double *uNW, double *uSE, double *uSW,
        double *x, double *r) {

    int i, j, k;

    double tmpO, tmpU, tmpD;

    MAT3(cc, *nx, *ny, *nz);
    MAT3(fc, *nx, *ny, *nz);
    MAT3(x, *nx, *ny, *nz);
    MAT3(r, *nx, *ny, *nz);

    MAT3(oC, *nx, *ny, *nz);
    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(oNE, *nx, *ny, *nz);
    MAT3(oNW, *nx, *ny, *nz);

    MAT3(uC, *nx, *ny, *nz);
    MAT3(uE, *nx, *ny, *nz);
    MAT3(uW, *nx, *ny, *nz);
    MAT3(uN, *nx, *ny, *nz);
    MAT3(uS, *nx, *ny, *nz);
    MAT3(uNE, *nx, *ny, *nz);
    MAT3(uNW, *nx, *ny, *nz);
    MAT3(uSE, *nx, *ny, *nz);
    MAT3(uSW, *nx, *ny, *nz);

    #pragma omp parallel for private(i, j, k, tmpO, tmpU, tmpD)
    for (k=2; k<=*nz-1; k++) {
        for (j=2; j<=*ny-1; j++) {
            for(i=2; i<=*nx-1; i++) {

                tmpO  =
                        + VAT3(  oN,   i,   j,   k) * VAT3(x,   i, j+1,   k)
                        + VAT3(  oN,   i, j-1,   k) * VAT3(x,   i, j-1,   k)
                        + VAT3(  oE,   i,   j,   k) * VAT3(x, i+1,   j,   k)
                        + VAT3(  oE, i-1,   j,   k) * VAT3(x, i-1,   j,   k)
                        + VAT3( oNE,   i,   j,   k) * VAT3(x, i+1, j+1,   k)
                        + VAT3( oNW,   i,   j,   k) * VAT3(x, i-1, j+1,   k)
                        + VAT3( oNW, i+1, j-1,   k) * VAT3(x, i+1, j-1,   k)
                        + VAT3( oNE, i-1, j-1,   k) * VAT3(x, i-1, j-1,   k);

                tmpU =
                        + VAT3(  uC,   i,   j,   k) * VAT3(x,   i,   j, k+1)
                        + VAT3(  uN,   i,   j,   k) * VAT3(x,   i, j+1, k+1)
                        + VAT3(  uS,   i,   j,   k) * VAT3(x,   i, j-1, k+1)
                        + VAT3(  uE,   i,   j,   k) * VAT3(x, i+1,   j, k+1)
                        + VAT3(  uW,   i,   j,   k) * VAT3(x, i-1,   j, k+1)
                        + VAT3( uNE,   i,   j,   k) * VAT3(x, i+1, j+1, k+1)
                        + VAT3( uNW,   i,   j,   k) * VAT3(x, i-1, j+1, k+1)
                        + VAT3( uSE,   i,   j,   k) * VAT3(x, i+1, j-1, k+1)
                        + VAT3( uSW,   i,   j,   k) * VAT3(x, i-1, j-1, k+1);

                tmpD =
                        + VAT3(  uC,   i,   j, k-1) * VAT3(x,   i,   j, k-1)
                        + VAT3(  uS,   i, j+1, k-1) * VAT3(x,   i, j+1, k-1)
                        + VAT3(  uN,   i, j-1, k-1) * VAT3(x,   i, j-1, k-1)
                        + VAT3(  uW, i+1,   j, k-1) * VAT3(x, i+1,   j, k-1)
                        + VAT3(  uE, i-1,   j, k-1) * VAT3(x, i-1,   j, k-1)
                        + VAT3( uSW, i+1, j+1, k-1) * VAT3(x, i+1, j+1, k-1)
                        + VAT3( uSE, i-1, j+1, k-1) * VAT3(x, i-1, j+1, k-1)
                        + VAT3( uNW, i+1, j-1, k-1) * VAT3(x, i+1, j-1, k-1)
                        + VAT3( uNE, i-1, j-1, k-1) * VAT3(x, i-1, j-1, k-1);

                VAT3(r, i, j, k) =  VAT3(fc, i, j, k) + tmpO + tmpU + tmpD
                           - (VAT3(oC, i, j, k) + VAT3(cc, i, j, k)) * VAT3(x, i, j, k);
            }
        }
    }
}



VPUBLIC void Vnmresid(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *r, double *w1) {

    int numdia;

    // Do in oNe step ***
    numdia = VAT(ipc, 11);
    if (numdia == 7) {
        Vnmresid7(nx, ny, nz, ipc, rpc, ac, cc, fc, x, r, w1);
    } else if (numdia == 27) {
        Vnmresid27(nx, ny, nz, ipc, rpc, ac, cc, fc, x, r, w1);
    } else {
        Vnm_print(2, "Vnmresid: invalid stencil type given...\n");
    }
}



VPUBLIC void Vnmresid7(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *r, double *w1) {

    MAT2(ac, *nx * *ny * *nz, 1);

    // Do in oNe step
    Vnmresid7_1s(nx, ny, nz,
            ipc, rpc,
            RAT2(ac, 1, 1), cc, fc,
            RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4),
            x, r, w1);
}

VPUBLIC void Vnmresid7_1s(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *oC, double *cc, double *fc,
        double *oE, double *oN, double *uC,
        double *x, double *r, double *w1) {

    int i, j, k;
    int ipkey;

    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(uC, *nx, *ny, *nz);
    MAT3(cc, *nx, *ny, *nz);
    MAT3(fc, *nx, *ny, *nz);
    MAT3(oC, *nx, *ny, *nz);
    MAT3( x, *nx, *ny, *nz);
    MAT3( r, *nx, *ny, *nz);
    MAT3(w1, *nx, *ny, *nz);

    // First get vector nonlinear term to avoid subroutine calls
    ipkey = VAT(ipc, 10);
    Vc_vec(cc, x, w1, nx, ny, nz, &ipkey);

    // The residual
    for (k=2; k<=*nz-1; k++) {
        for (j=2; j<=*ny-1; j++) {
            for (i=2; i<=*nx-1; i++) {
                VAT3(r, i, j, k) = VAT3(fc,   i,   j,   k)
                                  + VAT3(oN,   i,   j,   k) * VAT3(x,   i, j+1,   k)
                                  + VAT3(oN,   i, j-1,   k) * VAT3(x,   i, j-1,   k)
                                  + VAT3(oE,   i,   j,   k) * VAT3(x, i+1,   j,   k)
                                  + VAT3(oE, i-1,   j,   k) * VAT3(x, i-1,   j,   k)
                                  + VAT3(uC,   i,   j, k-1) * VAT3(x,   i,   j, k-1)
                                  + VAT3(uC,   i,   j,   k) * VAT3(x,   i,   j, k+1)
                                  - VAT3(oC,   i,   j,   k) * VAT3(x,   i,   j,   k)
                                  - VAT3(w1,   i,   j,   k);
            }
        }
    }
}



VPUBLIC void Vnmresid27(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *r, double *w1) {

    MAT2(ac, *nx * *ny * *nz, 1);

    // Do in oNe step
    Vnmresid27_1s(nx, ny, nz,
            ipc, rpc,
            RAT2(ac, 1,  1), cc, fc,
            RAT2(ac, 1,  2), RAT2(ac, 1,  3), RAT2(ac, 1,  4),
            RAT2(ac, 1,  5), RAT2(ac, 1,  6),
            RAT2(ac, 1,  7), RAT2(ac, 1,  8), RAT2(ac, 1,  9), RAT2(ac, 1, 10),
            RAT2(ac, 1, 11), RAT2(ac, 1, 12), RAT2(ac, 1, 13), RAT2(ac, 1, 14),
            x, r, w1);
}



VPUBLIC void Vnmresid27_1s(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *oC, double *cc, double *fc,
        double *oE, double *oN, double *uC,
        double *oNE, double *oNW,
        double *uE, double *uW, double *uN, double *uS,
        double *uNE, double *uNW, double *uSE, double *uSW,
        double *x, double *r, double *w1) {

    int i, j, k;
    int ipkey;
    double tmpO, tmpU, tmpD;

    MAT3( oC, *nx, *ny, *nz);
    MAT3( cc, *nx, *ny, *nz);
    MAT3( fc, *nx, *ny, *nz);
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
    MAT3(  x, *nx, *ny, *nz);
    MAT3(  r, *nx, *ny, *nz);
    MAT3( w1, *nx, *ny, *nz);

    // First get vector noNlinear term to avoid subroutine calls
    ipkey = VAT(ipc, 10);
    Vc_vec(cc, x, w1, nx, ny, nz, &ipkey);

    // The residual
    for (k=2; k<=*nz-1; k++) {
        for (j=2; j<=*ny-1; j++) {
            for (i=2; i<=*nx-1; i++) {

                tmpO =
                     +  VAT3( oN,   i,   j,   k) * VAT3(x,   i, j+1,   k)
                     +  VAT3( oN,   i, j-1,   k) * VAT3(x,   i, j-1,   k)
                     +  VAT3( oE,   i,   j,   k) * VAT3(x, i+1,   j,   k)
                     +  VAT3( oE, i-1,   j,   k) * VAT3(x, i-1,   j,   k)
                     +  VAT3(oNE,   i,   j,   k) * VAT3(x, i+1, j+1,   k)
                     +  VAT3(oNW,   i,   j,   k) * VAT3(x, i-1, j+1,   k)
                     +  VAT3(oNW, i+1, j-1,   k) * VAT3(x, i+1, j-1,   k)
                     +  VAT3(oNE, i-1, j-1,   k) * VAT3(x, i-1, j-1,   k);

                tmpU =
                     +  VAT3( uC,   i,   j,   k) * VAT3(x,   i,   j, k+1)
                     +  VAT3( uN,   i,   j,   k) * VAT3(x,   i, j+1, k+1)
                     +  VAT3( uS,   i,   j,   k) * VAT3(x,   i, j-1, k+1)
                     +  VAT3( uE,   i,   j,   k) * VAT3(x, i+1,   j, k+1)
                     +  VAT3( uW,   i,   j,   k) * VAT3(x, i-1,   j, k+1)
                     +  VAT3(uNE,   i,   j,   k) * VAT3(x, i+1, j+1, k+1)
                     +  VAT3(uNW,   i,   j,   k) * VAT3(x, i-1, j+1, k+1)
                     +  VAT3(uSE,   i,   j,   k) * VAT3(x, i+1, j-1, k+1)
                     +  VAT3(uSW,   i,   j,   k) * VAT3(x, i-1, j-1, k+1);

                tmpD =
                     +  VAT3( uC,   i,   j, k-1) * VAT3(x,   i,   j, k-1)
                     +  VAT3( uS,   i, j+1, k-1) * VAT3(x,   i, j+1, k-1)
                     +  VAT3( uN,   i, j-1, k-1) * VAT3(x,   i, j-1, k-1)
                     +  VAT3( uW, i+1,   j, k-1) * VAT3(x, i+1,   j, k-1)
                     +  VAT3( uE, i-1,   j, k-1) * VAT3(x, i-1,   j, k-1)
                     +  VAT3(uSW, i+1, j+1, k-1) * VAT3(x, i+1, j+1, k-1)
                     +  VAT3(uSE, i-1, j+1, k-1) * VAT3(x, i-1, j+1, k-1)
                     +  VAT3(uNW, i+1, j-1, k-1) * VAT3(x, i+1, j-1, k-1)
                     +  VAT3(uNE, i-1, j-1, k-1) * VAT3(x, i-1, j-1, k-1);

                VAT3(r, i, j, k) =
                        + tmpO + tmpU + tmpD
                        + VAT3(fc, i, j, k)
                        - VAT3(oC, i, j, k) * VAT3(x, i, j, k)
                        - VAT3(w1, i, j, k);
            }
        }
    }
}



VPUBLIC void Vrestrc(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        double *xin, double *xout, double *pc) {

    MAT2(pc, *nxc * *nyc * *nzc, 1 );

    Vrestrc2(nxf, nyf, nzf,
            nxc, nyc, nzc,
            xin, xout,
            RAT2(pc, 1, 1), RAT2(pc, 1, 2), RAT2(pc, 1, 3), RAT2(pc, 1, 4), RAT2(pc, 1, 5),
            RAT2(pc, 1, 6), RAT2(pc, 1, 7), RAT2(pc, 1, 8), RAT2(pc, 1, 9),
            RAT2(pc, 1,10), RAT2(pc, 1,11), RAT2(pc, 1,12), RAT2(pc, 1,13), RAT2(pc, 1,14),
            RAT2(pc, 1,15), RAT2(pc, 1,16), RAT2(pc, 1,17), RAT2(pc, 1,18),
            RAT2(pc, 1,19), RAT2(pc, 1,20), RAT2(pc, 1,21), RAT2(pc, 1,22), RAT2(pc, 1,23),
            RAT2(pc, 1,24), RAT2(pc, 1,25), RAT2(pc, 1,26), RAT2(pc, 1,27));
}



VEXTERNC void Vrestrc2(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        double  *xin, double *xout,
        double  *oPC, double  *oPN, double  *oPS, double  *oPE,  double *oPW,
        double *oPNE, double *oPNW, double *oPSE, double *oPSW,
        double  *uPC, double  *uPN, double  *uPS, double  *uPE,  double *uPW,
        double *uPNE, double *uPNW, double *uPSE, double *uPSW,
        double  *dPC, double  *dPN, double  *dPS, double  *dPE,  double *dPW,
        double *dPNE, double *dPNW, double *dPSE, double *dPSW) {

    int  i,  j,  k;
    int ii, jj, kk;
    int idimenshun = 3;

    double tmpO, tmpU, tmpD;
    double dimfac;

    MAT3(xin, *nxf, *nyf, *nzf);
    MAT3(xout, *nxc, *nyc, *nzc);

    MAT3(oPC, *nxc, *nyc, *nzc);
    MAT3(oPN, *nxc, *nyc, *nzc);
    MAT3(oPS, *nxc, *nyc, *nzc);
    MAT3(oPE, *nxc, *nyc, *nzc);
    MAT3(oPW, *nxc, *nyc, *nzc);

    MAT3(oPNE, *nxc, *nyc, *nzc);
    MAT3(oPNW, *nxc, *nyc, *nzc);
    MAT3(oPSE, *nxc, *nyc, *nzc);
    MAT3(oPSW, *nxc, *nyc, *nzc);

    MAT3(uPC, *nxc, *nyc, *nzc);
    MAT3(uPN, *nxc, *nyc, *nzc);
    MAT3(uPS, *nxc, *nyc, *nzc);
    MAT3(uPE, *nxc, *nyc, *nzc);
    MAT3(uPW, *nxc, *nyc, *nzc);

    MAT3(uPNE, *nxc, *nyc, *nzc);
    MAT3(uPNW, *nxc, *nyc, *nzc);
    MAT3(uPSE, *nxc, *nyc, *nzc);
    MAT3(uPSW, *nxc, *nyc, *nzc);

    MAT3(dPC, *nxc, *nyc, *nzc);
    MAT3(dPN, *nxc, *nyc, *nzc);
    MAT3(dPS, *nxc, *nyc, *nzc);
    MAT3(dPE, *nxc, *nyc, *nzc);
    MAT3(dPW, *nxc, *nyc, *nzc);

    MAT3(dPNE, *nxc, *nyc, *nzc);
    MAT3(dPNW, *nxc, *nyc, *nzc);
    MAT3(dPSE, *nxc, *nyc, *nzc);
    MAT3(dPSW, *nxc, *nyc, *nzc);

    // Verify correctness of the input boundary points
    VfboundPMG00(nxf, nyf, nzf, xin);

    dimfac = VPOW(2.0, idimenshun);

    // Handle the interior points as average of 5 finer grid pts ***
    #pragma omp parallel for private(k, kk, j, jj, i, ii, tmpO, tmpU, tmpD)
    for (k=2; k<=*nzc-1; k++) {
        kk = (k - 1) * 2 + 1;

        for (j=2; j<=*nyc-1; j++) {
            jj = (j - 1) * 2 + 1;

            for (i=2; i<=*nxc-1; i++) {
                ii = (i - 1) * 2 + 1;

                // Compute the restriction
                tmpO =
                     + VAT3( oPC, i, j, k) * VAT3(xin,   ii,   jj,   kk)
                     + VAT3( oPN, i, j, k) * VAT3(xin,   ii, jj+1,   kk)
                     + VAT3( oPS, i, j, k) * VAT3(xin,   ii, jj-1,   kk)
                     + VAT3( oPE, i, j, k) * VAT3(xin, ii+1,   jj,   kk)
                     + VAT3( oPW, i, j, k) * VAT3(xin, ii-1,   jj,   kk)
                     + VAT3(oPNE, i, j, k) * VAT3(xin, ii+1, jj+1,   kk)
                     + VAT3(oPNW, i, j, k) * VAT3(xin, ii-1, jj+1,   kk)
                     + VAT3(oPSE, i, j, k) * VAT3(xin, ii+1, jj-1,   kk)
                     + VAT3(oPSW, i, j, k) * VAT3(xin, ii-1, jj-1,   kk);

                tmpU =
                     + VAT3( uPC, i, j, k) * VAT3(xin,   ii,   jj, kk+1)
                     + VAT3( uPN, i, j, k) * VAT3(xin,   ii, jj+1, kk+1)
                     + VAT3( uPS, i, j, k) * VAT3(xin,   ii, jj-1, kk+1)
                     + VAT3( uPE, i, j, k) * VAT3(xin, ii+1,   jj, kk+1)
                     + VAT3( uPW, i, j, k) * VAT3(xin, ii-1,   jj, kk+1)
                     + VAT3(uPNE, i, j, k) * VAT3(xin, ii+1, jj+1, kk+1)
                     + VAT3(uPNW, i, j, k) * VAT3(xin, ii-1, jj+1, kk+1)
                     + VAT3(uPSE, i, j, k) * VAT3(xin, ii+1, jj-1, kk+1)
                     + VAT3(uPSW, i, j, k) * VAT3(xin, ii-1, jj-1, kk+1);

                tmpD =
                     + VAT3( dPC, i, j, k) * VAT3(xin,   ii,   jj, kk-1)
                     + VAT3( dPN, i, j, k) * VAT3(xin,   ii, jj+1, kk-1)
                     + VAT3( dPS, i, j, k) * VAT3(xin,   ii, jj-1, kk-1)
                     + VAT3( dPE, i, j, k) * VAT3(xin, ii+1,   jj, kk-1)
                     + VAT3( dPW, i, j, k) * VAT3(xin, ii-1,   jj, kk-1)
                     + VAT3(dPNE, i, j, k) * VAT3(xin, ii+1, jj+1, kk-1)
                     + VAT3(dPNW, i, j, k) * VAT3(xin, ii-1, jj+1, kk-1)
                     + VAT3(dPSE, i, j, k) * VAT3(xin, ii+1, jj-1, kk-1)
                     + VAT3(dPSW, i, j, k) * VAT3(xin, ii-1, jj-1, kk-1);

                VAT3(xout, i, j, k) = tmpO + tmpU + tmpD;
            }
        }
    }

    // Verify correctness of the output boundary points
    VfboundPMG00(nxc, nyc, nzc, xout);
}



VPUBLIC void VinterpPMG(int *nxc, int *nyc, int *nzc,
        int *nxf, int *nyf, int *nzf,
        double *xin, double *xout,
        double *pc) {

    MAT2(pc, *nxc * *nyc * *nzc, 1);

    VinterpPMG2(nxc, nyc, nzc,
            nxf, nyf, nzf,
            xin, xout,
            RAT2(pc, 1, 1), RAT2(pc, 1, 2), RAT2(pc, 1, 3), RAT2(pc, 1, 4), RAT2(pc, 1, 5),
            RAT2(pc, 1, 6), RAT2(pc, 1, 7), RAT2(pc, 1, 8), RAT2(pc, 1, 9),
            RAT2(pc, 1,10), RAT2(pc, 1,11), RAT2(pc, 1,12), RAT2(pc, 1,13), RAT2(pc, 1,14),
            RAT2(pc, 1,15), RAT2(pc, 1,16), RAT2(pc, 1,17), RAT2(pc, 1,18),
            RAT2(pc, 1,19), RAT2(pc, 1,20), RAT2(pc, 1,21), RAT2(pc, 1,22), RAT2(pc, 1,23),
            RAT2(pc, 1,24), RAT2(pc, 1,25), RAT2(pc, 1,26), RAT2(pc, 1,27));
}



VPUBLIC void VinterpPMG2(int *nxc, int *nyc, int *nzc,
        int *nxf, int *nyf, int *nzf,
        double *xin, double *xout,
        double  *oPC, double  *oPN, double  *oPS, double  *oPE, double  *oPW,
        double *oPNE, double *oPNW, double *oPSE, double *oPSW,
        double  *uPC, double  *uPN, double  *uPS, double  *uPE, double  *uPW,
        double *uPNE, double *uPNW, double *uPSE, double *uPSW,
        double  *dPC, double  *dPN, double  *dPS, double  *dPE, double  *dPW,
        double *dPNE, double *dPNW, double *dPSE, double *dPSW) {

    int  i,  j,  k;
    int ii, jj, kk;

    MAT3( xin, *nxc, *nyc, *nzc);
    MAT3(xout, *nxf, *nyf, *nzf);

    MAT3( oPC, *nxc, *nyc, *nzc);
    MAT3( oPN, *nxc, *nyc, *nzc);
    MAT3( oPS, *nxc, *nyc, *nzc);
    MAT3( oPE, *nxc, *nyc, *nzc);
    MAT3( oPW, *nxc, *nyc, *nzc);

    MAT3(oPNE, *nxc, *nyc, *nzc);
    MAT3(oPNW, *nxc, *nyc, *nzc);
    MAT3(oPSE, *nxc, *nyc, *nzc);
    MAT3(oPSW, *nxc, *nyc, *nzc);

    MAT3( uPC, *nxc, *nyc, *nzc);
    MAT3( uPN, *nxc, *nyc, *nzc);
    MAT3( uPS, *nxc, *nyc, *nzc);
    MAT3( uPE, *nxc, *nyc, *nzc);
    MAT3( uPW, *nxc, *nyc, *nzc);

    MAT3(uPNE, *nxc, *nyc, *nzc);
    MAT3(uPNW, *nxc, *nyc, *nzc);
    MAT3(uPSE, *nxc, *nyc, *nzc);
    MAT3(uPSW, *nxc, *nyc, *nzc);

    MAT3( dPC, *nxc, *nyc, *nzc);
    MAT3( dPN, *nxc, *nyc, *nzc);
    MAT3( dPS, *nxc, *nyc, *nzc);
    MAT3( dPE, *nxc, *nyc, *nzc);
    MAT3( dPW, *nxc, *nyc, *nzc);

    MAT3(dPNE, *nxc, *nyc, *nzc);
    MAT3(dPNW, *nxc, *nyc, *nzc);
    MAT3(dPSE, *nxc, *nyc, *nzc);
    MAT3(dPSW, *nxc, *nyc, *nzc);

    /* *********************************************************************
     * Setup
     * *********************************************************************/

    // Verify correctness of the input boundary points ***
    VfboundPMG00(nxc, nyc, nzc, xin);

    // Do it
    for (k=1; k<=*nzf-2; k+=2) {
        kk = (k - 1) / 2 + 1;

        for (j=1; j<=*nyf-2; j+=2) {
            jj = (j - 1) / 2 + 1;

            for (i=1; i<=*nxf-2; i+=2) {
                ii = (i - 1) / 2 + 1;

                /* ******************************************************** *
                 * Type 1 -- Fine grid points common to a coarse grid point *
                 * ******************************************************** */

                // Copy coinciding points from coarse grid to fine grid
                VAT3(xout, i, j, k) = VAT3(xin, ii, jj, kk);

                /* ******************************************************** *
                 * type 2 -- fine grid points common to a coarse grid plane *
                 * ******************************************************** */

                // Fine grid pts common only to y-z planes on coarse grid
                // (intermediate pts between 2 grid points on x-row)
                VAT3(xout, i+1, j, k) = VAT3(oPE,   ii, jj, kk) * VAT3(xin,   ii, jj, kk)
                                + VAT3(oPW, ii+1, jj, kk) * VAT3(xin, ii+1, jj, kk);

                // Fine grid pts common only to x-z planes on coarse grid
                // (intermediate pts between 2 grid points on a y-row)
                VAT3(xout, i, j+1, k) = VAT3(oPN, ii,   jj, kk) * VAT3(xin, ii,   jj, kk)
                                + VAT3(oPS, ii, jj+1, kk) * VAT3(xin, ii, jj+1, kk);

                // Fine grid pts common only to x-y planes on coarse grid
                // (intermediate pts between 2 grid points on a z-row)
                VAT3(xout, i, j, k+1) = VAT3(uPC, ii, jj,   kk) * VAT3(xin, ii, jj,   kk)
                                + VAT3(dPC, ii, jj, kk+1) * VAT3(xin, ii, jj, kk+1);

                /* ******************************************************* *
                 * type 3 -- fine grid points common to a coarse grid line *
                 * ******************************************************* */

                // Fine grid pts common only to z planes on coarse grid
                // (intermediate pts between 4 grid pts on the xy-plane

                VAT3(xout, i+1, j+1, k) = VAT3(oPNE,   ii,   jj, kk) * VAT3(xin,   ii,   jj, kk)
                                  + VAT3(oPNW, ii+1,   jj, kk) * VAT3(xin, ii+1,   jj, kk)
                                  + VAT3(oPSE,   ii, jj+1, kk) * VAT3(xin,   ii, jj+1, kk)
                                  + VAT3(oPSW, ii+1, jj+1, kk) * VAT3(xin, ii+1, jj+1, kk);

                // Fine grid pts common only to y planes on coarse grid
                // (intermediate pts between 4 grid pts on the xz-plane
                VAT3(xout, i+1, j, k+1) = VAT3(uPE,   ii, jj,   kk) * VAT3(xin,   ii, jj,   kk)
                                  + VAT3(uPW, ii+1, jj,   kk) * VAT3(xin, ii+1, jj,   kk)
                                  + VAT3(dPE,   ii, jj, kk+1) * VAT3(xin,   ii, jj, kk+1)
                                  + VAT3(dPW, ii+1, jj, kk+1) * VAT3(xin, ii+1, jj, kk+1);

                // Fine grid pts common only to x planes on coarse grid
                // (intermediate pts between 4 grid pts on the yz-plane***
                VAT3(xout, i, j+1, k+1) = VAT3(uPN, ii,   jj,  kk) * VAT3(xin, ii,   jj,   kk)
                                  + VAT3(uPS, ii, jj+1,  kk) * VAT3(xin, ii, jj+1,   kk)
                                  + VAT3(dPN, ii,   jj,kk+1) * VAT3(xin, ii,   jj, kk+1)
                                  + VAT3(dPS, ii, jj+1,kk+1) * VAT3(xin, ii, jj+1, kk+1);

                /* **************************************** *
                 * type 4 -- fine grid points not common to *
                 *           coarse grid pts/lines/planes   *
                 * **************************************** */

                // Completely interior points
                VAT3(xout, i+1,j+1,k+1) =
                     + VAT3(uPNE,   ii,   jj,   kk) * VAT3(xin,   ii,   jj,   kk)
                     + VAT3(uPNW, ii+1,   jj,   kk) * VAT3(xin, ii+1,   jj,   kk)
                     + VAT3(uPSE,   ii, jj+1,   kk) * VAT3(xin,   ii, jj+1,   kk)
                     + VAT3(uPSW, ii+1, jj+1,   kk) * VAT3(xin, ii+1, jj+1,   kk)
                     + VAT3(dPNE,   ii,   jj, kk+1) * VAT3(xin,   ii,   jj, kk+1)
                     + VAT3(dPNW, ii+1,   jj, kk+1) * VAT3(xin, ii+1,   jj, kk+1)
                     + VAT3(dPSE,   ii, jj+1, kk+1) * VAT3(xin,   ii, jj+1, kk+1)
                     + VAT3(dPSW, ii+1, jj+1, kk+1) * VAT3(xin, ii+1, jj+1, kk+1);
            }
        }
    }

    // Verify correctness of the output boundary points ***
    VfboundPMG00(nxf, nyf, nzf, xout);
}



VPUBLIC void Vextrac(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        double *xin, double *xout) {

    int  i,  j,  k;
    int ii, jj, kk;

    MAT3( xin, *nxf, *nyf, *nzf);
    MAT3(xout, *nxc, *nyc, *nzc);

    // Verify correctness of the input boundary points
    VfboundPMG00(nxf, nyf, nzf, xin);

    // Do it
    for (k=2; k<=*nzc-1; k++) {
       kk = (k - 1) * 2 + 1;

       for (j=2; j<=*nyc-1; j++) {
          jj = (j - 1) * 2 + 1;

          for (i=2; i<=*nxc-1; i++) {
             ii = (i - 1) * 2 + 1;

             // Compute the restriction
             VAT3(xout, i, j, k) = VAT3(xin, ii, jj, kk);
          }
       }
    }

    // Verify correctness of the output boundary points
    VfboundPMG00(nxc, nyc, nzc, xout);
}
