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

#include "buildBd.h"

VPUBLIC void Vbuildband(int *key, int *nx, int *ny, int *nz,
        int *ipc, double *rpc, double *ac,
        int *ipcB, double *rpcB, double *acB) {

    int numdia;
    int n, m;
    int lda, info;

    MAT2(ac, *nx * *ny * *nz, 1);

    // Do in one step
    numdia = VAT(ipc, 11);
    if (numdia == 7) {

       n   = (*nx - 2) * (*ny - 2) * (*nz - 2);
       m   = (*nx - 2) * (*ny - 2);
       lda = m + 1;

       Vbuildband1_7
              (nx, ny, nz,
               ipc, rpc,
               RAT2(ac, 1, 1), RAT2(ac, 1, 2), RAT2(ac, 1, 3), RAT2(ac, 1, 4),
               ipcB, rpcB, acB,
               &n, &m, &lda);

    } else if (numdia == 27) {

        n   = (*nx - 2) * (*ny - 2) * (*nz - 2);
        m   = (*nx - 2) * (*ny - 2) + (*nx - 2) + 1;
        lda = m + 1;

        Vbuildband1_27
               (nx, ny, nz,
                ipc, rpc,
                RAT2(ac, 1,  1), RAT2(ac, 1,  2), RAT2(ac, 1,  3), RAT2(ac, 1,  4),
                RAT2(ac, 1,  5), RAT2(ac, 1,  6),
                RAT2(ac, 1,  7), RAT2(ac, 1,  8), RAT2(ac, 1,  9), RAT2(ac, 1, 10),
                RAT2(ac, 1, 11), RAT2(ac, 1, 12), RAT2(ac, 1, 13), RAT2(ac, 1, 14),
                ipcB, rpcB, acB,
                &n, &m, &lda);
    } else {
        Vnm_print(2, "Vbuildband: invalid stencil type given...");
    }

    // Factor the system
    *key  = 0;
    info = 0;

    Vdpbfa(acB, &lda, &n, &m, &info);
    VAT(ipcB, 4) = 1;

    if (info != 0) {

        Vnm_print(2, "Vbuildband: dpbfa problem: %d\n", info);
        Vnm_print(2, "Vbuildband: leading principle minor not PD...\n");

        *key = 1;
    }
}



VPUBLIC void Vbuildband1_7(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *oC, double *oE, double *oN, double *uC,
        int *ipcB, double *rpcB, double *acB,
        int *n, int *m, int *lda) {

    int  i,  j,  k;
    int ii, jj, kk;

    MAT2(acB, *lda, *ny-1);

    MAT3(oC, *nx, *ny, *nz);
    MAT3(oE, *nx, *ny, *nz);
    MAT3(oN, *nx, *ny, *nz);
    MAT3(uC, *nx, *ny, *nz);

    WARN_UNTESTED;

    // Do it
    VAT(ipcB, 1) = *n;
    VAT(ipcB, 2) = *m;
    VAT(ipcB, 3) = *lda;
    VAT(ipcB, 4) = 0;

    jj = 0;

    //fprintf(data, "%s\n", PRINT_FUNC);

    for (k=2; k<=*nz-1; k++) {

        for (j=2; j<=*ny-1; j++) {

            for (i=2; i<=*nx-1; i++) {
                jj++;

                // Diagonal term
                ii = jj;
                kk = ii - jj + *m + 1;

                VAT2(acB, kk, jj) = VAT3(oC, i, j, k);

                // East neighbor
                ii = jj - 1;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(oE, i-1, j, k);

                // North neighbor
                ii = jj - (*nx - 2);
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(oN, i, j-1, k);

                // Up neighbor ***
                ii = jj - (*nx - 2) * (*ny - 2);
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(uC, i, j, k-1);

                //fprintf(data, "%19.12E\n", VAT2(acB, kk, jj));
            }
        }
    }
}



VPUBLIC void Vbuildband1_27(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double  *oC, double  *oE, double  *oN, double  *uC,
        double *oNE, double *oNW,
        double  *uE, double  *uW, double  *uN, double  *uS,
        double *uNE, double *uNW, double *uSE, double *uSW,
        int *ipcB, double *rpcB, double *acB,
        int *n, int *m, int *lda) {

    int  i,  j,  k;
    int ii, jj, kk;

    MAT2(acB, *lda, *ny-1);

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

    // Do it
    VAT(ipcB, 1) = *n;
    VAT(ipcB, 2) = *m;
    VAT(ipcB, 3) = *lda;
    VAT(ipcB, 4) = 0;

    jj = 0;

    //fprintf(data, "%s\n", PRINT_FUNC);

    for (k=2; k<=*nz-1; k++) {

        for (j=2; j<=*ny-1; j++) {

            for (i=2; i<=*nx-1; i++) {
                jj++;

                // Diagonal term
                ii = jj;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = VAT3(oC, i, j, k);

                // East neighbor
                ii = jj - 1;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(oE, i-1, j, k);

                // North neighbor
                ii = jj - (*nx - 2);
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(oN, i, j-1, k);

                // North-east neighbor
                ii = jj - (*nx - 2) + 1;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(oNE, i, j-1, k);

                // North-west neighbor
                ii = jj - (*nx - 2) - 1;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(oNW, i, j-1, k);

                // Up neighbor
                ii = jj - (*nx - 2) * (*ny - 2);
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(uC, i, j, k-1);

                // Up-east neighbor
                ii = jj - (*nx - 2) * (*ny - 2) +1;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(uE, i, j, k-1);

                // Up-west neighbor
                ii = jj - (*nx - 2) * (*ny - 2) - 1;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(uW, i, j, k-1);

                // Up-north neighbor
                ii = jj - (*nx - 2) * (*ny - 2) + (*nx - 2);
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(uN, i, j, k-1);

                // Up-south neighbor
                ii = jj - (*nx - 2) * (*ny - 2) - (*nx - 2);
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(uS, i, j, k-1);

                // Up-north-east neighbor
                ii = jj - (*nx - 2) * (*ny - 2) + (*nx - 2) + 1;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(uNE, i, j, k-1);

                // Up-north-west neighbor
                ii = jj - (*nx - 2) * (*ny - 2) + (*nx - 2) - 1;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(uNW, i, j, k-1);

                // Up-south-east neighbor
                ii = jj - (*nx - 2) * (*ny - 2) - (*nx - 2) + 1;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(uSE, i, j, k-1);

                // Up-south-west neighbor
                ii = jj - (*nx - 2) * (*ny - 2) - (*nx - 2) - 1;
                kk = ii - jj + *m + 1;
                VAT2(acB, kk, jj) = - VAT3(uSW, i, j, k-1);

                //fprintf(data, "%19.12E\n", VAT2(acB, kk, jj));
            }
        }
    }
}
