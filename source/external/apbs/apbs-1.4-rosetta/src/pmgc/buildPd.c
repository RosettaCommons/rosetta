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

#include "buildPd.h"

VPUBLIC void VbuildP(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        int *mgprol,
        int *ipc, double *rpc,
        double *pc, double *ac,
        double *xf, double *yf, double *zf) {

    int numdia;

    MAT2(pc, *nxc * *nyc * *nzc, 1);
    MAT2(ac, *nxf * *nyf * *nzf, 1);

    if (*mgprol == 0) {

            VbuildP_trilin(nxf, nyf, nzf,
                    nxc, nyc, nzc,
                    RAT2(pc, 1, 1),
                    xf, yf, zf);

    } else if (*mgprol == 1) {

        numdia = VAT(ipc, 11);

        if (numdia == 7) {
            VbuildP_op7(nxf, nyf, nzf,
                    nxc, nyc, nzc,
                    ipc, rpc,
                    RAT2(ac, 1, 1), RAT2(pc, 1, 1));
        } else if (numdia == 27) {
            VbuildP_op27(nxf, nyf, nzf,
                    nxc, nyc, nzc,
                    ipc, rpc,
                    RAT2(ac, 1, 1), RAT2(pc, 1, 1));
        } else {
            Vnm_print(2,"BUILDP: invalid stencil type given: %d\n", numdia);
        }
    }
}

VPUBLIC void VbuildP_trilin(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        double *pc,
        double *xf, double *yf, double *zf) {

    MAT2(pc, *nxc * *nyc * *nzc, 1);

    VbuildPb_trilin(nxf, nyf, nzf,
            nxc, nyc, nzc,
            RAT2(pc, 1,  1),RAT2(pc, 1,  2),RAT2(pc, 1,  3),RAT2(pc, 1,  4),RAT2(pc, 1,  5),
            RAT2(pc, 1,  6),RAT2(pc, 1,  7),RAT2(pc, 1,  8),RAT2(pc, 1,  9),
            RAT2(pc, 1, 10),RAT2(pc, 1, 11),RAT2(pc, 1, 12),RAT2(pc, 1, 13),RAT2(pc, 1, 14),
            RAT2(pc, 1, 15),RAT2(pc, 1, 16),RAT2(pc, 1, 17),RAT2(pc, 1, 18),
            RAT2(pc, 1, 19),RAT2(pc, 1, 20),RAT2(pc, 1, 21),RAT2(pc, 1, 22),RAT2(pc, 1, 23),
            RAT2(pc, 1, 24),RAT2(pc, 1, 25),RAT2(pc, 1, 26),RAT2(pc, 1, 27),
            xf, yf, zf);
}

VEXTERNC void VbuildPb_trilin(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        double *oPC, double *oPN, double *oPS, double *oPE, double *oPW,
        double *oPNE, double *oPNW, double *oPSE, double *oPSW,
        double *uPC, double *uPN, double *uPS, double *uPE, double *uPW,
        double *uPNE, double *uPNW, double *uPSE, double *uPSW,
        double *dPC, double *dPN, double *dPS, double *dPE, double *dPW,
        double *dPNE, double *dPNW, double *dPSE, double *dPSW,
        double *xf, double *yf, double *zf) {

    int i, j, k;

    /// @todo Document this function
    double won     = 1.0;
    double half    = 1.0 / 2.0;
    double quarter = 1.0 / 4.0;
    double eighth  = 1.0 / 8.0;

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

for(k=2; k<=*nzc-1; k++) {
    for(j=2; j<=*nyc-1; j++) {
            for(i=2; i<=*nxc-1; i++) {

                VAT3(oPC, i,j,k)  = won;

                VAT3(oPN, i,j,k)  = half;
                VAT3(oPS, i,j,k)  = half;
                VAT3(oPE, i,j,k)  = half;
                VAT3(oPW, i,j,k)  = half;
                VAT3(uPC, i,j,k)  = half;
                VAT3(dPC, i,j,k)  = half;

                VAT3(oPNE, i,j,k) = quarter;
                VAT3(oPNW, i,j,k) = quarter;
                VAT3(oPSE, i,j,k) = quarter;
                VAT3(oPSW, i,j,k) = quarter;
                VAT3(dPN, i,j,k)  = quarter;
                VAT3(dPS, i,j,k)  = quarter;
                VAT3(dPE, i,j,k)  = quarter;
                VAT3(dPW, i,j,k)  = quarter;
                VAT3(uPN, i,j,k)  = quarter;
                VAT3(uPS, i,j,k)  = quarter;
                VAT3(uPE, i,j,k)  = quarter;
                VAT3(uPW, i,j,k)  = quarter;

                VAT3(dPNE, i,j,k) = eighth;
                VAT3(dPNW, i,j,k) = eighth;
                VAT3(dPSE, i,j,k) = eighth;
                VAT3(dPSW, i,j,k) = eighth;
                VAT3(uPNE, i,j,k) = eighth;
                VAT3(uPNW, i,j,k) = eighth;
                VAT3(uPSE, i,j,k) = eighth;
                VAT3(uPSW, i,j,k) = eighth;
            }
        }
    }
}



VPUBLIC void VbuildP_op7(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        int *ipc, double *rpc,
        double *ac, double *pc) {

    MAT2(ac, *nxf * *nyf * *nzf, 1);
    MAT2(pc, *nxc * *nyc * *nzc, 1);

    WARN_UNTESTED;

    VbuildPb_op7(nxf, nyf, nzf,
            nxc, nyc, nzc,
            ipc, rpc,
            RAT2(ac, 1,  1), RAT2(ac, 1,  2), RAT2(ac, 1,  3),
            RAT2(ac, 1,  4),
            RAT2(pc, 1,  1), RAT2(pc, 1,  2), RAT2(pc, 1,  3), RAT2(pc, 1,  4), RAT2(pc, 1,  5),
            RAT2(pc, 1,  6), RAT2(pc, 1,  7), RAT2(pc, 1,  8), RAT2(pc, 1,  9),
            RAT2(pc, 1, 10), RAT2(pc, 1, 11), RAT2(pc, 1, 12), RAT2(pc, 1, 13), RAT2(pc, 1, 14),
            RAT2(pc, 1, 15), RAT2(pc, 1, 16), RAT2(pc, 1, 17), RAT2(pc, 1, 18),
            RAT2(pc, 1, 19), RAT2(pc, 1, 20), RAT2(pc, 1, 21), RAT2(pc, 1, 22), RAT2(pc, 1, 23),
            RAT2(pc, 1, 24), RAT2(pc, 1, 25), RAT2(pc, 1, 26), RAT2(pc, 1, 27));
}



VPUBLIC void VbuildPb_op7(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        int *ipc, double *rpc,
        double   *oC, double   *oE, double   *oN,
        double   *uC,
        double  *oPC, double  *oPN, double  *oPS, double  *oPE, double  *oPW,
        double *oPNE, double *oPNW, double *oPSE, double *oPSW,
        double  *uPC, double  *uPN, double  *uPS, double  *uPE, double  *uPW,
        double *uPNE, double *uPNW, double *uPSE, double *uPSW,
        double  *dPC, double  *dPN, double  *dPS, double  *dPE, double  *dPW,
        double *dPNE, double *dPNW, double *dPSE, double *dPSW) {

    int  i,  j,  k;
    int ii, jj, kk;
    int  im1,  ip1;
    int  im2,  ip2;
    int  jm1,  jp1;
    int  jm2,  jp2;
    int  km1,  kp1;
    int  km2,  kp2;
    int iim1, iip1;
    int jjm1, jjp1;
    int kkm1, kkp1;

    double won, half, quarter, eighth;

    MAT3(  oC, *nxf, *nyf, *nzf);
    MAT3(  oE, *nxf, *nyf, *nzf);
    MAT3(  oN, *nxf, *nyf, *nzf);
    MAT3(  uC, *nxf, *nyf, *nzf);
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

    WARN_UNTESTED;

    // interpolation stencil ***
    won     =  1.0;
    half    =  1.0 /  2.0;
    quarter =  1.0 /  4.0;
    eighth  =  1.0 /  8.0;

    //fprintf(data, "%s\n", PRINT_FUNC);

    for (kk = 2; kk < *nzc - 1; kk++) {
        k = 2 * kk - 1;

        for (jj = 2; jj < *nyc - 1; jj++) {
            j = 2 * jj - 1;

            for (ii = 2; ii < *nxc - 1; ii++) {
                i = 2 * ii - 1;

                // index computations ***
                im1  =  i - 1;
                ip1  =  i + 1;
                im2  =  i - 2;
                ip2  =  i + 2;
                jm1  =  j - 1;
                jp1  =  j + 1;
                jm2  =  j - 2;
                jp2  =  j + 2;
                km1  =  k - 1;
                kp1  =  k + 1;
                km2  =  k - 2;
                kp2  =  k + 2;
                iim1 = ii - 1;
                iip1 = ii + 1;
                jjm1 = jj - 1;
                jjp1 = jj + 1;
                kkm1 = kk - 1;
                kkp1 = kk + 1;

                // *************************************************************
                // *** > oPC;
                // *************************************************************

                VAT3( oPC, ii, jj, kk) = won;

                //fprintf(data, "%19.12E\n", VAT3(oPC, ii, jj, kk));

                // *************************************************************
                // *** > oPN;
                // *************************************************************

                VAT3( oPN, ii, jj, kk) =
                        VAT3(  oN, i, j, k) / ( VAT3(  oC,   i, jp1,   k)
                                               - VAT3(  oE, im1, jp1,   k)
                                               - VAT3(  oE,   i, jp1,   k)
                                               - VAT3(  uC,   i, jp1, km1)
                                               - VAT3(  uC,   i, jp1,   k));

                //fprintf(data, "%19.12E\n", VAT3(oPN, ii, jj, kk));

                // *************************************************************
                // *** > oPS;
                // *************************************************************

                VAT3( oPS, ii, jj, kk) =
                        VAT3(  oN, i, jm1, k) / ( VAT3(  oC,   i, jm1,   k)
                                                 - VAT3(  oE, im1, jm1,   k)
                                                 - VAT3(  oE,   i, jm1,   k)
                                                 - VAT3(  uC,   i, jm1, km1)
                                                 - VAT3(  uC,   i, jm1,   k));

                //fprintf(data, "%19.12E\n", VAT3(oPS, ii, jj, kk));

                // *************************************************************
                // *** > oPE;
                // *************************************************************

                VAT3( oPE, ii, jj, kk) =
                        VAT3(  oE, i, j, k) / ( VAT3(  oC, ip1,   j,   k)
                                               - VAT3(  uC, ip1,   j, km1)
                                               - VAT3(  uC, ip1,   j,   k)
                                               - VAT3(  oN, ip1,   j,   k)
                                               - VAT3(  oN, ip1, jm1,   k));

                //fprintf(data, "%19.12E\n", VAT3(oPE, ii, jj, kk));

                // *************************************************************
                // *** > oPW;
                // *************************************************************

                VAT3( oPW, ii, jj, kk) =
                        VAT3(  oE, im1, j, k) / ( VAT3(  oC, im1,   j,   k)
                                                 - VAT3(  uC, im1,   j, km1)
                                                 - VAT3(  uC, im1,   j,   k)
                                                 - VAT3(  oN, im1,   j,   k)
                                                 - VAT3(  oN, im1, jm1,   k));

                //fprintf(data, "%19.12E\n", VAT3(oPW, ii, jj, kk));

                // *************************************************************
                // *** > oPNE;
                // *************************************************************

                VAT3(oPNE, ii, jj, kk) =
                        (
                            VAT3(  oN, ip1,   j,   k) * VAT3( oPE, ii, jj, kk)
                          + VAT3(  oE,   i, jp1,   k) * VAT3( oPN, ii, jj, kk)
                        ) / (
                            VAT3(  oC, ip1, jp1,   k)
                          - VAT3(  uC, ip1, jp1, km1)
                          - VAT3(  uC, ip1, jp1,   k)
                        );

                //fprintf(data, "%19.12E\n", VAT3(oPNE, ii, jj, kk));

                // *************************************************************
                // *** > oPNW;
                // *************************************************************

                VAT3(oPNW, ii, jj, kk) =
                        (
                            VAT3(  oN, im1,   j,   k) * VAT3( oPW, ii, jj, kk)
                          + VAT3(  oE, im1, jp1,   k) * VAT3( oPN, ii, jj, kk)
                        ) / (
                            VAT3(  oC, im1, jp1,   k)
                          - VAT3(  uC, im1, jp1, km1)
                          - VAT3(  uC, im1, jp1,   k)
                        );

                //fprintf(data, "%19.12E\n", VAT3(oPNW, ii, jj, kk));

                // *************************************************************
                // *** > oPSE;
                // *************************************************************

                VAT3(oPSE, ii, jj, kk) =
                    (
                        VAT3(  oN, ip1, jm1,   k) * VAT3( oPE, ii, jj, kk)
                      + VAT3(  oE,   i, jm1,   k) * VAT3( oPS, ii, jj, kk)
                    ) / (
                        VAT3(  oC, ip1, jm1,   k)
                      - VAT3(  uC, ip1, jm1, km1)
                      - VAT3(  uC, ip1, jm1,   k)
                    );

                //fprintf(data, "%19.12E\n", VAT3(oPSE, ii, jj, kk));

                // *************************************************************
                // *** > oPSW;
                // *************************************************************

                VAT3(oPSW, ii, jj, kk) =
                    (
                        VAT3(  oN, im1, jm1,   k) * VAT3( oPW, ii, jj, kk)
                      + VAT3(  oE, im1, jm1,   k) * VAT3( oPS, ii, jj, kk)
                    ) / (
                        VAT3(  oC, im1, jm1,   k)
                      - VAT3(  uC, im1, jm1, km1)
                      - VAT3(  uC, im1, jm1,   k)
                    );

                //fprintf(data, "%19.12E\n", VAT3(oPSW, ii, jj, kk));

                // *************************************************************
                // *** > dPC;
                // *************************************************************

                VAT3( dPC, ii, jj, kk) =
                    VAT3(  uC, i, j, km1)
                    / (
                          VAT3(  oC,   i,   j, km1)
                        - VAT3(  oN,   i,   j, km1)
                        - VAT3(  oN,   i, jm1, km1)
                        - VAT3(  oE, im1,   j, km1)
                        - VAT3(  oE,   i,   j, km1)
                    );

                //fprintf(data, "%19.12E\n", VAT3(dPC, ii, jj, kk));

                // *************************************************************
                // *** > dPN;
                // *************************************************************

                VAT3( dPN, ii, jj, kk) =
                    (
                        VAT3(  oN,   i,   j, km1) * VAT3( dPC, ii, jj, kk)
                      + VAT3(  uC,   i, jp1, km1) * VAT3( oPN, ii, jj, kk)
                    ) / (
                        VAT3(  oC,   i, jp1, km1)
                      - VAT3(  oE, im1, jp1, km1)
                      - VAT3(  oE,   i, jp1, km1)
                    );

                //fprintf(data, "%19.12E\n", VAT3(dPN, ii, jj, kk));

                // *************************************************************
                // *** > dPS;
                // *************************************************************

                VAT3( dPS, ii, jj, kk) =
                    (
                        VAT3(  oN,   i, jm1, km1) * VAT3( dPC, ii, jj, kk)
                      + VAT3(  uC,   i, jm1, km1) * VAT3( oPS, ii, jj, kk)
                    ) / (
                        VAT3(  oC,   i, jm1, km1)
                      - VAT3(  oE, im1, jm1, km1)
                      - VAT3(  oE,   i, jm1, km1)
                    );

                //fprintf(data, "%19.12E\n", VAT3(dPS, ii, jj, kk));

                // *************************************************************
                // *** > dPE;
                // *************************************************************

                VAT3( dPE, ii, jj, kk) =
                    (
                        VAT3(  uC, ip1,   j, km1) * VAT3( oPE, ii, jj, kk)
                      + VAT3(  oE,   i,   j, km1) * VAT3( dPC, ii, jj, kk)
                    ) / (
                        VAT3(  oC, ip1,   j, km1)
                      - VAT3(  oN, ip1,   j, km1)
                      - VAT3(  oN, ip1, jm1, km1)
                    );

                //fprintf(data, "%19.12E\n", VAT3(dPE, ii, jj, kk));

                // *************************************************************
                // *** > dPW;
                // *************************************************************

                VAT3( dPW, ii, jj, kk) =
                    (
                        VAT3(  uC, im1,   j, km1) * VAT3( oPW, ii, jj, kk)
                      + VAT3(  oE, im1,   j, km1) * VAT3( dPC, ii, jj, kk)
                    ) / (
                        VAT3(  oC, im1,   j, km1)
                      - VAT3(  oN, im1,   j, km1)
                      - VAT3(  oN, im1, jm1, km1)
                    );

                //fprintf(data, "%19.12E\n", VAT3(dPW, ii, jj, kk));

                // *************************************************************
                // *** > dPNE;
                // *************************************************************

                VAT3(dPNE, ii, jj, kk) =
                    (
                        VAT3(  uC, ip1, jp1, km1) * VAT3(oPNE, ii, jj, kk)
                      + VAT3(  oE,   i, jp1, km1) * VAT3( dPN, ii, jj, kk)
                      + VAT3(  oN, ip1,   j, km1) * VAT3( dPE, ii, jj, kk)
                    ) / VAT3(  oC, ip1, jp1, km1);

                //fprintf(data, "%19.12E\n", VAT3(dPNE, ii, jj, kk));

                // *************************************************************
                // *** > dPNW;
                // *************************************************************

                VAT3(dPNW, ii, jj, kk) =
                    (
                        VAT3(  uC, im1, jp1, km1) * VAT3(oPNW, ii, jj, kk)
                      + VAT3(  oE, im1, jp1, km1) * VAT3( dPN, ii, jj, kk)
                      + VAT3(  oN, im1,   j, km1) * VAT3( dPW, ii, jj, kk)
                    ) / VAT3(  oC, im1, jp1, km1);

                //fprintf(data, "%19.12E\n", VAT3(dPNW, ii, jj, kk));

                // *************************************************************
                // *** > dPSE;
                // *************************************************************

                VAT3(dPSE, ii, jj, kk) =
                    (
                        VAT3(  uC, ip1, jm1, km1) * VAT3(oPSE, ii, jj, kk)
                      + VAT3(  oE,   i, jm1, km1) * VAT3( dPS, ii, jj, kk)
                      + VAT3(  oN, ip1, jm1, km1) * VAT3( dPE, ii, jj, kk)
                    ) / VAT3(  oC, ip1, jm1, km1);

                //fprintf(data, "%19.12E\n", VAT3(dPSE, ii, jj, kk));

                // *************************************************************
                // *** > dPSW;
                // *************************************************************

                VAT3(dPSW, ii, jj, kk) =
                    (
                        VAT3(  uC, im1, jm1, km1) * VAT3(oPSW, ii, jj, kk)
                      + VAT3(  oE, im1, jm1, km1) * VAT3( dPS, ii, jj, kk)
                      + VAT3(  oN, im1, jm1, km1) * VAT3( dPW, ii, jj, kk)
                    ) / VAT3(  oC, im1, jm1, km1);

                //fprintf(data, "%19.12E\n", VAT3(dPSW, ii, jj, kk));

                // *************************************************************
                // *** > uPC;
                // *************************************************************

                VAT3( uPC, ii, jj, kk) =
                    VAT3(  uC, i, j, k)
                    / ( VAT3(  oC,   i,   j, kp1)
                      - VAT3(  oN,   i,   j, kp1)
                      - VAT3(  oN,   i, jm1, kp1)
                      - VAT3(  oE, im1,   j, kp1)
                      - VAT3(  oE,   i,   j, kp1)
                    );

                //fprintf(data, "%19.12E\n", VAT3(uPC, ii, jj, kk));

                // *************************************************************
                // *** > uPN;
                // *************************************************************

                VAT3( uPN, ii, jj, kk) =
                    (
                        VAT3(  oN,   i,   j, kp1) * VAT3( uPC, ii, jj, kk)
                     +  VAT3(  uC,   i, jp1,   k) * VAT3( oPN, ii, jj, kk)
                    ) / (
                        VAT3(  oC,   i, jp1, kp1)
                      - VAT3(  oE, im1, jp1, kp1)
                      - VAT3(  oE,   i, jp1, kp1)
                    );

                //fprintf(data, "%19.12E\n", VAT3(uPN, ii, jj, kk));

                // *************************************************************
                // *** > uPS;
                // *************************************************************

                VAT3( uPS, ii, jj, kk) =
                    (
                        VAT3(  oN,   i, jm1, kp1) * VAT3( uPC, ii, jj, kk)
                     +  VAT3(  uC,   i, jm1,   k) * VAT3( oPS, ii, jj, kk)
                    ) / (
                        VAT3(  oC,   i, jm1, kp1)
                      - VAT3(  oE, im1, jm1, kp1)
                      - VAT3(  oE,   i, jm1, kp1)
                    );

                //fprintf(data, "%19.12E\n", VAT3(uPS, ii, jj, kk));

                // *************************************************************
                // *** > uPE;
                // *************************************************************

                VAT3( uPE, ii, jj, kk) =
                    (
                        VAT3(  uC, ip1,   j,   k) * VAT3( oPE, ii, jj, kk)
                      + VAT3(  oE,   i,   j, kp1) * VAT3( uPC, ii, jj, kk)
                    ) / (
                        VAT3(  oC, ip1,   j, kp1)
                      - VAT3(  oN, ip1,   j, kp1)
                      - VAT3(  oN, ip1, jm1, kp1)
                    );

                //fprintf(data, "%19.12E\n", VAT3(uPE, ii, jj, kk));

                // *************************************************************
                // *** > uPW;
                // *************************************************************

                VAT3( uPW, ii, jj, kk) =
                    (
                        VAT3(  uC, im1,   j,   k) * VAT3( oPW, ii, jj, kk)
                     +  VAT3(  oE, im1,   j, kp1) * VAT3( uPC, ii, jj, kk)
                    ) / (
                        VAT3(  oC, im1,   j, kp1)
                      - VAT3(  oN, im1,   j, kp1)
                      - VAT3(  oN, im1, jm1, kp1)
                    );

                //fprintf(data, "%19.12E\n", VAT3(uPW, ii, jj, kk));

                // *************************************************************
                // *** > uPNE;
                // *************************************************************

                VAT3(uPNE, ii, jj, kk) =
                    (
                        VAT3(  uC, ip1, jp1,   k) * VAT3(oPNE, ii, jj, kk)
                      + VAT3(  oE,   i, jp1, kp1) * VAT3( uPN, ii, jj, kk)
                      + VAT3(  oN, ip1,   j, kp1) * VAT3( uPE, ii, jj, kk)
                    ) / VAT3(  oC, ip1, jp1, kp1);

                //fprintf(data, "%19.12E\n", VAT3(uPNE, ii, jj, kk));

                // *************************************************************
                // *** > uPNW;
                // *************************************************************

                VAT3(uPNW, ii, jj, kk) =
                    (
                        VAT3(  uC, im1, jp1,   k) * VAT3(oPNW, ii, jj, kk)
                      + VAT3(  oE, im1, jp1, kp1) * VAT3( uPN, ii, jj, kk)
                      + VAT3(  oN, im1,   j, kp1) * VAT3( uPW, ii, jj, kk)
                    ) / VAT3(  oC, im1, jp1, kp1);

                //fprintf(data, "%19.12E\n", VAT3(uPNW, ii, jj, kk));

                // *************************************************************
                // *** > uPSE;
                // *************************************************************

                VAT3(uPSE, ii, jj, kk) =
                    (
                        VAT3(  uC, ip1, jm1,   k) * VAT3(oPSE, ii, jj, kk)
                      + VAT3(  oE,   i, jm1, kp1) * VAT3( uPS, ii, jj, kk)
                      + VAT3(  oN, ip1, jm1, kp1) * VAT3( uPE, ii, jj, kk)
                    ) / VAT3(  oC, ip1, jm1, kp1);

                //fprintf(data, "%19.12E\n", VAT3(uPSE, ii, jj, kk));

                // *************************************************************
                // *** > uPSW;
                // *************************************************************

                VAT3(uPSW, ii, jj, kk) =
                    (
                        VAT3(  uC, im1, jm1,   k) * VAT3(oPSW, ii, jj, kk)
                      + VAT3(  oE, im1, jm1, kp1) * VAT3( uPS, ii, jj, kk)
                      + VAT3(  oN, im1, jm1, kp1) * VAT3( uPW, ii, jj, kk)
                    ) / VAT3(  oC, im1, jm1, kp1);

                //fprintf(data, "%19.12E\n", VAT3(uPSW, ii, jj, kk));

            }
        }
    }
}



VPUBLIC void VbuildP_op27(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        int *ipc, double *rpc,
        double *ac, double *pc) {

    MAT2(ac, *nxf * *nyf * *nzf, 1);
    MAT2(pc, *nxc * *nyc * *nzc, 1);

    WARN_UNTESTED;

    VbuildPb_op27(nxf, nyf, nzf,
            nxc, nyc, nzc,
            ipc, rpc,

            RAT2(ac, 1,  1), RAT2(ac, 1,  2), RAT2(ac, 1,  3),
            RAT2(ac, 1,  4),
            RAT2(ac, 1,  5), RAT2(ac, 1,  6),
            RAT2(ac, 1,  7), RAT2(ac, 1,  8), RAT2(ac, 1,  9), RAT2(ac, 1, 10),
            RAT2(ac, 1, 11), RAT2(ac, 1, 12), RAT2(ac, 1, 13), RAT2(ac, 1, 14),
            RAT2(pc, 1,  1), RAT2(pc, 1,  2), RAT2(pc, 1,  3), RAT2(pc, 1,  4), RAT2(pc, 1,  5),
            RAT2(pc, 1,  6), RAT2(pc, 1,  7), RAT2(pc, 1,  8), RAT2(pc, 1,  9),
            RAT2(pc, 1, 10), RAT2(pc, 1, 11), RAT2(pc, 1, 12), RAT2(pc, 1, 13), RAT2(pc, 1, 14),
            RAT2(pc, 1, 15), RAT2(pc, 1, 16), RAT2(pc, 1, 17), RAT2(pc, 1, 18),
            RAT2(pc, 1, 19), RAT2(pc, 1, 20), RAT2(pc, 1, 21), RAT2(pc, 1, 22), RAT2(pc, 1, 23),
            RAT2(pc, 1, 24), RAT2(pc, 1, 25), RAT2(pc, 1, 26), RAT2(pc, 1, 27));
}

VPUBLIC void VbuildPb_op27(int *nxf, int *nyf, int *nzf,
        int *nxc, int *nyc, int *nzc,
        int *ipc, double *rpc,
        double   *oC, double   *oE, double *oN,
        double   *uC,
        double  *oNE, double  *oNW,
        double   *uE, double   *uW, double   *uN, double   *uS,
        double  *uNE, double  *uNW, double  *uSE, double  *uSW,
        double  *oPC, double  *oPN, double  *oPS, double  *oPE, double  *oPW,
        double *oPNE, double *oPNW, double *oPSE, double *oPSW,
        double  *uPC, double  *uPN, double  *uPS, double  *uPE, double  *uPW,
        double *uPNE, double *uPNW, double *uPSE, double *uPSW,
        double  *dPC, double  *dPN, double  *dPS, double  *dPE, double  *dPW,
        double *dPNE, double *dPNW, double *dPSE, double *dPSW) {

      int   i,     j,  k;
      int  ii,    jj, kk;
      int  im1,  ip1;
      int  im2,  ip2;
      int  jm1,  jp1;
      int  jm2,  jp2;
      int  km1,  kp1;
      int  km2,  kp2;
      int iim1, iip1;
      int jjm1, jjp1;
      int kkm1, kkp1;

      double won, half, quarter, eighth;

      MAT3(  oC, *nxf, *nyf, *nzf);
      MAT3(  oE, *nxf, *nyf, *nzf);
      MAT3(  oN, *nxf, *nyf, *nzf);
      MAT3(  uC, *nxf, *nyf, *nzf);
      MAT3( oNE, *nxf, *nyf, *nzf);
      MAT3( oNW, *nxf, *nyf, *nzf);
      MAT3(  uE, *nxf, *nyf, *nzf);
      MAT3(  uW, *nxf, *nyf, *nzf);
      MAT3(  uN, *nxf, *nyf, *nzf);
      MAT3(  uS, *nxf, *nyf, *nzf);
      MAT3( uNE, *nxf, *nyf, *nzf);
      MAT3( uNW, *nxf, *nyf, *nzf);
      MAT3( uSE, *nxf, *nyf, *nzf);
      MAT3( uSW, *nxf, *nyf, *nzf);
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

      WARN_UNTESTED;

      // Interpolation Stencil
      won     = 1.0;
      half    = 1.0 / 2.0;
      quarter = 1.0 / 4.0;
      eighth  = 1.0 / 8.0;

      //fprintf(data, "%s\n", PRINT_FUNC);

      for (kk = 2; kk <= *nzc - 1; kk++) {
          k = 2 * kk - 1;

          for (jj = 2; jj <= *nyc - 1; jj++) {
              j = 2 * jj - 1;

              for (ii = 2; ii <= *nxc - 1; ii++) {
                  i = 2 * ii - 1;

                  // Index Computations
                  im1  =  i - 1;
                  ip1  =  i + 1;
                  im2  =  i - 2;
                  ip2  =  i + 2;
                  jm1  =  j - 1;
                  jp1  =  j + 1;
                  jm2  =  j - 2;
                  jp2  =  j + 2;
                  km1  =  k - 1;
                  kp1  =  k + 1;
                  km2  =  k - 2;
                  kp2  =  k + 2;
                  iim1 = ii - 1;
                  iip1 = ii + 1;
                  jjm1 = jj - 1;
                  jjp1 = jj + 1;
                  kkm1 = kk - 1;
                  kkp1 = kk + 1;

                  //* **********************************************************
                  //* *** > oPC;
                  //* **********************************************************

                  VAT3( oPC,  ii,  jj,  kk) = won;

                  //fprintf(data, "%19.12E\n", VAT3(oPC, ii, jj, kk));

                  //* **********************************************************
                  //* *** > oPN;
                  //* **********************************************************

                  VAT3( oPN,  ii,  jj,  kk) =
                          (
                                VAT3( uNE, im1,   j, km1)
                              + VAT3(  uN,   i,   j, km1)
                              + VAT3( uNW, ip1,   j, km1)
                              + VAT3( oNE, im1,   j,   k)
                              + VAT3(  oN,   i,   j,   k)
                              + VAT3( oNW, ip1,   j,   k)
                              + VAT3( uSW,   i, jp1,   k)
                              + VAT3(  uS,   i, jp1,   k)
                              + VAT3( uSE,   i, jp1,   k)
                          ) / (
                                VAT3(  oC,   i, jp1,   k)
                              - VAT3(  oE, im1, jp1,   k)
                              - VAT3(  oE,   i, jp1,   k)
                              - VAT3(  uC,   i, jp1, km1)
                              - VAT3(  uE, im1, jp1, km1)
                              - VAT3(  uW, ip1, jp1, km1)
                              - VAT3(  uC,   i, jp1,   k)
                              - VAT3(  uW,   i, jp1,   k)
                              - VAT3(  uE,   i, jp1,   k)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(oPN, ii, jj, kk));

                  //* **********************************************************
                  //* *** > oPS;
                  //* **********************************************************

                  VAT3( oPS,  ii,  jj,  kk) =
                          (
                                VAT3( uSE, im1,   j, km1)
                              + VAT3(  uS,   i,   j, km1)
                              + VAT3( uSW, ip1,   j, km1)
                              + VAT3( oNW,   i, jm1,   k)
                              + VAT3(  oN,   i, jm1,   k)
                              + VAT3( oNE,   i, jm1,   k)
                              + VAT3( uNW,   i, jm1,   k)
                              + VAT3(  uN,   i, jm1,   k)
                              + VAT3( uNE,   i, jm1,   k)
                          ) / (
                                VAT3(  oC,   i, jm1,   k)
                              - VAT3(  oE, im1, jm1,   k)
                              - VAT3(  oE,   i, jm1,   k)
                              - VAT3(  uC,   i, jm1, km1)
                              - VAT3(  uE, im1, jm1, km1)
                              - VAT3(  uW, ip1, jm1, km1)
                              - VAT3(  uC,   i, jm1,   k)
                              - VAT3(  uW,   i, jm1,   k)
                              - VAT3(  uE,   i, jm1,   k)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(oPS, ii, jj, kk));

                  //* **********************************************************
                  //* *** > oPE;
                  //* **********************************************************

                  VAT3( oPE,  ii,  jj,  kk) =
                          (
                                VAT3( uSE,   i, jp1, km1)
                              + VAT3( oNW, ip1,   j,   k)
                              + VAT3( uNW, ip1,   j,   k)
                              + VAT3(  uE,   i,   j, km1)
                              + VAT3(  oE,   i,   j,   k)
                              + VAT3(  uW, ip1,   j,   k)
                              + VAT3( uNE,   i, jm1, km1)
                              + VAT3( oNE,   i, jm1,   k)
                              + VAT3( uSW, ip1,   j,   k)
                          ) / (
                                VAT3(  oC, ip1,   j,   k)
                              - VAT3(  uC, ip1,   j, km1)
                              - VAT3(  uC, ip1,   j,   k)
                              - VAT3(  oN, ip1,   j,   k)
                              - VAT3(  uS, ip1, jp1, km1)
                              - VAT3(  uN, ip1,   j,   k)
                              - VAT3(  oN, ip1, jm1,   k)
                              - VAT3(  uN, ip1, jm1, km1)
                              - VAT3(  uS, ip1,   j,   k)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(oPE, ii, jj, kk));

                  //* **********************************************************
                  //* *** > oPW;
                  //* **********************************************************

                  VAT3( oPW,  ii,  jj,  kk) =
                          (
                                VAT3( uSW,   i, jp1, km1)
                              + VAT3( oNE, im1,   j,   k)
                              + VAT3( uNE, im1,   j,   k)
                              + VAT3(  uW,   i,   j, km1)
                              + VAT3(  oE, im1,   j,   k)
                              + VAT3(  uE, im1,   j,   k)
                              + VAT3( uNW,   i, jm1, km1)
                              + VAT3( oNW,   i, jm1,   k)
                              + VAT3( uSE, im1,   j,   k)
                          ) / (
                                VAT3(  oC, im1,   j,   k)
                              - VAT3(  uC, im1,   j, km1)
                              - VAT3(  uC, im1,   j,   k)
                              - VAT3(  oN, im1,   j,   k)
                              - VAT3(  uS, im1, jp1, km1)
                              - VAT3(  uN, im1,   j,   k)
                              - VAT3(  oN, im1, jm1,   k)
                              - VAT3(  uN, im1, jm1, km1)
                              - VAT3(  uS, im1,   j,   k)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(oPW, ii, jj, kk));

                  //* **********************************************************
                  //* *** > oPNE;
                  //* **********************************************************

                  VAT3(oPNE,  ii,  jj,  kk) =
                          (
                                VAT3( uNE,   i,   j, km1)
                              + VAT3( oNE,   i,   j,   k)
                              + VAT3( uSW, ip1, jp1,   k)
                              + (
                                    VAT3(  uN, ip1,   j, km1)
                                  + VAT3(  oN, ip1,   j,   k)
                                  + VAT3(  uS, ip1, jp1,   k)
                                )
                              * VAT3( oPE,  ii,  jj,  kk)
                              + (
                                    VAT3(  uE,   i, jp1, km1)
                                  + VAT3(  oE,   i, jp1,   k)
                                  + VAT3(  uW, ip1, jp1,   k)
                                )
                              * VAT3( oPN,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC, ip1, jp1,   k)
                              - VAT3(  uC, ip1, jp1, km1)
                              - VAT3(  uC, ip1, jp1,   k)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(oPNE, ii, jj, kk));

                  //* **********************************************************
                  //* *** > oPNW;
                  //* **********************************************************

                  VAT3(oPNW,  ii,  jj,  kk) =
                          (
                                VAT3( uNW,   i,   j, km1)
                              + VAT3( oNW,   i,   j,   k)
                              + VAT3( uSE, im1, jp1,   k)
                              + (
                                    VAT3(  uN, im1,   j, km1)
                                  + VAT3(  oN, im1,   j,   k)
                                  + VAT3(  uS, im1, jp1,   k)
                                )
                              * VAT3( oPW,  ii,  jj,  kk)
                              + (
                                    VAT3(  uW,   i, jp1, km1)
                                  + VAT3(  oE, im1, jp1,   k)
                                  + VAT3(  uE, im1, jp1,   k)
                                )
                              * VAT3( oPN,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC, im1, jp1,   k)
                              - VAT3(  uC, im1, jp1, km1)
                              - VAT3(  uC, im1, jp1,   k)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(oPNW, ii, jj, kk));

                  //* **********************************************************
                  //* *** > oPSE;
                  //* **********************************************************

                  VAT3(oPSE,  ii,  jj,  kk) =
                          (
                                VAT3( uSE,   i,   j, km1)
                              + VAT3( oNW, ip1, jm1,   k)
                              + VAT3( uNW, ip1, jm1,   k)
                              + (
                                    VAT3(  uS, ip1,   j, km1)
                                  + VAT3(  oN, ip1, jm1,   k)
                                  + VAT3(  uN, ip1, jm1,   k)
                                )
                              * VAT3( oPE,  ii,  jj,  kk)
                              + (
                                    VAT3(  uE,   i, jm1, km1)
                                  + VAT3(  oE,   i, jm1,   k)
                                  + VAT3(  uW, ip1, jm1,   k)
                                )
                              * VAT3( oPS,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC, ip1, jm1,   k)
                              - VAT3(  uC, ip1, jm1, km1)
                              - VAT3(  uC, ip1, jm1,   k)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(oPSE, ii, jj, kk));

                  //* **********************************************************
                  //* *** > oPSW;
                  //* **********************************************************

                  VAT3(oPSW,  ii,  jj,  kk) =
                          (
                                VAT3( uSW,   i,   j, km1)
                              + VAT3( oNE, im1, jm1,   k)
                              + VAT3( uNE, im1, jm1,   k)
                              + (
                                    VAT3(  uS, im1,   j, km1)
                                  + VAT3(  oN, im1, jm1,   k)
                                  + VAT3(  uN, im1, jm1,   k)
                                )
                              * VAT3( oPW,  ii,  jj,  kk)
                              + (
                                    VAT3(  uW,   i, jm1, km1)
                                  + VAT3(  oE, im1, jm1,   k)
                                  + VAT3(  uE, im1, jm1,   k)
                                )
                              * VAT3( oPS,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC, im1, jm1,   k)
                              - VAT3(  uC, im1, jm1, km1)
                              - VAT3(  uC, im1, jm1,   k)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(oPSW, ii, jj, kk));

                  //* **********************************************************
                  //* *** > dPC;
                  //* **********************************************************

                  VAT3( dPC,  ii,  jj,  kk) =
                          (
                                VAT3( uNW,   i,   j, km1)
                              + VAT3(  uW,   i,   j, km1)
                              + VAT3( uSW,   i,   j, km1)
                              + VAT3(  uN,   i,   j, km1)
                              + VAT3(  uC,   i,   j, km1)
                              + VAT3(  uS,   i,   j, km1)
                              + VAT3( uNE,   i,   j, km1)
                              + VAT3(  uE,   i,   j, km1)
                              + VAT3( uSE,   i,   j, km1)
                          ) / (
                                VAT3(  oC,   i,   j, km1)
                              - VAT3(  oN,   i,   j, km1)
                              - VAT3(  oN,   i, jm1, km1)
                              - VAT3( oNW,   i,   j, km1)
                              - VAT3(  oE, im1,   j, km1)
                              - VAT3( oNE, im1, jm1, km1)
                              - VAT3( oNE,   i,   j, km1)
                              - VAT3(  oE,   i,   j, km1)
                              - VAT3( oNW, ip1, jm1, km1)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(dPC, ii, jj, kk));

                  //* **********************************************************
                  //* *** > dPN;
                  //* **********************************************************

                  VAT3( dPN,  ii,  jj,  kk) =
                          (
                                VAT3( uSW,   i, jp1, km1)
                              + VAT3(  uS,   i, jp1, km1)
                              + VAT3( uSE,   i, jp1, km1)
                              + (
                                    VAT3( oNE, im1,   j, km1)
                                  + VAT3(  oN,   i,   j, km1)
                                  + VAT3( oNW, ip1,   j, km1)
                                )
                              * VAT3( dPC,  ii,  jj,  kk)
                              + (
                                    VAT3(  uW,   i, jp1, km1)
                                  + VAT3(  uC,   i, jp1, km1)
                                  + VAT3(  uE,   i, jp1, km1)
                                )
                              * VAT3( oPN,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC,   i, jp1, km1)
                              - VAT3(  oE, im1, jp1, km1)
                              - VAT3(  oE,   i, jp1, km1)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(dPN, ii, jj, kk));

                  //* **********************************************************
                  //* *** > dPS;
                  //* **********************************************************

                  VAT3( dPS,  ii,  jj,  kk) =
                          (
                                VAT3( uNW,   i, jm1, km1)
                              + VAT3(  uN,   i, jm1, km1)
                              + VAT3( uNE,   i, jm1, km1)
                              + (
                                    VAT3( oNW,   i, jm1, km1)
                                  + VAT3(  oN,   i, jm1, km1)
                                  + VAT3( oNE,   i, jm1, km1)
                                )
                              * VAT3( dPC,  ii,  jj,  kk)
                              + (
                                    VAT3(  uW,   i, jm1, km1)
                                  + VAT3(  uC,   i, jm1, km1)
                                  + VAT3(  uE,   i, jm1, km1)
                                )
                              * VAT3( oPS,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC,   i, jm1, km1)
                              - VAT3(  oE, im1, jm1, km1)
                              - VAT3(  oE,   i, jm1, km1)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(dPS, ii, jj, kk));

                  //* **********************************************************
                  //* *** > dPE;
                  //* **********************************************************

                  VAT3( dPE,  ii,  jj,  kk) =
                          (
                                VAT3( uNW, ip1,   j, km1)
                              + VAT3(  uW, ip1,   j, km1)
                              + VAT3( uSW, ip1,   j, km1)
                              + (
                                    VAT3(  uN, ip1,   j, km1)
                                  + VAT3(  uC, ip1,   j, km1)
                                  + VAT3(  uS, ip1,   j, km1)
                                )
                              * VAT3( oPE,  ii,  jj,  kk)
                              + (
                                    VAT3( oNW, ip1,   j, km1)
                                  + VAT3(  oE,   i,   j, km1)
                                  + VAT3( oNE,   i, jm1, km1)
                                )
                              * VAT3( dPC,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC, ip1,   j, km1)
                              - VAT3(  oN, ip1,   j, km1)
                              - VAT3(  oN, ip1, jm1, km1)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(dPE, ii, jj, kk));

                  //* **********************************************************
                  //* *** > dPW;
                  //* **********************************************************

                  VAT3( dPW,  ii,  jj,  kk) =
                          (
                                VAT3( uNE, im1,   j, km1)
                              + VAT3(  uE, im1,   j, km1)
                              + VAT3( uSE, im1,   j, km1)
                              + (
                                    VAT3(  uN, im1,   j, km1)
                                  + VAT3(  uC, im1,   j, km1)
                                  + VAT3(  uS, im1,   j, km1)
                                )
                              * VAT3( oPW,  ii,  jj,  kk)
                              + (
                                    VAT3( oNE, im1,   j, km1)
                                  + VAT3(  oE, im1,   j, km1)
                                  + VAT3( oNW,   i, jm1, km1)
                                )
                              * VAT3( dPC,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC, im1,   j, km1)
                              - VAT3(  oN, im1,   j, km1)
                              - VAT3(  oN, im1, jm1, km1)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(dPW, ii, jj, kk));

                  //* **********************************************************
                  //* *** > dPNE;
                  //* **********************************************************

                  VAT3(dPNE,  ii,  jj,  kk) =
                          (
                                VAT3( uSW, ip1, jp1, km1)
                              + VAT3(  uW, ip1, jp1, km1)
                              * VAT3( oPN,  ii,  jj,  kk)
                              + VAT3(  uS, ip1, jp1, km1)
                              * VAT3( oPE,  ii,  jj,  kk)
                              + VAT3(  uC, ip1, jp1, km1)
                              * VAT3(oPNE,  ii,  jj,  kk)
                              + VAT3( oNE,   i,   j, km1)
                              * VAT3( dPC,  ii,  jj,  kk)
                              + VAT3(  oE,   i, jp1, km1)
                              * VAT3( dPN,  ii,  jj,  kk)
                              + VAT3(  oN, ip1,   j, km1)
                              * VAT3( dPE,  ii,  jj,  kk)
                          )
                          / VAT3(  oC, ip1, jp1, km1);

                  //fprintf(data, "%19.12E\n", VAT3(dPNE, ii, jj, kk));

                  //* **********************************************************
                  //* *** > dPNW;
                  //* **********************************************************

                  VAT3(dPNW,  ii,  jj,  kk) =
                          (
                                VAT3( uSE, im1, jp1, km1)
                              + VAT3(  uE, im1, jp1, km1)
                              * VAT3( oPN,  ii,  jj,  kk)
                              + VAT3(  uS, im1, jp1, km1)
                              * VAT3( oPW,  ii,  jj,  kk)
                              + VAT3(  uC, im1, jp1, km1)
                              * VAT3(oPNW,  ii,  jj,  kk)
                              + VAT3( oNW,   i,   j, km1)
                              * VAT3( dPC,  ii,  jj,  kk)
                              + VAT3(  oE, im1, jp1, km1)
                              * VAT3( dPN,  ii,  jj,  kk)
                              + VAT3(  oN, im1,   j, km1)
                              * VAT3( dPW,  ii,  jj,  kk)
                          )
                          / VAT3(  oC, im1, jp1, km1);

                  //fprintf(data, "%19.12E\n", VAT3(dPNW, ii, jj, kk));

                  //* **********************************************************
                  //* *** > dPSE;
                  //* **********************************************************

                  VAT3(dPSE,  ii,  jj,  kk) =
                          (
                                VAT3( uNW, ip1, jm1, km1)
                              + VAT3(  uW, ip1, jm1, km1)
                              * VAT3( oPS,  ii,  jj,  kk)
                              + VAT3(  uN, ip1, jm1, km1)
                              * VAT3( oPE,  ii,  jj,  kk)
                              + VAT3(  uC, ip1, jm1, km1)
                              * VAT3(oPSE,  ii,  jj,  kk)
                              + VAT3( oNW, ip1, jm1, km1)
                              * VAT3( dPC,  ii,  jj,  kk)
                              + VAT3(  oE,   i, jm1, km1)
                              * VAT3( dPS,  ii,  jj,  kk)
                              + VAT3(  oN, ip1, jm1, km1)
                              * VAT3( dPE,  ii,  jj,  kk)
                          )
                          / VAT3(  oC, ip1, jm1, km1);

                  //fprintf(data, "%19.12E\n", VAT3(dPSE, ii, jj, kk));

                  //* **********************************************************
                  //* *** > dPSW;
                  //* **********************************************************

                  VAT3(dPSW,  ii,  jj,  kk) =
                          (
                                VAT3( uNE, im1, jm1, km1)
                              + VAT3(  uE, im1, jm1, km1)
                              * VAT3( oPS,  ii,  jj,  kk)
                              + VAT3(  uN, im1, jm1, km1)
                              * VAT3( oPW,  ii,  jj,  kk)
                              + VAT3(  uC, im1, jm1, km1)
                              * VAT3(oPSW,  ii,  jj,  kk)
                              + VAT3( oNE, im1, jm1, km1)
                              * VAT3( dPC,  ii,  jj,  kk)
                              + VAT3(  oE, im1, jm1, km1)
                              * VAT3( dPS,  ii,  jj,  kk)
                              + VAT3(  oN, im1, jm1, km1)
                              * VAT3( dPW,  ii,  jj,  kk)
                          )
                          / VAT3(  oC, im1, jm1, km1);

                  //fprintf(data, "%19.12E\n", VAT3(dPSW, ii, jj, kk));

                  //* **********************************************************
                  //* *** > uPC;
                  //* **********************************************************

                  VAT3( uPC,  ii,  jj,  kk) =
                          (
                                VAT3( uSE, im1, jp1,   k)
                              + VAT3(  uE, im1,   j,   k)
                              + VAT3( uNE, im1, jm1,   k)
                              + VAT3(  uS,   i, jp1,   k)
                              + VAT3(  uC,   i,   j,   k)
                              + VAT3(  uN,   i, jm1,   k)
                              + VAT3( uSW, ip1, jp1,   k)
                              + VAT3(  uW, ip1,   j,   k)
                              + VAT3( uNW, ip1, jm1,   k)
                          ) / (
                                VAT3(  oC,   i,   j, kp1)
                              - VAT3(  oN,   i,   j, kp1)
                              - VAT3(  oN,   i, jm1, kp1)
                              - VAT3( oNW,   i,   j, kp1)
                              - VAT3(  oE, im1,   j, kp1)
                              - VAT3( oNE, im1, jm1, kp1)
                              - VAT3( oNE,   i,   j, kp1)
                              - VAT3(  oE,   i,   j, kp1)
                              - VAT3( oNW, ip1, jm1, kp1)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(uPC, ii, jj, kk));

                  //* **********************************************************
                  //* *** > uPN;
                  //* **********************************************************

                  VAT3( uPN,  ii,  jj,  kk) =
                          (
                                VAT3( uNE, im1,   j,   k)
                              + VAT3(  uN,   i,   j,   k)
                              + VAT3( uNW, ip1,   j,   k)
                              + (
                                    VAT3( oNE, im1,   j, kp1)
                                  + VAT3(  oN,   i,   j, kp1)
                                  + VAT3( oNW, ip1,   j, kp1)
                                )
                              * VAT3( uPC,  ii,  jj,  kk)
                              + (
                                    VAT3(  uE, im1, jp1,   k)
                                  + VAT3(  uC,   i, jp1,   k)
                                  + VAT3(  uW, ip1, jp1,   k)
                                )
                              * VAT3( oPN,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC,   i, jp1, kp1)
                              - VAT3(  oE, im1, jp1, kp1)
                              - VAT3(  oE,   i, jp1, kp1)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(uPN, ii, jj, kk));

                  //* **********************************************************
                  //* *** > uPS;
                  //* **********************************************************

                  VAT3( uPS,  ii,  jj,  kk) =
                          (
                                VAT3( uSE, im1,   j,   k)
                              + VAT3(  uS,   i,   j,   k)
                              + VAT3( uSW, ip1,   j,   k)
                              + (
                                    VAT3( oNW,   i, jm1, kp1)
                              + VAT3(  oN,   i, jm1, kp1)
                              + VAT3( oNE,   i, jm1, kp1)
                                )
                              * VAT3( uPC,  ii,  jj,  kk)
                              + (
                                    VAT3(  uE, im1, jm1,   k)
                              + VAT3(  uC,   i, jm1,   k)
                              + VAT3(  uW, ip1, jm1,   k)
                                )
                              * VAT3( oPS,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC,   i, jm1, kp1)
                              - VAT3(  oE, im1, jm1, kp1)
                              - VAT3(  oE,   i, jm1, kp1)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(uPS, ii, jj, kk));

                  //* **********************************************************
                  //* *** > uPE;
                  //* **********************************************************

                  VAT3( uPE,  ii,  jj,  kk) =
                          (
                                VAT3( uSE,   i, jp1,   k)
                              + VAT3(  uS, ip1, jp1,   k)
                              + VAT3( uNE,   i, jm1,   k)
                              + (
                                    VAT3(  uS, ip1, jp1,   k)
                              + VAT3(  uC, ip1,   j,   k)
                              + VAT3(  uN, ip1, jm1,   k)
                                )
                              * VAT3( oPE,  ii,  jj,  kk)
                              + (
                                    VAT3( oNW, ip1,   j, kp1)
                              + VAT3(  oE,   i,   j, kp1)
                              + VAT3( oNE,   i, jm1, kp1)
                                )
                              * VAT3( uPC,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC, ip1,   j, kp1)
                              - VAT3(  oN, ip1,   j, kp1)
                              - VAT3(  oN, ip1, jm1, kp1)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(uPE, ii, jj, kk));

                  //* **********************************************************
                  //* *** > uPW;
                  //* **********************************************************

                  VAT3( uPW,  ii,  jj,  kk) =
                          (
                                VAT3( uSW,   i, jp1,   k)
                              + VAT3(  uW,   i,   j,   k)
                              + VAT3( uNW,   i, jm1,   k)
                              + (
                                    VAT3(  uS, im1, jp1,   k)
                                  + VAT3(  uC, im1,   j,   k)
                                  + VAT3(  uN, im1, jm1,   k)
                                )
                              * VAT3( oPW,  ii,  jj,  kk)
                              + (
                                    VAT3( oNE, im1,   j, kp1)
                                  + VAT3(  oE, im1,   j, kp1)
                                  + VAT3( oNW,   i, jm1, kp1)
                                )
                              * VAT3( uPC,  ii,  jj,  kk)
                          ) / (
                                VAT3(  oC, im1,   j, kp1)
                              - VAT3(  oN, im1,   j, kp1)
                              - VAT3(  oN, im1, jm1, kp1)
                          );

                  //fprintf(data, "%19.12E\n", VAT3(uPW, ii, jj, kk));

                  //* **********************************************************
                  //* *** > uPNE;
                  //* **********************************************************

                  VAT3(uPNE,  ii,  jj,  kk) =
                          (
                                VAT3( uNE,   i,   j,   k)
                              + VAT3(  uE,   i, jp1,   k)
                              * VAT3( oPN,  ii,  jj,  kk)
                              + VAT3(  uN, ip1,   j,   k)
                              * VAT3( oPE,  ii,  jj,  kk)
                              + VAT3(  uC, ip1, jp1,   k)
                              * VAT3(oPNE,  ii,  jj,  kk)
                              + VAT3( oNE,   i,   j, kp1)
                              * VAT3( uPC,  ii,  jj,  kk)
                              + VAT3(  oE,   i, jp1, kp1)
                              * VAT3( uPN,  ii,  jj,  kk)
                              + VAT3(  oN, ip1,   j, kp1)
                              * VAT3( uPE,  ii,  jj,  kk)
                          )
                          / VAT3(  oC, ip1, jp1, kp1);

                  //fprintf(data, "%19.12E\n", VAT3(uPNE, ii, jj, kk));

                  //* **********************************************************
                  //* *** > uPNW;
                  //* **********************************************************

                  VAT3(uPNW,  ii,  jj,  kk) =
                          (
                                VAT3( uNW,   i,   j,   k)
                              + VAT3(  uW,   i, jp1,   k)
                              * VAT3( oPN,  ii,  jj,  kk)
                              + VAT3(  uN, im1,   j,   k)
                              * VAT3( oPW,  ii,  jj,  kk)
                              + VAT3(  uC, im1, jp1,   k)
                              * VAT3(oPNW,  ii,  jj,  kk)
                              + VAT3( oNW,   i,   j, kp1)
                              * VAT3( uPC,  ii,  jj,  kk)
                              + VAT3(  oE, im1, jp1, kp1)
                              * VAT3( uPN,  ii,  jj,  kk)
                              + VAT3(  oN, im1,   j, kp1)
                              * VAT3( uPW,  ii,  jj,  kk)
                          )
                          / VAT3(  oC, im1, jp1, kp1);

                  //fprintf(data, "%19.12E\n", VAT3(uPNW, ii, jj, kk));

                  //* **********************************************************
                  //* *** > uPSE;
                  //* **********************************************************

                  VAT3(uPSE,  ii,  jj,  kk) =
                          (
                                VAT3( uSE,   i,   j,   k)
                              + VAT3(  uE,   i, jm1,   k)
                              * VAT3( oPS,  ii,  jj,  kk)
                              + VAT3(  uS, ip1,   j,   k)
                              * VAT3( oPE,  ii,  jj,  kk)
                              + VAT3(  uC, ip1, jm1,   k)
                              * VAT3(oPSE,  ii,  jj,  kk)
                              + VAT3( oNW, ip1, jm1, kp1)
                              * VAT3( uPC,  ii,  jj,  kk)
                              + VAT3(  oE,   i, jm1, kp1)
                              * VAT3( uPS,  ii,  jj,  kk)
                              + VAT3(  oN, ip1, jm1, kp1)
                              * VAT3( uPE,  ii,  jj,  kk)
                          )
                          / VAT3(  oC, ip1, jm1, kp1);

                  //fprintf(data, "%19.12E\n", VAT3(uPSE, ii, jj, kk));

                  //* **********************************************************
                  //* *** > uPSW;
                  //* **********************************************************

                  VAT3(uPSW,  ii,  jj,  kk) =
                          (
                                VAT3( uSW,   i,   j,   k)
                              + VAT3(  uW,   i, jm1,   k)
                              * VAT3( oPS,  ii,  jj,  kk)
                              + VAT3(  uS, im1,   j,   k)
                              * VAT3( oPW,  ii,  jj,  kk)
                              + VAT3(  uC, im1, jm1,   k)
                              * VAT3(oPSW,  ii,  jj,  kk)
                              + VAT3( oNE, im1, jm1, kp1)
                              * VAT3( uPC,  ii,  jj,  kk)
                              + VAT3(  oE, im1, jm1, kp1)
                              * VAT3( uPS,  ii,  jj,  kk)
                              + VAT3(  oN, im1, jm1, kp1)
                              * VAT3( uPW,  ii,  jj,  kk)
                          )
                          / VAT3(  oC, im1, jm1, kp1);

                  //fprintf(data, "%19.12E\n", VAT3(uPSW, ii, jj, kk));

              }
          }
      }
}
