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

#include "mlinpckd.h"

VPUBLIC void Vdpbsl(double *abd, int *lda, int *n, int *m, double *b) {

    double t;
    int k, kb, la, lb, lm;

    MAT2(abd, *lda, 1);

    for (k=1; k<=*n; k++) {
        lm = VMIN2(k-1, *m);
        la = *m + 1 - lm;
        lb = k - lm;
        t = Vddot(lm, RAT2(abd, la, k ), 1, RAT(b, lb), 1 );
        VAT(b, k) = (VAT(b, k) - t) / VAT2(abd, *m+1, k);
    }

    // Solve R*X = Y
    for (kb=1; kb<=*n; kb++) {

        k = *n + 1 - kb;
        lm = VMIN2(k-1, *m);
        la = *m + 1 - lm;
        lb = k - lm;
        VAT(b, k) /= VAT2(abd, *m+1, k);
        t = -VAT(b, k);
        Vdaxpy(lm, t, RAT2(abd, la, k), 1, RAT(b, lb), 1);
    }
}

VPUBLIC void Vdaxpy(int n, double da,
        double *dx, int incx,
        double *dy, int incy) {

    int i, ix, iy, m, mp1;

    if (n <= 0)
        return;

    if (da == 0)
        return;

    if (incx == 1 && incy == 1) {

        m = n % 4;
        if (m != 0) {

            for (i=1; i<=m; i++)
                VAT(dy, i) += da * VAT(dx, i);
        }

        if (n < 4)
            return;

        mp1 = m + 1;

        for (i=mp1; i<=n; i+=4) {

            VAT(dy, i  ) += da * VAT(dx, i  );
            VAT(dy, i+1) += da * VAT(dx, i+1);
            VAT(dy, i+2) += da * VAT(dx, i+2);
            VAT(dy, i+3) += da * VAT(dx, i+3);
        }
    } else {

        ix = 1;
        if (incx < 0 )
            ix = (-n + 1) * incx + 1;

        iy = 1;
        if (incy < 0 )
            iy = (-n + 1) * incy + 1;

        for (i=1; i<=n; i++) {

            VAT(dy, iy) += da * VAT(dx, ix);
            ix += incx;
            iy += incy;
        }
    }
}


VPUBLIC double Vddot(int n, double *dx, int incx, double *dy, int incy) {

    double dtemp;
    int i, ix, iy, m, mp1;

    double ddot = 0.0;
    dtemp = 0.0;

    if (n <= 0)
        return ddot;

    if (incx == 1 && incy == 1) {

        m = n % 5;

        if (m != 0) {

            for (i=1; i<=m; i++)
                dtemp += VAT(dx, i) * VAT(dy, i);

            if (n < 5) {
                ddot = dtemp;
                return ddot;
            }
        }

        mp1 = m + 1;

        for (i=mp1; i<=n; i+=5)
            dtemp += VAT(dx,   i) * VAT(dy,   i)
                  +  VAT(dx, i+1) * VAT(dy, i+1)
                  +  VAT(dx, i+2) * VAT(dy, i+2)
                  +  VAT(dx, i+3) * VAT(dy, i+3)
                  +  VAT(dx, i+4) * VAT(dy, i+4);
    } else {

        ix = 1;
        if (incx < 0)
            ix = (-n + 1) * incx + 1;

        iy = 1;
        if (incy < 0)
            iy = (-n + 1) * incy + 1;

        for (i=1; i<=n; i++) {
            ix += incx;
            iy += incy;
        }
    }

    ddot = dtemp;
    return ddot;
}



VPUBLIC void Vdpbfa(double *abd, int *lda, int *n, int *m, int *info) {

    double t, s;
    int ik, j, jk, k, mu;

    MAT2(abd, *lda, 1);

    *info = 0;

    for(j = 1; j <= *n; j++) {

        s = 0.0;
        ik = *m + 1;
        jk = VMAX2(j - *m, 1);
        mu = VMAX2(*m + 2 - j, 1);

        if (*m >= mu ) {

            for(k = mu; k <= *m; k++) {
                t = VAT2(abd, k, j) - Vddot(k - mu,
                        RAT2(abd, ik, jk), 1,
                        RAT2(abd, mu,  j), 1);
                t /= VAT2(abd, *m + 1, jk);
                VAT2(abd, k, j) = t;
                s += t * t;
                ik--;
                jk++;
            }
        }

        s = VAT2(abd, *m + 1, j) - s;

        if (s <= 0.0) {
            *info = j;
            break;
        }

        VAT2(abd, *m + 1, j) = VSQRT(s);
    }
}
