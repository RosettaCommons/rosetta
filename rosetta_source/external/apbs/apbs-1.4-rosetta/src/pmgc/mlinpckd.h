/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief LINPACK interface
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

#ifndef MLINPCKD_H_
#define MLINPCKD_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vmatrix.h"

/** @brief   Solves the double precision symmetric positive definite band system
 *           A*X = B using the factors computed by dpbco or dpbfa
 *  @ingroup  PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    A division by zero will occur if the input factor contains
 *           a zero on the diagonal.  Technically this indicates singularity,
 *           but it is usually caused by imporper subroutine arugments.  It will
 *           not occur if the subroutines are called correctly and info == 0
 *  @note    Replaces dpbsl from mgsubd.f
 */
VEXTERNC void Vdpbsl(
        double *abd, ///< The output from dpbco or dpbfa
        int *lda,    ///< The leading dimension of the array abd
        int *n,      ///< The order of the matrix a
        int *m,      ///< The number of diagonals above the main diagonal
        double *b    ///< The right hand side vector
        );

/** Translation of LINPACK daxpy subroutine
 * *     jack dongarra, linpack, 3/11/78.
 * */
VEXTERNC void Vdaxpy(int n, double da,
        double *dx, int incx,
        double *dy, int incy);

/** Translation of LINPACK ddot subroutine
 * *     jack dongarra, linpack, 3/11/78.
 * */
VEXTERNC double Vddot(int n, double *dx, int incx, double *dy, int incy);

/** Translation of LINPACK dpbfa subroutine
 * *     jack dongarra, linpack, 3/11/78.
 * */
VEXTERNC void Vdpbfa(double *abd, int *lda, int *n, int *m, int *info);

#endif /* MLINPCKD_H_ */
