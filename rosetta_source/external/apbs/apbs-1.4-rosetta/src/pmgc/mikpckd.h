/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief   A collection of useful low-level routines (timing, etc).
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

#ifndef MIKPCKD_H_
#define MIKPCKD_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vmatrix.h"

/** @brief   Copy operation for a grid function with boundary values.
 *           Quite simply copies one 3d matrix to another
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces xcopy from mikpckd.f
 */
VEXTERNC void Vxcopy(
        int *nx,   ///< The size of the x dimension of the 3d matrix
        int *ny,   ///< The size of the y dimension of the 3d matrix
        int *nz,   ///< The size of the z dimension of the 3d matrix
        double *x, ///< The source matrix from which to copy data
        double *y  ///< The destination matrix to receive copied data
        );



/** @brief   Copy operation for a grid function with boundary values.
 *           Quite simply copies one 3d matrix to another
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces xcopy_small from mikpckd.f
 */
VEXTERNC void Vxcopy_small(
        int *nx,   ///< The size of the x dimension of the 3d matrix
        int *ny,   ///< The size of the y dimension of the 3d matrix
        int *nz,   ///< The size of the z dimension of the 3d matrix
        double *x, ///< The source matrix from which to copy data
        double *y  ///< The destination matrix to receive copied data
        );



/** @brief   Copy operation for a grid function with boundary values.
 *           Quite simply copies one 3d matrix to another
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces xcopy_large from mikpckd.f
 */
VEXTERNC void Vxcopy_large(
        int *nx,   ///< The size of the x dimension of the 3d matrix
        int *ny,   ///< The size of the y dimension of the 3d matrix
        int *nz,   ///< The size of the z dimension of the 3d matrix
        double *x, ///< The source matrix from which to copy data
        double *y  ///< The destination matrix to receive copied data
        );



/** @brief   saxpy operation for a grid function with boundary values.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces xaxpy from mikpckd.f
 */
VEXTERNC void Vxaxpy(
        int *nx,       ///< The size of the x dimension of the 3d matrix
        int *ny,       ///< The size of the y dimension of the 3d matrix
        int *nz,       ///< The size of the z dimension of the 3d matrix
        double *alpha, ///< @todo: Doc
        double *x,     ///< The source matrix from which to copy data
        double *y      ///< The destination matrix to receive copied data
        );

/** @brief   Norm operation for a grid function with boundary values.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces xnrm1 from mikpckd.f
 */
VEXTERNC double Vxnrm1(
        int *nx,  ///< The size of the x dimension of the 3d matrix
        int *ny,  ///< The size of the y dimension of the 3d matrix
        int *nz,  ///< The size of the z dimension of the 3d matrix
        double *x ///< The matrix to normalize
        );



/** @brief   Norm operation for a grid function with boundary values.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces xnrm2 from mikpckd.f
 */
VEXTERNC double Vxnrm2(
        int *nx,  ///< The size of the x dimension of the 3d matrix
        int *ny,  ///< The size of the y dimension of the 3d matrix
        int *nz,  ///< The size of the z dimension of the 3d matrix
        double *x ///< The matrix to normalize
        );

/** @brief   Inner product operation for a grid function with boundary values.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces xdot from mikpckd.f
 */
VEXTERNC double Vxdot(
        int *nx,   ///< The size of the x dimension of the 3d matrix
        int *ny,   ///< The size of the y dimension of the 3d matrix
        int *nz,   ///< The size of the z dimension of the 3d matrix
        double *x, ///< The first vector
        double *y  ///< The second vector
        );

/** @brief   Zero out operation for a grid function, including boundary values.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces azeros from mikpckd.f
 */
VEXTERNC void Vazeros(
        int *nx,  ///< The size of the x dimension of the 3d matrix
        int *ny,  ///< The size of the x dimension of the 3d matrix
        int *nz,  ///< The size of the x dimension of the 3d matrix
        double *x ///< The matrix to zero out
        );

/** @brief   Initialize a grid function to have a certain boundary value,
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces fboundPMG from mikpckd.f
 */
VEXTERNC void VfboundPMG(
        int *ibound, ///< @todo: Doc
        int *nx,     ///< @todo: Doc
        int *ny,     ///< @todo: Doc
        int *nz,     ///< @todo: Doc
        double *x,   ///< @todo: Doc
        double *gxc, ///< @todo: Doc
        double *gyc, ///< @todo: Doc
        double *gzc  ///< @todo: Doc
        );

/** @brief   Initialize a grid function to have a zero boundary value
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces fboundPMG00 from mikpckd.f
 */
VEXTERNC void VfboundPMG00(
        int    *nx, ///< The size of the x dimension of the 3d matrix
        int    *ny, ///< The size of the y dimension of the 3d matrix
        int    *nz, ///< The size of the z dimension of the 3d matrix
        double *x   ///< The 3d matrix to initialize
        );

/** @brief   Fill grid function with random values, including boundary values.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces axrand from mikpckd.f
 */
VEXTERNC void Vaxrand(
        int    *nx, ///< The size of the x dimension of the 3d matrix
        int    *ny, ///< The size of the y dimension of the 3d matrix
        int    *nz, ///< The size of the z dimension of the 3d matrix
        double *x   ///< The 3d matrix to fill
        );

/** @brief   Scale operation for a grid function with boundary values.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces xscal from mikpckd.f
 */
VEXTERNC void Vxscal(
        int    *nx,  ///< The size of the x dimension of the 3d matrix
        int    *ny,  ///< The size of the y dimension of the 3d matrix
        int    *nz,  ///< The size of the z dimension of the 3d matrix
        double *fac, ///< The scaling factor
        double *x    ///< The 3d matrix to scale
        );

/** @brief
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces prtmatd from mikpckd.f
 */
VEXTERNC void Vprtmatd(
        int    *nx,  ///< The size of the x dimension of the 3d matrix
        int    *ny,  ///< The size of the y dimension of the 3d matrix
        int    *nz,  ///< The size of the z dimension of the 3d matrix
        int    *ipc, ///< Integer parameters
        double *rpc, ///< Double parameters
        double *ac   ///< @todo  Document
        );

VEXTERNC void Vprtmatd7(
        int *nx,     ///< The size of the x dimension of the 3d matrix
        int *ny,     ///< The size of the y dimension of the 3d matrix
        int *nz,     ///< The size of the z dimension of the 3d matrix
        int *ipc,    ///< Integer parameters
        double *rpc, ///< Double parameters
        double *oC,  ///< @todo  Document
        double *oE,  ///< @todo  Document
        double *oN,  ///< @todo  Document
        double *uC   ///< @todo  Document
        );

VEXTERNC void Vprtmatd27(
        int    *nx,  ///< The size of the x dimension of the 3d matrix
        int    *ny,  ///< The size of the y dimension of the 3d matrix
        int    *nz,  ///< The size of the z dimension of the 3d matrix
        int    *ipc, ///< Integer parameters
        double *rpc, ///< Double parameters
        double *oC,  ///< @todo  Document
        double *oE,  ///< @todo  Document
        double *oN,  ///< @todo  Document
        double *uC,  ///< @todo  Document
        double *oNE, ///< @todo  Document
        double *oNW, ///< @todo  Document
        double *uE,  ///< @todo  Document
        double *uW,  ///< @todo  Document
        double *uN,  ///< @todo  Document
        double *uS,  ///< @todo  Document
        double *uNE, ///< @todo  Document
        double *uNW, ///< @todo  Document
        double *uSE, ///< @todo  Document
        double *uSW  ///< @todo  Document
        );

VEXTERNC void Vlinesearch(
        int    *nx,    ///< The size of the x dimension of the 3d matrix
        int    *ny,    ///< The size of the y dimension of the 3d matrix
        int    *nz,    ///< The size of the z dimension of the 3d matrix
        double *alpha, ///< @todo  Document
        int    *ipc,   ///< Integer parameters
        double *rpc,   ///< Double parameters
        double *ac,    ///< @todo  Document
        double *cc,    ///< @todo  Document
        double *fc,    ///< @todo  Document
        double *p,     ///< @todo  Document
        double *x,     ///< @todo  Document
        double *r,     ///< @todo  Document
        double *ap,    ///< @todo  Document
        double *zk,    ///< @todo  Document
        double *zkp1   ///< @todo  Document
        );

#endif /* MIKPCKD_H_ */
