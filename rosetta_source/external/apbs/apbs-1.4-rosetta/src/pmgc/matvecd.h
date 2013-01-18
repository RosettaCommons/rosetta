/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief  Matrix-vector multiplication routines
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

#ifndef _MATVECD_H_
#define _MATVECD_H_

#include "apbscfg.h"

#include "generic/vhal.h"
#include "generic/vmatrix.h"
#include "pmgc/mikpckd.h"
#include "pmgc/mypdec.h"

/** @brief   Break the matrix data-structure into diagonals and
 *           then call the matrix-vector routine.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces matvec from matvecd.f
 */
VEXTERNC void Vmatvec(
        int    *nx,  ///< @todo:  Doc
        int    *ny,  ///< @todo:  Doc
        int    *nz,  ///< @todo:  Doc
        int    *ipc, ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *ac,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *y    ///< @todo:  Doc
        );

VEXTERNC void Vmatvec7(
        int    *nx,  ///< @todo:  Doc
        int    *ny,  ///< @todo:  Doc
        int    *nz,  ///< @todo:  Doc
        int    *ipc, ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *ac,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *y    ///< @todo:  Doc
        );

VEXTERNC void Vmatvec7_1s(
        int    *nx,  ///< @todo:  Doc
        int    *ny,  ///< @todo:  Doc
        int    *nz,  ///< @todo:  Doc
        int    *ipc, ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *oC,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *oE,  ///< @todo:  Doc
        double *oN,  ///< @todo:  Doc
        double *uC,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *y    ///< @todo:  Doc
        );



VEXTERNC void Vmatvec27(
        int    *nx,  ///< @todo:  Doc
        int    *ny,  ///< @todo:  Doc
        int    *nz,  ///< @todo:  Doc
        int    *ipc, ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *ac,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *y    ///< @todo:  Doc
        );

VEXTERNC void Vmatvec27_1s(
        int    *nx,  ///< @todo:  Doc
        int    *ny,  ///< @todo:  Doc
        int    *nz,  ///< @todo:  Doc
        int    *ipc, ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *oC,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *oE,  ///< @todo:  Doc
        double *oN,  ///< @todo:  Doc
        double *uC,  ///< @todo:  Doc
        double *oNE, ///< @todo:  Doc
        double *oNW, ///< @todo:  Doc
        double *uE,  ///< @todo:  Doc
        double *uW,  ///< @todo:  Doc
        double *uN,  ///< @todo:  Doc
        double *uS,  ///< @todo:  Doc
        double *uNE, ///< @todo:  Doc
        double *uNW, ///< @todo:  Doc
        double *uSE, ///< @todo:  Doc
        double *uSW, ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *y    ///< @todo:  Doc
        );


/** @brief   Break the matrix data-structure into diagonals and
 *           then call the matrix-vector routine.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces nmatvec from matvecd.f
 */
VEXTERNC void Vnmatvec(
        int *nx,     ///< @todo  Document
        int *ny,     ///< @todo  Document
        int *nz,     ///< @todo  Document
        int *ipc,    ///< @todo  Document
        double *rpc, ///< @todo  Document
        double *ac,  ///< @todo  Document
        double *cc,  ///< @todo  Document
        double *x,   ///< @todo  Document
        double *y,   ///< @todo  Document
        double *w1   ///< @todo  Document
        );

VEXTERNC void Vnmatvec7(
        int    *nx,  ///< @todo  Document
        int    *ny,  ///< @todo  Document
        int    *nz,  ///< @todo  Document
        int    *ipc, ///< @todo  Document
        double *rpc, ///< @todo  Document
        double  *ac, ///< @todo  Document
        double  *cc, ///< @todo  Document
        double   *x, ///< @todo  Document
        double   *y, ///< @todo  Document
        double  *w1  ///< @todo  Document
        );

VEXTERNC void Vnmatvecd7_1s(
        int     *nx, ///< @todo  Document
        int     *ny, ///< @todo  Document
        int     *nz, ///< @todo  Document
        int    *ipc, ///< @todo  Document
        double *rpc, ///< @todo  Document
        double  *oC, ///< @todo  Document
        double  *cc, ///< @todo  Document
        double  *oE, ///< @todo  Document
        double  *oN, ///< @todo  Document
        double  *uC, ///< @todo  Document
        double   *x, ///< @todo  Document
        double   *y, ///< @todo  Document
        double  *w1  ///< @todo  Document
        );

VEXTERNC void Vnmatvec27(
        int    *nx,  ///< @todo  Document
        int    *ny,  ///< @todo  Document
        int    *nz,  ///< @todo  Document
        int    *ipc, ///< @todo  Document
        double *rpc, ///< @todo  Document
        double  *ac, ///< @todo  Document
        double  *cc, ///< @todo  Document
        double   *x, ///< @todo  Document
        double   *y, ///< @todo  Document
        double  *w1  ///< @todo  Document
        );

VEXTERNC void Vnmatvecd27_1s(
        int     *nx, ///< @todo  Document
        int     *ny, ///< @todo  Document
        int     *nz, ///< @todo  Document
        int    *ipc, ///< @todo  Document
        double *rpc, ///< @todo  Document
        double  *oC, ///< @todo  Document
        double  *cc, ///< @todo  Document
        double  *oE, ///< @todo  Document
        double  *oN, ///< @todo  Document
        double  *uC, ///< @todo  Document
        double *oNE, ///< @todo  Document
        double *oNW, ///< @todo  Document
        double  *uE, ///< @todo  Document
        double  *uW, ///< @todo  Document
        double  *uN, ///< @todo  Document
        double  *uS, ///< @todo  Document
        double *uNE, ///< @todo  Document
        double *uNW, ///< @todo  Document
        double *uSE, ///< @todo  Document
        double *uSW, ///< @todo  Document
        double   *x, ///< @todo  Document
        double   *y, ///< @todo  Document
        double  *w1  ///< @todo  Document
        );


/** @brief   Break the matrix data-structure into diagonals and
 *           then call the residual routine.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces mresid from matvecd.f
 */
VEXTERNC void Vmresid(
        int    *nx,  ///< @todo:  Doc
        int    *ny,  ///< @todo:  Doc
        int    *nz,  ///< @todo:  Doc
        int    *ipc, ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *ac,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *fc,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *r    ///< @todo:  Doc
        );

VEXTERNC void Vmresid7(
        int    *nx,  ///< @todo:  Doc
        int    *ny,  ///< @todo:  Doc
        int    *nz,  ///< @todo:  Doc
        int    *ipc, ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *ac,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *fc,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *r    ///< @todo:  Doc
        );

VEXTERNC void Vmresid7_1s(
        int *nx,     ///< @todo:  Doc
        int *ny,     ///< @todo:  Doc
        int *nz,     ///< @todo:  Doc
        int *ipc,    ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *oC,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *fc,  ///< @todo:  Doc
        double *oE,  ///< @todo:  Doc
        double *oN,  ///< @todo:  Doc
        double *uC,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *r    ///< @todo:  Doc
        );

VEXTERNC void Vmresid27(
        int    *nx,  ///< @todo:  Doc
        int    *ny,  ///< @todo:  Doc
        int    *nz,  ///< @todo:  Doc
        int    *ipc, ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *ac,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *fc,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *r    ///< @todo:  Doc
        );

VEXTERNC void Vmresid27_1s(
        int *nx,     ///< @todo:  Doc
        int *ny,     ///< @todo:  Doc
        int *nz,     ///< @todo:  Doc
        int *ipc,    ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *oC,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *fc,  ///< @todo:  Doc
        double *oE,  ///< @todo:  Doc
        double *oN,  ///< @todo:  Doc
        double *uC,  ///< @todo:  Doc
        double *oNE, ///< @todo:  Doc
        double *oNW, ///< @todo:  Doc
        double *uE,  ///< @todo:  Doc
        double *uW,  ///< @todo:  Doc
        double *uN,  ///< @todo:  Doc
        double *uS,  ///< @todo:  Doc
        double *uNE, ///< @todo:  Doc
        double *uNW, ///< @todo:  Doc
        double *uSE, ///< @todo:  Doc
        double *uSW, ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *r    ///< @todo:  Doc
        );



/** @brief   Break the matrix data-structure into diagonals and
 *           then call the residual routine.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces nmresid from matvecd.f
 */
VEXTERNC void Vnmresid(
        int    *nx,  ///< @todo:  Doc
        int    *ny,  ///< @todo:  Doc
        int    *nz,  ///< @todo:  Doc
        int    *ipc, ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *ac,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *fc,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *r,   ///< @todo:  Doc
        double *w1   ///< @todo:  Doc
        );

VEXTERNC void Vnmresid7(
        int *nx,     ///< @todo:  Doc
        int *ny,     ///< @todo:  Doc
        int *nz,     ///< @todo:  Doc
        int *ipc,    ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *ac,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *fc,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *r,   ///< @todo:  Doc
        double *w1   ///< @todo:  Doc
        );

VEXTERNC void Vnmresid7_1s(
        int *nx,     ///< @todo:  Doc
        int *ny,     ///< @todo:  Doc
        int *nz,     ///< @todo:  Doc
        int *ipc,    ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *oC,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *fc,  ///< @todo:  Doc
        double *oE,  ///< @todo:  Doc
        double *oN,  ///< @todo:  Doc
        double *uC,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *r,   ///< @todo:  Doc
        double *w1   ///< @todo:  Doc
        );

VEXTERNC void Vnmresid27(
        int *nx,     ///< @todo:  Doc
        int *ny,     ///< @todo:  Doc
        int *nz,     ///< @todo:  Doc
        int *ipc,    ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *ac,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *fc,  ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *r,   ///< @todo:  Doc
        double *w1   ///< @todo:  Doc
        );

VEXTERNC void Vnmresid27_1s(
        int *nx,     ///< @todo:  Doc
        int *ny,     ///< @todo:  Doc
        int *nz,     ///< @todo:  Doc
        int *ipc,    ///< @todo:  Doc
        double *rpc, ///< @todo:  Doc
        double *oC,  ///< @todo:  Doc
        double *cc,  ///< @todo:  Doc
        double *fc,  ///< @todo:  Doc
        double *oE,  ///< @todo:  Doc
        double *oN,  ///< @todo:  Doc
        double *uC,  ///< @todo:  Doc
        double *oNE, ///< @todo:  Doc
        double *oNW, ///< @todo:  Doc
        double *uE,  ///< @todo:  Doc
        double *uW,  ///< @todo:  Doc
        double *uN,  ///< @todo:  Doc
        double *uS,  ///< @todo:  Doc
        double *uNE, ///< @todo:  Doc
        double *uNW, ///< @todo:  Doc
        double *uSE, ///< @todo:  Doc
        double *uSW, ///< @todo:  Doc
        double *x,   ///< @todo:  Doc
        double *r,   ///< @todo:  Doc
        double *w1   ///< @todo:  Doc
        );



/** @brief   Apply the restriction operator
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces restrc from matvecd.f
 */
VEXTERNC void Vrestrc(
        int *nxf,     ///< @todo:  Doc
        int *nyf,     ///< @todo:  Doc
        int *nzf,     ///< @todo:  Doc
        int *nxc,     ///< @todo:  Doc
        int *nyc,     ///< @todo:  Doc
        int *nzc,     ///< @todo:  Doc
        double *xin,  ///< @todo:  Doc
        double *xout, ///< @todo:  Doc
        double *pc    ///< @todo:  Doc
        );

VEXTERNC void Vrestrc2(
        int    *nxf,///< @todo:  Doc
        int    *nyf,///< @todo:  Doc
        int    *nzf,///< @todo:  Doc
        int    *nxc,///< @todo:  Doc
        int    *nyc,///< @todo:  Doc
        int    *nzc,///< @todo:  Doc
        double *xin,///< @todo:  Doc
        double *xout,///< @todo:  Doc
        double *oPC,///< @todo:  Doc
        double *oPN,///< @todo:  Doc
        double *oPS,///< @todo:  Doc
        double *oPE,///< @todo:  Doc
        double *oPW,///< @todo:  Doc
        double *oPNE,///< @todo:  Doc
        double *oPNW,///< @todo:  Doc
        double *oPSE,///< @todo:  Doc
        double *oPSW,///< @todo:  Doc
        double *uPC,///< @todo:  Doc
        double *uPN,///< @todo:  Doc
        double *uPS,///< @todo:  Doc
        double *uPE,///< @todo:  Doc
        double *uPW,///< @todo:  Doc
        double *uPNE,///< @todo:  Doc
        double *uPNW,///< @todo:  Doc
        double *uPSE,///< @todo:  Doc
        double *uPSW,///< @todo:  Doc
        double *dPC,///< @todo:  Doc
        double *dPN,///< @todo:  Doc
        double *dPS,///< @todo:  Doc
        double *dPE,///< @todo:  Doc
        double *dPW,///< @todo:  Doc
        double *dPNE,///< @todo:  Doc
        double *dPNW,///< @todo:  Doc
        double *dPSE,///< @todo:  Doc
        double *dPSW///< @todo:  Doc
        );

/** @brief   Apply the prolongation operator
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces interpPMG from matvecd.f
 */
VEXTERNC void VinterpPMG(
        int    *nxc,  ///< @todo:  Doc
        int    *nyc,  ///< @todo:  Doc
        int    *nzc,  ///< @todo:  Doc
        int    *nxf,  ///< @todo:  Doc
        int    *nyf,  ///< @todo:  Doc
        int    *nzf,  ///< @todo:  Doc
        double *xin,  ///< @todo:  Doc
        double *xout, ///< @todo:  Doc
        double *pc    ///< @todo:  Doc
        );

VEXTERNC void VinterpPMG2(
        int     *nxc, ///< @todo:  Doc
        int     *nyc, ///< @todo:  Doc
        int     *nzc, ///< @todo:  Doc
        int     *nxf, ///< @todo:  Doc
        int     *nyf, ///< @todo:  Doc
        int     *nzf, ///< @todo:  Doc
        double  *xin, ///< @todo:  Doc
        double *xout, ///< @todo:  Doc
        double  *oPC, ///< @todo:  Doc
        double  *oPN, ///< @todo:  Doc
        double  *oPS, ///< @todo:  Doc
        double  *oPE, ///< @todo:  Doc
        double  *oPW, ///< @todo:  Doc
        double *oPNE, ///< @todo:  Doc
        double *oPNW, ///< @todo:  Doc
        double *oPSE, ///< @todo:  Doc
        double *oPSW, ///< @todo:  Doc
        double  *uPC, ///< @todo:  Doc
        double  *uPN, ///< @todo:  Doc
        double  *uPS, ///< @todo:  Doc
        double  *uPE, ///< @todo:  Doc
        double  *uPW, ///< @todo:  Doc
        double *uPNE, ///< @todo:  Doc
        double *uPNW, ///< @todo:  Doc
        double *uPSE, ///< @todo:  Doc
        double *uPSW, ///< @todo:  Doc
        double  *dPC, ///< @todo:  Doc
        double  *dPN, ///< @todo:  Doc
        double  *dPS, ///< @todo:  Doc
        double  *dPE, ///< @todo:  Doc
        double  *dPW, ///< @todo:  Doc
        double *dPNE, ///< @todo:  Doc
        double *dPNW, ///< @todo:  Doc
        double *dPSE, ///< @todo:  Doc
        double *dPSW  ///< @todo:  Doc
        );

/** @brief   Simple injection of a fine grid function into coarse grid.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces extrac from matvecd.f
 */
VEXTERNC void Vextrac(
        int *nxf,    ///< @todo:  Doc
        int *nyf,    ///< @todo:  Doc
        int *nzf,    ///< @todo:  Doc
        int *nxc,    ///< @todo:  Doc
        int *ny,     ///< @todo:  Doc
        int *nzc,    ///< @todo:  Doc
        double *xin, ///< @todo:  Doc
        double *xout ///< @todo:  Doc
        );

#endif /* _MATVECD_H_ */
