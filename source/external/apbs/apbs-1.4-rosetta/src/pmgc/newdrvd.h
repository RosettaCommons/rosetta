/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief   Driver for the Newton Solver
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

#ifndef _NEWDRVD_H_
#define _NEWDRVD_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "pmgc/mgsubd.h"
#include "pmgc/mikpckd.h"
#include "pmgc/newtond.h"
#include "pmgc/mgdrvd.h"

/** @brief   Driver for a screaming inexact-newton-multilevel solver.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces newdriv from newdrvd.f
 */
VEXTERNC void Vnewdriv(
        int    *iparm, ///< @todo:  Doc
        double *rparm, ///< @todo:  Doc
        int    *iwork, ///< @todo:  Doc
        double *rwork, ///< @todo:  Doc
        double *u,     ///< @todo:  Doc
        double *xf,    ///< @todo:  Doc
        double *yf,    ///< @todo:  Doc
        double *zf,    ///< @todo:  Doc
        double *gxcf,  ///< @todo:  Doc
        double *gycf,  ///< @todo:  Doc
        double *gzcf,  ///< @todo:  Doc
        double *a1cf,  ///< @todo:  Doc
        double *a2cf,  ///< @todo:  Doc
        double *a3cf,  ///< @todo:  Doc
        double *ccf,   ///< @todo:  Doc
        double *fcf,   ///< @todo:  Doc
        double *tcf    ///< @todo:  Doc
        );

/** @brief   Solves using Newton's Method
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *    This routine uses a newton's method, combined with a linear
 *    multigrid iteration, to solve the following three-dimensional,
 *    2nd order elliptic partial differential equation:
 *
 *         lu = f, u in omega
 *          u = g, u on boundary of omega
 *    where
 *
 *         omega = [xmin,xmax]x[ymin,ymax]x[zmin,zmax]
 *
 *    the multigrid code requires the operator in the form:
 *
 *         - \nabla \cdot (a \nabla u) + c(u) = f
 *
 *    with
 *
 *        a(x,y,z),f(x,y,z), scalar functions (possibly discontinuous)
 *        on omega.  (discontinuities must be along fine grid lines).
 *        boundary function g(x,y,z) is smooth on boundary of omega.
 *
 *        the function c(u) is a possibly nonlinear function of the
 *        unknown u, and varies (possibly discontinuously) with the
 *        spatial position also.
 *
 *   User inputs:
 *
 *    the user must provide the coefficients of the differential
 *    operator, some initial parameter settings in an integer and a
 *    real parameter array, and various work arrays.
 *
 *  @note    Replaces newdriv2 from newdrvd.f
 */
VEXTERNC void Vnewdriv2(
        int    *iparm, ///< @todo:  Doc
        double *rparm, ///< @todo:  Doc
        int    *nx,    ///< @todo:  Doc
        int    *ny,    ///< @todo:  Doc
        int    *nz,    ///< @todo:  Doc
        double *u,     ///< @todo:  Doc
        int    *iz,    ///< @todo:  Doc
        double *w1,    ///< @todo:  Doc
        double *w2,    ///< @todo:  Doc
        int    *ipc,   ///< @todo:  Doc
        double *rpc,   ///< @todo:  Doc
        double *pc,    ///< @todo:  Doc
        double *ac,    ///< @todo:  Doc
        double *cc,    ///< @todo:  Doc
        double *fc,    ///< @todo:  Doc
        double *xf,    ///< @todo:  Doc
        double *yf,    ///< @todo:  Doc
        double *zf,    ///< @todo:  Doc
        double *gxcf,  ///< @todo:  Doc
        double *gycf,  ///< @todo:  Doc
        double *gzcf,  ///< @todo:  Doc
        double *a1cf,  ///< @todo:  Doc
        double *a2cf,  ///< @todo:  Doc
        double *a3cf,  ///< @todo:  Doc
        double *ccf,   ///< @todo:  Doc
        double *fcf,   ///< @todo:  Doc
        double *tcf    ///< @todo:  Doc
        );

#endif /* _NEWDRVD_H_ */
