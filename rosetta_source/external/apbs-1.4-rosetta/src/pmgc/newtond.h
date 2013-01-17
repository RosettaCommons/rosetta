/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief Driver routines for the Newton method
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

#ifndef _NEWTOND_H_
#define _NEWTOND_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vmatrix.h"
#include "pmgc/matvecd.h"
#include "pmgc/mikpckd.h"
#include "pmgc/mgcsd.h"
#include "pmgc/mgsubd.h"
#include "pmgc/powerd.h"

/** @brief   Nested iteration for an inexact-newton-multilevel method.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces fnewton from newtond.f
 */
VPUBLIC void Vfnewton(
        int *nx,        ///< @todo: Doc
        int *ny,        ///< @todo: Doc
        int *nz,        ///< @todo: Doc
        double *x,      ///< @todo: Doc
        int *iz,        ///< @todo: Doc
        double *w0,     ///< @todo: Doc
        double *w1,     ///< @todo: Doc
        double *w2,     ///< @todo: Doc
        double *w3,     ///< @todo: Doc
        int *istop,     ///< @todo: Doc
        int *itmax,     ///< @todo: Doc
        int *iters,     ///< @todo: Doc
        int *ierror,    ///< @todo: Doc
        int *nlev,      ///< @todo: Doc
        int *ilev,      ///< @todo: Doc
        int *nlev_real, ///< @todo: Doc
        int *mgsolv,    ///< @todo: Doc
        int *iok,       ///< @todo: Doc
        int *iinfo,     ///< @todo: Doc
        double *epsiln, ///< @todo: Doc
        double *errtol, ///< @todo: Doc
        double *omega,  ///< @todo: Doc
        int *nu1,       ///< @todo: Doc
        int *nu2,       ///< @todo: Doc
        int *mgsmoo,    ///< @todo: Doc
        double *cprime, ///< @todo: Doc
        double *rhs,    ///< @todo: Doc
        double *xtmp,   ///< @todo: Doc
        int *ipc,       ///< @todo: Doc
        double *rpc,    ///< @todo: Doc
        double *pc,     ///< @todo: Doc
        double *ac,     ///< @todo: Doc
        double *cc,     ///< @todo: Doc
        double *fc,     ///< @todo: Doc
        double *tru     ///< @todo: Doc
        );

/** @brief   Inexact-newton-multilevel method.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces newton from newtond.f
 */
VEXTERNC void Vnewton(
        int *nx,        ///< @todo: Doc
        int *ny,        ///< @todo: Doc
        int *nz,        ///< @todo: Doc
        double *x,      ///< @todo: Doc
        int *iz,        ///< @todo: Doc
        double *w0,     ///< @todo: Doc
        double *w1,     ///< @todo: Doc
        double *w2,     ///< @todo: Doc
        double *w3,     ///< @todo: Doc
        int *istop,     ///< @todo: Doc
        int *itmax,     ///< @todo: Doc
        int *iters,     ///< @todo: Doc
        int *ierror,    ///< @todo: Doc
        int *nlev,      ///< @todo: Doc
        int *ilev,      ///< @todo: Doc
        int *nlev_real, ///< @todo: Doc
        int *mgsolv,    ///< @todo: Doc
        int *iok,       ///< @todo: Doc
        int *iinfo,     ///< @todo: Doc
        double *epsiln, ///< @todo: Doc
        double *errtol, ///< @todo: Doc
        double *omega,  ///< @todo: Doc
        int *nu1,       ///< @todo: Doc
        int *nu2,       ///< @todo: Doc
        int *mgsmoo,    ///< @todo: Doc
        double *cprime, ///< @todo: Doc
        double *rhs,    ///< @todo: Doc
        double *xtmp,   ///< @todo: Doc
        int *ipc,       ///< @todo: Doc
        double *rpc,    ///< @todo: Doc
        double *pc,     ///< @todo: Doc
        double *ac,     ///< @todo: Doc
        double *cc,     ///< @todo: Doc
        double *fc,     ///< @todo: Doc
        double *tru     ///< @todo: Doc
        );


/** @brief   Form the jacobian system.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces getjac from newtond.f
 */
VEXTERNC void Vgetjac(
        int *nx,        ///< @todo: Doc
        int *ny,        ///< @todo: Doc
        int *nz,        ///< @todo: Doc
        int *nlev_real, ///< @todo: Doc
        int *iz,        ///< @todo: Doc
        int *lev,       ///< @todo: Doc
        int *ipkey,     ///< @todo: Doc
        double *x,      ///< @todo: Doc
        double *r,      ///< @todo: Doc
        double *cprime, ///< @todo: Doc
        double *rhs,    ///< @todo: Doc
        double *cc,     ///< @todo: Doc
        double *pc      ///< @todo: Doc
        );

#endif /* _NEWTOND_H_ */
