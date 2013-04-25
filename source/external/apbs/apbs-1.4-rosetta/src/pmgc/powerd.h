/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief Power methods for eigenvalue estimation
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

#ifndef _POWERD_H_
#define _POWERD_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vmatrix.h"
#include "pmgc/mikpckd.h"
#include "pmgc/mgcsd.h"

/** @brief  Standard power method for maximum eigenvalue estimation of a matrix
c*
c*  @note   To test, note that the 3d laplacean has min/max eigenvalues:
c*
c*       lambda_min = 6 - 2*dcos(pi/(nx-1))
c*                      - 2*dcos(pi/(ny-1))
c*                      - 2*dcos(pi/(nz-1))
c*
c*       lambda_max = 6 - 2*dcos((nx-2)*pi/(nx-1))
c*                      - 2*dcos((ny-2)*pi/(ny-1))
c*                      - 2*dcos((nz-2)*pi/(nz-1))
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces power from powerd.f
 *  @note    Vpower is yet untested as a call stack including it hasn't been found
 */
VEXTERNC void Vpower(
        int *nx,              ///< @todo  Document
        int *ny,           	  ///< @todo  Document
        int *nz,              ///< @todo  Document
        int *iz,              ///< @todo  Document
        int *ilev,            ///< @todo  Document
        int *ipc,             ///< @todo  Document
        double *rpc,          ///< @todo  Document
        double *ac,           ///< @todo  Document
        double *cc,           ///< @todo  Document
        double *w1,           ///< @todo  Document
        double *w2,           ///< @todo  Document
        double *w3,           ///< @todo  Document
        double *w4,           ///< @todo  Document
        double *eigmax,       ///< @todo  Document
        double *eigmax_model, ///< @todo  Document
        double *tol,          ///< @todo  Document
        int *itmax,           ///< @todo  Document
        int *iters,           ///< @todo  Document
        int *iinfo            ///< @todo  Document
        );



/** @brief  Standard inverse power method for minimum eigenvalue estimation
 *
 *  @note   To test, note that the 3d laplacean has min/max eigenvalues:
 *
 *       lambda_min = 6 - 2*dcos(pi/(nx-1))
 *                      - 2*dcos(pi/(ny-1))
 *                      - 2*dcos(pi/(nz-1))
 *
 *       lambda_max = 6 - 2*dcos((nx-2)*pi/(nx-1))
 *                      - 2*dcos((ny-2)*pi/(ny-1))
 *                      - 2*dcos((nz-2)*pi/(nz-1))
 *
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces ipower from powerd.f
 */
VEXTERNC void Vipower(
        int *nx,              ///< @todo  Document
        int *ny,              ///< @todo  Document
        int *nz,              ///< @todo  Document
        double *u,            ///< @todo  Document
        int *iz,              ///< @todo  Document
        double *w0,           ///< @todo  Document
        double *w1,           ///< @todo  Document
        double *w2,           ///< @todo  Document
        double *w3,           ///< @todo  Document
        double *w4,           ///< @todo  Document
        double *eigmin,       ///< @todo  Document
        double *eigmin_model, ///< @todo  Document
        double *tol,          ///< @todo  Document
        int *itmax,           ///< @todo  Document
        int *iters,           ///< @todo  Document
        int *nlev,            ///< @todo  Document
        int *ilev,            ///< @todo  Document
        int *nlev_real,       ///< @todo  Document
        int *mgsolv,          ///< @todo  Document
        int *iok,             ///< @todo  Document
        int *iinfo,           ///< @todo  Document
        double *epsiln,       ///< @todo  Document
        double *errtol,       ///< @todo  Document
        double *omega,        ///< @todo  Document
        int *nu1,             ///< @todo  Document
        int *nu2,             ///< @todo  Document
        int *mgsmoo,          ///< @todo  Document
        int *ipc,             ///< @todo  Document
        double *rpc,          ///< @todo  Document
        double *pc,           ///< @todo  Document
        double *ac,           ///< @todo  Document
        double *cc,           ///< @todo  Document
        double *tru           ///< @todo  Document
        );



VEXTERNC void Vmpower(
        int    *nx,        ///< @todo  Document
        int    *ny,        ///< @todo  Document
        int    *nz,        ///< @todo  Document
        double *u,         ///< @todo  Document
        int    *iz,        ///< @todo  Document
        double *w0,        ///< @todo  Document
        double *w1,        ///< @todo  Document
        double *w2,        ///< @todo  Document
        double *w3,        ///< @todo  Document
        double *w4,        ///< @todo  Document
        double *eigmax,    ///< @todo  Document
        double *tol,       ///< @todo  Document
        int    *itmax,     ///< @todo  Document
        int    *iters,     ///< @todo  Document
        int    *nlev,      ///< @todo  Document
        int    *ilev,      ///< @todo  Document
        int    *nlev_real, ///< @todo  Document
        int    *mgsolv,    ///< @todo  Document
        int    *iok,       ///< @todo  Document
        int    *iinfo,     ///< @todo  Document
        double *epsiln,    ///< @todo  Document
        double *errtol,    ///< @todo  Document
        double *omega,     ///< @todo  Document
        int    *nu1,       ///< @todo  Document
        int    *nu2,       ///< @todo  Document
        int    *mgsmoo,    ///< @todo  Document
        int    *ipc,       ///< @todo  Document
        double *rpc,       ///< @todo  Document
        double *pc,        ///< @todo  Document
        double *ac,        ///< @todo  Document
        double *cc,        ///< @todo  Document
        double *fc,        ///< @todo  Document
        double *tru        ///< @todo  Document
        );

#endif /* _POWERD_H_ */
