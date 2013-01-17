/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief MG helper functions
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

#ifndef _MGCSD_H_
#define _MGCSD_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vmatrix.h"
#include "pmgc/mgsubd.h"
#include "pmgc/mikpckd.h"
#include "pmgc/smoothd.h"
#include "pmgc/gsd.h"
#include "pmgc/matvecd.h"
#include "pmgc/mlinpckd.h"

/** @brief   Screaming linear multilevel method.
 *
 *    algorithm:  linear multigrid iteration (cs)
 *
 *    multigrid v-cycle solver.
 *
 *    input:
 *       (1) fine and coarse grid discrete linear operators: L_h, L_H
 *       (2) fine grid source function: f_h
 *       (3) fine grid approximate solution: u_h
 *
 *    output:
 *       (1) fine grid improved solution: u_h
 *
 *    the two-grid algorithm is:
 *       (1) pre-smooth:               u1_h = smooth(L_h,f_h,u_h)
 *       (2) restrict defect:          d_H  = r(L_h(u1_h) - f_h)
 *       (3) solve for correction:     c_H  = L_H^{-1}(d_H)
 *       (4) prolongate and correct:   u2_h = u1_h - p(c_H)
 *       (5) post-smooth:              u_h  = smooth(L_h,f_h,u2_h)
 *
 *    (of course, c_H is determined with another two-grid algorithm)
 *
 *    implementation notes:
 *       (0) "u1_h" must be kept on each level until "c_H" is computed,
 *           and then both are used to compute "u2_h".
 *       (1) "u_h" (and then "u1_h") on all levels is stored in the "x" array.
 *       (2) "d_H" is identically "f_h" for f_h on the next coarser grid.
 *       (3) "c_h" is identically "u_h" for u_h on the next coarser grid.
 *       (4) "d_H" is stored in the "r" array (must be kept for post-smooth).
 *       (5) "f_h" is stored in the "fc" array.
 *       (6) "L_h" on all levels is stored in the "ac" array.
 *       (7) signs may be reveresed; i.e., residual is used in place
 *           of the defect in places, etc.
 *
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces mvcs from mgcsd.f
 */
VEXTERNC void Vmvcs(
        int    *nx,        ///< @todo: doc
        int    *ny,        ///< @todo: doc
        int    *nz,        ///< @todo: doc
        double *x,         ///< @todo: doc
        int    *iz,        ///< @todo: doc
        double *w0,        ///< @todo: doc
        double *w1,        ///< @todo: doc
        double *w2,        ///< @todo: doc
        double *w3,        ///< @todo: doc
        int    *istop,     ///< @todo: doc
        int    *itmax,     ///< @todo: doc
        int    *iters,     ///< @todo: doc
        int    *ierror,    ///< @todo: doc
        int    *nlev,      ///< @todo: doc
        int    *ilev,      ///< @todo: doc
        int    *nlev_real, ///< @todo: doc
        int    *mgsolv,    ///< @todo: doc
        int    *iok,       ///< @todo: doc
        int    *iinfo,     ///< @todo: doc
        double *epsiln,    ///< @todo: doc
        double *errtol,    ///< @todo: doc
        double *omega,     ///< @todo: doc
        int    *nu1,       ///< @todo: doc
        int    *nu2,       ///< @todo: doc
        int    *mgsmoo,    ///< @todo: doc
        int    *ipc,       ///< @todo: doc
        double *rpc,       ///< @todo: doc
        double *pc,        ///< @todo: doc
        double *ac,        ///< @todo: doc
        double *cc,        ///< @todo: doc
        double *fc,        ///< @todo: doc
        double *tru        ///< @todo: doc
        );

#endif /* _MGCSD_H_ */
