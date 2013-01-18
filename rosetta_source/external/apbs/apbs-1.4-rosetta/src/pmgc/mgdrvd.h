/** @defgroup PMGC C translation of Holst group PMG code
 *  @brief C translation of Holst group PMG code
 */

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

#ifndef _MGDRVD_H_
#define _MGDRVD_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vmatrix.h"
#include "pmgc/mgsubd.h"
#include "pmgc/mgcsd.h"
#include "pmgc/powerd.h"
#include "pmgc/mgfasd.h"

/** @brief   Multilevel solver driver
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  Replaces mgdriv from mgdrvd.f
 */
VEXTERNC void Vmgdriv(
        int* iparm,    ///< @todo: Doc
        double* rparm, ///< @todo: Doc
        int* iwork,    ///< @todo: Doc
        double* rwork, ///< @todo: Doc
        double* u,     ///< @todo: Doc
        double* xf,    ///< @todo: Doc
        double* yf,    ///< @todo: Doc
        double* zf,    ///< @todo: Doc
        double* gxcf,  ///< @todo: Doc
        double* gycf,  ///< @todo: Doc
        double* gzcf,  ///< @todo: Doc
        double* a1cf,  ///< @todo: Doc
        double* a2cf,  ///< @todo: Doc
        double* a3cf,  ///< @todo: Doc
        double* ccf,   ///< @todo: Doc
        double* fcf,   ///< @todo: Doc
        double* tcf    ///< @todo: Doc
        );

/** @brief   Solves the pde using the multi-grid method
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *	Replaces mgdriv2 from mgdrvd.f
 *
 *  This routine uses a multigrid method to solve the following
 *  three-dimensional, 2nd order elliptic partial differential
 *  equation:
 *
 *       lu = f, u in omega
 *        u = g, u on boundary of omega
 *  where
 *
 *       omega = [xmin,xmax]x[ymin,ymax]x[zmin,zmax]
 *
 *  the multigrid code requires the operator in the form:
 *
 *       - \nabla \cdot (a \nabla u) + c(u) = f
 *
 *  with
 *
 *      a(x,y,z),f(x,y,z), scalar functions (possibly discontinuous)
 *      on omega.  (discontinuities must be along fine grid lines).
 *      boundary function g(x,y,z) is smooth on boundary of omega.
 *
 *      the function c(u) is a possibly nonlinear function of the
 *      unknown u, and varies (possibly discontinuously) with the
 *      spatial position also.
 *
 *  user inputs:
 *
 *      the user must provide the coefficients of the differential
 *      operator, some initial parameter settings in an integer and a
 *      real parameter array, and various work arrays.
 */
VEXTERNC void Vmgdriv2(
        int *iparm,    ///< @todo: Doc
        double *rparm, ///< @todo: Doc
        int *nx,       ///< @todo: Doc
        int *ny,       ///< @todo: Doc
        int *nz,       ///< @todo: Doc
        double *u,     ///< @todo: Doc
        int *iz,       ///< @todo: Doc
        int *ipc,      ///< @todo: Doc
        double *rpc,   ///< @todo: Doc
        double *pc,    ///< @todo: Doc
        double *ac,    ///< @todo: Doc
        double *cc,    ///< @todo: Doc
        double *fc,    ///< @todo: Doc
        double *xf,    ///< @todo: Doc
        double *yf,    ///< @todo: Doc
        double *zf,    ///< @todo: Doc
        double *gxcf,  ///< @todo: Doc
        double *gycf,  ///< @todo: Doc
        double *gzcf,  ///< @todo: Doc
        double *a1cf,  ///< @todo: Doc
        double *a2cf,  ///< @todo: Doc
        double *a3cf,  ///< @todo: Doc
        double *ccf,   ///< @todo: Doc
        double *fcf,   ///< @todo: Doc
        double *tcf    ///< @todo: Doc
        );



/** @brief   This routine computes the required sizes of the real and integer
 *           work arrays for the multigrid code.  these two sizes are a
 *           (complicated) function of input parameters.
 *
 *   The work arrays must have been declared in the calling program as:
 *
 *       double precision rwork(iretot)
 *       integer          iwork(iintot)
 *
 *   where:
 *
 *       iretot   = function_of(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev)
 *       iintot   = function_of(mgcoar,mgdisc,mgsolv,nx,ny,nz,nlev)
 *
 *       mgcoar   = coarsening technique:
 *                  0=standard discretization
 *                  1=averaged coefficient + standard discretization
 *                  2=algebraic galerkin coarsening
 *
 *       mgdisc   = discretization technique:
 *                  0=box method
 *                  1=fem method
 *
 *       mgsolv   = coarse grid solver:
 *                  0=conjugate gradients
 *                  1=symmetric banded linpack solver
 *
 *       nx,ny,nz = grid dimensions in each direction,
 *                  including boundary points
 *
 *       nlev     = the number of multigrid levels desired for the
 *                  method.
 *
 *   other parameters:
 *
 *       nf       = number of unknowns on the finest mesh
 *                = nx * ny * nz
 *
 *       nc       = number of unknowns on the coarsest mesh
 *
 *       narr     = storage for one vector on all the meshes
 *
 *       narrc    = storage for one vector on all the meshes but the finest
 *
 *   the work arrays rwork and iwork will be chopped into smaller
 *   pieces according to:
 *
 *       double precision ac(STORE)         (system operators on all levels)
 *       double precision pc(27*narrc)      (prol. opers for coarse levels)
 *       double precision cc(narr),fc(narr) (helmholtz term, rhs -- all levels)
 *       double precision rpc(100*(nlev+1)) (real info for all levels)
 *       integer          ipc(100*(nlev+1)) (integer info for all levels)
 *       integer          iz(50,nlev+1),    (pointers into ac,pc,cc,fc,etc.)
 *
 *   where STORE depends on the discretization, coarsening, and coarse
 *   grid solver:
 *
 *       STORE =  4*nf +  4*narrc + NBAND*nc (mgdisc=box, mgcoar=stan/harm)
 *          or =  4*nf + 14*narrc + NBAND*nc (mgdisc=box, mgcoar=gal)
 *          or = 14*nf + 14*narrc + NBAND*nc (mgdisc=fem, mgcoar=stan/harm/gal)
 *
 *       NBAND = 0                           (mgsolv=iterative)
 *          or = 1+(nxc-2)*(nyc-2)           (mgsolv=7-pt banded linpack)
 *          or = 1+(nxc-2)*(nyc-2)+(nxc-2)+1 (mgsolv=27-pt banded linpack)
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  Replaces mgsz from mgdrvd.f
 */
VEXTERNC void Vmgsz(
        int *mgcoar, ///< @todo: Doc
        int *mgdisc, ///< @todo: Doc
        int *mgsolv, ///< @todo: Doc
        int *nx,     ///< @todo: Doc
        int *ny,     ///< @todo: Doc
        int *nz,     ///< @todo: Doc
        int *nlev,   ///< @todo: Doc
        int *nxc,    ///< @todo: Doc
        int *nyc,    ///< @todo: Doc
        int *nzc,    ///< @todo: Doc
        int *nf,     ///< @todo: Doc
        int *nc,     ///< @todo: Doc
        int *narr,   ///< @todo: Doc
        int *narrc,  ///< @todo: Doc
        int *n_rpc,  ///< @todo: Doc
        int *n_iz,   ///< @todo: Doc
        int *n_ipc,  ///< @todo: Doc
        int *iretot, ///< @todo: Doc
        int *iintot  ///< @todo: Doc
        );



#endif /* _MGDRVD_H_ */
