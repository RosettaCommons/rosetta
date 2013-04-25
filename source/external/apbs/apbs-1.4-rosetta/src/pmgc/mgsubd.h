/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief Multigrid subroutines
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

#ifndef _MGSUBD_H_
#define _MGSUBD_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vmatrix.h"
#include "pmgc/matvecd.h"
#include "pmgc/buildAd.h"
#include "pmgc/buildPd.h"
#include "pmgc/buildBd.h"
#include "pmgc/buildGd.h"

#define HARMO2(a, b)                   (2.0 * (a) * (b) / ((a) + (b)))
#define HARMO4(a, b, c, d)             (1.0 / ( 0.25 * ( 1.0/(a) + 1.0/(b) + 1.0/(c) + 1.0/(d))))
#define ARITH2(a, b)                   (((a) + (b)) / 2.0)
#define ARITH4(a, b, c, d)             (((a) + (b) + (c) + (d)) / 4.0)
#define ARITH6(a, b, c, d, e, f)       (((a) + (b) + (c) + (d) + (e) + (f)) / 6.0)
#define ARITH8(a, b, c, d, e, f, g, h) (((a) + (b) + (c) + (d) + (e) + (f) + (g) + (h)) / 8.0)

/** @brief   Build operators, boundary arrays, modify affine vectors
 *             ido==0: do only fine level
 *             ido==1: do only coarse levels (including second op at coarsest)
 *             ido==2: do all levels
 *             ido==3: rebuild the second operator at the coarsest level
 *
 *  @note    The fine level must be build before any coarse levels.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  Replaces buildops from mgsubd.f
 */
VEXTERNC void Vbuildops(
        int    *nx,     ///< @todo: doc
        int    *ny,     ///< @todo: doc
        int    *nz,     ///< @todo: doc
        int    *nlev,   ///< @todo: doc
        int    *ipkey,  ///< @todo: doc
        int    *iinfo,  ///< @todo: doc
        int    *ido,    ///< @todo: doc
        int    *iz,     ///< @todo: doc
        int    *mgprol, ///< @todo: doc
        int    *mgcoar, ///< @todo: doc
        int    *mgsolv, ///< @todo: doc
        int    *mgdisc, ///< @todo: doc
        int    *ipc,    ///< @todo: doc
        double *rpc,    ///< @todo: doc
        double *pc,     ///< @todo: doc
        double *ac,     ///< @todo: doc
        double *cc,     ///< @todo: doc
        double *fc,     ///< @todo: doc
        double *xf,     ///< @todo: doc
        double *yf,     ///< @todo: doc
        double *zf,     ///< @todo: doc
        double *gxcf,   ///< @todo: doc
        double *gycf,   ///< @todo: doc
        double *gzcf,   ///< @todo: doc
        double *a1cf,   ///< @todo: doc
        double *a2cf,   ///< @todo: doc
        double *a3cf,   ///< @todo: doc
        double *ccf,    ///< @todo: doc
        double *fcf,    ///< @todo: doc
        double *tcf     ///< @todo: doc
        );

/** @brief   Build the nexted operator framework in the array iz
 *  @note    iz(50,i) indexes into the gridfcn arrays
 *           for each level i=(1,...,nlev+1) as follows:
 *
 *           fun(i)    = fun (iz(1,i))
 *           bndx(i)   = bndx(iz(2,i))
 *           bndy(i)   = bndy(iz(3,i))
 *           bndz(i)   = bndz(iz(4,i))
 *           ipc(i)    = ipc(iz(5,i))
 *           rpc(i)    = rpc(iz(6,i))
 *           oper(i)   = oper(iz(7,i))
 *           grdx(i)   = brdx(iz(8,i))
 *           grdy(i)   = brdy(iz(9,i))
 *           grdz(i)   = brdz(iz(10,i))
 *
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces buildstr from mgsubd.f
 */
VEXTERNC void Vbuildstr(
        int *nx,   ///< @todo  Document
        int *ny,   ///< @todo  Document
        int *nz,   ///< @todo  Document
        int *nlev, ///< @todo  Document
        int *iz    ///< @todo  Document
        );

/** @brief   Form the Galerkin coarse grid system
 *  @note    Although the fine grid matrix may be 7 or 27 diagonal,
 *           the coarse grid matrix is always 27 diagonal.
 *           (only 14 stored due to symmetry.)
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces buildgaler0 from mgsubd.f
 */
VEXTERNC void Vbuildgaler0(
        int    *nxf,    ///< @todo: doc
        int    *nyf,    ///< @todo: doc
        int    *nzf,    ///< @todo: doc
        int    *nxc,    ///< @todo: doc
        int    *nyc,    ///< @todo: doc
        int    *nzc,    ///< @todo: doc
        int    *ipkey,  ///< @todo: doc
        int    *numdia, ///< @todo: doc
        double *pcFF,   ///< @todo: doc
        int    *ipcFF,  ///< @todo: doc
        double *rpcFF,  ///< @todo: doc
        double *acFF,   ///< @todo: doc
        double *ccFF,   ///< @todo: doc
        double *fcFF,   ///< @todo: doc
        int    *ipc,    ///< @todo: doc
        double *rpc,    ///< @todo: doc
        double *ac,     ///< @todo: doc
        double *cc,     ///< @todo: doc
        double *fc      ///< @todo: doc
        );



/** @brief   Coarsen a grid
 *           Compute the number of grid points in the coarser grid, given the
 *           number of grid points in a finger grid in each direction.
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces mkcors from mgsubd.f
 */
VEXTERNC void Vmkcors(
        int *numlev, ///< @todo: doc
        int *nxold,  ///< @todo: doc
        int *nyold,  ///< @todo: doc
        int *nzold,  ///< @todo: doc
        int *nxnew,  ///< @todo: doc
        int *nynew,  ///< @todo: doc
        int *nznew   ///< @todo: doc
        );

/** @brief   Coarsen a grid
 *           Compute the number of grid points in the coarser grid, given the
 *           number of grid points in a finger grid in each direction.
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces mkcors from mgsubd.f
 */
VEXTERNC void Vcorsr(
        int *nold, ///< @todo: doc
        int *nnew  ///< @todo: doc
        );



/** @brief   Refine a grid
 *           Compute the number of grid points in the finer grid, given the
 *           number of grid points in a coarser grid in each direction.
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces mkfine from mgsubd.f
 */
VEXTERNC void Vmkfine(
        int *numlev, ///< @todo: doc
        int *nxold,  ///< @todo: doc
        int *nyold,  ///< @todo: doc
        int *nzold,  ///< @todo: doc
        int *nxnew,  ///< @todo: doc
        int *nynew,  ///< @todo: doc
        int *nznew   ///< @todo: doc
        );

/** @brief   Refine a grid
 *           Compute the number of grid points in the finer grid, given the
 *           number of grid points in a coarser grid.
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces mkfine from mgsubd.f
 */
VEXTERNC void Vfiner(
        int *nold, ///< @todo: doc
        int *nnew  ///< @todo: doc
        );



/** @brief   Coarsen a single dimension of a grid
 *           Compute the number of grid points in the coarser grid, given the
 *           number of grid points in a finer grid.
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces corsr from mgsubd.f
 */
VEXTERNC int Vivariv(
        int *nu,   ///< @todo: doc
        int *level ///< @todo: doc
        );

/** @brief   Find maximum multigrid possible coarsenning
 *           common to three grid sizes
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces maxlev from mgsubd.f
 */
VEXTERNC int Vmaxlev(
        int n1, ///< The first grid size
        int n2, ///< The second grid size
        int n3  ///< The third grid size
        );


/// @todo  Get rid of these globals in refactor
double bf, oh, cputme;

/** @brief   This routine prints out some info and such from inside multigrid.
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces prtstp from mgsubd.f
 */
VEXTERNC void Vprtstp(
        int  iok,      ///< @todo  Document
        int  iters,    ///< @todo  Document
        double  rsnrm, ///< @todo  Document
        double  rsden, ///< @todo  Document
        double  orsnrm ///< @todo  Document
        );



/** @brief   Print out a column-compressed sparse matrix in Harwell-Boeing
 *           format.
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @bug  Can this path variable be replaced with a Vio socket?
 */
VEXTERNC void Vpackmg(
        int *iparm,
        double *rparm,
        int *nrwk,
        int *niwk,
        int *nx,
        int *ny,
        int *nz,
        int *nlev,
        int *nu1,
        int *nu2,
        int *mgkey,
        int *itmax,
        int *istop,
        int *ipcon,
        int *nonlin,
        int *mgsmoo,
        int *mgprol,
        int *mgcoar,
        int *mgsolv,
        int *mgdisc,
        int *iinfo,
        double *errtol,
        int *ipkey,
        double *omegal,
        double *omegan,
        int *irite,
        int *iperf
        );



/** @brief   Produce information for a coarser grid.
 *           Also harmonically average the problem coefficients.
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces buildharm0 from mgsubd.f
 */
VEXTERNC void Vbuildharm0(
        int *nx,      ///< @todo  Document
        int *ny,      ///< @todo  Document
        int *nz,      ///< @todo  Document
        int *nxf,     ///< @todo  Document
        int *nyf,     ///< @todo  Document
        int *nzf,     ///< @todo  Document
        double *xc,   ///< @todo  Document
        double *yc,   ///< @todo  Document
        double *zc,   ///< @todo  Document
        double *gxc,  ///< @todo  Document
        double *gyc,  ///< @todo  Document
        double *gzc,  ///< @todo  Document
        double *a1c,  ///< @todo  Document
        double *a2c,  ///< @todo  Document
        double *a3c,  ///< @todo  Document
        double *cc,   ///< @todo  Document
        double *fc,   ///< @todo  Document
        double *tc,   ///< @todo  Document
        double *xf,   ///< @todo  Document
        double *yf,   ///< @todo  Document
        double *zf,   ///< @todo  Document
        double *gxcf, ///< @todo  Document
        double *gycf, ///< @todo  Document
        double *gzcf, ///< @todo  Document
        double *a1cf, ///< @todo  Document
        double *a2cf, ///< @todo  Document
        double *a3cf, ///< @todo  Document
        double *ccf,  ///< @todo  Document
        double *fcf,  ///< @todo  Document
        double *tcf   ///< @todo  Document
        );



/** @brief   Produce information for a coarser grid.
 *           Also harmonically average the problem coefficients.
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces buildharm0 from mgsubd.f
 */
VEXTERNC void Vbuildcopy0(
        int    *nx,   ///< @todo  Document
        int    *ny,   ///< @todo  Document
        int    *nz,   ///< @todo  Document
        int    *nxf,  ///< @todo  Document
        int    *nyf,  ///< @todo  Document
        int    *nzf,  ///< @todo  Document
        double *xc,   ///< @todo  Document
        double *yc,   ///< @todo  Document
        double *zc,   ///< @todo  Document
        double *gxc,  ///< @todo  Document
        double *gyc,  ///< @todo  Document
        double *gzc,  ///< @todo  Document
        double *a1c,  ///< @todo  Document
        double *a2c,  ///< @todo  Document
        double *a3c,  ///< @todo  Document
        double *cc,   ///< @todo  Document
        double *fc,   ///< @todo  Document
        double *tc,   ///< @todo  Document
        double *xf,   ///< @todo  Document
        double *yf,   ///< @todo  Document
        double *zf,   ///< @todo  Document
        double *gxcf, ///< @todo  Document
        double *gycf, ///< @todo  Document
        double *gzcf, ///< @todo  Document
        double *a1cf, ///< @todo  Document
        double *a2cf, ///< @todo  Document
        double *a3cf, ///< @todo  Document
        double *ccf,  ///< @todo  Document
        double *fcf,  ///< @todo  Document
        double *tcf   ///< @todo  Document
        );



/** @brief   Build RHS algebraically for analysis purposes
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *  @note    Replaces buildALG from mgsubd.f
 *           The fine level must be built before any coarse levels
 */
VEXTERNC void Vbuildalg(
        int      *nx, ///< @todo  Document
        int      *ny, ///< @todo  Document
        int      *nz, ///< @todo  Document
        int    *mode, ///< @todo  Document
        int    *nlev, ///< @todo  Document
        int      *iz, ///< @todo  Document
        int     *ipc, ///< @todo  Document
        double  *rpc, ///< @todo  Document
        double   *ac, ///< @todo  Document
        double   *cc, ///< @todo  Document
        double   *fc, ///< @todo  Document
        double    *x, ///< @todo  Document
        double    *y, ///< @todo  Document
        double  *tmp  ///< @todo  Document
);



#endif // _MGSUBD_H_
