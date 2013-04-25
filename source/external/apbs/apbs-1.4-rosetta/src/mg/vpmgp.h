/** @defgroup Vpmgp Vpmgp class
 *  @brief  Parameter structure for Mike Holst's PMGP code
 *  @note   Variables and many default values taken directly from PMG
 */

/**
 *  @file     vpmgp.h
 *  @ingroup  Vpmgp
 *  @brief    Contains declarations for class Vpmgp
 *  @version  $Id$
 *  @author   Nathan A. Baker
 *  @note     Variables and many default values taken directly from PMG
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 *  Nathan A. Baker (nathan.baker@pnnl.gov)
 *  Pacific Northwest National Laboratory
 *
 *  Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2012 Battelle Memorial Institute. Developed at the
 * Pacific Northwest National Laboratory, operated by Battelle Memorial
 * Institute, Pacific Northwest Division for the U.S. Department of Energy.
 *
 * Portions Copyright (c) 2002-2010, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2010, Nathan A. Baker.
 * Portions Copyright (c) 1999-2002, The Regents of the University of
 * California.
 * Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of the developer nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#ifndef _VPMGP_H_
#define _VPMGP_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/mgparm.h"

/**
 *  @ingroup Vpmgp
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vpmgp class/module
 *  @bug     Value ipcon does not currently allow for preconditioning in PMG
 */
struct sVpmgp {

    /* ********** USER-SPECIFIED PARAMETERS ********** */
    int nx;  /**< Grid x dimensions [no default]  */
    int ny;  /**< Grid y dimensions [no default]  */
    int nz;  /**< Grid z dimensions [no default]  */
    int nlev;  /**< Number of mesh levels [no default] */
    double hx;  /**< Grid x spacings [no default]  */
    double hy;  /**< Grid y spacings [no default]  */
    double hzed;  /**< Grid z spacings [no default]  */
    int nonlin;  /**< Problem type [no default]
                  * \li 0: linear
                  * \li 1: nonlinear
                  * \li 2: linear then nonlinear */

    /* ********** DERIVED PARAMETERS ********** */
    int nxc;  /**< Coarse level grid x dimensions */
    int nyc;  /**< Coarse level grid y dimensions */
    int nzc;  /**< Coarse level grid z dimensions */
    int nf;  /**< Number of fine grid unknowns */
    int nc;  /**< Number of coarse grid unknowns */
    int narrc;  /**< Size of vector on coarse level */
    int n_rpc;  /**< Real info work array required storage */
    int n_iz;  /**< Integer storage parameter (index max) */
    int n_ipc;  /**< Integer info work array required storage */

    int nrwk;  /**< Real work storage */
    int niwk;  /**< Integer work storage */
    int narr;  /**< Array work storage */
    int ipkey;  /**< Toggles nonlinearity (set by nonlin)
                 * \li  -2: Size-Modified PBE
                 * \li  -1: Linearized PBE
                 * \li   0: Nonlinear PBE with capped sinh
                 *          term [default]
                 * \li  >1: Polynomial approximation to sinh,
                 *          note that ipkey must be odd  */

    /* ********** PARAMETERS WITH DEFAULT VALUES ********** */
    double xcent;  /**< Grid x center [0]  */
    double ycent;  /**< Grid y center [0]  */
    double zcent;  /**< Grid z center [0]  */
    double errtol;  /**< Desired error tolerance [default = 1e-9] */
    int itmax;  /**< Maximum number of iters [default = 100] */
    int istop;  /**< Stopping criterion [default = 1]
                 * \li 0: residual
                 * \li 1: relative residual
                 * \li 2: diff
                 * \li 3: errc
                 * \li 4: errd
                 * \li 5: aerrd */
    int iinfo;  /**< Runtime status messages [default = 1]
                 * \li 0: none
                 * \li 1: some
                 * \li 2: lots
                 * \li 3: more */
    Vbcfl bcfl;  /**< Boundary condition method [default = BCFL_SDH] */
    int key;  /**< Print solution to file [default = 0]
               * \li   0: no
               * \li   1: yes */
    int iperf;  /**< Analysis of the operator [default = 0]
                 * \li   0: no
                 * \li   1: condition number
                 * \li   2: spectral radius
                 * \li   3: cond. number & spectral radius */
    int meth;  /**< Solution method [default = 2]
                * \li   0: conjugate gradient multigrid
                * \li   1: newton
                * \li   2: multigrid
                * \li   3: conjugate gradient
                * \li   4: sucessive overrelaxation
                * \li   5: red-black gauss-seidel
                * \li   6: weighted jacobi
                * \li   7: richardson
                * \li   8: conjugate gradient multigrid aqua
                * \li   9: newton aqua */
    int mgkey;  /**< Multigrid method [default = 0]
                 * \li   0: variable v-cycle
                 * \li   1: nested iteration */
    int nu1;  /**< Number of pre-smoothings [default = 2] */
    int nu2;  /**< Number of post-smoothings [default = 2] */
    int mgsmoo;  /**< Smoothing method [default = 1]
                  * \li   0: weighted jacobi
                  * \li   1: gauss-seidel
                  * \li   2: SOR
                  * \li   3: richardson
                  * \li   4: cghs */
    int mgprol;  /**< Prolongation method [default = 0]
                  * \li   0: trilinear
                  * \li   1: operator-based
                  * \li   2: mod. operator-based */
    int mgcoar;  /**< Coarsening method [default = 2]
                  * \li   0: standard
                  * \li   1: harmonic
                  * \li   2: galerkin */
    int mgsolv;  /**< Coarse equation solve method [default = 1]
                  * \li   0: cghs
                  * \li   1: banded linpack */
    int mgdisc;  /**< Discretization method [default = 0]
                  * \li   0: finite volume
                  * \li   1: finite element */
    double omegal;  /**< Linear relax parameter [default = 8e-1] */
    double omegan;  /**< Nonlin relax parameter [default = 9e-1] */
    int irite;  /**< FORTRAN output unit [default = 8] */
    int ipcon;  /**< Preconditioning method [default = 3]
                 * \li   0: diagonal
                 * \li   1: ICCG
                 * \li   2: ICCGDW
                 * \li   3: MICCGDW
                 * \li   4: none */
    double xlen;  /**< Domain x length */
    double ylen;  /**< Domain y length */
    double zlen;  /**< Domain z length */
    double xmin;  /**< Domain lower x corner */
    double ymin;  /**< Domain lower y corner */
    double zmin;  /**< Domain lower z corner */
    double xmax;  /**< Domain upper x corner */
    double ymax;  /**< Domain upper y corner */
    double zmax;  /**< Domain upper z corner */
};

/**
 *  @ingroup Vpmgp
 *  @brief   Declaration of the Vpmgp class as the sVpmgp structure
 */
typedef struct sVpmgp Vpmgp;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Inlineable methods (vpmgp.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPMGP)
#else /* if defined(VINLINE_VPMGP) */
#endif /* if !defined(VINLINE_VPMGP) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpmgp: Non-Inlineable methods (vpmgp.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct PMG parameter object and initialize to default values
 *  @ingroup Vpmgp
 *  @author  Nathan Baker
 *  @param   mgparm	  MGParm object containing parameters to be used in setup
 *  @returns Newly allocated and initialized Vpmgp object
 */
VEXTERNC Vpmgp* Vpmgp_ctor(MGparm *mgparm);

/** @brief   FORTRAN stub to construct PMG parameter object and initialize to
 *           default values
 *  @ingroup Vpmgp
 *  @author  Nathan Baker
 *  @param   thee  Newly allocated PMG object
 *  @param   mgparm	  MGParm object containing parameters to be used in setup
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vpmgp_ctor2(Vpmgp *thee, MGparm *mgparm);

/** @brief   Object destructor
 *  @ingroup Vpmgp
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location for Vpmgp object
 */
VEXTERNC void Vpmgp_dtor(Vpmgp **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup Vpmgp
 *  @author  Nathan Baker
 *  @param   thee  Pointer to Vpmgp object
 */
VEXTERNC void Vpmgp_dtor2(Vpmgp *thee);

/**	@brief	Determine array sizes and parameters for multigrid solver
*	@ingroup Vpmgp
*	@author	Mike Holst and Nathan Baker
*/
VEXTERNC void Vpmgp_size(
    Vpmgp *thee	/**< Object to be sized */
    );

/**	@brief	Coarsen the grid by the desired number of levels and determine the resulting numbers of grid points.
*	@ingroup Vpmgp
*	@author	Mike Holst and Nathan Baker
*/
VEXTERNC void Vpmgp_makeCoarse(
    int numLevel,	/**< Number of levels to coarsen */
    int nxOld,	/**< Number of old grid points in this direction */
    int nyOld,	/**< Number of old grid points in this direction */
    int nzOld,	/**< Number of old grid points in this direction */
    int *nxNew,	/**< Number of new grid points in this direction */
    int *nyNew,	/**< Number of new grid points in this direction */
    int *nzNew	/**< Number of new grid points in this direction */
    );



#endif    /* ifndef _VPMGP_H_ */
