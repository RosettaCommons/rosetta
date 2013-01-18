/** @defgroup Vopot Vopot class
 *  @brief  Potential oracle for Cartesian mesh data
 */

/**
 *  @file    vopot.h
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @brief   Potential oracle for Cartesian mesh data
 *  @version $Id$
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

#ifndef _VOPOT_H_
#define _VOPOT_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/pbeparm.h"
#include "generic/vatom.h"
#include "generic/valist.h"
#include "generic/vunit.h"
#include "generic/vpbe.h"
#include "generic/pbeparm.h"
#include "mg/vmgrid.h"

/**
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @brief   Electrostatic potential oracle for Cartesian mesh data
 */
struct sVopot {

    Vmgrid *mgrid;  /**< Multiple grid object containing potential data (in
                     * units kT/e) */
    Vpbe   *pbe;  /**< Pointer to PBE object */
    Vbcfl bcfl;  /**< Boundary condition flag for returning potential
                  * values at points off the grid. */
};

/**
 *  @ingroup Vopot
 *  @brief   Declaration of the Vopot class as the Vopot structure
 */
typedef struct sVopot Vopot;

/** @brief   Construct Vopot object with values obtained from Vpmg_readDX (for
 *           example)
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   mgrid  Multiple grid object containing potential data (in units
 *                  kT/e)
 *  @param   pbe   Pointer to Vpbe object for parameters
 *  @param   bcfl  Boundary condition to use for potential values off the grid
 *  @returns Newly allocated and initialized Vopot object
 */
VEXTERNC Vopot*  Vopot_ctor(Vmgrid *mgrid, Vpbe *pbe, Vbcfl bcfl);

/** @brief   Initialize Vopot object with values obtained from Vpmg_readDX (for
 *           example)
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee  Pointer to newly allocated Vopot object
 *  @param   mgrid  Multiple grid object containing potential data (in units
 *                 kT/e)
 *  @param   pbe   Pointer to Vpbe object for parameters
 *  @param   bcfl  Boundary condition to use for potential values off the grid
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vopot_ctor2(Vopot *thee, Vmgrid *mgrid, Vpbe *pbe, Vbcfl bcfl);

/** @brief   Get potential value (from mesh or approximation) at a point
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee  Vopot obejct
 *  @param   x     Point at which to evaluate potential
 *  @param   pot   Set to dimensionless potential (units kT/e) at point x
 *  @returns        1 if successful, 0 otherwise
 */
VEXTERNC int Vopot_pot(Vopot *thee, double x[3], double *pot);

/** @brief   Object destructor
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void Vopot_dtor(Vopot **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void Vopot_dtor2(Vopot *thee);

/** @brief   Get second derivative values at a point
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee   Pointer to Vopot object
 *  @param   pt     Location to evaluate second derivative
 *  @param   cflag
 *             \li  0:  Reduced Maximal Curvature
 *             \li  1:  Mean Curvature (Laplace)
 *             \li  2:  Gauss Curvature
 *             \li  3:  True Maximal Curvature
 *  @param   curv   Set to specified curvature value
 *  @returns        1 if successful, 0 otherwise
 */
VEXTERNC int Vopot_curvature(Vopot *thee, double pt[3], int cflag, double
  *curv);

/** @brief   Get first derivative values at a point
 *  @ingroup Vopot
 *  @author  Nathan Baker
 *  @param   thee   Pointer to Vopot object
 *  @param   pt     Location to evaluate gradient
 *  @param   grad   Gradient
 *  @returns        1 if successful, 0 otherwise
 */
VEXTERNC int Vopot_gradient(Vopot *thee, double pt[3], double grad[3] );


#endif
