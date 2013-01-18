/** @defgroup Vmgrid Vmgrid class
 *  @brief    Oracle for Cartesian mesh data
 */

/**
 *  @file    vmgrid.h
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @brief   Multiresolution oracle for Cartesian mesh data
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

#ifndef _VMGRID_H_
#define _VMGRID_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "mg/vgrid.h"

/** @def VMGRIDMAX
 *  @brief The maximum number of levels in the grid hiearchy
 *  @ingroup Vmgrid
 */
#define VMGRIDMAX 20


/**
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @brief   Multiresoltion oracle for Cartesian mesh data
 */
struct sVmgrid {

    int ngrids;                /**< Number of grids in hiearchy */
    Vgrid *grids[VMGRIDMAX];   /**< Grids in hiearchy.  Our convention will be
                                *   to have the finest grid first, however,
                                *   this will not be enforced as it may be
                                *   useful to search multiple grids for
                                *   parallel datasets, etc. */
};

/**
 *  @ingroup Vmgrid
 *  @brief   Declaration of the Vmgrid class as the Vgmrid structure
 */
typedef struct sVmgrid Vmgrid;

/** @brief   Construct Vmgrid object
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @returns Newly allocated and initialized Vmgrid object
 */
VEXTERNC Vmgrid*  Vmgrid_ctor();

/** @brief   Initialize Vmgrid object
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee Newly allocated Vmgrid object
 *  @returns Newly allocated and initialized Vmgrid object
 */
VEXTERNC int Vmgrid_ctor2(Vmgrid *thee);

/** @brief   Get potential value (from mesh or approximation) at a point
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee  Vmgrid obejct
 *  @param   x     Point at which to evaluate potential
 *  @param   value Value of data at point x
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vmgrid_value(Vmgrid *thee, double x[3], double *value);

/** @brief   Object destructor
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void Vmgrid_dtor(Vmgrid **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void Vmgrid_dtor2(Vmgrid *thee);

/** @brief   Add a grid to the hierarchy
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 *  @param   grid   Grid to be added.  As mentioned above, we would prefer to
 *           have the finest grid added first, next-finest second, ...,
 *           coarsest last -- this is how the grid will be searched when
 *           looking up values for points.  However, this is not enforced to
 *           provide flexibility for cases where the dataset is decomposed into
 *           disjoint partitions, etc.
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vmgrid_addGrid(Vmgrid *thee, Vgrid *grid);


/** @brief   Get second derivative values at a point
 *  @ingroup Vmgrid
 *  @author  Nathan Baker (wrapper for Vgrid routine by Steve Bond)
 *  @param   thee   Pointer to Vmgrid object
 *  @param   pt     Location to evaluate second derivative
 *  @param   cflag
 *             \li  0:  Reduced Maximal Curvature
 *             \li  1:  Mean Curvature (Laplace)
 *             \li  2:  Gauss Curvature
 *             \li  3:  True Maximal Curvature
 *  @param   curv Specified curvature value
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vmgrid_curvature(Vmgrid *thee, double pt[3], int cflag,
  double *curv);

/** @brief   Get first derivative values at a point
 *  @ingroup Vmgrid
 *  @author  Nathan Baker and Steve Bond
 *  @param   thee   Pointer to Vmgrid object
 *  @param   pt     Location to evaluate gradient
 *  @param   grad   Gradient
 *  @return  1 if successful, 0 if off grid
 */
VEXTERNC int Vmgrid_gradient(Vmgrid *thee, double pt[3], double grad[3] );

/** @brief   Get specific grid in hiearchy
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to Vmgrid object
 *  @param   num    Number of grid in hiearchy
 *  @return  Pointer to specified grid
 */
VEXTERNC Vgrid* Vmgrid_getGridByNum(Vmgrid *thee, int num);

/** @brief   Get grid in hiearchy which contains specified point or VNULL
 *  @ingroup Vmgrid
 *  @author  Nathan Baker
 *  @param   thee   Pointer to Vmgrid object
 *  @param   pt     Point to check
 *  @return  Pointer to specified grid
 */
VEXTERNC Vgrid* Vmgrid_getGridByPoint(Vmgrid *thee, double pt[3]);

#endif

