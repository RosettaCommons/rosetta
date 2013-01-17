/** @defgroup Vgreen Vgreen class
 *  @brief    Provides capabilities for pointwise evaluation of free space
 *            Green's function for point charges in a uniform dielectric.
 *  @note     Right now, these are very slow methods without any fast multipole
 *            acceleration.
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

/**
 *  @file     vgreen.h
 *  @ingroup  Vgreen
 *  @brief    Contains declarations for class Vgreen
 *  @version  $Id$
 *  @author   Nathan A. Baker
 */

#ifndef _VGREEN_H_
#define _VGREEN_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vunit.h"
#include "generic/vatom.h"
#include "generic/valist.h"

/**
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vgreen class/module
 */
struct sVgreen {

  Valist *alist;  /**< Atom (charge) list for Green's function */
  Vmem *vmem;  /**< Memory management object */
  double *xp;  /**< Array of particle x-coordinates for use with
                * treecode routines */
  double *yp;  /**< Array of particle y-coordinates for use with
                * treecode routines */
  double *zp;  /**< Array of particle z-coordinates for use with
                * treecode routines */
  double *qp;  /**< Array of particle charges for use with
                * treecode routines */
  int np;  /**< Set to size of above arrays */
};

/**
 *  @ingroup Vgreen
 *  @brief   Declaration of the Vgreen class as the Vgreen structure
 */
typedef struct sVgreen Vgreen;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vgreen: Inlineable methods (vgreen.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VGREEN)

    /** @brief   Get the atom list associated with this Green's function object
     *  @ingroup Vgreen
     *  @author  Nathan Baker
     *  @param   thee  Vgreen object
     *  @return  Pointer to Valist object associated with this Green's function
     *           object
     */
    VEXTERNC Valist* Vgreen_getValist(Vgreen *thee);

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vgreen
     *  @author  Nathan Baker
     *  @param   thee  Vgreen object
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC unsigned long int Vgreen_memChk(Vgreen *thee);

#else /* if defined(VINLINE_VGREEN) */
#   define Vgreen_getValist(thee) ((thee)->alist)
#   define Vgreen_memChk(thee) (Vmem_bytes((thee)->vmem))
#endif /* if !defined(VINLINE_VGREEN) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vgreen: Non-Inlineable methods (vgreen.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct the Green's function oracle
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   alist  Atom (charge) list associated with object
 *  @return  Pointer to newly allocated Green's function oracle
 */
VEXTERNC Vgreen* Vgreen_ctor(Valist *alist);

/** @brief   FORTRAN stub to construct the Green's function oracle
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory allocated for object
 *  @param   alist  Atom (charge) list associated with object
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vgreen_ctor2(Vgreen *thee, Valist *alist);

/** @brief   Destruct the Green's function oracle
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee Pointer to memory location for object
 */
VEXTERNC void Vgreen_dtor(Vgreen **thee);

/** @brief   FORTRAN stub to destruct the Green's function oracle
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee Pointer to object
 */
VEXTERNC void Vgreen_dtor2(Vgreen *thee);

/** @brief   Get the Green's function for Helmholtz's equation integrated over
 *           the atomic point charges
 *
 *           Returns the potential \f$\phi\f$ defined by
 *           \f[ \phi(r) = \sum_i \frac{q_i e^{-\kappa r_i}}{r_i} \f]
 *
 *           where \f$\kappa\f$ is the inverse screening length (in &Aring;)
 *           \f$q_i\f$ is the atomic charge (in e), and \f$r_i\f$ r_i is the
 *           distance from atom \f$i\f$ to the observation point \f$r\f$.  The
 *           potential is scaled to units of V.
 *
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @bug     Not implemented yet
 *  @note    Not implemented yet
 *  @param   thee  Vgreen object
 *  @param   npos  Number of positions to evaluate
 *  @param   x  The npos x-coordinates
 *  @param   y  The npos y-coordinates
 *  @param   z  The npos z-coordinates
 *  @param   val  The npos values
 *  @param   kappa The value of \f$\kappa\f$ (see above)
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vgreen_helmholtz(Vgreen *thee, int npos, double *x, double *y,
  double *z, double *val, double kappa);

/** @brief   Get the gradient of Green's function for Helmholtz's equation
 *           integrated over the atomic point charges
 *
 *           Returns the field \f$\nabla \phi\f$ defined by
 *           \f[ \nabla \phi(r) = \nabla \sum_i \frac{q_i e^{-\kappa r_i}}{r_i}
 *           \f]
 *
 *           where \f$\kappa\f$ is the inverse screening length (in &Aring;).
 *           \f$q_i\f$ is the atomic charge (in e), and \f$r_i\f$ r_i is the
 *           distance from atom \f$i\f$ to the observation point \f$r\f$.  The
 *           potential is scaled to units of V/&Aring;.
 *
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @bug     Not implemented yet
 *  @note    Not implemented yet
 *  @param   thee  Vgreen object
 *  @param   npos  The number of positions to evaluate
 *  @param   x  The npos x-coordinates
 *  @param   y  The npos y-coordinates
 *  @param   z  The npos z-coordinates
 *  @param   gradx  The npos gradient x-components
 *  @param   grady  The npos gradient y-components
 *  @param   gradz  The npos gradient z-components
 *  @param   kappa The value of \f$\kappa\f$ (see above)
 *  @return  int 1 if sucessful, 0 otherwise
 */
VEXTERNC int Vgreen_helmholtzD(Vgreen *thee, int npos, double *x, double *y,
  double *z, double *gradx, double *grady, double *gradz, double kappa);

/** @brief   Get the Coulomb's Law Green's function (solution to Laplace's
 *           equation) integrated over the atomic point charges using direct
 *           summation
 *
 *           Returns the potential \f$\phi\f$ defined by
 *           \f[ \phi(r) = \sum_i \frac{q_i}{r_i} \f]
 *           where \f$q_i\f$ is the atomic charge (in e) and \f$r_i\f$ is the
 *           distance to the observation point \f$r\f$.  The potential is
 *           scaled to units of V.
 *
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee Vgreen object
 *  @param   npos  The number of positions to evaluate
 *  @param   x  The npos x-coordinates
 *  @param   y  The npos y-coordinates
 *  @param   z  The npos z-coordinates
 *  @param   val  The npos values
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vgreen_coulomb_direct(Vgreen *thee, int npos, double *x,
        double *y, double *z, double *val);

/** @brief   Get the Coulomb's Law Green's function (solution to Laplace's
 *           equation) integrated over the atomic point charges using direct
 *           summation or H. E. Johnston, R. Krasny FMM library (if available)
 *
 *           Returns the potential \f$\phi\f$ defined by
 *           \f[ \phi(r) = \sum_i \frac{q_i}{r_i} \f]
 *           where \f$q_i\f$ is the atomic charge (in e) and \f$r_i\f$ is the
 *           distance to the observation point \f$r\f$.  The potential is
 *           scaled to units of V.
 *
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee Vgreen object
 *  @param   npos  The number of positions to evaluate
 *  @param   x  The npos x-coordinates
 *  @param   y  The npos y-coordinates
 *  @param   z  The npos z-coordinates
 *  @param   val  The npos values
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vgreen_coulomb(Vgreen *thee, int npos, double *x, double *y,
  double *z, double *val);

/** @brief   Get gradient of the Coulomb's Law Green's function (solution to
 *           Laplace's equation) integrated over the atomic point charges using
 *           direct summation
 *
 *           Returns the field \f$\nabla \phi\f$ defined by
 *           \f[ \nabla \phi(r) = \sum_i \frac{q_i}{r_i} \f]
 *           where \f$q_i\f$ is the atomic charge (in e) and \f$r_i\f$ is the
 *           distance to the observation point \f$r\f$.  The field is
 *           scaled to units of V/&Aring;.
 *
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee Vgreen object
 *  @param   npos  The number of positions to evaluate
 *  @param   x  The npos x-coordinates
 *  @param   y  The npos y-coordinates
 *  @param   z  The npos z-coordinates
 *  @param   pot    The npos potential values
 *  @param   gradx  The npos gradient x-components
 *  @param   grady  The npos gradient y-components
 *  @param   gradz  The npos gradient z-components
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vgreen_coulombD_direct(Vgreen *thee, int npos, double *x,
        double *y, double *z, double *pot, double *gradx, double *grady, double
        *gradz);

/** @brief   Get gradient of the Coulomb's Law Green's function (solution to
 *           Laplace's equation) integrated over the atomic point charges using
 *           either direct summation or H. E. Johnston/R. Krasny FMM library
 *           (if available)
 *
 *           Returns the field \f$\nabla \phi\f$ defined by
 *           \f[ \nabla \phi(r) = \sum_i \frac{q_i}{r_i} \f]
 *           where \f$q_i\f$ is the atomic charge (in e) and \f$r_i\f$ is the
 *           distance to the observation point \f$r\f$.  The field is
 *           scaled to units of V/&Aring;.
 *
 *  @ingroup Vgreen
 *  @author  Nathan Baker
 *  @param   thee Vgreen object
 *  @param   npos  The number of positions to evaluate
 *  @param   x  The npos x-coordinates
 *  @param   y  The npos y-coordinates
 *  @param   z  The npos z-coordinates
 *  @param   pot    The npos potential values
 *  @param   gradx  The npos gradient x-components
 *  @param   grady  The npos gradient y-components
 *  @param   gradz  The npos gradient z-components
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vgreen_coulombD(Vgreen *thee, int npos, double *x, double *y,
  double *z, double *pot, double *gradx, double *grady, double *gradz);

#endif /* ifndef _VGREEN_H_ */
