/** @defgroup Vpee Vpee class
 *  @brief  This class provides some functionality for error esimation
 *          in parallel.
 *
 *    This class provides some functionality for error esimation in parallel.
 *    The purpose is to modulate the error returned by some external error
 *    estimator according to the partitioning of the mesh.  For example, the
 *    Bank/Holst parallel refinement routine essentially reduces the error
 *    outside the ``local" partition to zero.  However,  this leads to the need
 *    for a few final overlapping Schwarz solves to smooth out the errors near
 *    partition boundaries.  Supposedly, if the region in which we allow
 *    error-based refinement includes the ``local" partition and an external
 *    buffer zone approximately equal in size to the local region, then the
 *    solution will asymptotically approach the solution obtained via more
 *    typical methods.  This is essentially a more flexible parallel
 *    implementation of MC's AM_markRefine.
 */

/**
 *  @file     vpee.h
 *  @ingroup  Vpee
 *  @brief    Contains declarations for class Vpee
 *  @version  $Id$
 *  @author   Nathan A. Baker
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

#ifndef _VPEE_H
#define _VPEE_H

#include "apbscfg.h"

#include "maloc/maloc.h"
#include "mc/mc.h"

/**
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vpee class/module
 */
struct sVpee {

  Gem *gm;                     /**< Grid manager */
  int localPartID;             /**< The local partition ID: i.e. the partition
                                * whose boundary simplices we're keeping
                                * track of */
  double localPartCenter[3];   /**< The coordinates of the center of the local
                                * partition */
  double localPartRadius;      /**< The radius of the circle/sphere which
                                * circumscribes the local partition */
  int killFlag;                /**< A flag indicating the method we're using to
                                * artificially decrease the error esimate
                                * outside the local partition */
  double killParam;            /**< A parameter for the error estimate
                                * attenuation method */
  Vmem *mem;                   /**< Memory manager */

};

/**
 *  @ingroup Vpee
 *  @brief   Declaration of the Vpee class as the Vpee structure
 */
typedef struct sVpee Vpee;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee Inlineable methods
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPEE)
#else /* if defined(VINLINE_VPEE) */
#endif /* if !defined(VINLINE_VPEE) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vpee: Non-Inlineable methods (vpee.c)
/////////////////////////////////////////////////////////////////////////// */

/**
 * @brief   Construct the Vpee object
 * @ingroup Vpee
 * @author  Nathan Baker
 * @return   Newly constructed Vpee object
 */
VEXTERNC Vpee* Vpee_ctor(
        Gem *gm,  /**< FEtk geometry manager object */
        int localPartID,  /**< ID of the local partition (focus of refinement) */
        int killFlag,  /**< A flag to indicate how error estimates are to be
                         attenuated outside the local partition:
                         \li 0:  no attenuation
                         \li 1:  all error outside the local partition set to
                               zero
                         \li 2:  all error is set to zero outside a sphere of
                               radius (killParam*partRadius), where
                               partRadius is the radius of the sphere
                               circumscribing the local partition
                         \li 3:  all error is set to zero except for the local
                               partition and its immediate neighbors */
        double killParam /**< @see killFlag for usage */
        );

/**
 * @brief  FORTRAN stub to construct the Vpee object
 * @ingroup  Vpee
 * @author  Nathan Baker
 * @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vpee_ctor2(
        Vpee *thee,  /**< The Vpee object */
        Gem *gm,  /**< FEtk geometry manager object */
        int localPartID,  /**< ID of the local partition (focus of refinement) */
        int killFlag,  /**< A flag to indicate how error estimates are to be
                         attenuated outside the local partition:
                         \li 0:  no attenuation
                         \li 1:  all error outside the local partition set to
                               zero
                         \li 2:  all error is set to zero outside a sphere of
                               radius (killParam*partRadius), where
                               partRadius is the radius of the sphere
                               circumscribing the local partition
                         \li 3:  all error is set to zero except for the local
                               partition and its immediate neighbors */
        double killParam /**< @see killFlag for usage */
        );

/** @brief   Object destructor
 *  @ingroup Vpee
 *  @author  Nathan Baker
 */
VEXTERNC void Vpee_dtor(
        Vpee **thee /**< Pointer to memory location of the Vpee object */
        );

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vpee
 *  @author  Nathan Baker
 */
VEXTERNC void Vpee_dtor2(
        Vpee *thee /**< Pointer to object to be destroyed */
        );

/** @brief   Mark simplices for refinement based on attenuated error estimates.
 *
 *  A wrapper/reimplementation of AM_markRefine that allows for more flexible
 *  attenuation of error-based markings outside the local partition.  The error
 *  in each simplex is modified by the method (see killFlag) specified in the
 *  Vpee constructor.  This allows the user to confine refinement to an
 *  arbitrary area around the local partition.
 *
 *  @ingroup Vpee
 *  @author  Nathan Baker and Mike Holst
 *  @note  This routine borrows very heavily from FEtk routines by Mike Holst.
 *  @return The number of simplices marked for refinement.
 *  @bug  This function is no longer up-to-date with FEtk and may not function
 *  properly
 */
VEXTERNC int Vpee_markRefine(
        Vpee *thee,  /**< The Vpee object */
        AM *am,  /**< The FEtk algebra manager currently used to solve the PB */
        int level,  /**< The current level of the multigrid hierarchy */
        int akey,  /**< The marking method:
                      \li -1:  Reset markings  --> killFlag has no effect.
                      \li 0:  Uniform.
                      \li 1:  User defined (geometry-based).
                      \li >1:  A numerical estimate for the error has already been
                               set in am and should be attenuated according to
                               killFlag and used, in conjunction with etol, to mark
                               simplices for refinement. */
        int rcol, /**< The ID of the main parition on which to mark (or -1 if
                    all partitions should be marked).  NOte that we shouldhave
                    (rcol == thee->localPartID) for (thee->killFlag == 2 or 3) */
        double etol,  /**< The error tolerance criterion for marking */
        int bkey  /**< How the error tolerance is interpreted:
                     \li 0:  Simplex marked if error > etol.
                     \li 1:  Simplex marked if error >
                     sqrt(etol^2/L) where L$ is the number of simplices */
        );

/** @brief   Returns the number of simplices in the local partition
 *  @ingroup Vpee
 *  @author  Nathan Baker
 *  @return  Number of simplices in the local partition
 */
VEXTERNC int Vpee_numSS(
        Vpee *thee /**< The Vpee object */
        );

#endif    /* ifndef _VPEE_H_ */
