/** @defgroup Vclist Vclist class
 *  @brief    Atom cell list
 */

/**
 *  @file     vclist.h
 *  @ingroup  Vclist
 *  @brief    Contains declarations for class Vclist
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

#ifndef _VCLIST_H_
#define _VCLIST_H_

#include "apbscfg.h"

#include "maloc/maloc.h"
#if defined(HAVE_MC_H)
#include "mc/mc.h"
#endif

#include "generic/vhal.h"
#include "generic/valist.h"
#include "generic/vatom.h"
#include "generic/vunit.h"

/**
 *  @ingroup Vclist
 *  @author  Nathan Baker
 *  @brief   Atom cell list domain setup mode
 */
enum eVclist_DomainMode {
    CLIST_AUTO_DOMAIN,  /**< Setup the cell list domain automatically to
                         * encompass the entire molecule */
    CLIST_MANUAL_DOMAIN   /**< Specify the cell list domain manually through
                           * the constructor */
};

/**
 * @typedef Vclist_DomainMode
 * @ingroup Vclist
 * @brief Declaration of Vclist_DomainMode enumeration type
 */
typedef enum eVclist_DomainMode Vclist_DomainMode;

/**
 * @ingroup Vclist
 * @author Nathan Baker
 * @brief Atom cell list cell
 */
struct sVclistCell {
    Vatom **atoms;  /**< Array of atom objects associated with this cell */
    int natoms;  /**< Length of thee->atoms array */
};

/**
 *  @ingroup Vclist
 *  @brief   Declaration of the VclistCell class as the VclistCell structure
 */
typedef struct sVclistCell VclistCell;

/**
 *  @ingroup Vclist
 *  @author  Nathan Baker
 *  @brief   Atom cell list
 */
struct sVclist {

  Vmem *vmem;  /**< Memory management object for this class */
  Valist *alist;  /**< Original Valist structure for list of atoms */
  Vclist_DomainMode mode;  /**< How the cell list was constructed */
  int npts[VAPBS_DIM];  /**< Hash table grid dimensions */
  int n;  /**< n = nx*nz*ny */
  double max_radius;  /**< Maximum probe radius */
  VclistCell *cells;  /**< Cell array of length thee->n */
  double lower_corner[VAPBS_DIM]; /**< Hash table grid corner */
  double upper_corner[VAPBS_DIM]; /**< Hash table grid corner */
  double spacs[VAPBS_DIM];  /**< Hash table grid spacings */

};

/**
 *  @ingroup Vclist
 *  @brief   Declaration of the Vclist class as the Vclist structure
 */
typedef struct sVclist Vclist;

#if !defined(VINLINE_VCLIST)

    /** @brief   Get number of bytes in this object and its members
     *  @ingroup Vclist
     *  @author  Nathan Baker
     *  @returns Number of bytes allocated for object
     */
    VEXTERNC unsigned long int Vclist_memChk(
            Vclist *thee /**< Object for memory check */
            );

    /**
     * @brief  Get the max probe radius value (in A) the cell list was
     *         constructed with
     * @ingroup Vclist
     * @author Nathan Baker
     * @returns Max probe radius (in A)
     */
    VEXTERNC double Vclist_maxRadius(
            Vclist *thee /**< Cell list object */
            );

#else /* if defined(VINLINE_VCLIST) */

#   define Vclist_memChk(thee) (Vmem_bytes((thee)->vmem))
#   define Vclist_maxRadius(thee) ((thee)->max_radius)

#endif /* if !defined(VINLINE_VCLIST) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vclist: Non-Inlineable methods (vclist.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct the cell list object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 *  @returns Newly allocated Vclist object */
VEXTERNC Vclist* Vclist_ctor(
        Valist *alist, /**< Molecule for cell list queries */
        double max_radius, /**< Max probe radius (&Aring;) to be queried */
        int npts[VAPBS_DIM], /**< Number of in hash table points in each
                              * direction*/
        Vclist_DomainMode mode, /**< Mode to construct table */
        double lower_corner[VAPBS_DIM],  /**< Hash table lower corner for
                                           manual construction (see mode
                                           variable); ignored otherwise */
        double upper_corner[VAPBS_DIM]   /**< Hash table upper corner for
                                           manual construction (see mode
                                           variable); ignored otherwise */
        );

/** @brief   FORTRAN stub to construct the cell list object
 *  @ingroup Vclist
 *  @author  Nathan Baker, Yong Huang
 *  @returns Success enumeration */
VEXTERNC Vrc_Codes Vclist_ctor2(
        Vclist *thee, /**< Memory for Vclist objet */
        Valist *alist, /**< Molecule for cell list queries */
        double max_radius, /**< Max probe radius (&Aring;) to be queried */
        int npts[VAPBS_DIM], /**< Number of in hash table points in each
                              * direction*/
        Vclist_DomainMode mode, /**< Mode to construct table */
        double lower_corner[VAPBS_DIM],  /**< Hash table lower corner for
                                           manual construction (see mode
                                           variable); ignored otherwise */
        double upper_corner[VAPBS_DIM]   /**< Hash table upper corner for
                                           manual construction (see mode
                                           variable); ignored otherwise */
        );

/** @brief   Destroy object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 */
VEXTERNC void Vclist_dtor(
        Vclist **thee /**< Pointer to memory location of object */
        );

/** @brief   FORTRAN stub to destroy object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 */
VEXTERNC void Vclist_dtor2(
        Vclist *thee /**< Pointer to object */
        );

/**
 * @brief  Return cell corresponding to specified position or return VNULL.
 * @ingroup Vclist
 * @author Nathan Baker
 * @returns Pointer to VclistCell object or VNULL if no cell available (away
 * from molecule).
 */
VEXTERNC VclistCell* Vclist_getCell(
        Vclist *thee, /**< Pointer to Vclist cell list */
        double position[VAPBS_DIM] /**< Position to evaluate */
        );

/**
 * @brief  Allocate and construct a cell list cell object
 * @ingroup Vclist
 * @author Nathan Baker
 * @returns Pointer to newly-allocated and constructed object.
 */
VEXTERNC VclistCell* VclistCell_ctor(
        int natoms  /**< Number of atoms associated with this cell */
        );

/**
 * @brief  Construct a cell list object
 * @ingroup  Vclist
 * @author  Nathan Baker, Yong Huang
 * @returns Success enumeration
 */
VEXTERNC Vrc_Codes VclistCell_ctor2(
        VclistCell *thee,  /**< Memory location for object */
        int natoms  /**< Number of atoms associated with this cell */
        );

/** @brief   Destroy object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 */
VEXTERNC void VclistCell_dtor(
        VclistCell **thee /**< Pointer to memory location of object */
        );

/** @brief   FORTRAN stub to destroy object
 *  @ingroup Vclist
 *  @author  Nathan Baker
 */
VEXTERNC void VclistCell_dtor2(
        VclistCell *thee /**< Pointer to object */
        );

#endif    /* ifndef _VCLIST_H_ */
