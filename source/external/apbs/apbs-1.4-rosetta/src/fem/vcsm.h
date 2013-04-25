/** @defgroup Vcsm Vcsm class
 *  @brief  A charge-simplex map for evaluating integrals of delta functions
 *          in a finite element setting
 */

/**
 *  @file      vcsm.h
 *  @brief     Contains declarations for the Vcsm class
 *  @ingroup   Vcsm
 *  @version   $Id$
 *  @author    Nathan A. Baker
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

#ifndef _VCSM_H_
#define _VCSM_H_

#include "apbscfg.h"

#include "maloc/maloc.h"
#include "mc/mc.h"

#include "generic/vhal.h"
#include "generic/valist.h"

/** @brief   External function for FEtk Gem class to use during mesh refinement
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 */
VEXTERNC void Gem_setExternalUpdateFunction(
        Gem *thee, /**< The FEtk geometry manager */
        void (*externalUpdate)(SS **simps, int num) /**< Function pointer for
                                                      call during mesh
                                                      refinement */
        );

/** @brief   Charge-simplex map class
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 */
struct sVcsm {

  Valist *alist;      /**< Atom (charge) list */
  int natom;          /**< Size of thee->alist; redundant, but useful for
                       * convenience */
  Gem *gm;            /**< Grid manager (container class for master vertex
                       * and simplex lists as well as prolongation
                       * operator for updating after refinement ) */
  int **sqm;          /**< The map which gives the list charges associated with
                       * each simplex in gm->simplices.  The indices of
                       * the first dimension are associated with the
                       * simplex ID's in Vgm.  Each charge list (second
                       * dimension) contains entries corresponding to
                       * indicies in thee->alist with lengths given in
                       * thee->nsqm */
  int *nsqm;          /**< The length of the charge lists in thee->sqm */
  int nsimp;          /**< The _currently used) length of sqm, nsqm -- may not
                       * always be up-to-date with Gem */
  int msimp;          /**< The maximum number of entries that can be
                       * accomodated by sqm or nsqm  -- saves on realloc's */
  int **qsm;          /**< The inverse of sqm; the list of simplices
                       * associated with a given charge */
  int *nqsm;          /**< The length of the simplex lists in thee->qsm */
  int initFlag;       /**< Indicates whether the maps have been initialized
                       * yet */
  Vmem *vmem;         /**< Memory management object */

};

/**
 *  @ingroup Vcsm
 *  @brief   Declaration of the Vcsm class as the Vcsm structure
 */
typedef struct sVcsm Vcsm;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Inlineable methods (vcsm.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VCSM)

    /** @brief   Get atom list
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Pointer to Valist atom list
     */
    VEXTERNC Valist* Vcsm_getValist(
            Vcsm *thee /**< The Vcsm object */
            );

    /** @brief   Get number of atoms associated with a simplex
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Number of atoms associated with a simplex
     */
    VEXTERNC int Vcsm_getNumberAtoms(
            Vcsm *thee,  /**< The Vcsm object */
            int isimp  /**< Simplex ID */
            );

    /** @brief   Get particular atom associated with a simplex
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Array of atoms associated with a simplex
     */
    VEXTERNC Vatom* Vcsm_getAtom(
            Vcsm *thee,  /**< The Vcsm object */
            int iatom,  /**< Index of atom in Vcsm list ofr this simplex */
            int isimp  /**< Simplex ID */
            );

    /** @brief   Get ID of particular atom in a simplex
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Index of atom in Valist object
     */
    VEXTERNC int Vcsm_getAtomIndex(
            Vcsm *thee,  /**< The Vcsm object */
            int iatom,  /**< Index of atom in Vcsm list for this simplex */
            int isimp  /**< Simplex ID */
            );

    /** @brief   Get number of simplices associated with an atom
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Number of simplices associated with an atom
     */
    VEXTERNC int Vcsm_getNumberSimplices(
            Vcsm *thee,  /**< The Vcsm object */
            int iatom  /**< The Valist atom index */
            );

    /** @brief   Get particular simplex associated with an atom
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Pointer to simplex object
     */
    VEXTERNC SS* Vcsm_getSimplex(
            Vcsm *thee,  /**< The Vcsm object */
            int isimp,  /**< Index of simplex in Vcsm list */
            int iatom  /**< Valist atom index */
            );

    /** @brief   Get index particular simplex associated with an atom
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  Gem index of specified simplex
     */
    VEXTERNC int Vcsm_getSimplexIndex(
            Vcsm *thee,  /**< The Vcsm object */
            int isimp,  /**< Index of simplex in Vcsm list */
            int iatom  /**< Index of atom in Valist */
            );

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vcsm
     *  @author  Nathan Baker
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC unsigned long int Vcsm_memChk(
            Vcsm *thee /**< The Vcsm object */
            );

#else /* if defined(VINLINE_VCSM) */
#   define Vcsm_getValist(thee) ((thee)->alist)
#   define Vcsm_getNumberAtoms(thee, isimp) ((thee)->nsqm[isimp])
#   define Vcsm_getAtom(thee, iatom, isimp) (Valist_getAtom((thee)->alist, ((thee)->sqm)[isimp][iatom]))
#   define Vcsm_getAtomIndex(thee, iatom, isimp) (((thee)->sqm)[isimp][iatom])
#   define Vcsm_getNumberSimplices(thee, iatom) (((thee)->nqsm)[iatom])
#   define Vcsm_getSimplex(thee, isimp, iatom) (Gem_SS((thee)->gm, ((thee)->qsm)[iatom][isimp]))
#   define Vcsm_getSimplexIndex(thee, isimp, iatom) (((thee)->qsm)[iatom][isimp])
#   define Vcsm_memChk(thee) (Vmem_bytes((thee)->vmem))
#endif /* if !defined(VINLINE_VCSM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vcsm: Non-Inlineable methods (vcsm.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct Vcsm object
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @note    \li The initial mesh must be sufficiently coarse for the assignment
 *             procedures to be efficient
 *           \li The map is not built until Vcsm_init is called
 *  @return  Pointer to newly allocated Vcsm object
 */
VEXTERNC Vcsm* Vcsm_ctor(
        Valist *alist,  /**< List of atoms */
        Gem *gm  /**< FEtk geometry manager defining the mesh */
        );

/** @brief   FORTRAN stub to construct Vcsm object
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @note    \li The initial mesh must be sufficiently coarse for the assignment
 *             procedures to be efficient
 *           \li The map is not built until Vcsm_init is called
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vcsm_ctor2(
        Vcsm *thee,  /**< The Vcsm object */
        Valist *alist,  /**< The list of atoms */
        Gem *gm  /**< The FEtk geometry manager defining the mesh */
        );

/** @brief   Destroy Vcsm object
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 */
VEXTERNC void Vcsm_dtor(
        Vcsm **thee  /**< Pointer to memory location for Vcsm object */
        );

/** @brief   FORTRAN stub to destroy Vcsm object
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 */
VEXTERNC void Vcsm_dtor2(
        Vcsm *thee /**< Pointer to Vcsm object */
        );

/** @brief   Initialize charge-simplex map with mesh and atom data
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @note    The initial mesh must be sufficiently coarse for the assignment
 *            procedures to be efficient
 */
VEXTERNC void Vcsm_init(
        Vcsm *thee /**< The Vcsm object */
        );

/** @brief   Update the charge-simplex and simplex-charge maps after
 *           refinement
 *  @ingroup Vcsm
 *  @author  Nathan Baker
 *  @return  1 if successful, 0 otherwise
 */
VEXTERNC int Vcsm_update(
        Vcsm *thee, /**< The Vcsm object */
        SS **simps, /**< List of pointer to newly created (by refinement)
                      simplex objects.  The first simplex is expected to be
                      derived from the parent simplex and therefore have the
                      same ID.  The remaining simplices are the children and
                      should represent new entries in the charge-simplex map. */
        int num /**< Number of simplices in simps list */
        );

#endif /* ifndef _VCSM_H_ */
