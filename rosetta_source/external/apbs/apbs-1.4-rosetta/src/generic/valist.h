/** @defgroup Valist Valist class
 *  @brief    Container class for list of atom objects
 */

/**
 *  @file     valist.h
 *  @ingroup  Valist
 *  @brief    Contains declarations for class Valist
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

#ifndef _VALIST_H_
#define _VALIST_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vatom.h"
#include "generic/vparam.h"

/**
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @brief   Container class for list of atom objects
 */
struct sValist {

  int number;         /**< Number of atoms in list */
  double center[3];   /**< Molecule center (xmin - xmax)/2, etc.*/
  double mincrd[3];   /**< Minimum coordinates */
  double maxcrd[3];   /**< Maximum coordinates */
  double maxrad;      /**< Maximum radius */
  double charge;      /**< Net charge */
  Vatom *atoms;       /**< Atom list */
  Vmem *vmem;         /**< Memory management object */

};

/**
 *  @ingroup Valist
 *  @brief Declaration of the Valist class as the Valist structure
 */
typedef struct sValist Valist;

#if !defined(VINLINE_VATOM)

/**
 * @brief   Get actual array of atom objects from the list
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Array of atom objects
 */
VEXTERNC Vatom* Valist_getAtomList(
        Valist *thee /**< Atom list object */
        );

/** @brief   Get x-coordinate of molecule center
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  X-coordinate of molecule center
 */
VEXTERNC double Valist_getCenterX(
        Valist *thee /**< Atom list object */
        );

/** @brief   Get y-coordinate of molecule center
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Y-coordinate of molecule center
 */
VEXTERNC double Valist_getCenterY(
        Valist *thee /**< Atom list object */
        );

/** @brief   Get z-coordinate of molecule center
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Z-coordinate of molecule center
 */
VEXTERNC double Valist_getCenterZ(
        Valist *thee /**< Atom list object */
        );

/** @brief   Get number of atoms in the list
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Number of atoms in list
 */
VEXTERNC int Valist_getNumberAtoms(
        Valist *thee /**< Atom list object */
        );

/** @brief   Get pointer to particular atom in list
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Pointer to atom object i
 */
VEXTERNC Vatom* Valist_getAtom(
        Valist *thee, /**< Atom list object */
        int i /**< Index of atom in list */
        );

/** @brief   Get total memory allocated for this object and its members
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @return  Total memory in bytes
 */
VEXTERNC unsigned long int Valist_memChk(
        Valist *thee /**< Atom list object */
        );

#else /* if defined(VINLINE_VATOM) */
#   define Valist_getAtomList(thee) ((thee)->atoms)
#   define Valist_getNumberAtoms(thee) ((thee)->number)
#   define Valist_getAtom(thee, i) (&((thee)->atoms[i]))
#   define Valist_memChk(thee) (Vmem_bytes((thee)->vmem))
#   define Valist_getCenterX(thee) ((thee)->center[0])
#   define Valist_getCenterY(thee) ((thee)->center[1])
#   define Valist_getCenterZ(thee) ((thee)->center[2])
#endif /* if !defined(VINLINE_VATOM) */

/** @brief   Construct the atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 *  @returns Pointer to newly allocated (empty) atom list
 */
VEXTERNC Valist* Valist_ctor();

/** @brief   FORTRAN stub to construct the atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker, Yong Huang
 *  @returns Success enumeration
 */
VEXTERNC Vrc_Codes Valist_ctor2(
        Valist *thee /**< Storage for new atom list */
        );

/** @brief   Destroys atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 */
VEXTERNC void Valist_dtor(
        Valist **thee /**< Pointer to storage for atom list */
        );

/** @brief   FORTRAN stub to destroy atom list object
 *  @ingroup Valist
 *  @author  Nathan Baker
 */
VEXTERNC void Valist_dtor2(
        Valist *thee /**< Pointer to atom list object */
        );

/**
 * @brief  Fill atom list with information from a PQR file
 * @ingroup Valist
 * @author  Nathan Baker, Yong Huang
 * @returns	Success enumeration
 * @note  \li A PQR file has PDB structure with charge and radius in the last
 *            two columns instead of weight and occupancy
 *        \li We don't actually respect PDB format; instead recognize
 *            whitespace- or tab-delimited fields which allows us to deal with
 *            structures with coordinates > 999 or < -999.
 */
VEXTERNC Vrc_Codes Valist_readPQR(
        Valist *thee, /**< Atom list object */
        Vparam *param, /**< A pre-initialized parameter object */
        Vio *sock /**< Socket reading for reading PQR file */
        );

/**
 * @brief  Fill atom list with information from a PDB file
 * @ingroup Valist
 * @author  Nathan Baker, Todd Dolinsky, Yong Huang
 * @returns	Success enumeration
 * @note  We don't actually respect PDB format; instead recognize whitespace-
 * or tab-delimited fields which allows us to deal with structures with
 * coordinates > 999 or < -999.
 */
VEXTERNC Vrc_Codes Valist_readPDB(
        Valist *thee, /**< Atom list object */
        Vparam *param, /**< A pre-initialized parameter object */
        Vio *sock /**< Socket read for reading PDB file */
        );

/**
 * @brief  Fill atom list with information from an XML file
 * @ingroup Valist
 * @author  Todd Dolinsky, Yong Huang
 * @returns Success enumeration
 * @note  \li The XML file must adhere to some guidelines, notably the
 *            presence of an &lt;atom&gt; tag with all other useful information
 *            (x, y, z, charge, and radius) as nested elements.
 */
VEXTERNC Vrc_Codes Valist_readXML(
        Valist *thee, /**< Atom list object */
        Vparam *param, /**< A pre-initialized parameter object */
        Vio *sock /**< Socket reading for reading PQR file */
        );

/**
 * @brief   Load up Valist with various statistics
 * @ingroup Valist
 * @author  Nathan Baker, Yong Huang
 * @returns	Success enumeration
 */
VEXTERNC Vrc_Codes Valist_getStatistics(Valist *thee);


#endif /* ifndef _VALIST_H_ */
