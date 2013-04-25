/** @defgroup Vatom Vatom class
 *  @brief  Atom class for interfacing APBS with PDB files
 */

/**
 *  @file     vatom.h
 *  @ingroup  Vatom
 *  @brief    Contains declarations for class Vatom
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

#ifndef _VATOM_H_
#define _VATOM_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"

/**
 *  @ingroup Vatom
 *  @def VMAX_RECLEN
 *  @author  Nathan Baker, David Gohara, Mike Schneiders
 *  @brief   Residue name length
 */
#define VMAX_RECLEN		   64

/**
 *  @ingroup Vatom
 *  @author  Nathan Baker, David Gohara, Mike Schneiders
 *  @brief   Contains public data members for Vatom class/module
 */
struct sVatom {

    double position[3];  /**< Atomic position */
    double radius;  /**< Atomic radius   */
    double charge;  /**< Atomic charge   */
    double partID;  /**< Partition value for assigning atoms to particular
                     * processors and/or partitions   */
    double epsilon; /**< Epsilon value for WCA calculations */

    int id;  /**< Atomic ID; this should be a unique non-negative integer
              * assigned based on the index of the atom in a Valist atom
              * array */

    char resName[VMAX_RECLEN]; /**< Residue name from PDB/PQR file */
    char atomName[VMAX_RECLEN]; /**< Atom name from PDB/PDR file */

#if defined(WITH_TINKER)

    double dipole[3];          /**< Permanent dipole */
    double quadrupole[9];      /**< Permanent quadrupole */
    double inducedDipole[3];   /**< Induced dipole */
    double nlInducedDipole[3];  /**< Non-local induced dipole */

#endif /* if defined(WITH_TINKER) */
};

/**
 *  @ingroup Vatom
 *  @brief   Declaration of the Vatom class as the Vatom structure
 */
typedef struct sVatom Vatom;

#if !defined(VINLINE_VATOM)

    /** @brief   Get atomic position
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee  Vatom object
     *  @returns Pointer to 3*double array of atomic coordinates (in &Aring;)
     */
    VEXTERNC double* Vatom_getPosition(Vatom *thee);

    /** @brief   Set atomic radius
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @param   radius  Atomic radius (in &Aring;)
     */
    VEXTERNC void    Vatom_setRadius(Vatom *thee, double radius);

    /** @brief   Get atomic position
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee  Vatom object
     *  @returns Atomic radius (in &Aring;)
     */
    VEXTERNC double  Vatom_getRadius(Vatom *thee);

    /** @brief   Set partition ID
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @param   partID  Partition ID; a negative value means this atom is not
     *                   assigned to any partition
     */
    VEXTERNC void    Vatom_setPartID(Vatom *thee, int partID);

    /** @brief   Get partition ID
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @return  Partition ID; a negative value means this atom is not
     *           assigned to any partition
     */
    VEXTERNC double     Vatom_getPartID(Vatom *thee);

    /** @brief   Set atom ID
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @param   id  Unique non-negative number
     */
    VEXTERNC void Vatom_setAtomID(Vatom *thee, int id);

    /** @brief   Get atom ID
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @return  Unique non-negative number
     */
    VEXTERNC double Vatom_getAtomID(Vatom *thee);

    /** @brief   Set atomic charge
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @param   charge  Atom partial charge (in e)
     */
    VEXTERNC void    Vatom_setCharge(Vatom *thee, double charge);

    /** @brief   Get atomic charge
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee    Vatom object
     *  @return  Atom partial charge (in e)
     */
    VEXTERNC double  Vatom_getCharge(Vatom *thee);

    /** @brief   Set atomic epsilon
    *  @ingroup Vatom
    *  @author  David Gohara
    *  @param   thee    Vatom object
    *  @param   epsilon  Atomic epsilon (in &Aring;)
    */
    VEXTERNC void    Vatom_setEpsilon(Vatom *thee, double epsilon);

    /** @brief   Get atomic epsilon
    *  @ingroup Vatom
    *  @author  David Gohara
    *  @param   thee  Vatom object
    *  @returns Atomic epsilon (in &Aring;)
    */
    VEXTERNC double  Vatom_getEpsilon(Vatom *thee);

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vatom
     *  @author  Nathan Baker
     *  @param   thee  Vpmg object
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC unsigned long int Vatom_memChk(Vatom *thee);

#else /* if defined(VINLINE_VATOM) */
#   define Vatom_getPosition(thee) ((thee)->position)
#   define Vatom_setRadius(thee, tRadius) ((thee)->radius = (tRadius))
#   define Vatom_getRadius(thee) ((thee)->radius)
#   define Vatom_setPartID(thee, tpartID) ((thee)->partID = (double)(tpartID))
#   define Vatom_getPartID(thee) ((thee)->partID)
#   define Vatom_setAtomID(thee, tatomID) ((thee)->id = (tatomID))
#   define Vatom_getAtomID(thee) ((thee)->id)
#   define Vatom_setCharge(thee, tCharge) ((thee)->charge = (tCharge))
#   define Vatom_getCharge(thee) ((thee)->charge)
#   define Vatom_setEpsilon(thee, tEpsilon) ((thee)->epsilon = (tEpsilon))
#   define Vatom_getEpsilon(thee) ((thee)->epsilon)
#   define Vatom_memChk(thee) (sizeof(Vatom))
#endif /* if !defined(VINLINE_VATOM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vatom: Non-Inlineable methods (vatom.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Set residue name
*  @ingroup Vatom
*  @author  Jason Wagoner
*  @param   thee    Vatom object
*  @param   resName Residue Name
*/
VEXTERNC void    Vatom_setResName(Vatom *thee, char resName[VMAX_RECLEN]);

/** @brief   Set atom name
*  @ingroup Vatom
*  @author  Jason Wagoner
*/
VEXTERNC void    Vatom_setAtomName(
        Vatom *thee,  /**< Vatom object */
        char atomName[VMAX_RECLEN]  /**< Atom name */
        );

/** @brief   Retrieve residue name
*  @ingroup Vatom
*  @author  Jason Wagoner
*  @param   thee    Vatom object
*  @param   resName Residue Name
*/
VEXTERNC void    Vatom_getResName(Vatom *thee, char resName[VMAX_RECLEN]);

/** @brief   Retrieve atom name
*  @ingroup Vatom
*  @author  Jason Wagoner
*/
VEXTERNC void   Vatom_getAtomName(
        Vatom *thee, /**< Vatom object */
        char atomName[VMAX_RECLEN] /**< Atom name */
        );

/** @brief   Constructor for the Vatom class
 *  @author  Nathan Baker
 *  @ingroup Vatom
 *  @returns Pointer to newly allocated Vatom object
 */
VEXTERNC Vatom* Vatom_ctor();

/** @brief   FORTRAN stub constructor for the Vatom class
 *  @author  Nathan Baker
 *  @ingroup Vatom
 *  @param   thee Pointer to Vatom allocated memory location
 *  @returns 1 if succesful, 0 otherwise
 */
VEXTERNC int     Vatom_ctor2(Vatom *thee);

/** @brief   Object destructor
 *  @ingroup Vatom
 *  @author  Nathan Baker
 *  @param   thee   Pointer to memory location of object to be destroyed
 */
VEXTERNC void    Vatom_dtor(Vatom **thee);

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vatom
 *  @author  Nathan Baker
 *  @param   thee   Pointer to object to be destroyed
 */
VEXTERNC void    Vatom_dtor2(Vatom *thee);

/** @brief   Set the atomic position
 *  @ingroup Vatom
 *  @author  Nathan Baker
 *  @param   thee   Vatom object to be modified
 *  @param   position  Coordinates (in &Aring;)
 */
VEXTERNC void   Vatom_setPosition(Vatom *thee, double position[3]);

/**
 * @brief  Copy information to another atom
 * @ingroup  Vatom
 * @author  Nathan Baker
 * @param  thee Source for atom information
 * @param  dest Destination for atom information
 */
VEXTERNC void Vatom_copyTo(Vatom *thee, Vatom *dest);

/**
 * @brief  Copy information to another atom
 * @ingroup  Vatom
 * @author  Nathan Baker
 * @param  thee Destination for atom information
 * @param  src Source for atom information
 */
VEXTERNC void Vatom_copyFrom(Vatom *thee, Vatom *src);

#if defined(WITH_TINKER)

/** @brief   Set the induced dipole moment
 *  @ingroup Vatom
 *  @author  Michael Schnieders
 *  @param   thee   Vatom object to be modified
 *  @param   inducedDipole Induced dipole moment (e*A)
 */
VEXTERNC void   Vatom_setInducedDipole(Vatom *thee,
                                       double inducedDipole[3]);

/** @brief   Set the non-local induced dipole moment
 *  @ingroup Vatom
 *  @author  Michael Schnieders
 *  @param   thee   Vatom object to be modified
 *  @param   nlInducedDipole Induced dipole moment (e*A)
 */
VEXTERNC void   Vatom_setNLInducedDipole(Vatom *thee,
                                       double nlInducedDipole[3]);

/** @brief   Set the permanent dipole moment
 *  @ingroup Vatom
 *  @author  Michael Schnieders
 *  @param   thee   Vatom object to be modified
 *  @param   dipole Permanent dipole moment
 */
VEXTERNC void   Vatom_setDipole(Vatom *thee, double dipole[3]);

/** @brief   Set the permanent quadrupole moment
 *  @ingroup Vatom
 *  @author  Michael Schnieders
 *  @param   thee   Vatom object to be modified
 *  @param   quadrupole Permanent quadrupole moment
 */
VEXTERNC void   Vatom_setQuadrupole(Vatom *thee, double quadrupole[9]);

/** @brief   Get permanent dipole
 *  @ingroup Vatom
 *  @author  Michael Schnieders
 *  @param   thee    Vatom object
 */
VEXTERNC double*  Vatom_getDipole(Vatom *thee);

/** @brief   Get permanent quadrupole
 *  @ingroup Vatom
 *  @author  Michael Schnieders
 *  @param   thee    Vatom object
 */
VEXTERNC double*  Vatom_getQuadrupole(Vatom *thee);

/** @brief   Get induced dipole
 *  @ingroup Vatom
 *  @author  Michael Schnieders
 *  @param   thee    Vatom object
 */
VEXTERNC double*  Vatom_getInducedDipole(Vatom *thee);

/** @brief   Get non-local induced dipole
 *  @ingroup Vatom
 *  @author  Michael Schnieders
 *  @param   thee    Vatom object
 */
VEXTERNC double*  Vatom_getNLInducedDipole(Vatom *thee);
#endif /* if defined(WITH_TINKER) */

#endif /* ifndef _VATOM_H_ */
