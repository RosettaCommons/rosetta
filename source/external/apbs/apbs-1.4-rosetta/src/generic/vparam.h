/** @defgroup  Vparam Vparam class
 *  @brief  Reads and assigns charge/radii parameters
 */

/**
 *  @file     vparam.h
 *  @ingroup  Vparam
 *  @brief    Contains declarations for class Vparam
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

#ifndef _VPARAM_H_
#define _VPARAM_H_

#include "apbscfg.h"

#include "maloc/maloc.h"
#if defined(HAVE_MC_H)
#include "mc/mc.h"
#endif

#include "generic/vhal.h"
#include "generic/vunit.h"
#include "generic/vstring.h"

/**
 *  @ingroup  Vparam
 *  @author  Nathan Baker
 *  @brief  AtomData sub-class; stores atom data
 *  @note  The epsilon and radius members of this class refer use the following
 *  formula for calculating the van der Waals energy of atom \f$i\f$
 *  interacting with atom \f$j\f$:
 *  \f[  V_{ij}(r_{ij}) = \epsilon_{ij} \left[
 *       \left( \frac{\sigma_{ij}}{r_{ij}} \right)^{12} - 2
 *       \left( \frac{\sigma_{ij}}{r_{ij}} \right)^{6} \right]
 *  \f]
 *  where \f$\epsilon_{ij} = \sqrt{\epsilon_i \epsilon_j}\f$ is the well-depth
 *  (in the desired energy units), \f$r_{ij}\f$ is the distance between atoms
 *  \f$i\f$ and \f$j\f$, and \f$\sigma_{ij} = \sigma_i + \sigma_j\f$ is the sum
 *  of the van der Waals radii.
 */
struct sVparam_AtomData {
    char atomName[VMAX_ARGLEN];  /**< Atom name */
    char resName[VMAX_ARGLEN];  /**< Residue name */
    double charge;  /**< Atom charge (in e) */
    double radius;  /**< Atom VdW radius (\f$\sigma_i\f$ above; in &Aring;) */
    double epsilon;  /**< Atom VdW well depth (\f$\epsilon_i\f$ above; in
                      * kJ/mol) */
};

/**
 *  @ingroup Vparam
 *  @brief   Declaration of the Vparam_AtomData class as the sVparam_AtomData
 *           structure
 */
typedef struct sVparam_AtomData Vparam_AtomData;

/**
 *  @struct  Vparam_ResData
 *  @ingroup  Vparam
 *  @author  Nathan Baker
 *  @brief  ResData sub-class; stores residue data
 */
struct Vparam_ResData {
    Vmem *vmem;  /**<  Pointer to memory manager from Vparam master class */
    char name[VMAX_ARGLEN]; /**< Residue name */
    int nAtomData;  /**<  Number of Vparam_AtomData objects associated with
                     * this object */
    Vparam_AtomData *atomData;  /**<  Array of Vparam_AtomData natom objects */
};

/** @typedef Vparam_ResData
 *  @ingroup Vparam
 *  @brief   Declaration of the Vparam_ResData class as the Vparam_ResData
 *           structure
 */
typedef struct Vparam_ResData Vparam_ResData;

/**
 *  @struct  Vparam
 *  @ingroup  Vparam
 *  @author  Nathan Baker
 *  @brief  Reads and assigns charge/radii parameters
 */
struct Vparam {

  Vmem *vmem;  /**< Memory management object for this class */
  int nResData;  /**< Number of Vparam_ResData objects associated with
                  * this object */
  Vparam_ResData *resData;  /**< Array of nResData Vparam_ResData objects */
};

/** @typedef Vparam
 *  @ingroup Vparam
 *  @brief   Declaration of the Vparam class as the Vparam structure
 */
typedef struct Vparam Vparam;

/* ///////////////////////////////////////////////////////////////////////////
// Class Vparam: Inlineable methods (vparam.c)
/////////////////////////////////////////////////////////////////////////// */

#if !defined(VINLINE_VPARAM)

    /** @brief   Get number of bytes in this object and its members
     *  @ingroup Vparam
     *  @author  Nathan Baker
     *  @param   thee  Vparam object
     *  @returns Number of bytes allocated for object
     */
    VEXTERNC unsigned long int Vparam_memChk(Vparam *thee);

#else /* if defined(VINLINE_VPARAM) */

#   define Vparam_memChk(thee) (Vmem_bytes((thee)->vmem))

#endif /* if !defined(VINLINE_VPARAM) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vparam: Non-Inlineable methods (vparam.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct the object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @returns Newly allocated object */
VEXTERNC Vparam_AtomData* Vparam_AtomData_ctor();

/** @brief   FORTRAN stub to construct the object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Allocated memory
 *  @returns 1 if successful, 0 otherwise */
VEXTERNC int Vparam_AtomData_ctor2(Vparam_AtomData *thee);

/** @brief   Destroy object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of object */
VEXTERNC void Vparam_AtomData_dtor(Vparam_AtomData **thee);

/** @brief   FORTRAN stub to destroy object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Pointer to object */
VEXTERNC void Vparam_AtomData_dtor2(Vparam_AtomData *thee);

/**
 * @brief  Copy current atom object to destination
 * @ingroup  Vparam
 * @author  Nathan Baker
 * @param  thee  Pointer to source object
 * @param  dest  Pointer to destination object
 */
VEXTERNC void Vparam_AtomData_copyTo(Vparam_AtomData *thee,
  Vparam_AtomData *dest);

/**
 * @brief  Copy current residue object to destination
 * @ingroup  Vparam
 * @author  Todd Dolinsky
 * @param  thee  Pointer to source object
 * @param  dest  Pointer to destination object
 */
VEXTERNC void Vparam_ResData_copyTo(Vparam_ResData *thee,
  Vparam_ResData *dest);

/**
 * @brief  Copy current atom object from another
 * @ingroup  Vparam
 * @author  Nathan Baker
 * @param  thee  Pointer to destination object
 * @param  src  Pointer to source object
 */
VEXTERNC void Vparam_AtomData_copyFrom(Vparam_AtomData *thee,
  Vparam_AtomData *src);

/** @brief   Construct the object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param  mem  Memory object of Vparam master class
 *  @returns Newly allocated object */
VEXTERNC Vparam_ResData* Vparam_ResData_ctor(Vmem *mem);

/** @brief   FORTRAN stub to construct the object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Allocated memory
 *  @param  mem  Memory object of Vparam master class
 *  @returns 1 if successful, 0 otherwise */
VEXTERNC int Vparam_ResData_ctor2(Vparam_ResData *thee, Vmem *mem);

/** @brief   Destroy object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of object */
VEXTERNC void Vparam_ResData_dtor(Vparam_ResData **thee);

/** @brief   FORTRAN stub to destroy object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Pointer to object */
VEXTERNC void Vparam_ResData_dtor2(Vparam_ResData *thee);

/** @brief   Construct the object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @returns Newly allocated Vparam object */
VEXTERNC Vparam* Vparam_ctor();

/** @brief   FORTRAN stub to construct the object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Allocated Vparam memory
 *  @returns 1 if successful, 0 otherwise */
VEXTERNC int Vparam_ctor2(Vparam *thee);

/** @brief   Destroy object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of object */
VEXTERNC void Vparam_dtor(Vparam **thee);

/** @brief   FORTRAN stub to destroy object
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Pointer to object */
VEXTERNC void Vparam_dtor2(Vparam *thee);

/** @brief  Get residue data
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Vparam object
 *  @param   resName  Residue name
 *  @returns  Pointer to the desired residue object or VNULL if residue not
 *  found
 *  @note  Some method to initialize the database must be called before this
 *  method (e.g., @see Vparam_readFlatFile)
 */
VEXTERNC Vparam_ResData* Vparam_getResData(Vparam *thee,
  char resName[VMAX_ARGLEN]);

/** @brief  Get atom data
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @param   thee  Vparam object
 *  @param   resName  Residue name
 *  @param   atomName  Atom name
 *  @returns  Pointer to the desired atom object or VNULL if residue not
 *  found
 *  @note  Some method to initialize the database must be called before this
 *  method (e.g., @see Vparam_readFlatFile)
 */
VEXTERNC Vparam_AtomData* Vparam_getAtomData(Vparam *thee,
  char resName[VMAX_ARGLEN], char atomName[VMAX_ARGLEN]);

/** @brief  Read a flat-file format parameter database
 * @ingroup Vparam
 * @author  Nathan Baker
 * @param  thee Vparam object
 * @param   iodev  Input device type (FILE/BUFF/UNIX/INET)
 * @param   iofmt  Input device format (ASCII/XDR)
 * @param   thost  Input hostname (for sockets)
 * @param   fname  Input FILE/BUFF/UNIX/INET name
 * (see note below for format)
 * @returns 1 if successful, 0 otherwise
 * @note  The database file should have the following format:
 * <pre>
 * RESIDUE ATOM CHARGE RADIUS EPSILON
 * </pre>
 * where RESIDUE is the residue name string, ATOM is the atom name string,
 * CHARGE is the charge in e, RADIUS is the van der Waals radius
 * (\f$\sigma_i\f$) in &Aring;, and EPSILON is the van der Waals well-depth
 * (\f$\epsilon_i\f$) in kJ/mol.  See the Vparam structure documentation for
 * the precise definitions of \f$\sigma_i\f$ and \f$\epsilon_i\f$.
 *
 * ASCII-format flat files are provided with the APBS source code:
 * <dl>
 *   <dt> tools/conversion/vparam-amber-parm94.dat
 *   <dd> AMBER parm94 parameters
 *   <dt> tools/conversion/vparam-charmm-par_all27.dat
 *   <dd> CHARMM par_all27_prot_na parameters
 * </dl>
 * */
VEXTERNC int Vparam_readFlatFile(Vparam *thee, const char *iodev,
  const char *iofmt, const char *thost, const char *fname);

/** @brief  Read an XML format parameter database
 * @ingroup Vparam
 * @author  Todd Dolinsky
 * @param  thee Vparam object
 * @param   iodev  Input device type (FILE/BUFF/UNIX/INET)
 * @param   iofmt  Input device format (ASCII/XDR)
 * @param   thost  Input hostname (for sockets)
 * @param   fname  Input FILE/BUFF/UNIX/INET name
 * @returns 1 if successful, 0 otherwise
 * */
VEXTERNC int Vparam_readXMLFile(Vparam *thee, const char *iodev,
  const char *iofmt, const char *thost, const char *fname);

#endif    /* ifndef _VPARAM_H_ */
