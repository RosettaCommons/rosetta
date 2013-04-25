/** @defgroup PBEparm PBEparm class
 *  @brief    Parameter structure for PBE variables independent of solver
 */

/**
 *  @file     pbeparm.h
 *  @ingroup  PBEparm
 *  @brief    Contains declarations for class PBEparm
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

#ifndef _PBEPARM_H_
#define _PBEPARM_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vstring.h"

/** @brief   Number of things that can be written out in a single calculation
 *  @ingroup PBEparm
 */
#define PBEPARM_MAXWRITE 20

/**
 * @ingroup PBEparm
 * @brief  Define energy calculation enumeration
 */
enum ePBEparm_calcEnergy {
    PCE_NO=0, /**< Do not perform energy calculation */
    PCE_TOTAL=1, /**< Calculate total energy only */
    PCE_COMPS=2 /**< Calculate per-atom energy components */
};

/**
 * @ingroup PBEparm
 * @brief  Define ePBEparm_calcEnergy enumeration as PBEparm_calcEnergy
 */
typedef enum ePBEparm_calcEnergy PBEparm_calcEnergy;

/**
 * @ingroup PBEparm
 * @brief  Define force calculation enumeration
 */
enum ePBEparm_calcForce {
    PCF_NO=0, /**< Do not perform force calculation */
    PCF_TOTAL=1, /**< Calculate total force only */
    PCF_COMPS=2 /**< Calculate per-atom force components */
};

/**
 * @ingroup PBEparm
 * @brief  Define ePBEparm_calcForce enumeration as PBEparm_calcForce
 */
typedef enum ePBEparm_calcForce PBEparm_calcForce;

/**
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @brief   Parameter structure for PBE variables from input files
 *  @note    If you add/delete/change something in this class, the member
 *           functions -- especially PBEparm_copy -- must be modified
 *           accordingly
 */
struct sPBEparm {

    int molid;  /**< Molecule ID to perform calculation on */
    int setmolid;  /**< Flag, @see molid */
    int useDielMap;  /**< Indicates whether we use external
                      * dielectric maps (note plural) */
    int dielMapID;  /**< Dielectric map ID (if used) */
    int useKappaMap;  /**< Indicates whether we use an external
                       * kappa map */
    int kappaMapID;  /**< Kappa map ID (if used) */
    int usePotMap;  /**< Indicates whether we use an external
                       * kappa map */
    int potMapID;  /**< Kappa map ID (if used) */

    int useChargeMap;  /**< Indicates whether we use an external
                        * charge distribution map */
    int chargeMapID;  /**< Charge distribution map ID (if used) */
    Vhal_PBEType pbetype;  /**< Which version of the PBE are we solving? */
    int setpbetype;  /**< Flag, @see pbetype */
    Vbcfl bcfl;  /**< Boundary condition method */
    int setbcfl;  /**< Flag, @see bcfl */
    int nion;  /**< Number of counterion species */
    int setnion;  /**< Flag, @see nion */
    double ionq[MAXION];  /**< Counterion charges (in e) */
    double ionc[MAXION];  /**< Counterion concentrations (in M) */
    double ionr[MAXION];  /**< Counterion radii (in A) */
    int setion[MAXION];  /**< Flag, @see ionq */
    double pdie;  /**< Solute dielectric */
    int setpdie;  /**< Flag, @see pdie */
    double sdens; /**< Vacc sphere density */
    int setsdens; /**< Flag, @see sdens */
    double sdie;  /**< Solvent dielectric */
    int setsdie;  /**< Flag, @see sdie */
    Vsurf_Meth srfm;  /**< Surface calculation method */
    int setsrfm;  /**< Flag, @see srfm */
    double srad;  /**< Solvent radius */
    int setsrad;  /**< Flag, @see srad */
    double swin;  /**< Cubic spline window */
    int setswin;  /**< Flag, @see swin */
    double temp;  /**< Temperature (in K) */
    int settemp;  /**< Flag, @see temp */

    double smsize; /**< SMPBE size */
    int setsmsize; /**< Flag, @see temp */

    double smvolume; /**< SMPBE size */
    int setsmvolume; /**< Flag, @see temp */

    PBEparm_calcEnergy calcenergy;  /**< Energy calculation flag */
    int setcalcenergy;  /**< Flag, @see calcenergy */
    PBEparm_calcForce calcforce;  /**< Atomic forces calculation */
    int setcalcforce;  /**< Flag, @see calcforce */

    /*----------------------------------------------------------------*/
    /* Added by Michael Grabe                                         */
    /*----------------------------------------------------------------*/

    double zmem;               /**< z value of membrane bottom */
    int setzmem;               /**< Flag */
    double Lmem;               /**< membrane width */
    int setLmem;               /**< Flag */
    double mdie;               /**< membrane dielectric constant */
    int setmdie;               /**< Flag */
    double memv;               /**< Membrane potential */
    int setmemv;               /**< Flag */

    /*----------------------------------------------------------------*/

    int numwrite;  /**< Number of write statements encountered */
    char writestem[PBEPARM_MAXWRITE][VMAX_ARGLEN]; /**< File stem to write
                                                    * data to */
    Vdata_Type writetype[PBEPARM_MAXWRITE];  /**< What data to write */
    Vdata_Format writefmt[PBEPARM_MAXWRITE];  /**< File format to write data
                                               * in */
    int writemat;  /**< Write out the operator matrix?
                    * \li 0 => no
                    * \li 1 => yes */
    int setwritemat;  /**< Flag, @see writemat */
    char writematstem[VMAX_ARGLEN];  /**< File stem to write mat */
    int writematflag;  /**< What matrix should we write:
                        * \li 0 => Poisson (differential operator)
                        * \li 1 => Poisson-Boltzmann operator linearized around
                        * solution (if applicable) */

    int parsed;  /**< Has this been filled with anything other
                  * than the default values? */
};

/**
 *  @ingroup PBEparm
 *  @brief   Declaration of the PBEparm class as the PBEparm structure
 */
typedef struct sPBEparm PBEparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (mcsh.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Get charge (e) of specified ion species
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns Charge of ion species (e)
 */
VEXTERNC double PBEparm_getIonCharge(
        PBEparm *thee, /**< PBEparm object */
        int iion  /**< Ion species ID/index */
        );

/** @brief   Get concentration (M) of specified ion species
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns Concentration of ion species (M)
 */
VEXTERNC double PBEparm_getIonConc(
        PBEparm *thee, /**< PBEparm object */
        int iion /**< Ion species ID/index */
        );

/** @brief   Get radius (A) of specified ion species
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns Radius of ion species (A)
 */
VEXTERNC double PBEparm_getIonRadius(
        PBEparm *thee, /**< PBEparm object */
        int iion /**< Ion species ID/index */
        );


/** @brief   Construct PBEparm object
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns Newly allocated and initialized PBEparm object
 */
VEXTERNC PBEparm* PBEparm_ctor();

/** @brief   FORTRAN stub to construct PBEparm object
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns 1 if succesful, 0 otherwise
 */
VEXTERNC int PBEparm_ctor2(
        PBEparm *thee /**< Memory location for object */
        );

/** @brief   Object destructor
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 */
VEXTERNC void PBEparm_dtor(
        PBEparm **thee /**< Pointer to memory location of object */
        );

/** @brief   FORTRAN stub for object destructor
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 */
VEXTERNC void PBEparm_dtor2(
        PBEparm *thee /**< Pointer to object to be destroyed */
        );

/** @brief   Consistency check for parameter values stored in object
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @returns 1 if OK, 0 otherwise
 */
VEXTERNC int PBEparm_check(
        PBEparm *thee /**< Object to be checked */
        );

/** @brief   Copy PBEparm object into thee
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 */
VEXTERNC void PBEparm_copy(
        PBEparm *thee, /**< Target for copy */
        PBEparm *parm /**< Source for copy */
        );

/** @brief   Parse a keyword from an input file
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @return   1 if matched and assigned; -1 if matched, but there's some sort
 *            of error (i.e., too few args); 0 if not matched
 */
VEXTERNC int PBEparm_parseToken(
        PBEparm *thee, /**< Parsing object */
        char tok[VMAX_BUFSIZE], /**< Token to parse */
        Vio *sock /**< Socket for additional tokens */
        );


#endif

