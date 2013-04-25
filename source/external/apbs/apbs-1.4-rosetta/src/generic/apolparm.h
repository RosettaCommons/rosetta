/** @defgroup APOLparm APOLparm class
 *  @brief    Parameter structure for APOL-specific variables from input files
 */

/**
 *  @file     femparm.h
 *  @ingroup  APOLparm
 *  @brief    Contains declarations for class APOLparm
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


#ifndef _APOLPARM_H_
#define _APOLPARM_H_

/* Generic header files */
#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vstring.h"
#include "generic/vparam.h"

/**
* @ingroup APOLparm
 * @brief  Define energy calculation enumeration
 */
enum eAPOLparm_calcEnergy {
    ACE_NO=0, /**< Do not perform energy calculation */
    ACE_TOTAL=1, /**< Calculate total energy only */
    ACE_COMPS=2 /**< Calculate per-atom energy components */
};

/**
* @ingroup APOLparm
 * @brief  Define eAPOLparm_calcEnergy enumeration as APOLparm_calcEnergy
 */
typedef enum eAPOLparm_calcEnergy APOLparm_calcEnergy;

/**
* @ingroup APOLparm
 * @brief  Define force calculation enumeration
 */
enum eAPOLparm_calcForce {
    ACF_NO=0, /**< Do not perform force calculation */
    ACF_TOTAL=1, /**< Calculate total force only */
    ACF_COMPS=2 /**< Calculate per-atom force components */
};

/**
* @ingroup APOLparm
 * @brief  Define eAPOLparm_calcForce enumeration as APOLparm_calcForce
 */
typedef enum eAPOLparm_calcForce APOLparm_calcForce;

/**
* @ingroup APOLparm
 * @brief  Define force calculation enumeration
 */
enum eAPOLparm_doCalc {
    ACD_NO=0, /**< Do not perform calculation */
    ACD_YES=1, /**< Perform calculations */
    ACD_ERROR=2 /**< Error setting up calculation */
};

/**
* @ingroup APOLparm
 * @brief  Define eAPOLparm_calcForce enumeration as APOLparm_calcForce
 */
typedef enum eAPOLparm_doCalc APOLparm_doCalc;


/**
 *  @ingroup APOLparm
 *  @author  David Gohara
 *  @brief   Parameter structure for APOL-specific variables from input files
 */
struct sAPOLparm {

    int parsed;  /**< Flag:  Has this structure been filled with anything other than the default values? (0 = no, 1 = yes) */

    double grid[3];  /**< Grid spacing */
    int setgrid;  /**< Flag, @see grid */

    int molid;  /**< Molecule ID to perform calculation on */
    int setmolid;  /**< Flag, @see molid */

    double bconc; /**< Vacc sphere density */
    int setbconc; /**< Flag, @see bconc */

    double sdens; /**< Vacc sphere density */
    int setsdens; /**< Flag, @see sdens */

    double dpos; /**< Atom position offset */
    int setdpos; /**< Flag, @see dpos */

    double press; /**< Solvent pressure */
    int setpress; /**< Flag, @see press */

    Vsurf_Meth srfm;  /**< Surface calculation method */
    int setsrfm;  /**< Flag, @see srfm */

    double srad;  /**< Solvent radius */
    int setsrad;  /**< Flag, @see srad */

    double swin;  /**< Cubic spline window */
    int setswin;  /**< Flag, @see swin */

    double temp;  /**< Temperature (in K) */
    int settemp;  /**< Flag, @see temp */

    double gamma;  /**< Surface tension for apolar energies/forces
                    * (in kJ/mol/A^2) */
    int setgamma;  /**< Flag, @see gamma */

    APOLparm_calcEnergy calcenergy;  /**< Energy calculation flag */
    int setcalcenergy;  /**< Flag, @see calcenergy */

    APOLparm_calcForce calcforce;  /**< Atomic forces calculation */
    int setcalcforce;  /**< Flag, @see calcforce */

    double watsigma;  /**< Water oxygen Lennard-Jones radius (A) */
    double watepsilon;  /**< Water oxygen Lennard-Jones well depth (kJ/mol) */
    double sasa; /**< Solvent accessible surface area for this calculation */
    double sav;   /**< Solvent accessible volume for this calculation */
    double wcaEnergy; /**< wcaEnergy */
    double totForce[3]; /**< Total forces on x, y, z */

    int setwat; /**< Boolean for determining if a water parameter
                 *  is supplied. Yes = 1, No = 0 */
};

/** @typedef APOLparm
 *  @ingroup APOLparm
 *  @brief   Declaration of the APOLparm class as the APOLparm structure
 */
typedef struct sAPOLparm APOLparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (nosh.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct APOLparm
 *  @ingroup APOLparm
 *  @author  David Gohara
 *  @returns Newly allocated and initialized Vpmgp object
 */
VEXTERNC APOLparm* APOLparm_ctor();

/** @brief   FORTRAN stub to construct APOLparm
 *  @ingroup APOLparm
 *  @author  David Gohara, Yong Huang
 *  @param   thee Pointer to allocated APOLparm object
 *  @returns Success enumeration
 */
VEXTERNC Vrc_Codes APOLparm_ctor2(APOLparm *thee);

/** @brief   Object destructor
 *  @ingroup APOLparm
 *  @author  David Gohara
 *  @param   thee  Pointer to memory location of APOLparm object
 */
VEXTERNC void APOLparm_dtor(APOLparm **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup APOLparm
 *  @author  David Gohara
 *  @param   thee  Pointer to APOLparm object
 */
VEXTERNC void APOLparm_dtor2(APOLparm *thee);

/**
 * @brief   Consistency check for parameter values stored in object
 * @ingroup APOLparm
 * @author  David Gohara, Yong Huang
 * @param   thee   APOLparm object
 * @returns Success enumeration
 */
VEXTERNC Vrc_Codes APOLparm_check(APOLparm *thee);

/**	@brief	Copy target object into thee
    @ingroup	APOLparm
    @author	Nathan Baker
    @param	thee	Destination object
    @param	source	Source object
*/
VEXTERNC void APOLparm_copy(APOLparm *thee, APOLparm *source);

/**
 * @brief   Parse an MG keyword from an input file
 * @ingroup MGparm
 * @author  David Gohara
 * @param   thee   MGparm object
 * @param   tok    Token to parse
 * @param   sock   Stream for more tokens
 * @returns Success enumeration (1 if matched and assigned; -1 if matched, but there's
 * some sort of error (i.e., too few args); 0 if not matched)
 */
VEXTERNC Vrc_Codes APOLparm_parseToken(APOLparm *thee, char tok[VMAX_BUFSIZE],
  Vio *sock);

#endif

