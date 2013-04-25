/** @defgroup FEMparm FEMparm class
 *  @brief    Parameter structure for FEM-specific variables from input files
 */

/**
 *  @file     femparm.h
 *  @ingroup  FEMparm
 *  @brief    Contains declarations for class FEMparm
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


#ifndef _FEMPARM_H_
#define _FEMPARM_H_

/* Generic header files */
#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vstring.h"

/**
 * @brief  Adaptive refinment error estimate tolerance key
 * @ingroup FEMparm
 * @author  Nathan Baker
 */
enum eFEMparm_EtolType {
    FET_SIMP=0,  /**< per-simplex error tolerance */
    FET_GLOB=1,  /**< global error tolerance */
    FET_FRAC=2   /**< fraction of simplices we want to have refined */
};

/**
 * @brief  Declare FEparm_EtolType type
 * @ingroup  FEMparm
 * @author  Nathan Baker
 */
typedef enum eFEMparm_EtolType FEMparm_EtolType;

/**
 * @brief  Adaptive refinment error estimator method
 * @ingroup FEMparm
 * @note  Do not change these values; they correspond to settings in FEtk
 * @author  Nathan Baker
 */
enum eFEMparm_EstType {
    FRT_UNIF=0,  /**< Uniform refinement */
    FRT_GEOM=1,  /**< Geometry-based (i.e. surfaces and charges) refinement */
    FRT_RESI=2,  /**< Nonlinear residual estimate-based refinement */
    FRT_DUAL=3,  /**< Dual-solution weight nonlinear residual estimate-based
                  * refinement */
    FRT_LOCA=4  /**< Local problem error estimate-based refinement */
};

/**
 * @brief  Declare FEMparm_EstType type
 * @ingroup  FEMparm
 */
typedef enum eFEMparm_EstType FEMparm_EstType;

/**
 * @brief  Calculation type
 * @ingroup FEMparm
 */
enum eFEMparm_CalcType {
    FCT_MANUAL,  /**< fe-manual */
    FCT_NONE  /**< unspecified */
};

/**
 * @brief  Declare FEMparm_CalcType type
 * @ingroup  FEMparm
 */
typedef enum eFEMparm_CalcType FEMparm_CalcType;

/**
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @brief   Parameter structure for FEM-specific variables from input files
 */
struct sFEMparm {

    int parsed;  /**< Flag:  Has this structure been filled with
                  * anything other than * the default values? (0 = no,
                  * 1 = yes) */
    FEMparm_CalcType type;  /**<  Calculation type */
    int settype;  /**< Boolean */
    double glen[3];  /**< Domain side lengths (in &Aring;) */
    int setglen;  /**< Boolean */
    double etol;  /**< Error tolerance for refinement; interpretation depends
                   * on the adaptive refinement method chosen */
    int setetol;  /**< Boolean */
    FEMparm_EtolType ekey;  /**< Adaptive refinment interpretation of error
                            * tolerance */
    int setekey;  /**< Boolean */
    FEMparm_EstType akeyPRE;  /**< Adaptive refinment error estimator method
                               * for pre-solution refine.  Note, this should
                               * either be FRT_UNIF or FRT_GEOM.  */
    int setakeyPRE;  /**< Boolean */
    FEMparm_EstType akeySOLVE;  /**< Adaptive refinment error estimator method
                               * for a posteriori solution-based refinement. */
    int setakeySOLVE;  /**< Boolean */
    int targetNum;    /**< Initial mesh will continue to be marked and refined
                        * with the method specified by akeyPRE until the mesh
                        * contains this many vertices or until targetRes is
                        * reached. */
    int settargetNum;  /**< Boolean */
    double targetRes; /**< Initial mesh will continue to be marked and refined
                        * with the method specified by akeyPRE until the mesh
                        * contains no markable simplices with longest edges
                        * above this size or until targetNum is reached. */
    int settargetRes;  /**< Boolean */
    int maxsolve;  /**< Maximum number of solve-estimate-refine cycles */
    int setmaxsolve;  /**< Boolean */
    int maxvert;  /**< Maximum number of vertices in mesh (ignored if less
                   * than zero) */
    int setmaxvert;  /**< Boolean */
    int pkey;		/**< Boolean sets the pkey type for going into AM_Refine
                      * pkey = 0 for non-HB based methods
                      * pkey = 1 for HB based methods */
    int useMesh;  /**< Indicates whether we use external finite element mesh */
    int meshID;  /**< External finite element mesh ID (if used) */

};

/** @typedef FEMparm
 *  @ingroup FEMparm
 *  @brief   Declaration of the FEMparm class as the FEMparm structure
 */
typedef struct sFEMparm FEMparm;

/* ///////////////////////////////////////////////////////////////////////////
// Class NOsh: Non-inlineable methods (nosh.c)
/////////////////////////////////////////////////////////////////////////// */

/** @brief   Construct FEMparm
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param  type  FEM calculation type
 *  @returns Newly allocated and initialized Vpmgp object
 */
VEXTERNC FEMparm* FEMparm_ctor(FEMparm_CalcType type);

/** @brief   FORTRAN stub to construct FEMparm
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param   thee Pointer to allocated FEMparm object
 *  @param  type  FEM calculation type
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int FEMparm_ctor2(FEMparm *thee, FEMparm_CalcType type);

/** @brief   Object destructor
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of FEMparm object
 */
VEXTERNC void FEMparm_dtor(FEMparm **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup FEMparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to FEMparm object
 */
VEXTERNC void FEMparm_dtor2(FEMparm *thee);

/**
 * @brief   Consistency check for parameter values stored in object
 * @ingroup FEMparm
 * @author  Nathan Baker
 * @param   thee   FEMparm object
 * @returns 1 if OK, 0 otherwise
 */
VEXTERNC int FEMparm_check(FEMparm *thee);

/**	@brief	Copy target object into thee
    @ingroup	FEMparm
    @author	Nathan Baker
    @param	thee	Destination object
    @param	source	Source object
*/
VEXTERNC void FEMparm_copy(FEMparm *thee, FEMparm *source);

/**
 * @brief   Parse an MG keyword from an input file
 * @ingroup MGparm
 * @author  Nathan Baker
 * @param   thee   MGparm object
 * @param   tok    Token to parse
 * @param   sock   Stream for more tokens
 * @return   VRC_SUCCESS if matched and assigned; VRC_FAILURE if matched, but there's
 * some sort of error (i.e., too few args); VRC_WARNING if not matched
 */
VEXTERNC Vrc_Codes FEMparm_parseToken(FEMparm *thee, char tok[VMAX_BUFSIZE],
  Vio *sock);

#endif

