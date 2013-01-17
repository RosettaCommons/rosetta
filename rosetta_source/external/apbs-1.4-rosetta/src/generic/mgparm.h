/** @defgroup MGparm MGparm class
 *  @brief    Parameter which holds useful parameters for generic multigrid
 *            calculations
 */

/**
 *  @file     mgparm.h
 *  @ingroup  MGparm
 *  @brief    Contains declarations for class MGparm
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


#ifndef _MGPARM_H_
#define _MGPARM_H_

/* Generic header files */
#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vstring.h"

/**
 * @brief  Calculation type
 * @ingroup MGparm
 */
enum eMGparm_CalcType {
    MCT_MANUAL=0,  /**< mg-manual */
    MCT_AUTO=1,  /**< mg-auto */
    MCT_PARALLEL=2,  /**< mg-para */
    MCT_DUMMY=3,  /**< mg-dummy */
    MCT_NONE=4  /**< unspecified */
};

/**
 * @brief  Declare MGparm_CalcType type
 * @ingroup  MGparm
 */
typedef enum eMGparm_CalcType MGparm_CalcType;

/**
 * @brief  Centering method
 * @ingroup MGparm
 */
enum eMGparm_CentMeth {
    MCM_POINT=0, /**< Center on a point */
    MCM_MOLECULE=1,  /**< Center on a molecule */
    MCM_FOCUS=2  /**< Determined by focusing */
};

/**
 * @brief  Declare MGparm_CentMeth type
 * @ingroup  MGparm
 */
typedef enum eMGparm_CentMeth MGparm_CentMeth;
/**
 *  @ingroup MGparm
 *  @author  Nathan Baker and Todd Dolinsky
 *  @brief   Parameter structure for MG-specific variables from input files
 *  @note    If you add/delete/change something in this class, the member
 *           functions -- especially MGparm_copy -- must be modified
 *           accordingly
 */
struct sMGparm {

    MGparm_CalcType type;  /**< What type of MG calculation? */
    int parsed;  /**< Has this structure been filled? (0 = no, 1 = yes) */

    /* *** GENERIC PARAMETERS *** */
    int dime[3];  /**< Grid dimensions */
    int setdime;  /**< Flag, @see dime */
    Vchrg_Meth chgm;  /**< Charge discretization method */
    int setchgm;  /**< Flag, @see chgm */
    Vchrg_Src  chgs; /**< Charge source (Charge, Multipole, Induced Dipole,
                      * NL Induced */

    /* *** TYPE 0 PARAMETERS (SEQUENTIAL MANUAL) *** */
    int nlev;  /**< Levels in multigrid hierarchy
                *   @deprecated Just ignored now */
    int setnlev;  /**< Flag, @see nlev */
    double etol;  /**< User-defined error tolerance */
    int setetol;  /**< Flag, @see etol */
    double grid[3];  /**< Grid spacings */
    int setgrid;  /**< Flag, @see grid */
    double glen[3];  /**< Grid side lengths. */
    int setglen;  /**< Flag, @see glen */
    MGparm_CentMeth cmeth;  /**< Centering method */
    double center[3];  /**< Grid center. If ispart = 0, then this is
                        * only meaningful if cmeth = 0.  However, if
                        * ispart = 1 and cmeth = MCM_PNT, then this is the
                        * center of the non-disjoint (overlapping)
                        * partition.  If ispart = 1 and cmeth = MCM_MOL, then
                        * this is the vector that must be added to the
                        * center of the molecule to give the center of
                        * the non-disjoint partition.  */
    int centmol;  /**< Particular molecule on which we want to center the grid.
        This should be the appropriate index in an array of molecules, not the
        positive definite integer specified by the user. */
    int setgcent;  /**< Flag, @see cmeth */

    /* ******** TYPE 1 & 2 PARAMETERS (SEQUENTIAL & PARALLEL AUTO-FOCUS) *** */
    double cglen[3];  /**< Coarse grid side lengths */
    int setcglen;  /**< Flag, @see cglen */
    double fglen[3];  /**< Fine grid side lengths */
    int setfglen;  /**< Flag, @see fglen */
    MGparm_CentMeth ccmeth;  /**< Coarse grid centering method */
    double ccenter[3];  /**< Coarse grid center.  */
    int ccentmol;  /**< Particular molecule on which we want to center the grid.
        This should be the appropriate index in an array of molecules, not the
        positive definite integer specified by the user. */
    int setcgcent;  /**< Flag, @see ccmeth */
    MGparm_CentMeth fcmeth;  /**< Fine grid centering method */
    double fcenter[3];  /**< Fine grid center.  */
    int fcentmol; /**< Particular molecule on which we want to center the grid.
        This should be the appropriate index in an array of molecules, not the
        positive definite integer specified by the user. */
    int setfgcent;  /**< Flag, @see fcmeth */


    /* ********* TYPE 2 PARAMETERS (PARALLEL AUTO-FOCUS) ******** */
    double partDisjCenter[3];  /**< This gives the center
                                     of the disjoint partitions */
    double partDisjLength[3];  /**< This gives the lengths of the disjoint
                                * partitions */
    int partDisjOwnSide[6];  /**< Tells whether the boundary points are ours
                              * (1) or not (0) */

    int pdime[3];  /**< Grid of processors to be used in calculation */
    int setpdime;  /**< Flag, @see pdime */
    int proc_rank;  /**< Rank of this processor */
    int setrank;  /**< Flag, @see proc_rank */
    int proc_size;  /**< Total number of processors */
    int setsize;  /**< Flag, @see proc_size */
    double ofrac;  /**< Overlap fraction between procs */
    int setofrac;  /**< Flag, @see ofrac */
    int async; /**< Processor ID for asynchronous calculation */
    int setasync; /**< Flag, @see asynch */

    int nonlintype; /**< Linearity Type Method to be used */
    int setnonlintype; /**< Flag, @see nonlintype */

    int method;		/**< Solver Method */
    int setmethod; /**< Flag, @see method */

    int useAqua;  /**< Enable use of lpbe/aqua */
    int setUseAqua; /**< Flag, @see useAqua */
};

/** @typedef MGparm
 *  @ingroup MGparm
 *  @brief   Declaration of the MGparm class as the MGparm structure
 */
typedef struct sMGparm MGparm;

/** @brief   Get number of grid points in x direction
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Number of grid points in the x direction
 */
VEXTERNC int MGparm_getNx(MGparm *thee);

/** @brief   Get number of grid points in y direction
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Number of grid points in the y direction
 */
VEXTERNC int MGparm_getNy(MGparm *thee);

/** @brief   Get number of grid points in z direction
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Number of grid points in the z direction
 */
VEXTERNC int MGparm_getNz(MGparm *thee);

/** @brief   Get grid spacing in x direction (&Aring;)
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Grid spacing in the x direction
 */
VEXTERNC double MGparm_getHx(MGparm *thee);

/** @brief   Get grid spacing in y direction (&Aring;)
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Grid spacing in the y direction
 */
VEXTERNC double MGparm_getHy(MGparm *thee);

/** @brief   Get grid spacing in z direction (&Aring;)
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns Grid spacing in the z direction
 */
VEXTERNC double MGparm_getHz(MGparm *thee);

/** @brief   Set center x-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @param   x     x-coordinate
 */
VEXTERNC void MGparm_setCenterX(MGparm *thee, double x);

/** @brief   Set center y-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @param   y     y-coordinate
 */
VEXTERNC void MGparm_setCenterY(MGparm *thee, double y);

/** @brief   Set center z-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @param   z     z-coordinate
 */
VEXTERNC void MGparm_setCenterZ(MGparm *thee, double z);

/** @brief   Get center x-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  x-coordinate
 */
VEXTERNC double MGparm_getCenterX(MGparm *thee);

/** @brief   Get center y-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  y-coordinate
 */
VEXTERNC double MGparm_getCenterY(MGparm *thee);

/** @brief   Get center z-coordinate
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  MGparm object
 *  @returns  z-coordinate
 */
VEXTERNC double MGparm_getCenterZ(MGparm *thee);

/** @brief   Construct MGparm object
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   type Type of MG calculation
 *  @returns Newly allocated and initialized MGparm object
 */
VEXTERNC MGparm*  MGparm_ctor(MGparm_CalcType type);

/** @brief   FORTRAN stub to construct MGparm object
 *  @ingroup MGparm
 *  @author  Nathan Baker and Todd Dolinsky
 *  @param   thee Space for MGparm object
 *  @param   type Type of MG calculation
 *  @returns Success enumeration
 */
VEXTERNC Vrc_Codes      MGparm_ctor2(MGparm *thee, MGparm_CalcType type);

/** @brief   Object destructor
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to memory location of MGparm object
 */
VEXTERNC void     MGparm_dtor(MGparm **thee);

/** @brief   FORTRAN stub for object destructor
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee  Pointer to MGparm object
 */
VEXTERNC void     MGparm_dtor2(MGparm *thee);

/** @brief   Consistency check for parameter values stored in object
 *  @ingroup MGparm
 *  @author  Nathan Baker
 *  @param   thee   MGparm object
 *  @returns Success enumeration
 */
VEXTERNC Vrc_Codes      MGparm_check(MGparm *thee);

/** @brief   Copy MGparm object into thee
 *  @ingroup MGparm
 *  @author  Nathan Baker and Todd Dolinsky
 *  @param   thee   MGparm object (target for copy)
 *  @param   parm   MGparm object (source for copy)
 */
VEXTERNC void     MGparm_copy(MGparm *thee, MGparm *parm);

/** @brief   Parse an MG keyword from an input file
 *  @ingroup MGparm
 *  @author  Nathan Baker and Todd Dolinsky
 *  @param   thee   MGparm object
 *  @param   tok    Token to parse
 *  @param   sock   Stream for more tokens
 *  @returns Success enumeration (1 if matched and assigned; -1 if matched, but there's some sort
 *            of error (i.e., too few args); 0 if not matched)
 */
VEXTERNC Vrc_Codes      MGparm_parseToken(MGparm *thee, char tok[VMAX_BUFSIZE],
                    Vio *sock);

#endif

