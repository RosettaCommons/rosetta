/** @defgroup Vhal Vhal class
 *  @brief    A "class" which consists solely of macro definitions which are
 *            used by several other classes
 */

/**
 *  @file       vhal.h
 *  @ingroup    Vhal
 *  @brief      Contains generic macro definitions for APBS
 *  @version  $Id$
 *  @author   Nathan A. Baker
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (nathan.baker@pnl.gov)
 * Pacific Northwest National Laboratory
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2012 Battelle Memorial Institute. Developed at the Pacific Northwest National Laboratory, operated by Battelle Memorial Institute, Pacific Northwest Division for the U.S. Department Energy.  Portions Copyright (c) 2002-2010, Washington University in St. Louis.  Portions Copyright (c) 2002-2010, Nathan A. Baker.  Portions Copyright (c) 1999-2002, The Regents of the University of California. Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#ifndef _VAPBSHAL_H_
#define _VAPBSHAL_H_

#include "stdio.h"

FILE *data;

/**
 *    @ingroup Vhal
 *    @author     David Gohara
 *    @brief     Return code enumerations
 *  @note  Note that the enumerated values are opposite the standard for FAILURE and SUCCESS
 */
enum eVrc_Codes {

    VRC_WARNING=-1, /** A non-fatal error  */
    VRC_FAILURE=0,  /** A fatal error */
    VRC_SUCCESS=1   /** A successful execution */

};
typedef enum eVrc_Codes Vrc_Codes;

/**
 *    @ingroup Vhal
 *    @author     David Gohara
 *    @brief     Solution Method enumerations
 *  @note  Note that the enumerated values are opposite the standard for FAILURE and SUCCESS
 */
enum eVsol_Meth {

    VSOL_CGMG,        /* 0: conjugate gradient multigrid */
    VSOL_Newton,    /* 1: newton */
    VSOL_MG,        /* 2: multigrid */
    VSOL_CG,        /* 3: conjugate gradient */
    VSOL_SOR,        /* 4: sucessive overrelaxation */
    VSOL_RBGS,        /* 5: red-black gauss-seidel */
    VSOL_WJ,        /* 6: weighted jacobi */
    VSOL_Richardson,/* 7: richardson  */
    VSOL_CGMGAqua,    /* 8: conjugate gradient multigrid aqua */
    VSOL_NewtonAqua    /* 9: newton aqua */

};
typedef enum eVsol_Meth Vsol_Meth;

/**
 *  @ingroup Vhal
 *  @author  Nathan Baker
 *  @brief   Types of molecular surface definitions
 */
enum eVsurf_Meth {
    VSM_MOL=0, /**<  Ion accessibility is defined using inflated van der Waals
                *    radii, the dielectric coefficient ( ) is defined using the
                *    molecular (Conolly) surface definition without
                *    smoothing */
    VSM_MOLSMOOTH=1, /**<  As VSM_MOL but with a simple harmonic average
                      *    smoothing */
    VSM_SPLINE=2,    /**<  Spline-based surface definitions. This is primarily
                      * for use with force calculations, since it requires
                      * substantial reparameterization of radii. This is based
                      * on the work of Im et al, Comp. Phys.  Comm. 111 ,
                      * (1998) and uses a cubic spline to define a smoothly
                      * varying characteristic function for the surface-based
                      * parameters. Ion accessibility is defined using inflated
                      * van der Waals radii with the spline function and the
                      * dielectric coefficient is defined using the standard
                      * van der Waals radii with the spline function.  */
    VSM_SPLINE3=3,  /**<  A 5th order polynomial spline is used to create
                      *  a smoothly varying characteristic function
                      *  (continuity through 2nd derivatives) for surface
                      *  based paramters. */
    VSM_SPLINE4=4   /**<  A 7th order polynomial spline is used to create
                      *  a smoothly varying characteristic function
                      *  (continuity through 3rd derivatives) for surface
                      *  based paramters. */
};

/** @typedef Vsurf_Meth
 *  @ingroup Vhal
 *  @brief   Declaration of the Vsurf_Meth type as the Vsurf_Meth enum
 */
typedef enum eVsurf_Meth Vsurf_Meth;

/**
 * @brief  Version of PBE to solve
 * @ingroup  Vhal
 */
enum eVhal_PBEType {
    PBE_LPBE,  /**<  Traditional Poisson-Boltzmann equation, linearized */
    PBE_NPBE,  /**<  Traditional Poisson-Boltzmann equation, full */
    PBE_LRPBE,  /**<  Regularized Poisson-Boltzmann equation, linearized */
    PBE_NRPBE,  /** <  Regularized Poisson-Boltzmann equation, full */
    PBE_SMPBE    /**< SM PBE */
};

/**
 *  @brief   Declaration of the Vhal_PBEType type as the Vhal_PBEType enum
 *  @ingroup Vhal
 */
typedef enum eVhal_PBEType Vhal_PBEType;

/**
* @brief  Type of ipkey to use for MG methods
 * @ingroup  Vhal
 */
enum eVhal_IPKEYType {
    IPKEY_SMPBE = -2, /**<  SMPBE ipkey */
    IPKEY_LPBE,  /**<  LPBE ipkey */
    IPKEY_NPBE /**<  NPBE ipkey */
};

/**
*  @brief   Declaration of the Vhal_IPKEYType type as the Vhal_IPKEYType enum
*  @ingroup Vhal
*/
typedef enum eVhal_IPKEYType Vhal_IPKEYType;

/**
* @brief  Type of nonlinear to use for MG methods
 * @ingroup  Vhal
 */
enum eVhal_NONLINType {
    NONLIN_LPBE = 0,
    NONLIN_NPBE,
    NONLIN_SMPBE,
    NONLIN_LPBEAQUA,
    NONLIN_NPBEAQUA
};

/**
*  @brief   Declaration of the Vhal_NONLINType type as the Vhal_NONLINType enum
*  @ingroup Vhal
*/
typedef enum eVhal_NONLINType Vhal_NONLINType;

/**
 * @brief Output file format
 * @ingroup Vhal
 */
enum eVoutput_Format {
    OUTPUT_NULL,   /**< No output */
    OUTPUT_FLAT, /**< Output in flat-file format */
};

/**
 * @brief Declaration of the Voutput_Format type as the VOutput_Format enum
 * @ingroup Vhal
 */
typedef enum eVoutput_Format Voutput_Format;

/**
 * @ingroup  Vhal
 * @author  Nathan Baker
 * @brief  Types of boundary conditions
 */
enum eVbcfl {
    BCFL_ZERO=0,  /**< Zero Dirichlet boundary conditions */
    BCFL_SDH=1,  /**< Single-sphere Debye-Huckel Dirichlet boundary
                  * condition */
    BCFL_MDH=2,  /**< Multiple-sphere Debye-Huckel Dirichlet boundary
                  * condition */
    BCFL_UNUSED=3,  /**< Unused boundary condition method (placeholder) */
    BCFL_FOCUS=4,  /**< Focusing Dirichlet boundary condition */
    BCFL_MEM=5,  /**< Focusing membrane boundary condition */
    BCFL_MAP=6    /**< Skip first level of focusing use an external map */
};

/**
 * @brief  Declare Vbcfl type
 * @ingroup  Vhal
 */
typedef enum eVbcfl Vbcfl;

/**
 *  @ingroup Vhal
 *  @author  Nathan Baker
 *  @brief   Types of charge discretization methods
 */
enum eVchrg_Meth {
    VCM_TRIL=0,  /**< Trilinear interpolation of charge to 8 nearest grid
                  *   points.  The traditional method; not particularly good to
                  *   use with PBE forces. */
    VCM_BSPL2=1,  /**< Cubic B-spline across nearest- and
                  *   next-nearest-neighbors.  Mainly for use in grid-sensitive
                  *   applications (such as force calculations). */
    VCM_BSPL4=2  /**< 5th order B-spline for AMOEBA permanent multipoles. */
};

/** @typedef Vchrg_Meth
 *  @ingroup Vhal
 *  @brief   Declaration of the Vchrg_Meth type as the Vchrg_Meth enum
 */
typedef enum eVchrg_Meth Vchrg_Meth;

/**
 *  @ingroup Vhal
 *  @author  Michael Schnieders
 *  @brief   Charge source
 */
enum eVchrg_Src {
    VCM_CHARGE=0,     /**< Partial Charge source distribution */
    VCM_PERMANENT=1,  /**< Permanent Multipole source distribution */
    VCM_INDUCED=2,    /**< Induced Dipole source distribution */
    VCM_NLINDUCED=3   /**< NL Induced Dipole source distribution */
};

/** @typedef Vchrg_Src
 *  @ingroup Vhal
 *  @brief   Declaration of the Vchrg_Src type as the Vchrg_Meth enum
 */
typedef enum eVchrg_Src Vchrg_Src;

/**
 *  @ingroup Vhal
 *  @author  Nathan Baker
 *  @brief   Types of (scalar) data that can be written out of APBS
 */
enum eVdata_Type {
    VDT_CHARGE, /**< Charge distribution (e) */
    VDT_POT,    /**< Potential (kT/e) */
    VDT_ATOMPOT, /**< Atom potential (kT/e) */
    VDT_SMOL,   /**< Solvent accessibility defined by molecular/Connolly
                 * surface definition (1 = accessible, 0 = inaccessible) */
    VDT_SSPL,   /**< Spline-based solvent accessibility (1 = accessible, 0 =
                 * inaccessible) */
    VDT_VDW,    /**< van der Waals-based accessibility (1 = accessible, 0 =
                 * inaccessible) */
    VDT_IVDW,   /**< Ion accessibility/inflated van der Waals (1 =
                 * accessible, 0 = inaccessible) */
    VDT_LAP,    /**< Laplacian of potential (kT/e/A^2) */
    VDT_EDENS,  /**< Energy density \f$\epsilon (\nabla u)^2\f$, where \f$u\f$
                 * is potential (kT/e/A)^2 */
    VDT_NDENS,  /**< Ion number density \f$\sum c_i \exp (-q_i u)^2\f$,
             * where \f$u\f$ is potential (output in M) */
    VDT_QDENS,  /**< Ion charge density \f$\sum q_i c_i \exp (-q_i u)^2\f$,
             * where \f$u\f$ is potential (output in \f$e_c M\f$) */
    VDT_DIELX,  /**< Dielectric x-shifted map as calculated with the currently
                 * specified scheme (dimensionless) */
    VDT_DIELY,  /**< Dielectric y-shifted map as calculated with the currently
                 * specified scheme (dimensionless) */
    VDT_DIELZ,  /**< Dielectric y-shifted map as calculated with the currently
                 * specified scheme (dimensionless) */
    VDT_KAPPA   /**< Kappa map as calculated with the currently
                 * specified scheme (\f$\AA^{-3}\f$) */
};

/** @typedef Vdata_Type
 *  @ingroup Vhal
 *  @brief   Declaration of the Vdata_Type type as the Vdata_Type enum
 */
typedef enum eVdata_Type Vdata_Type;

/**
 *  @ingroup Vhal
 *  @author  Nathan Baker
 *  @brief   Format of data for APBS I/O
 */
enum eVdata_Format {
    VDF_DX=0,  /**< OpenDX (Data Explorer) format */
    VDF_UHBD=1, /**< UHBD format */
    VDF_AVS=2,  /**< AVS UCD format */
    VDF_MCSF=3,  /**< FEtk MC Simplex Format (MCSF) */
    VDF_GZ=4,    /**< Binary file (GZip) */
    VDF_FLAT=5  /**< Write flat file */
};

/** @typedef Vdata_Format
 *  @ingroup Vhal
 *  @brief   Declaration of the Vdata_Format type as the Vdata_Format enum
 */
typedef enum eVdata_Format Vdata_Format;

/**
 * @brief  APBS total execution timer ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_WALL_CLOCK 26

/**
 * @brief  APBS setup timer ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_SETUP 27

/**
 * @brief  APBS solver timer ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_SOLVER 28

/**
 * @brief  APBS energy timer ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_ENERGY 29

/**
 * @brief  APBS force timer ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_FORCE 30

/**
 * @brief  APBS temp timer #1 ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_TEMP1 31

/**
 * @brief  APBS temp timer #2 ID
 * @ingroup  Vhal
 */
#define APBS_TIMER_TEMP2 32

/** @brief The maximum number of molecules that can be involved in a single
 *         PBE calculation
 *  @ingroup Vhal
 */
#define MAXMOL 5

/** @brief The maximum number of ion species that can be involved in a single
 *         PBE calculation
 *  @ingroup Vhal
 */
#define MAXION 10

/** @brief The maximum number of times an MG calculation can be focused
 *  @ingroup Vhal
 */
#define MAXFOCUS 5

/** @brief   Minimum number of levels in a multigrid calculations
 *  @ingroup Vhal
 */
#define VMGNLEV 4

/** @brief   Maximum reduction of grid spacing during a focusing calculation
 *  @ingroup Vhal
 */
#define VREDFRAC 0.25

/** @brief  Number of vertices per simplex (hard-coded to 3D)
 *  @ingroup Vhal
 */
#define VAPBS_NVS 4

/** @brief  Our dimension
 * @ingroup Vhal
 */
#define VAPBS_DIM 3

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH
 *  @ingroup Vhal
 */
#define VAPBS_RIGHT 0

/** @brief    Maximum number of points on a sphere
    @note    Used by VaccSurf
    @ingroup    Vhal
    */
#define MAX_SPHERE_PTS 50000

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH
 *  @ingroup Vhal
 */
#define VAPBS_FRONT 1

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH
 *  @ingroup Vhal
 */
#define VAPBS_UP    2

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH
 *  @ingroup Vhal
 */
#define VAPBS_LEFT  3

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH
 *  @ingroup Vhal
 */
#define VAPBS_BACK  4

/** @brief   Face definition for a volume
 *  @note    Consistent with PMG if RIGHT = EAST, BACK = SOUTH
 *  @ingroup Vhal
 */
#define VAPBS_DOWN  5

/** @brief   A small number used in Vpmg to decide if points are on/off
 *           grid-lines or non-zer0 (etc.)
 *  @ingroup Vhal
 */
#define VPMGSMALL 1e-12

/** @brief   Used to set the min values acceptable for sinh chopping
 *  @def SINH_MIN
 *  @ingroup Vhal
 */
#define SINH_MIN -85.0

/** @brief   Used to set the max values acceptable for sinh chopping
 *  @def SINH_MAX
 *  @ingroup Vhal
 */
#define SINH_MAX 85.0

#define MAX_HASH_DIM 75

#if defined(VDEBUG)
#   if !defined(APBS_NOINLINE)
#       define APBS_NOINLINE 1
#   endif
#endif

#if !defined(APBS_NOINLINE)

/** @brief   Turns on inlining macros in Vacc class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VACC

/** @brief   Turns on inlining macros in Vatom class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VATOM

/** @brief   Turns on inlining macros in Vcsm class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VCSM

/** @brief   Turns on inlining macros in Vpbe class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VPBE

/** @brief   Turns on inlining macros in Vpee class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VPEE

/** @brief   Turns on inlining macros in Vgreen class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VGREEN

/** @brief   Turns on inlining macros in Vfetk class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VFETK

/** @brief   Turns on inlining macros in Vpmg class if defined
 *  @ingroup Vhal
 */
#   define VINLINE_VPMG

/** @brief   Defines the maximum hash table size in any direction
 *  @ingroup Vhal
 */

#endif

/* Fortran name mangling */
#if defined(VF77_UPPERCASE)
#   if defined(VF77_NOUNDERSCORE)
#       define VF77_MANGLE(name,NAME) NAME
#   elif defined(VF77_ONEUNDERSCORE)
#       define VF77_MANGLE(name,NAME) NAME ## _
#   else
#       define VF77_MANGLE(name,NAME) name
#   endif
#else
#   if defined(VF77_NOUNDERSCORE)
#       define VF77_MANGLE(name,NAME) name
#   elif defined(VF77_ONEUNDERSCORE)
#       define VF77_MANGLE(name,NAME) name ## _
#   else
        /** @brief  Name-mangling macro for using FORTRAN functions in C code
         *  @ingroup  Vhal
         */
#       define VF77_MANGLE(name,NAME) name
#   endif
#endif

/* Floating Point Error */
#if defined(FLOAT_EPSILON)
#     define VFLOOR(value) \
              ((floor(value) != floor(value + FLOAT_EPSILON)) ? \
              floor(value + FLOAT_EPSILON) : floor(value))
#else
      /** @brief  Wrapped floor to fix floating point issues in the Intel
       * compiler
       *  @author Todd Dolinksy
       *  @ingroup Vhal
       */
#     define VFLOOR(value) floor(value)
#endif

/* String embedding for ident */
#if defined(HAVE_EMBED)
/**
 * @brief  Allows embedding of RCS ID tags in object files.
 * @author Mike Holst
 * @ingroup Vhal */
#    define VEMBED(rctag) \
         VPRIVATE const char* rctag; \
         static void* use_rcsid=(0 ? &use_rcsid : (void**)&rcsid);
#else
/**
 * @brief  Allows embedding of RCS ID tags in object files.
 * @author Mike Holst
 * @ingroup Vhal */
#    define VEMBED(rctag)
#endif /* if defined(HAVE_EMBED) */



/** OS specific flags and etcetera */
#if !defined(WIN32) || defined(__MINGW32__)
#define PRINT_FUNC __PRETTY_FUNCTION__
#define OS_SEP_STR "/"
#define OS_SEP_CHAR '/'
#else
#define OS_SEP_STR "\\"
#define OS_SEP_CHAR '\\'
#define PRINT_FUNC __FUNCSIG__
#define snprintf sprintf_s
#endif

#ifdef VERGBOSE_DEBUG
#define ANNOUNCE_FUNCTION                                \
        do {                                             \
            Vnm_prrint(2, "%s() [%s:%d]\n",              \
                       PRINT_FUNC, __FILE__, __LINE__ ); \
        } while(0)

#define WARN_UNTESTED                                             \
        do {                                                      \
            Vnm_print(2, "%s() [%s:%d]: Untested Translation!\n", \
                      __FUNCTION__, __FILE__, __LINE__);          \
        } while(0)

#define WARN_PARTTESTED                                                   \
        do{                                                               \
            Vnm_print(2, "%s() [%s:%d]: Partially Tested Translation.\n", \
                      __FUNCTION__, __FILE__, __LINE__);                  \
        } while(0)
#else
#define ANNOUNCE_FUNCTION
#define WARN_UNTESTED
#define WARN_PARTTESTED
#endif



/* Utility messages.  Print out messages with location information */
#ifdef DEBUG
#define VCHANNELEDMESSAGE0(channel, msg)                      \
        do {                                                  \
            Vnm_print(channel, "%s:%d [%s()]: MESSAGE:\n"     \
                         "    %s\n\n",                        \
                      __FILE__, __LINE__, __FUNCTION__, msg); \
        } while(0)

#define VCHANNELEDMESSAGE1(channel, msg, arg)                  \
        do {                                                   \
            char buff[1000];                                   \
            snprintf( buff, 1000, msg, arg );                  \
            Vnm_print(channel, "%s:%d [%s()]: MESSAGE:\n"      \
                         "    %s\n\n",                         \
                      __FILE__, __LINE__, __FUNCTION__, buff); \
        } while(0)

#define VCHANNELEDMESSAGE2(channel, msg, arg0, arg1)           \
        do {                                                   \
            char buff[1000];                                   \
            snprintf( buff, 1000, msg, arg0, arg1 );           \
            Vnm_print(channel, "%s:%d [%s()]: MESSAGE:\n"      \
                         "    %s\n\n",                         \
                      __FILE__, __LINE__, __FUNCTION__, buff); \
        } while(0)

#define VCHANNELEDMESSAGE3(channel, msg, arg0, arg1, arg2)     \
        do {                                                   \
            char buff[1000];                                   \
            snprintf(buff, 1000, msg, arg0, arg1, arg2);       \
            Vnm_print(channel, "%s:%d [%s()]: MESSAGE:\n"      \
                         "    %s\n\n",                         \
                      __FILE__, __LINE__, __FUNCTION__, buff); \
        } while(0)

#define VMESSAGE0(msg)                   VCHANNELEDMESSAGE0(2, msg)
#define VMESSAGE1(msg, arg0)             VCHANNELEDMESSAGE1(2, msg, arg0)
#define VMESSAGE2(msg, arg0, arg1)       VCHANNELEDMESSAGE2(2, msg, arg0, arg1)
#define VMESSAGE3(msg, arg0, arg1, arg2) VCHANNELEDMESSAGE3(2, msg, arg0, arg1, arg2)

#define VERRMSG0(msg)                   VMESSAGE0(msg)
#define VERRMSG1(msg, arg0)             VMESSAGE1(msg, arg0)
#define VERRMSG2(msg, arg0, arg1)       VMESSAGE2(msg, arg0, arg1)
#define VERRMSG3(msg, arg0, arg1, arg2) VMESSAGE3(msg, arg0, arg1, arg2)
#else
#define VCHANNELEDMESSAGE0(channel, msg)                       \
        do {                                                   \
            Vnm_print(channel, "%s: %s\n", __FUNCTION__, msg); \
        } while(0)

#define VCHANNELEDMESSAGE1(channel, msg, arg0)                  \
        do {                                                    \
            char buff[1000];                                    \
            snprintf( buff, 1000, msg, arg0 );                  \
            Vnm_print(channel, "%s: %s\n", __FUNCTION__, buff); \
        } while(0)

#define VCHANNELEDMESSAGE2(channel, msg, arg0, arg1)            \
        do {                                                    \
            char buff[1000];                                    \
            snprintf( buff, 1000, msg, arg0, arg1 );            \
            Vnm_print(channel, "%s: %s\n", __FUNCTION__, buff); \
        } while(0)

#define VCHANNELEDMESSAGE3(channel, msg, arg0, arg1, arg2)       \
        do {                                                     \
            char buff[1000];                                     \
            snprintf(buff, 1000, msg, arg0, arg1, arg2);         \
            Vnm_print(channel, "%s: %s\n", __FUNCTION__, buff);  \
        } while(0)

#define VMESSAGE0(msg)                   VCHANNELEDMESSAGE0(0, msg)
#define VMESSAGE1(msg, arg0)             VCHANNELEDMESSAGE1(0, msg, arg0)
#define VMESSAGE2(msg, arg0, arg1)       VCHANNELEDMESSAGE2(0, msg, arg0, arg1)
#define VMESSAGE3(msg, arg0, arg1, arg2) VCHANNELEDMESSAGE3(0, msg, arg0, arg1, arg2)

#define VERRMSG0(msg)                    VCHANNELEDMESSAGE0(2, msg)
#define VERRMSG1(msg, arg0)              VCHANNELEDMESSAGE1(2, msg, arg0)
#define VERRMSG2(msg, arg0, arg1)        VCHANNELEDMESSAGE2(2, msg, arg0, arg1)
#define VERRMSG3(msg, arg0, arg1, arg2)  VCHANNELEDMESSAGE3(2, msg, arg0, arg1, arg2)
#endif



/* Utility assertions.  If they fail, they print out messages with possible
 * arguments and then abort
 * The do{...} while(0) simply enforces that a semicolon appears at the end
 */
#ifdef DEBUG
#define VASSERT_MSG0(cnd, msg)                                                                           \
    do {                                                            \
        if( (cnd) == 0 ) {                                          \
            Vnm_print(2, "%s:%d [%s()]: ERROR:\n"                   \
                         "    Assertion Failed (%s): %s\n\n",       \
                      __FILE__, __LINE__, __FUNCTION__, #cnd, msg); \
            abort();                                                \
        }                                                           \
    } while(0)

#define VASSERT_MSG1(cnd, msg, arg)                                  \
    do {                                                             \
        if( (cnd) == 0 ) {                                           \
            char buff[1000];                                         \
            snprintf( buff, 1000, msg, arg );                        \
            Vnm_print(2, "%s:%d [%s()]: ERROR:\n"                    \
                         "    Assertion Failed (%s): %s\n\n",        \
                      __FILE__, __LINE__, __FUNCTION__, #cnd, buff); \
            abort();                                                 \
        }                                                            \
    } while(0)

#define VASSERT_MSG2(cnd, msg, arg0, arg1)                           \
    do {                                                             \
        if( (cnd) == 0 ) {                                           \
            char buff[1000];                                         \
            snprintf( buff, 1000, msg, arg0, arg1 );                 \
            Vnm_print(2, "%s:%d [%s()]: ERROR:\n"                    \
                         "    Assertion Failed (%s): %s\n\n",        \
                      __FILE__, __LINE__, __FUNCTION__, #cnd, buff); \
            abort();                                                 \
        }                                                            \
    } while(0)
#else
#define VASSERT_MSG0(cnd, msg)                                                                           \
    do {                                                            \
        if( (cnd) == 0 ) {                                          \
            Vnm_print(2, "%[%s()]: ERROR:\n"                        \
                         "    Assertion Failed (%s): %s\n\n",       \
                      __FUNCTION__, #cnd, msg);                     \
            abort();                                                \
        }                                                           \
    } while(0)

#define VASSERT_MSG1(cnd, msg, arg)                                  \
    do {                                                             \
        if( (cnd) == 0 ) {                                           \
            char buff[1000];                                         \
            snprintf( buff, 1000, msg, arg );                        \
            Vnm_print(2, "[%s()]: ERROR:\n"                          \
                         "    Assertion Failed (%s): %s\n\n",        \
                      __FUNCTION__, #cnd, buff);                     \
            abort();                                                 \
        }                                                            \
    } while(0)

#define VASSERT_MSG2(cnd, msg, arg0, arg1)                           \
    do {                                                             \
        if( (cnd) == 0 ) {                                           \
            char buff[1000];                                         \
            snprintf( buff, 1000, msg, arg0, arg1 );                 \
            Vnm_print(2, "[%s()]: ERROR:\n"                          \
                         "    Assertion Failed (%s): %s\n\n",        \
                      __FUNCTION__, #cnd, buff);                     \
            abort();                                                 \
        }                                                            \
    } while(0)
#endif



/* Utility warning.  Tests a condition and if it fails prints out a message
 * with optional arguments
 * The do{...} while(0) simply enforces that a semicolon at the end
 */
#ifdef DEBUG
#define VWARN_MSG0(cnd, msg)                                              \
        do {                                                              \
                if( (cnd) == 0 ) {                                        \
                        Vnm_print(                                        \
                                2,                                        \
                                "%s:%d [%s()]: WARNING:\n"                \
                                "    Condition Failed (%s):\n    %s\n\n", \
                                __FILE__,                                 \
                                __LINE__,                                 \
                                __FUNCTION__,                             \
                                #cnd,                                     \
                                msg                                       \
                                );                                        \
                }                                                         \
        } while(0)

#define VWARN_MSG1(cnd, msg, arg0)                                        \
        do {                                                              \
                if( (cnd) == 0 ) {                                        \
                        char buff[1000];                                  \
                        snprintf(buff, 1000, msg, arg0);                  \
                        Vnm_print(                                        \
                                2,                                        \
                                "%s:%d [%s()]: WARNING:\n"                \
                                "    Condition Failed (%s):\n    %s\n\n", \
                                __FILE__,                                 \
                                __LINE__,                                 \
                                __FUNCTION__,                             \
                                #cnd,                                     \
                                buff                                      \
                                );                                        \
                }                                                         \
        } while(0)

#define VWARN_MSG2(cnd, msg, arg0, arg1)                                  \
        do {                                                              \
                if( (cnd) == 0 ) {                                        \
                        char buff[1000];                                  \
                        snprintf(buff, 1000, msg, arg0, arg1);            \
                        Vnm_print(                                        \
                                2,                                        \
                                "%s:%d [%s()]: WARNING:\n"                \
                                "    Condition Failed (%s):\n    %s\n\n", \
                                __FILE__,                                 \
                                __LINE__,                                 \
                                __FUNCTION__,                             \
                                #cnd,                                     \
                                buff                                      \
                                );                                        \
                }                                                         \
        } while(0)
#else
#define VWARN_MSG0(cnd, msg)                                              \
        do {                                                              \
                if( (cnd) == 0 ) {                                        \
                        Vnm_print(                                        \
                                2,                                        \
                                "[%s()]: WARNING:\n"                      \
                                "    %s\n\n",                             \
                                __FUNCTION__,                             \
                                msg                                       \
                                );                                        \
                }                                                         \
        } while(0)

#define VWARN_MSG1(cnd, msg, arg0)                                        \
        do {                                                              \
                if( (cnd) == 0 ) {                                        \
                        char buff[1000];                                  \
                        snprintf(buff, 1000, msg, arg0);                  \
                        Vnm_print(                                        \
                                2,                                        \
                                "[%s()]: WARNING:\n"                      \
                                "    %s\n\n",                             \
                                __FUNCTION__,                             \
                                buff                                      \
                                );                                        \
                }                                                         \
        } while(0)

#define VWARN_MSG2(cnd, msg, arg0, arg1)                                  \
        do {                                                              \
                if( (cnd) == 0 ) {                                        \
                        char buff[1000];                                  \
                        snprintf(buff, 1000, msg, arg0, arg1);            \
                        Vnm_print(                                        \
                                2,                                        \
                                "[%s()]: WARNING:\n"                      \
                                "    %s\n\n",                             \
                                __FUNCTION__,                             \
                                buff                                      \
                                );                                        \
                }                                                         \
        } while(0)
#endif

/* Utility Abort.  Prints a message with optional arugments and aborts */
#ifdef DEBUG
#define VABORT_MSG0(msg)                                      \
        do {                                                  \
            Vnm_print(2, "%s:%d [%s()]: ABORTING:\n"          \
                         "    %s\n\n",                        \
                      __FILE__, __LINE__, __FUNCTION__, msg); \
            abort();                                          \
        } while(0)

#define VABORT_MSG1(msg, arg)                                  \
        do {                                                   \
            char buff[1000];                                   \
            snprintf( buff, 1000, msg, arg );                  \
            Vnm_print(2, "%s:%d [%s()]: ABORTING:\n"           \
                         "    %s\n\n",                         \
                      __FILE__, __LINE__, __FUNCTION__, buff); \
            abort();                                           \
        } while(0)

#define VABORT_MSG2(msg, arg0, arg1)                           \
        do {                                                   \
            char buff[1000];                                   \
            snprintf( buff, 1000, msg, arg0, arg1);            \
            Vnm_print(2, "%s:%d [%s()]: ABORTING:\n"           \
                         "    %s\n\n",                         \
                      __FILE__, __LINE__, __FUNCTION__, buff); \
            abort();                                           \
        } while(0)
#else
#define VABORT_MSG0(msg)                        \
        do {                                    \
            Vnm_print(2, "%[%s()]: ABORTING:\n" \
                         "    %s\n\n",          \
                      __FUNCTION__, msg);       \
            abort();                            \
        } while(0)

#define VABORT_MSG1(msg, arg)                  \
        do {                                   \
            char buff[1000];                   \
            snprintf( buff, 1000, msg, arg );  \
            Vnm_print(2, "[%s()]: ABORTING:\n" \
                         "    %s\n\n",         \
                      __FUNCTION__, buff);     \
            abort();                           \
        } while(0)

#define VABORT_MSG2(msg, arg0, arg1)                \
        do {                                        \
            char buff[1000];                        \
            snprintf( buff, 1000, msg, arg0, arg1); \
            Vnm_print(2, "[%s()]: ABORTING:\n"      \
                         "    %s\n\n",              \
                      __FUNCTION__, buff);          \
            abort();                                \
        } while(0)
#endif



/* Utility expression printers.  Print the expression and its value */
#ifdef DEBUG
#define PRINT_INT(expr)                                             \
    do {                                                            \
        Vnm_print(2, "%s:%d [%s()]: %s == %d\n",                    \
                  __FILE__, __LINE__, __FUNCTION__, #expr, expr);   \
    } while(0)

#define PRINT_DBL(expr)                                           \
    do {                                                          \
        Vnm_print(2, "%s:%d [%s()]: %s == %f\n\n",                \
                  __FILE__, __LINE__, __FUNCTION__, #expr, expr); \
    } while(0)
#else
#define PRINT_INT(expr)
#define PRINT_DBL(expr)
#endif

#define VMALLOC(vmem, n, type) ((type*)Vmem_malloc(vmem, n, sizeof(type)))

#define VFREE(vmem, n, type, ptr) (Vmem_free(vmem, n, sizeof(type), (void **)&(ptr)))

#define VFILL(vec, n, val)                           \
    do {                                             \
        int fill_idx;                                \
        for (fill_idx = 0; fill_idx < n; fill_idx++) \
            vec[fill_idx] = val;                     \
    } while(0)

#define VCOPY(srcvec, dstvec, i, n) \
    do {                            \
        for (i = 0; i < n; i++)     \
            dstvec[i] = srcvec[i];  \
    } while(0)


char* wrap_text( char* str, int right_margin, int left_padding );

#define VAT(array, i) ((array)[(i) - 1])
#define RAT(array, i) ((array) + i - 1)

#endif /* #ifndef _VAPBSHAL_H_ */
