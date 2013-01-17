/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief Specifies the PDE definition for PMG to solve
 *  @version $Id:
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

#ifndef _MYPDE_H_
#define _MYPDE_H_

#include "math.h"

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "pmgc/mypdec.h"

#define MAXIONS    50
#define MAXPOLY   50
#define ZSMALL     1.0e-20
#define ZLARGE     1.0e20
#define SINH_MIN -85.0
#define SINH_MAX  85.0

/// @todo  Remove dependencies on global variables
double v1, v2, v3, conc1, conc2, conc3, vol, relSize;
int nion;
double charge[MAXIONS];
double sconc[MAXIONS];

#define Na 6.022045000e-04

/** @brief   Set up the ionic species to be used in later calculations.  This
 *           must be called before any other of the routines in this file.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Nathan Baker [Original]
 *
 *  @note    Replaces mypdefinitlpbe from mypde.f
 */
VEXTERNC void Vmypdefinitlpbe(
        int *tnion,       ///< The number if ionic species
        double *tcharge,  ///< The charge in electrons
        double *tsconc    /**< Prefactor for conterion Bolzmann distribution
                           *   terms.  Basically a scaled concentration
                           *     -(ion concentration/bulkIonicStrength)/2
                           */
        );



/** @brief   Set up the ionic species to be used in later calculations.  This
 *           must be called before any other of the routines in this file.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Nathan Baker [Original]
 *
 *  @note    Replaces mypdefinitnpbe from mypde.f
 */
VEXTERNC void Vmypdefinitnpbe(
        int *tnion,       ///< The number if ionic species
        double *tcharge,  ///< The charge in electrons
        double *tsconc    /**< Prefactor for conterion Bolzmann distribution
                           *   terms.  Basically a scaled concentration
                           *     -(ion concentration/bulkIonicStrength)/2
                           */
        );



/** @brief   Set up the ionic species to be used in later calculations.  This
 *           must be called before any other of the routines in this file.
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Nathan Baker [Original]
 *
 *  @note    Replaces mypdefinitsmpbe from mypde.f
 */
VEXTERNC void Vmypdefinitsmpbe(
        int *tnion,       ///< The number if ionic species
        double *tcharge,  ///< The charge in electrons
        double *tsconc,   /**< Prefactor for conterion Bolzmann distribution
                           *   terms.  Basically a scaled concentration
                           *     -(ion concentration/bulkIonicStrength)/2
                           */
        double *smvolume, ///< @todo: Doc
        double *smsize    ///< @todo: Doc
        );



/** @brief   Define the nonlinearity (vector version)
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces c_vec from mypde.f
 */
VEXTERNC void Vc_vec(
        double *coef, ///< @todo: Doc
        double *uin,  ///< @todo: Doc
        double *uout, ///< @todo: Doc
        int *nx,      ///< @todo: Doc
        int *ny,      ///< @todo: Doc
        int *nz,      ///< @todo: Doc
        int *ipkey    ///< @todo: Doc
        );



/** @brief   Define the derivative of the nonlinearity (vector version)
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces dc_vec from mypde.f
 */
VEXTERNC void Vdc_vec(
        double *coef, ///< @todo: Doc
        double *uin,  ///< @todo: Doc
        double *uout, ///< @todo: Doc
        int    *nx,   ///< @todo: Doc
        int    *ny,   ///< @todo: Doc
        int    *nz,   ///< @todo: Doc
        int    *ipkey ///< @todo: Doc
        );

VEXTERNC void Vdc_vecpmg(
        double *coef, ///< @todo: Doc
        double *uin,  ///< @todo: Doc
        double *uout, ///< @todo: Doc
        int    *nx,   ///< @todo: Doc
        int    *ny,   ///< @todo: Doc
        int    *nz,   ///< @todo: Doc
        int    *ipkey ///< @todo: Doc
        );

VEXTERNC void Vdc_vecsmpbe(
        double *coef, ///< @todo: Doc
        double *uin,  ///< @todo: Doc
        double *uout, ///< @todo: Doc
        int    *nx,   ///< @todo: Doc
        int    *ny,   ///< @todo: Doc
        int    *nz,   ///< @todo: Doc
        int    *ipkey ///< @todo: Doc
        );



/** @brief   Define the nonlinearity (vector version)
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces c_vecpmg from mypde.f
 */
VEXTERNC void Vc_vecpmg(
        double *coef, ///< @todo: Doc
        double *uin,  ///< @todo: Doc
        double *uout, ///< @todo: Doc
        int *nx,      ///< @todo: Doc
        int *ny,      ///< @todo: Doc
        int *nz,      ///< @todo: Doc
        int *ipkey    ///< @todo: Doc
        );



/** @brief   Define the nonlinearity (vector version)
 *  @ingroup PMGC
 *  @author  Tucker Beck [C Translation], Michael Holst [Original]
 *
 *  @note    Replaces c_vecpmg from mypde.f
 */
VEXTERNC void Vc_vecsmpbe(
        double *coef, ///< @todo: Doc
        double *uin,  ///< @todo: Doc
        double *uout, ///< @todo: Doc
        int *nx,      ///< @todo: Doc
        int *ny,      ///< @todo: Doc
        int *nz,      ///< @todo: Doc
        int *ipkey    ///< @todo: Doc
        );

#endif /* _MYPDE_H_ */
