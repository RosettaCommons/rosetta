/** @defgroup Vcap Vcap class
 *  @brief  Collection of routines which cap certain exponential and hyperbolic
 *          functions
 *  @note   These routines are based on FORTRAN code by Mike Holst
 */

/**
 *  @file     vcap.h
 *  @ingroup  Vcap
 *  @brief    Contains declarations for class Vcap
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

#ifndef _VCAP_H_
#define _VCAP_H_

#include "apbscfg.h"

/** @brief   Maximum argument for exp(), sinh(), or cosh()
 *  @ingroup Vcap
 */
#define EXPMAX  85.00

/** @brief   Minimum argument for exp(), sinh(), or cosh()
 *  @ingroup Vcap
 */
#define EXPMIN -85.00

#include "maloc/maloc.h"

/** @brief   Provide a capped exp() function
 *
 *           If the argument x of Vcap_exp() exceeds EXPMAX or EXPMIN, then we
 *           return exp(EXPMAX) or exp(EXPMIN) rather than exp(x).
 *
 *  @note    Original FORTRAN routine from PMG library by Mike Holst
 *           Original notes:
 *           to control overflow in the hyperbolic and exp functions, note
 *           that the following are the argument limits of the various
 *           functions on various machines after which overflow occurs:
 *           Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
 *           sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
 *           dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
 *
 *  @ingroup Vcap
 *  @author  Nathan Baker (based on FORTRAN code by Mike Holst)
 *  @return  exp(x) or capped equivalent
 */
VEXTERNC double Vcap_exp(
        double x, /**< Argument to exp() */
        int *ichop /**< Set to 1 if function capped, 0 otherwise */
        );


/** @brief   Provide a capped sinh() function
 *
 *           If the argument x of Vcap_sinh() exceeds EXPMAX or EXPMIN, then we
 *           return sinh(EXPMAX) or sinh(EXPMIN) rather than sinh(x).
 *
 *  @note    Original FORTRAN routine from PMG library by Mike Holst
 *           Original notes:
 *           to control overflow in the hyperbolic and exp functions, note
 *           that the following are the argument limits of the various
 *           functions on various machines after which overflow occurs:
 *           Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
 *           sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
 *           dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
 *
 *  @ingroup Vcap
 *  @author  Nathan Baker (based on FORTRAN code by Mike Holst)
 *  @return  sinh(x) or capped equivalent
 */
VEXTERNC double Vcap_sinh(
        double x, /**< Argument to sinh() */
        int *ichop /**< Set to 1 if function capped, 0 otherwise */
        );

/** @brief   Provide a capped cosh() function
 *
 *           If the argument x of Vcap_cosh() exceeds EXPMAX or EXPMIN, then we
 *           return cosh(EXPMAX) or cosh(EXPMIN) rather than cosh(x).
 *
 *  @note    Original FORTRAN routine from PMG library by Mike Holst
 *           Original notes:
 *           to control overflow in the hyperbolic and exp functions, note
 *           that the following are the argument limits of the various
 *           functions on various machines after which overflow occurs:
 *           Convex C240, Sun 3/60, Sun SPARC, IBM RS/6000:
 *           sinh, cosh, exp:     maximal argument (abs value) =  88.0d0
 *           dsinh, dcosh, dexp:  maximal argument (abs value) = 709.0d0
 *
 *  @ingroup Vcap
 *  @author  Nathan Baker (based on FORTRAN code by Mike Holst)
 *  @return  cosh(x) or capped equivalent
 */
VEXTERNC double Vcap_cosh(
        double x, /**< Argument to cosh() */
        int *ichop /**< Set to 1 if function capped, 0 otherwise */
        );

#endif    /* ifndef _VCAP_H_ */
