/**
 *  @file    vcap.c
 *  @ingroup Vcap
 *  @author  Nathan Baker
 *  @brief   Class Vcap methods
 *  @version $Id$
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

#include "vcap.h"

VPUBLIC double Vcap_exp(double x, int *ichop) {

    /* The two chopped arguments */
    if (x > EXPMAX) {
       (*ichop) = 1;
       return VEXP(EXPMAX);
    } else if (x < EXPMIN) {
       (*ichop) = 1;
       return VEXP(EXPMIN);
    }

    /* The normal EXP */
    (*ichop) = 0;
    return VEXP(x);
}

VPUBLIC double Vcap_sinh(double x, int *ichop) {

    /* The two chopped arguments */
    if (x > EXPMAX) {
       (*ichop) = 1;
       return VSINH(EXPMAX);
    } else if (x < EXPMIN) {
       (*ichop) = 1;
       return VSINH(EXPMIN);
    }

    /* The normal SINH */
    (*ichop) = 0;
    return VSINH(x);
}

VPUBLIC double Vcap_cosh(double x, int *ichop) {

    /* The two chopped arguments */
    if (x > EXPMAX) {
       (*ichop) = 1;
       return VCOSH(EXPMAX);
    } else if (x < EXPMIN) {
       (*ichop) = 1;
       return VCOSH(EXPMIN);
    }

    /* The normal COSH */
    (*ichop) = 0;
    return VCOSH(x);
}

