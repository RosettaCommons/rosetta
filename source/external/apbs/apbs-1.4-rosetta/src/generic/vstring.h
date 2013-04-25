/** @defgroup Vstring Vstring class
 *  @brief    Provides a collection of useful non-ANSI string functions
 */

/**
 *  @file     vstring.h
 *  @ingroup  Vstring
 *  @brief    Contains declarations for class Vstring
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

#ifndef _VSTRING_H_
#define _VSTRING_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"

/** @brief   Case-insensitive string comparison (BSD standard)
 *  @ingroup Vstring
 *  @author  Copyright (c) 1988-1993 The Regents of the University of
 *           California.  Copyright (c) 1995-1996 Sun Microsystems, Inc.
 *  @note    Copyright (c) 1988-1993 The Regents of the University of
 *           California.  Copyright (c) 1995-1996 Sun Microsystems, Inc.
 *  @param   s1   First string for comparison
 *  @param   s2   Second string for comparison
 *  @return  An integer less than, equal to, or greater than zero if s1 is
 *           found,  respectively,  to  be  less  than, to match, or be greater
 *           than s2. (Source:  Linux man pages)
 */
VEXTERNC int Vstring_strcasecmp(const char *s1, const char *s2);

/** @brief   A modified sscanf that examines the complete string
 *  @ingroup Vstring
 *  @author  Todd Dolinsky
 *  @param   tok   The string to examine
 *  @return  1 if the entire string is an integer, 0 if otherwise.
 */
VEXTERNC int Vstring_isdigit(const char *tok);

/** Creates a wrapped and indented string from an input string
 *  @author Tucker Beck
 *  @ingroup Vstring
 *  @note    This function allocates a new string, so be sure to free it!
 */
VEXTERNC char* Vstring_wrappedtext(
        const char* str,  /**< The input string to wrap and indent          */
        int right_margin, /**< The number of characters to the right margin */
        int left_padding  /**< The number of characters in the left indent  */
        );

#endif    /* ifndef _VSTRING_H_ */
