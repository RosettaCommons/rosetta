/**
 *  @file    vunit.h
 *  @ingroup Vunit
 *  @author  Nathan Baker
 *  @brief   Contains a collection of useful constants and conversion factors
 *  @author  Nathan A. Baker
 *  @version $Id$
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

/** @defgroup Vunit Vunit class
 *  @brief    Collection of constants and conversion factors
 */

#ifndef _VUNIT_H_
#define _VUNIT_H_

/** @brief   Multiply by this to convert J to cal 
 *  @ingroup Vunit */
#define Vunit_J_to_cal		4.1840000e+00

/** @brief   Multiply by this to convert cal to J
 *  @ingroup Vunit */
#define Vunit_cal_to_J		2.3900574e-01

/** @brief   Multiply by this to convert amu to kg
 *  @ingroup Vunit */
#define Vunit_amu_to_kg 	1.6605402e-27

/** @brief   Multiply by this to convert kg to amu
 *  @ingroup Vunit */
#define Vunit_kg_to_amu 	6.0221367e+26

/** @brief   Multiply by this to convert ec to C
 *  @ingroup Vunit */
#define Vunit_ec_to_C		1.6021773e-19

/** @brief   Multiply by this to convert C to ec
 *  @ingroup Vunit */
#define Vunit_C_to_ec		6.2415065e+18

/** @brief   Charge of an electron in C
 *  @ingroup Vunit */
#define Vunit_ec		1.6021773e-19

/** @brief   Boltzmann constant
 *  @ingroup Vunit */
#define Vunit_kb		1.3806581e-23

/** @brief   Avogadro's number
 *  @ingroup Vunit */
#define Vunit_Na		6.0221367e+23

/** @brief   Pi
 *  @ingroup Vunit */
#define Vunit_pi		VPI

/** @brief   Vacuum permittivity
 *  @ingroup Vunit */
#define Vunit_eps0		8.8541878e-12

/** @brief \f${e_c}^2/\AA\f$ in ESU units => kcal/mol 
 *  @ingroup Vunit */
#define Vunit_esu_ec2A		3.3206364e+02

/** @brief \f$k_b\f$ in ESU units => kcal/mol 
 *  @ingroup Vunit */
#define Vunit_esu_kb            1.9871913e-03

#endif /* ifndef _VUNIT_H_ */
