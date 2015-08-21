// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/sc/ShapeComplementarityCalculator.hh
/// @brief  Private headers for the Shape Complementarity Calculator
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)

#ifndef INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_Private_hh
#define INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_Private_hh

#include <string>
#include <stdio.h>
#include <stdarg.h>

#include <utility/excn/Exceptions.hh>

namespace core {
namespace scoring {
namespace sc {

////////////////////////////////////////////////////////////
// Defs from sc source

enum __ATTEN_ENUM__ {
	ATTEN_BLOCKER = 1,
	ATTEN_2  = 2,
	ATTEN_BURIED_FLAGGED = 5,
	ATTEN_6 = 6
};

#define MAX_SUBDIV 100

#if  (!defined MIN)
#define MIN(a,b) ((a) < (b) ? (a): (b))
#define MAX(a,b) ((a) > (b) ? (a): (b))
#endif

#define ABS(a) (((a) < 0) ? (-a) : (a))

#define PI numeric::NumericTraits< core::Real >::pi()

////////////////////////////////////////////////////////////

class ShapeComplementarityCalculatorException : public utility::excn::EXCN_Msg_Exception {

public:
	std::string error;

	ShapeComplementarityCalculatorException(const char *err, ...) :
		utility::excn::EXCN_Msg_Exception( std::string() )
	{
		va_list p;
		char buf[256];
		va_start(p, err);
		vsnprintf(buf, sizeof(buf), err, p);
		va_end(p);
		error = buf;
		add_msg(error);
	}
};

} //namespace sc
} //namespace filters
} //namespace protocols

#endif
