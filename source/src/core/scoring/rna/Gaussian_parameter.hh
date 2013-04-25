// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file   core/scoring/rna/Gaussian_parameter.hh
/// @brief  very simple class that records amplitude, center, and width of a Gaussian
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_Gaussian_parameter_HH
#define INCLUDED_core_scoring_rna_Gaussian_parameter_HH

#include <core/types.hh>

// Project headers
#include <utility/pointer/ReferenceCount.hh>

// C++
namespace core {
namespace scoring {
namespace rna {

class Gaussian_parameter {
public:
	Real amplitude, center, width;

	Gaussian_parameter ( Real const amplitude_in, Real const center_in, Real const width_in ):
		amplitude( amplitude_in ),
		center   ( center_in ),
		width    ( width_in )
	{}

	//
	Gaussian_parameter &
	operator=( Gaussian_parameter const & src )
	{
		amplitude = src.amplitude;
		center = src.center;
		width = src.width;
		return *this;
	}


};


}
}
}

#endif
