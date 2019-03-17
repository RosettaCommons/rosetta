// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

/// @file    core/io/rcsb/ExperimentalTechnique.hh
/// @brief   Enumeration definition for ExperimentalTechnique and helper function declarations.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @author  Matthew O'Meara


#ifndef INCLUDED_core_io_ExperimentalTechnique_HH
#define INCLUDED_core_io_ExperimentalTechnique_HH


// C++ headers
#include <list>
#include <string>


namespace core {
namespace io {
namespace rcsb {

enum ExperimentalTechnique {
	UNKNOWN_EXPDTA = 0,
	// Experimental Techniques for spec version 3.3
	X_RAY_DIFFRACTION = 1,
	FIBER_DIFFRACTION,
	NEUTRON_DIFFRACTION,
	ELECTRON_CRYSTALLOGRAPHY,
	ELECTRON_MICROSCOPY,
	SOLID_STATE_NMR,
	SOLUTION_NMR,
	SOLUTION_SCATTERING,

	THEORETICAL_MODEL,
	ExperimentalTechnique_max_current = THEORETICAL_MODEL,

	// Additional non-obsolete code.
	EPR,

	// Obsolete technique codes
	ELECTRON_DEFRACTION,
	CRYO_ELECTRON_MICROSCOPY,
	SOLUTION_SCATTERING_THEORETICAL_MODEL,
	FLORECENCE_TRANSFER,
	NMR, // Note the qualifying information is parsed not stored

	ExperimentalTechnique_max = NMR
};

typedef std::list< ExperimentalTechnique > ExperimentalTechniques;

std::string experimental_technique_to_string( ExperimentalTechnique technique );

ExperimentalTechnique string_to_experimental_technique( std::string const & technique );

}  // namespace rcsb
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_ExperimentalTechnique_HH
