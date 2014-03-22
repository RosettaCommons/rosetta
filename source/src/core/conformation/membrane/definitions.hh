// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/conformation/membrane/definitions.hh
///
/// @brief 		Definitions for membrane protein modeling data
/// @details 	Useful customd ata structures for implementing membrane protein
///				modeling framework
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_conformation_membrane_definitions_hh
#define INCLUDED_core_conformation_membrane_definitions_hh

// Unit headers
#include <core/conformation/membrane/definitions.fwd.hh>

// Package Headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>
#include <utility/vector0.hh>

// Platform headers
#include <platform/types.hh>

// C++ Headers
#include <cstddef>
#include <cstdlib>
#include <string>

namespace core {
namespace conformation {
namespace membrane {

/// @brief      Struct EmbedConfigInfo
/// @details    Stores data required for properly initializing membrane protein embedding
///             data either from user definition or external methods such as searc/score
struct EmbedConfigInfo : utility::pointer::ReferenceCount {

    // Tag for methods
    std::string normal_tag;
    std::string center_tag;

	// Normal and center with respect to membrane definition
	core::Vector normal;
	core::Vector center;

	// Chain depth for non-spanning chains
	core::Real depth;

}; // Struct EmbedConfigInfo

/// @brief      Struct Embedding Search Info
/// @details	Stores required parameters for performing a search and score for embedding (protein)
///             in a membrane
struct EmbedSearchParams : utility::pointer::ReferenceCount {

	// Search for center
	bool center_search;
    
	// Search for normal
	bool normal_search;
    
    // Apply Penalties
    bool penalties;
    
	// Specify a maximum deviation from initial center
	core::Size center_max_delta;

	// Specify a starting angle, max, and max deviation
	core::Size normal_start_angle;
	core::Size normal_max_angle;
	core::Size normal_delta_angle;

	// Membrane normal and center mag
	core::Size center_mag;
	core::Size normal_mag;

	// Normal cycles
	core::Size normal_cycles;
    
    // Interpolation in scoring
    bool no_interpolate_Mpair;

};

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_definitions_hh

