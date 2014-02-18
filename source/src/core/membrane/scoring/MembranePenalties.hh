// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file           MembranePenalties.hh
///
/// @brief          Membrane Penalties
/// @details		This class's methods evaluate scoring penalties for the membrane region
///                 This will later get instantiated in MembraneSearch and MembranePtoential
///
/// @note           Updated Documentation - Documented from Vladmir et al. 2006
///
/// @author         Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_scoring_MembranePenalties_hh
#define INCLUDED_core_membrane_scoring_MembranePenalties_hh

// Unit Headers
#include <core/membrane/scoring/MembranePenalties.fwd.hh>

// Project Headers
#include <core/membrane/util/definitions.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace membrane {
namespace scoring {

/// @brief 	Class: MembranePenalties
/// @detail This class contains methods for evaluating and applying penalties
///			for tm helix orientations relative to pose membrane center and normal.
class MembranePenalties : public utility::pointer::ReferenceCount {

public: // methods

	/// @brief      Constructor
	/// @details	Creates an instance of the Membrane Penalties object
	MembranePenalties();

	/// @brief      Evlauate TM Projection Penalty
	/// @details    Penalty for TM helixes that do not follow projection from
    ///             normal and center
	void
    tm_projection_penalty(
                         std::string desc,
                         std::string chain,
                         Vector const & normal,
                         Vector const & center,
                         core::Real & tm_proj
                        );

	/// @brief      Evaluate Non Membrane Helix Penalty
	/// @details    Penalty for Predicted TM Helixes that do not span the membrane
	void
    non_helix_in_membrane_penalty(
                       std::string desc,
                       std::string chain,
                       Vector const & normal,
                       Vector const & center,
                       core::Real & non_helix_pen
                       );

	/// @brief      Evaluate Termini Penalty
	/// @details    Penalty for termini positions contradictory to membrane
    ///             spanning
	void
    termini_penalty(
                         std::string desc,
                         std::string chain,
                         Vector const & normal,
                         Vector const & center,
                         Real & termini_pen
                         );

}; // class MembranePenalties

} // scoring
} // membrane
} // core

#endif // INCLUDED_core_membrane_scoring_MembranePenalties_hh