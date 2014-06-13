// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMLJScore.hh
/// @brief  Molecular mechanics lj score class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_scoring_mm_MMLJScore_hh
#define INCLUDED_core_scoring_mm_MMLJScore_hh

// Unit headers
#include <core/scoring/mm/MMLJScore.fwd.hh>
// AUTO-REMOVED #include <core/scoring/mm/MMLJLibrary.hh>

// Project headers
// AUTO-REMOVED #include <core/chemical/MMAtomTypeSet.fwd.hh>

// AUTO-REMOVED #include <core/scoring/EnergyMap.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>

#include <core/types.hh>

// Utility header
#include <utility/vector1.hh>

// C++ headers
// AUTO-REMOVED #include <string>
// AUTO-REMOVED #include <map>

#include <core/scoring/mm/MMLJLibrary.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>


namespace core {
namespace scoring {
namespace mm {

/// @brief Calculates scores of mm lj paramater sets given two mm atom types,
/// the path distance of the 2 atoms and actual distance between the two atoms
///
/// @details blah
///
///
class MMLJScore : public utility::pointer::ReferenceCount
{

public:

  /// @brief Default ctor
  MMLJScore();

  /// @brief Alternate ctor that inintalizes class with given MMLJLibrary
  MMLJScore( MMLJLibrary const & mmljl );

	virtual ~MMLJScore();

	/// @brief blah
	MMLJLibrary const &
	mm_lj_library() const
	{ return mm_lj_library_; }

  /// @brief blah
  Energy
	score( Size atom1, Size atom2, Size path_distance, Real distance ) const;

  /// @brief blah
	Energy
	deriv_score( Size atom1, Size atom2, Size path_distance, Real distance ) const;

  /// @brief blah
	Real
	min_dist( Size atom1, Size atom2, Size path_distance ) const;

private:

  /// @brief Local MMLJLibrary for looking up lj parameters
  MMLJLibrary const & mm_lj_library_;
};

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_mm_lj_score_HH
