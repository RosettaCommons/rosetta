// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mm/MMBondAngleScore.hh
/// @brief  Molecular mechanics bond angle score class
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_core_scoring_mm_MMBondAngleScore_hh
#define INCLUDED_core_scoring_mm_MMBondAngleScore_hh

// Unit headers
#include <core/scoring/mm/MMBondAngleScore.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.fwd.hh>

//#include <core/scoring/ScoringManager.hh>

#include <core/types.hh>

// Utility header
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

namespace core {
namespace scoring {
namespace mm {

/// @brief Calculates scores of mm bond angle paramater sets given an angle
///
/// @details
///
///
class MMBondAngleScore : public utility::pointer::ReferenceCount
{

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~MMBondAngleScore();

	/// @brief Default ctor
	MMBondAngleScore();

	/// @brief Alternate ctor that inintalizes class with given MMBondAngleLibrary
	MMBondAngleScore( MMBondAngleLibrary const & mmtl );

	/// @brief Returns energy given an mm_bondangle_atom_tri and an angle in radians
	Real score( Real Ktheta, Real theta0, Real angle ) const;

	/// @brief Returns energy given an mm_bondangle_atom_tri and an angle in radians
	Real score( mm_bondangle_atom_tri bondangle_atom_set, Real angle ) const;

	/// @brief Returns a derivative given an mm_bondangle_atom_trie and an angle in radians
	Real
	dscore( Real Ktheta, Real theta0, Real angle ) const;

	/// @brief Returns a derivative given an mm_bondangle_atom_trie and an angle in radians
	Real
	dscore( mm_bondangle_atom_tri mm_atomtype_set, Real angle ) const;

private:

	/// @brief Local MMBondAngleLibrary for looking up bond angle parameters
	MMBondAngleLibrary const & mm_bondangle_library_;

};

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_MMBondAngleScore_HH
