// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMTorsionScore.hh
/// @brief  Molecular mechanics torsion score class
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_scoring_mm_MMTorsionScore_hh
#define INCLUDED_core_scoring_mm_MMTorsionScore_hh

// Unit headers
#include <core/scoring/mm/MMTorsionScore.fwd.hh>
#include <core/scoring/mm/MMTorsionLibrary.fwd.hh>
#include <core/scoring/mm/MMTorsionLibrary.hh>

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

/// @brief Calculates scores of mm torsion paramater sets given an angle
///
/// @details
///
///
class MMTorsionScore : public utility::pointer::ReferenceCount
{

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~MMTorsionScore();

	/// @brief Default ctor
	MMTorsionScore();

	/// @brief Alternate ctor that inintalizes class with given MMTorsionLibrary
	MMTorsionScore( MMTorsionLibrary const & mmtl );

	/// @brief Returns energy given an mm_torsion_atom_quad and an angle in radians
	Real score( mm_torsion_atom_quad dihedral_atom_set, Real angle ) const;

	/// @brief Returns derivative of the energy given an mm_torsion_atom_quad and an angle in radians
	Real dscore( mm_torsion_atom_quad dihedral_atom_set, Real angle ) const;

private:

	/// @brief Local MMTorsionLibrary for looking up torsion parameters
	MMTorsionLibrary const & mm_torsion_library_;

};

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_mm_torsion_score_HH
