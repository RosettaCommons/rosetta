// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mm/MMBondLengthScore.hh
/// @brief  Molecular mechanics bond length score class
/// @author Frank DiMaio (based on Colin Smith's MMBondAngle potential)


#ifndef INCLUDED_core_scoring_mm_MMBondLengthScore_hh
#define INCLUDED_core_scoring_mm_MMBondLengthScore_hh

// Unit headers
/// you cannot #include yourself #include <core/scoring/mm/MMBondLengthScore.hh>
#include <core/scoring/mm/MMBondLengthLibrary.hh>

// Project headers
#include <core/chemical/MMAtomTypeSet.fwd.hh>

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


class MMBondLengthScore : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~MMBondLengthScore();

	MMBondLengthScore();
	MMBondLengthScore( MMBondLengthLibrary const & mmtl );

	Real score( Real Kb, Real b0, Real d ) const;
	Real score( mm_bondlength_atom_pair mm_atomtype_set, Real d ) const;

	Real dscore( Real Kb, Real b0, Real d ) const;
	Real dscore( mm_bondlength_atom_pair bondlength_atom_set, Real d ) const;

private:

	/// @brief Local MMBondLengthLibrary for looking up bond angle parameters
	MMBondLengthLibrary const & mm_bondlength_library_;

};

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_MMBondLengthScore_HH
