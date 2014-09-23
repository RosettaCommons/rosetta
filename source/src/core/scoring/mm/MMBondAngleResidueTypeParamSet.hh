// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondAngleResidueTypeParamSet.hh
/// @brief  Class to store bond angle parameters for a set of ResidueTypes
/// @author Colin A. Smith (colin.smith@ucsf.edu)


#ifndef INCLUDED_core_scoring_mm_MMBondAngleResidueTypeParamSet_hh
#define INCLUDED_core_scoring_mm_MMBondAngleResidueTypeParamSet_hh

// Unit headers
#include <core/scoring/mm/MMBondAngleResidueTypeParamSet.fwd.hh>

// Project headers
#include <core/id/AtomID.fwd.hh>
// AUTO-REMOVED #include <core/scoring/mm/MMBondAngleLibrary.hh>
// AUTO-REMOVED #include <core/scoring/mm/MMBondAngleResidueTypeParam.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.hh>

// Utility header
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

#include <core/chemical/ResidueType.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/scoring/mm/MMBondAngleResidueTypeParam.fwd.hh>
#include <utility/vector1_bool.hh>
#include <string>


namespace core {
namespace scoring {
namespace mm {

class MMBondAngleResidueTypeParamSet  : public utility::pointer::ReferenceCount
{

public:
	/// @brief ctor
	MMBondAngleResidueTypeParamSet();

	/// @brief copy ctor
	MMBondAngleResidueTypeParamSet(
		MMBondAngleResidueTypeParamSet const & src
	);

	virtual ~MMBondAngleResidueTypeParamSet();

	/// @brief get bond angle library for newly created ResidueTypeParam objects
	MMBondAngleLibrary const *
	mm_bondangle_library() const
	{
		return mm_bondangle_library_;
	}

	/// @brief set bond angle library for newly created ResidueTypeParam objects
	void
	mm_bondangle_library(
		MMBondAngleLibrary const * value
	)
	{
		mm_bondangle_library_ = value;
	}

	/// @brief get use of ResidueType theta0 in newly created ResidueTypeParam objects
	bool
	use_residue_type_theta0() const
	{
		return use_residue_type_theta0_;
	}

	/// @brief set use of ResidueType theta0 in newly created ResidueTypeParam objects
	void
	use_residue_type_theta0(
		bool value
	)
	{
		use_residue_type_theta0_ = value;
	}

	/// @brief get central atoms to score in newly created ResidueTypeParam objects
	utility::vector1<std::string> const &
	central_atoms_to_score() const
	{
		return central_atoms_to_score_;
	}

	/// @brief set central atoms to score in newly created ResidueTypeParam objects
	void
	central_atoms_to_score(
		utility::vector1<std::string> const & value
	)
	{
		central_atoms_to_score_ = value;
	}

	/// @brief lookup a param object for a given ResidueType
	MMBondAngleResidueTypeParam const &
	get(
		core::chemical::ResidueType const & residue_type
	);

	/// @brief lookup a param object for a given ResidueType, does not auto-create
	MMBondAngleResidueTypeParam const *
	get(
		core::chemical::ResidueType const & residue_type
	) const;

	/// @brief lookup Ktheta and theta0 for any bond angle in a conformation
	void
	lookup(
		core::conformation::Conformation const & conformation,
		core::id::AtomID const & atomid1,
		core::id::AtomID const & atomid2,
		core::id::AtomID const & atomid3,
		core::Real & Ktheta,
		core::Real & theta0
	);

	/// @brief lookup Ktheta and theta0 for any bond angle in a conformation
	void
	lookup(
		core::conformation::Conformation const & conformation,
		core::id::AtomID const & atomid1,
		core::id::AtomID const & atomid2,
		core::id::AtomID const & atomid3,
		core::Real & Ktheta,
		core::Real & theta0
	) const;

private:

	/// @brief lookup Ktheta and theta0 for any bond angle in a conformation
	void
	lookup(
		core::conformation::Conformation const & conformation,
		core::id::AtomID const & atomid1,
		core::id::AtomID const & atomid2,
		core::id::AtomID const & atomid3,
		MMBondAngleResidueTypeParam const & residue_type_param,
		core::Real & Ktheta,
		core::Real & theta0
	) const;

	/// @brief bond angle library for deriving new ResidueTypeParam objects
	MMBondAngleLibrary const * mm_bondangle_library_;

	/// @brief use ResidueType defined theta0 in newly created ResidueTypeParam objects
	bool use_residue_type_theta0_;

	/// @brief only score certain central atoms in newly created ResidueTypeParam objects
	utility::vector1<std::string> central_atoms_to_score_;

	/// @brief mapping from ResidueType name to ResidueTypeParam
	std::map<std::string, MMBondAngleResidueTypeParam> reside_type_param_map_;

};

/// @brief extract a MMBondAngleResidueTypeParamSet from a ScoreFunction, returning NULL if none exists
MMBondAngleResidueTypeParamSetCOP
mm_bond_angle_residue_type_param_set(
	core::scoring::ScoreFunction const & scorefxn
);

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_MMBondAngleResidueTypeParamSet_HH
