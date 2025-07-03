// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/modifications/Reprotonate.hh
/// @brief  Heuristically reprotonate a molecule to be in neutral pH
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_core_chemical_modifications_Reprotonate_hh
#define INCLUDED_core_chemical_modifications_Reprotonate_hh

#include <core/chemical/modifications/Reprotonate.fwd.hh>
#include <core/chemical/modifications/ChemistryBase.hh>

#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/types.hh>

namespace core {
namespace chemical {
namespace modifications {


class Reprotonate : public ChemistryBase
{
public:
	Reprotonate() :
		ChemistryBase(class_name())
	{}

	/// @brief Search the provided ResidueType for protons which should be added/removed
	/// in order to conform to the protonation state in neutral (pH 7-ish) aqueous solution.
	/// This is done heuristically, assuming that protons have already been added.
	/// The geometry of any existing hydrogen atoms will be adjusted, but not any heavy atoms.
	void apply(MutableResidueType & res ) override;

	core::chemical::VDVDMapping
	get_mapping() const override {
		return mapping_;
	}

	static std::string class_name();

private:
	/// @brief Find a doubly-bonded oxygen or sulfur attached to neighbor
	VD
	find_electron_sink(MutableResidueType const & restype, VD neighbor) const;

	/// @brief Is the provided able to participate in conjugation?
	bool
	is_conjugatable(MutableResidueType const & restype, VD atom) const;

	/// @brief Remove hydrogens from specified atom.
	/// Returns the number of hydrogens removed
	core::Size
	remove_hydrogens(MutableResidueType & restype, VD atom);

private:

	core::chemical::VDVDMapping mapping_;

};

}
}
}

#endif
