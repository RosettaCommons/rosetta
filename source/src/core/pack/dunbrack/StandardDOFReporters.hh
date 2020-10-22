// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/dunbrack/StandardDOFReporter.hh
/// @brief   Class to measure the Standard DOFs
/// @author  Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_dunbrack_StandardDOFReporter_HH
#define INCLUDED_core_pack_dunbrack_StandardDOFReporter_HH

// Package headers
#include <core/pack/dunbrack/ResidueDOFReporter.hh>

// Project headers
#include <core/id/PartialAtomID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// C++ headers
#include <set>

namespace core {
namespace pack {
namespace dunbrack {

class MainchainTorsionReporter : public ResidueDOFReporter {
public:
	MainchainTorsionReporter(
		Size tor_ind,
		Size upper_tor_ind,
		Real neutral_val,
		bool flip_neutral_for_d_aa,
		bool flip_neutral_for_mirrored_);

	~MainchainTorsionReporter() override;

	Real
	get_dof( conformation::Residue const & rsd, pose::Pose const & pose ) const override;

	void
	insert_atoms_defining_dof(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		std::set< id::PartialAtomID > & atoms
	) const override;

private:
	Size const tor_ind_;
	Size const upper_tor_ind_;
	Real const neutral_val_;
	bool const flip_neutral_for_d_aa_;
	bool const flip_neutral_for_mirrored_;
};

class PeptideTorsionReporter : public ResidueDOFReporter {
public:
	PeptideTorsionReporter(
		Size tor_ind,
		Size upper_tor_ind,
		Real neutral_val,
		bool flip_neutral_for_mirrored_);

	~PeptideTorsionReporter() override;

	Real
	get_dof( conformation::Residue const & rsd, pose::Pose const & pose ) const override;

	void
	insert_atoms_defining_dof(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		std::set< id::PartialAtomID > & atoms
	) const override;

private:
	int bb_torsion_index_for_rsd(conformation::Residue const & rsd) const;

private:
	Size const tor_ind_;
	Size const upper_tor_ind_;
	Real const neutral_val_;
	bool const flip_neutral_for_mirrored_;
};


}
}
}


#endif
