// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/dunbrack/PeptoidDOFReporter.hh
/// @brief   Class to measure the Peptoid DOFs
/// @author  Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_dunbrack_PeptoidDOFReporter_hh
#define INCLUDED_core_pack_dunbrack_PeptoidDOFReporter_hh

// Package headers
#include <core/pack/dunbrack/ResidueDOFReporter.hh>

// Project headers
#include <core/types.hh>
#include <core/id/PartialAtomID.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace core {
namespace pack {
namespace dunbrack {

class PeptoidOmegaReporter : public ResidueDOFReporter {
public:
	PeptoidOmegaReporter(Real const neutral_omega);
	~PeptoidOmegaReporter() override;

	Real
	get_dof(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const override;

	void
	insert_atoms_defining_dof(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		std::set< id::PartialAtomID > & atoms
	) const override;

private:
	Real const peptoid_neutral_omega_;

};


class PeptoidGeneralDOFReporter : public ResidueDOFReporter {
public:
	PeptoidGeneralDOFReporter(
		Size const mainchain_torsion,
		Size const n_tors,
		Real const neutral_val
	);
	~PeptoidGeneralDOFReporter() override;

	Real
	get_dof(
		conformation::Residue const & rsd,
		pose::Pose const & pose
	) const override;

	void
	insert_atoms_defining_dof(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		std::set< id::PartialAtomID > & atoms
	) const override;

private:
	Size const mainchain_torsion_;
	Size const n_tors_;
	Real const neutral_val_;
};

}
}
}


#endif
