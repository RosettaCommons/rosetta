// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/rotamer_set/AddResiduesRotamersOperation.cc
/// @brief Transform a vector of residues onto the backbone of the appropriate pose position and
/// add the transformed residues to the rotamer set
/// @author Tim Jacobs

//Unit
#include<core/pack/rotamer_set/AddResiduesRotamerSetOperation.hh>
#include <utility>

namespace core {
namespace pack {
namespace rotamer_set {

AddResiduesRotamerSetOperation::AddResiduesRotamerSetOperation(
	utility::vector1<core::conformation::ResidueOP> const & residues
):
	core::pack::rotamer_set::RotamerSetOperation(),
	residues_(residues)
{}

core::pack::rotamer_set::RotamerSetOperationOP
AddResiduesRotamerSetOperation::clone() const
{
	return core::pack::rotamer_set::RotamerSetOperationOP( new AddResiduesRotamerSetOperation( *this ) );
}

void AddResiduesRotamerSetOperation::alter_rotamer_set(
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & /*sfxn*/,
	core::pack::task::PackerTask const & /*ptask*/,
	utility::graph::GraphCOP /*packer_neighbor_graph*/,
	core::pack::rotamer_set::RotamerSet & rotamer_set
){
	auto seqnum = (core::Size) rotamer_set.resid();
	for ( core::Size i=1; i<=this->residues_.size(); ++i ) {
		core::conformation::ResidueOP cur_res = residues_[i]->clone();
		core::conformation::Residue const & existing_residue=pose.residue(seqnum);
		cur_res->place(existing_residue, pose.conformation());
		cur_res->seqpos(seqnum);
		cur_res->chain(existing_residue.chain());
		cur_res->copy_residue_connections_from(existing_residue);
		rotamer_set.add_rotamer_into_existing_group(*cur_res);
	}
}

} //rotamer_set
} //pack
} //core
