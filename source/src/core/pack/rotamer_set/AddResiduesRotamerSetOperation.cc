// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file core/pack/rotamer_set/AddResiduesRotamersOperation.cc
/// @brief Transform a vector of residues onto the backbone of the appropriate pose position and
/// add the transformed residues to the rotamer set
/// @author Tim Jacobs

//Unit
#include<core/pack/rotamer_set/AddResiduesRotamerSetOperation.hh>

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
	core::graph::GraphCOP /*packer_neighbor_graph*/,
	core::pack::rotamer_set::RotamerSet & rotamer_set
){
	core::Size seqnum = (core::Size) rotamer_set.resid();
	for(core::Size i=1; i<=this->residues_.size(); ++i){
		core::conformation::ResidueOP cur_res = residues_[i]->clone();
		core::conformation::Residue const & existing_residue=pose.residue(seqnum);
		cur_res->place(existing_residue, pose.conformation());
		cur_res->seqpos(seqnum);
		cur_res->chain(existing_residue.chain());
		cur_res->copy_residue_connections_from(existing_residue);
		rotamer_set.add_rotamer(*cur_res);
	}
}

} //rotamer_set
} //pack
} //core
