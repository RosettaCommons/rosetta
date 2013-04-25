// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/devel/helixAssembly/AddResiduesRotamersOperation.ccAddResiduesRotamersOperation.cc
/// @brief 
/// @author Tim Jacobs

//Unit
#include<devel/helixAssembly/AddResiduesRotamerSetOperation.hh>

AddResiduesRotamerSetOperation::AddResiduesRotamerSetOperation(utility::vector1<core::conformation::ResidueOP> residues):
core::pack::rotamer_set::RotamerSetOperation(),
residues_(residues)
{}

core::pack::rotamer_set::RotamerSetOperationOP
AddResiduesRotamerSetOperation::clone() const
{
	return new AddResiduesRotamerSetOperation( *this );
}

void AddResiduesRotamerSetOperation::alter_rotamer_set(
			core::pose::Pose const & pose,
			core::scoring::ScoreFunction const & sfxn,
			core::pack::task::PackerTask const & ptask,
			core::graph::GraphCOP packer_neighbor_graph,
			core::pack::rotamer_set::RotamerSet & rotamer_set
	){

	for(core::Size i=1; i<=this->residues_.size(); ++i){
		rotamer_set.add_rotamer(*residues_[i]);
	}
}
