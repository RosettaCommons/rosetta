// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/devel/helixAssembly/AddResiduesRotamersOperation.hhAddResiduesRotamersOperation.hh
/// @brief 
/// @author Tim Jacobs

#ifndef ADDRESIDUESROTAMERSOPERATION_HH_
#define ADDRESIDUESROTAMERSOPERATION_HH_

//Core
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/graph/Graph.hh>
#include <core/conformation/Residue.hh>

///@brief Adds in rotamers from a list of Residues,
class AddResiduesRotamerSetOperation : public core::pack::rotamer_set::RotamerSetOperation
{
public:

	AddResiduesRotamerSetOperation(utility::vector1<core::conformation::ResidueOP> residues);

	core::pack::rotamer_set::RotamerSetOperationOP
	clone() const;

	void
	alter_rotamer_set(
			core::pose::Pose const & pose,
			core::scoring::ScoreFunction const & sfxn,
			core::pack::task::PackerTask const & ptask,
			core::graph::GraphCOP packer_neighbor_graph,
			core::pack::rotamer_set::RotamerSet & rotamer_set
	);

private:
	utility::vector1<core::conformation::ResidueOP> residues_;
};

#endif /* ADDRESIDUESROTAMERSOPERATION_HH_ */
