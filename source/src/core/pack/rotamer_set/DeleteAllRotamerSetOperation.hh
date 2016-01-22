// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/core/pack/rotamer_set/DeleteAllRotamerSetOperation.hh
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_core_pack_rotamer_set_DeleteAllRotamerSetOperation_HH
#define INCLUDED_core_pack_rotamer_set_DeleteAllRotamerSetOperation_HH

//Unit headers
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/rotamer_set/DeleteAllRotamerSetOperation.fwd.hh>

//Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>

namespace core {
namespace pack {
namespace rotamer_set {

/// @brief Adds in rotamers from a list of Residues,
class DeleteAllRotamerSetOperation : public RotamerSetOperation
{
public:

	RotamerSetOperationOP
	clone() const;

	void
	alter_rotamer_set(
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		task::PackerTask const & ptask,
		graph::GraphCOP packer_neighbor_graph,
		RotamerSet & rotamer_set
	);
};

} //namespace rotamer_set
} //namespace pack
} //namespace core

#endif
