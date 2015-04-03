// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /rosetta/rosetta_source/src/devel/sewing/DeleteAllRotamerSetOperation.hh
/// @brief 
/// @author Tim Jacobs

#ifndef INCLUDED_devel_sewing_DeleteAllRotamerSetOperation_HH
#define INCLUDED_devel_sewing_DeleteAllRotamerSetOperation_HH

//Core
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/graph/Graph.hh>
#include <core/conformation/Residue.hh>

namespace devel {
namespace sewing {

/// @brief Adds in rotamers from a list of Residues,
class DeleteAllRotamerSetOperation : public core::pack::rotamer_set::RotamerSetOperation
{
public:

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
};

} //sewing namespace
} //devel namespace

#endif /* ADDRESIDUESROTAMERSOPERATION_HH_ */
