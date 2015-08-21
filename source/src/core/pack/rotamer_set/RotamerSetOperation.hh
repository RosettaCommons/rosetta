// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/RotamerSetOperation.hh
/// @brief  rotamer set operation class
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_RotamerSetOperation_hh
#define INCLUDED_core_pack_rotamer_set_RotamerSetOperation_hh

// Unit Headers
#include <core/pack/rotamer_set/RotamerSetOperation.fwd.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <core/types.hh>
#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {

/// @brief RotamerSetOperations are able to modify the contents of a RotamerSet
/// within the packer immediately after rotamer creation.
///
/// @li They are handed into the packer through a packer task; each ResidueLevelTask
///     keeps its own list of rotamer set operations.
///
/// @li Each RotamerSet, within the build_rotamers() method will iterate across
///     the RotamerSetOperation objects, and call alter_rotamer_set on each
///     element in the ResidueLevelTask's list of rotamer_set_operations.
///
/// @li RotamerSetOperations are visited in the order in which they are appended
///     to the ResidueLevelTasks's operation list.
///
/// @li RotamerSetOperations are unable to correlate changes to rotamer sets across
///     multiple rotamer sets -- they have access to only a single RotamerSet object
///     in their alter_rotamer_set operations.  For correlated alterations to rotamer
///     sets, the (as of yet undefined) RotamerSetsOperation class should be used.

class RotamerOperation : public utility::pointer::ReferenceCount
{
public:
	RotamerOperation();
	virtual ~RotamerOperation();

	virtual
	bool
	operator() (
		conformation::ResidueOP rotamer,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		pack::task::ResidueLevelTask const & rtask,
		graph::GraphCOP packer_neighbor_graph,
		pack::dunbrack::ChiSetOP chi_set
	) = 0;

};

class RotamerSetOperation : public utility::pointer::ReferenceCount
{
public:
	RotamerSetOperation();
	virtual ~RotamerSetOperation();

	virtual
	RotamerSetOperationOP
	clone() const = 0;

	virtual
	void
	alter_rotamer_set(
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		task::PackerTask const & ptask,
		graph::GraphCOP packer_neighbor_graph,
		RotamerSet & rotamer_set
	) = 0;

	virtual
	Real
	increase_packer_residue_radius(
		pose::Pose const & pose,
		task::PackerTaskCOP the_task,
		core::Size residue_in
	) const;

};

}
}
}

#endif
