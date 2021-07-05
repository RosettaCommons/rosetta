// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_set/PyRotamerEliminator.hh
/// @brief  Wrappers that allow the user to eliminate rotamers from pyrosetta
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_core_pack_rotamer_set_PyRotamerEliminator_hh
#define INCLUDED_core_pack_rotamer_set_PyRotamerEliminator_hh

// Unit Headers
#include <core/pack/rotamer_set/PyRotamerEliminator.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <utility/graph/Graph.fwd.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility Headers

#include <utility/vector1.hh>

#include <functional>

#if defined(PYROSETTA) && !defined(PYROSETTA_BINDER)
#include "pybind11/functional.h"
#endif //PYROSETTA

namespace core {
namespace pack {
namespace rotamer_set {

class PyRotamerEliminator : public RotamerSetsOperation
{
public:

	// Function takes all of the arguments for alter_rotamer_sets and returns a vector1<bool> of length rotamer_sets.nrotamers(). Elements with TRUE will be ELIMINATED
	using EliminatorFunction = std::function< utility::vector1< bool >( pose::Pose const &, scoring::ScoreFunctionCOP, task::PackerTaskCOP, utility::graph::GraphCOP, RotamerSets const & ) >;


	RotamerSetsOperationOP
	clone() const override;

	void
	alter_rotamer_sets(
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		task::PackerTask const & ptask,
		utility::graph::GraphCOP packer_neighbor_graph,
		RotamerSets & rotamer_sets
	) override;

	void
	set_eliminator_function( EliminatorFunction const & f ){
		function_ = f;
	}

private:

	EliminatorFunction function_; //inits as empty

};

class PyRotamerEliminatorTaskOperation : public task::operation::TaskOperation {
public:

	task::operation::TaskOperationOP clone() const override {
		return utility::pointer::make_shared< PyRotamerEliminatorTaskOperation >( *this );
	};

	void apply(
		core::pose::Pose const &,
		core::pack::task::PackerTask & task
	) const override {
		runtime_assert( eliminator_ != nullptr );
		task.append_rotamersets_operation( eliminator_ );
	}

	void set_eliminator( PyRotamerEliminatorOP e ){
		eliminator_ = e;
	}

private:
	PyRotamerEliminatorOP eliminator_ = nullptr;
};

}
}
}

#endif
