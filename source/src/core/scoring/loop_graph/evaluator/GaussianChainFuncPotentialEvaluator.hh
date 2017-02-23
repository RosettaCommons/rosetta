// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/loop_graph/evaluator/GaussianChainFuncPotentialEvaluator.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_loop_graph_GaussianChainFuncPotentialEvaluator_HH
#define INCLUDED_core_scoring_loop_graph_GaussianChainFuncPotentialEvaluator_HH

#include <core/scoring/loop_graph/evaluator/LoopClosePotentialEvaluator.hh>
#include <core/scoring/loop_graph/evaluator/GaussianChainFuncPotentialEvaluator.fwd.hh>
#include <core/scoring/loop_graph/LoopCycle.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>

namespace core {
namespace scoring {
namespace loop_graph {
namespace evaluator {

class GaussianChainFuncPotentialEvaluator: public LoopClosePotentialEvaluator {

public:

	//constructor
	GaussianChainFuncPotentialEvaluator( LoopCycle const & loop_cycle,
		core::Real const & loop_fixed_cost,
		pose::Pose const & pose );

	//destructor
	~GaussianChainFuncPotentialEvaluator();

public:

	virtual
	void
	get_f1_f2( Vector & f1, Vector & f2, bool const takeoff ) const;

private:

	void
	initialize( LoopCycle const & loop_cycle,
		core::pose::Pose const & pose );

private:

	core::Real const rna_gaussian_variance_per_residue_;
	core::Real const protein_gaussian_variance_per_residue_;
	core::scoring::func::FuncOP func_;

};

} //evaluator
} //loop_graph
} //scoring
} //core

#endif
