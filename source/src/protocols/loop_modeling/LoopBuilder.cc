// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/loop_modeling/types.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/LoopBuilder.hh>

// Protocol headers
#include <protocols/loop_modeling/samplers/KicSampler.hh>
#include <protocols/loop_modeling/refiners/RepackingRefiner.hh>
#include <protocols/loop_modeling/refiners/LocalMinimizationRefiner.hh>
#include <protocols/kinematic_closure/perturbers/RamaPerturber.hh>
#include <protocols/kinematic_closure/pivot_pickers/EndToEndPivots.hh>
#include <protocols/kinematic_closure/solution_pickers/FilteredSolutions.hh>

namespace protocols {
namespace loop_modeling {

LoopBuilder::LoopBuilder() : LoopMover() {}

LoopBuilder::LoopBuilder(Loop const & loop, ScoreFunctionOP score_function) 
		: LoopMover(loop, score_function) {

	using protocols::loop_modeling::samplers::KicSampler;
	using protocols::loop_modeling::samplers::KicSamplerOP;
	using protocols::loop_modeling::refiners::RepackingRefiner;
	using protocols::loop_modeling::refiners::LocalMinimizationRefiner;
	using protocols::kinematic_closure::perturbers::RamaPerturber;
	using protocols::kinematic_closure::pivot_pickers::EndToEndPivots;
	using protocols::kinematic_closure::solution_pickers::FilteredSolutions;

	KicSamplerOP sampler = new KicSampler;

	//sampler->add_opener(new IdealGeometryOpener);
	sampler->add_perturber(new RamaPerturber);
	sampler->set_pivot_picker(new EndToEndPivots);
	sampler->set_solution_picker(new FilteredSolutions(false)); 

	add_task(sampler);
	add_task(new RepackingRefiner);
	add_task(new LocalMinimizationRefiner);
}

}
}
