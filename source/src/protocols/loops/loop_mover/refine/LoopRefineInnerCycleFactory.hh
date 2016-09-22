// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/loops/loop_mover/refine/LoopRefineInnerCycleFactory.hh
/// @brief Factory for creating LoopRefineInnerCycle objects
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )


#ifndef INCLUDED_protocols_loops_LoopRefineInnerCycle_HH
#define INCLUDED_protocols_loops_LoopRefineInnerCycle_HH

// Unit headers
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycleFactory.fwd.hh>

// Package headers
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.fwd.hh>
#include <protocols/loops/loop_mover/refine/LoopRefineInnerCycle.fwd.hh>

// Project headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <map>
#include <string>

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

// Prefix all entries in this enum with "IC_" (InnerCycle) to avoid confusion with class names
enum LoopRefineInnerCycleName {
	// Individual LoopRefineInnerCycles
	IC_SmallMinCCDTrial = 1,
	IC_ShearMinCCDTrial,
	IC_RepackTrial,

	// Pre-made algorithms
	IC_RefineCCDStandard,
	number_of_loop_refine_inner_cycle_names = IC_RefineCCDStandard
};

/// Create LoopMover Reporters
class LoopRefineInnerCycleFactory : public utility::SingletonBase< LoopRefineInnerCycleFactory >
{
public:
	friend class utility::SingletonBase< LoopRefineInnerCycleFactory >;

public:
	// Warning this is not called because of the singleton pattern
	virtual ~LoopRefineInnerCycleFactory();

	/// @brief Create a LoopRefineInnerCycle giving it a pointer to the data it needs to function
	LoopRefineInnerCycleOP create_inner_cycle(
		LoopRefineInnerCycleName type_name,
		LoopMover_Refine_CCDAP loop_mover,
		moves::MonteCarloOP mc,
		core::scoring::ScoreFunctionOP scorefxn,
		core::pack::task::TaskFactoryOP tf
	) const;

private: // methods
	// Private constructor - singleton
	LoopRefineInnerCycleFactory();
	LoopRefineInnerCycleFactory(const LoopRefineInnerCycleFactory & src); // unimplemented

	LoopRefineInnerCycleFactory const &
	operator=( LoopRefineInnerCycleFactory const & ); // unimplemented

	LoopRefineInnerCycleOP make_inner_cycle_from_string_name( std::string const & name ) const;
	void setup_known_types();


private:

	utility::vector1< utility::vector1< std::string > > loop_refine_inner_cycle_name_to_string_;

	// TODO: Add a std::map< std::string, LoopRefineInnerCycleName > to allow for commandline or RosettaScripts based selection.
	//       core/scoring/ScoreTypeManager.cc has an example of setting something like this up

};

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols

#endif // INCLUDED_protocols_loops_LoopRefineInnerCycle_HH
