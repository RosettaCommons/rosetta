// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/RotamerDump/RotamerDumpMover.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_RotamerDump_RotamerDumpMOver_HH
#define INCLUDED_protocols_RotamerDump_RotamerDumpMOver_HH

#include <protocols/RotamerDump/RotamerDumpMover.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>

#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace RotamerDump {

/// @brief a mover for dumping interaction graph information.  The information is appended to the current job.
class RotamerDumpMover : public protocols::moves::Mover
{
public:

	RotamerDumpMover(core::pack::task::TaskFactoryOP task_factory, core::scoring::ScoreFunctionOP score_function);
	void apply(core::pose::Pose & pose) override;
	std::string get_name() const override ;

private:
	/// @brief dump the one body energy table
	std::string get_onebody_energy_table(core::pack::interaction_graph::InteractionGraphBaseOP ig,core::pack::rotamer_set::RotamerSetsOP rotamer_sets);

	/// @brief dump the two body energy table
	std::string get_twobody_energy_table(core::pack::interaction_graph::InteractionGraphBaseOP ig,core::pack::rotamer_set::RotamerSetsOP rotamer_sets);
	/// @brief dump the xyz coordinates of every atom in every rotamer
	std::string get_xyz_coord_table(core::pack::rotamer_set::RotamerSetsOP rotamer_sets);
	/// @brief return the interaction nodes selected by the annealer.  This function does not modify the pose.
	std::string get_annealer_pick_table(core::pack::interaction_graph::InteractionGraphBaseOP ig, core::pack::rotamer_set::RotamerSetsOP rotamer_sets, core::pose::Pose & pose, core::pack::task::PackerTaskCOP task);

	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP score_function_;

};

}
}

#endif
