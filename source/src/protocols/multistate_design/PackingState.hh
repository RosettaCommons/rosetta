// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PackingState.hh
/// @brief for multistate design that is based on (re)use of packer information (RotamerSets and InteractionGraph)
/// @author ashworth

#ifndef INCLUDED_protocols_multistate_design_PackingState_hh
#define INCLUDED_protocols_multistate_design_PackingState_hh

#include <protocols/multistate_design/PackingState.fwd.hh>
#include <protocols/multistate_design/SingleState.hh>

#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <utility/vector0.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace multistate_design {

class PackingState : public SingleState {

public:
	PackingState();
	PackingState( core::pose::Pose const & pose_in, bool is_positive );
	~PackingState() override;

	void
	create_packer_data(
		core::scoring::ScoreFunctionCOP scorefxn,
		core::pack::task::PackerTaskCOP ptask
	);

	void share_packer_data_from( PackingState & other );

	void run_packer( utility::vector0<int> const & rot_to_pack );

	core::pack::task::PackerTaskCOP ptask() const;
	core::pack::rotamer_set::RotamerSetsCOP rotamersets() const;
	core::pack::interaction_graph::InteractionGraphBaseCOP ig() const;

protected:
	core::pack::rotamer_set::RotamerSetsOP rotamersets();
	core::pack::interaction_graph::InteractionGraphBaseOP ig();

private:
	// forbidden copy constructor: to prevent unexpected behavior re: nonconst pointer data
	PackingState( PackingState const & other_state );
	core::pack::task::PackerTaskCOP ptask_p_;
	core::pack::rotamer_set::RotamerSetsOP rotamersets_p_;
	core::pack::interaction_graph::InteractionGraphBaseOP ig_p_;
};

} // namespace multistate_design
} // namespace protocols

#endif
