// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file PackingState.cc
/// @brief
/// @author ashworth

#include <protocols/multistate_design/PackingState.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility headers
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


namespace protocols {
namespace multistate_design {

using namespace core::pack::rotamer_set;
using namespace core::pack::task;
using namespace core::pack::interaction_graph;
using core::pack::pack_rotamers_setup;
using core::pack::pack_rotamers_run;

////////////////////////////////////////////////////////////////////////////////////////////////////
PackingState::PackingState()
: SingleState(),
	ptask_p_(/* 0 */),
	rotamersets_p_(/* 0 */),
	ig_p_(/* 0 */)
{}

PackingState::PackingState( core::pose::Pose const & pose, bool is_positive )
: SingleState( pose, is_positive ),
	ptask_p_(/* 0 */),
	rotamersets_p_(/* 0 */),
	ig_p_(/* 0 */)
{}

PackingState::~PackingState() = default;

PackerTaskCOP PackingState::ptask() const { return ptask_p_; }
RotamerSetsCOP PackingState::rotamersets() const { return rotamersets_p_; }
InteractionGraphBaseCOP PackingState::ig() const { return ig_p_; }
// protected (non-const)
RotamerSetsOP PackingState::rotamersets() { return rotamersets_p_; }
InteractionGraphBaseOP PackingState::ig() { return ig_p_; }

void
PackingState::create_packer_data(
	core::scoring::ScoreFunctionCOP scorefxn,
	PackerTaskCOP ptask
)
{
	ptask_p_ = ptask;
	rotamersets_p_ = core::pack::rotamer_set::RotamerSetsOP( new RotamerSets() );
	debug_assert( scorefxn && ptask_p_ && rotamersets_p_ );

	AnnealableGraphBaseOP ig;
	pack_rotamers_setup( nonconst_pose(), *scorefxn, ptask_p_, rotamersets_p_, ig );

	ig_p_ = utility::pointer::dynamic_pointer_cast<core::pack::interaction_graph::InteractionGraphBase>( ig );
	if ( ! ig_p_ ) {
		throw utility::excn::EXCN_Msg_Exception( "Interaction graph returned by pack_rotamers_setup is not a two-body interaction graph." );
	}
}

void
PackingState::share_packer_data_from( PackingState & other )
{
	debug_assert( other.ptask() && other.rotamersets() && other.ig() );
	ptask_p_ = other.ptask();
	rotamersets_p_ = other.rotamersets();
	ig_p_ = other.ig();
}

void
PackingState::run_packer( utility::vector0<int> const & rot_to_pack )
{
	debug_assert( ptask_p_ && rotamersets_p_ && ig_p_ );
	pack_rotamers_run( nonconst_pose(), ptask_p_, rotamersets_p_, ig_p_, rot_to_pack );
}

} // namespace multistate_design
} // namespace protocols
