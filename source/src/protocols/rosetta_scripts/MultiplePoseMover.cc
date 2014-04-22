// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rosetta_scripts/MultiplePoseMover.cc
/// @brief	This mover accepts multiple poses from a previous mover,
///   performs selection using a provided pose selector,
///   applies contained ROSETTASCRIPTS protocol (ParsedProtocol),
///   and output multiple poses to the next mover of JD2.
/// @author Luki Goldschmidt (lugo@uw.edu)

// Unit headers
#include <protocols/rosetta_scripts/MultiplePoseMover.hh>
#include <protocols/rosetta_scripts/MultiplePoseMoverCreator.hh>
#include <protocols/rosetta_scripts/PoseSelector.fwd.hh>
#include <protocols/rosetta_scripts/PoseSelector.hh>
#include <protocols/rosetta_scripts/PoseSelectorFactory.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <protocols/filters/Filter.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/util.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <boost/foreach.hpp>


namespace protocols {
namespace rosetta_scripts {

static basic::Tracer TR( "protocols.rosetta_scripts.MultiplePoseMover" );

////////////////////////////////////////////////////////////////////////

std::string MultiplePoseMoverCreator::mover_name()
{
	return "MultiplePoseMover";
}

std::string MultiplePoseMoverCreator::keyname() const
{
	return mover_name();
}

protocols::moves::MoverOP MultiplePoseMoverCreator::create_mover() const
{
	return new MultiplePoseMover();
}

////////////////////////////////////////////////////////////////////////

MultiplePoseMover::MultiplePoseMover() :
	protocols::moves::Mover( "MultiplePoseMover" ),
	max_poses_(0),
	rosetta_scripts_tag_(NULL),
	selector_tag_(NULL),
	previous_mover_(NULL),
	selected_poses_i_(0)
{
}

std::string MultiplePoseMover::get_name() const
{
	return MultiplePoseMoverCreator::mover_name();
}

///@brief Process all input poses (provided pose and from previous mover)
void MultiplePoseMover::apply(core::pose::Pose& pose)
{
	moves::MoverStatus status(protocols::moves::MS_SUCCESS);

	using namespace core;
	using namespace core::pose;

	utility::vector1 < core::pose::PoseOP > selected_poses_for_processing;

	selected_poses_.clear();
	poses_.clear();

	// Clone provided pose since we may be returning a different based on Selector
	poses_.push_back(pose.clone());

	// Collect all remaining input poses from previous mover
	if(previous_mover_) {
		TR << "Obtaining additional poses from previous mover " << previous_mover_->get_name() << std::endl;
		do {
			PoseOP next_pose( previous_mover_->get_additional_output() );
			if(!next_pose)
				break;
			poses_.push_back(next_pose);
		} while(poses_.size() < max_poses_ || !max_poses_);
	}

	TR << "Collected input poses: " << poses_.size() << std::endl;
	BOOST_FOREACH( PoseOP p, poses_ ) {
		TR << "\t" << p->sequence() << std::endl;
	}

	// Process selection criteria
	if(selectors_.size() > 0) {
		utility::vector1 < bool > selected_poses_by_selectors;

		// Apply selectors
		BOOST_FOREACH( PoseSelectorCOP const selector, selectors_ ) {
			selected_poses_by_selectors = selector->select_poses(poses_);
			// TODO: How do we handle multiple selectors? AND? OR?
			break;
		}

		if(selectors_.size() > 1) {
			TR << selectors_.size() << " selectors defined; only the first selector is currently used!" << std::endl;
		}

		// Make a new vector of PoseOP's for selected poses for easier handling
		core::Size i = 1;
		for(utility::vector1 < PoseOP >::iterator it = poses_.begin(); it != poses_.end(); ++it) {
			if(selected_poses_by_selectors[i])
				selected_poses_for_processing.push_back(*it);
			++i;
		}

		// Apply movers to selected poses
		TR << "Selected poses for processing: " << selected_poses_for_processing.size() << std::endl;
		BOOST_FOREACH( PoseOP p, selected_poses_for_processing ) {
			TR << "\t" << p->sequence() << std::endl;
		}

	} else {
		// No selector specified -- select all poses
		selected_poses_for_processing = poses_;
	}

	{
		if(rosetta_scripts_tag_) {
			// Run sub-protocol:

			// Collect additional poses first rather than adding them to selected_poses_ right away
			// as it may invalidate the iterator pointer when memory location changes
			utility::vector1 < PoseOP > additional_poses;

			BOOST_FOREACH( PoseOP p, selected_poses_for_processing ) {
				TR << "Applying mover to pose: " << p->sequence() << std::endl;
				if(process_pose(*p, additional_poses)) {
					// Processing successful (from sub-protocol and its filters)
					selected_poses_.push_back( p );
				}
			}

			if(additional_poses.size() > 0) {
				BOOST_FOREACH( PoseOP p, additional_poses ) {
					selected_poses_.push_back( p );
				}
				TR << additional_poses.size() << " additional poses obtained; total output poses: " << selected_poses_.size() << std::endl;
			}
		} else {
			// No sub-protocol provided, use all selected poses
			selected_poses_ = selected_poses_for_processing;
		}
	}

	// Return first selected pose, additional can be obtained via get_additional_output()
	selected_poses_i_ = 1;
	if(selected_poses_.size() > 0) {
		pose = *selected_poses_[selected_poses_i_];
		++selected_poses_i_;
	} else {
		TR << "No poses selected or processed succesfully by sub-mover. Setting status to FAIL_RETRY." << std::endl;
		status = protocols::moves::FAIL_RETRY;
	}

	protocols::moves::Mover::set_last_move_status(status);
}

///@brief Process a single input pose by the RosettaScripts mover
bool MultiplePoseMover::process_pose( core::pose::Pose & pose, utility::vector1 < core::pose::PoseOP > & additional_poses )
{
	if(!rosetta_scripts_tag_)
		return false;

	protocols::rosetta_scripts::RosettaScriptsParser parser;

	// rosetta_scripts_tag_ has been pre-parsed in parse_my_tag() so no parsing exception should be thrown here
	protocols::moves::MoverOP mover( parser.parse_protocol_tag( pose, rosetta_scripts_tag_ ) );
	if(!mover) {
		TR << "Failed to parse protocol? This should not happen. Not applying protocol to pose." << std::endl;
		return false;
	}

	mover->apply(pose);

	if(mover->get_last_move_status() != protocols::moves::MS_SUCCESS) {
		TR << "Sub-mover reported failure; passing status on." << std::endl;
		return false;
	}

	// Get additional poses from protocol mover;
	// these are selected for output by this mover
	core::pose::PoseOP p;
	while( ( p = mover->get_additional_output() ) ) {
		additional_poses.push_back(p);
	}

	return true;
}

///@brief Hook for multiple pose putput to JD2 or another mover
core::pose::PoseOP MultiplePoseMover::get_additional_output()
{
	TR.Debug << "get_additional_output, last index = " << selected_poses_i_ << std::endl;
	if(selected_poses_i_ >= 1 && selected_poses_i_ <= selected_poses_.size())
		return selected_poses_[selected_poses_i_++];
	return NULL;
}

void MultiplePoseMover::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {

	if(tag->hasOption("max_input_poses"))
		max_poses_ = tag->getOption<int>("max_input_poses", 0);

	try {

		// ROSETTASCRIPTS tag (optional)
		if(tag->hasTag("ROSETTASCRIPTS")) {
			rosetta_scripts_tag_ = tag->getTag("ROSETTASCRIPTS");
			// Try parsing the ROSETTASCRIPTS block to avoid tripping die_for_unaccessed_options()
			protocols::rosetta_scripts::RosettaScriptsParser parser;
			protocols::moves::MoverOP mover( parser.parse_protocol_tag( rosetta_scripts_tag_ ) );
		}

		// SELECT tag (optional)
		if(tag->hasTag("SELECT")) {
			selectors_.clear();
			TagCOP select_tag( tag->getTag("SELECT") );
			BOOST_FOREACH( TagCOP const curr_tag, select_tag->getTags() ) {
				PoseSelectorOP new_selector(
					PoseSelectorFactory::get_instance()->
						newPoseSelector( curr_tag, data, filters, movers, pose )
				);
				runtime_assert( new_selector );
				selectors_.push_back( new_selector );
				TR << "Defined selector of type " << curr_tag->getName() << std::endl;
			}
		}

		// Warn if no ROSETTASCRIPTS protocol and no SELECTOR (i.e. null mover)
		if(selectors_.size() < 1 && !rosetta_scripts_tag_) {
			std::string my_name( tag->getOption<std::string>("name") );
			TR.Warning << "Neither a ROSETTASCRIPTS protocol nor a SELECT statement specified in MultiplePoseMover with name \"" << my_name << "\". This mover has no effect. Are you sure this is what you intended?" << std::endl;
		}

		// TODO: Should we complain here is there are tags specified that we don't understand?

	} catch( utility::excn::EXCN_Msg_Exception const & e ) {
		std::string my_name( tag->getOption<std::string>("name") );
		throw utility::excn::EXCN_Msg_Exception("Exception in MultiplePoseMover with name \"" + my_name + "\": " + e.msg());
	}

	TR << "MultiplePoseMover\n";
	if(max_poses_ != 0)
		TR << "\tMax input poses: " << max_poses_ << std::endl;
	// TODO: other parsing summary
}

} //rosetta_scripts
} //protocols
