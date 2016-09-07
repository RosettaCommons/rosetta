// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/rosetta_scripts/AdditionalOutputWrapper.cc
/// @brief This mover wraps another mover to obtain additional output from it via regular call to
///   apply(). A new instance is created for each call to get_additional_output() on this wrapper.
/// @author Luki Goldschmidt (lugo@uw.edu)

// Unit headers
#include <protocols/rosetta_scripts/AdditionalOutputWrapper.hh>
#include <protocols/rosetta_scripts/AdditionalOutputWrapperCreator.hh>

// C/C++ headers
#include <iostream>
#include <string>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Project headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/rosetta_scripts/RosettaScriptsParser.hh>

// Package headers
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <boost/foreach.hpp>


namespace protocols {
namespace rosetta_scripts {

static THREAD_LOCAL basic::Tracer TR( "protocols.rosetta_scripts.AdditionalOutputWrapper" );

using namespace protocols::moves;

////////////////////////////////////////////////////////////////////////

std::string AdditionalOutputWrapperCreator::mover_name()
{
	return "AdditionalOutputWrapper";
}

std::string AdditionalOutputWrapperCreator::keyname() const
{
	return mover_name();
}

MoverOP AdditionalOutputWrapperCreator::create_mover() const
{
	return MoverOP( new AdditionalOutputWrapper() );
}

////////////////////////////////////////////////////////////////////////

AdditionalOutputWrapper::AdditionalOutputWrapper() :
	Mover( "AdditionalOutputWrapper" ),
	mover_tag_(/* NULL */),
	rosetta_scripts_tag_(/* NULL */),
	reference_pose_(/* NULL */),
	max_poses_(0),
	n_poses_(0)
{
}

std::string AdditionalOutputWrapper::get_name() const
{
	return AdditionalOutputWrapperCreator::mover_name();
}

/// @brief Process all input poses (provided pose and from previous mover)
void AdditionalOutputWrapper::apply(core::pose::Pose& pose)
{
	reference_pose_ = core::pose::PoseOP( new core::pose::Pose(pose) );
	generate_pose(pose);
	++n_poses_;
}

/// @brief Hook for multiple pose putput to JD2 or another mover
core::pose::PoseOP AdditionalOutputWrapper::get_additional_output()
{
	if ( !reference_pose_ ) {
		return nullptr;
	}
	if ( (max_poses_ > 0) && (n_poses_ >= max_poses_) ) {
		return nullptr;
	}

	core::pose::PoseOP new_pose( new core::pose::Pose(*reference_pose_) );
	generate_pose(*new_pose);
	++n_poses_;

	return new_pose;
}

void AdditionalOutputWrapper::generate_pose(core::pose::Pose & pose)
{
	// Empty objects... may not work...
	basic::datacache::DataMap data;
	protocols::filters::Filters_map filters;
	protocols::moves::Movers_map movers;
	MoverOP mover(nullptr);

	if ( !mover && rosetta_scripts_tag_ ) {
		protocols::rosetta_scripts::RosettaScriptsParser parser;
		mover = parser.parse_protocol_tag( pose, rosetta_scripts_tag_ );
	}

	if ( !mover && mover_tag_ ) {
		mover = MoverFactory::get_instance()->newMover(mover_tag_, data, filters, movers, pose );
	}

	runtime_assert( mover != nullptr );
	mover->apply(pose);

	if ( ! pose.data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		// Add new output tag
		using basic::datacache::DataCache_CacheableData;
		std::ostringstream tag;
		tag << name_ << "_" << (n_poses_+1);
		pose.data().set(
			core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG,
			DataCache_CacheableData::DataOP( new basic::datacache::CacheableString( tag.str() ) )
		);
	}

	protocols::moves::Mover::set_last_move_status(mover->get_last_move_status());
}

void AdditionalOutputWrapper::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {

	if ( tag->hasOption("name") ) {
		name_ = tag->getOption<std::string>("name");
	}

	if ( tag->hasOption("max_output_poses") ) {
		max_poses_ = tag->getOption<int>("max_output_poses", 0);
	}

	try {

		// Children of tag are movers
		BOOST_FOREACH ( utility::tag::TagCOP const curr_tag, tag->getTags() ) {
			// Try instantiating first mover from tag to test parsing
			if ( curr_tag->getName() == "ROSETTASCRIPTS" ) {
				// Treat subtag as a ROSETTASCRIPTS protocol
				protocols::rosetta_scripts::RosettaScriptsParser parser;
				protocols::moves::MoverOP mover( parser.parse_protocol_tag( curr_tag ) );
				rosetta_scripts_tag_ = curr_tag;
			} else {
				// Treat subtag as a regular mover tag
				std::string name = curr_tag->getOption<std::string>("name");
				protocols::moves::MoverOP new_mover(
					protocols::moves::MoverFactory::get_instance()->
					newMover(curr_tag, data, filters, movers, pose)
				);
				mover_tag_ = curr_tag;
			}
			// Only first mover used -- add warning when multiple defined?
			break;
		}

		if ( !mover_tag_ && !rosetta_scripts_tag_ ) {
			throw utility::excn::EXCN_Msg_Exception("No mover or ROSETTASCRIPTS tag found.");
		}

	} catch( utility::excn::EXCN_Msg_Exception const & e ) {
		std::string my_name( tag->getOption<std::string>("name") );
		throw utility::excn::EXCN_Msg_Exception("Exception in AdditionalOutputWrapper with name \"" + my_name + "\": " + e.msg());
	}
}

} //rosetta_scripts
} //protocols
