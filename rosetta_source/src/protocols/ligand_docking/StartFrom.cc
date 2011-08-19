// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com)
///

// Unit Headers
#include <protocols/ligand_docking/StartFrom.hh>
#include <protocols/ligand_docking/StartFromCreator.hh>

#include <protocols/geometry/RB_geometry.hh>
#include <protocols/moves/RigidBodyMover.hh>

// Utility Headers
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>

#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <core/types.hh>

//Auto Headers
#include <core/pose/util.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static basic::Tracer start_from_tracer("protocols.ligand_docking.ligand_options.Start_from", basic::t_debug);

std::string
StartFromCreator::keyname() const
{
	return StartFromCreator::mover_name();
}

protocols::moves::MoverOP
StartFromCreator::create_mover() const {
	return new StartFrom;
}

std::string
StartFromCreator::mover_name()
{
	return "StartFrom";
}

///@brief
StartFrom::StartFrom():
		//utility::pointer::ReferenceCount(),
		Mover("StartFrom"),
		chain_(""),
		starting_points_(){}

StartFrom::StartFrom(StartFrom const & that):
		//utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		chain_(that.chain_),
		starting_points_(that.starting_points_)
{}

StartFrom::~StartFrom() {}

protocols::moves::MoverOP StartFrom::clone() const {
	return new StartFrom( *this );
}

protocols::moves::MoverOP StartFrom::fresh_instance() const {
	return new StartFrom;
}

std::string StartFrom::get_name() const{
	return "StartFrom";
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
StartFrom::parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & /*datamap*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "StartFrom" ){
		utility_exit_with_message("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) utility_exit_with_message("'StartFrom' mover requires chain tag");

	chain_ = tag->getOption<std::string>("chain");

	utility::vector0< utility::tag::TagPtr >::const_iterator begin=tag->getTags().begin();
	utility::vector0< utility::tag::TagPtr >::const_iterator end=tag->getTags().end();

	for(; begin != end; ++begin){
		utility::tag::TagPtr child_tag= *begin;
		std::string name= child_tag->getName();
		if( name == "features"){
			std::cout << "found features tag with type '" << child_tag->getOption<std::string>("type") << "'" << std::endl;

		} else if( name != "Coordinates")
			utility_exit_with_message("StartFrom's children are always of type 'Coordinates'");
		if ( ! child_tag->hasOption("x") ) utility_exit_with_message("'StartFrom' mover requires 'x' coordinates option");
		if ( ! child_tag->hasOption("y") ) utility_exit_with_message("'StartFrom' mover requires 'y' coordinates option");
		if ( ! child_tag->hasOption("z") ) utility_exit_with_message("'StartFrom' mover requires 'z' coordinates option");

		core::Vector v(
				child_tag->getOption<core::Real>("x"),
				child_tag->getOption<core::Real>("y"),
				child_tag->getOption<core::Real>("z")
		);

		starting_points_.push_back(v);
	}
}

void StartFrom::apply(core::pose::Pose & pose){
	assert(!starting_points_.empty());
	int const starting_point_index= numeric::random::RG.random_range(1, starting_points_.size());
	core::Vector desired_centroid = starting_points_[starting_point_index];

	core::Size jump_id = core::pose::get_jump_id_from_chain(chain_, pose);
	move_ligand_to_desired_centroid(jump_id, desired_centroid, pose);
}

void
move_ligand_to_desired_centroid(
		core::Size const jump_id,
		core::Vector const desired_centroid,
		core::pose::Pose & pose
){
	core::Vector const ligand_centroid = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
	core::Vector const trans_vec = desired_centroid - ligand_centroid;
	core::Real const trans_len = trans_vec.length();
	if (trans_len > 1e-3) { // otherwise we get NaNs
		protocols::moves::RigidBodyTransMover mover(pose, jump_id);
		mover.step_size(trans_len);
		mover.trans_axis(trans_vec);
		mover.apply(pose);
	}
}


} //namespace ligand_docking
} //namespace protocols
