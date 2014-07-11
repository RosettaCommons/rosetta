// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/simple_movers/InsertPoseIntoPoseMover.hh
/// @brief  
/// @author  Jared Adolf-Bryfogle

#include <protocols/grafting/simple_movers/InsertPoseIntoPoseMover.hh>
#include <protocols/grafting/simple_movers/InsertPoseIntoPoseMoverCreator.hh>

#include <core/pose/selection.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/grafting/util.hh>

#include <basic/datacache/DataCache.hh>
#include <utility/PyAssert.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace grafting {
namespace simple_movers {

InsertPoseIntoPoseMover::InsertPoseIntoPoseMover(bool copy_pdbinfo /*false*/):
protocols::moves::Mover("InsertPoseIntoPoseMover"),
		src_pose_(NULL),
		start_(0),
		end_(0),
		copy_pdbinfo_(copy_pdbinfo),
		tag_(NULL)
{
	
}
	
InsertPoseIntoPoseMover::InsertPoseIntoPoseMover(
	const core::pose::Pose& src_pose,
	core::Size res_start,
	core::Size res_end,
	bool copy_pdbinfo /*false*/) :
	protocols::moves::Mover("InsertPoseIntoPoseMover"),
	start_(res_start),
	end_(res_end),
	copy_pdbinfo_(copy_pdbinfo),
	tag_(NULL)
{
	src_pose_ = new core::pose::Pose(src_pose);
}
	
	
InsertPoseIntoPoseMover::InsertPoseIntoPoseMover(const InsertPoseIntoPoseMover& src) :
	protocols::moves::Mover(src),
	
	start_(src.start_),
	end_(src.end_),
	copy_pdbinfo_(src.copy_pdbinfo_),
	src_pose_(src.src_pose_),
	tag_(src.tag_)
{
	
}

InsertPoseIntoPoseMover::~InsertPoseIntoPoseMover(){}

void
InsertPoseIntoPoseMover::src_pose(const core::pose::Pose& src_pose){
	src_pose_ = new core::pose::Pose(src_pose);
}

void
InsertPoseIntoPoseMover::start(core::Size res_start){
	start_ = res_start;
}

void
InsertPoseIntoPoseMover::end(core::Size res_end){
	end_ = res_end;
}

core::Size
InsertPoseIntoPoseMover::start() const{
	return start_;
}

core::Size
InsertPoseIntoPoseMover::end() const {
	return end_;
}

std::string
InsertPoseIntoPoseMover::get_name() const {
	return "InsertPoseIntoPoseMover";
}

protocols::moves::MoverOP
InsertPoseIntoPoseMover::clone() const {
	return new InsertPoseIntoPoseMover(*this);
}

void
InsertPoseIntoPoseMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap& data_map,
	const Filters_map& ,
	const moves::Movers_map& ,
	const Pose& )
{
	tag_ = tag->clone();
	
	//Protect from unused option crash.
	protocols::rosetta_scripts::parse_bogus_res_tag(tag, "start_");
	protocols::rosetta_scripts::parse_bogus_res_tag(tag, "end_");
	
	src_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag, data_map, "spm_reference_name");
	copy_pdbinfo_ = tag->getOption<bool>("copy_pdbinfo", false);
	
}

protocols::moves::MoverOP
InsertPoseIntoPoseMoverCreator::create_mover() const {
	return new InsertPoseIntoPoseMover;
}

std::string
InsertPoseIntoPoseMoverCreator::keyname() const {
	return InsertPoseIntoPoseMoverCreator::mover_name();
}

std::string
InsertPoseIntoPoseMoverCreator::mover_name(){
	return "InsertPoseIntoPoseMover";
}

void
InsertPoseIntoPoseMover::apply(core::pose::Pose& pose) {
	
	if (tag_){
		start_ = core::pose::get_resnum(tag_, pose, "start_");
		end_ = core::pose::get_resnum(tag_, pose, "end_");
	}
	
	PyAssert(start_ != 0, "Cannot insert region starting with 0 - make sure region is set for InsertPoseIntoPoseMover");
	PyAssert(end_ !=0, "Cannot insert region ending with 0 - make sure region is set for InsertPoseIntoPoseMover");
	PyAssert(end_ > start_, "Cannot insert into a region where end > start");
	PyAssert(end_ <= pose.total_residue(), "Cannot insert a region where end is > pose.total_residues of the scaffold");
	
	pose = protocols::grafting::insert_pose_into_pose(pose, *src_pose_, start_, end_);
	
}


}
}
}


