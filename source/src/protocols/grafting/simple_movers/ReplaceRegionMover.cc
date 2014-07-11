// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/simple_movers/ReplaceRegionMover.hh
/// @brief  
/// @author  Jared Adolf-Bryfogle

#include <protocols/grafting/simple_movers/ReplaceRegionMover.hh>
#include <protocols/grafting/simple_movers/ReplaceRegionMoverCreator.hh>

#include <core/pose/selection.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>



#include <protocols/grafting/util.hh>
#include <utility/PyAssert.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataCache.hh>

namespace protocols {
namespace grafting {
namespace simple_movers {

ReplaceRegionMover::ReplaceRegionMover(bool copy_pdbinfo /*false*/):
protocols::moves::Mover("ReplaceRegionMover"),
		src_pose_(NULL),
		src_pose_start_(0),
		target_pose_start_(0),
		span_(0),
		copy_pdbinfo_(copy_pdbinfo),
		tag_(NULL)
{
	
}
	
ReplaceRegionMover::ReplaceRegionMover(
		core::pose::Pose const & src_pose,
		Size const src_pose_start, 
		Size const target_pose_start,
		Size const span,
		bool copy_pdbinfo /*false*/) :
	protocols::moves::Mover("ReplaceRegionMover"),
	src_pose_start_(src_pose_start),
	target_pose_start_(target_pose_start),
	span_(span),
	copy_pdbinfo_(copy_pdbinfo),
	tag_(NULL)
{
	src_pose_ = new core::pose::Pose(src_pose);
}
	
	
ReplaceRegionMover::ReplaceRegionMover(const ReplaceRegionMover& src) :
	protocols::moves::Mover(src),
	src_pose_start_(src.src_pose_start_),
	target_pose_start_(src.target_pose_start_),
	span_(src.span_),
	copy_pdbinfo_(src.copy_pdbinfo_),
	src_pose_(src.src_pose_),
	tag_(src.tag_)
{
	
}

ReplaceRegionMover::~ReplaceRegionMover(){}

std::string
ReplaceRegionMover::get_name() const {
	return "ReplaceRegionMover";
}


void
ReplaceRegionMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap& data_map,
	const Filters_map& ,
	const moves::Movers_map& ,
	const Pose& )
{

	src_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag, data_map, "spm_reference_name");
	tag_ = tag->clone();
	span_ = tag->getOption<core::Size>("span");
	copy_pdbinfo_ = tag->getOption<bool>("copy_pdbinfo", false);
	
	//Protect from unused option crash.
	protocols::rosetta_scripts::parse_bogus_res_tag(tag, "src_pose_start_");
	protocols::rosetta_scripts::parse_bogus_res_tag(tag, "target_pose_start_");
	
	
}

protocols::moves::MoverOP
ReplaceRegionMoverCreator::create_mover() const {
	return new ReplaceRegionMover;
}

std::string
ReplaceRegionMoverCreator::keyname() const {
	return ReplaceRegionMoverCreator::mover_name();
}

std::string
ReplaceRegionMoverCreator::mover_name(){
	return "ReplaceRegionMover";
}

void ReplaceRegionMover::src_pose_start(core::Size start){ src_pose_start_ = start; }
core::Size ReplaceRegionMover::src_pose_start() const{ return src_pose_start_; }

void ReplaceRegionMover::target_pose_start(core::Size start) { target_pose_start_ = start; }
core::Size ReplaceRegionMover::target_pose_start() const { return target_pose_start_;}

void 
ReplaceRegionMover::src_pose(const core::pose::Pose& src_pose){
	src_pose_ = new core::pose::Pose(src_pose);
}

protocols::moves::MoverOP
ReplaceRegionMover::clone() const {
	return new ReplaceRegionMover(*this);
}

protocols::moves::MoverOP
ReplaceRegionMover::fresh_instance() const {
	return new ReplaceRegionMover;
}

void
ReplaceRegionMover::apply(core::pose::Pose& pose) {
	
	if (tag_){
		src_pose_start_ = core::pose::get_resnum(tag_, *src_pose_, "src_pose_start_");
		target_pose_start_= core::pose::get_resnum(tag_, pose, "target_pose_start_");
	}
	
	PyAssert(src_pose_start_ != 0, "Cannot copy from a region starting with 0 - make sure region is set for ReplaceRegionMover");
	PyAssert(target_pose_start_ != 0, "Cannot copy to a region starting with 0 - make sure region is set for ReplaceRegionMover");
	PyAssert(span_ <= src_pose_->total_residue(), "Cannot delete region ending with 0 - make sure region is set for ReplaceRegionMover");
	PyAssert(src_pose_start_ + span_ <= pose.total_residue(), "Not enough residues in pose to copy all of the span of the source pose");
	
	
	pose = protocols::grafting::replace_region (*src_pose_, pose, src_pose_start_, target_pose_start_, span_);
	 
}


}
}
}


