// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/simple_movers/DeleteRegionMover.hh
/// @brief  
/// @author  Jared Adolf-Bryfogle


#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMoverCreator.hh>

#include <core/pose/selection.hh>

#include <protocols/grafting/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/py/PyAssert.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataCache.hh>
#include <utility>

namespace protocols {
namespace grafting {
namespace simple_movers {
			
DeleteRegionMover::DeleteRegionMover():
	protocols::moves::Mover("DeleteRegionMover"),
	start_(0),
	end_(0),
	nter_overhang_(0),
	cter_overhang_(0),
	tag_(/* NULL */)
{
	
}

DeleteRegionMover::DeleteRegionMover(core::Size res_start, core::Size res_end):
	protocols::moves::Mover("DeleteRegionMover"),
	start_(res_start),
	end_(res_end),
	nter_overhang_(0),
	cter_overhang_(0),
	tag_(/* NULL */)
{
	
}
		
DeleteRegionMover::DeleteRegionMover(const DeleteRegionMover& src):
	protocols::moves::Mover(src),
	start_(src.start_),
	end_(src.end_),
	nter_overhang_(src.nter_overhang_),
	cter_overhang_(src.cter_overhang_),
	tag_(src.tag_)
{
	
}

DeleteRegionMover::~DeleteRegionMover(){}

std::string
DeleteRegionMover::get_name() const{
	return "DeleteRegionMover";
}

void
DeleteRegionMover::region(core::Size res_start, core::Size res_end){
	start_ = res_start;
	end_ = res_end;
}

std::pair<core::Size, core::Size>
DeleteRegionMover::region() const{
	return std::make_pair(start_, end_);
}

void
DeleteRegionMover::start(core::Size res_start){
	start_ = res_start;
}

void
DeleteRegionMover::end(core::Size res_end){
	end_ = res_end;
}

core::Size
DeleteRegionMover::start() const{
	return start_;
}

core::Size
DeleteRegionMover::end() const{
	return end_;
}

protocols::moves::MoverOP
DeleteRegionMover::clone() const{
	return protocols::moves::MoverOP( new DeleteRegionMover(*this) );
}

protocols::moves::MoverOP
DeleteRegionMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new DeleteRegionMover );
}

protocols::moves::MoverOP
DeleteRegionMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new DeleteRegionMover );
}

std::string
DeleteRegionMoverCreator::keyname() const {
	return DeleteRegionMoverCreator::mover_name();
}

std::string
DeleteRegionMoverCreator::mover_name(){
	return "DeleteRegionMover";
}

void
DeleteRegionMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap& ,
	const Filters_map& ,
	const moves::Movers_map& ,
	const Pose& )
{
	
	tag_ = tag->clone();
	//std::cout << "Parsing tag"<<std::endl;
	protocols::rosetta_scripts::parse_bogus_res_tag(tag_, "start_");
	protocols::rosetta_scripts::parse_bogus_res_tag(tag_, "end_");
	
	nter_overhang_ = tag_->getOption<core::Size>("nter_overhang", 0);
	cter_overhang_ = tag_->getOption<core::Size>("cter_overhang", 0);
	//std::cout << " N "<<nter_overhang_<< " C " << cter_overhang_<<std::endl;
	
	protocols::rosetta_scripts::parse_bogus_res_tag(tag, "start_");
	protocols::rosetta_scripts::parse_bogus_res_tag(tag, "end_");
	
	nter_overhang_ = tag->getOption<core::Size>("nter_overhang", 0);
	cter_overhang_ = tag->getOption<core::Size>("cter_overhang", 0);
	//std::cout << " N "<<nter_overhang_<< " C " << cter_overhang_<<std::endl;
}

void
DeleteRegionMover::apply(core::pose::Pose& pose){
	
	if (tag_){
		start_ = core::pose::get_resnum(tag_, pose, "start_");
		end_ = core::pose::get_resnum(tag_, pose, "end_");
	}
	
	//std::cout <<"Start: "<<start_ << " End: " << end_ << std::endl;
	PyAssert(start_ != 0, "Cannot delete region starting with 0 - make sure region is set for DeleteRegionMover");
	PyAssert(end_ !=0, "Cannot delete region ending with 0 - make sure region is set for DeleteRegionMover");
	PyAssert(end_ > start_, "Cannot delete region where end > start");
	PyAssert(end_ <= pose.total_residue(), "Cannot delete region where end is > pose total_residues");
	
	protocols::grafting::delete_region(pose, start_ -  nter_overhang_, end_ + cter_overhang_);
	
}

}
}
}
