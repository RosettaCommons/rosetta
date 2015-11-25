// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/grafting/simple_movers/KeepRegionMover.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/moves/Mover.hh>
#include <protocols/grafting/simple_movers/KeepRegionMover.hh>
#include <protocols/grafting/simple_movers/KeepRegionMoverCreator.hh>

#include <core/pose/selection.hh>
#include <protocols/grafting/util.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <utility/tag/Tag.hh>
#include <utility/py/PyAssert.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.grafting.simple_movers.KeepRegionMover" );

namespace protocols {
namespace grafting {
namespace simple_movers {


KeepRegionMover::KeepRegionMover() :
	protocols::moves::Mover("KeepRegionMover"),
	start_(0),
	end_(0),
	nter_overhang_(0),
	cter_overhang_(0),
	tag_(/* NULL */)
{

}

KeepRegionMover::KeepRegionMover(core::Size res_start, core::Size res_end):
	start_(res_start),
	end_(res_end),
	nter_overhang_(0),
	cter_overhang_(0),
	tag_(/* NULL */)
{

}

KeepRegionMover::~KeepRegionMover() {}

KeepRegionMover::KeepRegionMover(KeepRegionMover const & src) :
	protocols::moves::Mover(src),
	start_(src.start_),
	end_(src.end_),
	nter_overhang_(src.nter_overhang_),
	cter_overhang_(src.cter_overhang_),
	tag_(src.tag_)
{

}

protocols::moves::MoverOP
KeepRegionMover::clone() const {
	return protocols::moves::MoverOP( new KeepRegionMover(*this) );
}

protocols::moves::MoverOP
KeepRegionMover::fresh_instance() const {
	return protocols::moves::MoverOP( new KeepRegionMover );
}

std::string
KeepRegionMover::get_name() const {
	return "KeepRegionMover";
}

protocols::moves::MoverOP
KeepRegionMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new KeepRegionMover );
}

std::string
KeepRegionMoverCreator::keyname() const {
	return KeepRegionMoverCreator::mover_name();
}

std::string
KeepRegionMoverCreator::mover_name() {
	return "KeepRegionMover";
}

void
KeepRegionMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap&,
	const Filters_map&,
	const protocols::moves::Movers_map&,
	const Pose& )
{
	tag_ = tag->clone();
	//Protect from unused option crash.
	protocols::rosetta_scripts::parse_bogus_res_tag(tag, "start_");
	protocols::rosetta_scripts::parse_bogus_res_tag(tag, "end_");

	nter_overhang_ = tag->getOption<core::Size>("nter_overhang", 0);
	cter_overhang_ = tag->getOption<core::Size>("cter_overhang", 0);
}

void
KeepRegionMover::region(core::Size res_start, core::Size res_end){
	start_ = res_start;
	end_ = res_end;
}

std::pair<core::Size, core::Size>
KeepRegionMover::region() const{
	return std::make_pair(start_, end_);
}

void
KeepRegionMover::start(core::Size res_start){
	start_ = res_start;
}

void
KeepRegionMover::end(core::Size res_end){
	end_ = res_end;
}

core::Size
KeepRegionMover::start() const{
	return start_;
}

core::Size
KeepRegionMover::end() const{
	return end_;
}

void
KeepRegionMover::apply(core::pose::Pose& pose) {

	if ( tag_ ) {
		start_ = core::pose::get_resnum(tag_, pose, "start_");
		end_ = core::pose::get_resnum(tag_, pose, "end_");
	}

	if ( start_ == 1 && end_ == pose.total_residue() ) {
		return;
	}

	PyAssert(start_ != 0, "Cannot keep region starting with 0 - make sure region is set for KeepRegionMover");
	PyAssert(end_ !=0, "Cannot keep region ending with 0 - make sure region is set for KeepRegionMover");
	PyAssert(end_ > start_, "Cannot keep region where end > start");
	PyAssert(end_ < pose.total_residue(), "Cannot keep region where end is > pose total_residues");




	core::pose::Pose temp_pose = protocols::grafting::return_region(pose, start_- nter_overhang_, end_ + cter_overhang_);
	pose = temp_pose;

}

}
}
}
