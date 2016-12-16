// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

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

// XRW TEMP std::string
// XRW TEMP KeepRegionMover::get_name() const {
// XRW TEMP  return "KeepRegionMover";
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP KeepRegionMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new KeepRegionMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP KeepRegionMoverCreator::keyname() const {
// XRW TEMP  return KeepRegionMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP KeepRegionMover::mover_name() {
// XRW TEMP  return "KeepRegionMover";
// XRW TEMP }

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
	tag->getOption<std::string>( "start_" );
	tag->getOption<std::string>( "end_" );

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

	if ( start_ == 1 && end_ == pose.size() ) {
		return;
	}

	PyAssert(start_ != 0, "Cannot keep region starting with 0 - make sure region is set for KeepRegionMover");
	PyAssert(end_ !=0, "Cannot keep region ending with 0 - make sure region is set for KeepRegionMover");
	PyAssert(end_ > start_, "Cannot keep region where end > start");
	PyAssert(end_ < pose.size(), "Cannot keep region where end is > pose size");




	core::pose::Pose temp_pose = protocols::grafting::return_region(pose, start_- nter_overhang_, end_ + cter_overhang_);
	pose = temp_pose;

}

std::string KeepRegionMover::get_name() const {
	return mover_name();
}

std::string KeepRegionMover::mover_name() {
	return "KeepRegionMover";
}

void KeepRegionMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		//+ XMLSchemaAttribute::required_attribute( "start_", xsct_non_negative_integer, "First residue of region" )
		//+ XMLSchemaAttribute::required_attribute( "end_", xsct_non_negative_integer, "Last residue of region" )
		+ XMLSchemaAttribute::attribute_w_default( "nter_overhang", xsct_non_negative_integer, "Number of residues N terminal to start to include", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "cter_overhang", xsct_non_negative_integer, "Number of residues C terminal to end to include", "0");

	core::pose::attributes_for_get_resnum( attlist, "start_" );
	core::pose::attributes_for_get_resnum( attlist, "end_" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Keeps a specified region of the current pose and deletes the rest", attlist );
}

std::string KeepRegionMoverCreator::keyname() const {
	return KeepRegionMover::mover_name();
}

protocols::moves::MoverOP
KeepRegionMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new KeepRegionMover );
}

void KeepRegionMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	KeepRegionMover::provide_xml_schema( xsd );
}


}
}
}
