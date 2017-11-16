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
#include <utility/string_util.hh>
#include <utility/py/PyAssert.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.grafting.simple_movers.KeepRegionMover" );

namespace protocols {
namespace grafting {
namespace simple_movers {


KeepRegionMover::KeepRegionMover() :
	protocols::moves::Mover("KeepRegionMover"),
	nter_overhang_(0),
	cter_overhang_(0)
{

}

KeepRegionMover::KeepRegionMover(std::string const & res_start, std::string const & res_end):
	start_(res_start),
	end_(res_end),
	nter_overhang_(0),
	cter_overhang_(0)
{

}

KeepRegionMover::~KeepRegionMover() {}

KeepRegionMover::KeepRegionMover(KeepRegionMover const & ) = default;

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
	start_ = core::pose::get_resnum_string(tag, "start_");
	end_ = core::pose::get_resnum_string(tag, "end_");

	nter_overhang_ = tag->getOption<core::Size>("nter_overhang", 0);
	cter_overhang_ = tag->getOption<core::Size>("cter_overhang", 0);
}

void
KeepRegionMover::region(std::string const & res_start, std::string const & res_end){
	start_ = res_start;
	end_ = res_end;
}

std::pair<std::string, std::string>
KeepRegionMover::region() const{
	return std::make_pair(start_, end_);
}

void
KeepRegionMover::start(std::string const & res_start){
	start_ = res_start;
}

void
KeepRegionMover::end(std::string const & res_end){
	end_ = res_end;
}

void
KeepRegionMover::start(core::Size res_start){
	start_ = utility::to_string(res_start);
}

void
KeepRegionMover::end(core::Size res_end){
	end_ = utility::to_string(res_end);
}

std::string const &
KeepRegionMover::start() const{
	return start_;
}

std::string const &
KeepRegionMover::end() const{
	return end_;
}

void
KeepRegionMover::apply(core::pose::Pose& pose) {

	core::Size start = core::pose::parse_resnum(start_, pose);
	core::Size end = core::pose::parse_resnum(end_, pose);

	if ( start == 1 && end == pose.size() ) {
		return;
	}

	PyAssert(start != 0, "Cannot keep region starting with 0 - make sure region is set for KeepRegionMover");
	PyAssert(end !=0, "Cannot keep region ending with 0 - make sure region is set for KeepRegionMover");
	PyAssert(end > start, "Cannot keep region where end > start");
	PyAssert(end < pose.size(), "Cannot keep region where end is > pose size");


	core::pose::Pose temp_pose = protocols::grafting::return_region(pose, start - nter_overhang_, end + cter_overhang_);
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

	core::pose::attributes_for_get_resnum_string( attlist, "start_" );
	core::pose::attributes_for_get_resnum_string( attlist, "end_" );

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
