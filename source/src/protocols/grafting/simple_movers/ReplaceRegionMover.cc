// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/grafting/simple_movers/ReplaceRegionMover.hh
/// @brief
/// @author  Jared Adolf-Bryfogle

#include <protocols/grafting/simple_movers/ReplaceRegionMover.hh>
#include <protocols/grafting/simple_movers/ReplaceRegionMoverCreator.hh>

#include <core/pose/selection.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>


#include <protocols/grafting/util.hh>
#include <utility>
#include <utility/py/PyAssert.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataCache.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace grafting {
namespace simple_movers {

ReplaceRegionMover::ReplaceRegionMover(bool copy_pdbinfo /*false*/):
	protocols::moves::Mover("ReplaceRegionMover"),
	span_(0),
	copy_pdbinfo_(copy_pdbinfo),
	src_pose_(/* NULL */)
{

}

ReplaceRegionMover::ReplaceRegionMover(
	core::pose::Pose const & src_pose,
	std::string const & src_pose_start,
	std::string const & target_pose_start,
	Size const span,
	bool copy_pdbinfo /*false*/
) :
	protocols::moves::Mover("ReplaceRegionMover"),
	src_pose_start_(src_pose_start),
	target_pose_start_(target_pose_start),
	span_(span),
	copy_pdbinfo_(copy_pdbinfo)
{
	src_pose_ = core::pose::PoseOP( new core::pose::Pose(src_pose) );
}


ReplaceRegionMover::ReplaceRegionMover(const ReplaceRegionMover& ) = default;


ReplaceRegionMover::~ReplaceRegionMover()= default;

// XRW TEMP std::string
// XRW TEMP ReplaceRegionMover::get_name() const {
// XRW TEMP  return "ReplaceRegionMover";
// XRW TEMP }


void
ReplaceRegionMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap& data_map,
	const Filters_map& ,
	const moves::Movers_map& ,
	const Pose& )
{

	src_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag, data_map, "spm_reference_name");
	span_ = tag->getOption<core::Size>("span");
	copy_pdbinfo_ = tag->getOption<bool>("copy_pdbinfo", false);

	src_pose_start_ = tag->getOption<std::string>("src_pose_start");
	target_pose_start_= tag->getOption<std::string>("target_pose_start");
}

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ReplaceRegionMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ReplaceRegionMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ReplaceRegionMoverCreator::keyname() const {
// XRW TEMP  return ReplaceRegionMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ReplaceRegionMover::mover_name(){
// XRW TEMP  return "ReplaceRegionMover";
// XRW TEMP }

void ReplaceRegionMover::src_pose_start(std::string const & start){ src_pose_start_ = start; }
std::string const & ReplaceRegionMover::src_pose_start() const{ return src_pose_start_; }

void ReplaceRegionMover::target_pose_start(std::string const & start) { target_pose_start_ = start; }
std::string const &  ReplaceRegionMover::target_pose_start() const { return target_pose_start_;}

void
ReplaceRegionMover::src_pose(const core::pose::Pose& src_pose){
	src_pose_ = core::pose::PoseCOP( new core::pose::Pose(src_pose) );
}

protocols::moves::MoverOP
ReplaceRegionMover::clone() const {
	return protocols::moves::MoverOP( new ReplaceRegionMover(*this) );
}

protocols::moves::MoverOP
ReplaceRegionMover::fresh_instance() const {
	return protocols::moves::MoverOP( new ReplaceRegionMover );
}

void
ReplaceRegionMover::apply(core::pose::Pose& pose) {

	core::Size src_pose_start = core::pose::parse_resnum(src_pose_start_, *src_pose_);
	core::Size target_pose_start = core::pose::parse_resnum(target_pose_start_, pose);

	//std::cout << src_pose_start << std::endl;
	//std::cout << target_pose_start << std::endl;
	PyAssert(src_pose_start != 0, "Cannot copy from a region starting with 0 - make sure region is set for ReplaceRegionMover");
	PyAssert(target_pose_start != 0, "Cannot copy to a region starting with 0 - make sure region is set for ReplaceRegionMover");
	PyAssert(span_ <= src_pose_->size(), "Cannot delete region ending with 0 - make sure region is set for ReplaceRegionMover");
	PyAssert(src_pose_start + span_ <= pose.size(), "Not enough residues in pose to copy all of the span of the source pose");


	pose = protocols::grafting::replace_region (*src_pose_, pose, src_pose_start, target_pose_start, span_);

}

std::string ReplaceRegionMover::get_name() const {
	return mover_name();
}

std::string ReplaceRegionMover::mover_name() {
	return "ReplaceRegionMover";
}

void ReplaceRegionMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "span", xsct_non_negative_integer, "Length of the region to replace" )
		+ XMLSchemaAttribute::attribute_w_default( "copy_pdbinfo", xsct_rosetta_bool, "Copy PDBInfo to the new pose?", "false" )
		+ XMLSchemaAttribute::required_attribute( "src_pose_start", xs_string, "First residue to copy from the source pose" )
		+ XMLSchemaAttribute::required_attribute( "target_pose_start", xs_string, "First residue to replace in the target pose" );
	rosetta_scripts::attributes_for_saved_reference_pose( attlist, "spm_reference_name" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Replaces a region of length span in the target pose with a specified region of length span in the source pose.", attlist );
}

std::string ReplaceRegionMoverCreator::keyname() const {
	return ReplaceRegionMover::mover_name();
}

protocols::moves::MoverOP
ReplaceRegionMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ReplaceRegionMover );
}

void ReplaceRegionMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReplaceRegionMover::provide_xml_schema( xsd );
}



}
}
}
