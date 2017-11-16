// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/WriteFiltersToPose.cc
/// @brief Writes all filter results to the pose file.
/// @author Cody Krivacic (krivacic@berkeley.edu)


//Unit Headers
#include <protocols/simple_moves/WriteFiltersToPose.hh>
#include <protocols/simple_filters/ReportFilterCreator.hh>
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <protocols/simple_moves/WriteFiltersToPoseCreator.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>
#include <protocols/moves/mover_schemas.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
#include <map>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
namespace protocols {
namespace simple_moves {

using namespace core;
using namespace core::scoring;

//using namespace ObjexxFCL::format;

//static basic::Tracer TR( "protocols.simple_filters.ReportFilter" );

/// @brief Writes all filter results to the pose file.
/// @details This mover goes through all filters in the RosettaScripts file and writes their
/// user-defined name, (optionally) their filter type, and their numerical result to the
/// output .pdb file. The include_type option writes the filter type.
/// The user can also specify a prefix to append before the user-defined
/// filter name; for instance, let's say there is a PackStat filter, as well as several movers like
/// fastdesign and docking. The "prefix" option allows the user to have several versions of
/// WriteFiltersToPose which take the filter results from different steps in the process.
/// If the RosettaScripts XML file has these instances of WriteFiltersToPose:
///      <WriteFiltersToPose name="writer1" prefix="post_fastdesign_"/>
///      <WriteFiltersToPose name="writer2" prefix="post_docking_" include_type="true"/>
/// If the only filter is a PackStat filter named "packing",
/// the output pdb will have two lines at the end that look something like this:
///      post_docking_packing PackStat 0.6794
///      post_fastdesign_packing 0.6751
/// If there are multiple filters, each instance of this mover will write all
/// filters to the pose.
///
/// Note: If you only want to write a specific filter result to your pose,
/// check out FilterReportAsPoseExtraScoresMover

WriteFiltersToPose::WriteFiltersToPose(){
	//nothing
}


void WriteFiltersToPose::apply( core::pose::Pose &pose ) {
	std::string filtername;
	for ( auto const &f : filters_ ) {
		protocols::filters::FilterOP filter = f.second;
		if ( include_type_ ) {
			filtername = prefix_ + filter->get_user_defined_name() + " " + filter->get_type();
		} else {
			filtername = prefix_ + filter->get_user_defined_name();
		}
		if ( (filter->get_type()!="FalseFilter") && (filter->get_type()!="TrueFilter") ) {
			core::Real score = filter->report_sm( pose );
			core::pose::setPoseExtraScore( pose, filtername, score );
		}
	}
}

void WriteFiltersToPose::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &f,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	filters_=f;
	prefix_ = tag->getOption< std::string >( "prefix" , "" );
	include_type_= tag->getOption< bool > ( "include_type" ,false );
}

std::string WriteFiltersToPose::get_name() const{
	return "WriteFiltersToPose";
}

protocols::moves::MoverOP
WriteFiltersToPoseCreator::create_mover() const {
	return protocols::moves::MoverOP( new WriteFiltersToPose );
}
std::string WriteFiltersToPoseCreator::keyname() const {
	return WriteFiltersToPose::mover_name();
}

std::string WriteFiltersToPose::mover_name() {
	return "WriteFiltersToPose";
}

void WriteFiltersToPose::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//attlist
	attlist
		+ XMLSchemaAttribute::attribute_w_default(
		"prefix", xs_string,
		"Prefix to append to report output, so you can have multiple WriteFiltersToPose instances at different points in your RosettaScripts and not have them overwrite each other.","")
		+ XMLSchemaAttribute::attribute_w_default(
		"include_type",xsct_rosetta_bool,"Include the filter type in your output (whether you want this depends on how you are parsing the pose files).","false"
	);

	//+ XMLSchemaAttribute( "scorefxn", xs_string, "Score function to use when evaluating best amino acids at each position" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Takes results from any filters and writes them to the end of the pose file.", attlist );
}

void WriteFiltersToPoseCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	WriteFiltersToPose::provide_xml_schema( xsd );
}

}
}
