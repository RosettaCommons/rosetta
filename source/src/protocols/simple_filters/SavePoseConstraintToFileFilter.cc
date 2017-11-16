// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/SavePoseConstraintToFileFilter.cc
/// @brief Filter for looking at specific atom distances
/// @author Lei Shi(shilei@uw.edu)

#include <protocols/simple_filters/SavePoseConstraintToFileFilter.hh>
#include <protocols/simple_filters/SavePoseConstraintToFileFilterCreator.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

#include <fstream>
#include <iostream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.filters.SavePoseConstraintToFileFilter" );

/// @brief default ctor
SavePoseConstraintToFileFilter::SavePoseConstraintToFileFilter() :
	parent( "SavePoseConstraintToFile" )
{}

/// @return Print pose information and return true
bool SavePoseConstraintToFileFilter::apply(core::pose::Pose const & pose ) const
{
	compute( pose );
	return true;
}

core::Real
SavePoseConstraintToFileFilter::compute( core::pose::Pose const & pose ) const
{
	report(TR, pose);
	return 1.0;
}

/// @return Print pose information and return true
core::Real
SavePoseConstraintToFileFilter::report_sm( core::pose::Pose const & pose ) const
{
	compute( pose );
	return( 1 );
}

void SavePoseConstraintToFileFilter::report( std::ostream &, core::pose::Pose const & pose ) const
{
	std::ofstream outfile;
	outfile.open(filename_.c_str(),std::ios::out);
	if ( ! outfile.good() ) {
		utility_exit_with_message( "Unable to open file: " + filename_);
	}

	pose.constraint_set()->show_definition(outfile,pose);
}

void SavePoseConstraintToFileFilter::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/)
{
	overwrite_ = tag->getOption< bool >( "overwrite", false );

	if ( tag->hasOption("filename") ) {
		filename_ = tag->getOption<std::string>("filename");
		std::fstream file;
		file.open( filename_.c_str(), std::ios::out | std::ios::in );
		if ( file.is_open() ) {
			if ( overwrite_ ) {
				TR << filename_ << " exists and will overwrite" << std::endl;
				file.clear();
				file.close();
			} else {
				utility_exit_with_message(filename_+" exists, specify overwrite=1");
			}
		} //clear file contents and throw out warning.
	} else {
		utility_exit_with_message("Must specify a filename");
	}//require a filename to be specified
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP SavePoseConstraintToFileFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SavePoseConstraintToFileFilter ); }

// XRW TEMP std::string
// XRW TEMP SavePoseConstraintToFileFilterCreator::keyname() const { return "SavePoseConstraintToFile"; }

std::string SavePoseConstraintToFileFilter::name() const {
	return class_name();
}

std::string SavePoseConstraintToFileFilter::class_name() {
	return "SavePoseConstraintToFile";
}

void SavePoseConstraintToFileFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default("overwrite", xsct_rosetta_bool, "XRW TO DO", "false")
		+ XMLSchemaAttribute("filename", xs_string, "XRW TO DO");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string SavePoseConstraintToFileFilterCreator::keyname() const {
	return SavePoseConstraintToFileFilter::class_name();
}

protocols::filters::FilterOP
SavePoseConstraintToFileFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SavePoseConstraintToFileFilter );
}

void SavePoseConstraintToFileFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SavePoseConstraintToFileFilter::provide_xml_schema( xsd );
}



} // filters
} // protocols
