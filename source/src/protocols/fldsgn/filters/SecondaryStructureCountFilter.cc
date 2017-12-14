// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/SecondaryStructureCountFilter.cc
/// @brief filter structures by number of secondary structures
/// @details
/// @author Lei Shi ( shilei@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/SecondaryStructureCountFilter.hh>
#include <protocols/fldsgn/filters/SecondaryStructureCountFilterCreator.hh>

// Project Headers
#include <protocols/parser/BluePrint.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <basic/Tracer.hh>

// Core Headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>

// Parser headers
#include <protocols/filters/Filter.hh>
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/scoring/dssp/Dssp.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


//// C++ headers
static basic::Tracer tr( "protocols.fldsgn.filters.SecondaryStructureCountFilter" );

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
SecondaryStructureCountFilter::SecondaryStructureCountFilter( ):
	Filter( "SecondaryStructureCount" ),
	selector_( new core::select::residue_selector::TrueResidueSelector )
{}

// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool SecondaryStructureCountFilter::apply( core::pose::Pose const & pose ) const
{
	compute( pose );

	bool helix_filter=false;
	bool sheet_filter=false;
	bool loop_filter=false;
	bool helix_sheet_filter=false;

	if ( filter_helix_ ) {
		if  ( num_helix_pose_ >= num_helix_ ) {
			helix_filter=true;
		}
	} else {
		helix_filter=true;
	}

	if ( filter_sheet_ ) {
		if  ( num_sheet_pose_ >= num_sheet_ ) {
			sheet_filter=true;
		}
	} else {
		sheet_filter=true;
	}

	if ( filter_helix_sheet_ ) {
		if  ( num_helix_pose_ + num_sheet_pose_ >= num_helix_sheet_ ) {
			helix_sheet_filter=true;
		}
	} else {
		helix_sheet_filter=true;
	}

	if ( filter_loop_ ) {
		if  ( num_loop_pose_ >= num_loop_ ) {
			loop_filter=true;
		}
	} else {
		loop_filter=true;
	}

	if ( helix_filter && sheet_filter && helix_sheet_filter && loop_filter ) {
		return true;
	} else {
		return false;
	}
} // apply_filter

/// @brief parse xml
void
SecondaryStructureCountFilter::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	using protocols::parser::BluePrint;
	num_helix_ = tag->getOption<core::Size>( "num_helix", 0 );
	num_sheet_ = tag->getOption<core::Size>( "num_sheet", 0 );
	num_loop_ = tag->getOption<core::Size>( "num_loop", 0 );
	min_helix_length_ = tag->getOption<core::Size>( "min_helix_length", 4 );
	min_sheet_length_ = tag->getOption<core::Size>( "min_sheet_length", 3 );
	min_loop_length_ = tag->getOption<core::Size>( "min_loop_length", 0 );
	max_helix_length_ = tag->getOption<core::Size>( "max_helix_length", 9999 );
	max_sheet_length_ = tag->getOption<core::Size>( "max_sheet_length", 9999 );
	max_loop_length_ = tag->getOption<core::Size>( "max_loop_length", 9999 );
	num_helix_sheet_ = tag->getOption<core::Size>( "num_helix_sheet", 0);
	filter_helix_ = tag->getOption<bool>( "filter_helix", false );
	filter_sheet_ = tag->getOption<bool>( "filter_sheet", false );
	filter_loop_ = tag->getOption<bool>( "filter_loop", false );
	filter_helix_sheet_ = tag->getOption<bool>( "filter_helix_sheet", true );
	min_element_resis_ = tag->getOption<core::Size>( "min_element_resis", 1 );
	return_total_ = tag->getOption<bool>("return_total", false );

	std::string const selector_name = tag->getOption< std::string >( "residue_selector", "" );
	if ( !selector_name.empty() ) selector_ = core::select::residue_selector::get_residue_selector( selector_name, data );

	if ( filter_helix_ ) {
		tr << "filter on "<< num_helix_ << " helix with length: "<<min_helix_length_<<"-"<< max_helix_length_ << std::endl;
		runtime_assert( num_helix_ > 0 );
	}

	if ( filter_sheet_ ) {
		tr << "filter on "<< num_sheet_ << " sheet with length: "<<min_sheet_length_<<"-"<< max_sheet_length_ << std::endl;
		runtime_assert( num_sheet_ > 0 );
	}

	if ( filter_loop_ ) {
		tr << "filter on "<< num_loop_ << " loop with length: "<<min_loop_length_<<"-"<< max_loop_length_ << std::endl;
		runtime_assert( num_loop_ > 0 );
	}

	if ( filter_helix_sheet_ ) {
		tr << "filter on Sum of"<< num_helix_sheet_ << " helix with length: "<<min_helix_length_<<"-"<< max_helix_length_ << " AND sheet with length: "<<min_sheet_length_<<"-"<< max_sheet_length_ << std::endl;
		runtime_assert( num_helix_sheet_ > 0 );
	}

}

core::Size SecondaryStructureCountFilter::compute( core::pose::Pose const & pose ) const {
	core::pose::Pose pose_copy=pose;
	core::scoring::dssp::Dssp dssp( pose_copy );
	std::string dssp_ss=dssp.get_dssp_secstruct();
	tr << dssp_ss << std::endl;

	//create residue_selector subset
	using core::select::residue_selector::ResidueSubset;
	ResidueSubset const subset = selector_->apply( pose_copy );

	//Debugging residue_selector
	//ResidueSubset::const_iterator subs=subset.begin();
	//while ( subs< subset.end() ) {
	// tr << *subs << std::endl;
	// ++subs;
	//}

	num_helix_pose_=0;
	num_sheet_pose_=0;
	num_loop_pose_=0;

	core::Size tmp_count=0;
	core::Size resi_num=1;
	core::Size num_encounter=1;
	tr.Debug << "Only elements with >= " << min_element_resis_ << " residues will be counted." << std::endl;

	std::string::const_iterator iter=dssp_ss.begin();
	while ( iter< dssp_ss.end() ) {
		tr.Debug << "Looking at residue: " << resi_num << " with SS: " << *iter << std::endl;

		//if H
		if ( *iter=='H' ) {
			tmp_count=0;
			num_encounter=0;
			while ( *iter=='H' && iter< dssp_ss.end() ) {
				++tmp_count;
				++iter;
				++resi_num;
				if ( subset[ resi_num ] ) {
					++num_encounter;
				}
			}

			if ( tmp_count >= min_helix_length_ && tmp_count <= max_helix_length_ ) {
				if ( num_encounter >= min_element_resis_ ) {
					tr.Debug << "Terminal residue " << ( resi_num - 1 ) << " accepted as H SS element, encountered " << num_encounter << " times." << std::endl;
					num_helix_pose_+=1;
				}
			}

		} else if ( *iter=='E' ) { //end H
			//if E
			tmp_count=0;
			num_encounter=0;
			while ( *iter=='E' && iter< dssp_ss.end() ) {
				++tmp_count;
				++iter;
				++resi_num;
				if ( subset[ resi_num ] ) {
					++num_encounter;
				}
			}

			if ( tmp_count >= min_sheet_length_ && tmp_count <= max_sheet_length_ ) {
				if ( num_encounter >= min_element_resis_ ) {
					tr.Debug << "Terminal residue " << ( resi_num - 1 ) << " accepted as E SS element, encountered " << num_encounter << " times." << std::endl;
					num_sheet_pose_+=1;
				}
			}
		} else { //end E
			//if L
			tmp_count=0;
			num_encounter=0;
			while ( *iter=='L' && iter< dssp_ss.end() ) {
				++tmp_count;
				++iter;
				++resi_num;
				if ( subset[ resi_num ] ) {
					++num_encounter;
				}
			}

			if ( tmp_count >= min_loop_length_ && tmp_count <= max_loop_length_ ) {
				if ( num_encounter >= min_element_resis_ ) {
					tr.Debug << "Terminal residue " << ( resi_num - 1 ) << " accepted as L SS element, encountered " << num_encounter << " times." << std::endl;
					num_loop_pose_+=1;
				}
			}
			//increase iterator
			//++iter;
			//++resi_num;
		} //end L
	} //while iterator

	tr << " Pose has " << num_helix_pose_ << " helix, " << num_sheet_pose_  << " sheet, " << num_loop_pose_ << " loop (filtered elements), according to dssp_reduced definition." <<  std::endl;

	if ( return_total_ ) {
		tr << "Returning TOTAL number of filtered SS elements into score file." << std::endl;
		core::Size total=0;
		if ( filter_helix_ ) { total+=num_helix_pose_; }
		if ( filter_sheet_ ) { total+=num_sheet_pose_; }
		if ( filter_loop_ ) { total+=num_loop_pose_; }
		if ( filter_helix_sheet_ ) { total+=num_helix_pose_; total+=num_sheet_pose_; }
		return total;
	} else { return 0; }
}

core::Real SecondaryStructureCountFilter::report_sm( core::pose::Pose const & pose ) const {
	return compute( pose );
	//return 0;
}

void SecondaryStructureCountFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	compute( pose );
	out << " Pose has " << num_helix_pose_ << " helix " << num_sheet_pose_  << " sheet " << num_loop_pose_ << " loop (filtered elements) according to dssp_reduced definition." <<  std::endl;
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP SecondaryStructureCountFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SecondaryStructureCountFilter ); }

// XRW TEMP std::string
// XRW TEMP SecondaryStructureCountFilterCreator::keyname() const { return "SecondaryStructureCount"; }

std::string SecondaryStructureCountFilter::name() const {
	return class_name();
}

std::string SecondaryStructureCountFilter::class_name() {
	return "SecondaryStructureCount";
}

void SecondaryStructureCountFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "num_helix", xsct_non_negative_integer, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "num_sheet", xsct_non_negative_integer, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "num_helix_sheet", xsct_non_negative_integer, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "num_loop", xsct_non_negative_integer, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "min_helix_length", xsct_non_negative_integer, "XRW TO DO", "4" )
		+ XMLSchemaAttribute::attribute_w_default( "min_sheet_length", xsct_non_negative_integer, "XRW TO DO", "3" )
		+ XMLSchemaAttribute::attribute_w_default( "min_loop_length", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "max_helix_length", xsct_non_negative_integer, "XRW TO DO", "9999" )
		+ XMLSchemaAttribute::attribute_w_default( "max_sheet_length", xsct_non_negative_integer, "XRW TO DO", "9999" )
		+ XMLSchemaAttribute::attribute_w_default( "max_loop_length", xsct_non_negative_integer, "XRW TO DO", "9999" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_helix", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_sheet", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_loop", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_helix_sheet", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "return_total", xsct_rosetta_bool, "Returns total count to score file instead of 0.", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "min_element_resis", xsct_non_negative_integer, "Minimum number of residues on an element for it to be counted.", "1" )
		+ XMLSchemaAttribute( "residue_selector" , xs_string , "Explicitly set which SS elements to count using residue_selectors." );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string SecondaryStructureCountFilterCreator::keyname() const {
	return SecondaryStructureCountFilter::class_name();
}

protocols::filters::FilterOP
SecondaryStructureCountFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SecondaryStructureCountFilter );
}

void SecondaryStructureCountFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SecondaryStructureCountFilter::provide_xml_schema( xsd );
}



} // filters
} // fldsgn
} // protocols
