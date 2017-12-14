// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/filters/SecondaryStructureHasResidueFilter.cc
/// @brief filter structures by number of secondary structures
/// @details use this filter to insure that each secondary element has one or more of a given set of residues
///  motivated by observation that alpha,beta elements without at least 1 h'phobic residue often do not fold in simulation
///  reports fraction of 2' struct elements that contain N or more residues from given residue list (default VILMFYM)
/// @author Chris King ( chrisk1@uw.edu )

// Unit Headers
#include <protocols/fldsgn/filters/SecondaryStructureHasResidueFilter.hh>
#include <protocols/fldsgn/filters/SecondaryStructureHasResidueFilterCreator.hh>

// Project Headers
#include <protocols/parser/BluePrint.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/rosetta_scripts/util.hh>

// Utility headers
#include <basic/Tracer.hh>

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
static basic::Tracer tr("protocols.fldsgn.filters.SecondaryStructureHasResidueFilter");

namespace protocols {
namespace fldsgn {
namespace filters {

// @brief default constructor
SecondaryStructureHasResidueFilter::SecondaryStructureHasResidueFilter( ): Filter( "SecondaryStructureHasResidue" )
{}

// @brief returns true if the given pose passes the filter, false otherwise.
bool SecondaryStructureHasResidueFilter::apply( core::pose::Pose const & pose ) const
{

	if ( compute( pose ) >= threshold_ ) {
		return true;
	} else {
		return false;
	}
} // apply_filter

/// @brief parse xml
void
SecondaryStructureHasResidueFilter::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	using protocols::parser::BluePrint;
	min_helix_length_ = tag->getOption<core::Size>( "min_helix_length", 4 );
	min_sheet_length_ = tag->getOption<core::Size>( "min_sheet_length", 3 );
	min_loop_length_ = tag->getOption<core::Size>( "min_loop_length", 1 );
	max_helix_length_ = tag->getOption<core::Size>( "max_helix_length", 9999 );
	max_sheet_length_ = tag->getOption<core::Size>( "max_sheet_length", 9999 );
	max_loop_length_ = tag->getOption<core::Size>( "max_loop_length", 9999 );
	filter_helix_ = tag->getOption<bool>( "filter_helix", true );
	filter_sheet_ = tag->getOption<bool>( "filter_sheet", true );
	filter_loop_ = tag->getOption<bool>( "filter_loop", false );
	//chrisk get aa1 res list string here
	req_residue_str_ = tag->getOption< std::string >( "required_restypes", "VILMFYW" );
	threshold_ = tag->getOption< core::Real >( "secstruct_fraction_threshold", 1.0 );
	nres_req_per_secstruct_ = tag->getOption< core::Size >( "nres_required_per_secstruct", 1.0 );

	if ( filter_helix_ ) {
		tr << "filter on helix with length: "<<min_helix_length_<<"-"<< max_helix_length_ << std::endl;
	}

	if ( filter_sheet_ ) {
		tr << "filter on sheet with length: "<<min_sheet_length_<<"-"<< max_sheet_length_ << std::endl;
	}

	if ( filter_loop_ ) {
		tr << "filter on loop with length: "<<min_loop_length_<<"-"<< max_loop_length_ << std::endl;
	}

	std::string const res_taskop_str( tag->getOption< std::string >( "res_check_task_operations", "" ) );
	if ( res_taskop_str.empty() ) res_check_task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );
	else res_check_task_factory_ = protocols::rosetta_scripts::parse_task_operations( res_taskop_str, data );
	std::string const ss_taskop_str( tag->getOption< std::string >( "ss_select_task_operations", "" ) );
	if ( ss_taskop_str.empty() ) ss_select_task_factory_ = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );
	else ss_select_task_factory_ = protocols::rosetta_scripts::parse_task_operations( ss_taskop_str, data );
}

core::Size
SecondaryStructureHasResidueFilter::n_req_res_in_seq(
	std::string const & seq,
	utility::vector0< bool > const & is_checked
) const {
	debug_assert( seq.size() == is_checked.size() );
	if ( seq.size() == 0 || req_residue_str_.size() == 0 ) return 0;
	core::Size rescount = 0;
	//strings are indexed from 0, derpface!
	for ( core::Size i=0; i<=seq.size() - 1; ++i ) {
		//skip if we're ignoring this residue because taskop
		if ( !is_checked[ i ] ) continue;
		for ( core::Size j=0; j<=req_residue_str_.size() - 1; ++j ) {
			if ( seq[ i ] == req_residue_str_[ j ] ) {
				++rescount;
				break; //we can stop looking once we found a match
			}
		}
	}
	return rescount;
}

core::Real SecondaryStructureHasResidueFilter::compute( core::pose::Pose const & pose ) const {


	core::pose::Pose pose_copy=pose;
	core::scoring::dssp::Dssp dssp( pose_copy );
	std::string dssp_ss=dssp.get_dssp_secstruct();
	tr << dssp_ss << std::endl;

	core::Size num_helix_has_residue=0;
	core::Size num_sheet_has_residue=0;
	core::Size num_loop_has_residue=0;
	core::Size num_helix=0;
	core::Size num_sheet=0;
	core::Size num_loop=0;

	//this task selects which residues are to be checked for having the right restype
	core::pack::task::PackerTaskOP task( res_check_task_factory_->create_task_and_apply_taskoperations( pose_copy ) );
	//TODO need another task to select secstruct elements considered
	core::pack::task::PackerTaskOP ss_select_task( ss_select_task_factory_->create_task_and_apply_taskoperations( pose_copy ) );

	std::string ss_seq;
	core::Size ss_len( 0 );
	std::string::const_iterator iter=dssp_ss.begin();
	//boolean mask indexed from 0 to match seq which is a string
	utility::vector0< bool > is_checked;
	while ( iter< dssp_ss.end() ) {
		if ( *iter=='H' ) {
			ss_seq.clear();
			ss_len = 0;
			is_checked.clear();
			core::Size ss_chain( pose_copy.chain( iter - dssp_ss.begin() + 1 ) ); //get chain of this ss element
			//keep iterating until we fall off the end of the helix, accounts for chain endings
			while ( *iter=='H' && iter < dssp_ss.end() && ss_chain == pose_copy.chain( iter - dssp_ss.begin() + 1 ) ) {
				++ss_len;
				core::Size iter_index( iter - dssp_ss.begin() + 1 );
				//TODO skip this res if not in ss selector task
				if ( ss_select_task->residue_task( iter_index ).being_designed() && pose_copy.residue( iter_index ).is_protein() ) {
					//cache aa1 sequence of residue where iterator is pointing only if designable in taskop
					ss_seq += pose_copy.sequence()[ iter_index - 1 ];
					//push_back yes consider this res or no to ss seq bool mask based on taskop
					if ( task->residue_task( iter_index ).being_designed() && pose_copy.residue( iter_index ).is_protein() ) {
						is_checked.push_back( true );
					} else is_checked.push_back( false );
				}
				++iter;
			}
			if ( ss_len >= min_helix_length_ && ss_len <= max_helix_length_ && ss_seq.length() > 0 ) {
				++num_helix;
				//how many instances of requires restype(s) is in this sestruct element?
				if ( n_req_res_in_seq( ss_seq, is_checked ) >= nres_req_per_secstruct_ ) ++num_helix_has_residue;
			}
		} else if ( *iter=='E' ) {
			ss_seq.clear();
			ss_len = 0;
			is_checked.clear();
			core::Size ss_chain( pose_copy.chain( iter - dssp_ss.begin() + 1 ) ); //get chain of this ss element
			//keep iterating until we fall off the end of the sheet, accounts for chain endings
			while ( *iter=='E' && iter < dssp_ss.end() && ss_chain == pose_copy.chain( iter - dssp_ss.begin() + 1 ) ) {
				++ss_len;
				core::Size iter_index( iter - dssp_ss.begin() + 1 );
				//TODO skip this res if not in ss selector task
				if ( ss_select_task->residue_task( iter_index ).being_designed() && pose_copy.residue( iter_index ).is_protein() ) {
					//cache aa1 sequence of residue where iterator is pointing only if designable in taskop
					ss_seq += pose_copy.sequence()[ iter_index - 1 ];
					//push_back yes consider this res or no to ss seq bool mask based on taskop
					if ( task->residue_task( iter_index ).being_designed() && pose_copy.residue( iter_index ).is_protein() ) {
						is_checked.push_back( true );
					} else is_checked.push_back( false );
				}
				++iter;
			}
			if ( ss_len >= min_sheet_length_ && ss_len <= max_sheet_length_ && ss_seq.length() > 0 ) {
				++num_sheet;
				//how many instances of requires restype(s) is in this sestruct element?
				if ( n_req_res_in_seq( ss_seq, is_checked ) >= nres_req_per_secstruct_ ) ++num_sheet_has_residue;
			}
		} else if ( *iter=='L' ) {
			ss_seq.clear();
			ss_len = 0;
			is_checked.clear();
			core::Size ss_chain( pose_copy.chain( iter - dssp_ss.begin() + 1 ) ); //get chain of this ss element
			//keep iterating until we fall off the end of the loop, accounts for chain endings
			while ( *iter=='L' && iter < dssp_ss.end() && ss_chain == pose_copy.chain( iter - dssp_ss.begin() + 1 ) ) {
				++ss_len;
				core::Size iter_index( iter - dssp_ss.begin() + 1 );
				//TODO skip this res if not in ss selector task
				if ( ss_select_task->residue_task( iter_index ).being_designed() && pose_copy.residue( iter_index ).is_protein() ) {
					//cache aa1 sequence of residue where iterator is pointing only if designable in taskop
					ss_seq += pose_copy.sequence()[ iter_index - 1 ];
					//push_back yes consider this res or no to ss seq bool mask based on taskop
					if ( task->residue_task( iter_index ).being_designed() && pose_copy.residue( iter_index ).is_protein() ) {
						is_checked.push_back( true );
					} else is_checked.push_back( false );
				}
				++iter;
			}
			if ( ss_len >= min_loop_length_ && ss_len <= max_loop_length_ && ss_seq.length() > 0 ) {
				++num_loop;
				//how many instances of requires restype(s) is in this sestruct element?
				if ( n_req_res_in_seq( ss_seq, is_checked ) >= nres_req_per_secstruct_ ) ++num_loop_has_residue;
			}
		}
	}

	tr << " Pose has " << num_helix_has_residue << " / " << num_helix << " helix, " << num_sheet_has_residue  << " / " << num_sheet << " sheet, " << num_loop_has_residue << " / " << num_loop << " loop with enough required residues, according to dssp_reduced definition" <<  std::endl;
	//compute fraction of total ss elements with enough required residue
	Size num_secstruct( 0 );
	Size num_secstruct_has_residue( 0 );
	if ( filter_helix_ ) {
		num_secstruct_has_residue += num_helix_has_residue;
		num_secstruct += num_helix;
	}
	if ( filter_sheet_ ) {
		num_secstruct_has_residue += num_sheet_has_residue;
		num_secstruct += num_sheet;
	}
	if ( filter_loop_ ) {
		num_secstruct_has_residue += num_loop_has_residue;
		num_secstruct += num_loop;
	}
	return ( static_cast< core::Real >( num_secstruct_has_residue ) / static_cast< core::Real >( num_secstruct ) );
}

core::Real SecondaryStructureHasResidueFilter::report_sm( core::pose::Pose const & pose ) const {
	return compute( pose );
}

void SecondaryStructureHasResidueFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	//report message is inside compute function already
	compute( pose );
	out << "See SecondaryStructureHasResidueFilter tracer log for report" << std::endl;
}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP SecondaryStructureHasResidueFilterCreator::create_filter() const { return protocols::filters::FilterOP( new SecondaryStructureHasResidueFilter ); }

// XRW TEMP std::string
// XRW TEMP SecondaryStructureHasResidueFilterCreator::keyname() const { return "SecondaryStructureHasResidue"; }

std::string SecondaryStructureHasResidueFilter::name() const {
	return class_name();
}

std::string SecondaryStructureHasResidueFilter::class_name() {
	return "SecondaryStructureHasResidue";
}

void SecondaryStructureHasResidueFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "min_helix_length", xsct_non_negative_integer, "XRW TO DO", "4" )
		+ XMLSchemaAttribute::attribute_w_default( "min_sheet_length", xsct_non_negative_integer, "XRW TO DO", "3" )
		+ XMLSchemaAttribute::attribute_w_default( "min_loop_length", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "max_helix_length", xsct_non_negative_integer, "XRW TO DO", "9999" )
		+ XMLSchemaAttribute::attribute_w_default( "max_sheet_length", xsct_non_negative_integer, "XRW TO DO", "9999" )
		+ XMLSchemaAttribute::attribute_w_default( "max_loop_length", xsct_non_negative_integer, "XRW TO DO", "9999" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_helix", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_sheet", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "filter_loop", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "required_restypes", xs_string, "XRW TO DO", "VILMFYW" )
		+ XMLSchemaAttribute::attribute_w_default( "secstruct_fraction_threshold", xsct_real, "XRW TO DO", "1.0" )
		+ XMLSchemaAttribute::attribute_w_default( "nres_required_per_secstruct", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute( "res_check_task_operations", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "ss_select_task_operations", xs_string, "XRW TO DO" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string SecondaryStructureHasResidueFilterCreator::keyname() const {
	return SecondaryStructureHasResidueFilter::class_name();
}

protocols::filters::FilterOP
SecondaryStructureHasResidueFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SecondaryStructureHasResidueFilter );
}

void SecondaryStructureHasResidueFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SecondaryStructureHasResidueFilter::provide_xml_schema( xsd );
}



} // filters
} // fldsgn
} // protocols
