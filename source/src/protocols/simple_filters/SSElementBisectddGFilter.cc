// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SSElementBisectddGFilter
/// @brief filter structures by IntraRepeatContacts
/// @details Disconnects a protein by secondary structure elements and calculates the DDG between the elements
/// @author TJ Brunette

// Unit Headers
#include <protocols/simple_filters/SSElementBisectddGFilter.hh>
#include <protocols/simple_filters/SSElementBisectddGFilterCreator.hh>

#include <numeric/xyzVector.hh>

// Project Headers
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/pose/subpose_manipulation_util.hh>

#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>


#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/motif/reference_frames.hh>

#include <basic/options/option_macros.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <numeric/xyzTransform.hh>
// Parser headers
#include <protocols/filters/Filter.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_filters/DdgFilter.hh>

#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <set>
#include <map>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers
static THREAD_LOCAL basic::Tracer tr("protocols.filters.SSElementBisectddGFilter");

namespace protocols {
namespace simple_filters {

using utility::vector1;
using core::Size;
using core::Real;
using std::set;
using namespace core;
// @brief default constructor
SSElementBisectddGFilter::SSElementBisectddGFilter():
	Filter( "SSBisectddGFilter" ),
	filtered_value_( 99 )
{}

// @brief copy constructor
SSElementBisectddGFilter::SSElementBisectddGFilter( SSElementBisectddGFilter const & rval ):
	Super( rval ),
	filtered_value_( rval.filtered_value_ )
{}

// @brief destructor
SSElementBisectddGFilter::~SSElementBisectddGFilter() {}

// @brief set filtered value
void SSElementBisectddGFilter::filtered_value( Real const & value )
{
	filtered_value_ = value;
}

/// @brief
SSElementBisectddGFilter::Real
SSElementBisectddGFilter::report_sm( const Pose & pose ) const
{
	return  compute( pose );
}

/// @brief
void
SSElementBisectddGFilter::report( std::ostream & out, Pose const & pose ) const
{
	if ( report_avg_ ) {
		out << "Avg ddg at post SSE bisections" <<  compute( pose ) << std::endl;
	} else {
		out << "Worst ddg at post SSE bisections" <<  compute( pose ) << std::endl;
	}
}

/// @brief get ss_elements
protocols::loops::Loops SSElementBisectddGFilter::get_ss_elements(const Pose & pose) const{
	using protocols::loops::Loop;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	char lastSecStruct = dssp.get_dssp_secstruct( 1 );
	Size startSS = 0;
	Size endSS = 0;
	Size nres1 = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		nres1 = core::pose::symmetry::symmetry_info(pose)->num_total_residues_without_pseudo();
		if ( option[OptionKeys::score::motif_ignore_symmmetry]() ) {
			nres1= core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
		}
	}
	protocols::loops::Loops ss_elements;
	if ( dssp.get_dssp_secstruct(1)=='H' || dssp.get_dssp_secstruct(1)=='E' ) {
		startSS=  1;
	}
	for ( core::Size ii = 2; ii <= nres1; ++ii ) {
		if ( (dssp.get_dssp_secstruct(ii) == 'H' || dssp.get_dssp_secstruct(ii)) && lastSecStruct == 'L' ) {
			startSS = ii;
		}
		if ( dssp.get_dssp_secstruct(ii)!=lastSecStruct && (lastSecStruct ==  'H' || lastSecStruct ==  'E') ) {
			endSS = ii-1;
			if ( endSS-startSS >= 2 ) {
				ss_elements.add_loop(Loop(startSS,endSS));
			}
		}
		lastSecStruct = dssp.get_dssp_secstruct(ii);
	}
	return(ss_elements);
}





Real SSElementBisectddGFilter::get_ddg_bisect_score(Size element, protocols::loops::Loops ssElements, const Pose & pose) const{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	Size pose_length = pose.size();
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		pose_length= core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	}
	protocols::simple_filters::DdgFilter ddg_filter( 100/*ddg_threshold*/, scorefxn_, 1 /*jump*/, 1 /*repeats*/ );
	//ddg_filter.repack( repack() );
	core::pose::Pose nTerm_pose;
	core::pose::Pose cTerm_pose;
	utility::vector1<Size> nTerm_positions;
	utility::vector1<Size> cTerm_positions;
	for ( Size ii=1; ii<=ssElements[element].stop(); ++ii ) {
		nTerm_positions.push_back(ii);
	}
	for ( Size ii=ssElements[element+1].start(); ii<=pose_length; ++ii ) {
		cTerm_positions.push_back(ii);
	}
	pdbslice(nTerm_pose,pose,nTerm_positions);
	pdbslice(cTerm_pose,pose,cTerm_positions);
	append_pose_to_pose(nTerm_pose,cTerm_pose,true);
	nTerm_pose.dump_pdb("splitPose.pdb");
	Real ddG = ddg_filter.compute( nTerm_pose );
	return(ddG);
}


/// @brief
Real
SSElementBisectddGFilter::compute( const Pose & pose ) const
{
	using protocols::loops::Loop;
	protocols::loops::Loops ss_elements = get_ss_elements(pose);
	Size startElement = 1;
	Size endElement = ss_elements.size();
	if ( ignore_terminal_SS_>0 ) {
		startElement+=ignore_terminal_SS_;
		endElement-=ignore_terminal_SS_;
	}
	Real score = 0;
	if ( only_n_term_ ) {
		return(get_ddg_bisect_score(1,ss_elements,pose));
	}
	if ( only_c_term_ ) {
		return(get_ddg_bisect_score(ss_elements.size()-1,ss_elements,pose));
	}
	if ( report_avg_ ) {
		Real tmpScore = 0;
		Size ct = 0;
		for ( Size ii=startElement; ii<=endElement; ++ii ) {
			tmpScore+=get_ddg_bisect_score(ii,ss_elements,pose);
			ct+=1;
		}
		if ( ct!=0 ) {
			score = tmpScore/Real(ct);
		}
	} else {
		Real tmpScore = -999999;
		for ( Size ii=startElement; ii<=endElement; ++ii ) {
			Real tmp=get_ddg_bisect_score(ii,ss_elements,pose);
			if ( tmp>tmpScore ) {
				tmpScore=tmp;
			}
		}
		score= tmpScore;
	}
	return(score);
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool SSElementBisectddGFilter::apply(const Pose & pose ) const
{
	Real value = compute( pose );
	tr << "value" << value << "filtered_value_" << threshold_ << std::endl;
	if ( value >= threshold_ ) {
		tr << "Successfully filtered: " << value << std::endl;
		return true;
	} else {
		tr << "Filter failed current/threshold=" << value << "/" << threshold_ << std::endl;
		return false;
	}
} // apply_filter

/// @brief parse xml
void
SSElementBisectddGFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const &,
	Movers_map const &,
	Pose const & )
{
	// set threshold
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	threshold_ = tag->getOption<Real>("threshold",-999999);
	report_avg_ = tag->getOption<bool>( "report_avg", true );
	ignore_terminal_SS_ = tag->getOption<Size>("ignore_terminal_ss",0);
	only_n_term_ = tag->getOption<bool>("only_n_term",false);
	only_c_term_ = tag->getOption<bool>("only_c_term",false);
}

// XRW TEMP filters::FilterOP
// XRW TEMP SSElementBisectddGFilterCreator::create_filter() const { return protocols::filters::FilterOP(new SSElementBisectddGFilter); }

// XRW TEMP std::string
// XRW TEMP SSElementBisectddGFilterCreator::keyname() const { return "SSBisectddGFilter"; }

std::string SSElementBisectddGFilter::name() const {
	return class_name();
}

std::string SSElementBisectddGFilter::class_name() {
	return "SSBisectddGFilter";
}

void SSElementBisectddGFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default("threshold", xsct_real, "XRW TO DO", "-999999")
		+ XMLSchemaAttribute::attribute_w_default("report_avg", xsct_rosetta_bool, "XRW TO DO", "true")
		+ XMLSchemaAttribute::attribute_w_default("ignore_terminal_ss", xsct_non_negative_integer, "XRW TO DO", "0")
		+ XMLSchemaAttribute::attribute_w_default("only_n_term", xsct_rosetta_bool, "XRW TO DO", "false")
		+ XMLSchemaAttribute::attribute_w_default("only_c_term", xsct_rosetta_bool, "XRW TO DO", "false");

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string SSElementBisectddGFilterCreator::keyname() const {
	return SSElementBisectddGFilter::class_name();
}

protocols::filters::FilterOP
SSElementBisectddGFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new SSElementBisectddGFilter );
}

void SSElementBisectddGFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SSElementBisectddGFilter::provide_xml_schema( xsd );
}

} // filters
} // protocols
