// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/SSElementLengthFilter
/// @brief filter structures by longest,shortest or avg length of a given secondary structure type
/// @author TJ Brunette

// Unit Headers
#include <protocols/simple_filters/SSElementLengthFilter.hh>
#include <protocols/simple_filters/SSElementLengthFilterCreator.hh>

#include <numeric/xyzVector.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

//#include <core/pack/task/TaskFactory.hh>

#include <core/pose/util.hh>
#include <core/pose/subpose_manipulation_util.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

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
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_filters/InterfaceSasaFilter.hh>
#include <protocols/simple_moves/MutateResidue.hh>


#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <set>
#include <map>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

//// C++ headers
static basic::Tracer tr("protocols.filters.SSElementLengthFilter");

namespace protocols {
namespace simple_filters {

using utility::vector1;
using core::Size;
using core::Real;
using std::set;
using namespace core;
// @brief default constructor
SSElementLengthFilter::SSElementLengthFilter():
	Filter( "SSElementLengthFilter" ),
	selector_(utility::pointer::make_shared< core::select::residue_selector::TrueResidueSelector >())
{}

// @brief copy constructor
SSElementLengthFilter::SSElementLengthFilter( SSElementLengthFilter const & rval ):
	Super( rval ),
	report_avg_(rval.report_avg_),
	report_longest_(rval.report_longest_),
	report_shortest_(rval.report_shortest_),
	selector_(rval.selector_)
{}

// @brief destructor
SSElementLengthFilter::~SSElementLengthFilter() = default;


/// @brief
SSElementLengthFilter::Real
SSElementLengthFilter::report_sm( const Pose & pose ) const
{
	return  compute( pose );
}

/// @brief
void
SSElementLengthFilter::report( std::ostream & out, Pose const & pose ) const
{
	out << "Length stat:" <<  compute( pose ) << std::endl;
}

/// @brief get ss_elements
protocols::loops::Loops SSElementLengthFilter::get_ss_elements(const Pose & pose) const{
	using protocols::loops::Loop;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	char lastSecStruct = dssp.get_dssp_secstruct( 1 );
	Size startSS = 0;
	Size endSS = 0;

	Size first_residue = 1;
	Size last_residue = pose.size();
	core::select::residue_selector::ResidueSubset residues( pose.total_residue(), false );
	residues = selector_->apply( pose );
	bool first_found =false;
	bool last_found=false;
	for ( Size ii=1; ii<=residues.size(); ++ii ) {
		if ( residues[ii]==true && first_found==false ) {
			first_residue=ii;
			first_found=true;
		}
		if ( residues[ii]==false && first_found==true && last_found==false ) {
			last_residue=ii-1;
			last_found=true;
		}
	}
	protocols::loops::Loops ss_elements;
	if ( dssp.get_dssp_secstruct(first_residue)=='H' || dssp.get_dssp_secstruct(first_residue)=='E' ) {
		startSS=  first_residue;
	}
	for ( core::Size ii = first_residue+1; ii <= last_residue; ++ii ) {
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




/// @brief
Real
SSElementLengthFilter::compute( const Pose & pose ) const
{
	using protocols::loops::Loop;
	protocols::loops::Loops ss_elements = get_ss_elements(pose);
	Size startElement = 1;
	Size endElement = ss_elements.size();
	Size shortest_tmp= 999;
	Size longest_tmp=0;
	Size total_length_tmp=0;
	for ( Size ii=startElement; ii<=endElement; ++ii ) {
		Size length_tmp = ss_elements[ii].stop()-ss_elements[ii].start()+1;
		total_length_tmp+=length_tmp;
		if ( length_tmp<shortest_tmp ) {
			shortest_tmp=length_tmp;
		}
		if ( length_tmp>longest_tmp ) {
			longest_tmp = length_tmp;
		}
	}
	Real avg_length = ((Real)total_length_tmp)/(Real)(ss_elements.size());
	if ( report_longest_ ) {
		return(longest_tmp);
	}
	if ( report_shortest_ ) {
		return(shortest_tmp);
	}
	return(avg_length);
}


// @brief returns true if the given pose passes the filter, false otherwise.
// In this case, the test is whether the give pose is the topology we want.
bool SSElementLengthFilter::apply(const Pose & pose ) const
{
	Real value = compute( pose );
	tr << "value" << value << "threshold" <<  threshold_ << std::endl;
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
SSElementLengthFilter::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data ,
	filters::Filters_map const &,
	moves::Movers_map const &,
	Pose const & )
{
	// set threshold
	threshold_ = tag->getOption<Real>("threshold",-999999);
	report_avg_ = tag->getOption<bool>( "report_avg", true );
	report_shortest_ =  tag->getOption<bool>( "report_shortest", false );
	report_longest_ =  tag->getOption<bool>( "report_longest", false );
	selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data );
	if ( ! selector_ ) {
		selector_ = utility::pointer::make_shared< core::select::residue_selector::TrueResidueSelector >();
	}


}

// XRW TEMP filters::FilterOP
// XRW TEMP SSElementLengthFilterCreator::create_filter() const { return protocols::filters::FilterOP(new SSElementLengthFilter); }

// XRW TEMP std::string
// XRW TEMP SSElementLengthFilterCreator::keyname() const { return "SSBisectddGFilter"; }

std::string SSElementLengthFilter::name() const {
	return class_name();
}

std::string SSElementLengthFilter::class_name() {
	return "SSElementLengthFilter";
}

void SSElementLengthFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default("threshold", xsct_real, "cutoff length", "-999999")
		+ XMLSchemaAttribute::attribute_w_default("report_avg", xsct_rosetta_bool, "reports avg", "true")
		+ XMLSchemaAttribute::attribute_w_default("report_longest", xsct_rosetta_bool, "reports longest", "false")
		+ XMLSchemaAttribute::attribute_w_default("report_shortest", xsct_rosetta_bool, "reports shortest", "false");
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "parses residue selector" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "reports attributes of secondary structure. Currently SSElement is defined as helix OR sheet", attlist );

}

std::string SSElementLengthFilterCreator::keyname() const {
	return SSElementLengthFilter::class_name();
}

protocols::filters::FilterOP
SSElementLengthFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( utility::pointer::make_shared<SSElementLengthFilter>() );
}

void SSElementLengthFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SSElementLengthFilter::provide_xml_schema( xsd );
}

} // filters
} // protocols
