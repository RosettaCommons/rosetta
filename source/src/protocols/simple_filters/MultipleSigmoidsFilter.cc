// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/MultipleSigmoidsFilter.cc
/// @brief
/// @author Sarel Fleishman and Shira Warszawski (shiraw1@weizmann.ac.il)

//Unit Headers
#include <protocols/simple_filters/MultipleSigmoidsFilter.hh>
#include <protocols/simple_filters/MultipleSigmoidsFilterCreator.hh>
#include <protocols/simple_filters/RelativePoseFilter.hh>
#include <protocols/simple_filters/SigmoidFilter.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/string_util.hh>
#include <protocols/filters/BasicFilters.hh>
#include <limits>
#include <protocols/simple_filters/OperatorFilter.hh>
#include <protocols/simple_filters/RelativePoseFilter.hh>
#include <protocols/rosetta_scripts/util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.MultipleSigmoids" );
using namespace protocols::filters;

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP MultipleSigmoidsFilterCreator::create_filter() const { return protocols::filters::FilterOP( new MultipleSigmoids ); }

// XRW TEMP std::string
// XRW TEMP MultipleSigmoidsFilterCreator::keyname() const { return "MultipleSigmoids"; }

MultipleSigmoids::MultipleSigmoids() :
	protocols::filters::Filter( "MultipleSigmoids" ),
	file_names_ ( "" ), // dflt ""
	threshold_(0.0),
	r_pose_( /* NULL */ ),
	sig_( /* NULL */ ),
	operatorF_( /* NULL */ )
{
}

MultipleSigmoids::~MultipleSigmoids() = default;

void
MultipleSigmoids::reset_baseline( core::pose::Pose const & pose, bool const attempt_read_from_checkpoint ){
	operatorF_->reset_baseline( pose, attempt_read_from_checkpoint );
	TR<<"MultipleSigmoids: reset baseline"<<std::endl;
}

void
MultipleSigmoids::operator_filter( OperatorOP opt ){ operatorF_ = opt; }

OperatorOP
MultipleSigmoids::operator_filter() const{ return operatorF_; }

SigmoidOP
MultipleSigmoids::sigmoid_filter() const{ return sig_; }

void
MultipleSigmoids::relative_pose_filter( RelativePoseFilterOP rpose ){ r_pose_ = rpose; }

RelativePoseFilterOP
MultipleSigmoids::relative_pose_filter() const{ return r_pose_; }


void
MultipleSigmoids::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, filters::Filters_map const &filters, moves::Movers_map const &movers, core::pose::Pose const & pose )
{
	threshold( tag->getOption< core::Real >( "threshold", 0 ) );
	utility::vector1< std::string > const pdb_names( utility::string_split( tag->getOption< std::string >( "file_names" ), ',' ) ); //split file names
	operatorF_ = OperatorOP( new protocols::simple_filters::Operator );
	utility::vector1< utility::tag::TagCOP > const sub_tags( tag->getTags() ); //tags
	for ( utility::tag::TagCOP const sub_tag : sub_tags ) {
		if ( sub_tag->getName() == "Operator" ) {
			operatorF_->parse_my_tag( sub_tag, data, filters, movers, pose ) ;
		}
	}

	for ( std::string const & fname : pdb_names ) {
		for ( utility::tag::TagCOP const sub_tag : sub_tags ) {
			if ( sub_tag->getName() == "RelativePose" ) {
				r_pose_ = RelativePoseFilterOP( new RelativePoseFilter );
				TR<<"I'm now reading from RelativePose filter"<<std::endl;
				r_pose_->pdb_name(fname);
				r_pose_->parse_my_tag(sub_tag, data, filters, movers, pose);
			} else if ( sub_tag->getName() == "Sigmoid" ) {
				sig_ = SigmoidOP( new Sigmoid );
				TR<<"I'm now reading from Sigmoid filter for fname "<<fname<<std::endl;
				sig_->set_user_defined_name( fname );
				sig_->filter(r_pose_);
				sig_->parse_my_tag(sub_tag, data, filters, movers, pose);
			} else if ( sub_tag->getName() == "Operator" ) {
				operatorF_->add_filter(sig_);
				TR<<"Adding filter to the operator"<<std::endl;
			} else {
				utility_exit_with_message( "MultipleSigmoids subtag not recognized: " + sub_tag->getName() );
			}
		} //tags foreach
	} //pdbs foreach
} //parse_my_tag

bool MultipleSigmoids::apply( core::pose::Pose const & pose ) const
{
	core::Real const val ( compute( pose ) );
	return( val <= threshold() );
}
void
MultipleSigmoids::report( std::ostream &o, core::pose::Pose const & pose ) const {
	core::Real const val = compute( pose );
	o << "Multiplesigmoids returns "<<val<<std::endl;
}
core::Real
MultipleSigmoids::report_sm( core::pose::Pose const & pose ) const {
	operatorF_->report_sm( pose );
	return( compute( pose ) );
}
core::Real
MultipleSigmoids::compute(
	core::pose::Pose const & pose
) const
{
	core::Real const val( operatorF_->compute( pose ) );
	TR<<"filter MultipleSigmoids returns "<<val<<std::endl;
	return( val );
}

protocols::filters::FilterOP
MultipleSigmoids::clone() const{
	return protocols::filters::FilterOP( new MultipleSigmoids( *this ) );
}

protocols::filters::FilterOP
MultipleSigmoids::fresh_instance() const{
	return protocols::filters::FilterOP( new MultipleSigmoids() );
}

std::string MultipleSigmoids::name() const {
	return class_name();
}

std::string MultipleSigmoids::class_name() {
	return "MultipleSigmoids";
}

std::string subtag_for_multiple_sigmoids( std::string const & foo ) {
	return "multiple_sigmoids_subtag_" + foo + "_type";
}

void MultipleSigmoids::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// FUCK
	using namespace utility::tag;
	AttributeList attlist;
	// THIS
	attlist + XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Value above which the filter fails", "0" )
		+ XMLSchemaAttribute( "file_names", xs_string, "PDBs to apply to" );
	// NOISE

	AttributeList relative_pose_subtag_attributes, sigmoid_subtag_attributes, operator_subtag_attributes;

	RelativePoseFilter::attributes( relative_pose_subtag_attributes );
	Sigmoid::attributes( sigmoid_subtag_attributes );
	Operator::attributes( operator_subtag_attributes );

	utility::tag::XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "RelativePose", relative_pose_subtag_attributes, "Tags describing RelativePose filters to be applied"/*, 0 minoccurs*/ )
		.add_simple_subelement( "Sigmoid", sigmoid_subtag_attributes, "Tags describing Sigmoid filters to be applied"/*, 0 minoccurs*/ )
		.add_simple_subelement( "Operator", operator_subtag_attributes, "Tags describing Operators to be applied"/*, 0 minoccurs*/ )
		.complex_type_naming_func( & subtag_for_multiple_sigmoids );

	protocols::filters::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, class_name(), "Apply an architecture of sigmoid and relative pose filters to a set of PDBs", attlist, ssl );
}

std::string MultipleSigmoidsFilterCreator::keyname() const {
	return MultipleSigmoids::class_name();
}

protocols::filters::FilterOP
MultipleSigmoidsFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new MultipleSigmoids );
}

void MultipleSigmoidsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MultipleSigmoids::provide_xml_schema( xsd );
}


}
}
