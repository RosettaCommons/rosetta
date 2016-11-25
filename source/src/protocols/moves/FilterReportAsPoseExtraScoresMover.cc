// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/moves/FilterReportAsPoseExtraScoresMover.cc
/// @brief This Mover runs a Filter and dumps the report_sm value to Pose's extra scores (for later JD reporting)
/// @author Steven Lewis (smlewi@gmail.com)

// Unit headers
#include <protocols/moves/FilterReportAsPoseExtraScoresMover.hh>
#include <protocols/moves/FilterReportAsPoseExtraScoresMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <protocols/filters/Filter.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

#include <string>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.FilterReportAsPoseExtraScoresMover" );

namespace protocols {
namespace moves {

FilterReportAsPoseExtraScoresMover::FilterReportAsPoseExtraScoresMover():
	protocols::moves::Mover( FilterReportAsPoseExtraScoresMover::mover_name() )
{

}

FilterReportAsPoseExtraScoresMover::FilterReportAsPoseExtraScoresMover(
	protocols::filters::FilterOP filter,
	std::string report_as ):
	protocols::moves::Mover( FilterReportAsPoseExtraScoresMover::mover_name() ),
	filter_(std::move(filter)),
	report_as_(std::move(report_as))
{

}

FilterReportAsPoseExtraScoresMover::~FilterReportAsPoseExtraScoresMover()= default;

void
FilterReportAsPoseExtraScoresMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	//borrowed from GenericMonteCarloMover.cc:parse_my_tag
	std::string const filter_name( tag->getOption< std::string >( "filter_name", "true_filter" ) );
	auto find_filter( filters.find( filter_name ) );
	if ( find_filter == filters.end() ) {
		TR.Error << "ERROR !! filter not found in map: \n" << tag << std::endl;
		runtime_assert( find_filter != filters.end() );
	}
	set_filter(find_filter->second/*->clone()*/); // I haven't a clue if this is safe or if it should be a clone operation!
	if ( tag->hasOption( "report_as" ) ) {
		set_report_as(tag->getOption< std::string >( "report_as" ));
	} else {
		throw utility::excn::EXCN_BadInput(mover_name() + " requires report_as to know what to report its value to");
	}
}

protocols::moves::MoverOP
FilterReportAsPoseExtraScoresMover::clone() const
{
	return protocols::moves::MoverOP( new FilterReportAsPoseExtraScoresMover( *this ) );
}

protocols::moves::MoverOP
FilterReportAsPoseExtraScoresMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new FilterReportAsPoseExtraScoresMover );
}

// XRW TEMP std::string
// XRW TEMP FilterReportAsPoseExtraScoresMover::get_name() const
// XRW TEMP {
// XRW TEMP  return FilterReportAsPoseExtraScoresMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FilterReportAsPoseExtraScoresMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "FilterReportAsPoseExtraScoresMover";
// XRW TEMP }

void
FilterReportAsPoseExtraScoresMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, FilterReportAsPoseExtraScoresMover const & mover )
{
	mover.show(os);
	return os;
}

void
FilterReportAsPoseExtraScoresMover::apply( core::pose::Pose & pose)
{
	//ensure filter exists
	if ( !filter_ ) {
		throw utility::excn::EXCN_NullPointer(mover_name() + " has no Filter.");
	}

	//extract report_sm from filter to setPoseExtraScores
	core::pose::setPoseExtraScore(
		pose,
		get_report_as(), //ideally this would be the filter name, but that will get overwritten
		filter_->report_sm(pose)
	);
}

/// @brief set filter we report from
void FilterReportAsPoseExtraScoresMover::set_filter( protocols::filters::FilterOP filter ) { filter_ = filter; }

void FilterReportAsPoseExtraScoresMover::set_report_as( std::string const & report_as ) { report_as_ = report_as; }

/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP FilterReportAsPoseExtraScoresMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new FilterReportAsPoseExtraScoresMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP FilterReportAsPoseExtraScoresMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return FilterReportAsPoseExtraScoresMover::mover_name();
// XRW TEMP }

std::string FilterReportAsPoseExtraScoresMover::get_name() const {
	return mover_name();
}

std::string FilterReportAsPoseExtraScoresMover::mover_name() {
	return "FilterReportAsPoseExtraScoresMover";
}

void FilterReportAsPoseExtraScoresMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	attlist
		//TECHNICALLY has a default, but useless if not specified
		+ XMLSchemaAttribute::required_attribute( "filter_name", xs_string, "Filter (defined in FILTERS section) to report value of")
		+ XMLSchemaAttribute::required_attribute( "report_as", xs_string, "string to report value to (cannot use filter name because RosettaScripts is stupid and reports the filter's ENDING value there already)");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Reports the report_sm value of the provided Filter at the time this Mover runs (normally RosettaScripts reports at the end of the run, when the Filter's value will have changed)", attlist );
}

std::string FilterReportAsPoseExtraScoresMoverCreator::keyname() const {
	return FilterReportAsPoseExtraScoresMover::mover_name();
}

protocols::moves::MoverOP
FilterReportAsPoseExtraScoresMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new FilterReportAsPoseExtraScoresMover );
}

void FilterReportAsPoseExtraScoresMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FilterReportAsPoseExtraScoresMover::provide_xml_schema( xsd );
}


} //protocols
} //moves
