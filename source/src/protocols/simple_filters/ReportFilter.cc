// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ReportFilter.cc
/// @brief
/// @author Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/ReportFilter.hh>
#include <protocols/simple_filters/ReportFilterCreator.hh>
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <protocols/jd2/util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
namespace protocols {
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "protocols.simple_filters.ReportFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP ReportFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ReportFilter ); }

// XRW TEMP std::string
// XRW TEMP ReportFilterCreator::keyname() const { return "Report"; }

//default ctor
ReportFilter::ReportFilter() :
	protocols::filters::Filter( "Report" ),
	report_string_( /* NULL */ ),
	filter_( /* NULL */ ),
	report_filter_name_( "" ),
	filter_val_( -9999.9 ),
	checkpointing_file_("")
{}

ReportFilter::~ReportFilter() = default;

void
ReportFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & )
{
	report_filter_name_ = tag->getOption< std::string >( "name" );
	if ( tag->hasOption( "report_string" ) ) {
		report_string_ = data.get_ptr< basic::datacache::DataMapObj< std::string > >( "report_string", tag->getOption< std::string >( "report_string" )  );
	}
	if ( tag->hasOption( "filter" ) ) {
		filter( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "filter" ), filters ) );
	}
	checkpointing_file( tag->getOption< std::string >( "checkpointing_file", "" ) );
	TR<<"checkpointing_file: "<<checkpointing_file()<<std::endl;
}

void
ReportFilter::checkpoint_read() const{
	if ( checkpointing_file() == "" ) {
		return;
	}
	utility::io::izstream data_in( checkpointing_file_ );
	if ( !data_in ) {
		TR<<"checkpointing file not yet written. Skipping reading checkpoint"<<std::endl;
	}
	TR<<"Loading filter_val from checkpointing file: "<<checkpointing_file()<<std::endl;
	std::string line;
	getline( data_in, line );
	std::istringstream stream( line );
	stream >> filter_val_;
	TR<<"Filter_val read from checkpoint: "<<filter_val()<<std::endl;
}

void ReportFilter::checkpoint_write() const{
	if ( checkpointing_file() == "" ) {
		return;
	}
	TR<<"writing filter_val "<<filter_val()<<" to checkpointing file: "<<checkpointing_file()<<std::endl;
	utility::io::ozstream data_out( checkpointing_file_ );
	if ( !data_out ) {
		utility_exit_with_message( "Could not open "+checkpointing_file()+" for writing checkpoint information" );
	}
	data_out<<filter_val();
	TR<<"Wrote filter_val to checkpointing file"<<std::endl;

}

bool
ReportFilter::apply( core::pose::Pose const & pose ) const {
	filter_val_ = compute( pose );
	checkpoint_read();
	checkpoint_write();
	return true;/// does not filter; only reports
}

void
ReportFilter::report( std::ostream & out, core::pose::Pose const & /*pose*/ ) const {
	out<<"filter: "<<report_filter_name_;
	std::string job_name ( protocols::jd2::current_output_name() );
	if ( report_string_ != nullptr && report_string_->obj.length() > 0 ) {
		out<<"job name: "<<job_name<<" report_string: "<<report_string_->obj<<std::endl;
	}
	if ( filter_ != nullptr ) {
		out<<"job name: "<<job_name<<" reporting filter value: ";
	}
	out<<"Internal filter's value is: "<<filter_val()<<std::endl;
}

core::Real
ReportFilter::report_sm( core::pose::Pose const & /*pose*/ ) const {
	checkpoint_read();
	return( filter_val() );
}

core::Real
ReportFilter::compute(
	core::pose::Pose const & pose
) const {
	return( filter()->report_sm( pose ) );
}

void ReportFilter::report_string( std::string const & f )
{
	report_string_->obj = f;
}

std::string ReportFilter::report_string() const
{
	return report_string_->obj;
}

std::string ReportFilter::name() const {
	return class_name();
}

std::string ReportFilter::class_name() {
	return "Report";
}

void ReportFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute(
		"name", xs_string,
		"Name for tag the filter to be reported")
		+ XMLSchemaAttribute(
		"report_string", xs_string,
		"name of an object on the datamap that stores a value for reporting. "
		"This requires another mover/filter to be aware of this object and modify "
		"it. Currently no movers/filters use this functionality, "
		"but it could come in useful in future")
		+ XMLSchemaAttribute(
		"filter", xs_string,
		"name of a filter on the datamap that report will invoke")
		+ XMLSchemaAttribute::attribute_w_default(
		"checkpointing_file", xs_string,
		"If the protocol is checkpointed (e.g., through GenericMonteCarlo) "
		"this will make ReportFilter checkpoint its data. If the checkpointing "
		"file exists the value from the checkpointing file will be read into "
		"ReportFilter's internal value and will be reported at the end of the "
		"run. On apply, the filter's value will be written to the checkpointing file.",
		"");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"This filter reports the value of another filter with the current job "
		"name. Useful when running long trajectories where one wants to see "
		"intermediate values of successful trajectories",
		attlist );
}

std::string ReportFilterCreator::keyname() const {
	return ReportFilter::class_name();
}

protocols::filters::FilterOP
ReportFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ReportFilter );
}

void ReportFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ReportFilter::provide_xml_schema( xsd );
}


}
}
