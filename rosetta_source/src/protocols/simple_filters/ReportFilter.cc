// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ReportFilter.cc
/// @brief
/// @author Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/ReportFilter.hh>
#include <protocols/simple_filters/ReportFilterCreator.hh>
#include <protocols/moves/DataMapObj.hh>
#include <protocols/moves/DataMap.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/rosetta_scripts/util.hh>
namespace protocols{
namespace simple_filters {

using namespace core;
using namespace core::scoring;

static basic::Tracer TR( "protocols.simple_filters.ReportFilter" );

protocols::filters::FilterOP
ReportFilterCreator::create_filter() const { return new ReportFilter; }

std::string
ReportFilterCreator::keyname() const { return "Report"; }

//default ctor
ReportFilter::ReportFilter() :
protocols::filters::Filter( "Report" ),
report_string_( NULL ),
filter_( NULL ),
report_filter_name_( "" )
{}

ReportFilter::~ReportFilter() {}

void
ReportFilter::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap & data, filters::Filters_map const &filters, moves::Movers_map const &, core::pose::Pose const & )
{
	report_filter_name_ = tag->getOption< std::string >( "name" );
	if( tag->hasOption( "report_string" ) )
		report_string_ = data.get< protocols::moves::DataMapObj< std::string > * >( "report_string", tag->getOption< std::string >( "report_string" )  );
	if( tag->hasOption( "filter" ) )
		filter( protocols::rosetta_scripts::parse_filter( tag->getOption< std::string >( "filter" ), filters ) );
}

bool
ReportFilter::apply( core::pose::Pose const & pose ) const {
	return compute( pose );
}

void
ReportFilter::report( std::ostream & out, core::pose::Pose const & /*pose*/ ) const {
	using namespace protocols::jd2;
	out<<"filter: "<<report_filter_name_;
	protocols::jd2::JobOP job2 = jd2::JobDistributor::get_instance()->current_job();
	std::string job_name (JobDistributor::get_instance()->job_outputter()->output_name( job2 ) );
	if( report_string_->obj.length() > 0 )
		out<<"job name: "<<job_name<<" report_string: "<<report_string_->obj<<std::endl;
	if( filter_() != NULL )
		out<<"job name: "<<job_name<<" reporting filter value: ";
}

core::Real
ReportFilter::report_sm( core::pose::Pose const & /*pose*/ ) const {
	return( 0.0 );
}

core::Real
ReportFilter::compute(
	core::pose::Pose const &
) const {
	return( 0.0 );
}

void ReportFilter::report_string( std::string const f )
{
	report_string_->obj = f;
}

std::string ReportFilter::report_string() const
{
	return report_string_->obj;
}

}
}
