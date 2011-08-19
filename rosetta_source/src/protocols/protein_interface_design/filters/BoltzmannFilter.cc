// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 sw=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Sarel Fleishman (sarelf@uw.edu)
#include <protocols/protein_interface_design/filters/BoltzmannFilter.hh>
#include <protocols/protein_interface_design/filters/BoltzmannFilterCreator.hh>
#include <utility/string_util.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/moves/DataMap.hh>
#include <basic/Tracer.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH
#include <math.h>

namespace protocols {
namespace protein_interface_design{
namespace filters {

static basic::Tracer TR( "protocols.protein_interface_design.filters.BoltzmannFilter" );

///@brief default ctor
BoltzmannFilter::BoltzmannFilter() :
	temperature_( 0.6 ),
	fitness_threshold_( 0.0 )
{
	positive_filters_.clear();
	negative_filters_.clear();
}


core::Real
BoltzmannFilter::fitness_threshold() const{
	return fitness_threshold_;
}

void
BoltzmannFilter::fitness_threshold( core::Real const f ){
	fitness_threshold_ = f;
}

core::Real
BoltzmannFilter::temperature() const{
	return temperature_;
}

void
BoltzmannFilter::temperature( core::Real const t ){
	temperature_ = t;
}

utility::vector1< protocols::filters::FilterOP >
BoltzmannFilter::get_positive_filters() const{
	return positive_filters_;
}

utility::vector1< protocols::filters::FilterOP >
BoltzmannFilter::get_negative_filters() const{
	return negative_filters_;
}

void
BoltzmannFilter::add_positive_filter( protocols::filters::FilterOP f ){
	positive_filters_.push_back( f );
}

void
BoltzmannFilter::add_negative_filter( protocols::filters::FilterOP f ){
	negative_filters_.push_back( f );
}

bool
BoltzmannFilter::apply(core::pose::Pose const & pose ) const
{
	core::Real const fitness( compute( pose ) );
	return( fitness <= fitness_threshold() );
}

/// NOTICE that this returns -Fitness [-1:0] for use in optimization
core::Real
BoltzmannFilter::compute( core::pose::Pose const & pose ) const{
	using protocols::filters::FilterCOP;

	core::Real positive_sum( 0.0 ), negative_sum( 0.0 );

	foreach( FilterCOP filter, get_positive_filters() )
		positive_sum += exp( -filter->report_sm( pose ) / temperature() );

	foreach( FilterCOP filter, get_negative_filters() )
		negative_sum += exp( -filter->report_sm( pose ) / temperature() );

	return( -positive_sum / ( positive_sum + negative_sum ));
}

core::Real
BoltzmannFilter::report_sm( core::pose::Pose const & pose ) const
{
	return( compute( pose ) );
}

void
BoltzmannFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
}

void
BoltzmannFilter::parse_my_tag( utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & pose )
{
	TR << "BoltzmannFilter"<<std::endl;
	fitness_threshold( tag->getOption< core::Real >( "fitness_threshold", 0 ) );
	temperature( tag->getOption< core::Real >( "temperature", 0.6 ) );
	utility::vector1< std::string > const positive_filter_names( utility::string_split( tag->getOption< std::string >( "positive_filters" ), ',' ) );
	utility::vector1< std::string > const negative_filter_names( utility::string_split( tag->getOption< std::string >( "negative_filters" ), ',' ) );
	foreach( std::string const positive_filter_name, positive_filter_names )
		add_positive_filter( protocols::rosetta_scripts::parse_filter( positive_filter_name, filters ) );
	foreach( std::string const negative_filter_name, negative_filter_names )
		add_negative_filter( protocols::rosetta_scripts::parse_filter( negative_filter_name, filters ) );

	TR<<"with options temperature: "<<temperature()<<" fitness_threshold "<<fitness_threshold()<<"  "<<get_positive_filters().size()<<" positive and "<<get_negative_filters().size()<<" negative filters."<<std::endl;
}

protocols::filters::FilterOP
BoltzmannFilter::fresh_instance() const{
	return new BoltzmannFilter();
}

BoltzmannFilter::~BoltzmannFilter(){}

protocols::filters::FilterOP
BoltzmannFilter::clone() const{
	return new BoltzmannFilter( *this );
}

protocols::filters::FilterOP
BoltzmannFilterCreator::create_filter() const { return new BoltzmannFilter; }

std::string
BoltzmannFilterCreator::keyname() const { return "Boltzmann"; }

} // filters
} // protein_interface_design
} // protocols
