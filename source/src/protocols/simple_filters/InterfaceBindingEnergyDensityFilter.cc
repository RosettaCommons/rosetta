// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/InterfaceBindingEnergyDensityFilter.cc
/// @brief  Implementation of the binding-energy-density filter
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/simple_filters/InterfaceBindingEnergyDensityFilter.hh>
#include <protocols/simple_filters/InterfaceBindingEnergyDensityFilterCreator.hh>

// Package headers
#include <protocols/simple_filters/InterfaceSasaFilter.hh>
#include <protocols/simple_filters/DdgFilter.hh>
#include <basic/datacache/DataMap.hh>

// Project headers

//#include <core/pose/Pose.hh>
//#include <core/id/AtomID_Map.hh>
//#include <core/conformation/Residue.hh>
//#include <core/chemical/AtomType.hh>
//#include <basic/MetricValue.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// Basic headers
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace simple_filters {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_filters.InterfaceBindingEnergyDensityFilter" );

protocols::filters::FilterOP
InterfaceBindingEnergyDensityFilterCreator::create_filter() const { return protocols::filters::FilterOP( new InterfaceBindingEnergyDensityFilter ); }

std::string
InterfaceBindingEnergyDensityFilterCreator::keyname() const { return "InterfaceBindingEnergyDensityFilter"; }


InterfaceBindingEnergyDensityFilter::InterfaceBindingEnergyDensityFilter() :
	Filter( "InterfaceBindingEnergyDensityFilter" ),
	sasa_filter_( /* 0 */ ),
	ddG_filter_( /* 0 */ ),
	upper_threshold_( 0.0 )
{}

InterfaceBindingEnergyDensityFilter::InterfaceBindingEnergyDensityFilter(
	InterfaceSasaFilterOP sasa_filter,
	DdgFilterOP ddG_filter,
	core::Real threshold
):
	Filter( "InterfaceBindingEnergyDensityFilter" ),
	sasa_filter_(std::move( sasa_filter )),
	ddG_filter_(std::move( ddG_filter )),
	upper_threshold_( threshold )
{}


void InterfaceBindingEnergyDensityFilter::set_interface_sasa_filter( InterfaceSasaFilterOP sasa_filter ) { sasa_filter_ = sasa_filter; }
void InterfaceBindingEnergyDensityFilter::set_ddG_filter( DdgFilterOP ddG_filter ) { ddG_filter_ = ddG_filter; }
void InterfaceBindingEnergyDensityFilter::set_upper_threshold( core::Real threshold ) { upper_threshold_ = threshold; }

InterfaceBindingEnergyDensityFilter::~InterfaceBindingEnergyDensityFilter()= default;

filters::FilterOP
InterfaceBindingEnergyDensityFilter::clone() const{
	return filters::FilterOP( new InterfaceBindingEnergyDensityFilter( *this ) );
}

filters::FilterOP
InterfaceBindingEnergyDensityFilter::fresh_instance() const{
	return filters::FilterOP( new InterfaceSasaFilter );
}

void
InterfaceBindingEnergyDensityFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	filters::Filters_map const & filters_map,
	moves::Movers_map const &,
	core::pose::Pose const &
)
{
	if ( ! tag->hasOption("sasa_filter") ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "InterfaceBindingEnergyDensityFilter requires the sasa_filter option" );
	}
	if ( ! tag->hasOption("ddG_filter") ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "InterfaceBindingEnergyDensityFilter requires the ddG_filter option" );
	}

	std::string sasa_filter_name = tag->getOption< std::string> ("sasa_filter");
	auto sasaiter = filters_map.find( sasa_filter_name );
	if ( sasaiter == filters_map.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Could not locate requested sasa_filter with name " + sasa_filter_name + " in the Filters_map." );
	}
	filters::FilterOP sasafilter_baseptr = sasaiter->second;
	sasa_filter_ = utility::pointer::dynamic_pointer_cast< protocols::simple_filters::InterfaceSasaFilter > ( sasafilter_baseptr );
	if ( ! sasa_filter_ ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Dynamic cast of filter " + sasa_filter_name + " to type InterfaceSasaFilter failed" );
	}

	std::string ddG_filter_name = tag->getOption< std::string> ("ddG_filter");
	auto ddGiter = filters_map.find( ddG_filter_name );
	if ( ddGiter == filters_map.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Could not locate requested ddG_filter with name " + ddG_filter_name + " in the Filters_map." );
	}
	filters::FilterOP ddGfilter_baseptr = ddGiter->second;
	ddG_filter_ = utility::pointer::dynamic_pointer_cast< protocols::simple_filters::DdgFilter > ( ddGfilter_baseptr );
	if ( ! ddG_filter_ ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Dynamic cast of filter " + ddG_filter_name + " to type DdgFilter failed" );
	}

	upper_threshold_ = tag->getOption< core::Real > ("threshold", -0.015 );

}

bool
InterfaceBindingEnergyDensityFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const sasa( sasa_filter_->compute( pose ) );
	core::Real const ddG( ddG_filter_->compute( pose ) );

	core::Real const binding_energy_density( ddG / sasa );

	TR << "SASA is "<<sasa<<". ddG is "<<ddG <<". Ratio is " << binding_energy_density << ". ";
	if ( binding_energy_density <= upper_threshold_ ) {
		TR<<"passing." <<std::endl;
		return true;
	} else {
		TR<<"failing."<<std::endl;
		return false;
	}
}

void
InterfaceBindingEnergyDensityFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const sasa( sasa_filter_->compute( pose ) );
	core::Real const ddG( ddG_filter_->compute( pose ) );
	core::Real const binding_energy_density( ddG / sasa );

	out << "SASA is "<<sasa<<". ddG is "<<ddG <<". Ratio is " << binding_energy_density << ". ";
}

core::Real
InterfaceBindingEnergyDensityFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const sasa( sasa_filter_->compute( pose ) );
	core::Real const ddG( ddG_filter_->compute( pose ) );
	core::Real const binding_energy_density( ddG / sasa );

	return binding_energy_density;
}

core::Real
InterfaceBindingEnergyDensityFilter::compute( core::pose::Pose const & pose ) const {
	core::Real const sasa( sasa_filter_->compute( pose ) );
	core::Real const ddG( ddG_filter_->compute( pose ) );
	core::Real const binding_energy_density( ddG / sasa );

	return binding_energy_density;
}

}
}
