// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/filters/FilterFactory.cc
/// @brief
/// @author ashworth

#include <protocols/filters/FilterFactory.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/filter_schemas.hh>

// Package headers
#include <protocols/filters/BasicFilters.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace filters {


static THREAD_LOCAL basic::Tracer TR( "protocols.filters.FilterFactory" );

FilterFactory::FilterFactory()
{}

FilterFactory::~FilterFactory()= default;

/// @brief add a Filter prototype, using its default type name as the map key
void
FilterFactory::factory_register( FilterCreatorOP creator )
{
	runtime_assert( creator != nullptr );
	std::string const filter_type( creator->keyname() );
	if ( filter_type == "UNDEFINED NAME" ) {
		utility_exit_with_message("Can't map derived Filter with undefined type name.");
	}
	if ( filter_creator_map_.find( filter_type ) != filter_creator_map_.end() ) {
		utility_exit_with_message("FilterFactory::factory_register already has a filter creator with name \"" + filter_type + "\".  Conflicting Filter names" );
	}
	filter_creator_map_[ filter_type ] = creator;
}


/// @brief return new Filter by key lookup in filter_prototype_map_ (new Filter parses Tag if provided)
FilterOP
FilterFactory::newFilter( std::string const & filter_type )
{
	FilterMap::const_iterator iter( filter_creator_map_.find( filter_type ) );
	if ( iter != filter_creator_map_.end() ) {
		if ( ! iter->second ) {
			utility_exit_with_message( "Error: FilterCreatorOP prototype for " + filter_type + " is NULL!" );
		}
		// use of cloning method would be faithful to pre-initialized prototypes
		//return iter->second->clone();
		// fresh_instance prevents propagation of pre-initialized prototypes, which may be safer(?)
		return iter->second->create_filter();
	} else {
		TR<<"Available filters: ";
		for ( FilterMap::const_iterator filt_it = filter_creator_map_.begin(); filt_it != filter_creator_map_.end(); ++filt_it ) {
			TR<<filt_it->first<<", ";
		}
		TR<<std::endl;
		utility_exit_with_message( filter_type + " is not known to the FilterFactory. Was it registered via a FilterRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return nullptr;
	}
}

/// @brief return new Filter by Tag parsing
/*FilterOP
FilterFactory::newFilter(
TagCOP const tag,
basic::datacache::DataMap & data,
Filters_map const & filters,
moves::Movers_map const & movers,
Pose const & pose )
{
FilterOP filter( newFilter( tag->getName() ) );
runtime_assert( filter );
filter->parse_my_tag( tag, data, filters, movers, pose );
return filter;
}*/

/// @brief return new Filter by Tag parsing
FilterOP
FilterFactory::newFilter(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const & filters,
	moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	FilterOP filter( newFilter( tag->getName() ) );
	runtime_assert( filter != nullptr );
	if ( ! tag->hasOption("name") ) {
		utility_exit_with_message("Can't define unnamed Filter");
	}
	filter->set_user_defined_name( tag->getOption<std::string>("name") );
	filter->parse_my_tag( tag, data, filters, movers, pose );
	// if confidence specified, link to StochasticFilter and wrap inside CompoundFilter
	core::Real const confidence( tag->getOption< core::Real >( "confidence", 1.0 ) );
	if ( confidence < 0.999 ) { // fuzzy logic
		CompoundFilter::CompoundStatement fuzzy_statement;
		FilterCOP stochastic_filter( FilterOP( new StochasticFilter( confidence ) ) );
		fuzzy_statement.push_back( std::make_pair( stochastic_filter->clone(), OR ) );
		fuzzy_statement.push_back( std::make_pair( filter->clone(), OR ) );
		FilterOP compound_filter( new CompoundFilter( fuzzy_statement ) );
		compound_filter->set_user_defined_name( tag->getOption<std::string>("name") );
		return compound_filter;
	}
	return filter;
}


FilterFactory::FilterMap const & FilterFactory::filter_creator_map() const { return filter_creator_map_; }

void FilterFactory::define_filter_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	try {
		utility::tag::define_xml_schema_group(
			filter_creator_map_,
			filter_xml_schema_group_name(),
			& complex_type_name_for_filter,
			xsd );
	} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
		throw utility::excn::EXCN_Msg_Exception( "Could not generate an XML Schema for Filters from FilterFactory; offending class"
			" must call protocols::filters::complex_type_name_for_filter when defining"
			" its XML Schema\n" + e.msg() );
	}

}

std::string FilterFactory::filter_xml_schema_group_name()
{
	return "filter";
}

} //namespace filters
} //namespace protocols
