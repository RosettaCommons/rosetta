// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/filters/FilterFactory.cc
/// @brief
/// @author ashworth

#include <protocols/filters/FilterFactory.hh>
#include <protocols/filters/Filter.hh>

// Package headers
#include <protocols/filters/BasicFilters.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace protocols {
namespace filters {


static THREAD_LOCAL basic::Tracer TR( "protocols.filters.FilterFactory" );

#if defined MULTI_THREADED && defined CXX11
std::atomic< FilterFactory * > FilterFactory::instance_( 0 );
#else
FilterFactory * FilterFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex FilterFactory::singleton_mutex_;

std::mutex & FilterFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
FilterFactory * FilterFactory::get_instance()
{
	boost::function< FilterFactory * () > creator = boost::bind( &FilterFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

FilterFactory *
FilterFactory::create_singleton_instance()
{
	return new FilterFactory;
}


FilterFactory::FilterFactory()
{}

FilterFactory::~FilterFactory(){}

/// @brief add a Filter prototype, using its default type name as the map key
void
FilterFactory::factory_register( FilterCreatorOP creator )
{
	runtime_assert( creator != 0 );
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
		return NULL;
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
	runtime_assert( filter != 0 );
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


} //namespace filters
} //namespace protocols
