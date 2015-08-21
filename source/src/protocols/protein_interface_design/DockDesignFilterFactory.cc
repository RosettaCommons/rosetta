// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/protein_interface_design/DockDesignFilterFactory.cc
/// @brief
/// @author ashworth
#include <protocols/protein_interface_design/DockDesignFilterFactory.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/filters/BasicFilters.hh>


#include <utility/exit.hh> // runtime_assert, utility_exit_with_message
#include <utility/tag/Tag.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace protein_interface_design {
using namespace protocols::filters;

DockDesignFilterFactory::DockDesignFilterFactory()
: utility::pointer::ReferenceCount()
{
	// no Filters are registered by default
	// they must be registered using the add_type method
}

DockDesignFilterFactory::~DockDesignFilterFactory(){}

/// @brief add a Filter prototype, using its default type name as the map key
void
DockDesignFilterFactory::add_type( FilterOP dock_design_filter )
{
	runtime_assert( dock_design_filter != 0 );
	std::string const type( dock_design_filter->get_type() );
	if ( type == "UNDEFINED TYPE" ) {
		utility_exit_with_message("Can't map derived Filter with undefined type name.");
	}
	dock_design_filter_map_[ type ] = dock_design_filter;
}

/// @brief add a Filter prototype, using an arbitrary type name as the map key
void
DockDesignFilterFactory::add_type( std::string const & type, FilterOP dock_design_filter )
{
	runtime_assert( dock_design_filter != 0 );
	dock_design_filter_map_[ type ] = dock_design_filter;
}

/// @brief return new Filter by key lookup in dock_design_filter_map_ (new Filter parses Tag if provided)
FilterOP
DockDesignFilterFactory::newFilter( std::string const & type )
{
	filters::Filters_map::const_iterator iter( dock_design_filter_map_.find( type ) );
	if ( iter != dock_design_filter_map_.end() ) {
		if ( ! iter->second ) {
			utility_exit_with_message( "Error: FilterOP prototype for " + type + " is NULL!" );
		}
		return iter->second->fresh_instance();
	} else {
		utility_exit_with_message( type + " is not known to the FilterFactory. Was it registered via the add_type method?" );
		return NULL;
	}
}

/// @brief return new Filter by Tag parsing
FilterOP
DockDesignFilterFactory::newFilter(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & pose )
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

} //namespace protein_interface_design
} //namespace protocols
