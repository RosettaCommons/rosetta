// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/Filter.cc
/// @brief
/// @detailed
///	  Contains currently:
///
///
/// @author Florian Richter, Sarel Fleishman (sarelf@uw.edu)

// Unit Headers
#include <protocols/filters/Filter.hh>

// Project Headers
// AUTO-REMOVED #include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/pose/util.hh>


static thread_local basic::Tracer TR( "protocols.filters.Filter" );

namespace protocols {
namespace filters {

#ifdef USELUA
void lregister_Filter( lua_State * lstate ) {
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("filters")
		[
			luabind::class_<Filter>("Filter")
				.def("apply", ( void (Filter::*)( core::io::serialization::PipeMap & )) &Filter::apply)
				.def("score", ( void (Filter::*)( core::io::serialization::PipeMap & )) &Filter::score)
				.def("score", ( core::Real (Filter::*)( core::pose::Pose & )) &Filter::score)
		]
	];
}
#endif

using namespace core;
typedef std::pair< std::string const, FilterCOP > StringFilter_pair;
typedef utility::tag::TagCOP TagCOP;

FilterCollection::~FilterCollection() {}

bool
FilterCollection::apply( core::pose::Pose const & pose ) const
{
	BOOST_FOREACH(protocols::filters::FilterCOP filter, filters_){
		if( ! filter->apply( pose ) ){
			return false;
		}
	}

	return true;
}

void
FilterCollection::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	BOOST_FOREACH(protocols::filters::FilterCOP filter, filters_){
		filter->report( out, pose );
	}
}

Filter::Filter()
	: utility::pointer::ReferenceCount(),
		type_( "UNDEFINED TYPE" ),
		scorename_("defaultscorename")
{}

Filter::Filter( std::string const & type )
	: utility::pointer::ReferenceCount(),
		type_( type ),
		scorename_("defaultscorename")
{}

Filter::Filter( Filter const & init )
	:	utility::pointer::ReferenceCount(),
		type_( init.type_ ),
		user_defined_name_( init.user_defined_name_ ),
		scorename_("defaultscorename")
		
{}

Filter::~Filter() {}

void
Filter::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const & )
{}

// start mpr support
void Filter::apply( core::io::serialization::PipeMap & pmap ) {
	core::io::serialization::Pipe::iterator itr = pmap["input"]->begin();
	while( itr != pmap["input"]->end() ) {
		if( !apply ( **itr ) ) {
			itr = pmap["input"]->erase( itr );	
		} else {
			itr++;
		}
		clear();
	}
}
core::Real Filter::score( core::pose::Pose & pose ) {
	core::Real score = report_sm( pose );
	core::pose::setPoseExtraScore( pose, scorename_, score );
	return score;
}
void Filter::score( core::io::serialization::PipeMap & pmap ) {
	for( core::io::serialization::Pipe::iterator itr = pmap["input"]->begin(); itr != pmap["input"]->end(); itr++ ) {
		score( **itr );
		clear();
	}
}
void Filter::parse_def( utility::lua::LuaObject const & /*def*/,
				utility::lua::LuaObject const & /*score_fxns*/,
				utility::lua::LuaObject const & /*tasks*/ ){
	utility_exit_with_message("This Filter has not implemented parse_def()");
}
// end mpr support

} // filters
} // protocols
