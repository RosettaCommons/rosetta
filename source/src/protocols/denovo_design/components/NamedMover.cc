// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/NamedMover.cc
/// @brief Named segment of residues
/// @detailed
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/components/NamedMover.hh>

//Project Headers
#include <protocols/denovo_design/util.hh>

//Protocols Headers

//Core Headers

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

//Boost
#include <boost/algorithm/string.hpp>
//C++ Headers

static basic::Tracer TR("protocols.denovo_design.components.NamedMover");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace components {

NamedMover::NamedMover() :
	protocols::moves::Mover(),
	id_( "" ),
	parent_id_( "" )
{
}

NamedMover::NamedMover( std::string const & idval, std::string const & parent_id ):
	protocols::moves::Mover(),
	id_( idval ),
	parent_id_( parent_id )
{
}

/// @brief setup the parameters via an xml tag
void
NamedMover::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const & )
{
	if ( tag->hasOption( "name" ) ) {
		set_id( tag->getOption< std::string >( "name" ) );
	} else {
		throw utility::excn::EXCN_RosettaScriptsOption( "Name is required for " + tag->getName() );
	}
	if ( tag->hasOption( "parent" ) ) {
		set_parent_id( tag->getOption< std::string >( "parent" ) );
	}
}

/// @brief function that creates a name for mover data to be stored in SegmentData
std::string
NamedMover::data_name( std::string const & name_of_data ) const
{
	debug_assert( id() != "" );
	debug_assert( name_of_data != "" );
	return id() + DATA_DELIMETER + name_of_data;
}

/// @brief adds prefix if necessary, returns result
std::string
NamedMover::add_parent_prefix( std::string const & s ) const
{
	debug_assert( s != "" );
	std::string prefix = parent_id() + PARENT_DELIMETER;
	if ( !parent_id().empty() && !boost::starts_with(s, prefix) ) {
		return prefix + s;
	} else {
		return s;
	}
}

std::string const &
NamedMover::id() const
{
	return id_;
}

std::string const &
NamedMover::parent_id() const
{
	return parent_id_;
}

void
NamedMover::set_id( std::string const & idval )
{
	id_ = add_parent_prefix( idval );
}

void
NamedMover::set_parent_id( std::string const & parentid )
{
	parent_id_ = parentid;
	if ( ! id_.empty() )
		set_id( id() );
}

} // namespace components
} // namespace denovo_design
} // namespace protocols
