// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/denovo_design/architects/StructureArchitect.cc
/// @brief Designs topologies
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/architects/StructureArchitect.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "protocols.denovo_design.architects.StructureArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

// Defined to prevent pure virtual destructor error at run time.
StructureArchitect::~StructureArchitect()
{}

StructureArchitect::StructureArchitect( std::string const & id ):
	utility::pointer::ReferenceCount(),
	id_( id )
{}

/// @brief private constructor -- should never be called
StructureArchitect::StructureArchitect()
{}

void
StructureArchitect::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	id_ = tag->getOption< std::string >( "name" );
	parse_tag( tag, data );
}

void
StructureArchitect::attributes_for_parse_my_tag(utility::tag::AttributeList& attlist) {
	using namespace utility::tag;

	attlist + XMLSchemaAttribute::required_attribute(
		"name", xs_string,
		"name of the architect");
}

std::string const &
StructureArchitect::id() const
{
	return id_;
}

void
StructureArchitect::set_id( std::string const & new_id )
{
	id_ = new_id;
}

} //protocols
} //denovo_design
} //architects
