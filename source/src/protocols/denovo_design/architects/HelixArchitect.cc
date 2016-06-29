// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/architects/HelixArchitect.cc
/// @brief Architect for helices
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/architects/HelixArchitect.hh>
#include <protocols/denovo_design/architects/HelixArchitectCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.architects.HelixArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

HelixArchitect::HelixArchitect( std::string const & id_value ):
	DeNovoArchitect( id_value ),
	lengths_()
{
}

HelixArchitect::~HelixArchitect()
{}

DeNovoArchitectOP
HelixArchitect::clone() const
{
	return DeNovoArchitectOP( new HelixArchitect( *this ) );
}

std::string
HelixArchitect::type() const
{
	return HelixArchitect::class_name();
}

void
HelixArchitect::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	set_lengths( tag->getOption< std::string >( "length" ) );
}

StructureDataOP
HelixArchitect::design( core::pose::Pose const &, core::Real & random ) const
{
	if ( lengths_.empty() ) {
		utility_exit_with_message( type() + "Architect requires one or more motifs or lengths to be specified" );
	}
	core::Size const length_idx = extract_int( random, 1, lengths_.size() );

	std::stringstream ss;
	ss << 'L' << std::string( lengths_[ length_idx ], 'H' ) << 'L';
	std::stringstream abego;
	abego << 'X' << std::string( lengths_[ length_idx ], 'A' ) << 'X';

	components::Segment motif;
	motif.extend( ss.str(), abego.str() );

	TR << "Created new segment: " << motif << std::endl;

	StructureDataOP sd( new StructureData( id() ) );
	sd->add_segment( id(), motif );

	return sd;
}

void
HelixArchitect::set_lengths( std::string const & lengths_str )
{
	set_lengths( parse_length_str< core::Size >( lengths_str ) );
}

void
HelixArchitect::set_lengths( Lengths const & lengths )
{
	lengths_ = lengths;
}

///////////////////////////////////////////////////////////////////////////////
/// Creator
///////////////////////////////////////////////////////////////////////////////
DeNovoArchitectOP
HelixArchitectCreator::create_architect( std::string const & architect_id ) const
{
	return DeNovoArchitectOP( new HelixArchitect( architect_id ) );
}

std::string
HelixArchitectCreator::keyname() const
{
	return HelixArchitect::class_name();
}

} //protocols
} //denovo_design
} //architects
