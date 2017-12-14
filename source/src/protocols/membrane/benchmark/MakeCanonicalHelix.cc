// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/benchmark/MakeCanonicalHelix.cc
/// @brief Creates an ideal a-helix from sequence given a range of residues
/// @author Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <protocols/membrane/benchmark/MakeCanonicalHelix.hh>
#include <protocols/membrane/benchmark/MakeCanonicalHelixCreator.hh>

// Project Headers
#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.membrane.benchmark.MakeCanonicalHelix" );

namespace protocols {
namespace membrane {
namespace benchmark {

MakeCanonicalHelix::MakeCanonicalHelix():
	protocols::moves::Mover( "MakeCanonicalHelix" ),
	phi_( -57.0 ),
	psi_( -47.0 ),
	omega_( 175.0 ),
	helix_start_( 1 /* Start at the first residue */ ),
	helix_end_( 0 /* Set by mover or user */ )
{}

MakeCanonicalHelix::MakeCanonicalHelix( core::Size helix_start, core::Size helix_end ):
	protocols::moves::Mover( "MakeCanonicalHelix" ),
	phi_( -57.0 ),
	psi_( -47.0 ),
	omega_( 175.0 ),
	helix_start_( helix_start ),
	helix_end_( helix_end )
{}

MakeCanonicalHelix::~MakeCanonicalHelix()= default;

MakeCanonicalHelix::MakeCanonicalHelix( MakeCanonicalHelix const & /*src*/ ) = default;

void
MakeCanonicalHelix::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	// Read in helix start position
	if ( tag->hasOption( "helix_start" ) ) {
		helix_start_ = tag->getOption< core::Real >( "helix_start" );
	}

	// Read in helix end position
	if ( tag->hasOption( "helix_end" ) ) {
		helix_end_ = tag->getOption< core::Real >( "helix_end" );
	}
}

protocols::moves::MoverOP
MakeCanonicalHelix::clone() const{
	return protocols::moves::MoverOP( new MakeCanonicalHelix( *this ) );
}

moves::MoverOP
MakeCanonicalHelix::fresh_instance() const
{
	return protocols::moves::MoverOP( new MakeCanonicalHelix );
}

// XRW TEMP std::string
// XRW TEMP MakeCanonicalHelix::get_name() const {
// XRW TEMP  return "MakeCanonicalHelix";
// XRW TEMP }

void
MakeCanonicalHelix::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

void
MakeCanonicalHelix::apply( core::pose::Pose& pose )
{
	// Check start & end positions
	if ( helix_end_ == 0 ) helix_end_ = pose.size();
	is_valid( pose );

	// Override the dihedral angles according to the rules setup by this mover
	for ( core::Size i = helix_start_; i <= helix_end_; ++i ) {
		pose.set_phi( i, phi_ );
		pose.set_psi( i, psi_ );
		pose.set_omega( i, omega_ );
	}

}

void
MakeCanonicalHelix::is_valid( core::pose::Pose& pose ) {

	// Check the starting position
	if ( helix_start_ < 1 || helix_start_ > pose.size() ) {
		TR.Fatal << "Helix start set at " << helix_start_ << " is out of bounds for this pose" << std::endl;
		utility_exit();
	}

	// Check the end position
	if ( helix_end_ > pose.size() || helix_end_ < helix_start_ || helix_end_ < 0 ) {
		TR.Fatal << "helix end set at " << helix_end_ << " is out of bounds for this pose" << std::endl;
		utility_exit();
	}
}

/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP MakeCanonicalHelixCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new MakeCanonicalHelix );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP MakeCanonicalHelixCreator::keyname() const {
// XRW TEMP  return MakeCanonicalHelix::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP MakeCanonicalHelix::mover_name(){
// XRW TEMP  return "MakeCanonicalHelix";
// XRW TEMP }

std::string MakeCanonicalHelix::get_name() const {
	return mover_name();
}

std::string MakeCanonicalHelix::mover_name() {
	return "MakeCanonicalHelix";
}

void MakeCanonicalHelix::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute("helix_start", xsct_real, "Residue to start with; by default, residue 1")
		+ XMLSchemaAttribute("helix_end", xsct_real, "Residue to end with - user added");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Forms a canonical alpha helix between the specified residues", attlist );
}

std::string MakeCanonicalHelixCreator::keyname() const {
	return MakeCanonicalHelix::mover_name();
}

protocols::moves::MoverOP
MakeCanonicalHelixCreator::create_mover() const {
	return protocols::moves::MoverOP( new MakeCanonicalHelix );
}

void MakeCanonicalHelixCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MakeCanonicalHelix::provide_xml_schema( xsd );
}


} //protocols
} //membrane
} //benchmark


