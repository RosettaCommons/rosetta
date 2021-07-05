// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/movers/ExtendChainMover.cc
/// @brief Extends a chain by a structural motif
/// @details
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/movers/ExtendChainMover.hh>
#include <protocols/denovo_design/movers/ExtendChainMoverCreator.hh>

//Project Headers
#include <protocols/denovo_design/architects/CompoundArchitect.hh>
#include <protocols/denovo_design/architects/MotifArchitect.hh>
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/components/RandomTorsionPoseFolder.hh>
#include <protocols/denovo_design/components/RemodelLoopMoverPoseFolder.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>
#include <protocols/denovo_design/movers/FoldArchitectMover.hh>

//Protocol Headers

//Core Headers
#include <core/pose/Pose.fwd.hh>

//Basic Headers
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/string_util.hh>
#include <utility/stream_util.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

//ObjexxFCL Headers

//C++ Headers

static basic::Tracer TR( "protocols.denovo_design.movers.ExtendChainMover" );

///////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace movers {
///////////////////////////////////////////////////////////////////////////////



///  ---------------------------------------------------------------------------------
///  ExtendChainMover main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
ExtendChainMover::ExtendChainMover() :
	protocols::moves::Mover( ExtendChainMover::mover_name() ),
	architect_( utility::pointer::make_shared< architects::MotifArchitect >( "ExtendChain_Motif__" ) ),
	segment_names_(),
	chain_(),
	prepend_( false ),
	dry_run_( false )
{
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
ExtendChainMover::~ExtendChainMover() = default;


/// Return a copy of ourselves
protocols::moves::MoverOP
ExtendChainMover::clone() const
{
	return utility::pointer::make_shared< ExtendChainMover >( *this );
}

/// @brief return a fresh instance of ourselves
protocols::moves::MoverOP
ExtendChainMover::fresh_instance() const
{
	return utility::pointer::make_shared< ExtendChainMover >();
}

/// @brief return a fresh instance of ourselves

void
ExtendChainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
)
{
	architect_->parse_my_tag( tag, data );
	// scorefunction
	// overlap

	if ( tag->hasOption( "length" ) ) {
		std::stringstream msg;
		msg << architect().id() << ": The length option is not valid for ExtendChainMover -- please specify a motif using the \"motif\" option!";
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}
	if ( !tag->hasOption( "motif" ) ) {
		std::stringstream msg;
		msg << architect().id() << ": You must specify a motif (e.g. 1LG-1LB-10HA) to ExtendChainMover.";
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}

	std::string const & segment_id_str = tag->getOption< std::string >( "segment", "" );
	if ( !segment_id_str.empty() ) set_segment_names( segment_id_str );

	set_chain( tag->getOption< core::Size >( "chain", chain_ ) );

	set_prepend( tag->getOption< bool >( "prepend", prepend_ ) );

	if ( !segment_names_.empty() && chain_ ) {
		std::stringstream msg;
		msg << architect().id() << ": You cannot set both segment and chain in ExtendChainMover. Please specify one or the other.";
		msg << "segment1_ids = " << segment_names_ << std::endl;
		msg << "chain1 = " << chain_ << std::endl;
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}
	if ( segment_names_.empty() && !chain_ ) {
		std::stringstream msg;
		msg << architect().id() << ": You must set either segment or chain in ExtendChainMover, but you haven't specified either. Please specify either one or the other.";
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  msg.str() );
	}
}

void
ExtendChainMover::apply( core::pose::Pose & pose )
{
	architects::CompoundArchitect arch( "" );
	arch.add_architect( architects::PoseArchitect( "pose" ) );
	arch.add_architect( architect() );

	connection::ConnectionArchitect conn( "staple_extendchain__" );
	if ( prepend_ ) {
		conn.set_segment1_ids( architect().id() );
		conn.set_segment2_ids( "" );
	} else {
		conn.set_segment1_ids( "" );
		conn.set_segment2_ids( architect().id() );
	}
	arch.add_connection( conn );

	FoldArchitectMover assemble;
	assemble.set_architect( arch );

	assemble.set_build_overlap( 1 );
	if ( dry_run_ ) {
		assemble.set_folder( components::RandomTorsionPoseFolder() );
	} else {
		components::RemodelLoopMoverPoseFolder folder;
		assemble.set_folder( folder );
	}
	assemble.apply( pose );
}

architects::DeNovoArchitect const &
ExtendChainMover::architect() const
{
	return *architect_;
}

void
ExtendChainMover::set_segment_names( std::string const & segment_names_str )
{
	utility::vector1< std::string > const names = utility::string_split( segment_names_str, ',' );
	set_segment_names( SegmentNames( names.begin(), names.end() ) );
}

void
ExtendChainMover::set_segment_names( SegmentNames const & seg_names )
{
	segment_names_ = seg_names;
}

void
ExtendChainMover::set_chain( core::Size const chain )
{
	chain_ = chain;
}

void
ExtendChainMover::set_prepend( bool const prepend )
{
	prepend_ = prepend;
}

void
ExtendChainMover::set_dry_run( bool const dry_run )
{
	dry_run_ = dry_run;
}

std::string ExtendChainMover::get_name() const {
	return mover_name();
}

std::string ExtendChainMover::mover_name() {
	return "ExtendChain";
}

void ExtendChainMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// TO DO!
	using namespace utility::tag;
	AttributeList attlist; // TO DO: add attributes to this list

	attlist + XMLSchemaAttribute::required_attribute( "motif", xs_string, "Motif to add to chain (e.g. 1LG-1LB-10HA)")
		+ XMLSchemaAttribute( "segment", xs_string, "Segment to extend. Mutually exclusive with chain, but one is required.")
		+ XMLSchemaAttribute( "chain", xsct_non_negative_integer, "Chain to extend. Mutually exclusive with segment, but one is required.")
		+ XMLSchemaAttribute( "prepend", xsct_rosetta_bool, "Add motif to the beginning of the chain");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add a specified motif (with secondary structure and abego) to a chain", attlist );
}

std::string ExtendChainMoverCreator::keyname() const {
	return ExtendChainMover::mover_name();
}

protocols::moves::MoverOP
ExtendChainMoverCreator::create_mover() const {
	return utility::pointer::make_shared< ExtendChainMover >();
}

void ExtendChainMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ExtendChainMover::provide_xml_schema( xsd );
}


} // namespace connection
} // namespace denovo_design
} // namespace protocols
