// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/BridgeChainsMover.cc
/// @brief Creates a bridge connection between two chains using remodel
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/movers/BridgeChainsMover.hh>
#include <protocols/denovo_design/movers/BridgeChainsMoverCreator.hh>

// Protocol headers
#include <protocols/denovo_design/architects/CompoundArchitect.hh>
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/components/RandomTorsionPoseFolder.hh>
#include <protocols/denovo_design/components/RemodelLoopMoverPoseFolder.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>
#include <protocols/denovo_design/movers/BuildDeNovoBackboneMover.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.BridgeChainsMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

BridgeChainsMover::BridgeChainsMover():
	protocols::moves::Mover( BridgeChainsMover::mover_name() ),
	architect_( new connection::ConnectionArchitect( "BridgeChainsArchitect" ) ),
	scorefxn_(),
	overlap_( 1 ),
	iterations_( 0 ),
	dry_run_( false )
{
}

BridgeChainsMover::BridgeChainsMover( std::string const & mover_name ):
	protocols::moves::Mover( mover_name ),
	architect_( new connection::ConnectionArchitect( "BridgeChainsArchitect" ) ),
	scorefxn_(),
	overlap_( 1 ),
	iterations_( 0 ),
	dry_run_( false )
{
}

BridgeChainsMover::~BridgeChainsMover(){}

void
BridgeChainsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	// Build architect from tag data
	architect_->set_bridge( true );
	architect_->parse_my_tag( tag, data );
	set_overlap( tag->getOption< core::Size >( "overlap", overlap_ ) );
	set_dry_run( tag->getOption< bool >( "dry_run", dry_run_ ) );
	iterations_ = tag->getOption< core::Size >( "trials", iterations_ );

	core::scoring::ScoreFunctionCOP sfxn = protocols::rosetta_scripts::parse_score_function( tag, data );
	if ( sfxn ) set_scorefxn( *sfxn );
}

protocols::moves::MoverOP
BridgeChainsMover::clone() const
{
	return protocols::moves::MoverOP( new BridgeChainsMover( *this ) );
}

protocols::moves::MoverOP
BridgeChainsMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new BridgeChainsMover );
}

// XRW TEMP std::string
// XRW TEMP BridgeChainsMover::get_name() const
// XRW TEMP {
// XRW TEMP  return BridgeChainsMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP BridgeChainsMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "BridgeChainsMover";
// XRW TEMP }

void
BridgeChainsMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, BridgeChainsMover const & mover )
{
	mover.show(os);
	return os;
}

void
BridgeChainsMover::apply( core::pose::Pose & pose )
{
	architects::CompoundArchitect arch( "" );
	arch.add_architect( architects::PoseArchitect( "" ) );
	arch.add_connection( architect() );

	BuildDeNovoBackboneMover assemble;
	assemble.set_architect( arch );
	assemble.set_build_overlap( overlap() );
	assemble.set_iterations_per_phase( iterations_ );
	if ( dry_run() ) {
		assemble.set_folder( components::RandomTorsionPoseFolder() );
	} else {
		components::RemodelLoopMoverPoseFolder folder;
		if ( scorefxn_ ) {
			folder.set_scorefxn( scorefxn() );
		}
		assemble.set_folder( folder );
	}
	assemble.apply( pose );
}

connection::ConnectionArchitect const &
BridgeChainsMover::architect() const
{
	return *architect_;
}

void
BridgeChainsMover::set_id( std::string const & my_id )
{
	architect_->set_id( my_id );
}

void
BridgeChainsMover::set_segment1_ids( std::string const & segments )
{
	architect_->set_segment1_ids( segments );
}

void
BridgeChainsMover::set_segment2_ids( std::string const & segments )
{
	architect_->set_segment2_ids( segments );
}

void
BridgeChainsMover::set_motifs( std::string const & motifs_str, std::string const & cut_resis_str )
{
	architect_->set_motifs( motifs_str, cut_resis_str );
}

core::Size
BridgeChainsMover::overlap() const
{
	return overlap_;
}

void
BridgeChainsMover::set_overlap( core::Size const overlap_val )
{
	overlap_ = overlap_val;
}

core::scoring::ScoreFunction const &
BridgeChainsMover::scorefxn() const
{
	if ( !scorefxn_ ) {
		std::stringstream msg;
		msg << architect().id() << ": No score function set!" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return *scorefxn_;
}

void
BridgeChainsMover::set_scorefxn( core::scoring::ScoreFunction const & sfxn )
{
	scorefxn_ = sfxn.clone();
}

bool
BridgeChainsMover::dry_run() const
{
	return dry_run_;
}

void
BridgeChainsMover::set_dry_run( bool const dry_run )
{
	dry_run_ = dry_run;
}

/////////////// Creator ///////////////

std::string BridgeChainsMover::get_name() const {
	return mover_name();
}

std::string BridgeChainsMover::mover_name() {
	return "BridgeChains";
}

void BridgeChainsMover::setup_attlist_for_derived_classes( utility::tag::AttributeList & attlist ) {
	using namespace utility::tag;
	connection::ConnectionArchitect::attributes_for_parse_my_tag(attlist);

	attlist + XMLSchemaAttribute(
		"overlap", xsct_non_negative_integer,
		"Build overlap of nested BuildDeNovoBackboneMover");

	attlist + XMLSchemaAttribute(
		"dry_run", xsct_non_negative_integer,
		"Sets folder of BuildDeNovoBackboneMover to RandomTorsionPoseFolder "
		"instead of RemodelLoopMoverPoseFolde");

	attlist + XMLSchemaAttribute(
		"trials", xsct_non_negative_integer,
		"iterations per phase of nested BuildDeNovoBackboneMover");

	rosetta_scripts::attributes_for_parse_score_function(attlist);
}

void BridgeChainsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	BridgeChainsMover::setup_attlist_for_derived_classes( attlist );

	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Given a pose with a jump, this mover uses a fragment insertion monte carlo to connect "
		"the specified termini. The new fragment will connect the C-terminal residue of jump1 to "
		"the N-terminal residue of jump2, and will have secondary structure and ramachandran space "
		"given by \"motif.\" This mover uses the VarLengthBuild code. The input pose must have at least "
		"two chains (jumps) to connect, or it will fail.",
		attlist );
}

std::string BridgeChainsMoverCreator::keyname() const {
	return BridgeChainsMover::mover_name();
}

protocols::moves::MoverOP
BridgeChainsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new BridgeChainsMover );
}

void BridgeChainsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BridgeChainsMover::provide_xml_schema( xsd );
}


/// this one is to keep compatibility with Kuhlman, Jacobs, and Linsky 2016 book chapter
protocols::moves::MoverOP
BridgeChainsCreator::create_mover() const
{
	return protocols::moves::MoverOP( new BridgeChainsMover );
}

std::string
BridgeChainsCreator::keyname() const
{
	return "BridgeChains";
}

//std::string BridgeChainsMover::mover_name() {
// return "BridgeChains";
//}
//
//void BridgeChainsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
//{
//// TO DO!
// using namespace utility::tag;
// AttributeList attlist; // TO DO: add attributes to this list
// // TO DO: perhaps this is not the right function to call? -- also, delete this comment
// protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
//}
//
//std::string BridgeChainsCreator::keyname() const {
// return BridgeChainsMover::mover_name();
//}
//
//protocols::moves::MoverOP
//BridgeChainsCreator::create_mover() const {
// return protocols::moves::MoverOP( new BridgeChainsMover );
//}

void BridgeChainsCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BridgeChainsMover::provide_xml_schema( xsd );
}


} //protocols
} //denovo_design
} //movers
