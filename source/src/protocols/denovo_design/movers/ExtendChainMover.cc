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
#include <protocols/denovo_design/architects/PoseArchitect.hh>
#include <protocols/denovo_design/components/RandomTorsionPoseFolder.hh>
#include <protocols/denovo_design/components/RemodelLoopMoverPoseFolder.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>
#include <protocols/denovo_design/movers/BuildDeNovoBackboneMover.hh>

//Protocol Headers

//Core Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>

//Basic Headers
#include <basic/Tracer.hh>

//Utility Headers
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

// XRW TEMP std::string
// XRW TEMP ExtendChainMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ExtendChainMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ExtendChainMoverCreator::create_mover() const
// XRW TEMP {
// XRW TEMP  return protocols::moves::MoverOP( new ExtendChainMover() );
// XRW TEMP }

///  ---------------------------------------------------------------------------------
///  ExtendChainMover main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
ExtendChainMover::ExtendChainMover() :
	protocols::moves::Mover( ExtendChainMover::mover_name() ),
	architect_( new architects::DeNovoMotifArchitect( "ExtendChain_DeNovoMotif__" ) ),
	segment_names_(),
	chain_(),
	prepend_( false ),
	dry_run_( false )
{
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
ExtendChainMover::~ExtendChainMover() {}


/// Return a copy of ourselves
protocols::moves::MoverOP
ExtendChainMover::clone() const
{
	return protocols::moves::MoverOP( new ExtendChainMover( *this ) );
}

/// @brief return a fresh instance of ourselves
protocols::moves::MoverOP
ExtendChainMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new ExtendChainMover() );
}

/// @brief return a fresh instance of ourselves
// XRW TEMP std::string
// XRW TEMP ExtendChainMover::get_name() const
// XRW TEMP {
// XRW TEMP  return ExtendChainMover::mover_name();
// XRW TEMP }

void
ExtendChainMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
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

	BuildDeNovoBackboneMover assemble;
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
	return protocols::moves::MoverOP( new ExtendChainMover );
}

void ExtendChainMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ExtendChainMover::provide_xml_schema( xsd );
}


/*
/// @brief configures based on a permutation
/// @throw EXCN_Setup if no valid connection endpoints are found
void
ExtendChainMover::setup_permutation( components::StructureData & perm ) const
{
core::Real rand = numeric::random::rg().uniform();
core::Size const motifidx = extract_int( rand, 1, motifs().size() );
components::Motif const & motif = motifs()[ motifidx ];

std::string seg1 = "";
std::string seg2 = "";

if ( comp1_ids().size() ) {
utility::vector1< std::string > const local_comp1_ids = find_available_upper_termini( perm );
core::Size const s1_idx = extract_int( rand, 1, local_comp1_ids.size() );
debug_assert( seg1.empty() );
seg1 = local_comp1_ids[ s1_idx ];
}
if ( comp2_ids().size() ) {
utility::vector1< std::string > const local_comp2_ids = find_available_lower_termini( perm );
core::Size const s2_idx = extract_int( rand, 1, local_comp2_ids.size() );
debug_assert( seg2.empty() );
seg2 = local_comp2_ids[ s2_idx ];
}
if ( user_chain1() ) {
if ( !seg1.empty() ) {
throw CREATE_EXCEPTION(utility::excn::BadInput,  id() + ": you can only specify one segment/chain to ExtendChainMover!  You have specified more than one." );
}
utility::vector1< std::string > const local_comp1_ids = find_available_upper_termini( perm );
if ( user_chain1() > local_comp1_ids.size() ) {
std::stringstream msg;
msg << id() << ": The user-specified chain number (" << user_chain1()
<< ") is larger than the number of upper termini in the pose (" << local_comp1_ids.size()
<< "). Perm = " << perm << std::endl;
throw CREATE_EXCEPTION(utility::excn::BadInput,  msg.str() );
}
seg1 = local_comp1_ids[ user_chain1() ];
}
if ( user_chain2() ) {
if ( !seg2.empty() ) {
throw CREATE_EXCEPTION(utility::excn::BadInput,  id() + ": you can only specify one segment/chain to ExtendChainMover!  You have specified more than one." );
}
utility::vector1< std::string > const local_comp2_ids = find_available_lower_termini( perm );
if ( user_chain2() > local_comp2_ids.size() ) {
std::stringstream msg;
msg << id() << ": The user-specified chain2 number (" << user_chain2()
<< ") is larger than the number of lower termini in the pose (" << local_comp2_ids.size()
<< "). Perm = " << perm << std::endl;
throw CREATE_EXCEPTION(utility::excn::BadInput,  msg.str() );
}
seg2 = local_comp2_ids[ user_chain2() ];
}

if ( ( seg1.empty() && seg2.empty() ) || ( !seg1.empty() && !seg2.empty() ) ) {
std::stringstream msg;
msg << id() << ": you must specify exactly one valid segment/chain to ExtendChainMover!" << std::endl;
msg << "user_chain1() = " << user_chain1() << " user_chain2() = " << user_chain2() << std::endl;
msg << "segment1 = " << comp1_ids() << " segment2 = " << comp2_ids() << std::endl;
msg << "Perm = " << perm << std::endl;
throw CREATE_EXCEPTION(utility::excn::BadInput,  msg.str() );
}

std::string const segname = id();
std::stringstream complete_ss;
complete_ss << 'L' << motif.ss() << 'L';
std::stringstream complete_abego;
complete_abego << 'X' << motif.abego() << 'X';

components::StructureDataOP loop_sd( new components::StructureData( segname ) );
// WAS nterm + cterm included
components::Segment const seg( complete_ss.str(), abego_vector( complete_abego.str() ), 1, false, false );
loop_sd->add_segment( segname, seg );

setup_structuredata( *loop_sd,
motif.length(),
motif.ss(),
motif.abego(),
seg1,
seg2,
"",
"",
0 );

if ( !motif.length() ) {
TR.Warning << "length of motif given to ExtendChainMover is 0.  Doing nothing." << std::endl;
return;
}

perm.merge( *loop_sd );

if ( !seg1.empty() ) {
set_loop_upper( perm, "" );
set_loop_lower( perm, segname );
} else {
set_loop_upper( perm, segname );
set_loop_lower( perm, "" );
}

if ( !lower_segment_id( perm ).empty() ) {
debug_assert( perm.has_free_upper_terminus( lower_segment_id( perm ) ) );
}
if ( !upper_segment_id( perm ).empty() ) {
debug_assert( perm.has_free_lower_terminus( upper_segment_id( perm ) ) );
}

TR << "Going to extend " << lower_segment_id( perm ) << "-> "
<< loop_lower( perm ) << "-> " << loop_upper( perm )
<< "-> " << upper_segment_id( perm )
<< " len=" << build_len( perm ) << " cut=" << cut_resi( perm )
<< " ss=" << build_ss( perm ) << " abego=" << build_abego( perm ) << std::endl;

if ( !seg1.empty() ) {
perm.mark_connected( lower_segment_id( perm ), loop_lower( perm ) );
} else {
perm.mark_connected( loop_upper( perm ), upper_segment_id( perm ) );
}

}

void
ExtendChainMover::process_permutation( components::StructureDataWithPose & perm ) const
{
if ( !lower_segment_id( perm ).empty() ) {
debug_assert( !loop_lower( perm ).empty() );
perm.move_segment( lower_segment_id( perm ), lower_segment_id( perm ), loop_lower( perm ), loop_lower( perm ) );
} else if ( !upper_segment_id( perm ).empty() ) {
debug_assert( !loop_upper( perm ).empty() );
perm.move_segment( loop_upper( perm ), loop_upper( perm ), upper_segment_id( perm ), upper_segment_id( perm ) );
} else {
std::stringstream msg;
msg << id() << ": Lower segment id and upper segment id are both empty! Perm= " << perm << std::endl;
throw CREATE_EXCEPTION(utility::excn::BadInput,  msg.str() );
}

perm.chains_from_termini();

// do the work
apply_connection( perm );

for ( core::Size i=1; i<=perm.pose().size(); ++i ) {
TR << i << " : " << perm.pose().residue( i ).name() << " " << perm.pose().chain( i ) << std::endl;
}
}

/// @brief Does the work of remodeling the connection
void
ExtendChainMover::apply_connection( components::StructureDataWithPose & perm ) const
{
if ( !loop_lower( perm ).empty() ) {
connect_lower_loop( perm );
} else if ( !loop_upper( perm ).empty() ) {
connect_upper_loop( perm );
} else {
std::stringstream msg;
msg << id() << ": Neither loop_upper (" << loop_upper( perm ) << ") nor loop_lower("
<< loop_lower( perm ) << ") are the segment built by this connection. Perm=" << perm << std::endl;
throw CREATE_EXCEPTION(utility::excn::BadInput,  msg.str() );
}

// check sfxn
if ( !scorefxn() ) {
std::stringstream err;
err << "ExtendChain: You must set a valid scorefunction to "
<< id() << " before connecting" << std::endl;
throw CREATE_EXCEPTION(utility::excn::Exception,  err.str() );
}

// create loops
protocols::loops::LoopsOP loops( new protocols::loops::Loops() );
core::Size left, right;
if ( !loop_lower( perm ).empty() ) {
left = build_left( perm );
right = perm.segment( loop_lower( perm ) ).cterm_resi();
} else {
left = perm.segment( loop_upper( perm ) ).nterm_resi();
right = build_right( perm );
}
loops->add_loop( left, right, 0 );

protocols::moves::MoverOP remodel =
create_remodel_mover(
perm.pose(),
loops,
true,
perm.ss(),
abego_vector( perm.abego() ),
left,
right );

// store original in case there is a failure
components::StructureDataWithPoseOP orig = perm.attach_pose( perm.pose() );

// switch to centroid if necessary
bool const input_centroid = perm.pose().is_centroid();
core::pose::Pose stored;

if ( !input_centroid ) {
stored = perm.pose();
perm.switch_residue_type_set( "centroid" );
}

if ( do_remodel() ) {
apply_constraints( perm );
perm.apply_mover( remodel );
remove_constraints( perm );
}

if ( remodel->get_last_move_status() != protocols::moves::MS_SUCCESS ) {
throw EXCN_ConnectionFailed( "Failed to close loop during remodel in " + id() );
}

// switch back to FA if necessary
if ( !input_centroid ) {
perm.switch_residue_type_set( "fa_standard" );
copy_rotamers( perm, stored );
}

// Rebuild fold tree
utility::vector1< std::string > roots;
if ( !lower_segment_id( perm ).empty() ) {
roots.push_back( lower_segment_id( perm ) );
}
if ( !upper_segment_id( perm ).empty() ) {
roots.push_back( upper_segment_id( perm ) );
}
perm.consolidate_movable_groups( roots );
TR << "CLOSED THE LOOP. After remodel: " << perm.pose().fold_tree() << std::endl;

perm.chains_from_termini();
}
*/

} // namespace connection
} // namespace denovo_design
} // namespace protocols
