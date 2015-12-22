// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/connection/ExtendChain.cc
/// @brief The BridgeChains
/// @detailed
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/connection/ExtendChain.hh>
#include <protocols/denovo_design/connection/ExtendChainCreator.hh>

//Project Headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/loops/Loops.hh>

//Protocol Headers

//Core Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>

//Basic Headers
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/tag/Tag.hh>

// ObjexxFCL Headers

//C++ Headers

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.connection.ExtendChain" );

///////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace connection {
///////////////////////////////////////////////////////////////////////////////

std::string
ExtendChainCreator::keyname() const
{
	return ExtendChainCreator::mover_name();
}

protocols::moves::MoverOP
ExtendChainCreator::create_mover() const
{
	return protocols::moves::MoverOP( new ExtendChain() );
}

std::string
ExtendChainCreator::mover_name()
{
	return "ExtendChain";
}

///  ---------------------------------------------------------------------------------
///  ExtendChain main code:
///  ---------------------------------------------------------------------------------

/// @brief default constructor
ExtendChain::ExtendChain() :
	BridgeChains()
{
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
ExtendChain::~ExtendChain() {}


/// Return a copy of ourselves
protocols::moves::MoverOP
ExtendChain::clone() const
{
	return protocols::moves::MoverOP( new ExtendChain( *this ) );
}

/// @brief return a fresh instance of ourselves
protocols::moves::MoverOP
ExtendChain::fresh_instance() const
{
	return protocols::moves::MoverOP( new ExtendChain() );
}

/// @brief return a fresh instance of ourselves
std::string
ExtendChain::get_name() const
{
	return ExtendChainCreator::mover_name();
}

void
ExtendChain::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose )
{
	BridgeChains::parse_my_tag( tag, data, filters, movers, pose );

	if ( tag->hasOption( "length" ) ) {
		std::stringstream msg;
		msg << id() << ": The length option is not valid for ExtendChain -- please specify a motif using the \"motif\" option!";
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	if ( !tag->hasOption( "motif" ) ) {
		std::stringstream msg;
		msg << id() << ": You must specify a motif (e.g. 1LG-1LB-10HA) to ExtendChain.";
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	if ( ( comp1_ids().size() && comp2_ids().size() ) ||
			( comp1_ids().size() && user_chain2() ) ||
			( user_chain1() && comp2_ids().size() ) ||
			( user_chain1() && user_chain2() ) ) {
		std::stringstream msg;
		msg << id() << ": You cannot set both segment1/chain1 and segment2/chain2 in ExtendChain. Please specify one or the other.";
		msg << "segment1_ids = " << comp1_ids() << " segment2_ids = " << comp2_ids() << std::endl;
		msg << "chain1 = " << user_chain1() << " chain2 = " << user_chain2() << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	if ( !comp1_ids().size() && !comp2_ids().size() ) {
		if ( !user_chain1() && !user_chain2() ) {
			std::stringstream msg;
			msg << id() << ": You must set either segment1/chain1 or segment2/chain2 in ExtendChain, but you haven't specified either. Please specify either one or the other.";
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
	}
}

/// @brief configures based on a permutation
/// @throw EXCN_Setup if no valid connection endpoints are found
void
ExtendChain::setup_permutation( components::StructureData & perm ) const
{
	core::Real rand = numeric::random::rg().uniform();
	core::Size const motifidx = extract_int( rand, 1, motifs().size() );
	Motif const & motif = motifs()[ motifidx ];

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
			throw utility::excn::EXCN_BadInput( id() + ": you can only specify one segment/chain to ExtendChain!  You have specified more than one." );
		}
		utility::vector1< std::string > const local_comp1_ids = find_available_upper_termini( perm );
		if ( user_chain1() > local_comp1_ids.size() ) {
			std::stringstream msg;
			msg << id() << ": The user-specified chain number (" << user_chain1()
				<< ") is larger than the number of upper termini in the pose (" << local_comp1_ids.size()
				<< "). Perm = " << perm << std::endl;
			throw utility::excn::EXCN_BadInput( msg.str() );
		}
		seg1 = local_comp1_ids[ user_chain1() ];
	}
	if ( user_chain2() ) {
		if ( !seg2.empty() ) {
			throw utility::excn::EXCN_BadInput( id() + ": you can only specify one segment/chain to ExtendChain!  You have specified more than one." );
		}
		utility::vector1< std::string > const local_comp2_ids = find_available_lower_termini( perm );
		if ( user_chain2() > local_comp2_ids.size() ) {
			std::stringstream msg;
			msg << id() << ": The user-specified chain2 number (" << user_chain2()
				<< ") is larger than the number of lower termini in the pose (" << local_comp2_ids.size()
				<< "). Perm = " << perm << std::endl;
			throw utility::excn::EXCN_BadInput( msg.str() );
		}
		seg2 = local_comp2_ids[ user_chain2() ];
	}

	if ( ( seg1.empty() && seg2.empty() ) || ( !seg1.empty() && !seg2.empty() ) ) {
		std::stringstream msg;
		msg << id() << ": you must specify exactly one valid segment/chain to ExtendChain!" << std::endl;
		msg << "user_chain1() = " << user_chain1() << " user_chain2() = " << user_chain2() << std::endl;
		msg << "segment1 = " << comp1_ids() << " segment2 = " << comp2_ids() << std::endl;
		msg << "Perm = " << perm << std::endl;
		throw utility::excn::EXCN_BadInput( msg.str() );
	}

	std::string const segname = id();
	std::stringstream complete_ss;
	complete_ss << 'L' << motif.ss << 'L';
	std::stringstream complete_abego;
	complete_abego << 'X' << motif.abego << 'X';

	components::StructureDataOP loop_sd(
		new components::SingleChainStructureData(
		segname,
		motif.len,
		motif.len + 2,
		false,
		complete_ss.str(),
		abego_vector( complete_abego.str() ) ) );

	setup_structuredata( *loop_sd,
		motif.len,
		motif.ss,
		motif.abego,
		seg1,
		seg2,
		"",
		"",
		0 );

	if ( !motif.len ) {
		TR << "Warning: length of motif given to ExtendChain is 0.  Doing nothing." << std::endl;
		return;
	}

	if ( perm.pose() ) {
		core::pose::PoseOP newpose = build_pose( *loop_sd );
		loop_sd->set_pose( newpose );
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
ExtendChain::process_permutation( components::StructureData & perm ) const
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
		throw utility::excn::EXCN_BadInput( msg.str() );
	}

	perm.chains_from_termini();

	// do the work
	try {
		apply_connection( perm );
	} catch ( EXCN_ConnectionFailed const & e ) {
		throw components::EXCN_Process( e.message() );
	}

	for ( core::Size i=1; i<=perm.pose()->total_residue(); ++i ) {
		TR << i << " : " << perm.pose()->residue( i ).name() << " " << perm.pose()->chain( i ) << std::endl;
	}
}

/// @brief Does the work of remodeling the connection
void
ExtendChain::apply_connection( components::StructureData & perm ) const
{
	if ( !loop_lower( perm ).empty() ) {
		connect_lower_loop( perm );
	} else if ( !loop_upper( perm ).empty() ) {
		connect_upper_loop( perm );
	} else {
		std::stringstream msg;
		msg << id() << ": Neither loop_upper (" << loop_upper( perm ) << ") nor loop_lower("
			<< loop_lower( perm ) << ") are the segment built by this connection. Perm=" << perm << std::endl;
		throw utility::excn::EXCN_BadInput( msg.str() );
	}

	// check sfxn
	if ( !scorefxn() ) {
		std::stringstream err;
		err << "ExtendChain: You must set a valid scorefunction to "
			<< id() << " before connecting" << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err.str() );
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

	debug_assert( perm.pose() );
	protocols::moves::MoverOP remodel =
		create_remodel_mover(
		*perm.pose(),
		loops,
		true,
		perm.ss(),
		perm.abego(),
		left,
		right );

	// store original in case there is a failure
	components::StructureDataOP orig = perm.clone();

	// switch to centroid if necessary
	bool const input_centroid = perm.pose()->is_centroid();
	core::pose::Pose stored;

	if ( !input_centroid ) {
		stored = *(perm.pose());
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

	// Rebuild fold tree
	utility::vector1< std::string > roots;
	if ( !lower_segment_id( perm ).empty() ) {
		roots.push_back( lower_segment_id( perm ) );
	}
	if ( !upper_segment_id( perm ).empty() ) {
		roots.push_back( upper_segment_id( perm ) );
	}
	perm.consolidate_movable_groups( roots );
	TR << "CLOSED THE LOOP. After remodel: " << perm.pose()->fold_tree() << std::endl;

	// switch back to FA if necessary
	if ( !input_centroid ) {
		perm.switch_residue_type_set( "fa_standard" );
		copy_rotamers( perm, stored );
	}
	perm.chains_from_termini();
}

} // namespace connection
} // namespace denovo_design
} // namespace protocols
