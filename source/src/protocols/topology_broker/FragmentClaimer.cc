// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @detailed responsibilities:
/// @author Oliver Lange

// Unit Headers
#include <protocols/topology_broker/FragmentClaimer.hh>

// Package Headers
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/SequenceClaim.hh>
#include <protocols/topology_broker/claims/BBClaim.hh>
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/SequenceNumberResolver.hh>


// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh> //visualize
#include <utility/exit.hh>
#include <core/kinematics/Exceptions.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIterator.hh>
#include <protocols/simple_moves/FragmentMover.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/LoopsFileIO.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <sstream>

//#include <basic/options/option.hh>

//// C++ headers

// option key includes


static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;

FragmentClaimer::FragmentClaimer() :
	TopologyClaimer(),
	mover_( NULL ),
	mover_tag_( "NoTag" ),
	bInitDofs_( false ),
	claim_right_( claims::DofClaim::CAN_INIT )
{
	movemap_ = new kinematics::MoveMap;
}

FragmentClaimer::FragmentClaimer( simple_moves::FragmentMoverOP mover, std::string tag, weights::AbinitioMoverWeightOP weight ) :
	TopologyClaimer( weight ),
	mover_( mover ),
	mover_tag_( tag ),
	bInitDofs_( false ),
	claim_right_( claims::DofClaim::CAN_INIT )
{
	movemap_ = new kinematics::MoveMap;
	runtime_assert( fragments() );
}

FragmentClaimer::FragmentClaimer( simple_moves::FragmentMoverOP mover, std::string tag, weights::AbinitioMoverWeightOP weight, std::string label, core::fragment::FragSetOP frags ) :
	TopologyClaimer( weight ),
	mover_( mover ),
	mover_tag_( tag ),
	bInitDofs_( false ),
	claim_right_( claims::DofClaim::CAN_INIT )
{
	movemap_ = new kinematics::MoveMap;
	set_label( label );
	set_fragments( frags );
	runtime_assert( fragments() );
}


FragmentClaimer::FragmentClaimer( simple_moves::FragmentMoverOP mover ) :
	TopologyClaimer(),
	mover_( mover ),
	mover_tag_( "NoTag" ),
	bInitDofs_( false ),
	claim_right_( claims::DofClaim::CAN_INIT )
{
	if ( mover_) mover_tag_ = mover_->type();
	movemap_ = new kinematics::MoveMap;
	runtime_assert( fragments() );
}


/// @details  SHALLOW COPY to duplicate the pre-9/7/2009 behavior provided
/// by the compiler.
FragmentClaimer::FragmentClaimer( FragmentClaimer const & src ) :
	Parent( src )
{
	movemap_ = src.movemap_;
	mover_ = src.mover_;
	mover_tag_ = src.mover_tag_;
	bInitDofs_ = src.bInitDofs_;
	claim_right_ = src.claim_right_;
	region_ = src.region_;
}


FragmentClaimer::~FragmentClaimer() {}

void FragmentClaimer::get_sequence_region( std::set< Size >& start_region ) const {
	//TODO: this should probably use local positions rather than absolute positions, but AIN'T NOBODY GOT TIME FOR REFACTORING
	//runtime_assert( false );

	start_region.clear();
	tr.Trace << "FragmentClaimer::get_sequence_region" << std::endl;
	for ( core::Size i = 1; i<= active_sequence_labels_.size(); ++i ) {
		tr.Trace << " look for label : " << active_sequence_labels_[ i ] << std::endl;
		claims::SequenceClaim const& seq_claim = broker().resolve_sequence_label( active_sequence_labels_[ i ] );
		core::Size start_pos = broker().sequence_number_resolver().find_global_pose_number( active_sequence_labels_[ i ] );
		core::Size end_pos = start_pos + seq_claim.length() - 1 ;
		for ( core::Size pos = start_pos; pos <= end_pos; ++pos ) {
			start_region.insert( pos );
		}
	}
}


void FragmentClaimer::generate_claims( claims::DofClaims& new_claims ) {

	// At this point, the global sequence is set and the fragment positions have to be updated.
	core::Size fragment_offset = broker().sequence_number_resolver().offset( label() );

	core::Size seq_length = broker().resolve_sequence_label( label() ).length();
	core::Size frag_seq_length = fragments()->max_pos() - fragments()->min_pos() + 1;

	if ( seq_length != frag_seq_length ){
		std::ostringstream msg;
		msg << " Sequence length of SequenceClaim with label " << label() << "(length: "<< seq_length
				<< ") and sequence length of corresponding FragmentClaimer (length: " << frag_seq_length << ") do not match." << std::endl;
		throw utility::excn::EXCN_BadInput( msg.str() );
	}
  tr.Debug << "sequence length: " << seq_length << "; fragment sequence length: " << frag_seq_length << std::endl;

	tr.Debug << "Adapt fragment positions of FragmentClaimer " << label() << "-" << mover_tag_
			<< " to offset of " << fragment_offset <<std::endl;

	core::fragment::FragSetOP shifted_fragments = fragments()->clone_shifted( fragment_offset );
  assert( ( shifted_fragments->max_pos() - shifted_fragments->min_pos() ) ==
          ( fragments()->max_pos() - fragments()->min_pos() ) );
	set_fragments( shifted_fragments );

	if ( region_.size() ) {
		region_.switch_movemap( *movemap_, core::id::BB, true );
	} else {
		movemap_->set_bb( true ); //Support only BB Claims right now
	}
	movemap_->set_jump( true );

	core::fragment::InsertMap insert_map;
	core::fragment::InsertSize insert_size;
	if ( !mover_ || !fragments() ) return;

	fragments()->generate_insert_map( *movemap_, insert_map, insert_size );

	if ( insert_size.size() == 0 ){
		std::cerr << "Insert Size is 0 in  FragmentClaimer::generate_claims( .. ). File: " << __FILE__ << " Line: " << __LINE__ << std::endl;
		utility_exit_with_message("Cannot continue - Error generating insert map"	);
	}

	Size last_insert( insert_size.back() );
	for ( Size i = 1; i<=last_insert; ++i ) {
		insert_size.push_back( 0 );
	}

	Size const total_insert ( insert_map.size() );
	tr.Trace << "insert_size: ";
	if ( tr.Trace.visible() && total_insert ) for ( Size i = 1; i<=insert_map[ total_insert ]; i++ ) tr.Trace << " " << insert_size[ i ];

	for ( core::fragment::InsertMap::const_iterator it = insert_map.begin(), eit = insert_map.end();
				it != eit; ++it ) {
		Size const start ( *it - fragment_offset );
		Size const length( insert_size[ *it ] );
		new_claims.push_back( new claims::BBClaim( this, std::make_pair( label(), start ), claim_right_) );
		//new_claims.push_back( new claims::BBClaim( this, *it, claim_right_ ) );
		for ( Size i = start + 1; i < start+length && insert_size[ i + fragment_offset ] == 0; i++ ) {
			new_claims.push_back( new claims::BBClaim( this, std::make_pair( label(), i), claim_right_ ));
			//new_claims.push_back( new claims::BBClaim( this, i, claim_right_ ) );
		}
	}
	tr.Trace << std::endl;
}

bool FragmentClaimer::accept_declined_claim( claims::DofClaim const& was_declined ) {
	was_declined.toggle( *movemap_, false ); //this could even be BaseClass --- or BaseClass with MoveMap...
	tr.Debug << "OK: FragmentClaimer couldn't get " << was_declined << std::endl;
	return true; // full tolerance here ---
}

void FragmentClaimer::set_mover( simple_moves::FragmentMoverOP mover ) {
	mover_ = mover;
}

void FragmentClaimer::set_mover_tag( std::string const& str ) {
	mover_tag_ = str;
	if ( mover_ ) mover_->type( str );
}


simple_moves::FragmentMoverOP FragmentClaimer::get_frag_mover_ptr() {
	return mover_;
}


void FragmentClaimer::initialize_dofs( core::pose::Pose& pose, claims::DofClaims const& init_dofs, claims::DofClaims& /*failed_to_init*/ ) {
	//need to copy coords and jumps --- if chunks were idealized no problem .... but non-idealized stuff ?
	//also take care of fullatom vs centroid...

	core::fragment::InsertMap insert_map;
	core::fragment::InsertSize insert_size;
	kinematics::MoveMap test_map;
	test_map.set_bb( true );
	test_map.set_jump( true );
	if ( mover_ ) {
		fragments()->generate_insert_map( test_map, insert_map, insert_size );
	}

	Size const total_insert ( insert_map.size() );
	tr.Trace << "size of insertmap: " << total_insert << " -- ";
	for ( Size i = 1; i<=total_insert; i++ ) tr.Trace << " " << insert_map[ i ];
	tr.Trace << "insert_size: ";
	if ( total_insert ) for ( Size i = 1; i<=insert_map[ total_insert ]; i++ ) tr.Trace << " " << insert_size[ i ];
	tr.Trace << std::endl;

	kinematics::MoveMapOP init_map = new kinematics::MoveMap;
	init_map->set_bb( false );
	init_map->set_jump( false );

	for ( claims::DofClaims::const_iterator it = init_dofs.begin(), eit = init_dofs.end();
				it != eit; ++it ) {
		if ( (*it)->owner()==this ) {
			(*it)->toggle( *init_map, true );
			//don't really know how this looks for jumps
			//			Size pos( (*it)->pos( 1 ) );
			//if ( pos >= insert_size.size() + insert_size.back() || ( pos<= insert_size.size() && !insert_size[ pos ] ) ) {
				//				failed_to_init.push_back( *it );
			//}
		}
	}

	if ( mover_ && bInitDofs_ ) {
		mover_->set_movemap( init_map );
		mover_->apply_at_all_positions( pose );
		if( tr.Debug.visible() ){
			tr.Debug << type() << " " << label() << " init-dof-map: ";
			core::kinematics::simple_visualize_fold_tree_and_movemap( pose.fold_tree(), *movemap_, tr.Debug );
			tr.Debug << std::endl;
		}
	}
	if ( mover_ ) mover_->set_movemap( movemap_ );

	if( tr.Debug.visible() ){
		tr.Info << type() << " " << label() << " movemap: ";
		core::kinematics::simple_visualize_fold_tree_and_movemap( pose.fold_tree(), *movemap_, tr.Info );
		tr.Info << std::endl;
	}
}//initialize dofs

bool FragmentClaimer::read_tag( std::string tag, std::istream& is ) {
	loops::PoseNumberedLoopFileReader loop_file_reader;
	loop_file_reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
	if ( tag == "REGION" ) {
		loops::SerializedLoopList loops = loop_file_reader.read_pose_numbered_loops_file( is, type(), false /*no strict checking */ );
		region_ = loops::Loops( loops );
	} else if ( tag == "region_file" || tag == "REGION_FILE" ) {
		std::string file;
		is >> file;
		std::ifstream infile( file.c_str() );

		if (!infile.good()) {
			utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + file + "'" );
		}
		loops::SerializedLoopList loops = loop_file_reader.read_pose_numbered_loops_file( infile, file, false /*no strict checking */ );
		region_ = loops::Loops( loops ); // <==
	} else if ( tag == "SEQUENCE_REGION" ) {
		std::string label;
		is >> label;
		active_sequence_labels_.push_back( label );
		tr.Trace << type() << " reads SEQUENCE_REGION with label " << label;
	} else return Parent::read_tag( tag, is );
	return true;
}

void FragmentClaimer::init_after_reading() {
	if ( active_sequence_labels_.size() == 0 ) {
		active_sequence_labels_.push_back( label() );
	}
}

moves::MoverOP FragmentClaimer::get_mover( core::pose::Pose const& /*pose*/ ) const {
	return mover_;
}

void FragmentClaimer::set_fragments( core::fragment::FragSetOP frags ){
	mover_->set_fragments( frags );
}

core::fragment::FragSetCOP FragmentClaimer::fragments() {
	return mover_->fragments();
}

} //topology_broker
} //protocols
