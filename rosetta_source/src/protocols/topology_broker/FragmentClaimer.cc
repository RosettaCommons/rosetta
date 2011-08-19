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
#include <protocols/topology_broker/DofClaim.hh>
#include <protocols/topology_broker/weights/AbinitioMoverWeight.hh>
#include <protocols/topology_broker/TopologyBroker.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh> //visualize
#include <utility/exit.hh>

#include <core/fragment/FragSet.hh>
#include <protocols/basic_moves/FragmentMover.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>

//#include <basic/options/option.hh>

//// C++ headers
// AUTO-REMOVED #include <fstream>

// option key includes


static basic::Tracer tr("protocols.topo_broker",basic::t_info);
//static numeric::random::RandomGenerator RG(181134);

namespace protocols {
namespace topology_broker {

using namespace core;

FragmentClaimer::FragmentClaimer() :
	mover_( NULL ),
	mover_tag_( "NoTag" ),
	bInitDofs_( false ),
	claim_right_( DofClaim::CAN_INIT )
{
	movemap_ = new kinematics::MoveMap;
}

FragmentClaimer::FragmentClaimer( basic_moves::FragmentMoverOP mover, std::string tag, weights::AbinitioMoverWeightOP weight ) :
	TopologyClaimer( weight ),
	mover_( mover ),
	mover_tag_( tag ),
	bInitDofs_( false ),
	claim_right_( DofClaim::CAN_INIT )
{
	movemap_ = new kinematics::MoveMap;
}

FragmentClaimer::FragmentClaimer( basic_moves::FragmentMoverOP mover ) :
	mover_( mover )
{
	if ( mover_) mover_tag_ = mover_->type();
	movemap_ = new kinematics::MoveMap;
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
	start_region.clear();
	tr.Trace << "FragmentClaimer::get_sequence_region" << std::endl;
	for ( core::Size i = 1; i<= active_sequence_labels_.size(); ++i ) {
		tr.Trace << " look for label : " << active_sequence_labels_[ i ] << std::endl;
		SequenceClaim const& seq_claim = broker().resolve_sequence_label( active_sequence_labels_[ i ] );
		for ( core::Size pos = seq_claim.offset(); pos <= seq_claim.last_residue(); ++pos ) {
			start_region.insert( pos );
		}
	}
}


void FragmentClaimer::generate_claims( DofClaims& new_claims ) {

	if ( region_.size() ) {
		region_.switch_movemap( *movemap_, core::id::BB, true );
	} else {
		movemap_->set_bb( true ); //Support only BB Claims right now
	}
	movemap_->set_jump( true );

	core::fragment::InsertMap insert_map;
	core::fragment::InsertSize insert_size;
	if ( !mover_ || !mover_->fragments() ) return;

	mover_->fragments()->generate_insert_map( *movemap_, insert_map, insert_size );

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
		Size const start ( *it );
		Size const length( insert_size[ start ] );
		new_claims.push_back( new BBClaim( this, *it, claim_right_ ) );
		for ( Size i = start + 1; i < start+length && insert_size[ i ] == 0; i++ ) {
			new_claims.push_back( new BBClaim( this, i, claim_right_ ) );
		}
	}
	tr.Trace << std::endl;
}

bool FragmentClaimer::accept_declined_claim( DofClaim const& was_declined ) {
	was_declined.toggle( *movemap_, false ); //this could even be BaseClass --- or BaseClass with MoveMap...
	tr.Debug << "OK: FragmentClaimer couldn't get " << was_declined << std::endl;
	return true; // full tolerance here ---
}

void FragmentClaimer::set_mover( basic_moves::FragmentMoverOP mover ) {
	mover_ = mover;
}

void FragmentClaimer::set_mover_tag( std::string const& str ) {
	mover_tag_ = str;
	if ( mover_ ) mover_->type( str );
}


basic_moves::FragmentMoverOP FragmentClaimer::get_frag_mover_ptr() {
	return mover_;
}


void FragmentClaimer::initialize_dofs( core::pose::Pose& pose, DofClaims const& init_dofs, DofClaims& /*failed_to_init*/ ) {
	//need to copy coords and jumps --- if chunks were idealized no problem .... but non-idealized stuff ?
	//also take care of fullatom vs centroid...

	core::fragment::InsertMap insert_map;
	core::fragment::InsertSize insert_size;
	kinematics::MoveMap test_map;
	test_map.set_bb( true );
	test_map.set_jump( true );
	if ( mover_ ) {
		mover_->fragments()->generate_insert_map( test_map, insert_map, insert_size );
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

	for ( DofClaims::const_iterator it = init_dofs.begin(), eit = init_dofs.end();
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
		tr.Info << type() << " " << label() << " init-dof-map: ";
		core::kinematics::simple_visualize_fold_tree_and_movemap( pose.fold_tree(), *movemap_, tr.Info );
		tr.Info << std::endl;
	}
	if ( mover_ ) mover_->set_movemap( movemap_ );
	tr.Info << type() << " " << label() << " movemap: ";
	core::kinematics::simple_visualize_fold_tree_and_movemap( pose.fold_tree(), *movemap_, tr.Info );
	tr.Info << std::endl;

}//initialize dofs

bool FragmentClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "REGION" ) {
		region_.read_stream_to_END( is, false /*no strict checking */, type(), "RIGID" );
	} else if ( tag == "region_file" || tag == "REGION_FILE" ) {
		std::string file;
		is >> file;
		region_.read_loop_file( file, false /*no strict looprlx checking*/, "RIGID" );  // <==
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
		active_sequence_labels_.push_back("main");
	}
}

moves::MoverOP FragmentClaimer::get_mover( core::pose::Pose const& /*pose*/ ) const {
	return mover_;
}

} //topology_broker
} //protocols
