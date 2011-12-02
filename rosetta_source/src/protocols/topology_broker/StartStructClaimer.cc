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
#include <protocols/topology_broker/StartStructClaimer.hh>

// Package Headers
#include <protocols/topology_broker/DofClaim.hh>
#include <protocols/topology_broker/Exceptions.hh>
// Project Headers
#include <core/pose/Pose.hh>

#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/util.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/BBTorsionSRFD.hh>

#include <protocols/simple_moves/FragmentMover.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/excn/Exceptions.hh>
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>
//#include <basic/options/option.hh>
#include <numeric/random/random.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


//// C++ headers

static numeric::random::RandomGenerator RG(199234234);


static basic::Tracer tr("protocols.topo_broker",basic::t_info);
//static numeric::random::RandomGenerator RG(181134);

namespace protocols {
namespace topology_broker {

using namespace core;

StartStructClaimer::StartStructClaimer() : bUseInputPose_( true ), perturb_( 0.0 )
{
	set_claim_right( DofClaim::INIT );
	set_bInitDofs( true );
}

StartStructClaimer::StartStructClaimer( core::pose::Pose const& pose ) : bUseInputPose_( true ), perturb_( 0.0 )
{
	set_claim_right( DofClaim::INIT );
	//	generate_init_frags( pose );
	start_pose_ = pose;
	set_bInitDofs( true );
}

void StartStructClaimer::new_decoy() {
	//	start_pose_.clear();
}

void StartStructClaimer::new_decoy( core::pose::Pose const& pose ) {
	if ( bUseInputPose_ ) {
		start_pose_ = pose;
	}
	//CANNOT do this here since SEQUENCE_CLAIMS need to be resolved first to know the sequence-offset
	// generate_init_frags( start_pose_ );

	//	start_pose_.clear();
}

void StartStructClaimer::generate_sequence_claims( DofClaims& new_claims ) {
	new_claims.push_back( new SequenceClaim( this, 0, 0, label(), DofClaim::NEED_TO_KNOW /* for now... eventually CAN_INIT ? */ ) );
}

void StartStructClaimer::generate_init_frags( core::pose::Pose const& pose ) {
	std::set< Size > start_region;
	get_sequence_region( start_region );
 	//for ( Size i = 1; i<= pose.total_residue(); ++i ) {
	// 		start_region.insert( i );
	// 	}
	if (tr.Trace.visible() ) {
		tr.Trace << " start region for StartStructClaimer "<< std::endl;
		for ( std::set< Size >::iterator it = start_region.begin(); it != start_region.end(); ++it ) {
			tr.Trace << *it << " ";
		}
		tr.Trace << std::endl;
	}


	using namespace fragment;
	ConstantLengthFragSetOP fragset = new ConstantLengthFragSet( 1 );
	steal_frag_set_from_pose( pose, *fragset, new FragData( new BBTorsionSRFD, 1 ), start_region );
	simple_moves::ClassicFragmentMoverOP mover = new simple_moves::ClassicFragmentMover( fragset );
	mover->set_check_ss( false ); /* not good if we want to initialize from 1mer fragments */
	set_mover( mover );

}

void StartStructClaimer::initialize_residues( core::pose::Pose&, SequenceClaimOP , DofClaims&  ) {
	//now the sequence is known:
	if ( start_pose_.total_residue() > 0 ) {
		generate_init_frags( start_pose_ ); //requires sequence definition --- had to wait until now.
	}
	//start_pose_.clear(); // save the space..
}


void StartStructClaimer::initialize_dofs(
		 core::pose::Pose& pose,
		 DofClaims const& init_dofs,
		 DofClaims& failed_to_init
) {

	try{
		mover();
	} catch( utility::excn::EXCN_NullPointer excn ) {
		throw( EXCN_Input( "StartStructureClaimer needs JobInputter or FILE <pdb-file> entry in broker-setup"));
	}

	FragmentClaimer::initialize_dofs( pose, init_dofs, failed_to_init );
	if ( perturb_ == 0.0 ) return;
	for ( DofClaims::const_iterator it = init_dofs.begin(), eit = init_dofs.end();
				it != eit; ++it ) {
		//don't really know how this looks for jumps
		if ( (*it)->type() == DofClaim::BB ) {
			Size pos( (*it)->pos( 1 ) );
			pose.set_phi( pos, pose.phi( pos ) + RG.gaussian()*perturb_ );
			pose.set_psi( pos, pose.psi( pos ) + RG.gaussian()*perturb_ );
		}
	}
}

bool StartStructClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "file" || tag == "FILE" ) {
		std::string filename;
		is >> filename;
		tr.Debug << type() << " initialized with file " << filename << std::endl;
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( start_pose_,
			*core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID ),
			filename );
	} else if ( tag == "PERTURB" ) {
		is >> perturb_;
	} else if ( tag =="NO_USE_INPUT_POSE" ) {
		bUseInputPose_ = false;
	} else return FragmentClaimer::read_tag( tag, is );
	return true;
}



} //topology_broker
} //protocols
