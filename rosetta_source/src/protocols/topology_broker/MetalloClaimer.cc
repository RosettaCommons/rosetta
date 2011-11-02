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
#include <protocols/topology_broker/MetalloClaimer.hh>

// Package Headers
#include <protocols/topology_broker/DofClaim.hh>
#include <protocols/topology_broker/Exceptions.hh>
#include <protocols/topology_broker/TopologyBroker.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>

#include <protocols/jumping/ResiduePairJumpSetup.hh>
#include <protocols/jumping/ResiduePairJump.hh>
#include <protocols/jumping/JumpSetup.hh>

// ObjexxFCL Headers

// Utility headers
// AUTO-REMOVED #include <core/sequence/util.hh>


//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
//#include <basic/options/option.hh>

//// C++ headers
// AUTO-REMOVED #include <fstream>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector1.hh>


static basic::Tracer tr("protocols.topo_broker",basic::t_info);

namespace protocols {
namespace topology_broker {

using namespace core;
MetalloClaimer::MetalloClaimer() : residue_pair_jump_( NULL )
{}

//bool MetalloClaimer::allow_claim( DofClaim const& foreign_claim ) {
// 	if ( foreign_claim.owner() == this ) return true; // always allow your own claims!
// 	if ( foreign_claim.type() == DofClaim::SEQUENCE ) {
// 		if ( foreign_claim.pos( 1 ) == 1  && foreing_claim.right() == DofClaim::EXCLUSIVE ) {
// 			return false; // don't accept any other fixed positions right now --- maybe never ?!
// 		}
// 	}
//	return true;
//} // MetalloClaimer::allow_claim()
void MetalloClaimer::initialize_residues( core::pose::Pose& pose, SequenceClaimOP my_claim, DofClaims& failed_to_init ) {
	SequenceClaimer::initialize_residues( pose, my_claim, failed_to_init );
	// I know this guy isn't doing anything JumpClaimer::initialize_residues( pose, my_claim, failed_to_init );

	resolved_anchor_residue_ = broker().resolve_residue( anchor_chain_, anchor_residue_ );
	tr.Trace << "MetalloClaimer: setup jump between " << anchor_residue_ << " " << my_claim->offset() << std::endl;
	jump_setup_->add_jump(
			jumping::Interval( anchor_residue_, my_claim->offset()), //jump
			jumping::Interval( my_claim->offset() - 1, my_claim->offset() - 1 ) //cutpoint interval
	);
	new_decoy();
}

void MetalloClaimer::add_constraints( core::pose::Pose& ) {
}

void  MetalloClaimer::set_defaults() {
	SequenceClaimer::set_defaults();
	JumpClaimer::set_defaults();
	anchor_chain_ = ""; //usually anchored to DEFAULT chain
	anchor_residue_ = 0;
	residue_pair_jump_ = new jumping::ResiduePairJump;
}

bool MetalloClaimer::read_tag( std::string tag, std::istream& is ) {
	using namespace jumping;
	if ( tag == "ligand" ) {
		is >> ligand_;
	} else if ( tag == "anchor" ) {
		is >> anchor_residue_;
	} else if ( tag == "aa" ) {
		for ( int i = 1; i <= 2; ++i ) {
			std::string name;
			core::chemical::ResidueTypeSetCAP residue_type_set(
					core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );
			is >> name;
			if ( residue_type_set->name_map(name).is_protein() )
				name = name + "_p:CtermProteinFull_p:NtermProteinFull";
			core::chemical::ResidueType const & res_type( residue_type_set->name_map(name) );
			residue_pair_jump_->add_residue_single( res_type );
		}
	} else if ( tag == "cst_atoms" ) {
		for ( int i = 1; i <= 2; ++i ) {
			std::string name;
			for ( int j = 1; j <= 3; ++j ) {
				is >> name;
				residue_pair_jump_->set_cstAtoms(i,j,name);
			}
		}
	} else if ( tag == "jump_atoms" ) {
		for ( int i = 1; i <= 2; ++i ) {
			std::string name;
			for ( int j = 1; j <= 3; ++j ) {
				is >> name;
				residue_pair_jump_->set_jumpAtoms(i,j,name);
			}
		}
	} else if ( tag == "disAB"
			|| tag == "angleA" || tag == "angleB"
			|| tag == "dihedralA" || tag == "dihedralB" || tag == "dihedralAB" ) {
		//these are all multi valued read until line ends
		std::string line;
		getline( is, line );
		std::istringstream in( line );
		Real value;
		if ( tag == "disAB" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( disAB, value );
			}
		} else if ( tag == "angleA" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( angleA, value );
			}
		} else if ( tag == "angleB" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( angleB, value );
			}
		} else if ( tag == "dihedralA" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( dihedralA, value );
			}
		} else if ( tag == "dihedralB" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( dihedralB, value );
			}
		} else if ( tag == "dihedralAB" ) {
			while (in >> value) {
				residue_pair_jump_->set_cstInfo( dihedralAB, value );
			}
		} else runtime_assert( 0 ); //if you are here you missed a tag in the statement above.
	} else if ( SequenceClaimer::read_tag( tag, is ) ) {
		//noop
	} else if ( JumpClaimer::read_tag( tag, is ) ) {
		//noop
	} else return false;
	return true;
}

void MetalloClaimer::init_after_reading() {
	if ( !anchor_residue_ ) {
		throw EXCN_Input( "need to specify anchor residue for MetalloLigand "+ligand_ );
	}

	set_label( ligand_ + ObjexxFCL::string_of( anchor_residue_ ) + anchor_chain_ );
	if ( !anchor_chain_.size() ) anchor_chain_ = "DEFAULT"; //do this after labelling so we don't have "DEFAULT" appearing in the label

	set_sequence( "Z["+ligand_+"]" );

	residue_pair_jump_->init_mini_pose();
	jump_setup_ = new	jumping::ResiduePairJumpSetup;
	jump_setup_->add_residue_pair_jump( residue_pair_jump_ );
	set_jump_def( jump_setup_ ); //tell the underlying JumpClaimer
}



} //topology_broker
} //protocols
