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
#include <protocols/topology_broker/SequenceClaimer.hh>

// Package Headers
#include <protocols/topology_broker/DofClaim.hh>
#include <protocols/topology_broker/Exceptions.hh>
// Project Headers
#include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// ObjexxFCL Headers

// Utility headers
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>


// utility
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/pose/annotated_sequence.hh>

//// C++ headers
#include <string>
#include <iterator>

static basic::Tracer tr("protocols.topo_broker",basic::t_info);

namespace protocols {
namespace topology_broker {

using namespace core;
SequenceClaimer::SequenceClaimer() :
	sequence_( "NO_SEQUENCE" ),
	rsd_type_set_( core::chemical::CENTROID ),
	offset_( 0 ),
	nr_res_( 0 )
{}

SequenceClaimer::SequenceClaimer( std::string const& sequence, std::string const& rsd_type_set_identifier, std::string label ) :
	rsd_type_set_( rsd_type_set_identifier ),
	offset_( 0 )
{
	set_sequence( sequence );
	set_label( label );
}

TopologyClaimerOP SequenceClaimer::clone() const
{
	return new SequenceClaimer( *this );
}


void SequenceClaimer::set_sequence( std::string const& str ) {
	sequence_ = str;
	pose::Pose my_pose;
	tr.Info << "make pose from sequence: " << str << std::endl;
	core::pose::make_pose_from_sequence(
				my_pose,
				sequence_,
				*( chemical::ChemicalManager::get_instance()->residue_type_set( rsd_type_set_ ))
	);
	annotated_sequence_ = my_pose.annotated_sequence();
	nr_res_ = my_pose.total_residue();
}

void SequenceClaimer::generate_sequence_claims( DofClaims& new_claims ) {

	new_claims.push_back( new SequenceClaim( this, 1, nr_res_, label(), DofClaim::INIT /* for now... eventually CAN_INIT ? */ ) );
}

void SequenceClaimer::generate_claims( DofClaims& new_claims ) {
	// that doesn't seem necessary. the atom-tree has them as lower and upper termini anyway
 	if ( offset() > 1 ) {
 		new_claims.push_back( new CutClaim( this, offset() - 1, DofClaim::INIT /* for now... eventually CAN_INIT ? */ ) );
 	}

	//special --- if only 1 residue chain... the torsion will be irrelvant and probably unclaimed
	// ... make broker happy but don't do anything...
	if ( nr_res_ == 1 ) new_claims.push_back( new BBClaim( this, offset() ) );
}

void SequenceClaimer::initialize_dofs( core::pose::Pose& pose, DofClaims const& init_claims, DofClaims& failed_to_init ) {
	if ( nr_res_ > 1 ) {
		TopologyClaimer::initialize_dofs( pose, init_claims, failed_to_init );
	} else { //special case the BBTorsion claim for position offset() is ours and is left unitialized...
		for ( DofClaims::const_iterator it = init_claims.begin(), eit = init_claims.end();
					it != eit; ++it ) {
			if ( (*it)->owner()==this ) {
				if ( (*it)->type() == DofClaim::BB  && (*it)->pos( 1 ) == offset() ) {
					tr.Trace << "SequenceClaimer: 1-residue chain --- no need to init: " << **it << std::endl;
					continue;
				}	else {
					tr.Trace << "SequenceClaimer: can't handle " << **it << std::endl;
					failed_to_init.push_back( *it );
				}
			} // our claim
		} //loop
	} //nr_res_ == 1
}

//bool SequenceClaimer::allow_claim( DofClaim const& foreign_claim ) {
// 	if ( foreign_claim.owner() == this ) return true; // always allow your own claims!
// 	if ( foreign_claim.type() == DofClaim::SEQUENCE ) {
// 		if ( foreign_claim.pos( 1 ) == 1  && foreing_claim.right() == DofClaim::EXCLUSIVE ) {
// 			return false; // don't accept any other fixed positions right now --- maybe never ?!
// 		}
// 	}
//	return true;
//} // SequenceClaimer::allow_claim()


void SequenceClaimer::initialize_residues( core::pose::Pose& pose, SequenceClaimOP my_claim, DofClaims& /*failed_to_init*/ ) {
	//need to copy coords and jumps --- if chunks were idealized no problem .... but non-idealized stuff ?
	//also take care of fullatom vs centroid...
	tr.Debug << "add sequence " << sequence_ << " to " << pose.annotated_sequence() << std::endl;
	offset_ = my_claim->offset();
	runtime_assert( offset_ == pose.total_residue() + 1 );
	runtime_assert( nr_res_ > 0 );
	core::pose::make_pose_from_sequence(
				 pose,
				 pose.annotated_sequence()+annotated_sequence_,
				 *( chemical::ChemicalManager::get_instance()->residue_type_set( rsd_type_set_ ))
	);
	// make extended chain
	for ( Size pos = 1; pos <= pose.total_residue(); pos++ ) {
		if ( ! pose.residue(pos).is_protein() ) continue;
		pose.set_phi( pos, -150 );
		pose.set_psi( pos, 150);
		pose.set_omega( pos, 180 );
	}

}

bool SequenceClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "file" || tag == "FILE" || tag == "file:" ) {
		is >> tag;
		set_sequence( core::sequence::read_fasta_file( tag )[1]->sequence() );
	} else if ( tag == "DEF" ) {
		while ( is >> tag && tag != "END_DEF" ) {
			if ( tag[0]=='#' ) {
				getline( is, tag );
				continue;
			}
			copy( tag.begin(), tag.end(), std::back_inserter( sequence_ ) );
			set_sequence( sequence_ ); //to get the dependent values updated!
		}
		if ( tag != "END_DEF" ) {
			throw EXCN_Input( "END_DEF expected after DEF while reading sequence" );
		}
	} else if ( tag == "CMD_FLAG" ) {
			using namespace basic::options;
			using namespace basic::options::OptionKeys;
			if ( option[ in::file::fasta ].user() ) {
				set_sequence( core::sequence::read_fasta_file( option[ in::file::fasta ]()[1] )[1]->sequence() );
				tr.Info << "read fasta sequence: " << sequence_.size() << " residues\n"  << sequence_ << std::endl;
			} else {
				throw EXCN_Input( "SequenceClaimer found tag CMD_FLAG but no -in:file:fasta on command-line");
			}
	} else return Parent::read_tag( tag, is );
	return true;
}

void SequenceClaimer::init_after_reading() {
	if ( nr_res_ == 0 ) {
		throw EXCN_Input( "no sequence found when reading " +type() );
	}
}

} //topology_broker
} //protocols
