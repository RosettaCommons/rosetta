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
#include <protocols/topology_broker/claims/DofClaim.hh>
#include <protocols/topology_broker/claims/BBClaim.hh>
#include <protocols/topology_broker/claims/SequenceClaim.hh>
#include <protocols/topology_broker/claims/CutClaim.hh>
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

#include <core/pose/annotated_sequence.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <utility>

#ifdef WIN32
#include <iterator>
#endif

static thread_local basic::Tracer tr( "protocols.topo_broker", basic::t_info );

namespace protocols {
namespace topology_broker {

using namespace core;
SequenceClaimer::SequenceClaimer() :
	rsd_type_set_( core::chemical::CENTROID ),
	priority_( 0.0 ),
	input_sequence_( "" ),
	sequence_claim_( NULL )
{}

SequenceClaimer::SequenceClaimer( std::string const& sequence, std::string const& label,
	std::string const& rsd_type_set_identifier = core::chemical::CENTROID ) :
	rsd_type_set_( rsd_type_set_identifier ),
	priority_( 0.0 ),
	input_sequence_( sequence ),
	sequence_claim_( NULL )
{
	Parent::set_label(label);
}


TopologyClaimerOP SequenceClaimer::clone() const
{
	return new SequenceClaimer( *this );
}

void SequenceClaimer::make_sequence_claim() {
	tr.Debug << "SequenceClaimer " << label() << " making pose from sequence '"
            << input_sequence_ << "' to get annotated sequence." << std::endl;
	pose::Pose my_pose;
	core::pose::make_pose_from_sequence(
		my_pose,
		input_sequence_,
		*( chemical::ChemicalManager::get_instance()->residue_type_set( rsd_type_set_ ))
	);
	sequence_claim_ =
		new claims::SequenceClaim(
				this,
				my_pose.annotated_sequence(),
				label(),
				priority_
		);
	runtime_assert( my_pose.total_residue() == sequence_claim_->length() );
}

void SequenceClaimer::generate_sequence_claims( claims::DofClaims& new_claims ) {
	if ( !sequence_claim_ ) {
		make_sequence_claim();
	}
	new_claims.push_back( sequence_claim_ );
}

void SequenceClaimer::read_fasta_file( std::string file ) {
	using namespace core::sequence;
	utility::vector1< SequenceOP > sequences = core::sequence::read_fasta_file( file );
	if ( sequences.size() > 1 ) {
		throw EXCN_Input( "SequenceClaimer found multiple sequences in fasta file '"+file+
			"'\nTo resolve this, split fasta file into multiple files, and claim using multiple SequenceClaimers.");
	}
	input_sequence_ = sequences[1]->sequence();
    tr.Debug << "SequenceClaimer read sequence " << input_sequence_
            << " from fasta file '" << file << "'" << std::endl;
}

bool SequenceClaimer::read_tag( std::string tag, std::istream& is ) {
	if ( tag == "file" || tag == "FILE" || tag == "file:" ) {
		is >> tag;
		read_fasta_file(tag);
	} else if ( tag == "PRIORITY" || tag == "Priority" || tag == "priority" ) {
		is >> priority_;
	}	else if ( tag == "DEF" ) {
		//Expects input of the form:
		//DEF
		//HHHHHHAG GMPADFADF ADFAKDFLAK ADKFJ
		//FFFFF ADKFAJLDFJ ADFKAJ LFJ
		//END_DEF

		if ( input_sequence_ != "" ) {
			throw EXCN_Input( "SequenceClaimer found multiple definitions for its sequence. Ambiguity cannot be resolved.");
		}
		while ( is >> tag && tag != "END_DEF" ) {
			if ( tag[0]=='#' ) {
				getline( is, tag );
				continue;
			}
			input_sequence_ += tag;
		}
		if ( tag != "END_DEF" ) {
			throw EXCN_Input( "END_DEF expected after DEF while reading sequence" );
		}
	} else if ( tag == "CMD_FLAG" ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		if ( option[ in::file::fasta ].user() ) {
			read_fasta_file( option[ in::file::fasta ]()[1] );
		} else {
			throw EXCN_Input( "SequenceClaimer found tag CMD_FLAG but option '-in:file:fasta' is not set on cmd-line" );
		}
	} else return Parent::read_tag( tag, is );
	return true;
}

void SequenceClaimer::init_after_reading() {
	if ( input_sequence_.size() == 0 ) {
		throw EXCN_Input( "No sequence found when reading " +type() );
	}
	if ( label() == "NO_LABEL" ) {
		throw EXCN_Input( "SequenceClaimer left unlabeled: use the LABEL command to specify a label for each SequenceClaimer" );
	}
}

void SequenceClaimer::generate_claims( claims::DofClaims& new_claims ) {
	//Make a new cut at the end of this sequence.
	//TODO: Make TopologyBroker get rid of cuts outside valid sequence.
	new_claims.push_back( new claims::CutClaim( this, std::make_pair( label(), sequence_claim_->length() ) ) );

	//special --- if only 1 residue chain... the torsion will be irrelvant and probably unclaimed
	// ... make broker happy but don't do anything...
	if ( sequence_claim_->length() == 1 ) {

		new_claims.push_back( new claims::BBClaim( this, std::make_pair( label(), 1 ) ) );
	}
}


// void SequenceClaimer::initialize_residues( core::pose::Pose& pose, claims::SequenceClaimOP my_claim, claims::DofClaims& /*failed_to_init*/ ) {
// 	//need to copy coords and jumps --- if chunks were idealized no problem .... but non-idealized stuff ?
// 	//also take care of fullatom vs centroid...

// 	//    core::sequence::SequenceOP seq = NULL;
// 	//    for( Size i=1; i <= sequences_.size(); ++i){
// 	//        if( sequences_.at(i)->id() == my_claim->label() ){
// 	//            seq = sequences_.at(i);
// 	//        }
// 	//    }
// 	//    runtime_assert(seq != NULL);

// }


} //topology_broker
} //protocols




//bool SequenceClaimer::allow_claim( claims::DofClaim const& foreign_claim ) {
// 	if ( foreign_claim.owner() == this ) return true; // always allow your own claims!
// 	if ( foreign_claim.type() == claims::DofClaim::SEQUENCE ) {
// 		if ( foreign_claim.pos( 1 ) == 1  && foreing_claim.right() == claims::DofClaim::EXCLUSIVE ) {
// 			return false; // don't accept any other fixed positions right now --- maybe never ?!
// 		}
// 	}
//	return true;
//} // SequenceClaimer::allow_claim()

//void SequenceClaimer::initialize_dofs( core::pose::Pose& pose, claims::DofClaims const& init_claims,	claims::DofClaims& failed_to_init ) {
	//special case the BBTorsion claim for position offset() is ours and is left unitialized...
	//	for ( claims::DofClaims::const_iterator it = init_claims.begin(), eit = init_claims.end();
	//		it != eit; ++it ) {
		//		claims::BBClaimOP bb_ptr( dynamic_cast< claims::BBClaim* >( it->get() ) );

