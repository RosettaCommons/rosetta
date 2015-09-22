// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file IterativeAbrelax
/// @brief iterative protocol starting with ab initio and getting progressively more concerned with full-atom relaxed structures
/// @details
/// @author Oliver Lange

// Unit Headers
#include <protocols/abinitio/HedgeArchive.hh>

#include <protocols/jd2/archive/ArchiveManager.hh>

// Package Headers
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <ObjexxFCL/string.functions.hh>
#include <utility/file/file_sys_util.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.iterative.HedgeArchive" );
using basic::mem_tr;


using core::Real;


namespace protocols {
namespace abinitio {


HedgeArchive::HedgeArchive( std::string name ) :
	score_cut_per_batch_( 0.1 ),
	add_fuzzy_( 0.1 ) //0 for strictly score based, 1 for totally random
{
	set_name( name );
	set_max_nstruct( static_cast<Size>(1e6) ); //more than this and we have clearly too much file system load.
	set_evaluate_local( false ); //never re-evaluate decoys in this Archive
}

void HedgeArchive::incorporate_batch( core::Size batch_id ) {
	tr.Debug << "batch " << batch_id << " has finished... incorporating into hedge archive " << std::endl;
	SilentStructs& sorted_decoys( incoming_structures_[ batch_id ] );
	sorted_decoys.sort();
	Size ind_max( static_cast< Size > ( sorted_decoys.size()*score_cut_per_batch_ ) );
	for ( SilentStructs::const_iterator sit = sorted_decoys.begin(); sit != sorted_decoys.end() && ind_max>0; ++sit ) {
		if ( numeric::random::rg().uniform() < (1-add_fuzzy_) ) {
			set_max_nstruct( static_cast<Size>(1e6) );
			tr.Debug << "add decoy from batch " << batch_id << std::endl;
			add_structure_at_position( decoys().end(), sit->second, NULL );
			--ind_max;
		}
	}
	incoming_structures_.erase( batch_id );
	remove_pending_decoys( batch_id );
	old_batches_.insert( batch_id );
}

bool HedgeArchive::add_evaluated_structure(
	core::io::silent::SilentStructOP evaluated_decoy,
	core::io::silent::SilentStructOP /*alternative_decoy*/,
	jd2::archive::Batch const& batch
) {
	core::io::silent::SilentStructOP incoming_decoy = evaluated_decoy->clone(); //clone because that makes it more consistent after restart...
	if ( old_batches_.find( batch.id() ) != old_batches_.end() ) return false;
	incoming_structures_[ batch.id() ].push_back( std::make_pair( select_score( incoming_decoy ), incoming_decoy ) );
	tr.Debug << "incoming from batch: " << batch.id() << " " << incoming_structures_[ batch.id() ].size() << " " << batch.decoys_returned() << " " << batch.has_finished() << std::endl;
	if ( batch.has_finished() && incoming_structures_[ batch.id() ].size() >= batch.decoys_returned() ) { //90% of structures have arrived
		incorporate_batch( batch.id() );
	}
	return true;
}

std::string filename( core::Size batch_id ) {
	return std::string("pending_")+ObjexxFCL::lead_zero_string_of( batch_id , 4 );
}

void HedgeArchive::save_pending_decoys( SilentStructs const& decoys, core::Size batch_id ) const {
	std::string const& dirname( name() );
	std::string const ffilename ( dirname + "/" + filename( batch_id ) );
	std::string const backup_filename ( ffilename+".backup" );
	std::string const tmp_filename ( ffilename+".tmp" );
	//handle output myself... so it keeps the order of decoys.

	core::io::silent::SilentFileData sfd;
	//  utility::io::ozstream output( tmp_filename );
	//  if ( decoys.begin() != decoys.end() ) (*decoys.begin())->print_header( output );
	for ( SilentStructs::const_iterator it = decoys.begin(); it != decoys.end(); ++it ) {
		sfd.add_structure( *(it->second) ); //only add OP to sfd
	}
	sfd.write_all( tmp_filename );

	//rename to final
	if ( utility::file::file_exists( ffilename ) ) {
		rename( ffilename.c_str(), backup_filename.c_str() );
	}
	rename( tmp_filename.c_str(), ffilename.c_str() );
}

void HedgeArchive::remove_pending_decoys( core::Size batch_id ) const {
	std::string const& dirname( name() );
	std::string const ffilename ( dirname + "/" + filename( batch_id ) );
	std::string const backup_filename ( ffilename+".backup" );
	std::string const tmp_filename ( ffilename+".tmp" );
	utility::file::file_delete( ffilename );
	utility::file::file_delete( backup_filename );
	utility::file::file_delete( tmp_filename );
}


void HedgeArchive::collect( jd2::archive::Batch const& batch, core::io::silent::SilentStructOPs& start_decoys ) const {
	utility::vector1< core::Size > indices( decoys().size() );
	for ( core::Size i=1; i<=decoys().size(); ++i ) {
		indices[i]=i<=batch.nstruct();
	}
	numeric::random::random_permutation(indices, numeric::random::rg());
	numeric::random::random_permutation(indices, numeric::random::rg());
	core::Size i=1;
	start_decoys.reserve( batch.nstruct() );
	for ( Parent::SilentStructs::const_iterator it = decoys().begin(); it != decoys().end(); ++it, ++i ) {
		if ( indices[i] ) {
			start_decoys.push_back( *it );
		}
	}
}

void HedgeArchive::save_status( std::ostream& os ) const {
	Parent::save_status( os );
	os << "OPEN_BATCHES" << std::endl;
	for ( BatchStructuresMap::const_iterator it=incoming_structures_.begin(); it!=incoming_structures_.end(); ++it ) {
		if ( !it->second.size() ) continue;
		os << it->first << std::endl;
		HedgeArchive::save_pending_decoys( it->second, it->first );
	}
	os << "END_BATCHES" << std::endl;
}

void HedgeArchive::restore_status( std::istream& is ) {
	Parent::restore_status( is );
	std::string tag;
	is >> tag;
	runtime_assert( tag == "OPEN_BATCHES" );
	while ( is >> tag ) {
		if ( tag == "END_BATCHES" ) break;
		Size batch_id = ObjexxFCL::int_of( tag );
		std::string const& dirname( name() );
		std::string const ffilename ( dirname + "/" + filename( batch_id ) );
		//std::string const backup_filename ( ffilename+".backup" );
		//std::string const tmp_filename ( ffilename+".tmp" );
		if ( utility::file::file_exists( ffilename ) ) {
			using namespace core::io::silent;
			SilentFileData sfd;
			if ( !sfd.read_file( ffilename ) ) throw ( utility::excn::EXCN_BadInput( "problem reading silent file"+ffilename ) );
			for ( SilentFileData::iterator it=sfd.begin(), eit=sfd.end(); it!=eit; ++it ) {
				incoming_structures_[ batch_id ].push_back( std::make_pair( select_score( *it ), *it ) );
			}
		}
	}
}


} //abinitio
} //protocols
