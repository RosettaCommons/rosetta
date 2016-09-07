// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file IterativeAbrelax
/// @brief iterative protocol starting with abinitio and getting progressively more concerned with full-atom relaxed structures
/// @details
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/abinitio/IterativeCentroid.hh>

// Package Headers

// Project Headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/Tracer.hh>

// Option Headers

//// C++ headers
#include <string>


static THREAD_LOCAL basic::Tracer tr( "protocols.iterative" );

using core::Real;
using namespace core;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;


namespace protocols {
namespace abinitio {
using namespace jd2::archive;

void IterativeCentroid::gen_diversity_pool( jd2::archive::Batch& batch, bool fullatom ) {
	if ( fullatom ) fullatom_pool_ptr_->gen_diversity_pool( batch, fullatom );
	else Parent::gen_diversity_pool( batch, fullatom );
}

void IterativeCentroid::update_noesy_filter_files(
	std::string const& current,
	bool fullatom
) {
	fullatom_pool_ptr_->update_noesy_filter_files( current, fullatom );
	Parent::update_noesy_filter_files( current, fullatom );
}


class OrderSortPredicate {
public:
	OrderSortPredicate() = default;
	bool operator() ( core::io::silent::SilentStructOP const& pss1, core::io::silent::SilentStructOP const& pss2 ) {
		return pss1->get_silent_energy( "_archive_index" ).value() < pss2->get_silent_energy( "_archive_index" ).value();
	}
};

void IterativeCentroid::erase_decoy( std::string const&  tag ) {
	Parent::erase_decoy( tag );
	SilentStructs::iterator it;
	for ( it = stage2_decoys_.begin(); it != stage2_decoys_.end(); ++it ) {
		if ( (*it)->decoy_tag() == tag ) break;
	}
	stage2_decoys_.erase( it );
}

void IterativeCentroid::collect_alternative_decoys( SilentStructs /*primary_decoys*/, std::string /*alternative_decoy_file*/, SilentStructVector& output_decoys ) {

	tr.Info << "resample_stage2 \n";

	for (auto & stage2_decoy : stage2_decoys_) {
		std::string batch_prefix( stage2_decoy->decoy_tag() );
		batch_prefix='f'+batch_prefix.substr(9,3);
		core::io::silent::SilentStructOP new_decoy=stage2_decoy->clone();
		new_decoy->set_decoy_tag( batch_prefix+"_"+stage2_decoy->decoy_tag().substr(14,5) );
		output_decoys.push_back( new_decoy );
	}

	//  typedef std::map< std::string, utility::vector1< std::string > > SourceFiles;
	//  typedef std::map< std::string, utility::vector1< core::io::silent::SilentStructOP > > AlternativeDecoys;

	//  SourceFiles sources;
	//  AlternativeDecoys alternative_decoys;
	//  Size ct_in( 0 );

	//  //to find the stage2 structures collect first all tags for a specific file
	//  for ( const_decoy_iterator it = primary_decoys.begin(); it != primary_decoys.end(); ++it ) {
	//   runtime_assert( (*it)->has_comment( TAG_IN_FILE ) );
	//   std::string tag( (*it)->get_comment( TAG_IN_FILE ) );
	//   utility::file::FileName file( (*it)->get_comment( SOURCE_FILE ) );
	//   std::string stage2_file( file.path()+"/"+alternative_decoy_file );

	//   //creates map <filename> <list of tags>
	//   sources[ stage2_file ].push_back( tag );
	//   alternative_decoys[ stage2_file ].push_back( (*it) );
	//   ++ct_in;
	//  }

	//  //read selected structures from each file
	//  Size ct_read( 0 );
	//  for ( SourceFiles::const_iterator it = sources.begin(); it != sources.end(); ++it ) {
	//   /// it->first is filename, it->second are all tags collected for this file
	//   io::silent::SilentFileData sfd;
	//   try { //read structures
	//    sfd._read_file( it->first, it->second, true /*throw exceptions */ );
	//    if ( sfd.size() > it->second.size() ) {
	//     tr.Warning << "[WARNING] multiple decoys with same tag detected in file " << it->first << std::endl;
	//    }
	//    //copy( sfd.begin(), sfd.end(), std::back_inserter( output_decoys ) );
	//    for ( core::io::silent::SilentFileData::iterator sit = sfd.begin(); sit != sfd.end(); ++sit ) {
	//     std::string batch_prefix( it->first );
	//      batch_prefix='f'+batch_prefix.substr(9,3);
	//     sit->set_decoy_tag( batch_prefix+"_"+sit->decoy_tag() );
	//     output_decoys.push_back( *sit );
	//    }
	//    ct_read += sfd.size();
	//   } catch ( utility::excn::EXCN_IO& excn ) { //ERROR
	//    tr.Warning << "[WARNING] Problem reading silent-file " << it->first << " for " << it->second.size() << " structures " << std::endl;
	//    excn.show( tr.Warning );
	//    tr.Warning << std::endl;
	//    tr.Warning << "use the respective structures in the pool as starting structure instead" << std::endl;
	//    copy( alternative_decoys[ it->first ].begin(), alternative_decoys[ it->first ].end(), std::back_inserter( output_decoys ) );
	//    ct_read += alternative_decoys[ it->first ].size();
	//   }
	//  }

	//  tr.Debug << "structures from pool" << ct_in << " structure retrieved from " << alternative_decoy_file << "-files "
	//       << ct_read << " start structs: " << output_decoys.size() << std::endl;
	//  if ( output_decoys.size() != primary_decoys.size() ) {
	//   tr.Warning << "[WARNING] why do we have a different number of decoys in pool and start_decoys ? " << std::endl;
	//  }
}

void IterativeCentroid::save_to_file( std::string suffix ) {
	Parent::save_to_file( suffix );
	if ( stage2_decoys_.size() ) {
		//add running number to stage2_decoys so that we can sort them later.
		core::Size ct( 1 );
		for ( SilentStructs::const_iterator it = stage2_decoys_.begin(); it != stage2_decoys_.end(); ++it ) {
			(*it)->add_energy( "_archive_index", ct++, 1.0 );
		}
		std::string const dirname( name() + suffix );
		save_decoys( dirname, "stage2_decoys", stage2_decoys_ );
	}
}

bool IterativeCentroid::restore_from_file() {
	std::string const& dirname( name() );
	std::string const filename ( dirname + "/stage2_decoys.out" );
	load_decoys( filename, stage2_decoys_ );
	stage2_decoys_.sort( OrderSortPredicate() );
	return Parent::restore_from_file();
}

/// @brief call to insert structure at position given by iterator
void IterativeCentroid::add_structure_at_position (
	SilentStructs::iterator iss,
	core::io::silent::SilentStructOP new_decoy,
	core::io::silent::SilentStructOP alternative_decoy
) {

	if ( alternative_decoy ) {
		auto alt_iss = stage2_decoys_.begin();
		if ( iss != decoys().end() ) {
			std::string const& tag( (*iss)->decoy_tag() );
			std::cerr << "try to add alternative decoy with tag " << alternative_decoy->decoy_tag() << " before decoy with tag " << tag << std::endl;
			//find position in sorted list to insert
			while ( alt_iss != stage2_decoys_.end() ) {
				if ( (*alt_iss)->decoy_tag() == tag ) break;
				++alt_iss;
			}
			if ( alt_iss == stage2_decoys_.end() ) {
				std::cerr << "[ERROR] decoy " << tag << " not found in alternative decoys " << std::endl;
			} else {
				runtime_assert( (*alt_iss)->decoy_tag() == tag );
			}
		} else {
			alt_iss = stage2_decoys_.end();
		}
		//if we are at the end this decoy has a worse score than all others
		//  if ( alt_iss != stage2_decoys_.end() || stage2_decoys_.size() < nstruct() ) {
		std::cerr << "equivalent decoy tag has been found in stage2_decoys, inserting now" << std::endl;
		if ( alt_iss == stage2_decoys_.end() ) tr.Debug << "inserting at end..." << std::endl;
		stage2_decoys_.insert( alt_iss, alternative_decoy );
		if ( stage2_decoys_.size() > max_nstruct() ) { //take all decoys until full
			stage2_decoys_.pop_back();
		}
		runtime_assert( new_decoy->decoy_tag() == alternative_decoy->decoy_tag() );
	} else {
		std::cerr << "[ERROR] no alternative decoy when inserting " << new_decoy->decoy_tag() << std::endl;
	}
	Parent::add_structure_at_position( iss, new_decoy, alternative_decoy );

	//check list is consistent
	SilentStructs::const_iterator it1, it2, eit1, eit2;
	it1=decoys().begin();
	eit1=decoys().end();
	it2=stage2_decoys_.begin();
	eit2=stage2_decoys_.end();
	while ( it1 != eit1 && it2 !=eit2 ) {
		if ( (*it1)->decoy_tag() != (*it2)->decoy_tag() ) {
			std::cerr <<"[ERROR] inconsistency in list" << std::endl;
			it1=decoys().begin();
			eit1=decoys().end();
			it2=stage2_decoys_.begin();
			eit2=stage2_decoys_.end();
			while ( it1 != eit1 && it2 !=eit2 ) {
				std::cerr << (*it1)->decoy_tag() << "   " << (*it2)->decoy_tag() <<std::endl;
				++it1;
				++it2;
			}
			runtime_assert( false );
		}
		++it1;
		++it2;
	}
	runtime_assert( it1==eit1 && it2==eit2 );
}

} //abinitio
} //protocols
