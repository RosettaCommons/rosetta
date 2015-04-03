// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file IterativeAbrelax
/// @brief iterative protocol starting with abinitio and getting progressively more concerned with full-atom relaxed structures
/// @details
/// @author Oliver Lange

// Unit Headers
#include <protocols/jd2/archive/VarianceStatisticsArchive.hh>

// Package Headers
#include <core/io/silent/SilentStruct.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <basic/MemTracer.hh>
#include <numeric/random/random.hh>


static thread_local basic::Tracer tr( "protocols.iterative.VarianceStatistics" );
using basic::mem_tr;


using core::Real;


namespace protocols {
namespace jd2 {
namespace archive {


VarianceStatisticsArchive::VarianceStatisticsArchive( std::string name )
	: insertion_prob_( 0.1 )
{
	set_name( name );
}


bool VarianceStatisticsArchive::add_evaluated_structure(
  core::io::silent::SilentStructOP evaluated_decoy,
	core::io::silent::SilentStructOP /*alternative_decoy*/,
	Batch const&
) {
	if ( decoys().size() < nstruct() ) {
		tr.Debug << "added " << evaluated_decoy->decoy_tag() << " to " << name() << std::endl;
		decoys().insert( decoys().begin(), evaluated_decoy );
		invalidate_score_variations();
		return true;
	}

	if ( numeric::random::rg().uniform() < insertion_prob_ ) { //keep or not ?
		//replace with random element
		Size rg_pos( static_cast< int >( numeric::random::rg().uniform() * decoys().size() ) );
		runtime_assert( rg_pos < decoys().size() );
		SilentStructs::iterator it=decoys().begin();
		while ( rg_pos-- > 0 ) {
			++it;
		}
		runtime_assert( it != decoys().end() );
		*it=evaluated_decoy;
		invalidate_score_variations();
		return true;
	}

	return false;
}

}
} //abinitio
} //protocols
