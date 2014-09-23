// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSamplerComb.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSamplerComb.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSampler.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.sampler.copy_dofs.ResidueAlternativeStepWiseSamplerComb" );

using namespace core;
using core::conformation::Residue;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace copy_dofs {

	//Constructor
	ResidueAlternativeStepWiseSamplerComb::ResidueAlternativeStepWiseSamplerComb()
	{}

	//Destructor
	ResidueAlternativeStepWiseSamplerComb::~ResidueAlternativeStepWiseSamplerComb()
	{}

	/// @brief Add one more rotamer sampler to this sampler
	void
	ResidueAlternativeStepWiseSamplerComb::add_residue_alternative_rotamer( ResidueAlternativeStepWiseSamplerOP const & rotamer ) {
		residue_alternative_rotamer_map_[ rotamer->representative_seqpos() ] = rotamer;
		residue_alternative_rotamer_list_.push_back( rotamer );
		rotamer_list_.push_back( rotamer );
		set_init( false );
	}

	bool
	ResidueAlternativeStepWiseSamplerComb::has_resnum( Size const seqpos ){
		return ( residue_alternative_rotamer_map_.find( seqpos ) != residue_alternative_rotamer_map_.end() );
	}

	Size
	ResidueAlternativeStepWiseSamplerComb::find_resnum( Size const seqpos ){
		for ( Size i = 1; i <= residue_alternative_rotamer_list_.size(); i++ ){
			if ( residue_alternative_rotamer_list_[i]->representative_seqpos() == seqpos ) return i;
		}
		return 0;
	}

	Size
	ResidueAlternativeStepWiseSamplerComb::id_for_resnum( Size const seqpos ){
		runtime_assert( has_resnum( seqpos ) );
		return residue_alternative_rotamer_map_[ seqpos ]->id();
	}


	Residue const &
	ResidueAlternativeStepWiseSamplerComb::get_residue_at_origin( Size const seqpos ){
		return residue_alternative_rotamer_map_[ seqpos ]->get_residue_at_origin();
	}

	Residue const &
	ResidueAlternativeStepWiseSamplerComb::get_residue_at_origin_with_matching_type( Size const seqpos, Residue const & rsd_in ){
		return residue_alternative_rotamer_map_[ seqpos ]->get_residue_at_origin_with_matching_type( rsd_in );
	}

	// fast-forward inner loops so that we can get to the next residue pair.
	void
	ResidueAlternativeStepWiseSamplerComb::fast_forward_to_next_residue_pair( Size const i, Size const j){
		Size const which_rotamer_i = find_resnum( i );
		Size const which_rotamer_j = find_resnum( j );
		//		TR << "fast_forward_to_next_residue_pair " << i << "--" << j << ": " << which_rotamer_i << "--" << which_rotamer_j  << std::endl;
		if ( which_rotamer_i == 0 && which_rotamer_j == 0 ){
			fast_forward(); // fast forward to next rigid body.
		} else if ( which_rotamer_i == 0 ){
			fast_forward( which_rotamer_j - 1 );
		} else if ( which_rotamer_j == 0 ){
			fast_forward( which_rotamer_i - 1 );
		} else {
			Size which_rotamer = ( which_rotamer_i < which_rotamer_j ) ? which_rotamer_i : which_rotamer_j; // inner-most loop
			runtime_assert( which_rotamer > 0 );
			if ( which_rotamer > 1 ) 	fast_forward( which_rotamer - 1 ); // next ++ will give us another residue pair.
		}
	}

	// fast-forward inner loops so that we can get to the next residue pair.
	void
	ResidueAlternativeStepWiseSamplerComb::fast_forward_to_next_residue( Size const i ){
		Size const which_rotamer = find_resnum( i );
		//		TR << "fast_forward_to_next_residue " << i << ": " << which_rotamer << std::endl;
		if ( which_rotamer == 0 ){
			fast_forward(); // fast forward to next rigid body.
		} else if ( which_rotamer > 1 ) {
			fast_forward( which_rotamer - 1 ); // next ++ will give us another residue pair.
		}
	}

} //copy_dofs
} //sampler
} //stepwise
} //protocols
