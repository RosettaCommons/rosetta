// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/scoring/RsdEnergiesBreakdown.fwd.hh
///
///	@brief		Residue Energies Breakdown
///	@details	Returns individual one body, two-body and whole strucutre energies after a pose is
///				sored. Useful for unit testing. Rocco wrote an awesome app for this in
///				public/analysis/residue_energy_breakdown.cc. I just made it into a class.
///				Last Modified: 4/24/14
///
///	@author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_RsdEnergiesBreakdown_fwd_hh
#define INCLUDED_core_scoring_RsdEnergiesBreakdown_fwd_hh

// Unit Headers
#include <core/scoring/RsdEnergiesBreakdown.fwd.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>

// Utility headers


// C++ headers
#include <cstdlib>

namespace core {
namespace scoring {

	/// @brief Constructor
	RsdEnergiesBreakdown::RsdEnergiesBreakdown(
		core::pose::PoseOP pose,
		core::scoring:ScoreFunctionOP sfxn
	) :
		utility::pointer::ReferenceCount(),
		pose_( pose ),
		sfxn_( sfxn )
	{
	

	}
	
	/// @brief Destructor
	RsdEnergiesBreakdown::~RsdEnergiesBreakdown() {}
	
	/// @brief Get One-Body Energies for ScoreType
	core::Real one_body( ScoreType type, core::Size seqpos ) {
		
	}
	
	/// @brief Get Two-Body Energies for ScoreType
	core::Real two_body( ScoreType type, core::Size seqpos1, core::Size seqpos2 ) {
	
	}
	
	
	/// @brief Get SUm Two-Body
	core::Real two_body_sum( ScoreType type ) {
		
	}
	
	/// @brief Get Sum One-Body
	core::Real one_body_sum( ScoreType type ) {
		
	}
	
	/// @brief Initialize Pair Energies Map
	void initialize_energies( core::pose::PoseOP pose, ScoreFunctionOP sfxn ) {
		
		using namespace core::scoring;
		
		// Get ScoreTypes and Weights
		ScoreTypes const scoretypes( scorefxn->get_nonzero_weighted_scoretypes() );
		EnergyMap weights( scorefxn->weights() );
		
		// Score Pose with the current scoring function
		(*sfxn)(*pose);
		Energies const & pose_energies( pose.energies() );
		EnergyGraph const & egraph( pose_energies.energy_graph() );
		
		// Compute One_Body Energies
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		
			EnergyNode const * node(egraph.get_energy_node(ii));
			runtime_assert( node != 0 );
			EnergyMap unwt_residue1b;
		
			scorefxn->eval_ci_1b(current_pose.residue(ii), current_pose, unwt_residue1b);
			scorefxn->eval_cd_1b(current_pose.residue(ii), current_pose, unwt_residue1b);
			scorefxn->eval_intrares_energy(current_pose.residue(ii), current_pose, unwt_residue1b);
			
			EnergyMap residue1b( unwt_residue1b * weights );
			
			// Create a hashmap of resnum => energies
			std::map< core::Size, core::Real > residue1b[ type ];
			for ( ScoreTypes::const_iterator iter( scoretypes.begin() ); iter != scoretypes.end(); ++iter ) {
				ss->add_energy(  name_from_score_type( *iter ), residue1b[ *iter ] );
			} // for non-zero score types
			ss->add_energy( "total", residue1b.sum() );
			sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
		}
		
		
		// Short Range
		EnergyEdge const * edge(egraph.find_energy_edge(ii, jj));
		if ( edge != 0 ) {
			EnergyMap unwt_pair_energies( edge->fill_energy_map() );
			pair_energies += unwt_pair_energies * weights;
			output = true;
		}
		
		// Long Range
		for( Size lr = 1; lr <= core::scoring::methods::n_long_range_types; lr++){
			LREnergyContainerCOP lrec = pose_energies.long_range_container( core::scoring::methods::LongRangeEnergyType( lr ) );
			if( !lrec || lrec->empty()) { continue; }
			for ( ResidueNeighborConstIteratorOP rni( lrec->const_upper_neighbor_iterator_begin( ii ) ),
				 end( lrec->const_upper_neighbor_iterator_end( ii ) ); *rni != *end; ++(*rni) ) {
				if( rni->upper_neighbor_id() != jj ) { continue; }
				rni->accumulate_energy( pair_energies );
				output = true;
			}
		}
		
		//Output
		if( !output ) { continue; }
		
		SilentStructOP ss( new ScoreFileSilentStruct );
		ss->decoy_tag( tag + "_" + string_of(ii) + "_" + string_of(jj) );
		ss->add_string_value( "pose_id", tag );
		ss->add_string_value( "resi1", string_of(ii) );
		ss->add_string_value( "restype1", current_pose.residue_type(ii).name3() );
		ss->add_string_value( "resi2", string_of(jj) );
		ss->add_string_value( "restype2", current_pose.residue_type(jj).name3() );
		
		bool nonzero(false);
		for ( ScoreTypes::const_iterator iter( scoretypes.begin() ); iter != scoretypes.end(); ++iter ) {
			ss->add_energy(  name_from_score_type( *iter ), pair_energies[ *iter ] );
			if ( pair_energies[ *iter ] <= -0.001 || 0.001 <= pair_energies[ *iter ] ) {
				nonzero = true;
			}
		} // for non-zero score types
		ss->add_energy( "total", pair_energies.sum() );
		if ( pair_energies.sum() <= -0.001 || 0.001 <= pair_energies.sum() ) {
			nonzero = true;
		}
		if (nonzero) { // Ignore pairs which would print as zero for all components
			sfd.write_silent_struct( *ss, option[ out::file::silent ]() );
		}
	}
	

} // scoring
} // core

#endif // ICNLUDED_core_scoring_RsdEnergiesBreakdown_cc

