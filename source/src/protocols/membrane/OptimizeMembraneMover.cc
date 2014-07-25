// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/OptimizeMembraneMover.cc
///
/// @brief      Monte Carlo Search for Membrane Position (based on pose)
/// @details	Uses monte carlo embedding search protocol described in the
///				old MembranePotential class described in Yarov-Yaravy et al.
///				Last Modified (6/16/14)
///
/// @author     Rebecca Alford (rfalford12@gmail.com)

/*

Notes: 

This optimization mover will be different than the one previously used in MembranePotential. It is going to address a few issues

 - radians vs. degree conversion mis happs occuring in search memb norm and rot matrix_degrees (no idea why this is)
  - uses score functions instead of introducing loops
  - should just set & trial (modified mc oh this is so exciting???)
  
  also - they did span checks and would throw out any membrane set that did not respect the spanning. this makes sense to me. could be a break step?
  
  lots of variables still to initialize here, but working on it...
  
  if center/normal search were false and fixed_membrane was false, they did this crazy full fledged inefficent random search. I think I'm going to see how the initial search protocol works but keep in mind I might want to extend to that later.

*/

#ifndef INCLUDED_protocols_membrane_OptimizeMembraneMover_cc
#define INCLUDED_protocols_membrane_OptimizeMembraneMover_cc

// Unit Headers
#include <protocols/membrane/OptimizeMembraneMover.hh>

// Project Headers
#include <protocols/moves/Mover.hh>

// Package Headers
#include <core/pose/Pose.hh> 
#include <core/types.hh> 

#include <core/scoring/ScoreFunction.hh> 
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/moves/MonteCarlo.hh>

// Utility Headers
#include <utility/vector1.hh> 
// ...cont'd...


namespace protocols {
namespace membrane {

////////////////////////
/// Constructors	 ///
////////////////////////

/// @brief Default Constructor
OptimizeMembraneMover::OptimizeMembraneMover();

/// @brief Custom Choose Normal or Center Constructor
OptimizeMembraneMover::OptimizeMembraneMover( bool search_center, bool search_normal );

/// @brief Copy Constructor
OptimizeMembraneMover::OptimizeMembraneMover( OptimizeMembraneMover & src );

/// @brief Virtual Destructor
OptimizeMembraneMover::~OptimizeMembraneMover();


/// @brief Register Options
void
OptimizeMembraneMover::register_options() {
	
}

/// @brief Initialize Options from Command Line
void
OptimizeMembraneMover::init_from_cmd() {
}

/// @brief Setup sfxn
void
OptimizeMembraneMover::setup_sfxns() {
	
	using namespace core::scoring;
	
	// Setup fullatom vs. centroid energy method
	if ( fullatom_ ) {
		sfxn_ = ScoreFunctionFactory::create_score_function( "cen_membrane_2014" );
	} else {
		sfxn_ = ScoreFunctionFactory::create_score_function( "fa_membrane_2014" );
	}
	
	// Add penalties if user specified (manually patching them in)
	if ( include_penalties_ ) {
	
		sfxn_->set_weight( MPTermini, 1.0 );
		sfxn_->set_weight( MPNonHelix, 1.0 );
		sfxn_->set_weight( MPTMProj, 1.0 );
	
	}
	
	// scorefunction setup!
}
		


//////////////////////
/// Mover Methods  ///
//////////////////////

/// @brief Alternate Optimize Mmebrane Protocol
////	 definitely more expensive search, possibly more correct?
void
OptimizeMembraneMover::apply( Pose & pose ) {
	
	using namespace core::scoring;
	using namespace protocols::moves;

	// Before you proceed, check that you are in fact a membrane pose
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Warning! Cannot optimize the membrane on a non membrane pose. Exiting from Optimize membrane mover" );
	}

	// Create new monte carlo object
	MonteCarlo mc( pose, *sfxn_, 1.0 ); // last arg is temp, can adjust?
	
	// Store the current center/normal
	Vector center = pose.conformation().membrane_center();
	Vector normal = pose.conformation().membrane_normal();
	
	// currently, overwriting normal and center too much. if mc was not accepted,
	// needs to revert to the old center/normal. will update it in the pose for us, but
	// not in the parameters held in the apply method.
	
	// Cycle through (center normal might not be in the right order)
	for ( Size i = 1; i <= ncycles_; ++i ) {
		
		// If allowed, search for a new center
		if ( center_search_ ) {
		
			for ( Size delta_= -max_delta_; delta <= max_delta_; ++delta ) {
				
				// Search center, update move, mc.accept?
				search_center( center, normal, delta );
				pose.conformation().update_embedding( center, normal );
				mc.boltzmann( pose );
				
				if ( !mc.accepted() ) {
					center = prev_center;
				}
				 
			}
		}
		
		// If allowed, search for a new normal
		if ( normal_search_ ) {
		
			// still goin gwith this is a ridiculous amount of sampling...???
			for ( Size alpha_ = alpha_start_; alpha_ <= max_alpha_; alpha_ += delta_alpha_ ) {
				for ( Size theta = 0; theta < 360; theta += 60 ) {
					
					// Search normal
					search_normal( normal, alpha, theta );
					pose.conformation().update_embedding( center, normal );
					mc.boltzman( pose );
					
					if ( !mc.accepted() ) {
						normal = prev_normal;
					}
				}
				
			}
			
		}
		
		
	}



	
}

/// @brief Obligatory get string method
std::string OptimizeMembraneMover::get_name() const { return "OptimizeMembraneMover"; }

/// @brief Score Relevant Center and Normal Parameters
void
OptimizeMembraneMover::score_normal_center(
	core::pose::Pose & pose,
	core::Vector & normal,
	core::Vector & center,
	core::Real & score
	) {
	
	using namespace core::scoring;
	
	// Grab total residue in the pose and topology by chain
	Size const nres = pose.total_residue();
	utility::vector1< SpanningTopology > topology = pose.conformation().membrane()->spanning_topology();
	
	// Create a new scoring function and turn on env weight
	ScoreFunctionOP sfxn = new ScoreFunctionOP;
	sfxn->set_weight( MPEnv, idk );
	
	// Include penalties in the energy function if specified by the user
	if ( include_penalties_ ) {
	
		sfxn->set_weight( MPTermini, idk );
		sfxn->set_weight( MPTMProj, idk );
		sfxn->set_weight( MPNonHelix, idk ); // well, this was untested
		// and maybe some others for fun?
		
	}
	
	// Update parameters, holding onto old parameters
	// maybe we should actually use mc instead of just a straight score comparison
	
	// Initialize Start Score Parameters (which I think will be replaced by an efunc
		
				
	// Initialize Scoring Parameters
	score = 0;
	Real residue_score(0);
	Real tm_projection(0);
	Real non_helix_pen(0);
	Real termini_pen(0);
	
	// For every residue in the pose, evaluate membrane environment based on center and normal given. Makes corrections
	// for symmetry
	for ( Size i = 1; i <= nres; ++i ) {
		
		Size rsdSeq(i);

		
		// If any of these conditions apply, skip to the end of the loop
		if (rsdSeq == 0 ) continue;
		if ( pose.residue(rsdSeq).aa() == core::chemical::aa_vrt ) continue;
		if(!topology.allow_scoring(rsdSeq)) continue;
		
		// Grab CA coords, compute depth, score based on env and append score
		Vector const & xyz( pose.residue( i ).atom( 2 ).xyz());
		core::Real depth = dot( xyz-center, normal ) + 30;
		evaluate_env( pose, pose.residue(i), depth, residue_score );
		score+=residue_score;
		
	}
	
	// if the user specified to apply mp penalties, append the score
	if(Menv_penalties_) {
		tm_projection_penalty( pose, normal, center, tm_projection );
		non_helix_in_membrane_penalty( pose, normal, center, non_helix_pen );
		termini_penalty( pose, normal, center, termini_pen );
		score+=tm_projection+non_helix_pen+termini_pen;
	}
}
	
}

/// @brief Search for Membrane Normal thorugh Random Perturbation
void
OptimizeMembraneMover::search_memb_normal(
	core::Vector & normal,
	core::Real & alpha,
	core::Real & theta
) {

	using namespace std;
	using namespace numeric::conversions;
	
	// Convert to radians
	Real r_alpha = radians( alpha );
	Real r_theta = radians( theta );
	
	// Compute the dcos vector
	Vector u( sin( r_alpha ) * cos( r_theta ),
			  sin( r_alpha ) * sin( r_theta ),
			  cos( r_alpha )
			 );
	
	// Rotate and apply
	n = rotation_matrix_degrees( u, alpha ) * n;
	
}

/// @brief Search for Membrane Center through Random Perturbation
void
OptimizeMembraneMover::search_memb_center(
	core::Vector & center,
	core::Vector & normal,
	core::Real & delta
) {
	c = c + delta*n;
}

/// @brief Randomnly Rotate and Perturb Vector
void
OptimizeMembraneMover::rot_perturb_vector( core::Vector & v, Real & std_dev ) {

	using namespace numeric::random;
	
	Vector u( gaussian(), gaussian(), gaussian() );
	Real alpha( gaussian() * std_dev );
	v = rotation_matrix( u, alpha ) * v;
	
}

/// @brief Randomnly Perturb Vector (Rigid Random Gaussian Perturbation)
void
OptimizeMembraneMover::rigid_pertrub_vector( core::Vector & v, Real & std_dev ) {

	using namespace numeric::random;
	
	Vector u( gaussian(), gaussian(), gaussian() );
	u.normalize();
	v = v * std_dev * u;

}
	
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_OptimizeMembraneMover_cc