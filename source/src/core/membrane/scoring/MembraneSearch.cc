// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   	core/membrane/scoring/MembraneSearch.cc
///
/// @brief  	Membrane Search
/// @detail		Performs a search method for optimal membrane normal and center from scoring
///				this is used in both conformation and scoring
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_scoring_MembraneSearch_cc
#define INCLUDED_core_membrane_scoring_MembraneSearch_cc

// Unit Headers
#include <core/membrane/scoring/MembraneSearch.hh>

// Project Headers
#include <core/membrane/scoring/MembranePenalties.hh>
#include <core/membrane/scoring/MembraneScoring.hh>

#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/Exceptions.hh>

// Package Headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>

#include <core/pose/Pose.hh>

#include <core/types.hh>

#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

#include <basic/Tracer.hh>

#include <basic/resource_manager/ResourceManager.hh>
#include <basic/resource_manager/util.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.membrane.scoring.MembraneSearch");
static numeric::random::RandomGenerator RG(280628);  // <- Magic number, do not change it! <- :P bad magic number

namespace core {
namespace membrane {
namespace scoring {

	/// @brief 	Empty Constructor
	/// @detail	instantiate a membrane search object
	///			always internally instantiates a membrane penalties
	///			and scoring methods object
	///
	/// @param 	desc
    ///             resource description
    /// @param  chain
    ///             chain identifier
	/// @return MembraneSearch Object
    MembraneSearch::MembraneSearch(
           std::string desc,
           std::string chain
           )
           : utility::pointer::ReferenceCount(),
            desc_(desc),
            chain_(chain)
    {
        // Initialize
		penalty_ = new core::membrane::scoring::MembranePenalties();
		scoring_ = new core::membrane::scoring::MembraneScoring();
        config_ = core::membrane::util::init_embedConfigInfo();

	}

    /// @brief Get Membrane Config
    /// @details Returns configuration for membrane spanning
    ///
    /// @throws Membrane Exception
    core::membrane::util::EmbedConfigInfoOP
    MembraneSearch::getConfig() { return config_; }
    
	/// @brief  Compute Embed config
	/// @detail Store accepted parameters in an embedding configuraiton object
	///
	/// Precondition: Embedding config definition is somewhat initialized in the constructor
	///
	/// @param  normal
	///				normal to set
	/// @param 	center
	///				center to set
	/// @param  spanning
	///				spanning works correctly?
	///
	/// @return [none]
    void
    MembraneSearch::compute_embed_config( core::Vector normal, core::Vector center, bool spanning ) {

        using namespace core::membrane::util;
        
        config_->normal = normal;
        config_->center = center;
        config_->spanning = spanning;

	}

	/// @brief 	Compute Normal and Center
	/// @detail Computes Membrane Normal and Center via score and search
	///
	/// Precondition: Normal and center are already initialized to their desired
	///	initial values
	///
	/// @param pose
	///			pose of interest
	/// @param normal
	///			pre-initialized normal vector to search for
	/// @param center
	/// 		pre-initialized center vector to search for
	///
	/// @return - still deciding upon search type
	void
    MembraneSearch::search_normal_center(
                    core::pose::Pose & pose,
                    core::Vector normal,
                    core::Vector center
        )
	{
        using namespace core::membrane::util;
        using namespace core::membrane::geometry;
        
        // About to get used a lot
        MembraneBoundsChecking mbc;

        // Get Membrane Option Definitions and Resource Management
	 	SpanningTopologyOP topology_ = load_topology(desc_, chain_);
	 	EmbedSearchParamsOP embed_params = load_embed_params(desc_, chain_);

	 	// Initialize Some starting values
	 	Real score(0), best_score(999999), accepted_score(999999);
	 	Real temperature(2.0);

	 	// Initialize trial, original, best, and accepted vars
	 	Vector trial_normal(normal), trial_center(center);
	 	Vector orig_trial_normal(normal), orig_trial_center(center);
	 	Vector best_normal(normal), best_center(center);
	 	Vector accepted_normal(normal), accepted_center(center);

	 	// Score Normal and Center
	 	score_normal_center( pose, trial_normal, trial_center, best_score );
	 	accepted_score = best_score;

	 	// Confirm that normal and center align with spanning
	 	bool spanning = mbc.check_spanning( pose, trial_normal, trial_center, desc_, chain_ );

	 	// Set some local vars
	 	Real normal_mag(embed_params->normal_mag);
	 	Real center_mag(embed_params->center_mag);

	 	// They are using an int, I am going to use a size
	 	core::Size center_max_delta(embed_params->center_max_delta);
	 	core::Size alpha_start(embed_params->normal_start_angle);
	 	core::Size delta_alpha(embed_params->normal_delta_angle);
	 	core::Size max_alpha(embed_params->normal_max_angle);
        
        // Init Best SPanning Param
        bool best_spanning(false);

	 	// scoring vars
	 	core::Size nres(pose.total_residue());
	 	core::Size counter(0), accepted(0), thermally_accepted(0);

	 	// I am going to leave the symmetry param in here but there is something
	 	// that seriously bothers me about it...?? @@FIX
	 	if ( embed_params->center_search || core::pose::symmetry::is_symmetric(pose) ) {

	 		for ( int delta_center = -center_max_delta; delta_center <= center_max_delta; ++delta_center ) {

	 			// Search for potential centers
	 			trial_center = orig_trial_center;
	 			search_memb_center( trial_center, trial_normal, delta_center);

	 			// Check against spanning - if it violates the spanning move back
	 			if ( !mbc.check_spanning( pose, trial_normal, trial_center, desc_, chain_ ) ) {

	 				trial_center = accepted_center;
	 				trial_normal = accepted_normal;
	 				continue;
	 			}

	 			// Score normal center
	 			score_normal_center( pose, trial_normal, trial_center, score );

	 			// If score is lower than accepted score and best score, reassign
	 			if ( score < accepted_score ) {

	 				if ( score < best_score ) {
	 					best_score = score;
						best_center = trial_center;
						best_normal = trial_normal;
						best_spanning = true;
	 				}

	 				accepted_score = score;
	 				accepted_center = trial_center;
	 				accepted_normal = trial_normal;
	 			}
	 		}

	 		// Save best projection
	 		compute_embed_config(best_normal, best_center, best_spanning);
	 	}

	 	// Search for membrane normal
	 	if ( embed_params->normal_search ) {

	 		for( Size alpha = alpha_start; alpha <= max_alpha; alpha += delta_alpha ) {
				for( Size theta = 0; theta < 360; theta += 60 ) {

					// Search for membrane normal
					trial_normal=orig_trial_normal;
					search_memb_normal(trial_normal, alpha, theta);

					// Check the normal works with spanning topology
					if(! mbc.check_spanning(pose, trial_normal, trial_center, desc_, chain_)){
						trial_center=accepted_center;
						trial_normal=accepted_normal;
						continue;
					}

					// Score normal and center
					score_normal_center(pose,trial_normal,trial_center,score);

						if(score<accepted_score) {

							if(score<best_score) {

								best_score=score;
								best_center=trial_center;
								best_normal=trial_normal;
								best_spanning=true;
							}

							accepted_score=score;
							accepted_center=trial_center;
							accepted_normal=trial_normal;
						}
					}
				}

			// Accept
			compute_embed_config( best_normal, best_center, best_spanning );

	 	}


	 	// If no membrane center search or normal search, optimize???? (WTF Is this!!!!)
		if ( !embed_params->center_search && !embed_params->normal_search ) {

			for( Size cycles = 1; cycles <= embed_params->normal_cycles; ++cycles ) {

				// Set initial temperature for search
				temperature = 2.0/cycles;

				// Randomnly search
				if(RG.uniform()<0.5) {
					rigid_perturb_vector(trial_center,center_mag);
				} else {
					rot_perturb_vector(trial_normal,normal_mag);
				}

				// Check that everything is ok with spanning
				if(!mbc.check_spanning( pose, trial_normal, trial_center, desc_, chain_) ) {

					trial_center=accepted_center;
					trial_normal=accepted_normal;
					continue;
				}

				// Score mmebrane normal and center
				score_normal_center(pose,trial_normal,trial_center,score);

				// Score is better, accept the score
				if(score<accepted_score) {

					if(score<best_score) {

						best_score = score;
						best_center = trial_center;
						best_normal = trial_normal;
						best_spanning = true; //bw if you are here it is spanning....
					}

					accepted_score = score;
					accepted_center = trial_center;
					accepted_normal = trial_normal;
					++accepted;

				// Cool - I would like to add some better comments for this
				// because I think this if block does some important stuff
				} else {

					++counter;
					Real const boltz_factor=(accepted_score-score)/temperature;
					Real const probability = std::exp( std::min ((core::Real)40.0, std::max((core::Real)-40.0,boltz_factor)) );

					if(RG.uniform()<probability)
					{
						accepted_score=score;
						accepted_center=trial_center;
						accepted_normal=trial_normal;
						++thermally_accepted;
						++accepted;

					} else {

						trial_center=accepted_center;
					    trial_normal=accepted_normal;
					}
				}
			}


			// Save the final best projection
			compute_embed_config( best_normal, best_center, best_spanning );
		}
	}

	/// @brief Search Membrane Normal
	/// @details
	///
	/// @param vector const v
	void
	MembraneSearch::search_memb_normal( Vector & v, Real const & alpha, Real const & theta) const
	{
		Real r_alpha = numeric::conversions::radians(alpha);
		Real r_theta = numeric::conversions::radians(theta);
		Vector u(std::sin(r_alpha) * std::cos(r_theta), std::sin(r_alpha) * std::sin(r_theta), std::cos(r_alpha));
		v = rotation_matrix_degrees( u, alpha) * v;
	}

	/// @brief Search Membrane Center
	/// @detail
	///
	/// @param vector center
	/// @param vector normal
	void
	MembraneSearch::search_memb_center(Vector & c, Vector & n, Real const & delta) const
	{
		c = c + delta * n;
	}

	/// @brief Rotate and perturb vector
	/// @detail
	///
	/// @param vector to perturb
	void
	MembraneSearch::rot_perturb_vector(Vector & v, Real const & std_dev) const
	{
		Vector u( numeric::random::gaussian(), numeric::random::gaussian(), numeric::random::gaussian());
		Real alpha(numeric::random::gaussian()*std_dev);
		v = rotation_matrix(u, alpha) * v;
	}

	/// @brief Rigid perturb vector
	void
	MembraneSearch::rigid_perturb_vector(Vector & v, Real const & std_dev) const
	{
		Vector u(numeric::random::gaussian(),numeric::random::gaussian(),numeric::random::gaussian());
		u.normalize();
		v = v + std_dev * u;
	}

	/// @brief Score Normal and Center Vector
	void
	MembraneSearch::score_normal_center( pose::Pose const & pose,
                                         Vector const & normal,
                                         Vector const & center,
                                         Real & score
                                        ) const
	{

		using namespace core::conformation::symmetry;
		using namespace core::membrane::util;

		// Get total residue and membrane topology data
		Size const nres = pose.total_residue();
        SpanningTopologyOP topology = load_topology(desc_, chain_);

        // Initialize score to 0
        score = 0;
        
		// Initialize some variables for scoring
		Real residue_score(0), tm_projection(0), non_helix_pen(0), termini_pen(0);

		// Loop through the pose
		for ( Size i = 1; i <= nres; ++i ) {

			Size rsdSeq(i);

			// If the pose is symmetric,
			if ( core::pose::symmetry::is_symmetric( pose ) ) {

				SymmetricConformation const & symm_conf (
					dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
					SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

			if (! symm_info->bb_is_independent(pose.residue(i).seqpos()) ) {
				rsdSeq = symm_info->bb_follows(pose.residue(i).seqpos());
			}

			if ( symm_info->is_virtual(i) ) {
				rsdSeq = 0;
			}

		}

		if (rsdSeq ==0 ) continue;

		//CA coords
		if ( pose.residue(rsdSeq).aa() == core::chemical::aa_vrt ) continue;
		if(!topology->allow_scoring[rsdSeq]) continue;

		Vector const & xyz( pose.residue( i ).atom( 2 ).xyz());
		core::Real depth = dot(xyz-center,normal)+30; // thickness??
		scoring_->evaluate_env( pose, pose.residue(i), depth, residue_score );
		score+=residue_score;

        // Get embedding parameters
        EmbedSearchParamsOP params = load_embed_params( desc_, chain_ );

		// if user specifies to evlaute membrane environment penalties
		if( params->penalties ) {

			// Evaluate penalties
			penalty_->tm_projection_penalty( desc_, chain_, normal, center, tm_projection);
			penalty_->non_helix_in_membrane_penalty( desc_, chain_, normal, center, non_helix_pen );
			penalty_->termini_penalty( desc_, chain_, normal, center, termini_pen);
			score += tm_projection + non_helix_pen + termini_pen; // bw skipping term_penalty+50.0*term_penalty; //bw 0.5*c++ version.
		}
	}
    }
    
    /// @brief Load Required Resources
    /// @details Load required resources based on default description (startstruct base description)
    ///
    /// @return [none]
    /// @throws EXCN_Resource_Definition
    void load_required_resources() {
        
        // load in search parameters
        
        
    }
    
    /// @brief Modify Search resources
    /// @details Search for a specific chain by resource and not the whole pose
    void modify_required_resources() {
        
        // modify search parameters
    }

} // scoring
} // membrane
} // core


#endif // INCLUDED_core_membrane_scoring_MembraneSearch_cc