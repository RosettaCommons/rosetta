// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   	core/membrane/scoring/MembraneSearch.hh
///
/// @brief  	Membrane Search
/// @detail		Performs a search method for optimal membrane normal and center from scoring
///				this is used in both conformation and scoring
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_scoring_MembraneSearch_hh
#define INCLUDED_core_membrane_scoring_MembraneSearch_hh

// Unit Headers
#include <core/membrane/scoring/MembraneSearch.fwd.hh>

// Project Headers
#include <core/membrane/scoring/MembranePenalties.hh>
#include <core/membrane/scoring/MembraneScoring.hh>

#include <core/membrane/util/definitions.hh>
#include <core/membrane/util/definitions_util.hh>

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

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>

namespace core {
namespace membrane {
namespace scoring {

/// @brief  	Membrane Search
/// @detail		Performs a search method for optimal membrane normal and center from scoring
///				this is used in both conformation and scoring
class MembraneSearch : utility::pointer::ReferenceCount {

public: // methods

	/// @brief 	Standard Constructor
	/// @detail	instantiate a membrane search object
	///			always internally instantiates a membrane penalties
	///			and scoring methods object
	///
	/// @param 	desc
    ///             resource description
    /// @param  chain
    ///             chain identifier
	/// @return MembraneSearch Object
	MembraneSearch( std::string desc, std::string chain );
    
    /// @brief Get Membrane Config
    /// @details Returns configuration for membrane spanning
    ///
    /// @throws Membrane Exception
    core::membrane::util::EmbedConfigInfoOP getConfig();

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
	/// @return - Returns the embedding config info object stored here
	void
    search_normal_center(
                         core::pose::Pose & pose,
                         core::Vector normal,
                         core::Vector center
                         );
    
private: // methods

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
	void compute_embed_config(
		core::Vector normal,
		core::Vector center,
		bool spanning
		);
    
    /// @brief Score Normal and Center Pair
    /// @details Score membrane normal and center with respect to membrane
    ///
    /// @throws <none>
    void
	score_normal_center(
                        pose::Pose const & pose,
                        Vector const & normal,
                        Vector const & center,
                        Real & score
                        ) const;

	/// @brief Search Membrane Normal
	/// @detail
	///
	/// @param vector const v
	void
	search_memb_normal( Vector & v, Real const & alpha, Real const & theta) const;

	/// @brief Search Membrane Center
	/// @detail
	///
	/// @param vector center
	/// @param vector normal
	void
	search_memb_center(Vector & c, Vector & n, Real const & delta) const;

	/// @brief Rotate and perturb vector
	/// @detail
	///
	/// @param vector to perturb
	void
	rot_perturb_vector(Vector & v, Real const & std_dev) const;

	/// @brief Rigid perturb vector
	void
	rigid_perturb_vector(Vector & v, Real const & std_dev) const;


private: // data
    
    // Store job description and chain id
    std::string desc_;
    std::string chain_;

	// Membrane Scoring Evals
	MembraneScoringOP scoring_;

	// Membrane Penalties Methods
	MembranePenaltiesOP penalty_;
    
    // Store a Config Object
    core::membrane::util::EmbedConfigInfoOP config_;

}; // class MembraneSearch

} // scoring
} // membrane
} // core



#endif // INCLUDED_core_membrane_scoring_MembraneSearch_hh