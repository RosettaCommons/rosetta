// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/carbohydrates/CHIEnergyFunction.hh
/// @brief   Method declarations for CHIEnergyFunction.
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_scoring_carbohydrates_CHIEnergyFunction_HH
#define INCLUDED_core_scoring_carbohydrates_CHIEnergyFunction_HH

// Unit header
#include <core/scoring/carbohydrates/CHIEnergyFunction.fwd.hh>
#include <core/chemical/carbohydrates/LinkageType.hh>

// Project Headers
#include <core/types.hh>
#include <core/id/types.hh>

// Utity header
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <map>


namespace core {
namespace scoring {
namespace carbohydrates {


	///@brief Struct for CHI-based dihedral sampling.
	struct CHIDihedralSamplingData{
		chemical::carbohydrates::LinkageType linkage_type;
		Real step_size;
		
		utility::vector1< Real > probabilities;
		utility::vector1< Angle > angles;
	};


/// @details  This is an implementation of the "CarboHydrate Intrinsic" (CHI) energy function developed by Woods Lab.\n
/// The Gaussian parameters for the function depend on whether the glycosidic bond in question is a phi or psi angle.\n
/// The parameters further depend on if the phi angles are at alpha or beta linkages and on if the psi angles are at
/// ->2-axial, ->3-equatorial, or ->4-axial OR ->2-equatorial, ->3-axial, or ->4-equatorial linkages.\n
/// The function has not been developed for ->6 linkages (with omega angles).
/// @ref      A.K. Nivedha et al. J. Comput. Chem. 2014, 35, 526-39
class CHIEnergyFunction : public utility::pointer::ReferenceCount {
public:  // Standard Methods //////////////////////////////////////////////////
	/// @brief  Default constructor
	CHIEnergyFunction();

	virtual ~CHIEnergyFunction();


public:  // Other Public Methods //////////////////////////////////////////////
	Energy operator()( chemical::carbohydrates::LinkageType type, Angle x ) const;

	Real evaluate_derivative( chemical::carbohydrates::LinkageType type, Angle x ) const;

public: // Dihedral Sampling Methods //////////////////////////////////////////

	//CHIDihedralSamplingData const &
	//get_chi_sampling_data(id::MainchainTorsionType torsion, LikageType type);
	
	CHIDihedralSamplingData const &
	get_chi_sampling_data(chemical::carbohydrates::LinkageType type) const;
	
	///@brief Sets up CHI sampling data
	void
	setup_for_sampling(core::Real step_size = 0.1);
	
	bool
	sampling_data_setup() const;
	
	bool
	sampling_data_setup( chemical::carbohydrates::LinkageType type ) const;
	
private:  // Private methods //////////////////////////////////////////////////
	// Return single CHI energy function term, ae^-((x-b)^2/c), for the given type and index.
	Energy evaluate_term( chemical::carbohydrates::LinkageType type, uint i, Angle x ) const;

	// Sum the individual terms.
	Energy evaluate_function( chemical::carbohydrates::LinkageType type, Angle x ) const;

	void init();


private:  // Private Data /////////////////////////////////////////////////////
	// Gaussian function parameters for the CHI energy function as defined in:
	// A.K. Nivedha et al. J. Comput. Chem. 2014, 35, 526-39.
	// The outer vector is indexed by the chemical::carbohydrates::LinkageType enum.
	utility::vector1< utility::vector1 < Real > > a_;  // magnitude of the Gaussian distribution
	utility::vector1< utility::vector1 < Real > > b_;  // midpoint of the Gaussian distribution
	utility::vector1< utility::vector1 < Real > > c_;  // twice the square of the width of the Gaussian distribution
	utility::vector1< Real > d_;  // the intercept (coefficient of the zeroth term)
	
private:
	//Formatted Sampling Data
	
	std::map< chemical::carbohydrates::LinkageType, CHIDihedralSamplingData > dihedral_sampling_data_;
};  // class CHIEnergyFunction


// Helper methods /////////////////////////////////////////////////////////////
// This allows one to use a for loop with LinkageType enum values.
chemical::carbohydrates::LinkageType & operator++( chemical::carbohydrates::LinkageType & type );


}  // namespace carbohydrates
}  // namespace scoring
}  // namespace core

#endif // INCLUDED_core_scoring_carbohydrates_CHIEnergyFunction_HH
