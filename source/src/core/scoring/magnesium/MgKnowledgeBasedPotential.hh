// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/magnesium/MgKnowledgeBasedPotential.hh
/// @brief  Information on which atoms to use for computing clashes
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_MgKnowledgeBasedPotential_hh
#define INCLUDED_core_scoring_rna_MgKnowledgeBasedPotential_hh

// Unit Headers
#include <core/scoring/magnesium/MgKnowledgeBasedPotential.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace scoring {
namespace magnesium {

class MgKnowledgeBasedPotential : public utility::pointer::ReferenceCount {

public:

	/// @brief ctor, reads data file
	MgKnowledgeBasedPotential();

	core::chemical::rna::GaussianParameter
	get_mg_potential_gaussian_parameter( core::conformation::Residue const & rsd, Size const j, bool & is_phosphate_oxygen ) const;

	core::chemical::rna::GaussianParameter
	get_mg_potential_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const;

	core::chemical::rna::GaussianParameter
	get_mg_potential_indirect_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const;

	void setup_info_for_mg_calculation( core::pose::Pose & pose ) const;

	core::chemical::rna::GaussianParameter
	get_mg_potential_costheta_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const;

	core::chemical::rna::GaussianParameter
	get_mg_potential_costheta_indirect_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const;

	core::Real
	v_angle_width() const { return v_angle_width_; }

private: //data

	core::chemical::rna::GaussianParameter const gaussian_parameter_phosphate_oxygen_;
	core::chemical::rna::GaussianParameter const gaussian_parameter_imine_           ;
	core::chemical::rna::GaussianParameter const gaussian_parameter_exocyclic_oxygen_;
	core::chemical::rna::GaussianParameter const gaussian_parameter_o2prime_         ;
	core::chemical::rna::GaussianParameter const gaussian_parameter_water_oxygen_    ;

	core::chemical::rna::GaussianParameter const gaussian_parameter_phosphate_p_     ;
	core::chemical::rna::GaussianParameter const gaussian_parameter_polar_H_         ;
	core::chemical::rna::GaussianParameter const gaussian_parameter_nonpolar_H_         ;

	core::chemical::rna::GaussianParameter const gaussian_parameter_phosphate_oxygen_indirect_;
	core::chemical::rna::GaussianParameter const gaussian_parameter_imine_indirect_           ;
	core::chemical::rna::GaussianParameter const gaussian_parameter_exocyclic_oxygen_indirect_;
	core::chemical::rna::GaussianParameter const gaussian_parameter_o2prime_indirect_          ;

	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_phosphate_oxygen_;
	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_imine_           ;
	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_exocyclic_oxygen_;
	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_o2prime_         ;
	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_water_oxygen_    ;
	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_polar_H_         ;
	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_nonpolar_H_      ;

	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_phosphate_oxygen_indirect_;
	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_imine_indirect_           ;
	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_exocyclic_oxygen_indirect_;
	core::chemical::rna::GaussianParameter const gaussian_parameter_costheta_o2prime_indirect_         ;

	Real const v_angle_width_;
};

} //magnesium
} //scoring
} //core

#endif // INCLUDED_core_scoring_ScoreFunction_HH
