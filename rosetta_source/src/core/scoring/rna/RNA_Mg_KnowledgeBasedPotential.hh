// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_Mg_KnowledgeBasedPotential.hh
/// @brief  Information on which atoms to use for computing clashes
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_Mg_KnowledgeBasedPotential_hh
#define INCLUDED_core_scoring_rna_RNA_Mg_KnowledgeBasedPotential_hh

// Unit Headers
#include <core/scoring/rna/RNA_Mg_KnowledgeBasedPotential.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/rna/Gaussian_parameter.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace scoring {
namespace rna {

class RNA_Mg_KnowledgeBasedPotential : public utility::pointer::ReferenceCount {

public:

	/// @brief ctor, reads data file
	RNA_Mg_KnowledgeBasedPotential();

	Gaussian_parameter
	get_mg_potential_gaussian_parameter( core::conformation::Residue const & rsd, Size const j, bool & is_phosphate_oxygen ) const;

	Gaussian_parameter
	get_mg_potential_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const;

	Gaussian_parameter
	get_mg_potential_indirect_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const;

	void setup_info_for_mg_calculation( core::pose::Pose & pose ) const;

	Gaussian_parameter
	get_mg_potential_costheta_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const;

	Gaussian_parameter
	get_mg_potential_costheta_indirect_gaussian_parameter( core::conformation::Residue const & rsd, Size const j ) const;

private: //data

	Gaussian_parameter const gaussian_parameter_phosphate_oxygen_;
	Gaussian_parameter const gaussian_parameter_imine_           ;
	Gaussian_parameter const gaussian_parameter_exocyclic_oxygen_;
	Gaussian_parameter const gaussian_parameter_o2star_          ;

	Gaussian_parameter const gaussian_parameter_phosphate_p_     ;
	Gaussian_parameter const gaussian_parameter_polar_H_         ;
	Gaussian_parameter const gaussian_parameter_nonpolar_H_         ;

	Gaussian_parameter const gaussian_parameter_phosphate_oxygen_indirect_;
	Gaussian_parameter const gaussian_parameter_imine_indirect_           ;
	Gaussian_parameter const gaussian_parameter_exocyclic_oxygen_indirect_;
	Gaussian_parameter const gaussian_parameter_o2star_indirect_          ;

	Gaussian_parameter const gaussian_parameter_costheta_phosphate_oxygen_;
	Gaussian_parameter const gaussian_parameter_costheta_imine_           ;
	Gaussian_parameter const gaussian_parameter_costheta_exocyclic_oxygen_;
	Gaussian_parameter const gaussian_parameter_costheta_o2star_          ;
	Gaussian_parameter const gaussian_parameter_costheta_polar_H_          ;
	Gaussian_parameter const gaussian_parameter_costheta_nonpolar_H_          ;

	Gaussian_parameter const gaussian_parameter_costheta_phosphate_oxygen_indirect_;
	Gaussian_parameter const gaussian_parameter_costheta_imine_indirect_           ;
	Gaussian_parameter const gaussian_parameter_costheta_exocyclic_oxygen_indirect_;
	Gaussian_parameter const gaussian_parameter_costheta_o2star_indirect_          ;

	};



}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
