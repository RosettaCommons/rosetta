// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/PcsEnergyParameter.hh
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////


#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsEnergyParameter_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsEnergyParameter_hh

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers

// Objexx headers

// C++ headers
#include <string>

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

class PcsEnergyParameter{
public:

private:
	core::Size include_only_start_;
	core::Size include_only_end_;
	core::Size n_trial_min_;
	core::Real pcs_weight_;
	core::Real individual_scale_;
	utility::vector1<std::string> vec_filename_;
	utility::vector1<core::Real> vec_individual_weight_;

public:

	PcsEnergyParameter(); //Construct

	~PcsEnergyParameter(); //Destruct

	PcsEnergyParameter & // =
	operator=(PcsEnergyParameter const & other);

	PcsEnergyParameter(PcsEnergyParameter const & other); //copy


	/// @Set the vector an weight name for this PcsEnergyParameter
	void
	set_vector_name_and_weight(utility::vector1<std::string> const vec_filename,
		utility::vector1<core::Real> const vec_individual_weight);

	/// @Set the grid parameter for this PcsEnergyParameter
	void
	set_grid_param(
		core::Size const include_only_start,
		core::Size const include_only_end,
		core::Size const n_trial_min,
		core::Real const pcs_weight,
		core::Real const individual_scale
	);


	/// @brief Output myself on the stream
	friend std::ostream &
	operator<<(std::ostream& out, const PcsEnergyParameter &me);

	/// @brief
	core::Size
	get_include_only_start() const;

	/// @brief
	core::Size
	get_include_only_end() const;

	/// @brief
	core::Size
	get_n_trial_min() const;

	/*
	/// @brief Give me the gride_large_cutoff value
	core::Real
	get_grid_large_cutoff() const;

	/// @brief Give me the gride_cone_angle_cutoff value
	core::Real
	get_grid_cone_angle_cutoff() const;

	/// @brief Give me the gride_atom_name_1 value
	std::string
	get_grid_atom_name_1() const;

	/// @brief Give me the gride_atom_name_2 value
	std::string
	get_grid_atom_name_2() const;

	/// @brief Give me the gride_residue_num_1 value
	core::Size
	get_grid_residue_num_1() const;

	/// @brief Give me the gride_residue_num_2 value
	core::Size
	get_grid_residue_num_2() const;

	/// @brief Give me the gride_k_vector value
	core::Real
	get_grid_k_vector() const;

	/// @brief Give me the bool if the minimization of the best tensor is active
	bool
	get_minimize_best_tensor() const;
	*/
	/// @brief Give me the pcs_weight value
	core::Real
	get_pcs_weight() const;

	core::Real
	get_individual_scale() const;


	/// @brief Give me the vector of the name
	utility::vector1<std::string> const &
	get_vector_filename() const;

	/// @brief Give me the vector of the weight
	utility::vector1<core::Real> const &
	get_vector_weight() const;
};

} //PCS
} //methods
} //scoring
} //core

#endif
