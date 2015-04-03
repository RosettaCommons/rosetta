// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 /// @file protocols/scoring/methods/pcs2/PcsGridSearchParameter.hh
 ///
 /// @authorv Christophe Schmitz
 ///
 ////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsGridSearchParameter_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsGridSearchParameter_hh

// Unit Headers

// Package Headers

// Project Headers
#include <core/types.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{


class PcsGridSearchParameter{


public:

	core::Size include_only_start_stage1_;
	core::Size include_only_start_stage2_;
	core::Size include_only_start_stage3_;
	core::Size include_only_start_stage4_;

	core::Size include_only_end_stage1_;
	core::Size include_only_end_stage2_;
	core::Size include_only_end_stage3_;
	core::Size include_only_end_stage4_;

	core::Size n_trial_min_stage1_;
	core::Size n_trial_min_stage2_;
	core::Size n_trial_min_stage3_;
	core::Size n_trial_min_stage4_;

	core::Real pcs_weight_stage1_;
	core::Real pcs_weight_stage2_;
	core::Real pcs_weight_stage3_;
	core::Real pcs_weight_stage4_;

	core::Real individual_scale_stage1_;
	core::Real individual_scale_stage2_;
	core::Real individual_scale_stage3_;
	core::Real individual_scale_stage4_;

	utility::vector1<std::string> filenames_;
	utility::vector1<core::Real> individual_weights_;

public:
	PcsGridSearchParameter(); //Construct

	~PcsGridSearchParameter(); //Destruct

	PcsGridSearchParameter & // =
	operator=(PcsGridSearchParameter const & other);

	PcsGridSearchParameter(PcsGridSearchParameter const & other); //copy

	/// @brief Print me
	friend
	std::ostream &
	operator << ( std::ostream& out, const PcsGridSearchParameter &me );

	/// @brief This is to make sure the grid is correctly set up
	void
	control_grid_param();

};


}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
