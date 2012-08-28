// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file   core/scoring/rna/RNA_FittedTorsionInfo.hh
/// @brief  Statistically derived torsion information for RNA
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_RNA_FittedTorsionInfo_HH
#define INCLUDED_core_scoring_rna_RNA_FittedTorsionInfo_HH

#include <core/types.hh>

// Project headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++

namespace core {
namespace scoring {
namespace rna {

enum _RNA_FittedTorsionInfo_ { WHATEVER, ALPHA, BETA, GAMMA, DELTA, EPSILON, ZETA, CHI, NU2, NU1, O2H};

class Gaussian_parameter {
public:
	Real amplitude, center, width;

	Gaussian_parameter ( Real const amplitude_in, Real const center_in, Real const width_in ):
		amplitude( amplitude_in ),
		center   ( center_in ),
		width    ( width_in )
	{}

	//
	Gaussian_parameter &
	operator=( Gaussian_parameter const & src )
	{
		amplitude = src.amplitude;
		center = src.center;
		width = src.width;
		return *this;
	}


};

typedef utility::vector1< Gaussian_parameter > Gaussian_parameter_set;

class RNA_FittedTorsionInfo : public utility::pointer::ReferenceCount  {

public:

	RNA_FittedTorsionInfo();

	virtual ~RNA_FittedTorsionInfo();

	Real delta_cutoff() const { return DELTA_CUTOFF_; }

	Gaussian_parameter_set gaussian_parameter_set_alpha() const{ return gaussian_parameter_set_alpha_; }
	Gaussian_parameter_set gaussian_parameter_set_beta() const{ return gaussian_parameter_set_beta_; }
	Gaussian_parameter_set gaussian_parameter_set_gamma() const{ return gaussian_parameter_set_gamma_; }
	Gaussian_parameter_set gaussian_parameter_set_delta_north() const{ return gaussian_parameter_set_delta_north_; }
	Gaussian_parameter_set gaussian_parameter_set_delta_south() const{ return gaussian_parameter_set_delta_south_; }
	Gaussian_parameter_set gaussian_parameter_set_epsilon_north() const{ return gaussian_parameter_set_epsilon_north_; }
	Gaussian_parameter_set gaussian_parameter_set_epsilon_south() const{ return gaussian_parameter_set_epsilon_south_; }
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_sc_minus() const{ return gaussian_parameter_set_zeta_alpha_sc_minus_; }
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_sc_plus() const{ return gaussian_parameter_set_zeta_alpha_sc_plus_; }
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_ap() const{ return gaussian_parameter_set_zeta_alpha_ap_; }
	Gaussian_parameter_set gaussian_parameter_set_chi_north() const{ return gaussian_parameter_set_chi_north_; }
	Gaussian_parameter_set gaussian_parameter_set_chi_south() const{ return gaussian_parameter_set_chi_south_; }
	Gaussian_parameter_set gaussian_parameter_set_nu2_north() const{ return gaussian_parameter_set_nu2_north_; }
	Gaussian_parameter_set gaussian_parameter_set_nu2_south() const{ return gaussian_parameter_set_nu2_south_; }
	Gaussian_parameter_set gaussian_parameter_set_nu1_north() const{ return gaussian_parameter_set_nu1_north_; }
	Gaussian_parameter_set gaussian_parameter_set_nu1_south() const{ return gaussian_parameter_set_nu1_south_; }

	Real ideal_delta_north() const{ return ideal_delta_north_; }
	Real ideal_nu2_north() const{ return ideal_nu2_north_; }
	Real ideal_nu1_north() const{ return ideal_nu1_north_; }

	Real ideal_delta_south() const{ return ideal_delta_south_; }
	Real ideal_nu2_south() const{ return ideal_nu2_south_; }
	Real ideal_nu1_south() const{ return ideal_nu1_south_; }

private:


	void
	init_rna_torsion_gaussian_parameters();

	bool rna_tight_torsions_;

	Gaussian_parameter_set gaussian_parameter_set_alpha_;
	Gaussian_parameter_set gaussian_parameter_set_beta_;
	Gaussian_parameter_set gaussian_parameter_set_gamma_;
	Gaussian_parameter_set gaussian_parameter_set_delta_north_;
	Gaussian_parameter_set gaussian_parameter_set_delta_south_;
	Gaussian_parameter_set gaussian_parameter_set_epsilon_north_;
	Gaussian_parameter_set gaussian_parameter_set_epsilon_south_;
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_sc_minus_;
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_sc_plus_;
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_ap_;
	Gaussian_parameter_set gaussian_parameter_set_chi_north_;
	Gaussian_parameter_set gaussian_parameter_set_chi_south_;
	Gaussian_parameter_set gaussian_parameter_set_nu2_north_;
	Gaussian_parameter_set gaussian_parameter_set_nu2_south_;
	Gaussian_parameter_set gaussian_parameter_set_nu1_north_;
	Gaussian_parameter_set gaussian_parameter_set_nu1_south_;

	Real const DELTA_CUTOFF_;

	Real const ideal_delta_north_;
	Real const ideal_nu2_north_;
	Real const ideal_nu1_north_;

	Real const ideal_delta_south_;
	Real const ideal_nu2_south_;
	Real const ideal_nu1_south_;

};


}
}
}

#endif
