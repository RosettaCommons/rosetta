// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/rna/RNA_FittedTorsionInfo.cc
/// @brief  Statistically derived torsion information for RNA
/// @author Rhiju Das


// Unit headers
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>

///////////////////////////////////////////////////////
// Keep track of information from, e.g., chemical
// accessibility experiments -- useful for chemical.
///////////////////////////////////////////////////////

namespace core {
namespace chemical {
namespace rna {

RNA_FittedTorsionInfo::RNA_FittedTorsionInfo():
	utility::pointer::ReferenceCount(),
	rna_tight_torsions_( true ),
	//// Following are torsions that are optimized to minimize rna_torsion + sugar_close.
	////  [ Due to frustration, impossible to satisfy ideal bond lengths+angles for arbitrary delta. ]
	ideal_delta_north_( 85.4184 ),
	ideal_nu2_north_( 37.3643 ),
	ideal_nu1_north_( 92.2042 ),
	ideal_delta_south_( 152.467 ),
	ideal_nu2_south_( -38.085 ),
	ideal_nu1_south_( 150.133 ),
	delta_cutoff_(115.0)
{
	init_rna_torsion_gaussian_parameters();
}

RNA_FittedTorsionInfo::~RNA_FittedTorsionInfo(){}

///////////////////////////////////////////////////////////////////////////////////////
void
RNA_FittedTorsionInfo::init_rna_torsion_gaussian_parameters()
{
	// These numbers refer to amplitude, mean, and width of fitted Gaussians.

	////////////////////////////////////////////////////////////////////////////
	// Parameters for alpha,beta,gamma,delta, etc. torsion constraints...
	////////////////////////////////////////////////////////////////////////////
	// "Tight torsion" means use the width of the major peak in the torsion
	// histogram to determine how strong the harmonic constraint will be. Otherwise, be less restrictive,
	// and use width of the minor fatter peak that lies under the sharp peak.

	if ( rna_tight_torsions_ ) {
		gaussian_parameter_set_alpha_.emplace_back( 222.03, -64.11,  9.64 );
	} else {
		gaussian_parameter_set_alpha_.emplace_back( 13.12, -73.74,  35.97 );
	}
	gaussian_parameter_set_alpha_.emplace_back( 10.51, 66.01,  18.09 );
	gaussian_parameter_set_alpha_.emplace_back( 18.40, 161.80,  18.12 );

	if ( rna_tight_torsions_ ) {
		gaussian_parameter_set_beta_.emplace_back( 181.33, 176.33,  11.54 );
	} else {
		gaussian_parameter_set_beta_.emplace_back( 32.30, 174.52,  43.56 );
	}

	if ( rna_tight_torsions_ ) {
		gaussian_parameter_set_gamma_.emplace_back( 366.90, 53.08,  6.64 );
	} else {
		gaussian_parameter_set_gamma_.emplace_back( 18.26, 56.59,  20.57 );
	}
	gaussian_parameter_set_gamma_.emplace_back( 21.61, 178.19,  13.61 );
	gaussian_parameter_set_gamma_.emplace_back( 3.98, -64.02,  17.76 );

	gaussian_parameter_set_delta_north_.emplace_back( 687.92, 82.90,  3.99 );
	gaussian_parameter_set_delta_south_.emplace_back( 53.18, 145.25,  6.35 );

	gaussian_parameter_set_epsilon_north_.emplace_back( 178.08, -150.17,  14.64 );
	gaussian_parameter_set_epsilon_north_.emplace_back( 2.52, 68.28,  32.29 );

	gaussian_parameter_set_epsilon_south_.emplace_back( 11.95, -98.45, 26.80 );
	gaussian_parameter_set_epsilon_south_.emplace_back( 0.58, 159.70, 103.86 );

	if ( rna_tight_torsions_ ) {
		gaussian_parameter_set_zeta_alpha_sc_minus_.emplace_back( 143.97, -71.45, 7.91 );
	} else {
		gaussian_parameter_set_zeta_alpha_sc_minus_.emplace_back( 78.74, -68.60, 16.19 );
	}
	// gaussian_parameter_set_zeta_alpha_sc_minus_.emplace_back(  2.43, 178.84, 114.82 );
	gaussian_parameter_set_zeta_alpha_sc_minus_.emplace_back(  2.43, 178.84, 40.00 );

	gaussian_parameter_set_zeta_alpha_sc_plus_.emplace_back( 2.08, -137.28,  63.12 );
	gaussian_parameter_set_zeta_alpha_sc_plus_.emplace_back( 3.37, 87.07,  32.69 );

	gaussian_parameter_set_zeta_alpha_ap_.emplace_back( 13.65, -69.74,  15.28 );
	gaussian_parameter_set_zeta_alpha_ap_.emplace_back( 2.69, 63.03,  33.61 );

	//These need to be flipped by -180.0, due to change in chi def. from r++ to mini.
	// gaussian_parameter_set_chi_north_.emplace_back( 228.92, -100.57, 11.43 );
	// gaussian_parameter_set_chi_north_.emplace_back( 1.07, 129.24, 26.11 );
	gaussian_parameter_set_chi_north_.emplace_back( 228.92, 79.43, 11.43 );
	gaussian_parameter_set_chi_north_.emplace_back( 1.07, -50.76, 26.11 );

	// gaussian_parameter_set_chi_south_.emplace_back( 12.53, -63.40, 25.23 );
	// gaussian_parameter_set_chi_south_.emplace_back( 0.64, 130.09, 28.10 );
	gaussian_parameter_set_chi_south_.emplace_back( 12.53, 116.60, 25.23 );
	gaussian_parameter_set_chi_south_.emplace_back( 0.64, -49.91, 28.10 );

	gaussian_parameter_set_nu2_north_.emplace_back( 1148.01, 38.82, -2.77 );
	gaussian_parameter_set_nu2_south_.emplace_back( 173.15, -37.22, 3.00 );

	gaussian_parameter_set_nu1_north_.emplace_back( 631.56, 95.34,  4.20 );
	gaussian_parameter_set_nu1_south_.emplace_back( 57.04, 155.51,  6.00 );

}

}
}
}
