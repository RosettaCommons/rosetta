// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/carbon_hbonds/CarbonHBondPotential.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_carbon_hbonds_CarbonHBondPotential_hh
#define INCLUDED_core_scoring_carbon_hbonds_CarbonHBondPotential_hh

// Unit Headers
#include <core/scoring/carbon_hbonds/CarbonHBondPotential.fwd.hh>

// Package headers

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <ObjexxFCL/FArray1D.hh>


#include <utility/vector1_bool.hh>


namespace core {
namespace scoring {
namespace carbon_hbonds {

class CarbonHBondPotential : public utility::pointer::ReferenceCount {

public:

	/// @brief ctor, reads data file
	CarbonHBondPotential();

	/// @brief Calculate chbond energies for non-rna atom pairs.
	Real get_potential(
		Size const & atom_type,
		Vector const & H_A_vector, //vector of hydrogen to acceptor
		Vector const & D_H_vector, //vector of donor hv atom to hydrogen
		Vector const & B_A_vector, //vector of acceptor's base atom to acceptor
		bool calculate_deriv, // early exit if this is false
		Vector & deriv_vector
	) const;

	/// @brief Calculate the rna-specific chbond energies.  The derivative vector returned is the
	/// force vector on the acceptor atom.  Multiply by -1 to get the force vector on the donor atom
	Real
	get_potential_RNA(
		Vector const & r_H_A,
		Vector const & z_i /*unit vector pointing from donor to its connected hydrogen*/,
		bool const & update_deriv,
		Vector & deriv
	) const;

	/// @brief second declaration to allow skipping deriv; gcc 4.1.3 does not like setting default parameters for a pass-by-reference parameter
	Real get_potential(
		Size const & atom_type,
		Vector const & H_A_vector, //vector of hydrogen to acceptor
		Vector const & D_H_vector, //vector of donor hv atom to hydrogen
		Vector const & B_A_vector  //vector of acceptor's base atom to acceptor
	) const;

	Real max_dis() const {
		return max_dist_;
	}

private: //data

	void
	read_potential();

	//Which atoms to loop over during VDW check?
	Real bin_width_;
	Real max_dist_;
	Real const aroC_scale_factor_;

	Size const num_carbon_donor_atoms_;
	Size const num_bins_;
	utility::vector1<  ObjexxFCL::FArray1D <Real>  > carbon_hbond_parameter_;
	// ObjexxFCL::FArray2D < Real > carbon_hbond_deriv_ /*(5, 45)*/;

	utility::vector1< std::string > carbon_donors_;
	ObjexxFCL::FArray1D < Size > standard_atomtype_to_carbon_donor_index_;

	Real const rna_cos_theta_cutoff_;
	Real const rna_cos_theta_fade_zone_;
	Real const rna_ch_o_bond_distance_;

};

}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
