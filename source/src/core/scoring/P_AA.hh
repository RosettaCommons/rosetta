// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/P_AA.hh
/// @brief  Amino acid probability arrays and functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Andrew Leaver-Fay -- porting Stuart's code

#ifndef INCLUDED_core_scoring_P_AA_hh
#define INCLUDED_core_scoring_P_AA_hh

// Unit Headers
#include <core/scoring/P_AA.fwd.hh>

// Package headers
#include <core/scoring/types.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/TorsionID.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/interpolation/spline/Bicubic_spline.hh>

namespace core {
namespace scoring {

class P_AA : public utility::pointer::ReferenceCount
{

private:

	/// @brief Amino acid probability array: P(aa)
	//typedef  utility::keys::SmallKeyVector< conformation::amino::AminoAcidKey, Probability >  Probability_AA;
	//extern Probability_AA P_AA;
	typedef utility::vector1< Probability > Probability_AA;
	Probability_AA P_AA_;

	/// @brief Amino acid conditional probability wrt number of neighbors array: P(aa|neighbors)
	//typedef  utility::keys::SmallKeyVector< conformation::amino::AminoAcidKey, FArray1D_Probability >  Probability_AA_n;
	//extern Probability_AA_n P_AA_n;
	typedef utility::vector1< utility::vector1< Probability > > Probability_AA_n;
	Probability_AA_n P_AA_n_;

	/// @brief Amino acid conditional probability wrt (phi,psi) array: P(aa|phi,psi)
	//typedef  utility::keys::SmallKeyVector< conformation::amino::AminoAcidKey, FArray2D_Probability >  Probability_AA_pp;
	//extern Probability_AA_pp P_AA_pp;
	typedef utility::vector1< FArray2D_Probability > Probability_AA_pp;
	Probability_AA_pp P_AA_pp_;

	utility::vector1< numeric::interpolation::spline::BicubicSpline > P_AA_pp_energy_splines_;

public:
	P_AA();

	~P_AA() override;
private:

	/// @brief Gets whether this amino acid enum type is a canonical D-amino acid.
	bool
	is_canonical_d_aminoacid(
		core::chemical::AA const res_aa
	) const;

	core::chemical::AA
	get_l_equivalent(
		core::chemical::AA const d_aa
	) const;

	/// @brief Initialize the amino acid probability data structures
	void
	P_AA_initialize();


	/// @brief Read the amino acid probability file into P_AA
	void
	read_P_AA();


	/// @brief Read the amino acid conditional probability wrt (neighbors) file into P_AA_n
	void
	read_P_AA_n();


	/// @brief Read the amino acid conditional probability wrt (phi,psi) file into P_AA_pp
	void
	read_P_AA_pp();

	/// @brief Symmetrize the glyceine P_AA_pp table, if the user has used the -symmetric_gly_tables option.
	/// @details The gly table must already be loaded before this is called.  Also, this should be called before
	/// the bicubic splines are set up for the energy table.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void symmetrize_gly_table();

public:
	/// @brief Probability energies from P(aa|phi,psi)
	//Energy
	//P_AA_pp_energy( conformation::amino::AminoAcid const & aa );
	Energy
	P_AA_pp_energy( conformation::Residue const & ) const;

	/// @brief Low-level probability energies from P(aa|phi,psi)
	/// @brief Probability energies from P(aa|phi,psi): Low level calculation for non-terminus position
	//Energy
	//P_AA_pp_energy( conformation::amino::AminoAcidKey const & key, Angle const phi, Angle const psi );
	Energy
	P_AA_pp_energy( chemical::AA aa, Angle const phi, Angle const psi ) const;

	EnergyDerivative
	get_Paa_pp_deriv(
		conformation::Residue const & res,
		id::TorsionID const & tor_id
	) const;

	/// @brief Probability energies for P(aa)
	Energy
	P_AA_energy( conformation::Residue const & ) const;


};

} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_P_AA_HH
