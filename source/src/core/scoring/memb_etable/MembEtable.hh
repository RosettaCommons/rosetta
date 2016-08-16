// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/scoring/memb_etable/MembEtable.hh
/// @brief Table of pre-computed LK membrane solvation energies
/// @author Patrick Barth (original)
/// @author Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_scoring_memb_etable_MembEtable_hh
#define INCLUDED_core_scoring_memb_etable_MembEtable_hh

// Unit Headers
#include <core/scoring/memb_etable/MembEtable.fwd.hh>

// Package Headers
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/etable/EtableOptions.hh>

// Project Headers
#include <core/chemical/AtomTypeSet.hh>
#include <utility/pointer/access_ptr.hh>

// Utility Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray2.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

namespace core {
namespace scoring {
namespace etable {

/// @brief Table of pre-computed LK membrane solvation energies
class MembEtable : public Etable {

public:

	// Constructors

	MembEtable(
		chemical::AtomTypeSetCAP atom_set_in,
		EtableOptions const & options,
		std::string const & alternate_parameter_set = ""
	);

	MembEtable( MembEtable const & src );

	virtual ~MembEtable();

	// Accessors

	/// @brief Provide Constnat Access to Arrays
	ObjexxFCL::FArray3D< Real > const &
	solv1() const
	{
		return solv1_;
	}

	ObjexxFCL::FArray3D< Real > const &
	solv2() const
	{
		return solv2_;
	}

	ObjexxFCL::FArray3D< Real > const &
	memb_solv1() const
	{
		return memb_solv1_;
	}

	ObjexxFCL::FArray3D< Real > const &
	memb_solv2() const
	{
		return memb_solv2_;
	}

	/// @brief Return the solvation derivative table for the desolvation of atom1 by atom2
	ObjexxFCL::FArray3D< Real > const &
	dsolv1() const
	{
		return dsolv1_;
	}

	/// @brief Return the solvation derivative table that combines atom1 and atom2's desolvations
	ObjexxFCL::FArray3D< Real > const &
	dsolv2() const
	{
		return dsolv2_;
	}

	/// @brief Return the solvation derivative table for the desolvation of atom1 by atom2
	ObjexxFCL::FArray3D< Real > const &
	memb_dsolv1() const
	{
		return memb_dsolv1_;
	}

	/// @brief return the solvation derivative table that combines atom1 and atom2's desolvations
	ObjexxFCL::FArray3D< Real > const &
	memb_dsolv2() const
	{
		return memb_dsolv2_;
	}

	/// @brief Returns the maximum lj radius for any non-hydrogen
	/// atom as defined by the atom-type-set used to create this Etable.
	virtual
	Real
	max_non_hydrogen_lj_radius() const {
		return max_non_hydrogen_lj_radius_;
	}

	/// @brief Returns the maximum lj radius for any hydrogen atom as
	/// defined by the input atom-type-set used to create this Etable.
	virtual
	Real
	max_hydrogen_lj_radius() const {
		return max_hydrogen_lj_radius_;
	}

	/// set these up in the ctor
	virtual
	Real
	lj_radius( int const i ) const
	{
		return lj_radius_[i];
	}

	virtual
	Real
	lk_dgfree( int const i ) const
	{
		return lk_dgfree_[i];
	}

	virtual
	Real
	lk_volume( int const i ) const
	{
		return lk_volume_[i];
	}

	virtual
	Real
	lk_lambda( int const i ) const
	{
		return lk_lambda_[i];
	}

	Real
	memb_lk_dgfree( int const i ) const
	{
		return memb_lk_dgfree_[i];
	}

	ObjexxFCL::FArray1D< Real > const &
	lk_dgrefce() const
	{
		return lk_dgrefce_;
	}

	ObjexxFCL::FArray1D< Real > const &
	memb_lk_dgrefce() const
	{
		return memb_lk_dgrefce_;
	}

public: // Interfaces for convenient IO

	virtual
	void
	output_etable(
		ObjexxFCL::FArray3D<Real> & etable,
		std::string label,
		std::ostream & out
	);

	virtual
	void
	input_etable(
		ObjexxFCL::FArray3D<Real> & etable,
		std::string const & label,
		std::istream & in
	);


private: // methods

	// Initialization functions
	void dimension_memb_etable_arrays();
	void initialize_from_atomset( chemical::AtomTypeSetCAP atom_set_in );

	// Etable Calculations
	void smooth_etables();
	void modify_pot();
	void make_pairenergy_table();

	void
	precalc_etable_coefficients(
		ObjexxFCL::FArray2< Real > & lj_sigma,
		ObjexxFCL::FArray1< Real > & lk_inv_lambda2,
		ObjexxFCL::FArray2< Real > & lk_coeff,
		ObjexxFCL::FArray2< Real > & memb_lk_coeff,
		ObjexxFCL::FArray2< Real > & lk_min_dis2sigma_value,
		ObjexxFCL::FArray2< Real > & memb_lk_min_dis2sigma_value
	);

	void
	calc_etable_value(
		Real & dis2,
		int & atype1,
		int & atype2,
		Real & solvE1,
		Real & solvE2,
		Real & dsolvE1,
		Real & dsolvE2,
		ObjexxFCL::FArray2< Real > & lj_sigma,
		ObjexxFCL::FArray1< Real > & lk_inv_lambda2,
		ObjexxFCL::FArray2< Real > & lk_coeff,
		ObjexxFCL::FArray2< Real > & lk_min_dis2sigma_value,
		Real & memb_solvE1,
		Real & memb_solvE2,
		ObjexxFCL::FArray2< Real > & memb_lk_coeff,
		ObjexxFCL::FArray2< Real > & memb_lk_min_dis2sigma_value,
		Real & memb_dsolvE1,
		Real & memb_dsolvE2
	);

private: // data

	Real max_non_hydrogen_lj_radius_;
	Real max_hydrogen_lj_radius_;

	utility::vector1< Real > lj_radius_;
	utility::vector1< Real > lk_dgfree_;
	utility::vector1< Real > lk_volume_;
	utility::vector1< Real > lk_lambda_;
	utility::vector1< Real > memb_lk_dgfree_;
	ObjexxFCL::FArray1D< Real > lk_dgrefce_;
	ObjexxFCL::FArray1D< Real > memb_lk_dgrefce_;

	ObjexxFCL::FArray3D< Real > solv1_;
	ObjexxFCL::FArray3D< Real > solv2_;
	ObjexxFCL::FArray3D< Real > dsolv1_;
	ObjexxFCL::FArray3D< Real > dsolv2_;
	ObjexxFCL::FArray3D< Real > memb_solv1_;
	ObjexxFCL::FArray3D< Real > memb_solv2_;
	ObjexxFCL::FArray3D< Real > memb_dsolv1_;
	ObjexxFCL::FArray3D< Real > memb_dsolv2_;

};

} // etable
} // scoring
} // core

#endif // INCLUDED_core_scoring_memb_etable_MembEtable_hh
