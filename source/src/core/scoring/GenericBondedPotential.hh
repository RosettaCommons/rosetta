// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/GenericBondedPotential.hh
/// @brief  A "generic" (atom-type-only-based) torsional potential
/// @author Hahnbeom Park and Frank DiMaio


#ifndef INCLUDED_core_scoring_GenericBondedPotential_hh
#define INCLUDED_core_scoring_GenericBondedPotential_hh

// Unit headers
#include <core/scoring/GenericBondedPotential.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/EnergyMap.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {

// Bond-length/Bond-angle/Improper-tors params
//    mult_ describes the specificity/generality of a constraint
class SpringParams {
public:
	SpringParams() : k_(0), delta_(0), mult_(1) {}
	SpringParams( core::Real k_in, core::Real delta_in, core::Size mult_in ) : k_(k_in), delta_(delta_in), mult_(mult_in) {}

	Real energy ( core::Real value ) const;
	Real deriv ( core::Real value ) const;
	bool is_null () const { return (k_==0); }

	Size multiplicity() const { return mult_; }

private:
	Real k_, delta_;
	Size mult_;
};

// Torsional parameter
class GenTorsionParams {
public:
	GenTorsionParams( ) : k1_(0), k2_(0), k3_(0), f1_(0), f2_(0), f3_(0), mult_(1) {}
	GenTorsionParams(
		core::Real k1_in, core::Real k2_in, core::Real k3_in,
		core::Real f1_in, core::Real f2_in, core::Real f3_in, core::Size mult_in
	) : k1_(k1_in), k2_(k2_in), k3_(k3_in), f1_(f1_in), f2_(f2_in), f3_(f3_in), mult_(mult_in) {}

	Real energy ( core::Real value ) const;
	Real deriv ( core::Real value ) const;
	bool is_null () const { return (k1_==0 && k2_==0 && k3_==0); }

	Size multiplicity() const { return mult_; }
	Real get_params(std::string keyword) const;

private:
	Real k1_,k2_,k3_; // magnitude
	Real f1_,f2_,f3_; // phase
	Size mult_;
};

// "Special" torsional parameter (high-order components)
class SpecialGenTorsionParams {
public:
	SpecialGenTorsionParams( ) : k1_(0), k2_(0), k3_(0), k4_(0), k8_(0), f1_(0), f2_(0), f3_(0), f4_(0), f8_(0), mult_(1) {}
	SpecialGenTorsionParams(
		core::Real k1_in, core::Real k2_in, core::Real k3_in, core::Real k4_in, core::Real k8_in,
		core::Real f1_in, core::Real f2_in, core::Real f3_in, core::Real f4_in, core::Real f8_in, core::Size mult_in
	) : k1_(k1_in), k2_(k2_in), k3_(k3_in), k4_(k4_in), k8_(k8_in), f1_(f1_in), f2_(f2_in), f3_(f3_in), f4_(f4_in), f8_(f8_in), mult_(mult_in) {}

	Real energy ( core::Real value ) const;
	Real deriv ( core::Real value ) const;
	bool is_null () const { return (k1_==0 && k2_==0 && k3_==0 && k4_==0 && k8_==0); }

	Size multiplicity() { return mult_; }
	Real get_params(std::string keyword) const;

private:
	Real k1_,k2_,k3_,k4_,k8_; // magnitude
	Real f1_,f2_,f3_,f4_,f8_; // phase
	Size mult_;
};



///////////////////////////////////

class GenericBondedPotential : public utility::pointer::ReferenceCount  {
public:
	GenericBondedPotential();

	void
	setup_for_scoring( pose::Pose & pose ) const;

	void
	residue_energy(
		conformation::Residue const & rsd,
		EnergyMap & emap
	) const;

	void
	residue_derivatives(
		conformation::Residue const & rsd,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		EnergyMap & emap
	) const;

	void
	residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs_r1,
		utility::vector1< DerivVectorPair > & atom_derivs_r2
	) const;

private:
	void
	modify_tors_params(
		std::string const atm1, std::string const atm2,
		std::string const atm3, std::string const atm4,
		core::Real k1_in, core::Real k2_in, core::Real k3_in,
		core::Real f1_in=0.0, core::Real f2_in=0.0,
		core::Real f3_in=0.0
	);

	void
	modify_special_tors_params(
		std::string const atm1, std::string const atm2,
		std::string const atm3, std::string const atm4,
		core::Real k1_in, core::Real k2_in, core::Real k3_in,
		core::Real k4_in, core::Real k8_in,
		core::Real f1_in=0.0, core::Real f2_in=0.0,
		core::Real f3_in=0.0, core::Real f4_in=0.0,
		core::Real f8_in=0.0
	);

private:
	void
	read_database( std::string filename );

	void
	modify_torsion_params_from_cmd_line();

	void
	modify_special_torsion_params_from_cmd_line();

	inline uint64_t
	get_parameter_hash( Size type1, Size type2, Size type3=0, Size type4=0 ) const;

	// fast database lookup
	SpringParams const &lookup_bond_params( Size type1, Size type2 ) const;
	SpringParams const &lookup_angle_params( Size type1, Size type2, Size type3 ) const;
	GenTorsionParams const &lookup_tors_params( Size type1, Size type2, Size type3, Size type4 ) const;
	SpecialGenTorsionParams const &lookup_special_tors_params( Size type1, Size type2, Size type3, Size type4 ) const;
	SpringParams const &lookup_improper_params(
		Size type1, Size type2, Size type3, Size type4, int &idx2, int &idx3, int &idx4
	) const;

	// evaluate each component (1-body)
	Real
	eval_residue_energy_and_derivative_bond(
		conformation::Residue const & rsd,
		utility::vector1< DerivVectorPair > & atom_derivs,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	Real
	eval_residue_energy_and_derivative_angle(
		conformation::Residue const & rsd,
		utility::vector1< DerivVectorPair > & atom_derivs,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	Real
	eval_residue_energy_and_derivative_torsion(
		conformation::Residue const & rsd,
		utility::vector1< DerivVectorPair > & atom_derivs,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	Real
	eval_residue_energy_and_derivative_improper(
		conformation::Residue const & rsd,
		utility::vector1< DerivVectorPair > & atom_derivs,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	// evaluate each component (2-body)
	Real
	eval_residue_pair_energy_and_derivative_bond(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		utility::vector1< DerivVectorPair > & atom_derivs_r1,
		utility::vector1< DerivVectorPair > & atom_derivs_r2,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	Real
	eval_residue_pair_energy_and_derivative_angle(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		utility::vector1< DerivVectorPair > & atom_derivs_r1,
		utility::vector1< DerivVectorPair > & atom_derivs_r2,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	Real
	eval_residue_pair_energy_and_derivative_torsion(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		utility::vector1< DerivVectorPair > & atom_derivs_r1,
		utility::vector1< DerivVectorPair > & atom_derivs_r2,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	Real
	eval_residue_pair_energy_and_derivative_improper(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		utility::vector1< DerivVectorPair > & atom_derivs_r1,
		utility::vector1< DerivVectorPair > & atom_derivs_r2,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

private:
	// all the atom types known to this potential
	utility::vector1< int > defined_atom_types_;
	// make all the atom name to index map known to this potential
	std::map< std::string, utility::vector1< core::Size > > name_index_map;
	core::Size natomtypes;
	// the potentials
	utility::vector1< SpringParams > bond_pot_;
	utility::vector1< SpringParams > angle_pot_;
	utility::vector1< SpringParams > improper_pot_;
	utility::vector1< GenTorsionParams > tors_pot_;
	utility::vector1< SpecialGenTorsionParams > special_tors_pot_;

	// fast lookup: maps a hash to an index in the vectors above
	boost::unordered_map< uint64_t, Size > bond_lookup_;
	boost::unordered_map< uint64_t, Size > angle_lookup_;
	boost::unordered_map< uint64_t, Size > improper_lookup_;
	boost::unordered_map< uint64_t, Size > tors_lookup_;
	boost::unordered_map< uint64_t, Size > special_tors_lookup_;

	static const SpringParams null_sp_param;
	static const GenTorsionParams null_tors_param;
	static const SpecialGenTorsionParams null_sp_tors_param;
};

//
uint64_t
GenericBondedPotential::get_parameter_hash( Size type1, Size type2, Size type3, Size type4 ) const {
	uint64_t retval=type4;
	retval = (retval<<16) + type3;
	retval = (retval<<16) + type2;
	retval = (retval<<16) + type1;
	return retval;
}

} // scoring
} // core


#endif // INCLUDED_core_scoring_GenericBondedPotential_hh

