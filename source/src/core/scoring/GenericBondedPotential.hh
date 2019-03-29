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
#include <core/chemical/Bond.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/EnergyMap.hh>
#include <basic/datacache/CacheableData.hh>

#include <utility/vector1.hh>

namespace core {
namespace scoring {

/// @brief helper function to convert string specification of bondorders to indices
utility::vector1< core::Size >
bondorders_map(std::string bt);

/// @brief convert a bond type to a bin index
core::Size
bin_from_bond(core::chemical::Bond const &b);

/// @brief Parameter set for one torsion angle
/// @details Stores a set of constants required for enumerating a harmonic
///    function given length/angle as variable.
///    This set of parameters are required for bond-length/-angle/improper torsion
///    (gen_bonded_bond/_angle/_improper) scoring. More descriptions below:
///
///    E = k * (x-delta_)**2
///
///    E: Energy value
///    x: variable of the function (e.g. bond length, angle, improper torsion angle)
///    mult_: the specificity/generality of a constraint; when there is ambiguity
///           in multiple parameter set, one with lower mult_ is selected
///
class SpringParams {
public:
	/// @brief Default Constructor
	SpringParams() : k_(0), delta_(0), mult_(1) {}

	/// @brief Constructor with input params
	SpringParams( core::Real k_in, core::Real delta_in, core::Size mult_in ) : k_(k_in), delta_(delta_in), mult_(mult_in) {}

	/// getters
	Real energy ( core::Real value ) const;
	Real deriv ( core::Real value ) const;
	bool is_null () const { return (k_==0); }

	Real k () const { return k_; }
	Real delta () const { return delta_; }

	Size multiplicity() const { return mult_; }

private:
	Real k_, delta_;
	Size mult_;
};

/// @brief Parameter set for one torsion angle
/// @details Stores a set of constants required for enumerating a Karplus
///    cosine function given torsion angle as variable.
///    This set of parameters are required for regular torsion term
///    (gen_bonded_torsion) scoring. More descriptions below:
///
///    E = sum_over_n { k_n * cos(n*x - f_n) } + offset
///
///    E: Energy value
///    n: order of Karplus equation, from 1 to 4
///    x: variable of the function (torsion angle in this case)
///    k1_ ~ k4_: coefficients k_n for each order of n
///    f1_ ~ f4_: phases f_n for each order of n
///    k6_: sixth order k; used only in special case (not included in default params)
///    mult_: the specificity/generality of a constraint; when there is ambiguity
///           in multiple parameter set, one with lower mult_ is selected
///    offset: constant offset to E
///    torsion_type: string tag (just for labeling)
///
class GenTorsionParams {
public:
	/// @brief Default Constructor
	GenTorsionParams( ) : k1_(0), k2_(0), k3_(0), k4_(0), f1_(0), f2_(0), f3_(0), f4_(0), k6_(0),  mult_(1), offset_(0), torsion_type_("") {
		calculate_offset();
	}

	/// @brief Constructor with input params
	GenTorsionParams(
		core::Real k1_in, core::Real k2_in, core::Real k3_in, core::Real k4_in,
		core::Real f1_in, core::Real f2_in, core::Real f3_in, core::Real f4_in,
		core::Size mult_in
	) : k1_(k1_in), k2_(k2_in), k3_(k3_in), k4_(k4_in),
		f1_(f1_in), f2_(f2_in), f3_(f3_in), f4_(f4_in), mult_(mult_in), offset_(0), torsion_type_("") {
		calculate_offset();
	}

	/// @brief Constructor with input params and torsion tag
	GenTorsionParams(
		core::Real k1_in, core::Real k2_in, core::Real k3_in, core::Real k4_in,
		core::Real f1_in, core::Real f2_in, core::Real f3_in, core::Real f4_in,
		core::Size mult_in, std::string torsion_type_in
	) : k1_(k1_in), k2_(k2_in), k3_(k3_in), k4_(k4_in),
		f1_(f1_in), f2_(f2_in), f3_(f3_in), f4_(f4_in), mult_(mult_in), offset_(0), torsion_type_(torsion_type_in) {
		calculate_offset();
	}

	/// getters
	Real energy ( core::Real value ) const;
	Real deriv ( core::Real value ) const;
	bool is_null () const { return (k1_==0 && k2_==0 && k3_==0 && k4_==0 && k6_==0); }

	Size multiplicity() const { return mult_; }
	Real get_params(std::string keyword) const;

	Real k(core::Size idx) const {
		if ( idx==1 ) return k1_;
		if ( idx==2 ) return k2_;
		if ( idx==3 ) return k3_;
		if ( idx==4 ) return k4_;
		if ( idx==6 ) return k6_;
		return 0.0;
	}

	Real offset() const { return offset_; }

	void k6( core::Real value ) { k6_ = value; }

	std::string torsion_type() const {return torsion_type_;}

private:
	/// @brief calculate offset to given torsion param set that makes minimum 0
	void calculate_offset();

private:
	Real k1_,k2_,k3_,k4_; // magnitude
	Real f1_,f2_,f3_,f4_; // phase
	Real k6_; // special case
	Size mult_;
	Real offset_; // constant to make minimum value zero
	std::string torsion_type_;
};

/// @brief Potential for core/scoring/methods/GenericBondedEnergy method.
/// @details Main class calculating gen_bonded energy term. Consists of energy terms
///   of bond-length, -angle, torsion, and improper torsions. Reads in
///   database/scoring/score_function/generic_bonded/generic_bonded.XX.txt and
///   stores it as SpringParams or GenTorsionParams as local variables for
///   score enumeration.
class GenericBondedPotential : public utility::pointer::ReferenceCount  {
public:
	/// @brief Default Constructor
	GenericBondedPotential();

	/// @brief interface with same named function in GenericBondedEnergy
	void
	setup_for_scoring( pose::Pose & pose, scoring::ScoreFunction const &sfxn ) const;

	/// @brief interface with same named function in GenericBondedEnergy
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &pose,
		EnergyMap & emap
	) const;

	/// @brief interface with same named function in GenericBondedEnergy
	void
	residue_derivatives(
		conformation::Residue const & rsd,
		pose::Pose const &pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;

	/// @brief interface with same named function in GenericBondedEnergy
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const &pose,
		EnergyMap & emap
	) const;

	/// @brief interface with same named function in GenericBondedEnergy
	void
	residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const &pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs_r1,
		utility::vector1< DerivVectorPair > & atom_derivs_r2
	) const;

private:

	/// @brief modify torsional parameters "after the fact"
	void
	modify_tors_params(
		std::string const atm1, std::string const atm2,
		std::string bondtype,
		std::string const atm3, std::string const atm4,
		core::Real k1_in, core::Real k2_in, core::Real k3_in, core::Real k4_in,
		core::Real f1_in=0.0, core::Real f2_in=0.0, core::Real f3_in=0.0, core::Real f4_in=0.0
	);

	/// @brief read database file located at database/scoring/score_function/generic_potential
	void
	read_database( std::string filename );

	/// @brief read flags or command line and apply modify_tors_params funciton
	void
	modify_torsion_params_from_cmd_line();

	/// @brief get a uint64_t hash from 2~4 atomtype indices & bond type
	inline uint64_t
	get_parameter_hash( Size bondtype, Size type1, Size type2, Size type3=0, Size type4=0 ) const;

	/// @brief fast database lookup; bond/angles (currently) keyed on type only
	SpringParams const &lookup_bond_params( Size type1, Size type2 ) const;
	SpringParams const &lookup_angle_params( Size type1, Size type2, Size type3 ) const;

	/// @brief fast database lookup; torsions conditioned on types + 2-3 bond
	GenTorsionParams const &lookup_tors_params(
		core::chemical::Bond const &bt, Size type1, Size type2, Size type3, Size type4 ) const;

	/// @brief fast database lookup; improper keyed on types only
	SpringParams const &lookup_improper_params(
		Size type1, Size type2, Size type3, Size type4, int &idx2, int &idx3, int &idx4
	) const;

	/// @brief evaluate each component (1-body) for bond-length term
	Real
	eval_residue_energy_and_derivative_bond(
		conformation::Residue const & rsd,
		utility::vector1< DerivVectorPair > & atom_derivs,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	/// @brief evaluate each component (1-body) for bond-angle term
	Real
	eval_residue_energy_and_derivative_angle(
		conformation::Residue const & rsd,
		utility::vector1< DerivVectorPair > & atom_derivs,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	/// @brief evaluate each component (1-body) for regular torsion term
	Real
	eval_residue_energy_and_derivative_torsion(
		conformation::Residue const & rsd,
		utility::vector1< DerivVectorPair > & atom_derivs,
		//ResidueExclParamsCOP excl_params,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	/// @brief evaluate each component (1-body) for improper-torsion term
	Real
	eval_residue_energy_and_derivative_improper(
		conformation::Residue const & rsd,
		utility::vector1< DerivVectorPair > & atom_derivs,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	/// @brief evaluate each component (2-body) for bond-length term
	Real
	eval_residue_pair_energy_and_derivative_bond(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		utility::vector1< DerivVectorPair > & atom_derivs_r1,
		utility::vector1< DerivVectorPair > & atom_derivs_r2,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	/// @brief evaluate each component (2-body) for bond-angle term
	Real
	eval_residue_pair_energy_and_derivative_angle(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		utility::vector1< DerivVectorPair > & atom_derivs_r1,
		utility::vector1< DerivVectorPair > & atom_derivs_r2,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	/// @brief evaluate each component (2-body) for regular torsion term
	Real
	eval_residue_pair_energy_and_derivative_torsion(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		utility::vector1< DerivVectorPair > & atom_derivs_r1,
		utility::vector1< DerivVectorPair > & atom_derivs_r2,
		//utility::vector1< ResidueExclParamsCOP > const &excl_params,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	/// @brief evaluate each component (2-body) for improper-torsion term
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

	// fast lookup: maps a hash to an index in the vectors above
	boost::unordered_map< uint64_t, Size > bond_lookup_;
	boost::unordered_map< uint64_t, Size > angle_lookup_;
	boost::unordered_map< uint64_t, Size > improper_lookup_;
	boost::unordered_map< uint64_t, Size > tors_lookup_;

	static const SpringParams null_sp_param;
	static const GenTorsionParams null_tors_param;

};

/// @brief compress 5 values into one unsigned int; use 12 bits for each
uint64_t
GenericBondedPotential::get_parameter_hash( Size bondtypr, Size type1, Size type2, Size type3, Size type4 ) const {
	uint64_t retval=type4;
	retval = (retval<<12) + type3;
	retval = (retval<<12) + type2;
	retval = (retval<<12) + type1;
	retval = (retval<<12) + bondtypr;
	return retval;
}

} // scoring
} // core


#endif // INCLUDED_core_scoring_GenericBondedPotential_hh

