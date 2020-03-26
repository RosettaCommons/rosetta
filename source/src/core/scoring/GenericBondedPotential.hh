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
#include <core/chemical/Bond.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/scoring/DerivVectorPair.hh>
#include <core/scoring/EnergyMap.hh>
#include <basic/datacache/CacheableData.hh>

#include <utility/vector1.hh>

//#include <boost/unordered_map.hpp>
#include <unordered_map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

/// @brief compress 5 values into one unsigned int; use 12 bits for each
uint64_t
get_parameter_hash( Size bondtypr, Size type1, Size type2, Size type3=0, Size type4=0 );

/// @brief helper function to convert string specification of bondorders to indices
utility::vector1< core::Size >
bondorders_map(std::string bt);

/// @brief convert a bond type to a bin index
core::Size
bin_from_bond(core::chemical::BondName bn, core::chemical::BondRingness br);

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

	void set_offset( core::Real value ) { offset_ = value; }

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
class GenericBondedPotential : public utility::VirtualBase  {
public:
	/// @brief Default Constructor
	GenericBondedPotential();

	/// @brief Constructor with exclusion behaviors
	//GenericBondedPotential( bool score_full, bool score_hybrid );

	/// @brief interface with same named function in GenericBondedEnergy
	void
	setup_for_scoring( pose::Pose & pose,
		scoring::ScoreFunction const &sfxn,
		bool const &score_full=false,
		bool const &score_hybrid=true
	) const;

	/// @brief interface with same named function in GenericBondedEnergy
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const &pose,
		EnergyMap & emap,
		bool const &score_full,
		bool const &score_hybrid
	) const;

	/// @brief interface with same named function in GenericBondedEnergy
	void
	residue_derivatives(
		conformation::Residue const & rsd,
		pose::Pose const &pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs,
		bool const &score_full,
		bool const &score_hybrid
	) const;

	/// @brief interface with same named function in GenericBondedEnergy
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const &pose,
		EnergyMap & emap,
		bool const &score_full,
		bool const &score_hybrid
	) const;

	/// @brief interface with same named function in GenericBondedEnergy
	void
	residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const &pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs_r1,
		utility::vector1< DerivVectorPair > & atom_derivs_r2,
		bool const &score_full,
		bool const &score_hybrid
	) const;

	/*
	void
	score_full( bool setting ){ score_full_ = setting; } // score all angle/tors

	void
	score_hybrid( bool setting ){ score_hybrid_ = setting; } // score only those are missing in consideration by other terms
	*/

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

	/*
	/// @brief get a uint64_t hash from 2~4 atomtype indices & bond type
	inline uint64_t
	get_parameter_hash( Size bondtype, Size type1, Size type2, Size type3=0, Size type4=0 ) const;
	*/

	/// @brief fast database lookup; bond/angles (currently) keyed on type only
	SpringParams const &lookup_bond_params( Size type1, Size type2 ) const;
	SpringParams const &lookup_angle_params( Size type1, Size type2, Size type3 ) const;

	/// @brief fast database lookup; torsions conditioned on types + 2-3 bond
	GenTorsionParams const & lookup_tors_params(
		core::chemical::BondName bn, core::chemical::BondRingness br,
		Size type1, Size type2, Size type3, Size type4 ) const;

	/// @brief fast database lookup; improper keyed on types only
	SpringParams const & lookup_improper_params(
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
		ResidueExclParamsCOP excl_params,
		Real const weight = 0.0,
		bool const calc_deriv = false
	) const;

	/// @brief evaluate each component (1-body) for improper-torsion term
	Real
	eval_residue_energy_and_derivative_improper(
		conformation::Residue const & rsd,
		utility::vector1< DerivVectorPair > & atom_derivs,
		ResidueExclParamsCOP excl_params,
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
		//utility::vector1< ResidueExclParamsCOP > excl_params,
		ResidueExclParamsCOP excl_info_res1,
		ResidueExclParamsCOP excl_info_res2,
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
		ResidueExclParamsCOP excl_info_res1,
		ResidueExclParamsCOP excl_info_res2,
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
	std::unordered_map< uint64_t, Size > bond_lookup_;
	std::unordered_map< uint64_t, Size > angle_lookup_;
	std::unordered_map< uint64_t, Size > improper_lookup_;
	std::unordered_map< uint64_t, Size > tors_lookup_;

	static const SpringParams null_sp_param;
	static const GenTorsionParams null_tors_param;

};

class ResidueExclParams : public utility::VirtualBase {

public:
	ResidueExclParams( core::chemical::ResidueTypeCOP rsd_type,
		bool const score_full,
		bool const score_hybrid );

	void create_excl_info( bool const score_full,
		bool const score_hybrid );

	void fully_counted( bool const value ) { fully_counted_ = value; }
	bool fully_counted() const { return fully_counted_; }

	void fully_excluded_1b( bool const value ) { fully_excluded_1b_ = value; }
	bool fully_excluded_1b() const { return fully_excluded_1b_; }

	// atm index
	void atm_excluded( Size i, bool value ) { atm_excluded_[i] = value; }

	// atm index
	bool atm_excluded( Size i ) const { return atm_excluded_[i]; }

	// either dihedral or bond index
	void store_by_bondno( Size i ) { bondnos_[i] = true; }

	// either dihedral or bond index
	bool find_by_bondno( Size i ) const { return bondnos_[i]; }

	// central 2 atms
	void store_by_atmpairno( Size i, Size j );

	// central 2 atms
	bool find_by_atmpairno( Size i, Size j ) const;

	void rsd1type( core::chemical::ResidueType const &rsdtype ){ rsdtype1_ = rsdtype.name(); }
	void rsd2type( core::chemical::ResidueType const &rsdtype ){ rsdtype2_ = rsdtype.name(); }

	//bool empty() const { return ; }
	/*
	core::chemical::ResidueType rsd1type() const { return rsdtype1_; }

	core::chemical::ResidueType rsd2type() const { return rsdtype2_; }
	*/
	bool same_rsdpairs( core::conformation::Residue const &rsd1,
		core::conformation::Residue const &rsd2 ) const
	{ return (rsd1.type().name() == rsdtype1_) &&
		(rsd2.type().name() == rsdtype2_); }

private:
	bool fully_counted_;
	bool fully_excluded_1b_;
	std::string rsdtype1_, rsdtype2_;
	utility::vector1< bool > atm_excluded_;
	utility::vector1< bool > bondnos_;
	std::unordered_map< uint64_t, bool > atmpairnos_;
	core::chemical::ResidueTypeCOP rsd_type_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ResidueExclParams();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class GenBondedExclInfo : public basic::datacache::CacheableData {

public:
	// constructor
	GenBondedExclInfo( bool const score_full,
		bool const score_hybrid );

	// virtual
	basic::datacache::CacheableDataOP clone() const override
	{ return basic::datacache::CacheableDataOP( new GenBondedExclInfo( *this ) );}

	//void update_poseinfo( pose::Pose const &pose );

	void
	add_residue_exclude_torsions( chemical::ResidueType const &rsd_type
	);

	/* // unused
	void
	add_residue_pair_exclude_torsions( core::conformation::Residue const &rsd1,
	core::conformation::Residue const &rsd2
	);
	*/

	ResidueExclParamsCOP
	get_residue_data( core::chemical::ResidueType const &rsd_type ) const;

	ResidueExclParamsCOP
	get_residue_pair_data( core::Size seqpos1,
		core::Size seqpos2 ) const;


private:
	bool score_full_;
	bool score_hybrid_;
	std::unordered_map< std::string, ResidueExclParamsOP > perres_excls_;
	std::unordered_map< uint64_t, ResidueExclParamsOP > respair_excls_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	GenBondedExclInfo();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_ResidueExclParams )
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_GenBondedExclInfo )
#endif // SERIALIZATION

#endif // INCLUDED_core_scoring_GenericBondedPotential_hh

