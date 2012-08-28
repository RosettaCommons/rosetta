// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CartesianBondedEnergy.hh
/// @brief
/// @author

#ifndef INCLUDED_core_scoring_methods_CartesianBondedEnergy_hh
#define INCLUDED_core_scoring_methods_CartesianBondedEnergy_hh

// Unit headers
#include <core/scoring/methods/CartesianBondedEnergy.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/mm/MMBondAngleLibrary.fwd.hh>
#include <core/scoring/mm/MMBondLengthLibrary.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/types.hh>

// boost
#include <boost/unordered_map.hpp>
#include <boost/tuple/tuple.hpp>

// C++ headers
#include <iostream>

#include <utility/vector1.hh>

//#include <map>


typedef boost::tuples::tuple< std::string, int, int, int, int > residx_atm_quad;
typedef boost::tuples::tuple< std::string, int, int, int > residx_atm_triple;
typedef boost::tuples::tuple< std::string, int, int > residx_atm_pair;
//typedef boost::tuples::tuple< std::string, std::string, std::string > atm_name_triple;

namespace boost {
namespace tuples {

std::size_t hash_value(residx_atm_quad const& e); 
std::size_t hash_value(residx_atm_triple const& e); 
std::size_t hash_value(residx_atm_pair const& e); 
//std::size_t hash_value(atm_name_triple const& e); 
bool operator==(residx_atm_quad const& a,residx_atm_quad const& b); 
bool operator==(residx_atm_triple const& a,residx_atm_triple const& b); 
bool operator==(residx_atm_pair const& a,residx_atm_pair const& b); 
//bool operator==(atm_name_triple const& a,atm_name_triple const& b); 

}
}

namespace core {
namespace scoring {
namespace methods {

///////////////////
//fpd  cache ideal torsions for fast lookup
class TorsionDatabase  : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~TorsionDatabase();
	TorsionDatabase(Real k_tors, Real k_tors_prot);

	// lookup ideal torsion; insert in DB if not there
	void													
	lookup( core::chemical::ResidueType const & restype,
	        int atm1, int atm2, int atm3, int atm4, Real &Kphi, Real &phi0, Real &phi_step );

	Real
	ktors() { return k_torsion; }

private:
	//fpd involves a string comparison (restype) .. could be faster
	//std::map< residx_atm_triple, core::Real > bondangles_;
	boost::unordered_map< residx_atm_quad, core::Real > torsions_;
	boost::unordered_map< residx_atm_quad, core::Real > torsion_steps_;
	boost::unordered_map< residx_atm_quad, core::Real > Kphis_;

	Real k_torsion, k_torsion_proton;   // can only be set in constructor
};


///////////////////
//fpd  cache ideal bond angles for fast lookup
class BondAngleDatabase  : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~BondAngleDatabase();
	BondAngleDatabase(Real k_ang);

	// lookup ideal bondangle; insert in DB if not there
	void
	lookup( core::pose::Pose const & pose, core::conformation::Residue const & res,
		   int atm1, int atm2, int atm3, Real &Ktheta, Real &d0);


private:
	//fpd involves a string comparison (restype) .. could be faster
	//std::map< residx_atm_triple, core::Real > bondangles_;
	boost::unordered_map< residx_atm_triple, core::Real > bondangles_;
	boost::unordered_map< residx_atm_triple, core::Real > Kangles_;

	Real k_angle;  // can only be set in constructor
	mm::MMBondAngleLibrary const & mm_bondangle_library_;
};


///////////////////
//fpd  cache ideal bond lengths for fast lookup
class BondLengthDatabase  : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~BondLengthDatabase();
	BondLengthDatabase(Real k_len);

	// lookup ideal bondlength; insert in DB if not there
	void
	lookup( core::pose::Pose const & pose, core::conformation::Residue const & res, int atm1, int atm2, Real &Kd, Real &d0 );

private:
	//fpd involves a string comparison (restype) .. could be faster
	//std::map< residx_atm_pair, core::Real > bondlengths_;
	boost::unordered_map< residx_atm_pair, core::Real > bondlengths_;
	boost::unordered_map< residx_atm_pair, core::Real > Klengths_;

	Real k_length;  // can only be set in constructor
	mm::MMBondLengthLibrary const & mm_bondlength_library_;
};


///////////////////
class CartesianBondedEnergy : public ContextIndependentLRTwoBodyEnergy {
public:
	typedef ContextIndependentLRTwoBodyEnergy  parent;

public:
	CartesianBondedEnergy( methods::EnergyMethodOptions const & options );

	CartesianBondedEnergy( CartesianBondedEnergy const & src );

	~CartesianBondedEnergy();

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	virtual	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;

	///
	virtual	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	virtual	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const ;

	// helper functions to handle improper torsions
	// intra-res impropers
	void
	eval_improper_torsions(
	  conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	// inter-res impropers
	void
	eval_improper_torsions(
	  conformation::Residue const & rsd1,
	  conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	// deriv impropers
 	virtual	void
 	eval_improper_torsions_derivative(
 		id::AtomID const & id,
 		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
 		ScoreFunction const & sfxn,
 		EnergyMap const & emap,
 		Vector & F1,
		Vector & F2
 	) const;


	virtual	void
	eval_intrares_energy(
	  conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap
	) const;

	virtual	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & emap,
		Vector & F1,
		Vector & F2
	) const;


	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

	methods::LongRangeEnergyType
	long_range_type() const;

private:

	mutable BondAngleDatabaseOP db_angle_;
	mutable BondLengthDatabaseOP db_length_;
	mutable TorsionDatabaseOP db_torsion_;

	// options
	bool linear_bonded_potential_;
	core::Real cartbonded_len_, cartbonded_ang_, cartbonded_tors_, cartbonded_proton_;

private:

	virtual
	core::Size version() const;

};

} // namespace methods
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_methods_CartesianBondedEnergy_HH
