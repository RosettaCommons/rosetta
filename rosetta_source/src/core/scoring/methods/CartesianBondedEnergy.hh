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
#include <core/scoring/methods/CartBondedParameters.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

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


typedef boost::tuples::tuple< std::string, std::string, std::string, std::string, std::string > atm_name_quad;
typedef boost::tuples::tuple< std::string, std::string, std::string, std::string > atm_name_triple;
typedef boost::tuples::tuple< std::string, std::string, std::string > atm_name_pair;
typedef boost::tuples::tuple< std::string, std::string > atm_name_single;

namespace boost {
namespace tuples {

std::size_t hash_value(atm_name_quad const& e); 
std::size_t hash_value(atm_name_triple const& e); 
std::size_t hash_value(atm_name_pair const& e); 
std::size_t hash_value(atm_name_single const& e); 
bool operator==(atm_name_quad const& a,atm_name_quad const& b); 
bool operator==(atm_name_triple const& a,atm_name_triple const& b); 
bool operator==(atm_name_pair const& a,atm_name_pair const& b); 
bool operator==(atm_name_single const& a,atm_name_single const& b); 

}
}

namespace core {
namespace scoring {
namespace methods {




////////////////////
//fpd
//  Database stores all ideal parameters
class IdealParametersDatabase  : public utility::pointer::ReferenceCount {
public:
	IdealParametersDatabase(Real k_len, Real k_ang, Real k_tors, Real k_tors_prot, Real k_tors_improper);

	void init(Real k_len, Real k_ang, Real k_tors, Real k_tors_prot, Real k_tors_improper);

	CartBondedParametersCAP
	lookup_torsion(
		core::conformation::Residue const & res,
		std::string atm1_name, std::string atm2_name, std::string atm3_name, std::string atm4_name );

	// needs both names (for keying off databases) and indices (for building those not found from ideal)
	CartBondedParametersCAP													
	lookup_angle(
		core::conformation::Residue const & res,
		bool pre_proline,
		std::string atm1_name, std::string atm2_name, std::string atm3_name,
		int atm1idx, int atm2idx, int atm3idx);

	// needs both names (for keying off databases) and indices (for building those not found from ideal)
	CartBondedParametersCAP													
	lookup_length(
		core::conformation::Residue const & res,
		bool pre_proline,
		std::string atm1_name, std::string atm2_name,
		int atm1idx, int atm2idx);

	// old-style interface to database
	void													
	lookup_torsion_legacy( core::chemical::ResidueType const & restype,
	        int atm1, int atm2, int atm3, int atm4, Real &Kphi, Real &phi0, Real &phi_step );

	// old-style interface to database
	void
	lookup_angle_legacy( core::pose::Pose const & pose, core::conformation::Residue const & res,
		   int atm1, int atm2, int atm3, Real &Ktheta, Real &d0);

	// old-style interface to database
	void
	lookup_length_legacy( core::pose::Pose const & pose, core::conformation::Residue const & res, int atm1, int atm2, Real &Kd, Real &d0 );

	bool bbdep_bond_params() { return bbdep_bond_params_; }
	bool bbdep_bond_devs() { return bbdep_bond_devs_; }

private:
	// helper functions: find the ideal values by constructing from Rosetta's params file
	void			
	lookup_bondangle_buildideal(
		core::conformation::Residue const & res,
		int atm1, int atm2, int atm3, Real &Ktheta, Real &theta0 );

	void
	lookup_bondlength_buildideal(
		core::conformation::Residue const & res,
		int atm1, int atm2, Real &Kd, Real &d0 );

	// another helper function: read backbone dependent db files
	void
	read_bbdep_table(
		std::string filename,
		boost::unordered_map< atm_name_single, CartBondedParametersOP > &bondlengths,
		boost::unordered_map< atm_name_pair, CartBondedParametersOP > &bondangles,
		std::string res );

	
	// defaults (they should be rarely used as everything should be in the DB now)
	Real k_length_, k_angle_, k_torsion_, k_torsion_proton_, k_torsion_improper_;

	// backbone-independent parameters (keyed on atom names)
	boost::unordered_map< atm_name_quad, 	CartBondedParametersOP > torsions_indep_;
	boost::unordered_map< atm_name_triple, CartBondedParametersOP > bondangles_indep_;
	boost::unordered_map< atm_name_pair, CartBondedParametersOP > bondlengths_indep_;

	// backbone-dependent parameter sets
	boost::unordered_map< atm_name_pair, CartBondedParametersOP >
		bondangles_bbdep_def_, bondangles_bbdep_pro_, bondangles_bbdep_valile_, bondangles_bbdep_prepro_, bondangles_bbdep_gly_;
	boost::unordered_map< atm_name_single, CartBondedParametersOP >
		bondlengths_bbdep_def_, bondlengths_bbdep_pro_, bondlengths_bbdep_valile_, bondlengths_bbdep_prepro_, bondlengths_bbdep_gly_;

	// options
	bool bbdep_bond_params_, bbdep_bond_devs_;
};




///////////////////
///
/// the energy method
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

	/// here, rsd1.seqpos() < rsd2.seqpos()
	void
	residue_pair_energy_sorted(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sf,
		EnergyMap & emap
	) const;

	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const;

	//fpd because of separate pre-proline distributions
	//    this function must be called from residue_pair_energy
	//    rather than defined in intrares energy
	virtual	void
	eval_singleres_energy(
	  conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap & emap,
		bool preproline
	) const;

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

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResSingleMinimizationData const &,
		ResSingleMinimizationData const &,
		ResPairMinimizationData const & min_data,
		pose::Pose const &,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	// improper derivatives
 	virtual	void
 	eval_improper_torsions_derivative(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
 	) const;

	// improper derivatives
 	virtual	void
 	eval_improper_torsions_derivative(
		conformation::Residue const & rsd,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r_atom_derivs
 	) const;

	// single-residue derivatives
	void
	eval_singleres_derivatives(
		conformation::Residue const & rsd,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r_atom_derivs,
		bool preproline) const;

	// dof (bbdep) derivatives
	virtual	Real
	eval_intraresidue_dof_derivative(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
	) const;

	virtual	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual	bool
	defines_intrares_dof_derivatives( pose::Pose const & p ) const { return true; }

	//fpd  use the new minimizer interface
	virtual bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

	methods::LongRangeEnergyType
	long_range_type() const;

private:

	// the ideal parameter database
	static IdealParametersDatabaseOP db_;

	// option
	bool linear_bonded_potential_;

	virtual
	core::Size version() const;

};

} // namespace methods
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_methods_CartesianBondedEnergy_HH
