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
/// @author Frank DiMaio
/// @author Andrew Leaver-Fay optimized the code a bit

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
#include <map>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/fixedsizearray1.hh>
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


class ResidueCartBondedParameters : public utility::pointer::ReferenceCount {
public:
	typedef utility::fixedsizearray1< Size, 2 > Size2;
	typedef utility::fixedsizearray1< Size, 3 > Size3;
	typedef utility::fixedsizearray1< Size, 4 > Size4;
	typedef std::pair< Size2, CartBondedParametersCOP > length_parameter;
	typedef std::pair< Size3, CartBondedParametersCOP > angle_parameter;
	typedef std::pair< Size4, CartBondedParametersCOP > torsion_parameter;

public:
	ResidueCartBondedParameters();
	virtual ~ResidueCartBondedParameters();

	void add_length_parameter(  Size2 atom_inds, CartBondedParametersCOP );
	void add_angle_parameter(   Size3 atom_inds, CartBondedParametersCOP );
	void add_torsion_parameter( Size4 atom_inds, CartBondedParametersCOP );
	void add_improper_torsion_parameter( Size4 atom_inds, CartBondedParametersCOP );
	void add_bbdep_length_parameter(  Size2 atom_inds, CartBondedParametersCOP );
	void add_bbdep_angle_parameter(   Size3 atom_inds, CartBondedParametersCOP );
	void add_lower_connect_angle_params( Size3 atom_inds, CartBondedParametersCOP );
	void add_upper_connect_angle_params( Size3 atom_inds, CartBondedParametersCOP );

	void bb_N_index( Size index );
	void bb_CA_index( Size index );
	void bb_C_index( Size index );
	void bb_O_index( Size index );
	void bb_H_index( Size index );
	void pro_CD_index( Size index );

	void ca_cprev_n_h_interres_torsion_params( CartBondedParametersCOP );
	void oprev_cprev_n_h_interres_torsion_params( CartBondedParametersCOP );
	void ca_nnext_c_o_interres_torsion_params( CartBondedParametersCOP );
	void pro_cd_cprev_n_ca_interres_torsion_params( CartBondedParametersCOP );
	void cprev_n_bond_length_params( CartBondedParametersCOP );

	utility::vector1< length_parameter > const &
	length_parameters() const {
		return length_params_;
	}

	utility::vector1< angle_parameter > const &
	angle_parameters() const {
		return angle_params_;
	}

	utility::vector1< torsion_parameter > const &
	torsion_parameters() const {
		return torsion_params_;
	}

	/// @brief Exactly the same as proper torsion parameters, but parceled out
	/// into their own section so that debugging information can be given for
	/// these torsions in particular.
	utility::vector1< torsion_parameter > const &
	improper_torsion_parameters() const {
		return improper_torsion_params_;
	}

	/// @brief just the list of length parameters that are dependent on phi and psi; used for calculating dE/dphi and dE/dpsi
	utility::vector1< length_parameter > const &
	bbdep_length_parameters() const {
		return bbdep_length_params_;
	}

	/// @brief just the list of angle parameters that are dependent on phi and psi; used for calculating dE/dphi and dE/dpsi
	utility::vector1< angle_parameter > const &
	bbdep_angle_parameters() const {
		return bbdep_angle_params_;
	}

	utility::vector1< angle_parameter > const &
	lower_connect_angle_params() const {
		return lower_connect_angle_params_;
	}

	utility::vector1< angle_parameter > const &
	upper_connect_angle_params() const {
		return upper_connect_angle_params_;
	}

	Size bb_N_index()  const { return bb_N_index_;  }
	Size bb_CA_index() const { return bb_CA_index_; }
	Size bb_C_index()  const { return bb_C_index_;  }
	Size bb_O_index()  const { return bb_O_index_;  }
	Size bb_H_index()  const { return bb_H_index_;  }
	Size pro_CD_index( )  const { return pro_CD_index_;  }

	CartBondedParametersCOP
	ca_cprev_n_h_interres_torsion_params() const {
		return ca_cprev_n_h_interres_torsion_params_;
	}

	CartBondedParametersCOP
	oprev_cprev_n_h_interres_torsion_params() const {
		return oprev_cprev_n_h_interres_torsion_params_;
	}

	CartBondedParametersCOP
	ca_nnext_c_o_interres_torsion_params() const {
		return ca_nnext_c_o_interres_torsion_params_;
	}

	CartBondedParametersCOP
	pro_cd_cprev_n_ca_interres_torsion_params() const {
		return pro_cd_cprev_n_ca_interres_torsion_params_;
	}

	CartBondedParametersCOP
	cprev_n_bond_length_params() const {
		return cprev_n_bond_length_params_;
	}


private:
	utility::vector1< length_parameter  > length_params_;
	utility::vector1< angle_parameter   > angle_params_;
	utility::vector1< torsion_parameter > torsion_params_;
	utility::vector1< torsion_parameter > improper_torsion_params_;
	utility::vector1< length_parameter  > bbdep_length_params_;
	utility::vector1< angle_parameter   > bbdep_angle_params_;

	/// For amino acids only: if they have a lower connection,
	/// then what are the angle parameters for Cprev-at1-at2 for all
	/// atoms at2 bonded to lower-connect-atom at1?
	utility::vector1< angle_parameter   > lower_connect_angle_params_;

	/// For amino acids only: if they have an upper connection,
	/// then what are the angle parameters for Nnext-at1-at2 for all
	/// atoms at2 bonded to upper-connect-atom at1?
	utility::vector1< angle_parameter   > upper_connect_angle_params_;


	Size bb_N_index_;
	Size bb_CA_index_;
	Size bb_C_index_;
	Size bb_O_index_;
	Size bb_H_index_;
	Size pro_CD_index_;

	CartBondedParametersCOP ca_cprev_n_h_interres_torsion_params_;
	CartBondedParametersCOP oprev_cprev_n_h_interres_torsion_params_;
	CartBondedParametersCOP ca_nnext_c_o_interres_torsion_params_;
	CartBondedParametersCOP pro_cd_cprev_n_ca_interres_torsion_params_;
	CartBondedParametersCOP cprev_n_bond_length_params_;

};

////////////////////
//fpd
//  Database stores all ideal parameters
class IdealParametersDatabase  : public utility::pointer::ReferenceCount {
public:
	IdealParametersDatabase(Real k_len, Real k_ang, Real k_tors, Real k_tors_prot, Real k_tors_improper);

	void init(Real k_len, Real k_ang, Real k_tors, Real k_tors_prot, Real k_tors_improper);

	CartBondedParametersCOP
	lookup_torsion(
		core::chemical::ResidueType const & rsd_type,
		std::string const & atm1_name,
		std::string const & atm2_name,
		std::string const & atm3_name,
		std::string const & atm4_name
	);

	// needs both names (for keying off databases) and indices (for building those not found from ideal)
	CartBondedParametersCOP
	lookup_angle(
		core::chemical::ResidueType const & rsd_type,
		bool pre_proline,
		std::string const & atm1_name,
		std::string const & atm2_name,
		std::string const & atm3_name,
		int atm1idx,
		int atm2idx,
		int atm3idx
	);

	// needs both names (for keying off databases) and indices (for building those not found from ideal)
	CartBondedParametersCOP
	lookup_length(
		core::chemical::ResidueType const & rsd_type,
		bool pre_proline,
		std::string const & atm1_name,
		std::string const & atm2_name,
		int atm1idx,
		int atm2idx
	);

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

	/// @brief Return a list of all the bond lengths, bond angles, and bond torsions
	/// for a single residue type.  This list is constructed lazily as required.
	/// (This may cause thread safety issues!).
	ResidueCartBondedParameters const &
	parameters_for_restype(
		core::chemical::ResidueType const & restype,
		bool prepro
	);

private:
	// helper functions: find the ideal values by constructing from Rosetta's params file
	void
	lookup_bondangle_buildideal(
		core::chemical::ResidueType const & restype,
		int atm1,
		int atm2,
		int atm3,
		Real &Ktheta,
		Real &theta0
	);

	void
	lookup_bondlength_buildideal(
		core::chemical::ResidueType const & restype,
		int atm1,
		int atm2,
		Real &Kd,
		Real &d0
	);

	// read bb indep tables
	void read_length_database(std::string infile);
	void read_angle_database(std::string infile);
	void read_torsion_database(std::string infile);

	// another helper function: read backbone dependent db files
	void
	read_bbdep_table(
		std::string filename,
		boost::unordered_map< atm_name_single, CartBondedParametersOP > &bondlengths,
		boost::unordered_map< atm_name_pair, CartBondedParametersOP > &bondangles,
		std::string res );

	void
	create_parameters_for_restype(
		core::chemical::ResidueType const & restype,
		bool prepro
	);

private:

	// defaults (they should be rarely used as everything should be in the DB now)
	Real k_length_, k_angle_, k_torsion_, k_torsion_proton_, k_torsion_improper_;

	// backbone-independent parameters (keyed on atom names)
	boost::unordered_map< atm_name_quad,  CartBondedParametersOP > torsions_indep_;
	boost::unordered_map< atm_name_triple, CartBondedParametersOP > bondangles_indep_;
	boost::unordered_map< atm_name_pair, CartBondedParametersOP > bondlengths_indep_;

	// backbone-dependent parameter sets
	boost::unordered_map< atm_name_pair, CartBondedParametersOP >
		bondangles_bbdep_def_, bondangles_bbdep_pro_, bondangles_bbdep_valile_, bondangles_bbdep_prepro_, bondangles_bbdep_gly_;
	boost::unordered_map< atm_name_single, CartBondedParametersOP >
		bondlengths_bbdep_def_, bondlengths_bbdep_pro_, bondlengths_bbdep_valile_, bondlengths_bbdep_prepro_, bondlengths_bbdep_gly_;

	// options
	bool bbdep_bond_params_, bbdep_bond_devs_;

	// per residue-type data
	std::map< chemical::ResidueType const *, ResidueCartBondedParametersOP > prepro_restype_data_;
	std::map< chemical::ResidueType const *, ResidueCartBondedParametersOP > nonprepro_restype_data_;

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

	virtual void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & sfxn ) const;

	/// @brief Idealize the virtual NV atom of every proline in the pose. This
	///prevents innacurate pro-close scores when switching between cartesian
	///and non-cartesian score functions.
	void
	idealize_proline_nvs(
		pose::Pose & pose
	) const;

	virtual
	bool
	defines_residue_pair_energy(
		pose::Pose const & pose,
		Size res1,
		Size res2
	) const;


	virtual void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const;

	void
	eval_intrares_derivatives(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & res_data_cache,
		pose::Pose const & pose,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & atom_derivs
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

	// dof (bbdep) derivatives
	virtual
	Real
	eval_intraresidue_dof_derivative(
		conformation::Residue const & rsd,
		ResSingleMinimizationData const & min_data,
		id::DOF_ID const & dof_id,
		id::TorsionID const & torsion_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return true; }

	virtual
	bool
	defines_intrares_dof_derivatives( pose::Pose const & ) const { return true; }

	//fpd  use the new minimizer interface
	virtual
	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

	methods::LongRangeEnergyType
	long_range_type() const;

private:

	//////////////////////////////
	/// Score evaluation methods
	//////////////////////////////

	/// @brief here, rsd1.seqpos() < rsd2.seqpos()
	void
	residue_pair_energy_sorted(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const & sf,
		EnergyMap & emap
	) const;

	/// Methods for intra-residue energies

	//@brief because of separate pre-proline distributions
	//    this function must be called from residue_pair_energy
	//    rather than defined in intrares energy
	void
	eval_singleres_energy(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & rsdparams,
		Real phi,
		Real psi,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief helper function to handle intrares improper torsions
	void
	eval_singleres_improper_torsion_energies(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;


	/// @brief helper function to handle intrares bond torsions
	void
	eval_singleres_torsion_energies(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief helper function to handle intrares bond angles
	void
	eval_singleres_angle_energies(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief helper function to handle intrares bond lengths
	void
	eval_singleres_length_energies(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// Methods for inter-residue energies

	/// @brief Evaluate all the inter-residue components
	void
	eval_residue_pair_energies(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi1,
		Real psi1,
		Real phi2,
		Real psi2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/// @brief Evaluate all inter
	void
	eval_interresidue_angle_energies_two_from_rsd1(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi1,
		Real psi1,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	void
	eval_interresidue_angle_energies_two_from_rsd2(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi2,
		Real psi2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	void
	eval_interresidue_bond_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi2,
		Real psi2,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;


	void
	eval_improper_torsions(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	/////////////////////////////////
	/// Derivative evaluation methods
	/////////////////////////////////

	// single-residue derivatives

	/// @brief evaluate all intra-residue derivatives
	void
	eval_singleres_derivatives(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real phi,
		Real psi,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r_atom_derivs
	) const;

	/// @brief evaluate intra-residue (proper) torsion derivatives
	void
	eval_singleres_torsion_derivatives(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r_atom_derivs
	) const;

	/// @brief evaluate intra-residue angle derivatives
	void
	eval_singleres_angle_derivatives(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r_atom_derivs
	) const;

	/// @brief evaluate intra-residue bond-length derivatives
	void
	eval_singleres_length_derivatives(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		Real const phi,
		Real const psi,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r_atom_derivs
	) const;


	/// @brief evaluate intra-residue improper torsion derivatives
	void
	eval_singleres_improper_torsions_derivatives(
		conformation::Residue const & rsd,
		ResidueCartBondedParameters const & resparams,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r_atom_derivs
	) const;

	// residue-pair derivatives

	/// @brief evaluate inter-residue angle derivatives where
	/// two of the atoms defining the angle are from rsd1
	void
	eval_interresidue_angle_derivs_two_from_rsd1(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi1,
		Real psi1,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	/// @brief evaluate inter-residue angle derivatives where
	/// two of the atoms defining the angle are from rsd2
	void
	eval_interresidue_angle_derivs_two_from_rsd2(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi2,
		Real psi2,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	/// @brief evaluate inter-residue bond-length derivatives
	void
	eval_interresidue_bond_length_derivs(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & rsd1params,
		ResidueCartBondedParameters const & rsd2params,
		Real phi2,
		Real psi2,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	/// @brief evaluate inter-residue improper torsion derivatives
	void
	eval_improper_torsion_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		ResidueCartBondedParameters const & res1params,
		ResidueCartBondedParameters const & res2params,
		EnergyMap const & weights,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	////////////////////////////////////////////////////////////////////
	/// Common to evaluating the score for torsions, angles, and lengths
	////////////////////////////////////////////////////////////////////


	/// @brief Evaluate either the harmonic or linearized-harmonic energy
	/// given by either:
	/// score = 0.5 * K * (val-val0)^2
	/// or
	/// score = 0.5 * K * (val-val0)^2 if std::abs(val-val0) < 1
	///       = 0.5 * K * std::abs(val-val0) otherwise
	Real eval_score( Real val, Real K, Real val0 ) const;

	/// @brief Evaluate the derivative for a

private:

	// the ideal parameter database
	static IdealParametersDatabaseOP db_;

	// option
	bool linear_bonded_potential_;

	virtual
	core::Size version() const;

	std::string pro_nv_;

};

} // namespace methods
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_methods_CartesianBondedEnergy_HH
