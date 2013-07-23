// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// @file:   core/scoring/facts/FACTSPotential.hh
// @brief:  Header-file for FACTSPotential class declaration
//          FACTS: Fast Analytical Continuum Treatment of Solvation by URS HABERTHUR and AMEDEO CAFLISCH
// @author: Hahnbeom Park

#ifndef INCLUDED_core_scoring_facts_FACTSPotential_HH
#define INCLUDED_core_scoring_facts_FACTSPotential_HH

// Unit headers
#include <core/scoring/facts/FACTSResidue.hh>
#include <core/scoring/facts/FACTSPose.hh>
#include <core/scoring/facts/FACTSPotential.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/RotamerSetBase.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <iostream>

using namespace std;

namespace core {
namespace scoring {

/**************************************************************************************************/
/*                                                                                                */
/*    @brief: The FACTSPotential class provides all the functions, constants, and parameters      */
/*            common to all atoms required to calculate the free energy of solvation of a         */
/* (macro)molecule embedded in a continuum solvent using FACTS method               */
/*                                                                                                */
/**************************************************************************************************/

class FACTSPotential: public utility::pointer::ReferenceCount {

public:
	typedef conformation::Residue Residue;

public:
	FACTSPotential();

	void set_default();

	void setup_for_scoring(pose::Pose & pose, bool const & packing ) const;

	void setup_for_derivatives(pose::Pose & pose) const;

	void setup_for_packing(
												 pose::Pose & pose,
												 utility::vector1< bool > const & repacking_residues ) const;

	void update_residue_for_packing( pose::Pose & pose,
																	 Size const seqpos
																	 ) const;

	void get_rotamers_born_radii(pose::Pose const & pose, conformation::RotamerSetBase & rotamer_set) const;

	void evaluate_polar_energy( Residue const & rsd1,
															FACTSResidueInfo const & facts1,
															Residue const & rsd2,
															Real & E_elec,
															Real & E_solv_self,
															Real & E_solv_pair
															) const;

	Real evaluate_nonpolar_energy( Residue const & rsd1,
																 FACTSResidueInfo const & facts1,
																 Residue const & rsd2
																 ) const;

	void evaluate_context_change_for_packing(
						 Residue const & rsd1_ref,
						 Residue const & rsd1,
						 FACTSResidueInfo const & facts1,
						 Residue const & rsd2_ref,
						 Residue const & rsd2,
						 FACTSResidueInfo const & facts2,
						 utility::vector1< Real > & dBRi1,
						 utility::vector1< Real > & dBRi2,
						 utility::vector1< Real > & dSAi1,
						 utility::vector1< Real > & dSAi2
						 ) const;

	void evaluate_polar_otf_energy(Residue const & rsd1,
																 FACTSResidueInfo const & facts1,
																 Residue const & rsd2,
																 FACTSResidueInfo const & facts2,
																 Real & E_elec,
																 Real & E_solv_self,
																 Real & E_solv_pair
																 ) const;

	void eval_atom_polar_derivative(
					id::AtomID const & id,
					Real const weight_elec,
					Real const weight_solv,
					pose::Pose const & pose,
					kinematics::DomainMap const &, //domain_map,
					bool const, //exclude_DNA_DNA,
					Vector & F1,
					Vector & F2
					) const;

	void eval_atom_nonpolar_derivative(
					id::AtomID const & id,
					Real const weight,
					pose::Pose const & pose,
					kinematics::DomainMap const &, //domain_map,
					bool const, //exclude_DNA_DNA,
					Vector & F1,
					Vector & F2
					) const;

	void get_single_rotamer_born_radii(
																	Residue const & rsd1,
																	pose::Pose const & pose,
																	FACTSPoseInfo const & facts_info,
																	FACTSResidueInfo & facts1
																	) const;

	Real polar_energy_pack_corrector(
																	 Residue const & ref_rsd,
																	 Residue const & rsd,
																	 FACTSResidueInfo const & facts_info
																	 ) const;

private:
	void res_res_burial(
											Residue const & rsd1,
											FACTSResidueInfo & facts1,
											Residue const & rsd2,
											FACTSResidueInfo const & facts2
											) const;

	void res_res_burial_for_scoring(
												Residue const & rsd1,
												FACTSResidueInfo & facts1,
												Residue const & rsd2,
												FACTSResidueInfo & facts2
												) const;

	void get_self_terms( FACTSRsdTypeInfoCOP factstype1,
											 FACTSResidueInfo & facts1,
											 bool const packing
											 ) const;

	void atompair_scale( FACTSRsdTypeInfoCOP factstype1, 
											 scoring::etable::count_pair::CountPairFunctionCOP cpfxn14,
											 conformation::Residue const &rsd1,
											 conformation::Residue const &rsd2,
											 Size const atm1,
											 Size const atm2,
											 Real &scale_solv,
											 Real &scale_elec,
											 bool &self_pair,
											 bool const same_res,
											 bool const adjacent ) const;

	void calculate_GBpair_fast(
													Residue const & rsd1,
													Residue const & rsd2,
													FACTSResidueInfo & facts1,
													FACTSResidueInfo & facts2
													) const;

	void calculate_GBpair_exact(
													Residue const & rsd1,
													Residue const & rsd2,
													FACTSResidueInfo & facts1,
													FACTSResidueInfo & facts2
													) const;


	void atom_atom_context_derivative( FACTSResidueInfo & facts1,
																		 FACTSResidueInfo & facts2,
																		 Size const & atm1,
																		 Size const & atm2,
																		 Vector const & dxyz,
																		 bool const full_update
																		 ) const;

	void get_template_born_radii(
														pose::Pose const & pose,
														FACTSPoseInfo & facts_info
														) const;

	// Accessors
	inline Real Tau() const{	return Tau_; }
	inline Real inv_die() const{	return inv_die_; }
	inline Real Kappa() const { return Kappa_; }
	inline Real MultiplicitiveFactor() const { return MultiplicitiveFactor_; };
	inline Real GBPair_cut() const { return GBpair_cut_; };
	inline Real adjbb_elec_scale( Size const i ) const { return adjbb_elec_scale_[i]; }
	inline Real adjbb_solv_scale( Size const i ) const { return adjbb_solv_scale_[i]; }
	inline Real adjsc_elec_scale( Size const i ) const { return adjsc_elec_scale_[i]; }
	inline Real adjsc_solv_scale( Size const i ) const { return adjsc_solv_scale_[i]; }

private: //list of private variables and parameters for the FACTS method common to all atoms

	// Map storing parameters for residue types
	mutable FACTSRsdTypeMap FACTSrsdtypemap_;

	bool options_registered_;

	Real MultiplicitiveFactor_;
	Real inv_die_;
	Real Tau_;
	Real Kappa_;
	Real GBpair_cut_;
	bool do_apprx_;

	utility::vector1< Real > adjbb_elec_scale_;
	utility::vector1< Real > adjbb_solv_scale_;
	utility::vector1< Real > adjsc_elec_scale_;
	utility::vector1< Real > adjsc_solv_scale_;

	Real saltbridge_correction_;
	Real dshift2_;
	Real dshift2_saltbridge_;

	// Deprecated
	//Real elec_sh_exponent_;
	//Real selfenergy_scale_;
	//Real intrares_scale_;
	//Real min_dis_;

	// Below are not being used
	Real cut_off_born_radius_; //The cut off used for calculating born radius
	Real extra_cut_off_self_;
	Real extra_cut_off_interaction_;
	Real dummy_radius_;
	Real dummy_scale_;
	Real dummy_distance_; // also implicitly defined by the gb placeholder params file

};

} // scoring
} // core

#endif
