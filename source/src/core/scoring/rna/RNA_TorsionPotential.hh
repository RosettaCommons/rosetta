// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/RNA_TorsionPotential.hh
/// @brief  RNA_TorsionPotential potential class delcaration
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_RNA_TorsionPotential_HH
#define INCLUDED_core_scoring_rna_RNA_TorsionPotential_HH

// Unit Headers
#include <core/scoring/rna/RNA_TorsionPotential.fwd.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.fwd.hh>

// Project Headers
#include <core/chemical/rna/RNA_FittedTorsionInfo.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/SumFunc.fwd.hh>
#include <core/scoring/func/SumFunc.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>


#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <utility/vector1.hh>
#include <string>


namespace core {
namespace scoring {
namespace rna {

class RNA_TorsionPotential : public utility::pointer::ReferenceCount
{

public:

	RNA_TorsionPotential( RNA_EnergyMethodOptions const & options );

	virtual ~RNA_TorsionPotential() ; // auto-removing definition from header{}


	//	void update_constraints( pose::Pose & pose ) const;
	Real
	eval_intrares_energy(
											 core::conformation::Residue const & rsd,
											 pose::Pose const & pose
											 );


	Real
	residue_pair_energy(
											core::conformation::Residue const & rsd1,
											core::conformation::Residue const & rsd2,
											pose::Pose const & pose ) const;

	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


	//Real
	//	compute_torsion_potential( Size const & torsion_number, Real const & value, Real const & delta, Real const & next_alpha ) const;

	Real
	intrares_side_chain_score() const{ return intrares_side_chain_score_; }

	void set_verbose( bool const setting ){ verbose_ = setting; }

private:

	bool
	check_intra_residue( id::TorsionID const & torsion_id, pose::Pose const & pose, Size const seqpos ) const;

	void
	init_potentials_from_rna_torsion_database_files();

	void
	initialize_potential_from_file( core::scoring::func::FuncOP & func,
																	std::string const & filename );

	void
	init_fade_functions();

	bool
	get_f1_f2( core::id::TorsionID const & torsion_id,
					 core::pose::Pose const & pose, core::id::AtomID const & id, Vector & f1, Vector & f2 ) const;

	std::string path_to_torsion_files_;

	// alpha, beta, gamma, delta, epsilon, zeta
	// unused bool const rna_tight_torsions_;
	Real const delta_fade_;
	Real const alpha_fade_;

	core::scoring::func::FuncOP alpha_potential_, beta_potential_, gamma_potential_, delta_north_potential_,
	delta_south_potential_, epsilon_north_potential_, epsilon_south_potential_, zeta_alpha_sc_minus_potential_,
	zeta_alpha_sc_plus_potential_, zeta_alpha_ap_potential_, nu2_north_potential_, nu2_south_potential_,
	nu1_north_potential_, nu1_south_potential_, chi_north_potential_others_, chi_south_potential_others_,
	chi_north_potential_guanosine_, chi_south_potential_guanosine_, chi_purine_north_potential_, chi_purine_south_potential_,
		chi_pyrimidine_north_potential_, chi_pyrimidine_south_potential_, o2h_north_potential_, o2h_south_potential_, chi_potential_syn_guanosine_bonus_;

	core::scoring::func::FuncOP fade_delta_north_, fade_delta_south_;
	core::scoring::func::FuncOP fade_alpha_sc_minus_, fade_alpha_sc_plus_;
	core::scoring::func::SumFuncOP fade_alpha_ap_;

	bool const skip_chainbreak_torsions_;
	bool verbose_;
	bool use_new_potential_;
	bool const use_2prime_OH_potential_;
	Real const syn_G_potential_bonus_;
	chemical::rna::RNA_FittedTorsionInfoOP rna_fitted_torsion_info_;

	Real intrares_side_chain_score_;

};

} //rna
} //scoring
} //core

#endif
