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
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_rna_RNA_TorsionPotential_hh
#define INCLUDED_core_scoring_rna_RNA_TorsionPotential_hh

// Unit Headers
#include <core/scoring/rna/RNA_TorsionPotential.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray4D.hh>

//Auto Headers
#include <core/conformation/Residue.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <utility/vector1_bool.hh>
#include <map>


namespace core {
namespace scoring {
namespace rna {

typedef std::pair< id::TorsionID, utility::vector1< constraints::FuncOP > > RNA_SideChainTorsionTether;
typedef std::map< Size, utility::vector1< RNA_SideChainTorsionTether >  > RNA_SideChainTorsionTethers;

enum __TorsionPotential__ { WHATEVER, ALPHA, BETA, GAMMA, DELTA, EPSILON, ZETA, CHI, NU2, NU1, O2H};

class Gaussian_parameter {
public:
	Real amplitude, center, width;

	Gaussian_parameter ( Real const amplitude_in, Real const center_in, Real const width_in ):
		amplitude( amplitude_in ),
		center   ( center_in ),
		width    ( width_in )
	{}

	//
	Gaussian_parameter &
	operator=( Gaussian_parameter const & src )
	{
		amplitude = src.amplitude;
		center = src.center;
		width = src.width;
		return *this;
	}


};

typedef utility::vector1< Gaussian_parameter > Gaussian_parameter_set;

class RNA_TorsionPotential : public utility::pointer::ReferenceCount
{

public:
	RNA_TorsionPotential();
	~RNA_TorsionPotential() {}


    /* Undefinded, comented out to make python bindings complile
	void
	eval_rna_torsion_score_residue(
	  conformation::Residue const & rsd,
		Size const rna_torsion_number,
		Real & score,
		Real & deriv ) const;


	void
	eval_ribose_closure_score( conformation::Residue const & rsd, Real & score ) const;


	void
	eval_ribose_closure_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		Vector & f1,
		Vector & f2 ) const;
		*/

	//	void update_constraints( pose::Pose & pose ) const;

	void
	setup_constraints( pose::Pose & pose, constraints::ConstraintSetOP & rna_torsion_constraints,
	 constraints::ConstraintSetOP & rna_sugar_close_constraints,
	 rna::RNA_SideChainTorsionTethers & rna_side_chain_torsion_tethers ) const;

	Real delta_cutoff() const { return DELTA_CUTOFF_; }

	Gaussian_parameter_set gaussian_parameter_set_alpha() const{ return gaussian_parameter_set_alpha_; }
	Gaussian_parameter_set gaussian_parameter_set_beta() const{ return gaussian_parameter_set_beta_; }
	Gaussian_parameter_set gaussian_parameter_set_gamma() const{ return gaussian_parameter_set_gamma_; }
	Gaussian_parameter_set gaussian_parameter_set_delta_north() const{ return gaussian_parameter_set_delta_north_; }
	Gaussian_parameter_set gaussian_parameter_set_delta_south() const{ return gaussian_parameter_set_delta_south_; }
	Gaussian_parameter_set gaussian_parameter_set_epsilon_north() const{ return gaussian_parameter_set_epsilon_north_; }
	Gaussian_parameter_set gaussian_parameter_set_epsilon_south() const{ return gaussian_parameter_set_epsilon_south_; }
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_sc_minus() const{ return gaussian_parameter_set_zeta_alpha_sc_minus_; }
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_sc_plus() const{ return gaussian_parameter_set_zeta_alpha_sc_plus_; }
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_ap() const{ return gaussian_parameter_set_zeta_alpha_ap_; }
	Gaussian_parameter_set gaussian_parameter_set_chi_north() const{ return gaussian_parameter_set_chi_north_; }
	Gaussian_parameter_set gaussian_parameter_set_chi_south() const{ return gaussian_parameter_set_chi_south_; }
	Gaussian_parameter_set gaussian_parameter_set_nu2_north() const{ return gaussian_parameter_set_nu2_north_; }
	Gaussian_parameter_set gaussian_parameter_set_nu2_south() const{ return gaussian_parameter_set_nu2_south_; }
	Gaussian_parameter_set gaussian_parameter_set_nu1_north() const{ return gaussian_parameter_set_nu1_north_; }
	Gaussian_parameter_set gaussian_parameter_set_nu1_south() const{ return gaussian_parameter_set_nu1_south_; }

	void
	add_sugar_ring_closure_constraints( conformation::Residue const & rsd, constraints::ConstraintSet & cst_set ) const;

private:

	void
	add_sugar_ring_closure_constraints( pose::Pose & pose, constraints::ConstraintSet & cst_set ) const;

	void
	add_rna_torsion_tethers(
	  pose::Pose & pose,
		constraints::ConstraintSet & cst_set,
		rna::RNA_SideChainTorsionTethers & rna_side_chain_torsion_tethers ) const;

	void
	add_o2star_torsion_constraints(
	  pose::Pose & pose,
		constraints::ConstraintSet & cst_set) const;

	void
	add_RNA_torsion_constraint(
			 pose::Pose & pose,
			 Size const i,
			 constraints::ConstraintSet & cst_set,
			 rna::RNA_SideChainTorsionTethers & rna_side_chain_torsion_tethers,
			 Size const rna_torsion_number,
			 Gaussian_parameter_set const & gaussian_parameter_set ) const;


	void
	get_o2star_potential( Real const & torsion_value, Real & score, Real & deriv ) const;

	void
	init_rna_torsion_gaussian_parameters();

	Gaussian_parameter_set gaussian_parameter_set_alpha_;
	Gaussian_parameter_set gaussian_parameter_set_beta_;
	Gaussian_parameter_set gaussian_parameter_set_gamma_;
	Gaussian_parameter_set gaussian_parameter_set_delta_north_;
	Gaussian_parameter_set gaussian_parameter_set_delta_south_;
	Gaussian_parameter_set gaussian_parameter_set_epsilon_north_;
	Gaussian_parameter_set gaussian_parameter_set_epsilon_south_;
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_sc_minus_;
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_sc_plus_;
	Gaussian_parameter_set gaussian_parameter_set_zeta_alpha_ap_;
	Gaussian_parameter_set gaussian_parameter_set_chi_north_;
	Gaussian_parameter_set gaussian_parameter_set_chi_south_;
	Gaussian_parameter_set gaussian_parameter_set_nu2_north_;
	Gaussian_parameter_set gaussian_parameter_set_nu2_south_;
	Gaussian_parameter_set gaussian_parameter_set_nu1_north_;
	Gaussian_parameter_set gaussian_parameter_set_nu1_south_;

	// alpha, beta, gamma, delta, epsilon, zeta
	bool const rna_tight_torsions_;
	Real const DELTA_CUTOFF_;
	Real const scale_rna_torsion_tether_;
	Real const scale_rna_torsion_sd_;

	// Ribose closure
	Size const o4star_index_;
	Size const c1star_index_;
	Size const c2star_index_;
	Size const c4star_index_;
	Distance const o4star_c1star_bond_length_;
	Distance const o4star_c1star_sd_;
	constraints::HarmonicFuncOP o4star_c1star_dist_harm_func_;

	Real const angle_sd_;
	Real const o4star_c1star_c2star_bond_angle_;
	constraints::HarmonicFuncOP o4star_c1star_c2star_angle_harm_func_;
	Real const o4star_c1star_first_base_bond_angle_;
	constraints::HarmonicFuncOP o4star_c1star_first_base_angle_harm_func_;
	Real const c4star_o4star_c1star_bond_angle_;
	constraints::HarmonicFuncOP c4star_o4star_c1star_angle_harm_func_;

	// 2'-OH proton torsion
	Real const o2star_potential_weight_;
	Real const o2star_best_torsion_;

	constraints::FuncOP o2star_dihedral_constraint_func1_;
	constraints::FuncOP o2star_dihedral_constraint_func2_;

	bool const verbose_;

};

}
}
}

#endif
