// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/dna/DNATorsionPotential.hh
/// @brief  DNATorsionPotential potential class delcaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_dna_DNATorsionPotential_HH
#define INCLUDED_core_scoring_dna_DNATorsionPotential_HH

// Unit Headers
#include <core/scoring/dna/DNATorsionPotential.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/func/AmberPeriodicFunc.hh>
#include <core/scoring/func/HarmonicFunc.fwd.hh>
#include <core/scoring/func/HarmonicFunc.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray4D.hh>

namespace core {
namespace scoring {
namespace dna {

//enum{ WHATEVER, ALPHA, BETA, GAMMA, DELTA, EPSILON, ZETA, CHI, NU0, NU1, NU2, NU3, NU4 };

class DNATorsionPotential : public utility::pointer::ReferenceCount
{

public:
	DNATorsionPotential();
	~DNATorsionPotential() {}


    /* Undefinded, comented out to make python bindings complile
	void
	eval_dna_torsion_score_residue(
	  conformation::Residue const & rsd,
		Size const dna_torsion_number,
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
	setup_constraints( pose::Pose & pose, constraints::ConstraintSetOP & dna_torsion_constraints,
	 constraints::ConstraintSetOP & dna_sugar_close_constraints, constraints::ConstraintSetOP & dna_base_distance_constraints ) const;

	Real delta_cutoff() const { return DELTA_CUTOFF_; }

	utility::vector1< constraints::AmberPeriodicFuncOP > & alpha_components() { return alpha_components_; }
	utility::vector1< constraints::AmberPeriodicFuncOP > & beta_components() { return beta_components_; }
	utility::vector1< constraints::AmberPeriodicFuncOP > & gamma_components() { return gamma_components_; }
	utility::vector1< constraints::AmberPeriodicFuncOP > & delta_components() { return delta_components_; }
	utility::vector1< constraints::AmberPeriodicFuncOP > & epsilon_components() { return epsilon_components_; }
	utility::vector1< constraints::AmberPeriodicFuncOP > & zeta_components() { return zeta_components_; }

	utility::vector1< constraints::AmberPeriodicFuncOP > & nu0_components() { return nu0_components_; }
	utility::vector1< constraints::AmberPeriodicFuncOP > & nu1_components() { return nu1_components_; }
	utility::vector1< constraints::AmberPeriodicFuncOP > & nu2_components() { return nu2_components_; }
	utility::vector1< constraints::AmberPeriodicFuncOP > & nu3_components() { return nu3_components_; }
	utility::vector1< constraints::AmberPeriodicFuncOP > & nu4_components() { return nu4_components_; }

	void
	add_sugar_ring_closure_constraints( conformation::Residue const & rsd, constraints::ConstraintSet & cst_set ) const;

private:

	void
	add_sugar_ring_closure_constraints( pose::Pose & pose, constraints::ConstraintSet & cst_set ) const;

	void
	add_dna_base_distance_constraints( pose::Pose & pose, constraints::ConstraintSet & cst_set ) const;

	void
	add_dna_torsion_tethers(
	  pose::Pose & pose,
		constraints::ConstraintSet & cst_set ) const;

	void
	add_DNA_torsion_constraint(
			 pose::Pose & pose,
			 Size const i,
			 constraints::ConstraintSet & cst_set,
			 Size const dna_torsion_number,
			 utility::vector1< constraints::AmberPeriodicFuncOP > const & gaussian_parameter_set ) const;


	bool
	get_atom_ids_by_torsion(
		Size const dna_torsion_number,
		pose::Pose & pose,
		Size const resid,
		id::AtomID & id1,
		id::AtomID & id2,
		id::AtomID & id3,
		id::AtomID & id4 ) const;

	void
	init_dna_torsion_parameters();

	utility::vector1< constraints::AmberPeriodicFuncOP > alpha_components_;
	utility::vector1< constraints::AmberPeriodicFuncOP > beta_components_;
	utility::vector1< constraints::AmberPeriodicFuncOP > gamma_components_;
	utility::vector1< constraints::AmberPeriodicFuncOP > delta_components_;
	utility::vector1< constraints::AmberPeriodicFuncOP > epsilon_components_;
	utility::vector1< constraints::AmberPeriodicFuncOP > zeta_components_;

	utility::vector1< constraints::AmberPeriodicFuncOP > nu0_components_;
	utility::vector1< constraints::AmberPeriodicFuncOP > nu1_components_;
	utility::vector1< constraints::AmberPeriodicFuncOP > nu2_components_;
	utility::vector1< constraints::AmberPeriodicFuncOP > nu3_components_;
	utility::vector1< constraints::AmberPeriodicFuncOP > nu4_components_;

	utility::vector1< std::string > alpha_atom_names_;
	utility::vector1< std::string > beta_atom_names_;
	utility::vector1< std::string > gamma_atom_names_;
	utility::vector1< std::string > delta_atom_names_;
	utility::vector1< std::string > epsilon_atom_names_;
	utility::vector1< std::string > zeta_atom_names_;

	utility::vector1< std::string > nu0_atom_names_;
	utility::vector1< std::string > nu1_atom_names_;
	utility::vector1< std::string > nu2_atom_names_;
	utility::vector1< std::string > nu3_atom_names_;
	utility::vector1< std::string > nu4_atom_names_;

	// alpha, beta, gamma, delta, epsilon, zeta
	Real const DELTA_CUTOFF_;
	Real const scale_dna_torsion_tether_;
	Real const scale_dna_torsion_sd_;

	// Ribose closure
	Distance const c2prime_c3prime_bond_length_;
	Distance const c2prime_c3prime_sd_;
	constraints::HarmonicFuncOP c2prime_c3prime_dist_harm_func_;

	Real const c4prime_c3prime_c2prime_bond_angle_;
	constraints::HarmonicFuncOP c4prime_c3prime_c2prime_angle_harm_func_;
	Real const o3prime_c3prime_c2prime_bond_angle_;
	constraints::HarmonicFuncOP o3prime_c3prime_c2prime_angle_harm_func_;
	Real const c3prime_c2prime_c1prime_bond_angle_;
	constraints::HarmonicFuncOP c3prime_c2prime_c1prime_angle_harm_func_;

	bool const verbose_;

};

}
}
}

#endif
