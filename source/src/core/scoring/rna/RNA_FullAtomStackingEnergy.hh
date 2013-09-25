// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_FullAtomStacking.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Phil Bradley
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_FullAtomStackingEnergy_HH
#define INCLUDED_core_scoring_rna_RNA_FullAtomStackingEnergy_HH

// Unit Headers
#include <core/scoring/rna/RNA_FullAtomStackingEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.fwd.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.fwd.hh>



namespace core {
namespace scoring {
namespace rna {

typedef  numeric::xyzMatrix< Real > Matrix;

///

class RNA_FullAtomStackingEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;

public:

	///
	RNA_FullAtomStackingEnergy();


	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & scfxn) const;
    
    void
    setup_for_minimizing(
        pose::Pose & pose,
        ScoreFunction const & sfxn,
        kinematics::MinimizerMapBase const & min_map
    ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & scfxn) const;

	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	virtual
	void
	eval_intrares_energy(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &,
		EnergyMap &
	) const {}


 	virtual
 	void
 	eval_atom_derivative(
 		id::AtomID const & atom_id,
 		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
 		ScoreFunction const &,
 		EnergyMap const & weights,
 		Vector & F1,
 		Vector & F2
 	) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &// totals
	) const;

	virtual
	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {}

	/// @brief Interface function for class NeighborList.
	etable::count_pair::CountPairFunctionCOP
	get_intrares_countpair(
		conformation::Residue const &,
		pose::Pose const &,
		ScoreFunction const &
	) const;

	/// @brief Interface function for class NeighborList.
	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		Size const,
		Size const,
		pose::Pose const &,
		ScoreFunction const &
	) const;

	etable::count_pair::CountPairFunctionCOP
	get_count_pair_function(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2
	) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	//	Real
	//	get_fa_stack_score( Distance const dist, Real const cos_kappa ) const;

	Real
	get_fa_stack_score( Vector const r_vec, Matrix const M_i ) const;

	Vector
	get_fa_stack_deriv( Vector const r_vec, Matrix const M_i ) const;

	Real
	residue_pair_energy_one_way(
															conformation::Residue const & rsd1,
															conformation::Residue const & rsd2,
															pose::Pose const & pose,
															Real & score_aro
															) const;

	bool
	check_base_base_OK(
	   conformation::Residue const & rsd1,
		 conformation::Residue const & rsd2,
		 Size const & m, Size const & n ) const;

	bool
	is_aro(
				 conformation::Residue const & rsd1,
				 Size const & m ) const;

	virtual
	core::Size version() const;

  Real const prefactor_;
  Distance const full_stack_cutoff_;
  Distance const dist_cutoff_;
  Real const dist_cutoff2_;
	bool const base_base_only_;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
