// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/LK_PolarNonPolarEnergy.hh
/// @brief  LK Solvation using hemisphere culling class declaration
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_methods_LK_PolarNonPolarEnergy_hh
#define INCLUDED_core_scoring_methods_LK_PolarNonPolarEnergy_hh

// Unit Headers
#include <core/scoring/methods/LK_PolarNonPolarEnergy.fwd.hh>

// Package headers
#include <core/conformation/Atom.fwd.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
#include <core/scoring/etable/Etable.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray3D.fwd.hh>


namespace core {
namespace scoring {
namespace methods {


class LK_PolarNonPolarEnergy : public ContextIndependentTwoBodyEnergy  {
public:
	typedef ContextIndependentTwoBodyEnergy  parent;

public:

	LK_PolarNonPolarEnergy( etable::Etable const & etable_in,
										 bool const analytic_etable_evaluation );

	/// clone
	virtual
	EnergyMethodOP
	clone() const;

	LK_PolarNonPolarEnergy( LK_PolarNonPolarEnergy const & src );

	virtual
 	void
	setup_for_derivatives(
		pose::Pose & pose,
		ScoreFunction const & scfxn
	) const;

    virtual
    void
    setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;
    
    virtual
    void
    setup_for_minimizing(
        pose::Pose & pose,
        ScoreFunction const & sfxn,
        kinematics::MinimizerMapBase const & min_map
    ) const;
    
    virtual
    bool
    defines_score_for_residue_pair(
        conformation::Residue const & rsd1,
        conformation::Residue const & rsd2,
        bool res_moving_wrt_eachother
    ) const;
    
    virtual
    etable::count_pair::CountPairFunctionCOP
    get_count_pair_function(
        Size const res1,
        Size const res2,
        pose::Pose const & pose,
        ScoreFunction const &
    ) const;
    
    virtual
    etable::count_pair::CountPairFunctionCOP
    get_count_pair_function(
        conformation::Residue const & rsd1,
        conformation::Residue const & rsd2
    ) const;
    
    virtual
    etable::count_pair::CountPairFunctionCOP
    get_intrares_countpair(
        conformation::Residue const & rsd1,
        pose::Pose const &,
        ScoreFunction const &
    ) const;
    
    virtual
    bool
    use_extended_residue_pair_energy_interface() const;
    
    virtual
    void
    setup_for_minimizing_for_residue_pair(
        conformation::Residue const & rsd1,
        conformation::Residue const & rsd2,
        pose::Pose const & pose,
        ScoreFunction const &,
        kinematics::MinimizerMapBase const &,
        ResSingleMinimizationData const &,
        ResSingleMinimizationData const &,
        ResPairMinimizationData & pair_data
    ) const;
    
    virtual
    void
    residue_pair_energy_ext(
        conformation::Residue const & rsd1,
        conformation::Residue const & rsd2,
        ResPairMinimizationData const & min_data,
        pose::Pose const & pose,
        ScoreFunction const & scorefxn,
        EnergyMap & emap
    ) const;

    
	void
	eval_atom_derivative_intra_RNA( //Called by eval_atom_derivative, specific case for RNA intra_res. Parin Sripakdeevong June 27, 2011.
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


	/// called during gradient-based minimization inside dfunc
	/**
		 F1 and F2 are not zeroed -- contributions from this atom are
		 just summed in
	**/
	virtual
	void
	eval_atom_derivative(
		id::AtomID const & id,
		pose::Pose const & pose,
		kinematics::DomainMap const & domain_map,
		ScoreFunction const & sfxn,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;


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
	bool
	defines_intrares_energy( EnergyMap const & weights ) const;


	virtual
	void
	eval_intrares_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		ScoreFunction const & ,
		EnergyMap & emap
	) const;


	virtual
	Distance
	atomic_interaction_cutoff() const;


	void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;


	virtual
	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap &// totals
	) const;

private:

	Vector
	get_base_vector( conformation::Residue const & rsd1, Size const i, pose::Pose const & pose ) const;

	void
	get_residue_energy_RNA_intra(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		Real & lk_polar_intra_RNA_score,
		Real & lk_nonpolar_intra_RNA_score,
		Real & lk_costheta_intra_RNA_score,
		bool const compute_polar,
		bool const compute_nonpolar
	) const;

	void
	get_residue_pair_energy_one_way(
		 conformation::Residue const & rsd1,
		 conformation::Residue const & rsd2,
		 pose::Pose const & pose,
		 Real & lk_polar_score,
		 Real & lk_nonpolar_score,
		 Real & lk_costheta_score,
		 bool const compute_polar,
		 bool const compute_nonpolar
	 ) const;

	Real
	eval_lk(
	conformation::Atom const & atom1,
	conformation::Atom const & atom2,
	Real & deriv,
    bool const & eval_deriv
    ) const;
    
    std::pair<Real,Real>
    eval_lk_efficient(
        conformation::Atom const & atom1,
        conformation::Atom const & atom2,
        Real & deriv,
        bool const & eval_deriv
    ) const;

	void
	distribute_pseudo_base_atom_derivatives( pose::Pose const & pose ) const;

	virtual
	core::Size version() const;


/////////////////////////////////////////////////////////////////////////////
// data
/////////////////////////////////////////////////////////////////////////////
private:

	etable::EtableEvaluatorOP etable_evaluator_;
	Real const safe_max_dis2_;
	Real const max_dis_;
	bool const verbose_;

};

}
}
}

#endif // INCLUDED_core_scoring_methods_LK_PolarNonPolarEnergy_HH
