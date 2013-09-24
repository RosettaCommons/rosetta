// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/GeometricSolEnergy.hh
/// @brief  Geometric solvation energy (context dependent)
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_geometric_solvation_GeometricSolEnergy_hh
#define INCLUDED_core_scoring_geometric_solvation_GeometricSolEnergy_hh

// Unit Headers
#include <core/scoring/geometric_solvation/GeometricSolEnergy.fwd.hh>
#include <core/scoring/geometric_solvation/GeometricSolEnergyEvaluator.fwd.hh>

#include <core/types.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/methods/ContextDependentTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/scoring/NeighborList.tmpl.hh>
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/MinimizationData.hh>
#include <core/kinematics/MinimizerMapBase.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
#include <core/scoring/etable/count_pair/CountPairNone.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/types.hh>


namespace core {
namespace scoring {
namespace geometric_solvation {

///
class GeometricSolEnergy : public methods::ContextDependentTwoBodyEnergy  {
public:
	typedef methods::ContextDependentTwoBodyEnergy  parent;
public:

	///
	GeometricSolEnergy( methods::EnergyMethodOptions const & options );

	///@brief copy c-tor
	GeometricSolEnergy( GeometricSolEnergy const & src );

	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	///
	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	///
	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;
    
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
        conformation::Residue const &,
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
        ScoreFunction const &,
        EnergyMap & emap
    ) const;
    
    virtual
    void
    finalize_total_energy(
        pose::Pose & pose,
        ScoreFunction const &,
        EnergyMap & totals
    ) const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	/// This evaluates everything for now,
	/// but eventually may want to split this
	/// based on backbone/backbone vs. others,
	/// as is carried out in HBondEnergy.cc
	virtual
	void
	residue_pair_energy(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		ScoreFunction const &,
		EnergyMap & emap
	) const;

	// Undefined, commenting out to fix PyRosetta build
	/* void
	eval_atom_derivative_intra_RNA(
		 id::AtomID const & atom_id,
		 pose::Pose const & pose,
		 EnergyMap const & weights,
		 Vector & F1,
		 Vector & F2
	) const;
	*/

	/// f1 and f2 are zeroed
	virtual
	void
	eval_atom_derivative(
		id::AtomID const & atom_id,
		pose::Pose const & pose,
		kinematics::DomainMap const &,
		ScoreFunction const &,
		EnergyMap const & weights,
		Vector & F1,
		Vector & F2
	) const;

	// Undefined, commenting out to fix PyRosetta build
	/* Real
	eval_atom_energy(
		id::AtomID const & atom_id,
		pose::Pose const & pose
		) const; */

	///
	//	virtual
	//	void
	//	finalize_total_energy(
	//		pose::Pose & pose,
	//		ScoreFunction const &,
	//		EnergyMap & totals
	//	) const;

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

	//Real
	//hydrogen_interaction_cutoff2() const;

	///@brief GeometricSolEnergy is context sensitive
	virtual
	void indicate_required_context_graphs(
		utility::vector1< bool > & context_graphs_required ) const;

private:

	virtual
	core::Size version() const;

private:

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
	methods::EnergyMethodOptionsOP options_;

	GeometricSolEnergyEvaluatorOP evaluator_;

};

} // hbonds
} // scoring
} // core

#endif
