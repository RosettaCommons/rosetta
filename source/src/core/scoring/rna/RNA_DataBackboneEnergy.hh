// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_DataBackboneEnergy.hh
/// @brief  Statistically derived rotamer pair potential class declaration
/// @author Phil Bradley
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_rna_RNA_DataBackboneEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_DataBackboneEnergy_hh

// Unit Headers
#include <core/scoring/rna/RNA_DataBackboneEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>

#include <core/scoring/func/Func.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
// AUTO-REMOVED #include <numeric/xyzVector.hh>

#include <utility/vector1.hh>



namespace core {
namespace scoring {
namespace rna {

///

typedef numeric::xyzVector< core::Real > Vector;

class RNA_DataBackboneEnergy : public methods::ContextIndependentTwoBodyEnergy  {
public:
	typedef methods::ContextIndependentTwoBodyEnergy  parent;
public:

	///
	RNA_DataBackboneEnergy();


	/// clone
	virtual
	methods::EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// scoring
	/////////////////////////////////////////////////////////////////////////////

	virtual
	void
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	setup_for_packing( pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

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
											 ScoreFunction const & scorefxn,
											 EnergyMap const & weights,
											 Vector & F1,
											 Vector & F2
											 ) const;

	virtual
	bool
	defines_intrares_energy( EnergyMap const & /*weights*/ ) const { return false; }

	//	virtual
	//	void
	//	finalize_total_energy(
	//		pose::Pose & pose,
	//		ScoreFunction const &,
	//		EnergyMap &// totals
	//	) const;

	virtual
 	Distance
	atomic_interaction_cutoff() const;

	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const {}

	void
	initialize_atom_numbers_sugar();

	//Vector
	//	get_mean_sugar_pos( core::conformation::Residue const & rsd ) const;

	Real
	get_sugar_env_score( core::conformation::Residue const & rsd_buried, core::conformation::Residue const & rsd_other ) const;

	bool
	check_sugar_atom( Size const & n ) const;

	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////

private:

	Real const dist_cutoff_;
	Real const dist_fade_;
	Real const well_depth_burial_;
	Real const well_depth_exposed_;
	utility::vector1< Size > atom_numbers_sugar_;
	utility::vector1< Size > atom_numbers_sugar_coarse_;
	constraints::FuncOP burial_function_;
virtual
core::Size version() const;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
