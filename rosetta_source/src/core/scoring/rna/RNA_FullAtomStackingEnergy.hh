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


#ifndef INCLUDED_core_scoring_rna_RNA_FullAtomStackingEnergy_hh
#define INCLUDED_core_scoring_rna_RNA_FullAtomStackingEnergy_hh

// Unit Headers
#include <core/scoring/rna/RNA_FullAtomStackingEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_LowResolutionPotential.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
// AUTO-REMOVED #include <numeric/xyzMatrix.hh>

#include <platform/types.hh>
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/RotamerSetBase.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/scoring/DerivVectorPair.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/MinimizationData.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/ContextIndependentTwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethod.fwd.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/methods/EnergyMethodCreator.fwd.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/TwoBodyEnergy.fwd.hh>
#include <core/scoring/methods/TwoBodyEnergy.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <string>
#include <vector>



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
	setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual
	void
	setup_for_derivatives( pose::Pose & pose, ScoreFunction const & ) const;

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
															pose::Pose const & pose
															) const;

	bool
	check_base_base_OK(
	   conformation::Residue const & rsd1,
		 conformation::Residue const & rsd2,
		 Size const & m, Size const & n ) const;

  Real const prefactor_;
  Distance const full_stack_cutoff_;
  Distance const dist_cutoff_;
  Real const dist_cutoff2_;
	bool const base_base_only_;
virtual
core::Size version() const;

};


}
}
}

#endif // INCLUDED_core_scoring_ScoreFunction_HH
