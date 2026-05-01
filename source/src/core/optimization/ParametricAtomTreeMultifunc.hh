// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/ParametricAtomTreeMultifunc.hh
/// @brief  Multifunc that co-minimizes parametric DOFs (Crick parameters) alongside standard
///         atom-tree DOFs (chi angles, backbone torsions, jumps).
/// @details Parametric DOFs are appended to the end of the standard DOF vector. The operator()
/// rebuilds parametric backbone coordinates via the Crick equations before applying standard DOFs
/// and scoring. The dfunc() computes standard derivatives via atom_tree_dfunc, then adds
/// parametric derivatives via the analytical Crick Jacobian chain rule.
/// @author Andy Watkins

#ifndef INCLUDED_core_optimization_ParametricAtomTreeMultifunc_hh
#define INCLUDED_core_optimization_ParametricAtomTreeMultifunc_hh

#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/MinimizerMap.fwd.hh>
#include <core/optimization/parametric_minimize_util.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>

namespace core {
namespace optimization {

class ParametricAtomTreeMultifunc : public Multifunc {

public:

	ParametricAtomTreeMultifunc(
		pose::Pose & pose_in,
		MinimizerMap & min_map_in,
		scoring::ScoreFunction const & scorefxn_in,
		utility::vector1< ParametricDOFInfo > const & parametric_dofs_in,
		bool deriv_check_in = false,
		bool deriv_check_verbose_in = false
	);

	~ParametricAtomTreeMultifunc() override;

	Real operator()( Multivec const & vars ) const override;

	void dfunc( Multivec const & vars, Multivec & dE_dvars ) const override;

	void dump( Multivec const & vars, Multivec const & vars2 ) const override;

	Size n_standard_dofs() const { return n_standard_dofs_; }
	Size n_parametric_dofs() const { return parametric_dofs_.size(); }
	Size total_dofs() const { return n_standard_dofs_ + parametric_dofs_.size(); }

private:

	void apply_parametric_dofs( Multivec const & vars ) const;

	pose::Pose & pose_;
	MinimizerMap & min_map_;
	scoring::ScoreFunction const & score_function_;

	utility::vector1< ParametricDOFInfo > parametric_dofs_;
	Size n_standard_dofs_;

	bool deriv_check_;
	bool deriv_check_verbose_;

}; // ParametricAtomTreeMultifunc

} // namespace optimization
} // namespace core

#endif
