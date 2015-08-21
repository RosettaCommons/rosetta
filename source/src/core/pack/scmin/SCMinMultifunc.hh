// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/scmin/SCMinMultifunc.hh
/// @brief  Class for interfacing with the minimizer during sidechain minimization.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_scmin_SCMinMultifunc_hh
#define INCLUDED_core_pack_scmin_SCMinMultifunc_hh

// Unit headers
#include <core/pack/scmin/SCMinMultifunc.fwd.hh>

// Package Headers
#include <core/pack/scmin/SCMinMinimizerMap.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/optimization/Multifunc.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/scoring/ScoreType.hh>


namespace core {
namespace pack {
namespace scmin {

/// @brief
class SCMinMultifunc : public optimization::Multifunc
{
public:
	typedef optimization::Multivec Multivec;
public:
	/// @brief Constructor.  The pose is only modified during setup-for-scoring calls.  Its
	/// residue objects are untouched (though they are accessed).  The SCMinMultifunc modifies
	/// the residues that are contained in the SCMinMinimizerMap's AtomTreeCollection.  The
	/// minmap also describes which degrees of freedom are free.
	/// The minimization graph should have already been setup so that its edges connect all
	/// neighboring residues which will be involved in the minimization:
	/// any residue that's being minimized must have all of its neighbors represented by
	/// edges in the graph.  Active edges and nodes must already have had
	/// "setup_for_minimizing_for_{residue/residue_pair}" invoked: basically,
	/// the SCMinMultifunc is absolved of all responsibility for setting up.
	SCMinMultifunc(
		pose::Pose & p,
		utility::vector1< conformation::ResidueCOP > const & bg_residues,
		scoring::ScoreFunction const & sfxn,
		scoring::MinimizationGraph & mingraph,
		SCMinMinimizerMap & scminmap
	);

	virtual ~SCMinMultifunc();

	virtual
	Real
	operator ()( Multivec const & chi ) const;


	virtual
	void
	dfunc( Multivec const & chi, Multivec & dE_dchi ) const;

	virtual
	bool
	abort_min( Multivec const & ) const;

	using optimization::Multifunc::dump;

	/// @brief Error state reached; dump out something corresponding to the
	/// var assignment.  Default base class implementation: no_op();
	virtual
	void
	dump( Multivec const & /*vars*/, Multivec const & /*vars2*/ ) const;

	void
	scmin_numerical_derivative_check( Multivec const & start_vars, Multivec & dE_dvars ) const;

private:

	/*void
	eval_atom_deriv(
	id::AtomID const & atom,
	Vector & F1,
	Vector & F2
	) const;*/


private:
	pose::Pose & pose_;
	utility::vector1< conformation::ResidueCOP > const & bg_residues_;
	scoring::ScoreFunction const & sfxn_;
	scoring::MinimizationGraph & g_;
	SCMinMinimizerMap & scminmap_;

	//mutable scoring::EnergyMap emap_;
	scoring::ScoreTypes scoretypes_;
};


} // namespace scmin
} // namespace pack
} // namespace core

#endif
