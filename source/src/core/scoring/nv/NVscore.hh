// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/NV/NVscore.hh
/// @brief  Neighbor Vector algorithm for solvation approximation class declaration
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nv_NVscore_hh
#define INCLUDED_core_scoring_nv_NVscore_hh

//unit headers
#include <core/scoring/nv/NVscore.fwd.hh>

//project headers
#include <core/scoring/methods/ContextDependentOneBodyEnergy.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/nv/NVlookup.hh>

#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


//STL headers

namespace core {
namespace scoring {
namespace nv {

class NVscore : public methods::ContextDependentOneBodyEnergy  {
public:
	typedef methods::ContextDependentOneBodyEnergy  parent;

public:

	NVscore();

	/// clone
	virtual methods::EnergyMethodOP clone() const;

	virtual void setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const;

	virtual void setup_for_packing(pose::Pose & pose, utility::vector1< bool > const &, utility::vector1< bool > const & ) const;

	virtual void setup_for_derivatives( pose::Pose &pose, ScoreFunction const &  ) const;

	virtual void setup_for_minimizing(pose::Pose & pose, ScoreFunction const & ,kinematics::MinimizerMapBase const &) const;

	virtual void residue_energy(
			conformation::Residue const & rsd,
					pose::Pose const & pose,
					EnergyMap & emap
	) const;

	virtual void indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const;

	Real neighbor_weight( Vector::Value const & distance) const;

private:
		NVlookup const &lookup_table_;
		core::Real lower_bound_;
		core::Real lower_bound_squared_;
		core::Real upper_bound_;
		core::Real upper_bound_squared_;

virtual
core::Size version() const;
};

} //NV
} //scoring
} //core

#endif
