// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/SingleResidueMultifunc.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_core_optimization_SingleResidueMultifunc_hh
#define INCLUDED_core_optimization_SingleResidueMultifunc_hh

#include <core/optimization/SingleResidueMultifunc.fwd.hh>

#include <utility/graph/Graph.fwd.hh>
#include <core/optimization/types.hh>

#include <core/optimization/AtomTreeMultifunc.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace optimization {

/// @brief A streamlined AtomTreeMultifunc designed specifically for RTMIN.
///
/// @details Evaluates only the energies between the specified residue and the rest
/// of the Pose, assuming the nbr_atoms do not move (as in rotamer trials and repacking).
/// Could probably be sped up further with a customized dfunc().
/// DFPMIN seems to spend most of its time in func() rather than dfunc(),
/// so there's not as much to gain there anyway.
///
class SingleResidueMultifunc : public AtomTreeMultifunc
{
	typedef AtomTreeMultifunc parent;

public:

	SingleResidueMultifunc(
		pose::Pose & pose_in,
		Size const rsd_id_in,
		MinimizerMap & min_map_in,
		scoring::ScoreFunction const & scorefxn_in,
		utility::graph::GraphCOP packer_neighbor_graph_in,
		bool const deriv_check_in = false,
		bool const deriv_check_verbose_in = false
	);

	~SingleResidueMultifunc() override;

	// func

	Real
	operator ()( Multivec const & vars ) const override;

private:
	Size rsd_id_;
	utility::graph::GraphCOP packer_neighbor_graph_;

}; // SingleResidueMultifunc


} // namespace optimization
} // namespace core

#endif // INCLUDED_core_optimization_SingleResidueMultifunc_HH
