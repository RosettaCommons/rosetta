// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/nmrspinlabel_util.hh
/// @brief   utility functions for working with NMRSpinlabel class used both in PCSEnergy and PREEnergy
/// @details last Modified: 09/24/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_nmrspinlabel_util_HH
#define INCLUDED_protocols_nmr_nmrspinlabel_util_HH

#include <utility/graph/Graph.hh>
#include <core/types.hh>
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace nmr {

typedef core::scoring::nmr::NMRSpinlabel::WeightCoordVector WeightCoordVector;

/// @brief energetic filter for spinlabel ensemble
WeightCoordVector
filter_spinlabel_ensemble_by_packerenergy(
	core::pose::Pose const & pose,
	core::scoring::nmr::NMRSpinlabel & spinlabel,
	core::Size const spinlabel_position
);

/// @brief computes the bump energy of a spinlabel conformer
core::PackerEnergy
bump_check(
	core::conformation::ResidueCOP rotamer,
	core::Size resid,
	core::scoring::ScoreFunction const & sf,
	core::pose::Pose const & pose,
	core::pack::task::PackerTask const & task,
	utility::graph::GraphCOP packer_neighbor_graph
);

} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_nmrspinlabel_util_HH
