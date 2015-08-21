// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rtmin.hh
/// @brief  rotamer trials with minimization module header
/// @author Ian W. Davis (ian.w.davis@gmail.com)

#ifndef INCLUDED_core_pack_rtmin_hh
#define INCLUDED_core_pack_rtmin_hh

// pack headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/scmin/SCMinMinimizerMap.fwd.hh>

// conformation headers
#include <core/conformation/Residue.fwd.hh>

// pose headers
#include <core/pose/Pose.fwd.hh>

// scoring headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>

// utility headers
#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {

class RTMin {

public:

	RTMin();
	RTMin(
		bool minimize_ligand_chis,
		bool minimize_ligand_jumps
	);

	~RTMin();

public:
	void set_nonideal(bool nonideal_in) { nonideal_ = nonideal_in; }
	void set_cartesian(bool cartesian_in) { cartesian_ = cartesian_in; }

	void
	rtmin(
		pose::Pose & pose,
		scoring::ScoreFunction const & sfxn,
		task::PackerTaskOP input_task
	) const;

private:
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool minimize_ligand_chis_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool minimize_ligand_jumps_;
	bool nonideal_;
	bool cartesian_;
};

void
reinitialize_mingraph_neighborhood_for_residue(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn,
	utility::vector1< conformation::ResidueCOP > const & bgres,
	pack::scmin::SCMinMinimizerMap const & scminmap,
	conformation::Residue const & rsd,
	scoring::MinimizationGraph & mingraph
);


}
}

#endif
