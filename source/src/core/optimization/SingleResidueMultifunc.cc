// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/optimization/SingleResidueMultifunc.cc
///
/// @brief
/// @author Ian W. Davis

#include <core/optimization/SingleResidueMultifunc.hh>

#include <utility/graph/Graph.hh>
#include <core/optimization/types.hh>
#include <core/optimization/AtomTreeMultifunc.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <basic/prof.hh>

#include <core/optimization/NumericalDerivCheckResult.fwd.hh>

#include <utility>
#include <utility/vector1.hh>


namespace core {
namespace optimization {

SingleResidueMultifunc::SingleResidueMultifunc(
	pose::Pose & pose_in,
	Size const rsd_id_in,
	MinimizerMap & min_map_in,
	scoring::ScoreFunction const & scorefxn_in,
	utility::graph::GraphCOP packer_neighbor_graph_in,
	bool const deriv_check_in,
	bool const deriv_check_verbose_in
):
	AtomTreeMultifunc(pose_in, min_map_in, scorefxn_in, deriv_check_in, deriv_check_verbose_in),
	rsd_id_(rsd_id_in),
	packer_neighbor_graph_(std::move(packer_neighbor_graph_in))
{}

SingleResidueMultifunc::~SingleResidueMultifunc() = default;

// func
Real
SingleResidueMultifunc::operator ()( Multivec const & vars ) const {
	using namespace conformation;
	using namespace scoring;

	pose::Pose & pose = parent::pose();
	MinimizerMap const & min_map = parent::min_map();
	scoring::ScoreFunction const & score_function = parent::score_function();

	PROF_START( basic::FUNC );
	min_map.copy_dofs_to_pose( pose, vars );
	conformation::Residue const & rsd = pose.residue(rsd_id_);

	// Code adapted from RotamerSet_::compute_onebody_energies()
	EnergyMap emap;
	EnergyMap emap2b;

	pose.scoring_begin( score_function );
	// Setup (particularly hbonds) is expensive and so is done once, in RTMIN.

	score_function.eval_ci_1b( rsd, pose, emap );
	score_function.eval_cd_1b( rsd, pose, emap );
	score_function.eval_intrares_energy( rsd, pose, emap );

	for ( utility::graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph_->get_node( rsd_id_ )->const_edge_list_begin(),
			ire = packer_neighbor_graph_->get_node( rsd_id_ )->const_edge_list_end();
			ir != ire; ++ir ) {
		int const neighbor_id( (*ir)->get_other_ind( rsd_id_ ) );
		//if ( task.pack_residue( neighbor_id ) ) continue;
		Residue const & neighbor( pose.residue( neighbor_id ) );

		emap2b.zero();
		score_function.eval_ci_2b( rsd, neighbor, pose, emap2b );
		emap += emap2b;

		emap2b.zero();
		score_function.eval_cd_2b( rsd, neighbor, pose, emap2b );
		emap += emap2b;
	}

	// long-range energy interactions
	// Iterate across the long range energy functions and use the iterators generated
	// by the LRnergy container object
	for ( auto
			lr_iter = score_function.long_range_energies_begin(),
			lr_end  = score_function.long_range_energies_end();
			lr_iter != lr_end; ++lr_iter ) {
		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*lr_iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue; // only score non-emtpy energies.

		// Potentially O(N) operation leading to O(N^2) behavior
		for ( ResidueNeighborConstIteratorOP
				rni = lrec->const_neighbor_iterator_begin( rsd_id_ ),
				rniend = lrec->const_neighbor_iterator_end( rsd_id_ );
				(*rni) != (*rniend); ++(*rni) ) {
			Size const neighbor_id = rni->neighbor_id();
			debug_assert( neighbor_id != rsd_id_ );

			(*lr_iter)->residue_pair_energy( rsd, pose.residue( neighbor_id ), pose, score_function, emap );

		} // (potentially) long-range neighbors of rsd
	} // long-range energy functions

	// give energyfunctions a chance update/finalize energies
	// etable nblist calculation is performed here (?)
	for ( auto it=score_function.all_energies_begin(),
			it_end = score_function.all_energies_end(); it != it_end; ++it ) {
		(*it)->finalize_total_energy( pose, score_function, emap );
	}

	//emap.show_weighted(std::cout, score_function.weights());

	Real const score = score_function.weights().dot( emap );
	pose.scoring_end( score_function );
	PROF_STOP( basic::FUNC );
	return score;
}

} // namespace optimization
} // namespace core

