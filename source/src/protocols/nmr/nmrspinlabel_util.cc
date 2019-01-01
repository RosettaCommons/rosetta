// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/nmrspinlabel_util.cc
/// @brief   utility functions for working with NMRSpinlabel class used both in PCSEnergy and PREEnergy
/// @details last Modified: 09/24/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <protocols/nmr/nmrspinlabel_util.hh>

// Project headers
#include <core/scoring/nmr/NMRSpinlabel.hh>
#include <core/scoring/nmr/NMRDummySpinlabelEnsemble.hh>

// Package headers
#include <core/types.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pack/make_symmetric_task.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/BumpSelector.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/graph/Graph.hh>

// Numeric headers
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

// C++ headers
#include <string>
#include <iostream>
#include <algorithm>

namespace protocols {
namespace nmr {

static basic::Tracer TR( "protocols.nmr.nmrspinlabel_util" );

/// @brief energetic filter for spinlabel ensemble
WeightCoordVector
filter_spinlabel_ensemble_by_packerenergy(
	core::pose::Pose const & pose,
	core::scoring::nmr::NMRSpinlabel & spinlabel,
	core::Size const spinlabel_position
)
{
	using namespace core::scoring;
	using namespace core::scoring::nmr;
	using namespace core::pack;
	using namespace core::pack::task;
	using namespace core::pack::rotamer_set;
	using NMRDummySpinlabelConformerTableIter = core::scoring::nmr::NMRDummySpinlabelEnsemble::NMRDummySpinlabelConformerTableIter;

	// Setup some default scorefunction
	std::string default_weight_set = basic::options::option[ basic::options::OptionKeys::score::weights ].default_value();
	ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( default_weight_set );

	// Setup packer task
	PackerTaskOP packertask = TaskFactory::create_packer_task( pose );
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		PackerTaskOP symm_packertask = make_new_symmetric_PackerTask_by_requested_method( pose, packertask );
		packertask = symm_packertask;
	}
	utility::vector1<bool> residues_to_pack(pose.total_residue(), false);
	residues_to_pack[ spinlabel_position ] = true;
	packertask->restrict_to_residues(residues_to_pack);
	packertask->restrict_to_repacking();
	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, *sfxn, packertask );

	// Fill rotamer vector with NMRDummySpinlabelConformers
	// First, define pairs of atoms in spin-label and target residue
	// which should be aligned
	utility::vector1< ResidueOP > spinlabel_conf_vector;
	spinlabel_conf_vector.reserve(spinlabel.get_dummy_ensemble()->get_ensemble_size());
	utility::vector1< std::pair< std::string, std::string > > atm2aln;
	atm2aln.push_back(std::make_pair("CA", "CA"));
	if ( pose.residue(spinlabel_position).has("CB") ) {
		atm2aln.push_back(std::make_pair("CB","CB"));
	} else if ( pose.residue(spinlabel_position).has("2HA") ) {
		atm2aln.push_back(std::make_pair("2HA","CB"));
	} else {
		utility_exit_with_message("ERROR when trying to orient NMRDummySpinlabel onto protein. Target residue has no CB or 2HA atom.");
	}
	atm2aln.push_back(std::make_pair("C","C"));
	// Second, align every NMRDummySpinlabelConformer and put it in vector
	for ( NMRDummySpinlabelConformerTableIter ndslct_iter = spinlabel.get_dummy_ensemble()->get_conformer_table().begin(),
			ndslct_end = spinlabel.get_dummy_ensemble()->get_conformer_table().end(); ndslct_iter != ndslct_end; ndslct_iter++ ) {
		ResidueOP rsd = (*ndslct_iter)->get_residue();
		rsd->orient_onto_residue(pose.residue(spinlabel_position), atm2aln);
		spinlabel_conf_vector.push_back(rsd);
	}

	// Perform bump check of every NMRDummySpinlabelConformer
	BumpSelector bump_selector( packertask->max_rotbump_energy() );
	utility::vector1<bool> passing(spinlabel_conf_vector.size(),false);
	utility::vector1<core::Real> conf_scores(spinlabel_conf_vector.size());
	TR.Debug << "Filtering NMRDummySpinlabel conformers of residue " << spinlabel.get_code() << " on bump energy for site at " << spinlabel_position << std::endl;
	for ( core::Size i(1); i <= spinlabel_conf_vector.size(); ++i ) {
		ResidueOP conf = spinlabel_conf_vector[ i ];
		TR.Debug << " Suggested conformer " << i << ": ";
		for ( core::Size j(1); j <= conf->nchi(); ++j ) TR.Debug << conf->chi()[j] << ' ';
		core::PackerEnergy bumpenergy = bump_check( conf, spinlabel_position, *sfxn, pose, *packertask, packer_neighbor_graph );
		conf_scores[i] = static_cast<core::Real>(bumpenergy);
		TR.Debug << " Bump energy: " << bumpenergy;
		BumpSelectorDecision decision =  bump_selector.iterate_bump_selector( bumpenergy );
		switch ( decision ) {
		case KEEP_ROTAMER :
			TR.Debug << " ... added" << std::endl;
			passing[i] = true;
			break;
		case DELETE_PREVIOUS_ROTAMER :
			TR.Debug << " ... replace previous" << std::endl;
			if ( passing.size() == 0 ) {
				utility_exit_with_message("Internal consistency error: cannot replace non-existant previous residue.");
			}
			passing[i-1] = false;
			passing[i] = true;
			break;
		case DELETE_ROTAMER : // Do nothing.
			TR.Debug << " ... deleted" << std::endl;
			break;
		}
	}

	// Emergency conformer
	if ( std::count(passing.begin(), passing.end(), true ) == 0 ) {
		TR.Debug << "No NMRDummySpinlabel conformer passed bump energy filter. Picking random emergency conformer." << std::endl;
		passing[numeric::random::random_range(1, spinlabel_conf_vector.size())] = true;
	}

	WeightCoordVector sl_wghts_coords = spinlabel.filter_spinlabel_ensemble_by_mask(pose,spinlabel_position,passing,conf_scores);
	return sl_wghts_coords;
}

/// @brief computes the bump energy of a spinlabel conformer
core::PackerEnergy
bump_check(
	core::conformation::ResidueCOP rotamer,
	core::Size resid,
	core::scoring::ScoreFunction const & sf,
	core::pose::Pose const & pose,
	core::pack::task::PackerTask const & task,
	utility::graph::GraphCOP packer_neighbor_graph
)
{
	using namespace core::scoring;
	using namespace core::conformation;

	EnergyMap emap;

	for ( utility::graph::Graph::EdgeListConstIter
			ir  = packer_neighbor_graph->get_node( resid )->const_edge_list_begin(),
			ire = packer_neighbor_graph->get_node( resid )->const_edge_list_end();
			ir != ire; ++ir ) {
		int const neighbor_id( (*ir)->get_other_ind( resid ) );
		Residue const & neighbor( pose.residue( neighbor_id ) );

		if ( ! task.pack_residue( neighbor_id ) ) {
			sf.bump_check_full( *rotamer, neighbor, pose, emap);
		} else {
			sf.bump_check_backbone( *rotamer, neighbor, pose, emap);
		}
	}
	return static_cast< core::PackerEnergy > (sf.weights().dot( emap ));
}

} // namespace nmr
} // namespace protocols
