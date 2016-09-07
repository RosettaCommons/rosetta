// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file HelixNCapperMover.cc
/// @brief Mover to find and fix poorly N-capped helices.
/// @author ben bborgo@genetics.wustl.edu

#include <protocols/helix_capper/HelixNCapperMover.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cp.OptionKeys.gen.hh>
#include <basic/database/open.hh>

#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/methods/ContextIndependentLRTwoBodyEnergy.hh>
#include <core/scoring/methods/ContextDependentLRTwoBodyEnergy.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/PDBInfo.hh>

#include <protocols/simple_moves/ScoreMover.hh>

#include <utility/io/izstream.hh>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>

using namespace core;
using namespace core::conformation;
using namespace core::scoring;
using namespace core::scoring::hbonds;
using namespace utility::libsvm;

namespace protocols {
namespace helix_capper {

static THREAD_LOCAL basic::Tracer TR( "HelixNCapperMover" );

HelixNCapperMover::HelixNCapperMover() { setup_svm(); }
HelixNCapperMover::HelixNCapperMover(
	core::pose::Pose p
){
	start_pose_ = p;
	setup_svm();
}

HelixNCapperMover::~HelixNCapperMover()= default;

void HelixNCapperMover::dump_pdb_to_file( core::pose::Pose & posey, std::string filename ){
	posey.dump_pdb(filename);
}

void HelixNCapperMover::setup_svm() {

	string ncap_svm_filename = basic::database::full_name( "external/svm_models/helix_cap/Ncap_svm_model" );
	utility::io::izstream is( ncap_svm_filename );
	if ( !is.good() ) {
		utility_exit_with_message("Error: missing svm model for ncap classification in database!");
	}

	ncap_model_ = Svm_rosettaOP( new Svm_rosetta( ncap_svm_filename.c_str() ) );

	return;
}

void HelixNCapperMover::get_start_positions() {

	core::Size min_helix_length( 7 );

	helix_start_positions_.clear();
	core::scoring::dssp::Dssp all_dssp( start_pose_ );

	for ( core::Size i = 5 ; i <= start_pose_.total_residue()-5 ; ++i ) {

		if ( !start_pose_.residue(i).is_protein() ) continue;

		bool is_previous_helix(
			all_dssp.get_dssp_secstruct(i-1) == 'G' ||
			all_dssp.get_dssp_secstruct(i-1) == 'I' ||
			all_dssp.get_dssp_secstruct(i-1) == 'H'
		);

		bool is_this_helix(
			all_dssp.get_dssp_secstruct(i) == 'G' ||
			all_dssp.get_dssp_secstruct(i) == 'I' ||
			all_dssp.get_dssp_secstruct(i) == 'H'
		);

		if ( !is_previous_helix && is_this_helix ) {
			core::Size work_res( i+1 );
			core::Size helix_length( 1 );
			while ( work_res <= start_pose_.total_residue() &&
					( all_dssp.get_dssp_secstruct( work_res ) == 'G' ||
					all_dssp.get_dssp_secstruct( work_res ) == 'H' ||
					all_dssp.get_dssp_secstruct( work_res ) == 'I' ) ) {
				++work_res;
				++helix_length;
			}

			if ( helix_length >= min_helix_length ) {
				helix_start_positions_.push_back( i );
				//TR << "Found helix Ncap at " << i << std::endl;
			}
		}
	}
}

void HelixNCapperMover::get_Ncap_scores() {

	ncap_scores_.clear();

	// The features that are used are the one_body, two_body, hbond reporters
	core::pose::Pose pose_copy( start_pose_ );
	core::scoring::ScoreFunctionOP sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "talaris2013" );
	sfxn->setup_for_scoring( pose_copy );
	(*sfxn)( pose_copy );

	HBondSet hbond_set;
	hbond_set.setup_for_residue_pair_energies( pose_copy, false, false );

	// retrieve cached energies object
	Energies const & energies( pose_copy.energies() );
	assert(energies.energies_updated());
	// the neighbor/energy links
	EnergyGraph const & energy_graph( energies.energy_graph() );

	for ( core::Size ihelix = 1 ; ihelix <= helix_start_positions_.size() ; ++ihelix ) {
		core::Size icap( helix_start_positions_[ ihelix ] );

		// For each cap calculate the scores that are used as features in the Ncap
		// SVM classifier.

		utility::vector1< core::Real > features;
		core::scoring::EnergyMap emap;

		for ( core::Size resid = icap, end_resid = icap + 4 ; resid <= end_resid ; ++resid ) {

			core::conformation::Residue const & rsd( pose_copy.residue( resid ) );
			// Get the onebody energies
			emap.clear();
			sfxn->eval_ci_1b(rsd, pose_copy, emap);
			features.push_back( emap[ p_aa_pp ] );
			features.push_back( emap[ fa_dun ] );

			// Now the twobody energies
			// Need:  fa_atr, fa_rep, fa_sol, hack_elec_bb_bb, hack_elec_bb_sc, hack_elec_sc_sc
			// these are all treated as short-range ci2b energies

			core::Real fa_atr_sum( 0.0 );
			core::Real fa_rep_sum( 0.0 );
			core::Real fa_sol_sum( 0.0 );
			core::Real fa_elec_bb_bb_sum( 0.0 );
			core::Real fa_elec_bb_sc_sum( 0.0 );
			core::Real fa_elec_sc_sc_sum( 0.0 );

			for ( core::Size ires = 1 ; ires <= resid ; ++ires ) {
				for ( core::graph::Graph::EdgeListConstIter
						iru  = energy_graph.get_node(ires)->const_upper_edge_list_begin(),
						irue = energy_graph.get_node(ires)->const_upper_edge_list_end();
						iru != irue; ++iru ) {
					EnergyEdge const & edge( static_cast< EnergyEdge const &> (**iru) );
					Size const other_resid( edge.get_second_node_ind() );

					if ( (ires != resid) && (other_resid != resid) ) continue;

					emap.clear();

					//core::Size res1(0);
					//core::Size res2(0);
					if ( ires == resid ) {
						//res1 = resid;
						//res2 = other_resid;
						Residue const & other_rsd( pose_copy.residue(other_resid) );
						sfxn->eval_ci_2b(rsd, other_rsd, pose_copy, emap);
					} else {
						//res1 = ires;
						//res2 = resid;
						Residue const & other_rsd( pose_copy.residue(ires) );
						sfxn->eval_ci_2b(other_rsd, rsd, pose_copy, emap);
					}

					fa_atr_sum += emap[ fa_atr ];
					fa_rep_sum += emap[ fa_rep ];
					fa_sol_sum += emap[ fa_sol ];
					fa_elec_bb_bb_sum += emap[ fa_elec_bb_bb ];
					fa_elec_bb_sc_sum += emap[ fa_elec_bb_sc ];
					fa_elec_sc_sc_sum += emap[ fa_elec_sc_sc ];

					//TR << "Twobody fa_atr interaction between " << res1 << " and " << res2 << " is " << emap[ fa_atr ] << std::endl;
				}
			}

			features.push_back( fa_atr_sum );
			features.push_back( fa_rep_sum );
			features.push_back( fa_sol_sum );
			features.push_back( fa_elec_bb_bb_sum );
			features.push_back( fa_elec_bb_sc_sum );
			features.push_back( fa_elec_sc_sc_sum );

			// find backbone NH hbond for resid

			core::Real nh_hbond_sum( 0.0 );

			for ( Size i = 1; i<= hbond_set.nhbonds(); i++ ) {
				HBond const & hbond = hbond_set.hbond(i);
				if ( hbond.don_res() != resid ) continue;
				if ( !hbond.don_hatm_is_protein_backbone() ) continue;

				nh_hbond_sum += hbond.energy()/2;
			}

			features.push_back( nh_hbond_sum );
		}

		//   TR << "Features:  ";
		//   for( core::Size i = 1 ; i <= features.size() ; ++i ) {
		//    TR << features[i] << " ";
		//   }
		//   TR << std::endl;

		core::Real this_ncap_native_prob( ncap_prob_from_svm( features ) );
		ncap_scores_.push_back( this_ncap_native_prob );

		TR << "Helix Ncap at " << icap << " has prob score of " << this_ncap_native_prob << std::endl;
	}
}

core::Real HelixNCapperMover::ncap_prob_from_svm( utility::vector1< core::Real > & features ) {

	utility::vector1< core::Real > return_prob_vector;

	// Make the svm nodes
	utility::vector1< Svm_node_rosettaOP > ncap_feature_nodes;
	for ( Size i = 1 ; i <= features.size() ; ++i ) {
		Svm_node_rosettaOP new_node = Svm_node_rosettaOP( new Svm_node_rosetta( i, features[i] ) );
		ncap_feature_nodes.push_back( new_node );
	}

	return_prob_vector = ncap_model_->predict_probability( ncap_feature_nodes );

	//TR << "The size of the returned vector is " << return_prob_vector.size() << std::endl;

	//  for( Size i = 1 ; i <= return_prob_vector.size() ; ++i ) {
	//   TR << "Element " << i << " is " << return_prob_vector[i] << std::endl;
	//  }

	return return_prob_vector[1];
}


void HelixNCapperMover::apply(){
	get_start_positions();
	get_Ncap_scores();
}

void
HelixNCapperMover::set_excluded_positions() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	excluded_positions_.clear();

	if ( option[ cp::exclude_file ].active() ) {
		std::string exclude_file_name( option[ cp::exclude_file ] );

		// Probe for file
		std::ifstream exc_file( exclude_file_name.c_str() );

		if ( !exc_file ) {
			TR << "Exclude_file " << exclude_file_name << " not found." << std::endl;
			TR << "No positions will be excluded." << std::endl;
			return;
		}

		// Process one at a time

		core::Size exc_pos;
		char exc_chain;

		exc_file >> exc_pos;
		while ( !exc_file.eof() ) {
			exc_file >> exc_chain;
			TR << "Adding position " << exc_pos << " chain " << exc_chain << " to exclude list." << std::endl;
			excluded_positions_.push_back( start_pose_.pdb_info()->pdb2pose( exc_chain, exc_pos ) );
			exc_file >> exc_pos;
		}
		exc_file.close();
	}

	// Exclude any non-amino acid residues from mutation
	for ( core::Size i(1), ei( start_pose_.total_residue() ) ; i <= ei ; ++i ) {
		if ( !start_pose_.residue( i ).is_protein() &&
				( std::find(excluded_positions_.begin(), excluded_positions_.end(), i ) ==
				excluded_positions_.end() ) ) {
			excluded_positions_.push_back( i );
		}
	}

	if ( excluded_positions_.size() > 0 ) {
		TR << "Found " << excluded_positions_.size() << " positions excluded from mutation" << std::endl;
	} else {
		TR << "No positions will be excluded." << std::endl;
	}
	return;
}

}}
