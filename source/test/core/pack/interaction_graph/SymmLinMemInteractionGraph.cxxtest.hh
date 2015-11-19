// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/SymmLinMemInteractionGraph.cxxtest.hh
/// @brief  test suite for the symmetric, minimalist on-the-fly interaction graph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/interaction_graph/SymmLinMemInteractionGraph.hh>

#include <core/chemical/AA.hh>
#include <core/graph/Graph.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/make_symmetric_task.hh>

#include <core/io/silent/SilentFileData.hh>

// Test headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <utility/vector1.hh>


class SymmLinMemInteractionGraphTests : public CxxTest::TestSuite {
public:

	void setUp() {
		core_init_with_additional_options( "-no_optH -symmetry:symmetry_definition core/scoring/symmetry/sym_def.dat" );
	}

	core::PackerEnergy
	get_bb_bbE(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2
	)
	{
		core::scoring::EnergyMap emap;
		sfxn.eval_ci_2b_bb_bb( res1, res2, pose, emap );
		sfxn.eval_cd_2b_bb_bb( res1, res2, pose, emap );
		return static_cast< core::PackerEnergy > ( sfxn.weights().dot( emap ) );
	}

	core::PackerEnergy
	get_sc_bbE(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2
	)
	{
		core::scoring::EnergyMap emap;
		sfxn.eval_ci_2b_bb_sc( res2, res1, pose, emap );
		sfxn.eval_cd_2b_bb_sc( res2, res1, pose, emap );
		return static_cast< core::PackerEnergy > ( sfxn.weights().dot( emap ) );
	}

	core::PackerEnergy
	get_sc_scE(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2
	)
	{
		core::scoring::EnergyMap emap;
		sfxn.eval_ci_2b_sc_sc( res2, res1, pose, emap );
		sfxn.eval_cd_2b_sc_sc( res2, res1, pose, emap );
		return static_cast< core::PackerEnergy > ( sfxn.weights().dot( emap ) );
	}

	core::PackerEnergy
	get_residue_pairE(
		core::pose::Pose const & pose,
		core::scoring::ScoreFunction const & sfxn,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2
	)
	{
		core::scoring::EnergyMap emap;
		sfxn.eval_ci_2b( res2, res1, pose, emap );
		sfxn.eval_cd_2b( res2, res1, pose, emap );
		return static_cast< core::PackerEnergy > ( sfxn.weights().dot( emap ) );
	}

	void verify_that_rotamer_coordinates_are_right(
		core::pack::interaction_graph::SymmLinearMemoryInteractionGraphOP symlinmem_ig,
		core::conformation::symmetry::SymmetryInfoCOP symm_info,
		core::pose::Pose & pose,
		core::pack::task::PackerTaskOP /*task*/,
		core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets

	)
	{
		using namespace core::pack::rotamer_set;
		using namespace core::pack::interaction_graph;

		Size const nres_asu = symm_info->num_independent_residues();

		for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			SymmLinearMemNode * ii_node = static_cast< SymmLinearMemNode * > ( symlinmem_ig->get_node( ii ) );
			RotamerSetCOP ii_rotset = rotsets->rotamer_set_for_moltenresidue( ii );
			Size ii_residue1 =  ( rotsets->moltenres_2_resid( ii ) - 1 ) % nres_asu + 1;

			// std::cout << "Residue1 " << ii_residue1 << " ii " << ii << " iiresid " << rotsets->moltenres_2_resid(ii) << " nrots " << ii_rotset->num_rotamers() << std::endl;

			for ( Size jj = 1; jj <= ii_rotset->num_rotamers(); ++jj ) {
				//std::cout << "Rotamer " << jj << std::endl;
				pose.replace_residue( ii_rotset->resid(), *ii_rotset->rotamer(jj), false );
				for ( Size kk = 1; kk <= symm_info->subunits(); ++kk ) {
					Size kkresid = (kk-1)*nres_asu + ii_residue1;
					// std::cout << "Subunit " << kk << " kkresid " << kkresid << std::endl;
					// std::cout << "transform: " << symlinmem_ig->symmetric_transform( kk ) << std::endl;
					core::conformation::Residue const & ig_res( ii_node->get_rotamer( jj, kk ));
					core::conformation::Residue const & pose_res( pose.residue( kkresid ));
					// std::cout << "resids -- pose: " << pose_res.seqpos() << " vs ig: " << ig_res.seqpos() << std::endl;
					TS_ASSERT( & ig_res.type() == & pose_res.type() );
					if ( & ig_res.type() != & pose_res.type() ) continue;
					for ( Size ll = 1; ll <= 1 /*ig_res.natoms()*/; ++ll ) {
						TS_ASSERT_DELTA( ig_res.xyz(ll).x(), pose_res.xyz(ll).x(), 5e-5 );
						TS_ASSERT_DELTA( ig_res.xyz(ll).y(), pose_res.xyz(ll).y(), 5e-5 );
						TS_ASSERT_DELTA( ig_res.xyz(ll).z(), pose_res.xyz(ll).z(), 5e-5 );
					}
				}
			}
		}
	}

	void verify_proper_edges_exist(
		core::pack::interaction_graph::SymmLinearMemoryInteractionGraphOP symlinmem_ig,
		core::conformation::symmetry::SymmetryInfoCOP symm_info,
		core::pose::Pose const & pose,
		core::scoring::ScoreFunctionOP sfxn,
		core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets
	)
	{
		using namespace core;
		using namespace core::conformation;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::interaction_graph;

		Size const nres_asu = symm_info->num_independent_residues();
		Size const nsubunits = symm_info->subunits();
		Real const sfxn_reach = sfxn->max_atomic_interaction_cutoff();

		utility::vector1< Real > max_radius( rotsets->nmoltenres(), 0.0 );
		for ( core::Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			RotamerSetCOP iirotset = rotsets->rotamer_set_for_moltenresidue( ii );
			for ( core::Size jj = 1; jj <= iirotset->get_n_residue_groups(); ++jj ) {
				Real jj_radius = iirotset->rotamer( iirotset->get_residue_group_begin(jj) )->nbr_radius();
				if ( max_radius[ ii ] < jj_radius ) { max_radius[ ii ] = jj_radius; }
			}
		}

		for ( core::Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			Size ii_resid  = (rotsets->moltenres_2_resid(ii)-1) % nres_asu + 1;
			Size ii_resid1 = (ii_resid-1) % nres_asu + 1;
			for ( core::Size jj = ii+1; jj <= rotsets->nmoltenres(); ++jj ) {
				Size jj_resid = rotsets->moltenres_2_resid(jj);
				Size jj_resid1 = (jj_resid-1) % nres_asu + 1;

				Real contact_thresh = max_radius[ii] + max_radius[jj] + sfxn_reach;
				contact_thresh = contact_thresh * contact_thresh;
				bool in_contact = false;
				for ( core::Size kk = 1; kk <= nsubunits; ++kk ) {
					Size iiresid_on_kk = (kk-1)*nres_asu + ii_resid1;
					Size jjresid_on_kk = (kk-1)*nres_asu + jj_resid1;
					{
						Residue const & iires = pose.residue( ii_resid );
						Residue const & jjres = pose.residue( jjresid_on_kk );
						if ( iires.xyz( iires.nbr_atom() ).distance_squared( jjres.xyz( jjres.nbr_atom() )) <= contact_thresh ) {
							in_contact = true;
							break;
						}
					}
					{
						Residue const & iires = pose.residue( iiresid_on_kk );
						Residue const & jjres = pose.residue( jj_resid );
						if ( iires.xyz( iires.nbr_atom() ).distance_squared( jjres.xyz( jjres.nbr_atom() )) <= contact_thresh ) {
							in_contact = true;
							break;
						}
					}
				}

				TS_ASSERT_EQUALS( symlinmem_ig->get_edge_exists( ii, jj ), in_contact );
			}
		}
	}

	void verify_edge_adjacency_matrices_are_accurate(
		core::pack::interaction_graph::SymmLinearMemoryInteractionGraphOP symlinmem_ig,
		core::conformation::symmetry::SymmetryInfoCOP symm_info,
		core::scoring::ScoreFunctionOP sfxn,
		core::pack::task::PackerTaskOP /*task*/,
		core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets

	) {
		using namespace core;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::interaction_graph;

		Size const nres_asu = symm_info->num_independent_residues();
		Size const nsubunits = symm_info->subunits();
		Size const asymmetric_unit = symlinmem_ig->asymmetric_unit();
		Real const sfxnreach = sfxn->max_atomic_interaction_cutoff();
		for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			int const iiresid = rotsets->moltenres_2_resid(ii);
			int const iiresid1 = (iiresid-1)%nres_asu + 1;
			RotamerSetCOP iirotset = rotsets->rotamer_set_for_residue( iiresid );
			Size const iinresgroups = iirotset->get_n_residue_groups();
			SymmLinearMemNode * iinode = static_cast< SymmLinearMemNode * > ( symlinmem_ig->get_node( ii ) );

			for ( int jj = 1; jj <= iinode->get_num_incident_edges(); ++jj ) {
				SymmOnTheFlyEdge * iijjedge = iinode->get_incident_otf_edge(jj);

				//Size which_node = 1;
				//if ( iijjedge->get_other_ind( ii ) < int( ii ) ) which_node = 2;

				if ( iijjedge->get_other_ind( ii ) < (int) ii ) continue; // only look at upper edges

				SymmOnTheFlyNode * jjnode = iinode->get_adjacent_otf_node(jj);
				Size const jjmoltenres = iijjedge->get_second_node_ind();
				Size const jjresid = rotsets->moltenres_2_resid(jjmoltenres);
				Size const jjresid1 = (jjresid-1)%nres_asu + 1;
				//std::cout << "examining iimoltenres " << ii << " jjmoltenres " << jjmoltenres << std::endl;
				Size const jjnresgroups = rotsets->rotamer_set_for_moltenresidue( jjmoltenres )->get_n_residue_groups();

				// 1. think of ii as having a rotamer coming from the asymmetric unit, iterate across all subunits for
				// jjmoltenres and the check every pair of residue types
				for ( Size kk = 1; kk <= nsubunits; ++kk ) {
					core::Size kkresid = jjresid1 + (kk-1)*nres_asu;
					core::Size kk_score_multiply = symm_info->score_multiply( iiresid, kkresid );
					//std::cout << " 1 kk= " << kk << " score multiply: " << kk_score_multiply << std::endl;
					for ( Size ll = 1; ll <= iinresgroups; ++ll ) {
						core::conformation::Residue llrot = iinode->get_asu_rotamer( iinode->get_state_offset_for_restype_group( ll ) + 1 );
						for ( Size mm = 1; mm <= jjnresgroups; ++mm ) {
							core::conformation::Residue mmrot = jjnode->get_rotamer( jjnode->get_state_offset_for_restype_group( mm ) + 1, kk );
							Real llmm_cutoff = llrot.nbr_radius() + mmrot.nbr_radius() + sfxnreach;
							Real d2cutoff = llmm_cutoff * llmm_cutoff;
							Real d2 = llrot.xyz( llrot.nbr_atom() ).distance_squared( mmrot.xyz( mmrot.nbr_atom() ) );
							Real expected_val = d2 < d2cutoff ? kk_score_multiply : 0;
							// Real actual_val = iijjedge->residues_adjacent_for_subunit_pair( which_node, kk, which_node == 1 ? ll : mm , which_node == 1 ? mm : ll );
							Real actual_val = iijjedge->residues_adjacent_for_subunit_pair( 1, kk, ll, mm );
							TS_ASSERT_EQUALS( actual_val, expected_val );
							if ( actual_val !=  expected_val ) {
								std::cout << "ii: " << ii << " jj: " << jjmoltenres << " kk " << kk << " ll " << ll << " mm " << mm << " d2: " << d2 << " cutoff: " << d2cutoff << std::endl;
							}
						}
					}
				}

				// 2. think of jjmoltenres as having a rotamer coming from the asymmetric unit, iterate across all subunits for
				// iimoltenres and the check every pair of residue types
				for ( Size kk = 1; kk <= nsubunits; ++kk ) {
					core::Size kkresid = iiresid1 + (kk-1)*nres_asu;
					core::Size kk_score_multiply = symm_info->score_multiply( jjresid, kkresid );
					if ( kk == asymmetric_unit ) {
						// the interactions within the asymmetric unit should only be counted once; and by convention
						// they are counted with the lower node is the one considered to originate from the asymmetric unit
						kk_score_multiply = 0;
					}
					//std::cout << " 2 kk= " << kk << " score multiply: " << kk_score_multiply << std::endl;
					for ( Size ll = 1; ll <= jjnresgroups; ++ll ) {
						core::conformation::Residue llrot = jjnode->get_asu_rotamer( jjnode->get_state_offset_for_restype_group( ll ) + 1 );
						for ( Size mm = 1; mm <= iinresgroups; ++mm ) {
							core::conformation::Residue mmrot = iinode->get_rotamer( iinode->get_state_offset_for_restype_group( mm ) + 1, kk );
							Real llmm_cutoff = llrot.nbr_radius() + mmrot.nbr_radius() + sfxnreach;
							Real d2cutoff = llmm_cutoff * llmm_cutoff;
							Real d2 = llrot.xyz( llrot.nbr_atom() ).distance_squared( mmrot.xyz( mmrot.nbr_atom() ) );
							// NOTE: do not count asu/asu interactions twice, so if kk==1, then all interactions should be marked with a 0
							Real expected_val = d2 < d2cutoff ? kk_score_multiply : 0;
							Real actual_val = iijjedge->residues_adjacent_for_subunit_pair( 2, kk, ll, mm );
							TS_ASSERT_EQUALS( actual_val, expected_val );
							if ( actual_val !=  expected_val ) {
								std::cout << "iimoltenres: " << ii << " jjmoltenres: " << jjmoltenres << " kk " << kk << " ll " << ll << " mm " << mm << " d2: " << d2 << " cutoff: " << d2cutoff << std::endl;
							}
						}
					}
				}
			}
		}

	}

	void test_symm_minimalist_ig() {
		using namespace core::chemical;
		using namespace core::conformation::symmetry;
		using namespace core::graph;
		using namespace core::pack;
		using namespace core::pack::interaction_graph;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;
		using namespace core::pose;
		using namespace core::scoring;
		using core::Size;
		using core::Real;

		Pose pose;
		core::import_pose::pose_from_pdb( pose, "core/scoring/symmetry/test_in.pdb" );
		core::pose::symmetry::make_symmetric_pose( pose );
		SymmetricConformation const & symconf( static_cast< SymmetricConformation const & > (pose.conformation() ));
		SymmetryInfoCOP symm_info = symconf.Symmetry_Info();
		Size const nres_asu = symm_info->num_independent_residues();

		PackerTaskOP task = TaskFactory::create_packer_task( pose );

		utility::vector1< bool > allowed_aas( num_canonical_aas, false );
		allowed_aas[ aa_ala ] = allowed_aas[ aa_gly ] = allowed_aas[ aa_phe ] = allowed_aas[ aa_pro ] = allowed_aas[ aa_arg ] = true;

		for ( Size ii = 1; ii <= 103; ++ii ) {
			if ( ii == 11 || ii == 12 || ii == 13 ) {
				task->nonconst_residue_task( ii ).restrict_absent_canonical_aas( allowed_aas );
			} else {
				task->nonconst_residue_task( ii ).prevent_repacking();
			}
		}

		make_symmetric_PackerTask_by_truncation( pose, task );

		ScoreFunctionOP sfxn( new core::scoring::symmetry::SymmetricScoreFunction ); //get_score_function();
		sfxn->set_weight( fa_atr, 0.8 );
		sfxn->set_weight( fa_rep, 0.44 );
		sfxn->set_weight( fa_sol, 0.65 );
		methods::EnergyMethodOptionsOP emopts( new methods::EnergyMethodOptions( sfxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		sfxn->set_energy_method_options( *emopts );

		//core::Real const initial_score = (*sfxn)( pose ); // score the pose first;
		(*sfxn)( pose ); // score the pose first;
		sfxn->setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( pose, *sfxn, task );

		core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets( new core::pack::rotamer_set::symmetry::SymmetricRotamerSets() );
		rotsets->set_task( task );
		rotsets->build_rotamers( pose, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( pose, *sfxn );


		PDInteractionGraphOP regular_ig( new PDInteractionGraph( 3 ) );
		SymmLinearMemoryInteractionGraphOP symlinmem_ig( new SymmLinearMemoryInteractionGraph( 3 ) );
		symlinmem_ig->set_pose( pose );
		symlinmem_ig->set_score_function( *sfxn );

		rotsets->compute_energies( pose, *sfxn, packer_neighbor_graph, symlinmem_ig );
		rotsets->compute_energies( pose, *sfxn, packer_neighbor_graph, regular_ig );

		regular_ig->prepare_for_simulated_annealing();
		symlinmem_ig->prepare_for_simulated_annealing();

		// OK: let's make sure that when we ask for a rotamer from a particular node on a non-asu,
		// we get back the right set of coordinates for that rotamer.
		verify_that_rotamer_coordinates_are_right( symlinmem_ig, symm_info, pose, task, rotsets );
		verify_proper_edges_exist( symlinmem_ig, symm_info, pose, sfxn, rotsets );
		verify_edge_adjacency_matrices_are_accurate( symlinmem_ig, symm_info, sfxn, task, rotsets );


		// now let's make sure that the on-the-fly interaction graph is correctly calculating the
		// change in energies for a series of rotamer substitutions.
		{
			core::PackerEnergy deltaE, prevnode_energy;
			core::PackerEnergy regular_deltaE, regular_prevnode_energy;
			for ( Size ii = 11; ii <= 13; ++ii ) {
				pose.replace_residue( ii, *rotsets->rotamer_set_for_residue( ii )->rotamer( 1 ), false );
				symlinmem_ig->consider_substitution( ii-10, 1, deltaE, prevnode_energy );
				symlinmem_ig->commit_considered_substitution();
				regular_ig->consider_substitution( ii-10, 1, deltaE, prevnode_energy );
				regular_ig->commit_considered_substitution();
			}
			core::Real state1_energy = (*sfxn)( pose );
			for ( Size ii = 11; ii <= 13; ++ii ) {
				Size iimoltenresid = ii-10;
				//std::cout << "testing rotamers from residue " << ii << std::endl;
				for ( Size jj = 1; jj <= rotsets->rotamer_set_for_residue( ii )->num_rotamers(); ++jj ) {
					//std::cout << "  rotamer #" << jj << std::endl;
					pose.replace_residue( ii, *rotsets->rotamer_set_for_residue( ii )->rotamer( jj ), false );
					core::Real jjenergy = (*sfxn)( pose );
					//if ( ii == 11 && jj == 2 ) {
					// for ( core::graph::Node::EdgeListConstIter
					//   eiter = pose.energies().energy_graph().get_node(11)->const_edge_list_begin(),
					//   eiter_end = pose.energies().energy_graph().get_node(11)->const_edge_list_end();
					//   eiter != eiter_end; ++eiter ) {
					//  EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*eiter);
					//  Size const other_node_ind = ( eedge->get_other_ind( 11 ) - 1 ) % 32 + 1;
					//  if ( other_node_ind == 12 || other_node_ind == 13 ) continue;
					//  std::cout << "pose energies: 11 rotamer 1 w/ " << eedge->get_other_ind( 11 ) << "= " << eedge->dot( sfxn->weights() ) << std::endl;
					// }
					//}
					//utility::vector1< core::Real > two_body_energy_sum( 3, 0.0 );
					//for ( Size kk = 11; kk <= 13; ++kk ) {
					// for ( core::graph::Node::EdgeListConstIter
					//   eiter = pose.energies().energy_graph().get_node(kk)->const_upper_edge_list_begin(),
					//   eiter_end = pose.energies().energy_graph().get_node(kk)->const_upper_edge_list_end();
					//   eiter != eiter_end; ++eiter ) {
					//  EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*eiter);
					//  Size const other_node_ind = ( eedge->get_other_ind( kk ) - 1 ) % 32 + 1;
					//  if ( other_node_ind == kk ) continue;
					//  if ( kk != ii && other_node_ind != ii ) continue;
					//  if ( other_node_ind == 11 || other_node_ind == 12 || other_node_ind == 13 ) {
					//   Size notii = other_node_ind == ii ? kk : other_node_ind;
					//   two_body_energy_sum[ notii-10 ] += eedge->dot( sfxn->weights() );
					//
					//   //std::cout << "edge between " << eedge->get_first_node_ind() << " " << eedge->get_second_node_ind() << " energy " << eedge->dot( sfxn->weights() ) << "      ";
					//   //std::cout << "kkcbeta: " << pose.residue( kk ).xyz( "CB" ).x() << " "  << pose.residue( kk ).xyz( "CB" ).y() << " " << pose.residue( kk ).xyz( "CB" ).z();
					//   //std::cout << "; othercbeta: " << pose.residue( eedge->get_second_node_ind() ).xyz( "CB" ).x() << " "  << pose.residue( eedge->get_second_node_ind() ).xyz( "CB" ).y() << " " << pose.residue( eedge->get_second_node_ind() ).xyz( "CB" ).z() << std::endl;
					//  }
					// }
					//}

					symlinmem_ig->consider_substitution( iimoltenresid, jj, deltaE, prevnode_energy );
					regular_ig->consider_substitution( iimoltenresid, jj, regular_deltaE, regular_prevnode_energy );


					//SymmLinearMemNode * iinode = static_cast< SymmLinearMemNode * > ( symlinmem_ig->get_node( ii-10 ) );
					//for ( int kk = 1; kk <= iinode->get_num_incident_edges(); ++kk ) {
					// Size other_node = iinode->get_index_of_adjacent_node( kk );
					// TS_ASSERT_DELTA( two_body_energy_sum[ other_node ],  jj == 1 ? iinode->get_incident_symmlinmem_edge(kk)->alt_state_energy() : iinode->get_incident_symmlinmem_edge(kk)->curr_state_energy(), 1e-5 );
					// //std::cout << "IG: edge to node " << iinode->get_index_of_adjacent_node( kk ) << " alt energy: " << iinode->get_incident_symmlinmem_edge(kk)->alt_state_energy() << std::endl;
					//}

					if ( ! (std::abs( ( deltaE - (jjenergy - state1_energy)) / ( jjenergy - state1_energy + 1e-6 ) ) < 2e-4 || ( std::abs( jjenergy - state1_energy ) < 1.0 && std::abs( deltaE - ( jjenergy - state1_energy )) < 2e-4 ) ) ) {
						if ( rotsets->rotamer_set_for_residue( ii )->rotamer( jj )->aa() == aa_gly ) {
							// iterate across edges; subtract half of the difference between the alanine and glycine bb/bb interaction energy from deltaE
							// note: there is error here in the IG which we're tolerating for the sake of speed
							for ( Size kk = 11; kk <= 13; ++kk ) {
								Size kkmoltenresid = kk-10;
								core::conformation::Residue const & kkres = symlinmem_ig->get_on_the_fly_node( kkmoltenresid )->get_rotamer( kk == ii ? jj : 1, 1 );
								for ( core::graph::Node::EdgeListConstIter
										eiter = pose.energies().energy_graph().get_node(kk)->const_upper_edge_list_begin(),
										eiter_end = pose.energies().energy_graph().get_node(kk)->const_upper_edge_list_end();
										eiter != eiter_end; ++eiter ) {
									EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*eiter);
									Size const other_node_ind = ( eedge->get_other_ind( kk ) - 1 ) % 32 + 1;
									if ( other_node_ind == kk ) continue;
									if ( other_node_ind != ii && kk != ii ) continue; // make sure one of the two nodes is node ii
									if ( other_node_ind == 11 || other_node_ind == 12 || other_node_ind == 13 ) {
										Size othernodemoltenresid = other_node_ind-10;
										Size other_subunit = (eedge->get_other_ind(kk)-1)/nres_asu + 1;
										Size ii_subunit = other_node_ind == ii ? other_subunit : 1;
										Size scale = symm_info->score_multiply( kk, eedge->get_other_ind(kk) );
										if ( scale == 0 ) continue;

										core::conformation::Residue const & otherres = symlinmem_ig->get_on_the_fly_node( othernodemoltenresid )->get_rotamer( other_node_ind == ii ? jj : 1, other_subunit );
										core::conformation::Residue const & iires = other_node_ind == ii ? otherres : kkres;
										core::conformation::Residue const & nbrres = other_node_ind == ii ? kkres : otherres;

										Real gly_bb_bbE = get_bb_bbE( pose, *sfxn, nbrres, iires ) * scale;

										core::conformation::Residue const & ii_ala_res = symlinmem_ig->get_on_the_fly_node( iimoltenresid )->get_rotamer( 1, ii_subunit );
										Real ala_bb_bbE = get_bb_bbE( pose, *sfxn, nbrres, ii_ala_res ) * scale;
										//std::cout << "adding in " << 0.5 * (gly_bb_bbE - ala_bb_bbE) << std::endl;
										deltaE += 0.5 * (gly_bb_bbE - ala_bb_bbE);
									}
								}
							}
						}
					}

					TS_ASSERT( std::abs( ( deltaE - (jjenergy - state1_energy)) / ( jjenergy - state1_energy + 1e-6 ) ) < 2e-4 || ( std::abs( jjenergy - state1_energy ) < 1.0 && std::abs( deltaE - ( jjenergy - state1_energy )) < 2e-4 ));
					if ( ! (std::abs( ( deltaE - (jjenergy - state1_energy)) / ( jjenergy - state1_energy + 1e-6 ) ) < 2e-4 || ( std::abs( jjenergy - state1_energy ) < 1.0 && std::abs( deltaE - ( jjenergy - state1_energy )) < 2e-4 ) ) ) {
						std::cout << ii << " " << jj << " bad deltaE for " << rotsets->rotamer_set_for_residue( ii )->rotamer( jj )->name() << ": " << deltaE << " vs " << jjenergy - state1_energy << " relative: " << std::abs( ( deltaE - (jjenergy - state1_energy)) / ( jjenergy - state1_energy + 1e-6 ) ) << " diff: " << std::abs( ( deltaE - (jjenergy - state1_energy))) << std::endl;


						SymmLinearMemNode * iinode = static_cast< SymmLinearMemNode * > ( symlinmem_ig->get_node( ii-10 ) );
						for ( int kk = 1; kk <= iinode->get_num_incident_edges(); ++kk ) {
							Size other_node = iinode->get_index_of_adjacent_node( kk );
							std::cout << "IG: edge to node " << iinode->get_index_of_adjacent_node( kk ) << " pro correction: " << iinode->get_incident_symmlinmem_edge(kk)->get_proline_correction_for_node( other_node, 1 ) << std::endl;
						}

						for ( Size kk = 11; kk <= 13; ++kk ) {
							Size kkmoltenresid = kk-10;
							core::conformation::Residue const & kkres = symlinmem_ig->get_on_the_fly_node( kkmoltenresid )->get_rotamer( kk == ii ? jj : 1, 1 );
							for ( core::graph::Node::EdgeListConstIter
									eiter = pose.energies().energy_graph().get_node(kk)->const_upper_edge_list_begin(),
									eiter_end = pose.energies().energy_graph().get_node(kk)->const_upper_edge_list_end();
									eiter != eiter_end; ++eiter ) {
								EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*eiter);
								Size const other_node_ind = ( eedge->get_other_ind( kk ) - 1 ) % 32 + 1;
								if ( other_node_ind == kk ) continue;
								if ( other_node_ind != ii && kk != ii ) continue; // make sure one of the two nodes is node ii
								if ( other_node_ind == 11 || other_node_ind == 12 || other_node_ind == 13 ) {
									Size othernodemoltenresid = other_node_ind-10;
									Size other_subunit = (eedge->get_other_ind(kk)-1)/nres_asu + 1;
									Size ii_subunit = other_node_ind == ii ? other_subunit : 1;
									Size scale = symm_info->score_multiply( kk, eedge->get_other_ind(kk) );
									if ( scale == 0 ) continue;

									core::conformation::Residue const & otherres = symlinmem_ig->get_on_the_fly_node( othernodemoltenresid )->get_rotamer( other_node_ind == ii ? jj : 1, other_subunit );
									core::conformation::Residue const & iires = other_node_ind == ii ? otherres : kkres;
									core::conformation::Residue const & nbrres = other_node_ind == ii ? kkres : otherres;

									std::cout << "iires: " << iires.name() << std::endl;
									Real pro_bb_bbE = get_bb_bbE( pose, *sfxn, nbrres, iires ) * scale;
									Real pro_bb_scE = get_sc_bbE( pose, *sfxn, nbrres, iires ) * scale;

									core::conformation::Residue const & ii_ala_res = symlinmem_ig->get_on_the_fly_node( iimoltenresid )->get_rotamer( 1, ii_subunit );
									std::cout << "ii_ala_res: " << ii_ala_res.name() << std::endl;
									Real nonpro_bb_bbE = get_bb_bbE( pose, *sfxn, nbrres, ii_ala_res ) * scale;
									Real nonpro_bb_scE = get_sc_bbE( pose, *sfxn, nbrres, ii_ala_res ) * scale;
									std::cout << "edge between " << eedge->get_first_node_ind() << " " << eedge->get_second_node_ind();
									std::cout << " pro_bb_bbE: " << pro_bb_bbE << " pro_bb_scE " << pro_bb_scE;
									std::cout << " nonpro_bb_bbE: " << nonpro_bb_bbE << " nonpro_bb_scE " << nonpro_bb_scE;
									std::cout << " diff: bb/bb " << pro_bb_bbE - nonpro_bb_bbE << " sc/bb " << pro_bb_scE - nonpro_bb_scE << std::endl;

									//std::cout << "edge between " << eedge->get_first_node_ind() << " " << eedge->get_second_node_ind() << " energy " << eedge->dot( sfxn->weights() ) << "      ";
									//std::cout << "kkcbeta: " << pose.residue( kk ).xyz( "CB" ).x() << " "  << pose.residue( kk ).xyz( "CB" ).y() << " " << pose.residue( kk ).xyz( "CB" ).z();
									//std::cout << "; othercbeta: " << pose.residue( eedge->get_second_node_ind() ).xyz( "CB" ).x() << " "  << pose.residue( eedge->get_second_node_ind() ).xyz( "CB" ).y() << " " << pose.residue( eedge->get_second_node_ind() ).xyz( "CB" ).z() << std::endl;
								}
							}
						}
					}
					TS_ASSERT( std::abs( ( regular_deltaE - deltaE ) / ( regular_deltaE + 1e-6 ) ) < 2e-4 || ( std::abs( regular_deltaE ) < 1.0 && std::abs( deltaE - regular_deltaE ) < 2e-4 ));
					if ( ! (std::abs( ( regular_deltaE - deltaE ) / ( regular_deltaE + 1e-6 ) ) < 2e-4 || ( std::abs( regular_deltaE ) < 1.0 && std::abs( deltaE - regular_deltaE ) < 2e-4 ) ) ) {
						std::cout << ii << " " << jj << " bad deltaE for " << rotsets->rotamer_set_for_residue( ii )->rotamer( jj )->name() << " vs regular ig: " << deltaE << " vs " << regular_deltaE << " relative: " << std::abs( ( regular_deltaE - deltaE ) / ( regular_deltaE + 1e-6 ) ) << " diff " << std::abs( ( regular_deltaE - deltaE )) << std::endl;
					}
				}
				// restore pose to state1 assigment
				pose.replace_residue( ii, *rotsets->rotamer_set_for_residue( ii )->rotamer( 1 ), false );
			}
		}

	}


	void test_sym_linmem_ig_with_asu_in_middle() {

		using namespace core::chemical;
		using namespace core::conformation::symmetry;
		using namespace core::graph;
		using namespace core::pack;
		using namespace core::pack::interaction_graph;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;
		using namespace core::pose;
		using namespace core::scoring;

		TS_ASSERT( true );
		core::io::silent::SilentFileData sfd;
		sfd.read_file( "core/pack/interaction_graph/fibril.silent" );
		std::string tag = sfd.begin()->decoy_tag();
		Pose pose; sfd.get_structure( tag ).fill_pose( pose );
		//std::cout << "sym linmem ig, asu in center; " << pose.total_residue() << std::endl;

		PackerTaskOP task = TaskFactory::create_packer_task( pose );

		SymmetricConformation & symm_conf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		SymmetryInfoCOP symm_info = symm_conf.Symmetry_Info();
		// Size const nres_asu = symm_info->num_independent_residues();

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			//if ( symm_info->bb_follows( ii ) != 0 ) continue;
			//if ( pose.residue(ii).is_virtual_residue() || ii != (79*3+50) ) {
			if ( ii < 79*3  + 1 || ii > 79*3+10  ) {
				task->nonconst_residue_task( ii ).prevent_repacking();
			} else {
				task->nonconst_residue_task( ii ).restrict_to_repacking();
				task->nonconst_residue_task( ii ).or_include_current( true );
				//std::cout << "asymm residue " << ii << std::endl;
			}
		}
		make_symmetric_PackerTask_by_truncation( pose, task );

		ScoreFunctionOP sfxn( new core::scoring::symmetry::SymmetricScoreFunction ); //get_score_function();
		sfxn->set_weight( fa_atr, 0.8 );
		sfxn->set_weight( fa_rep, 0.44 );
		sfxn->set_weight( fa_sol, 0.65 );
		methods::EnergyMethodOptionsOP emopts( new methods::EnergyMethodOptions( sfxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		sfxn->set_energy_method_options( *emopts );
		//Real const sfxnreach = sfxn->max_atomic_interaction_cutoff();

		(*sfxn)( pose ); // score the pose first;

		sfxn->setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( pose, *sfxn, task );

		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( symm_info->bb_follows( ii ) != 0 ) continue;
			for ( core::Size jj = 1; jj <= pose.total_residue(); ++jj ) {
				if ( symm_info->score_multiply( ii, jj ) == 0 ) continue;
				if ( ii == jj ) continue;
				Real dis( pose.residue( ii ).xyz( pose.residue( ii ).nbr_atom() ).distance( pose.residue( jj ).xyz( pose.residue( jj ).nbr_atom() ) ));
				Real dis_cut = pose.residue(ii).nbr_radius() + pose.residue(jj).nbr_radius() + sfxn->max_atomic_interaction_cutoff();
				if ( dis <= dis_cut ) {
					TS_ASSERT( pose.energies().energy_graph().get_edge_exists( ii, jj ) );
					TS_ASSERT( packer_neighbor_graph->get_edge_exists( ii, jj ) );
				} else {
					TS_ASSERT( ! pose.energies().energy_graph().get_edge_exists( ii, jj ) );
					TS_ASSERT( ! packer_neighbor_graph->get_edge_exists( ii, jj ) );
				}
				// if ( (ii-1) % nres_asu == 78 && (jj-1) % nres_asu == 72 ) {
				//  std::cout << "dist " << ii << " " << jj << " " << dis << " (" << dis*dis << ") cut: " << dis_cut << " (" << dis_cut*dis_cut << ")" << std::endl;
				//  std::cout << " ii radius: " << pose.residue(ii).name() << " " << pose.residue(ii).nbr_radius();
				//  std::cout << " jj radius: " << pose.residue(jj).name() << " " << pose.residue(jj).nbr_radius();
				//  std::cout << " sfxn cutoff: " << sfxn->max_atomic_interaction_cutoff() << std::endl;
				// }
			}
		}


		// Size res79_inASU = 79 + 3*nres_asu;
		// for ( Node::EdgeListConstIter
		//   iter = packer_neighbor_graph->get_node( res79_inASU )->const_edge_list_begin(),
		//   iter_end = packer_neighbor_graph->get_node( res79_inASU )->const_edge_list_end();
		//   iter != iter_end; ++iter ) {
		//  Size other = (*iter)->get_other_ind( res79_inASU );
		//  if ( (other-1) % nres_asu+1 == 73 ) {
		//   std::cout << "Edge between 79 and 73 in packer neighbor graph: " << res79_inASU << " " << other << std::endl;
		//  }
		// }
		// for ( core::Size ii = 1; ii <= 7; ++ii ) {
		//  Size ii_73 = 73 + ( ii-1 ) * nres_asu;
		//  core::conformation::Residue const & r79 = pose.residue( res79_inASU );
		//  core::conformation::Residue const & r73 = pose.residue( ii_73 );
		//  std::cout << "res " << res79_inASU << " vs res " << ii_73 << " on subunit " << ii << " sqrdist: "
		//       << r79.xyz( r79.nbr_atom() ).distance_squared( r73.xyz( r73.nbr_atom() ) ) << std::endl;
		// }
		//
		// Size res73_inASU = 73 + 3*nres_asu;
		// for ( Node::EdgeListConstIter
		//   iter = packer_neighbor_graph->get_node( res73_inASU )->const_edge_list_begin(),
		//   iter_end = packer_neighbor_graph->get_node( res73_inASU )->const_edge_list_end();
		//   iter != iter_end; ++iter ) {
		//  Size other = (*iter)->get_other_ind( res73_inASU );
		//  if ( (other-1) % nres_asu+1 == 79 ) {
		//   std::cout << "Edge between 73 and 79 in packer neighbor graph: " << res73_inASU << " " << other << std::endl;
		//  }
		// }
		// for ( core::Size ii = 1; ii <= 7; ++ii ) {
		//  Size ii_79 = 79 + ( ii-1 ) * nres_asu;
		//  core::conformation::Residue const & r73 = pose.residue( res73_inASU );
		//  core::conformation::Residue const & r79 = pose.residue( ii_79 );
		//
		//  std::cout << "res " << res73_inASU << " vs res " << ii_79 << " on subunit " << ii << " sqrdist: "
		//       << r79.xyz( r79.nbr_atom() ).distance_squared( r73.xyz( r73.nbr_atom() ) ) << std::endl;
		// }

		core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets( new core::pack::rotamer_set::symmetry::SymmetricRotamerSets() );
		rotsets->set_task( task );
		rotsets->build_rotamers( pose, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( pose, *sfxn );

		std::cout << "Creating interaction graphs with " << rotsets->nmoltenres() << " nodes" << std::endl;

		PDInteractionGraphOP regular_ig( new PDInteractionGraph( rotsets->nmoltenres() ) );
		SymmLinearMemoryInteractionGraphOP symlinmem_ig( new SymmLinearMemoryInteractionGraph( rotsets->nmoltenres() ) );
		symlinmem_ig->set_pose( pose );
		symlinmem_ig->set_score_function( *sfxn );

		rotsets->compute_energies( pose, *sfxn, packer_neighbor_graph, symlinmem_ig );
		rotsets->compute_energies( pose, *sfxn, packer_neighbor_graph, regular_ig );

		regular_ig->prepare_for_simulated_annealing();
		symlinmem_ig->prepare_for_simulated_annealing();

		verify_that_rotamer_coordinates_are_right( symlinmem_ig, symm_info, pose, task, rotsets );
		verify_proper_edges_exist( symlinmem_ig, symm_info, pose, sfxn, rotsets );
		verify_edge_adjacency_matrices_are_accurate( symlinmem_ig, symm_info, sfxn, task, rotsets );

		// let's ensure that the energies computed by the regular ig and the symlinmem ig are equivalent
		core::PackerEnergy otf_deltaE, otf_prevnode_energy;
		core::PackerEnergy reg_deltaE, reg_prevnode_energy;
		for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			symlinmem_ig->consider_substitution( ii, 1, otf_deltaE, otf_prevnode_energy );
			symlinmem_ig->commit_considered_substitution();

			regular_ig->consider_substitution( ii, 1, reg_deltaE, reg_prevnode_energy );
			regular_ig->commit_considered_substitution();
		}
		for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			Size iicurrentrot = rotsets->rotamer_set_for_moltenresidue( ii )->id_for_current_rotamer();

			symlinmem_ig->consider_substitution( ii, iicurrentrot, otf_deltaE, otf_prevnode_energy );
			symlinmem_ig->commit_considered_substitution();

			regular_ig->consider_substitution( ii, iicurrentrot, reg_deltaE, reg_prevnode_energy );
			regular_ig->commit_considered_substitution();
			//pose.replace_residue( rotsets->moltenres_2_resid(ii), *rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( 1 ), false );
			//std::cout << "current rotamer for " << ii << " is " << iicurrentrot << std::endl;
		}
		// Real ref_score = (*sfxn)( pose );

		//Size ii = 1;
		//Size jj = 1;
		//
		//symlinmem_ig->consider_substitution( ii, jj, otf_deltaE, otf_prevnode_energy );
		//regular_ig->consider_substitution(   ii, jj, reg_deltaE, reg_prevnode_energy );
		//
		//TS_ASSERT_DELTA( otf_deltaE, reg_deltaE, 1e-5 );
		//std::cout << "Delta otf_deltaE-reg_deltaE: " << otf_deltaE - reg_deltaE << std::endl;
		//
		//pose.replace_residue( rotsets->moltenres_2_resid(ii), *rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( jj ), false );
		//(*sfxn)(pose);
		//
		//for ( Size kk = 1; kk <= 7; ++kk ) {
		// Size kkresid = (rotsets->moltenres_2_resid( ii )-1) % nres_asu + 1 + (kk-1)*nres_asu;
		// core::conformation::Residue const & kkres( pose.residue( kkresid ) );
		// for ( core::graph::Node::EdgeListConstIter
		//     eiter = pose.energies().energy_graph().get_node(kkresid)->const_edge_list_begin(),
		//     eiter_end = pose.energies().energy_graph().get_node(kkresid)->const_edge_list_end();
		//    eiter != eiter_end; ++eiter ) {
		//  EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*eiter);
		//  // Size const other_node_ind = ( eedge->get_other_ind( kkresid ) - 1 ) % nres_asu + 1;
		//
		//  // Size othernodemoltenresid = other_node_ind-10;
		//  // Size other_subunit = (eedge->get_other_ind(kkresid)-1)/nres_asu + 1;
		//  // Size kk_subunit = other_node_ind == kk ? other_subunit : 1;
		//  Size scale = symm_info->score_multiply( kkresid, eedge->get_other_ind(kkresid) );
		//  if ( scale == 0 ) continue;
		//
		//  core::conformation::Residue const & otherres = pose.residue( eedge->get_other_ind( kkresid ) );
		//  Real edge_energy = eedge->dot( sfxn->weights() );
		//  Real bb_bbE = get_bb_bbE( pose, *sfxn, otherres, kkres ) * scale;
		//  Real bb_scE = get_sc_bbE( pose, *sfxn, otherres, kkres ) * scale;
		//  Real sc_bbE = get_sc_bbE( pose, *sfxn, kkres, otherres ) * scale;
		//  Real sc_scE = get_sc_scE( pose, *sfxn, otherres, kkres ) * scale;
		//  Real rpe = get_residue_pairE( pose, *sfxn, otherres, kkres ) * scale;
		//
		//  if ( edge_energy == 0 && bb_bbE == 0 && bb_scE == 0 && sc_bbE == 0 && sc_scE == 0 && rpe == 0 ) continue;
		//
		//  // core::conformation::Residue const & kk_ala_res = symlinmem_ig->get_on_the_fly_node( kkmoltenresid )->get_rotamer( 1, kk_subunit );
		//  //Real ala_bb_bbE = get_bb_bbE( pose, *sfxn, nbrres, kk_ala_res ) * scale;
		//  //std::cout << "adding in " << 0.5 * (gly_bb_bbE - ala_bb_bbE) << std::endl;
		//  //deltaE += 0.5 * (gly_bb_bbE - ala_bb_bbE);
		//  std::cout << "edge between " << eedge->get_first_node_ind() << " " << eedge->get_second_node_ind() << " energy " << edge_energy << "; ";
		//  std::cout << " bb_bb " << bb_bbE;
		//  std::cout << " bb_sc " << bb_scE;
		//  std::cout << " sc_bb " << sc_bbE;
		//  std::cout << " sc_sc " << sc_scE << " total " << bb_bbE + bb_scE + sc_bbE + sc_scE << " vs rpe " << rpe <<  std::endl;
		//  std::cout << "EnergyMap: ";
		//  EnergyMap emap = eedge->fill_energy_map();
		//  emap.show_weighted( std::cout, sfxn->weights() );
		//  std::cout << std::endl;
		//  Vector xyz = kkres.xyz(1);
		//  std::cout << "kkres       : " << xyz.x() << " " << xyz.y() << " " << xyz.z() << std::endl;
		//  xyz = otherres.xyz(1);
		//  std::cout << "otherres    : " << xyz.x() << " " << xyz.y() << " " << xyz.z() << std::endl;
		//
		// }
		//}
		//
		////pose.dump_pdb( "fibril_met1_rot1.pdb" );
		//// Size iicurrentrot = rotsets->rotamer_set_for_moltenresidue( ii )->id_for_current_rotamer();
		//
		//symlinmem_ig->commit_considered_substitution();
		//regular_ig->commit_considered_substitution();
		//
		//jj = 9;
		//symlinmem_ig->consider_substitution( ii, jj, otf_deltaE, otf_prevnode_energy );
		//regular_ig->consider_substitution(   ii, jj, reg_deltaE, reg_prevnode_energy );
		//
		//TS_ASSERT_DELTA( otf_deltaE, reg_deltaE, 1e-5 );
		//std::cout << "Delta otf_deltaE-reg_deltaE: " << otf_deltaE - reg_deltaE << std::endl;
		//
		//pose.replace_residue( rotsets->moltenres_2_resid(ii), *rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( jj ), false );
		//(*sfxn)(pose);
		//
		//for ( Size kk = 1; kk <= 7; ++kk ) {
		// Size kkresid = (rotsets->moltenres_2_resid( ii )-1) % nres_asu + 1 + (kk-1)*nres_asu;
		// core::conformation::Residue const & kkres( pose.residue( kkresid ) );
		// for ( core::graph::Node::EdgeListConstIter
		//     eiter = pose.energies().energy_graph().get_node(kkresid)->const_edge_list_begin(),
		//     eiter_end = pose.energies().energy_graph().get_node(kkresid)->const_edge_list_end();
		//    eiter != eiter_end; ++eiter ) {
		//  EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*eiter);
		//  // Size const other_node_ind = ( eedge->get_other_ind( kkresid ) - 1 ) % nres_asu + 1;
		//
		//  // Size othernodemoltenresid = other_node_ind-10;
		//  // Size other_subunit = (eedge->get_other_ind(kkresid)-1)/nres_asu + 1;
		//  // Size kk_subunit = other_node_ind == kk ? other_subunit : 1;
		//  Size scale = symm_info->score_multiply( kkresid, eedge->get_other_ind(kkresid) );
		//  if ( scale == 0 ) continue;
		//
		//  core::conformation::Residue const & otherres = pose.residue( eedge->get_other_ind( kkresid ) );
		//  Real edge_energy = eedge->dot( sfxn->weights() );
		//  Real bb_bbE = get_bb_bbE( pose, *sfxn, otherres, kkres ) * scale;
		//  Real bb_scE = get_sc_bbE( pose, *sfxn, otherres, kkres ) * scale;
		//  Real sc_bbE = get_sc_bbE( pose, *sfxn, kkres, otherres ) * scale;
		//  Real sc_scE = get_sc_scE( pose, *sfxn, otherres, kkres ) * scale;
		//  Real rpe = get_residue_pairE( pose, *sfxn, otherres, kkres ) * scale;
		//
		//  if ( edge_energy == 0 && bb_bbE == 0 && bb_scE == 0 && sc_bbE == 0 && sc_scE == 0 && rpe == 0 ) continue;
		//
		//  // core::conformation::Residue const & kk_ala_res = symlinmem_ig->get_on_the_fly_node( kkmoltenresid )->get_rotamer( 1, kk_subunit );
		//  //Real ala_bb_bbE = get_bb_bbE( pose, *sfxn, nbrres, kk_ala_res ) * scale;
		//  //std::cout << "adding in " << 0.5 * (gly_bb_bbE - ala_bb_bbE) << std::endl;
		//  //deltaE += 0.5 * (gly_bb_bbE - ala_bb_bbE);
		//  std::cout << "edge between " << eedge->get_first_node_ind() << " " << eedge->get_second_node_ind() << " energy " << edge_energy << "; ";
		//  std::cout << " bb_bb " << bb_bbE;
		//  std::cout << " bb_sc " << bb_scE;
		//  std::cout << " sc_bb " << sc_bbE;
		//  std::cout << " sc_sc " << sc_scE << " total " << bb_bbE + bb_scE + sc_bbE + sc_scE << " vs rpe " << rpe <<  std::endl;
		//  std::cout << "EnergyMap: ";
		//  EnergyMap emap = eedge->fill_energy_map();
		//  emap.show_weighted( std::cout, sfxn->weights() );
		//  std::cout << std::endl;
		//  Vector xyz = kkres.xyz(1);
		//  std::cout << "kkres       : " << xyz.x() << " " << xyz.y() << " " << xyz.z() << std::endl;
		//  xyz = otherres.xyz(1);
		//  std::cout << "otherres    : " << xyz.x() << " " << xyz.y() << " " << xyz.z() << std::endl;
		//
		// }
		//}

		for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
			for ( Size jj = 1; jj <= rotsets->rotamer_set_for_moltenresidue( ii )->num_rotamers(); ++jj ) {
				// if ( ii == 1 && jj == 9 ) {
				//  symlinmem_ig->consider_substitution( ii, 1, otf_deltaE, otf_prevnode_energy );
				//  regular_ig->consider_substitution(   ii, 1, reg_deltaE, reg_prevnode_energy );
				//  symlinmem_ig->commit_considered_substitution();
				//  regular_ig->commit_considered_substitution();
				// }

				symlinmem_ig->consider_substitution( ii, jj, otf_deltaE, otf_prevnode_energy );
				regular_ig->consider_substitution(   ii, jj, reg_deltaE, reg_prevnode_energy );

				//if ( ii == 1 && (jj == 1||jj==9) ) {
				// pose.replace_residue( rotsets->moltenres_2_resid(ii), *rotsets->rotamer_set_for_moltenresidue(ii)->rotamer( jj ), false );
				// (*sfxn)(pose);
				//
				// regular_ig->print_vertices();
				// symlinmem_ig->print_vertices();
				// std::cout << "Difference: " << otf_deltaE - reg_deltaE << std::endl;
				//
				// for ( Size kk = 1; kk <= 7; ++kk ) {
				//  Size kkresid = (rotsets->moltenres_2_resid( ii )-1) % nres_asu + 1 + (kk-1)*nres_asu;
				//  core::conformation::Residue const & kkres( pose.residue( kkresid ) );
				//  for ( core::graph::Node::EdgeListConstIter
				//      eiter = pose.energies().energy_graph().get_node(kkresid)->const_edge_list_begin(),
				//      eiter_end = pose.energies().energy_graph().get_node(kkresid)->const_edge_list_end();
				//     eiter != eiter_end; ++eiter ) {
				//   EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*eiter);
				//   // Size const other_node_ind = ( eedge->get_other_ind( kkresid ) - 1 ) % nres_asu + 1;
				//
				//   // Size othernodemoltenresid = other_node_ind-10;
				//   // Size other_subunit = (eedge->get_other_ind(kkresid)-1)/nres_asu + 1;
				//   // Size kk_subunit = other_node_ind == kk ? other_subunit : 1;
				//   Size scale = symm_info->score_multiply( kkresid, eedge->get_other_ind(kkresid) );
				//   if ( scale == 0 ) continue;
				//
				//   core::conformation::Residue const & otherres = pose.residue( eedge->get_other_ind( kkresid ) );
				//   Real edge_energy = eedge->dot( sfxn->weights() );
				//   Real bb_bbE = get_bb_bbE( pose, *sfxn, otherres, kkres ) * scale;
				//   Real bb_scE = get_sc_bbE( pose, *sfxn, otherres, kkres ) * scale;
				//   Real sc_bbE = get_sc_bbE( pose, *sfxn, kkres, otherres ) * scale;
				//   Real sc_scE = get_sc_scE( pose, *sfxn, otherres, kkres ) * scale;
				//   Real rpe = get_residue_pairE( pose, *sfxn, otherres, kkres ) * scale;
				//
				//   if ( edge_energy == 0 && bb_bbE == 0 && bb_scE == 0 && sc_bbE == 0 && sc_scE == 0 && rpe == 0 ) continue;
				//
				//   // core::conformation::Residue const & kk_ala_res = symlinmem_ig->get_on_the_fly_node( kkmoltenresid )->get_rotamer( 1, kk_subunit );
				//   //Real ala_bb_bbE = get_bb_bbE( pose, *sfxn, nbrres, kk_ala_res ) * scale;
				//   //std::cout << "adding in " << 0.5 * (gly_bb_bbE - ala_bb_bbE) << std::endl;
				//   //deltaE += 0.5 * (gly_bb_bbE - ala_bb_bbE);
				//   std::cout << "edge between " << eedge->get_first_node_ind() << " " << eedge->get_second_node_ind() << " energy " << edge_energy << "; ";
				//   std::cout << " bb_bb " << bb_bbE;
				//   std::cout << " bb_sc " << bb_scE;
				//   std::cout << " sc_bb " << sc_bbE;
				//   std::cout << " sc_sc " << sc_scE << " total " << bb_bbE + bb_scE + sc_bbE + sc_scE << " vs rpe " << rpe <<  std::endl;
				//   Vector xyz = kkres.xyz(1);
				//   std::cout << "kkres       : " << xyz.x() << " " << xyz.y() << " " << xyz.z() << std::endl;
				//   xyz = otherres.xyz(1);
				//   std::cout << "otherres    : " << xyz.x() << " " << xyz.y() << " " << xyz.z() << std::endl;
				//
				//  }
				// }
				//
				// //pose.dump_pdb( "fibril_met1_rot1.pdb" );
				// Size iicurrentrot = rotsets->rotamer_set_for_moltenresidue( ii )->id_for_current_rotamer();
				// pose.replace_residue( rotsets->moltenres_2_resid(ii), *rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( iicurrentrot ), false );
				//}

				// TS_ASSERT_DELTA( otf_deltaE, reg_deltaE, 1e-5 );
				TS_ASSERT( std::abs( ( otf_deltaE - reg_deltaE) / ( reg_deltaE + 1e-6 ) ) < 2e-4 || ( std::abs( reg_deltaE ) < 1.0 && std::abs( otf_deltaE - reg_deltaE ) < 2e-4 ));


				//std::cout << "Delta otf_deltaE-reg_deltaE: " << otf_deltaE - reg_deltaE << std::endl;
				// if ( std::abs( otf_deltaE - reg_deltaE ) > 1e-5 ) {
				//  pose.replace_residue( rotsets->moltenres_2_resid(ii), *rotsets->rotamer_set_for_moltenresidue(ii)->rotamer( jj ), false );
				//  Real new_score = (*sfxn )( pose );
				//  std::cout << "Score discrepancy: " << ii << " " << jj << " " << rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( jj )->name() << " reg= " << reg_deltaE << " otf: " << otf_deltaE << " pose: " << (new_score - ref_score ) << " ig disc: " << otf_deltaE - reg_deltaE << std::endl;
				//  Size iicurrentrot = rotsets->rotamer_set_for_moltenresidue( ii )->id_for_current_rotamer();
				//  pose.replace_residue( rotsets->moltenres_2_resid(ii), *rotsets->rotamer_set_for_moltenresidue( ii )->rotamer( iicurrentrot ), false );
				// }
			}
		}

	}

};
