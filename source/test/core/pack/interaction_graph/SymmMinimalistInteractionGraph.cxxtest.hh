// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/SymmMinimalistInteractionGraph.cxxtest.hh
/// @brief  test suite for the symmetric, minimalist on-the-fly interaction graph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test framework headers
#include <cxxtest/TestSuite.h>

// Core Headers
#include <core/pack/interaction_graph/SymmMinimalistInteractionGraph.hh>

#include <core/chemical/AA.hh>

#include <core/graph/Graph.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/make_symmetric_task.hh>

// Test headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

//Auto Headers
#include <utility/vector1.hh>


class SymmMinimalistInteractionGraphTests : public CxxTest::TestSuite {
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
		sfxn->set_weight( atom_pair_constraint, 0.5 );

		methods::EnergyMethodOptionsOP emopts( new methods::EnergyMethodOptions( sfxn->energy_method_options() ) );
    emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
    sfxn->set_energy_method_options( *emopts );

		pose.add_constraint( core::scoring::constraints::ConstraintCOP( new core::scoring::constraints::AtomPairConstraint( core::id::AtomID( 4, 11 ), core::id::AtomID( 4, 12 ), core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc( 5.0, 1.0 ) )) ));
		//core::Real const initial_score = (*sfxn)( pose ); // score the pose first;
		(*sfxn)( pose ); // score the pose first;
		sfxn->setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );
		GraphOP packer_neighbor_graph = create_packer_graph( pose, *sfxn, task );

		core::pack::rotamer_set::symmetry::SymmetricRotamerSetsOP rotsets( new core::pack::rotamer_set::symmetry::SymmetricRotamerSets() );
		rotsets->set_task( task );
		rotsets->build_rotamers( pose, *sfxn, packer_neighbor_graph );
		rotsets->prepare_sets_for_packing( pose, *sfxn );


		PDInteractionGraphOP regular_ig( new PDInteractionGraph( 3 ) );
		SymmMinimalistInteractionGraphOP symmin_ig( new SymmMinimalistInteractionGraph( 3 ) );
		symmin_ig->set_pose( pose );
		symmin_ig->set_score_function( *sfxn );

		rotsets->compute_energies( pose, *sfxn, packer_neighbor_graph, symmin_ig );
		rotsets->compute_energies( pose, *sfxn, packer_neighbor_graph, regular_ig );

		regular_ig->prepare_for_simulated_annealing();
		symmin_ig->prepare_for_simulated_annealing();

		// OK: let's make sure that when we ask for a rotamer from a particular node on a non-asu,
		// we get back the right set of coordinates for that rotamer.
		{ // scope
			SymmMinimalistNode * node1 = static_cast< SymmMinimalistNode * > ( symmin_ig->get_node( 1 ) );
			RotamerSetCOP rotset1 = rotsets->rotamer_set_for_moltenresidue( 1 );
			for ( Size ii = 1; ii <= rotset1->num_rotamers(); ++ii ) {
				pose.replace_residue( 11, *rotset1->rotamer(ii), false );
				for ( Size jj = 1; jj <= 3; ++jj ) {
					core::conformation::Residue const & ig_res( node1->get_rotamer( ii, jj ));
					core::conformation::Residue const & pose_res( pose.residue( (jj-1)*32 + 11 ) );
					TS_ASSERT( & ig_res.type() == & pose_res.type() );
					if ( & ig_res.type() != & pose_res.type() ) continue;
					for ( Size kk = 1; kk <= ig_res.natoms(); ++kk ) {
						TS_ASSERT_DELTA( ig_res.xyz(kk).x(), pose_res.xyz(kk).x(), 1e-5 );
						TS_ASSERT_DELTA( ig_res.xyz(kk).y(), pose_res.xyz(kk).y(), 1e-5 );
						TS_ASSERT_DELTA( ig_res.xyz(kk).z(), pose_res.xyz(kk).z(), 1e-5 );
					}
				}
			}
		}

		//// Now let's make sure that the edge adjacency matrix is accurate
		{
			Size const nsubunits = 3;
			Real const sfxnreach = sfxn->max_atomic_interaction_cutoff();
			for ( int ii = 11; ii <= 13; ++ii ) {
				int const iimoltenres = ii-10;
				RotamerSetCOP iirotset = rotsets->rotamer_set_for_residue( ii );
				Size const iinrestypes = iirotset->get_n_residue_types();
				SymmMinimalistNode * iinode = static_cast< SymmMinimalistNode * > ( symmin_ig->get_node( iimoltenres ) );

				for ( int jj = 1; jj <= iinode->get_num_incident_edges(); ++jj ) {
					SymmOnTheFlyEdge * iijjedge = iinode->get_incident_otf_edge(jj);
					if ( iijjedge->get_other_ind( iimoltenres ) < iimoltenres ) continue; // only look at upper edges
					SymmOnTheFlyNode * jjnode = iinode->get_adjacent_otf_node(jj);
					Size const jjmoltenres = iijjedge->get_second_node_ind();
					//std::cout << "examining iimoltenres " << iimoltenres << " jjmoltenres " << jjmoltenres << std::endl;
					Size const jjnrestypes = rotsets->rotamer_set_for_moltenresidue( jjmoltenres )->get_n_residue_types();

					/// 1. think of iimoltenres as having a rotamer coming from the asymmetric unit, iterate across all subunits for
					/// jjmoltenres and the check every pair of residue types
					for ( Size kk = 1; kk <= nsubunits; ++kk ) {
						core::Size kk_score_multiply = symm_info->score_multiply( ii, jjmoltenres+10 + (kk-1)*nres_asu );
						//std::cout << " 1 kk= " << kk << " score multiply: " << kk_score_multiply << std::endl;
						for ( Size ll = 1; ll <= iinrestypes; ++ll ) {
							core::conformation::Residue llrot = iinode->get_rotamer( iinode->get_state_offset_for_restype_group( ll ) + 1, 1 );
							for ( Size mm = 1; mm <= jjnrestypes; ++mm ) {
								core::conformation::Residue mmrot = jjnode->get_rotamer( jjnode->get_state_offset_for_restype_group( mm ) + 1, kk );
								Real llmm_cutoff = llrot.nbr_radius() + mmrot.nbr_radius() + sfxnreach;
								Real d2cutoff = llmm_cutoff * llmm_cutoff;
								Real d2 = llrot.xyz( llrot.nbr_atom() ).distance_squared( mmrot.xyz( mmrot.nbr_atom() ) );
								Size expected_val = d2 < d2cutoff ? kk_score_multiply : 0;
								TS_ASSERT_EQUALS( iijjedge->residues_adjacent_for_subunit_pair( 1, kk, ll, mm ), expected_val );
								if ( iijjedge->residues_adjacent_for_subunit_pair( 1, kk, ll, mm ) !=  expected_val ) {
									std::cout << "iimoltenres: " << iimoltenres << " jjmoltenres: " << jjmoltenres << " kk " << kk << " ll " << ll << " mm " << mm << " d2: " << d2 << " cutoff: " << d2cutoff << std::endl;
								}
							}
						}
					}

					/// 2. think of jjmoltenres as having a rotamer coming from the asymmetric unit, iterate across all subunits for
					/// iimoltenres and the check every pair of residue types
					for ( Size kk = 1; kk <= nsubunits; ++kk ) {
						core::Size kk_score_multiply = symm_info->score_multiply( jjmoltenres+10, ii + (kk-1) * nres_asu );
						//std::cout << " 2 kk= " << kk << " score multiply: " << kk_score_multiply << std::endl;
						for ( Size ll = 1; ll <= jjnrestypes; ++ll ) {
							core::conformation::Residue llrot = jjnode->get_rotamer( jjnode->get_state_offset_for_restype_group( ll ) + 1, 1 );
							for ( Size mm = 1; mm <= iinrestypes; ++mm ) {
								core::conformation::Residue mmrot = iinode->get_rotamer( iinode->get_state_offset_for_restype_group( mm ) + 1, kk );
								Real llmm_cutoff = llrot.nbr_radius() + mmrot.nbr_radius() + sfxnreach;
								Real d2cutoff = llmm_cutoff * llmm_cutoff;
								Real d2 = llrot.xyz( llrot.nbr_atom() ).distance_squared( mmrot.xyz( mmrot.nbr_atom() ) );
								/// NOTE: do not count asu/asu interactions twice, so if kk==1, then all interactions should be marked with a 0
								Size expected_val = kk != 1 && d2 < d2cutoff ? kk_score_multiply : 0;
								TS_ASSERT_EQUALS( iijjedge->residues_adjacent_for_subunit_pair( 2, kk, ll, mm ), expected_val );
								if ( iijjedge->residues_adjacent_for_subunit_pair( 2, kk, ll, mm ) !=  expected_val ) {
									std::cout << "iimoltenres: " << iimoltenres << " jjmoltenres: " << jjmoltenres << " kk " << kk << " ll " << ll << " mm " << mm << " d2: " << d2 << " cutoff: " << d2cutoff << std::endl;
								}
							}
						}
					}
				}
			}
		}

		/// now let's make sure that the on-the-fly interaction graph is correctly calculating the
		/// change in energies for a series of rotamer substitutions.
		{
			core::PackerEnergy deltaE, prevnode_energy;
			core::PackerEnergy regular_deltaE, regular_prevnode_energy;
			for ( Size ii = 11; ii <= 13; ++ii ) {
				pose.replace_residue( ii, *rotsets->rotamer_set_for_residue( ii )->rotamer( 1 ), false );
				symmin_ig->consider_substitution( ii-10, 1, deltaE, prevnode_energy );
				symmin_ig->commit_considered_substitution();
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
					//	for ( core::graph::Node::EdgeListConstIter
					//			eiter = pose.energies().energy_graph().get_node(11)->const_edge_list_begin(),
					//			eiter_end = pose.energies().energy_graph().get_node(11)->const_edge_list_end();
					//			eiter != eiter_end; ++eiter ) {
					//		EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*eiter);
					//		Size const other_node_ind = ( eedge->get_other_ind( 11 ) - 1 ) % 32 + 1;
					//		if ( other_node_ind == 12 || other_node_ind == 13 ) continue;
					//		std::cout << "pose energies: 11 rotamer 1 w/ " << eedge->get_other_ind( 11 ) << "= " << eedge->dot( sfxn->weights() ) << std::endl;
					//	}
					//}
					//utility::vector1< core::Real > two_body_energy_sum( 3, 0.0 );
					//for ( Size kk = 11; kk <= 13; ++kk ) {
					//	for ( core::graph::Node::EdgeListConstIter
					//			eiter = pose.energies().energy_graph().get_node(kk)->const_upper_edge_list_begin(),
					//			eiter_end = pose.energies().energy_graph().get_node(kk)->const_upper_edge_list_end();
					//			eiter != eiter_end; ++eiter ) {
					//		EnergyEdge const * eedge = static_cast< EnergyEdge const * > (*eiter);
					//		Size const other_node_ind = ( eedge->get_other_ind( kk ) - 1 ) % 32 + 1;
					//		if ( other_node_ind == kk ) continue;
					//		if ( kk != ii && other_node_ind != ii ) continue;
					//		if ( other_node_ind == 11 || other_node_ind == 12 || other_node_ind == 13 ) {
					//			Size notii = other_node_ind == ii ? kk : other_node_ind;
					//			two_body_energy_sum[ notii-10 ] += eedge->dot( sfxn->weights() );
					//
					//			//std::cout << "edge between " << eedge->get_first_node_ind() << " " << eedge->get_second_node_ind() << " energy " << eedge->dot( sfxn->weights() ) << "      ";
					//			//std::cout << "kkcbeta: " << pose.residue( kk ).xyz( "CB" ).x() << " "  << pose.residue( kk ).xyz( "CB" ).y() << " " << pose.residue( kk ).xyz( "CB" ).z();
					//			//std::cout << "; othercbeta: " << pose.residue( eedge->get_second_node_ind() ).xyz( "CB" ).x() << " "  << pose.residue( eedge->get_second_node_ind() ).xyz( "CB" ).y() << " " << pose.residue( eedge->get_second_node_ind() ).xyz( "CB" ).z() << std::endl;
					//		}
					//	}
					//}

					symmin_ig->consider_substitution( iimoltenresid, jj, deltaE, prevnode_energy );
					regular_ig->consider_substitution( iimoltenresid, jj, regular_deltaE, regular_prevnode_energy );


					//SymmMinimalistNode * iinode = static_cast< SymmMinimalistNode * > ( symmin_ig->get_node( ii-10 ) );
					//for ( int kk = 1; kk <= iinode->get_num_incident_edges(); ++kk ) {
					//	Size other_node = iinode->get_index_of_adjacent_node( kk );
					//	TS_ASSERT_DELTA( two_body_energy_sum[ other_node ],  jj == 1 ? iinode->get_incident_symmin_edge(kk)->alt_state_energy() : iinode->get_incident_symmin_edge(kk)->curr_state_energy(), 1e-5 );
					//	//std::cout << "IG: edge to node " << iinode->get_index_of_adjacent_node( kk ) << " alt energy: " << iinode->get_incident_symmin_edge(kk)->alt_state_energy() << std::endl;
					//}

					if ( ! (std::abs( ( deltaE - (jjenergy - state1_energy)) / ( jjenergy - state1_energy + 1e-6 ) ) < 2e-4 || ( std::abs( jjenergy - state1_energy ) < 1.0 && std::abs( deltaE - ( jjenergy - state1_energy )) < 2e-4 ) )) {
						if ( rotsets->rotamer_set_for_residue( ii )->rotamer( jj )->aa() == aa_gly ) {
							// iterate across edges; subtract half of the difference between the alanine and glycine bb/bb interaction energy from deltaE
							// note: there is error here in the IG which we're tolerating for the sake of speed
							for ( Size kk = 11; kk <= 13; ++kk ) {
								Size kkmoltenresid = kk-10;
								core::conformation::Residue const & kkres = symmin_ig->get_on_the_fly_node( kkmoltenresid )->get_rotamer( kk == ii ? jj : 1, 1 );
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

										core::conformation::Residue const & otherres = symmin_ig->get_on_the_fly_node( othernodemoltenresid )->get_rotamer( other_node_ind == ii ? jj : 1, other_subunit );
										core::conformation::Residue const & iires = other_node_ind == ii ? otherres : kkres;
										core::conformation::Residue const & nbrres = other_node_ind == ii ? kkres : otherres;

										Real gly_bb_bbE = get_bb_bbE( pose, *sfxn, nbrres, iires ) * scale;

										core::conformation::Residue const & ii_ala_res = symmin_ig->get_on_the_fly_node( iimoltenresid )->get_rotamer( 1, ii_subunit );
										Real ala_bb_bbE = get_bb_bbE( pose, *sfxn, nbrres, ii_ala_res ) * scale;
										//std::cout << "adding in " << 0.5 * (gly_bb_bbE - ala_bb_bbE) << std::endl;
										deltaE += 0.5 * (gly_bb_bbE - ala_bb_bbE);
									}
								}
							}
						}
					}

					TS_ASSERT( std::abs( ( deltaE - (jjenergy - state1_energy)) / ( jjenergy - state1_energy + 1e-6 ) ) < 2e-4 || ( std::abs( jjenergy - state1_energy ) < 1.0 && std::abs( deltaE - ( jjenergy - state1_energy )) < 2e-4 ));
					if ( ! (std::abs( ( deltaE - (jjenergy - state1_energy)) / ( jjenergy - state1_energy + 1e-6 ) ) < 2e-4 || ( std::abs( jjenergy - state1_energy ) < 1.0 && std::abs( deltaE - ( jjenergy - state1_energy )) < 2e-4 ) )) {
						std::cout << ii << " " << jj << " bad deltaE for " << rotsets->rotamer_set_for_residue( ii )->rotamer( jj )->name() << ": " << deltaE << " vs " << jjenergy - state1_energy << " relative: " << std::abs( ( deltaE - (jjenergy - state1_energy)) / ( jjenergy - state1_energy + 1e-6 ) ) << " diff: " << std::abs( ( deltaE - (jjenergy - state1_energy))) << std::endl;


						SymmMinimalistNode * iinode = static_cast< SymmMinimalistNode * > ( symmin_ig->get_node( ii-10 ) );
						for ( int kk = 1; kk <= iinode->get_num_incident_edges(); ++kk ) {
							Size other_node = iinode->get_index_of_adjacent_node( kk );
							std::cout << "IG: edge to node " << iinode->get_index_of_adjacent_node( kk ) << " pro correction: " << iinode->get_incident_symmin_edge(kk)->get_proline_correction_for_node( other_node, 1 ) << std::endl;
						}

						for ( Size kk = 11; kk <= 13; ++kk ) {
							Size kkmoltenresid = kk-10;
							core::conformation::Residue const & kkres = symmin_ig->get_on_the_fly_node( kkmoltenresid )->get_rotamer( kk == ii ? jj : 1, 1 );
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

									core::conformation::Residue const & otherres = symmin_ig->get_on_the_fly_node( othernodemoltenresid )->get_rotamer( other_node_ind == ii ? jj : 1, other_subunit );
									core::conformation::Residue const & iires = other_node_ind == ii ? otherres : kkres;
									core::conformation::Residue const & nbrres = other_node_ind == ii ? kkres : otherres;

									std::cout << "iires: " << iires.name() << std::endl;
									Real pro_bb_bbE = get_bb_bbE( pose, *sfxn, nbrres, iires ) * scale;
									Real pro_bb_scE = get_sc_bbE( pose, *sfxn, nbrres, iires ) * scale;

									core::conformation::Residue const & ii_ala_res = symmin_ig->get_on_the_fly_node( iimoltenresid )->get_rotamer( 1, ii_subunit );
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
					if ( ! (std::abs( ( regular_deltaE - deltaE ) / ( regular_deltaE + 1e-6 ) ) < 2e-4 || ( std::abs( regular_deltaE ) < 1.0 && std::abs( deltaE - regular_deltaE ) < 2e-4 ) )) {
						std::cout << ii << " " << jj << " bad deltaE for " << rotsets->rotamer_set_for_residue( ii )->rotamer( jj )->name() << " vs regular ig: " << deltaE << " vs " << regular_deltaE << " relative: " << std::abs( ( regular_deltaE - deltaE ) / ( regular_deltaE + 1e-6 ) ) << " diff " << std::abs( ( regular_deltaE - deltaE )) << std::endl;
					}
				}
				// restore pose to state1 assigment
				pose.replace_residue( ii, *rotsets->rotamer_set_for_residue( ii )->rotamer( 1 ), false );
			}
		}
	}


};
