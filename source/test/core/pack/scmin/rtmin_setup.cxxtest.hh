// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/SCMinMultifunc.cxxtest.hh
/// @brief  Sidechain minimization multifunc class tests
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers

// Package headers
#include <core/pack/scmin/SCMinMultifunc.hh>
#include <core/pack/scmin/SCMinMinimizerMap.hh>

#include <core/pack/rtmin.hh>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

#include <platform/types.hh>

#include <core/graph/Graph.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/DOF_Node.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>

//#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>


#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>


#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

/// Etable
#include <core/scoring/ResidueNeighborList.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>

#include <core/optimization/AtomTreeMultifunc.hh>

#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/atom_tree_minimize.hh>

#include <core/types.hh>


#include <test/UTracer.hh>

// Numeric headers
// Auto-header: duplicate removed #include <numeric/random/random.hh>

//Auto Headers
#include <core/kinematics/DomainMap.hh>
#include <core/optimization/MinimizerMap.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("core.pack.scmin.SCMinMultifunc.cxxtest");

using namespace core;


class rtmin_setup_Tests : public CxxTest::TestSuite
{

public:
	rtmin_setup_Tests() {}

	// Shared initialization goes here.
	void setUp() {
		core_init_with_additional_options( "-no_optH -ignore_unrecognized_res" );
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_rtmin_setup()
	{
		using namespace conformation;
		using namespace chemical;
		using namespace pack;
		using namespace pack::rotamer_set;
		using namespace pack::scmin;
		using namespace pose;
		using namespace scoring;
		using namespace scoring::methods;
		using namespace optimization;
		using namespace graph;


		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		//scoring::ScoreFunctionOP scorefxn = new core::scoring::ScoreFunction();

		/// 1. test etable energies -- bad?
		//scorefxn->set_weight( fa_rep, 0.44 );
		//scorefxn->set_weight( fa_atr, 0.8 );
		//scorefxn->set_weight( fa_sol, 0.6 );

		/// 2. test rama -- good
		/// scorefxn->set_weight( rama, 0.2 );

		/// 3. test dun -- good
		//scorefxn->set_weight( fa_dun, 0.4 );

		/// 4. test fa_pair -- good
		/// scorefxn->set_weight( fa_pair, 0.45 );

		/// 5. test hbonds
		//scorefxn->set_weight( hbond_sr_bb, 1.17 );
		//scorefxn->set_weight( hbond_lr_bb, 1.17 );
		//scorefxn->set_weight( hbond_bb_sc, 1.1 );
		//scorefxn->set_weight( hbond_sc, 1.1 );

		EnergyMethodOptionsOP emopts( new EnergyMethodOptions( scorefxn->energy_method_options() ) );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		scorefxn->set_energy_method_options( *emopts );

		// read in pose
		Pose pose = create_trpcage_ideal_pose();
		(*scorefxn)( pose );
		rtmin_run( pose, scorefxn );
	}

	void rtmin_run(
		pose::Pose & pose,
		scoring::ScoreFunctionOP scorefxn
	)
	{
		using namespace conformation;
		using namespace chemical;
		using namespace pack;
		using namespace pack::rotamer_set;
		using namespace pack::scmin;
		using namespace pose;
		using namespace scoring;
		using namespace scoring::methods;
		using namespace optimization;
		using namespace graph;
		//typedef utility::vector1< core::conformation::ResidueCOP > ResidueCOPs;

		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			//if ( ii % 2 == 0 || ii == 1 ) { task->nonconst_residue_task( ii ).prevent_repacking(); }
			task->nonconst_residue_task( ii ).restrict_to_repacking();
		}
		Real score_begin( (*scorefxn)(pose));
		//clock_t starttime = clock();
		RTMin rtminer;
		rtminer.rtmin( pose, *scorefxn, task );
		//clock_t stoptime = clock();
		Real score_end( (*scorefxn)(pose));
		//std::cout << "Rtmin took " << ((double) stoptime - starttime ) / CLOCKS_PER_SEC << " score begin: " << score_begin << " score_end: " << score_end << std::endl;

		TS_ASSERT( score_end < score_begin );

		/*utility::vector1< Size > inactive_neighbors; inactive_neighbors.reserve( pose.total_residue() );
		utility::vector1< bool > residue_is_inactive_neighbor( pose.total_residue(), false );
		utility::vector1< bool > active_residue_has_been_visited( pose.total_residue(), false );

		utility::vector1< Size > active_residues = pack::repackable_residues( *task );
		numeric::random::random_permutation( active_residues, numeric::random::rg() );

		utility::vector1< conformation::ResidueCOP > bgres( pose.total_residue() );
		graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( pose, *scorefxn, task );
		scoring::MinimizationGraph mingraph( pose.total_residue() );

		SCMinMinimizerMap scminmap;
		scminmap.set_total_residue( pose.total_residue() );

		EnergyMap emap_dummy;

		optimization::MinimizerOptions min_options( "dfpmin", 0.1, true, false, false );

		for ( Size ii = 1; ii <= task->num_to_be_packed(); ++ii ) {
		Size iires = active_residues[ ii ];
		for ( graph::Node::EdgeListConstIter
		eiter = packer_neighbor_graph->get_node( iires )->const_edge_list_begin(),
		eiter_end = packer_neighbor_graph->get_node( iires )->const_edge_list_end();
		eiter != eiter_end; ++eiter ) {
		Size jjres = (*eiter)->get_other_ind( iires );
		if ( ! bgres[ jjres ] && ! task->being_packed( jjres )) {
		inactive_neighbors.push_back( jjres );
		residue_is_inactive_neighbor[ jjres ] = true;
		bgres[ jjres ] = new Residue( pose.residue( jjres ) );
		scminmap.set_natoms_for_residue( jjres, bgres[ jjres ]->natoms() );
		/// Do setup_for_minimizing for background nodes once and leave them alone for
		/// the rest of the trajectory
		scorefxn->setup_for_minimizing_for_node(
		* mingraph.get_minimization_node( jjres ), pose.residue( jjres ),
		scminmap, pose, false, emap_dummy );
		}
		if ( ! task->being_packed( jjres ) || iires < jjres ) {
		mingraph.add_edge( iires, jjres ); // add edges, but don't bother calling setup_for_minimization yet
		}
		}
		}

		task->set_bump_check( false );
		task->or_include_current( true );
		task->temporarily_fix_everything();

		/// in real rtmin, the active residues will be examined in a random order;
		/// random_shuffle( active_residues );

		optimization::Multivec chi( 4 ); // guess -- resized smaller

		for ( Size ii = 1; ii <= active_residues.size(); ++ii ) {
		/// Now, build rotamers, prep the nodes and edges of the minimization graph
		/// and build the AtomTreeCollection for this residue;
		Size iiresid = active_residues[ ii ];
		conformation::Residue const & trial_res = pose.residue( iiresid );
		scminmap.activate_residue_chi( iiresid );

		//pretend this is a repacking and only this residue is being repacked
		//while all other residues are being held fixed.
		task->temporarily_set_pack_residue( iiresid, true );

		RotamerSetFactory rsf;
		rotamer_set::RotamerSetOP iirotset = rsf.create_rotamer_set( trial_res );
		iirotset->set_resid( iiresid );
		iirotset->build_rotamers( pose, *scorefxn, *task, packer_neighbor_graph );
		Size const ii_curr_rot = iirotset->id_for_current_rotamer();
		TS_ASSERT( ii_curr_rot != 0 );

		AtomTreeCollectionOP ii_atc = new AtomTreeCollection( iirotset, false );
		ii_atc->residue_atomtree_collection( iiresid ).set_active_restype_index( 1 ); // start at the beginning.
		ii_atc->residue_atomtree_collection( iiresid ).set_rescoords( * iirotset->rotamer( 1 ) );
		ii_atc->residue_atomtree_collection( iiresid ).update_atom_tree();
		scminmap.setup( ii_atc );

		{/// SCOPE -- make sure the minimization graph is ready to optimize this residue
		Residue const & iirsd( ii_atc->residue_atomtree_collection( iiresid ).active_residue() );
		if ( ! bgres[ iiresid ] ) {
		// we have not ever done setup for scoring for this residue
		scorefxn->setup_for_minimizing_for_node(
		* mingraph.get_minimization_node( iiresid ), iirsd,
		scminmap, pose, false, emap_dummy );
		} else {
		scorefxn->reinitialize_minnode_for_residue(
		* mingraph.get_minimization_node( iiresid ), iirsd,
		scminmap, pose );
		}
		for ( graph::Node::EdgeListIter
		eiter = mingraph.get_node( iiresid )->edge_list_begin(),
		eiter_end = mingraph.get_node( iiresid )->edge_list_end();
		eiter != eiter_end; ++eiter ) {

		Size jjresid = (*eiter)->get_other_ind( iiresid );
		if ( ! bgres[ jjresid ] ) {
		/// we have an active residue which we have not yet visited in the rtmin traversal
		bgres[ jjresid ] = new Residue( pose.residue( jjresid ) );
		scorefxn->setup_for_minimizing_for_node(
		* mingraph.get_minimization_node( jjresid ),
		* bgres[ jjresid ],
		scminmap, pose, false, emap_dummy );
		scminmap.set_natoms_for_residue( jjresid, bgres[ jjresid ]->natoms() );
		}
		Residue const & jjrsd( * bgres[ jjresid ] );
		MinimizationEdge & min_edge( static_cast< MinimizationEdge & > ( **eiter ));
		//std::cout << "Minedge " << iiresid << " " << jjresid << std::endl;
		if ( jjresid < iiresid ) {
		if ( residue_is_inactive_neighbor[ jjresid ] || ! active_residue_has_been_visited[ jjresid ] ) {
		scorefxn->setup_for_minimizing_sr2b_enmeths_for_minedge(
		jjrsd, iirsd, min_edge, scminmap, pose, true, false, ( EnergyEdge * ) 0, emap_dummy );
		} else {
		min_edge.reinitialize_active_energy_methods( iirsd, jjrsd, pose, true);
		}
		/// TEMP!
		min_edge.setup_for_minimizing( jjrsd, iirsd, pose, *scorefxn, scminmap );

		} else {
		if ( residue_is_inactive_neighbor[ jjresid ]  || ! active_residue_has_been_visited[ jjresid ] ) {
		scorefxn->setup_for_minimizing_sr2b_enmeths_for_minedge(
		iirsd, jjrsd, min_edge, scminmap, pose, true, false, ( EnergyEdge * ) 0, emap_dummy );
		} else {
		min_edge.reinitialize_active_energy_methods( jjrsd, iirsd, pose, true);
		}
		/// TEMP!
		min_edge.setup_for_minimizing( iirsd, jjrsd, pose, *scorefxn, scminmap );
		}
		}
		//// LONG RANGE SETUP
		for ( ScoreFunction::LR_2B_MethodIterator
		iter = scorefxn->long_range_energies_begin(),
		iter_end = scorefxn->long_range_energies_end();
		iter != iter_end; ++iter ) {

		if ( (*iter)->minimize_in_whole_structure_context( pose ) ) continue;

		LREnergyContainerCOP lrec = pose.energies().long_range_container( (*iter)->long_range_type() );
		if ( !lrec || lrec->empty() ) continue;

		EnergyMap dummy_emap;

		// Potentially O(N) operation...
		for ( ResidueNeighborConstIteratorOP
		rni = lrec->const_neighbor_iterator_begin( iiresid ), // traverse both upper and lower neighbors
		rniend = lrec->const_neighbor_iterator_end( iiresid );
		(*rni) != (*rniend); ++(*rni) ) {
		Size const r1 = rni->lower_neighbor_id();
		Size const r2 = rni->upper_neighbor_id();
		Size const jjresid = ( r1 == iiresid ? r2 : r1 );
		bool const res_moving_wrt_eachother( true );

		/// We've already set up the long-range energy methods for this edge if
		/// jjresid is an active residue that has already had its conformation optimized
		if ( active_residue_has_been_visited[ jjresid ] ) continue;
		conformation::Residue const & lower_res( r1 == iiresid ? iirsd : *bgres[ jjresid ] );
		conformation::Residue const & upper_res( r1 == iiresid ? *bgres[ jjresid ] : iirsd );
		scorefxn->setup_for_lr2benmeth_minimization_for_respair(
		lower_res, upper_res, *iter, mingraph, scminmap, pose,
		res_moving_wrt_eachother, false, rni, dummy_emap );
		}
		}

		} /// END MinimizationGraph initialization SCOPE

		/// OK: now start iterating across rotamers, setting up for minimization when the residue type changes,
		/// initializing the scminmultifunc
		Size next_restype_index( 2 );

		Real best_score( 0.0 ); bool first_pass( true );
#ifdef APL_FULL_DEBUG
		Real best_real_score( 0.0 );
#endif
		ResidueAtomTreeCollectionMomento momento;
		scminmap.set_natoms_for_residue( iiresid, ii_atc->residue_atomtree_collection( iiresid ).active_residue().natoms()  );
		for ( Size jj = 1, jj_end = iirotset->num_rotamers(); jj <= jj_end; ++jj ) {
		if (next_restype_index <= iirotset->get_n_residue_types() &&
		iirotset->get_residue_type_begin( next_restype_index ) == jj ) {
		ii_atc->residue_atomtree_collection( iiresid ).set_active_restype_index( next_restype_index );
		++next_restype_index;

		ii_atc->residue_atomtree_collection( iiresid ).set_rescoords( * iirotset->rotamer( jj ));
		ii_atc->residue_atomtree_collection( iiresid ).update_atom_tree();
		scminmap.set_natoms_for_residue( iiresid, iirotset->rotamer( jj )->natoms() );

		scminmap.setup( ii_atc ); // traverse the atom tree and identify dofs
		}

		ii_atc->residue_atomtree_collection( iiresid ).set_rescoords( * iirotset->rotamer( jj ));
		ii_atc->residue_atomtree_collection( iiresid ).update_atom_tree();
		chi = iirotset->rotamer( jj )->chi();

		reinitialize_mingraph_neighborhood_for_residue( pose, scorefxn, bgres, scminmap, scminmap.residue( iiresid ), mingraph );

#ifdef APL_FULL_DEBUG
		pose.replace_residue( iiresid, ii_atc->residue_atomtree_collection( iiresid ).active_residue(), false );
		Real const real_start_score( (*scorefxn)( pose ) );
#endif
		//pose.dump_pdb( "rtmin_before_" + utility::to_string( iiresid ) + "_" + utility::to_string( jj ) + ".pdb" );
		/// OK: Minimization graph is initialized.  Now setup the SCMinMultifunc
		SCMinMultifunc scmin_multifunc( pose, bgres, *scorefxn, mingraph, scminmap );
		//Real const start_score( scmin_multifunc( chi ) );

		//std::cout << "Starting comparison: " << iiresid << " " << start_score  << " " << iirotset->rotamer( jj )->name() << std::endl;
#ifdef APL_FULL_DEBUG
		deriv_check_for_residue( iiresid, jj, scmin_multifunc, chi );
		compare_mingraph_and_energy_graph( iiresid, pose, *scorefxn, mingraph );
#endif

		Minimizer minimizer( scmin_multifunc, min_options );
		//Real const start_func = scmin_multifunc( chi );
		//Real const end_func =
		minimizer.run( chi );
		/// Note: our neighborlist may have gone out-of-date.  Update now to make sure the best rotamer is placed in the pose
		reinitialize_mingraph_neighborhood_for_residue( pose, scorefxn, bgres, scminmap, scminmap.residue( iiresid ), mingraph );
		Real const end_score = scmin_multifunc( chi );
		//std::cout << "start func: " << start_func << " end func: " << end_func << " end score: " << end_score << std::endl;

		//for ( Size kk = 1; kk <= chi.size(); ++kk ) {
		// std::cout << "chi " << kk << " " << chi[ kk ] << " vs "
		//  << ii_atc->residue_atomtree_collection( iiresid ).active_residue().chi()[ kk ]
		//  << " ";
		//}
		//std::cout << std::endl;

#ifdef APL_FULL_DEBUG
		pose.replace_residue( iiresid, ii_atc->residue_atomtree_collection( iiresid ).active_residue(), false );
		Real const real_end_score( (*scorefxn)( pose ) );
		//std::cout << "Ending comparison: " << iiresid  << " " << end_score << " " << real_end_score << " " << iirotset->rotamer( jj )->name() << std::endl;
		deriv_check_for_residue( iiresid, jj, scmin_multifunc, chi );
		compare_mingraph_and_energy_graph( iiresid, pose, *scorefxn, mingraph );
		//if ( iiresid == 14 && jj == 7 ) {
		// atom_tree_multifunc_dump( pose, *scorefxn, chi, ii );
		//}
#endif

		if ( first_pass || end_score <= best_score ) {
		best_score = end_score;
#ifdef APL_FULL_DEBUG
		best_real_score = real_end_score;
#endif
		first_pass = false;
		ii_atc->residue_atomtree_collection( iiresid ).save_momento( momento );
		}

		// ok -- lets get here
		//std::cout << "iiresid " << iiresid << " rot: " << jj << " start score: "
		// << start_score << " end score: " << end_score << " real start: " << real_start_score
		// << " real end:" << real_end_score << " ddScore " << ( end_score - start_score ) - ( real_end_score - real_start_score )
		// << std::endl;

		//pose.dump_pdb( "rtmin_after_" + utility::to_string( iiresid ) + "_" + utility::to_string( jj ) + ".pdb" );

		}
		ii_atc->residue_atomtree_collection( iiresid ).update_from_momento( momento );
		bgres[ iiresid ] = new Residue( ii_atc->residue_atomtree_collection( iiresid ).active_residue() );

		/// NOW, we must call setup_for_scoring_for_residue for this residue we've just replaced, and
		/// for the edges adjacent to this residue and to other non-background residues so that the guarantee
		/// that setup_for_scoring_for_residue has been called on a residue before the next time its score is
		/// evaluated as a two-body energy

		//scorefxn->reinitialize_minnode_for_residue( * mingraph.get_minimization_node( iiresid ),
		// *bgres[ iiresid ], scminmap, pose );
		reinitialize_mingraph_neighborhood_for_residue( pose, scorefxn, bgres, scminmap, *bgres[ iiresid ], mingraph );

		//for ( graph::Graph::EdgeListIter
		//  edgeit = mingraph.get_node( iiresid )->edge_list_begin(),
		//  edgeit_end = mingraph.get_node( iiresid )->edge_list_end();
		//  edgeit != edgeit_end; ++edgeit ) {
		// Size const jjresid = (*edgeit)->get_other_ind( iiresid );
		// if ( residue_is_inactive_neighbor[ jjresid ] ) continue;
		// MinimizationEdge & min_edge = static_cast< MinimizationEdge & > ( (**edgeit) );
		// if ( iiresid < jjresid ) {
		//  min_edge.reinitialize_active_energy_methods( *bgres[ iiresid ], *bgres[ jjresid ], pose, true);
		//  min_edge.setup_for_minimizing( *bgres[ iiresid ], *bgres[ jjresid ], pose, *scorefxn, scminmap );
		// } else {
		//  min_edge.reinitialize_active_energy_methods( *bgres[ jjresid ], *bgres[ iiresid ], pose, true);
		//  min_edge.setup_for_minimizing( *bgres[ jjresid ], *bgres[ iiresid ], pose, *scorefxn, scminmap );
		// }
		//}

		active_residue_has_been_visited[ iiresid ] = true;
		scminmap.clear_active_chi();
#ifdef APL_FULL_DEBUG
		pose.replace_residue( iiresid, *bgres[ iiresid ], false );
		for ( Size jj = 1; jj <= bgres[ iiresid ]->natoms(); ++jj ) {
		TS_ASSERT( bgres[ iiresid ]->xyz( jj ).distance( pose.residue( iiresid ).xyz( jj ) ) < 1e-5 );
		}
		Real const ii_final_score( (*scorefxn)( pose ) );
		TS_ASSERT_DELTA( best_real_score, ii_final_score, 1e-13 );
#endif
		//pose.dump_pdb( "rtmin_selected_" + utility::to_string( iiresid ) + ".pdb" );
		//std::cout << "Round " << ii << " final score: " << (*scorefxn)(pose) << std::endl;
		}*/

	}

	void compare_mingraph_and_energy_graph(
		Size resid,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		scoring::MinimizationGraph const & mingraph
	)
	{
		using namespace pose;
		using namespace scoring;
		using namespace graph;


		EnergyMap const & one_body_emap( pose.energies().onebody_energies( resid ));
		EnergyMap min_node_1b;
		eval_res_onebody_energies_for_minnode( * mingraph.get_minimization_node( resid ), pose.residue( resid ), pose, sfxn, min_node_1b );
		for ( Size kk = 1; kk <= total_score; ++kk ) {
			ScoreType kkst = ScoreType(kk);
			if ( sfxn.weights()[ kkst ] != 0.0 ) {
				if ( std::abs( one_body_emap[ kkst ] - min_node_1b[ kkst ] ) > 1e-10 ) {
					std::cout << "   one body discrepancy " << kkst << ": " << one_body_emap[ kkst ] << " " << min_node_1b[ kkst ] << std::endl;
				}
			}
		}

		EnergyGraph const & eg( pose.energies().energy_graph() );
		for ( Node::EdgeListConstIter iter = eg.get_node(resid)->const_edge_list_begin(),
				iter_end = eg.get_node(resid)->const_edge_list_end();
				iter != iter_end; ++iter ) {
			Size ii( (*iter)->get_first_node_ind() );
			Size jj( (*iter)->get_second_node_ind() );
			if ( ii != resid && jj != resid ) continue; // only compare nodes that are involved in the current optimization
			MinimizationEdge const * minedge = static_cast< MinimizationEdge const * > ( mingraph.find_edge( ii, jj ) );
			EnergyEdge const * eedge = static_cast< EnergyEdge const * > ( *iter );
			EnergyMap emap = eedge->fill_energy_map();;

			if ( ! minedge ) {
				std::cout << "Minimization edge " << ii << " " << jj << " missing from minimization graph" << std::endl;
				emap.show_if_nonzero_weight( std::cout, pose.energies().weights() );
				continue;
			}
			bool etab_discrepancy( false );
			EnergyMap emap2;
			eval_res_pair_energy_for_minedge( *minedge, pose.residue(ii), pose.residue(jj), pose, sfxn, emap2 );
			for ( Size kk = 1; kk <= total_score; ++kk ) {
				ScoreType kkst = ScoreType(kk);
				if ( sfxn.weights()[ kkst ] != 0.0 ) {
					if ( std::abs( emap[ kkst ] - emap2[ kkst ] ) > 1e-10 ) {
						std::cout << "   " << ii << " " << jj << " " << kkst << " discrepancy: " << emap[ kkst ] << " " << emap2[ kkst ] << std::endl;
						if ( kkst == fa_atr || kkst == fa_rep || kkst == fa_sol ) {
							etab_discrepancy = true;
						}
					}
				}
			}
			if ( etab_discrepancy ) {
				using namespace scoring;
				using namespace scoring::etable;
				methods::EnergyMethodOptions options; // default is fine
				TableLookupEtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() ).lock()), options );
				ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > (minedge->res_pair_min_data().get_data_ref( etab_pair_nblist )) );
				debug_nblist_for_respair( pose.residue(ii), pose.residue(jj), pose, sfxn, etab_energy, nblist );
			}
		}
	}

	void deriv_check_for_residue(
		Size, // resid
		Size, // rotid
		pack::scmin::SCMinMultifunc & scmin_multifunc,
		utility::vector1< Real > const & chi
	)
	{
		using namespace optimization;

		//std::cout << "RTMin deriv check for residue " << resid << " rotamer #" << rotid << std::endl;
		utility::vector1< Real > dEdchi( chi );
		scmin_multifunc.dfunc( chi, dEdchi );

		bool write_output( false );//resid == 11 && rotid == 3 );
		if ( write_output ) {
			std::cout << "chi:";
			for ( Size ii = 1; ii <= chi.size(); ++ii ) {
				std::cout << " " << chi[ ii ];
			}
			std::cout << std::endl;
			std::cout << "dEdchi:";
			for ( Size ii = 1; ii <= dEdchi.size(); ++ii ) {
				std::cout << " " << dEdchi[ ii ];
			}
			std::cout << std::endl;
		}

		SimpleDerivCheckResult deriv_check = simple_numeric_deriv_check( scmin_multifunc, chi, dEdchi, write_output, write_output );
		for ( Size ii = 1; ii <= deriv_check.nangles(); ++ii ) {
			for ( Size jj = 1; jj <= 1; ++jj ) {
				//for ( Size jj = 1; jj <= deriv_check.nsteps(); ++jj ) {
				if ( deriv_check.step_data(ii,jj).ratio() != 0.0 ) {
					TS_ASSERT_DELTA( deriv_check.step_data(ii,jj).ratio(), 1.0, 0.5 );
				}
			}
			if ( deriv_check.abs_deriv_dev( ii ) != 0.0 ) {
				TS_ASSERT_DELTA( deriv_check.abs_deriv_dev( ii ), 0.0, 2e-1 );
				TS_ASSERT_DELTA( deriv_check.rel_deriv_dev( ii ), 0.0, 5e-1 );
			}
		}
		if ( deriv_check.best_abs_log_norm_ratio() != 1.0
				&& deriv_check.best_abs_log_norm_ratio() !=  100.0
				&& deriv_check.best_abs_log_norm_ratio() != -100.0  ) {
			TS_ASSERT_DELTA( deriv_check.best_abs_log_norm_ratio(), 0, 0.5 );
			TS_ASSERT_DELTA( deriv_check.best_cos_theta(), 1, 3e-1 );
		}
	}

	void atom_tree_multifunc_dump(
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		optimization::Multivec const &,// chivals,
		Size rsdno
	)
	{
		using namespace pose;
		using namespace scoring;
		using namespace kinematics;
		using namespace optimization;
		using namespace graph;

		Pose pose_copy( pose );
		sfxn( pose_copy );

		MoveMap movemap;
		movemap.set_chi( rsdno, true );

		MinimizerOptions min_options( "dfpmin_armijo", 0.01, true, false, false );

		MinimizerMap min_map;
		min_map.setup( pose_copy, movemap );

		pose_copy.energies().set_use_nblist( pose_copy, min_map.domain_map(), min_options.nblist_auto_update() );

		sfxn.setup_for_minimizing( pose_copy, min_map );

		Multivec dofs( min_map.nangles() );
		min_map.copy_dofs_from_pose( pose_copy, dofs );

		// setup the function that we will pass to the low-level minimizer
		AtomTreeMultifunc f( pose_copy, min_map, sfxn,
			min_options.deriv_check(), min_options.deriv_check_verbose() );
		Multivec dEdchi( dofs );

		f.dfunc( dofs, dEdchi );
		//std::cout << "chi (dofs):";
		for ( Size ii = 1; ii <= dofs.size(); ++ii ) {
			//std::cout << " " << dofs[ ii ];
		}
		//std::cout << "dEdchi:";
		for ( Size ii = 1; ii <= dEdchi.size(); ++ii ) {
			//std::cout << " " << dEdchi[ ii ];
		}
		//std::cout << std::endl;

		SimpleDerivCheckResult deriv_check = simple_numeric_deriv_check( f, dofs, dEdchi, true, true );

		f.dump( dofs, dofs );
	}

	void test_noop_test()
	{
		TS_ASSERT( true );
	}

	void dont_test_rtmin_for_several_pdbs() {
		utility::vector1< std::string > pdbs;
		pdbs.push_back( "1a1x_.pdb" );
		pdbs.push_back( "1a62_.pdb" );
		pdbs.push_back( "1aho_.pdb" );
		pdbs.push_back( "1arb_.pdb" );
		pdbs.push_back( "1bkb_.pdb" );
		pdbs.push_back( "1byi_.pdb" );
		pdbs.push_back( "1c52_.pdb" );
		pdbs.push_back( "1hyp_.pdb" );
		pdbs.push_back( "1koe_.pdb" );
		pdbs.push_back( "1mrj_.pdb" );
		pdbs.push_back( "1pne_.pdb" );
		pdbs.push_back( "2tgi_.pdb" );
		pdbs.push_back( "3vub_.pdb" );
		pdbs.push_back( "4mt2_.pdb" );

		Size niterations = 10;
		for ( core::Size ii = 1; ii <= pdbs.size(); ++ii ) {
			pose::Pose orig_pose;
			scoring::ScoreFunctionOP sfxn = scoring::get_score_function();
			core::import_pose::pose_from_file( orig_pose, "core/pack/" + pdbs[ ii ] , core::import_pose::PDB_file);
			(*sfxn)( orig_pose );
			//clock_t start_time = clock();
			for ( core::Size jj = 1; jj <= niterations; ++jj ) {
				pose::Pose copy_pose( orig_pose );
				//pack::task::PackerTaskOP task
				// ( pack::task::TaskFactory::create_packer_task( copy_pose ));
				//task->initialize_from_command_line().restrict_to_repacking();

				//pack::RTMin rtmin;
				//rtmin.rtmin( copy_pose, *sfxn, task );
				rtmin_run( copy_pose, sfxn );
			}
			//clock_t stop_time = clock();
			//std::cout << "RTMIN TIMING: " << pdbs[ii] << " nres " << orig_pose.total_residue() << " avg time: " << ((double) stop_time - start_time ) / ( niterations * CLOCKS_PER_SEC ) << std::endl;
		}

	}

	void reinitialize_mingraph_neighborhood_for_residue(
		pose::Pose & pose,
		scoring::ScoreFunctionCOP scorefxn,
		utility::vector1< conformation::ResidueCOP > const & bgres,
		pack::scmin::SCMinMinimizerMap const & scminmap,
		conformation::Residue const & rsd,
		scoring::MinimizationGraph & mingraph
	)
	{
		//Residue const & iires( ii_atc->residue_atomtree_collection( iiresid ).active_residue() );
		Size const resid = rsd.seqpos();
		//std::cout << "reinitialize_mingraph_neighborhood_for_residue: " << resid << std::endl;

		/// Setup the minimization graph for this new restype
		scorefxn->reinitialize_minnode_for_residue(
			* mingraph.get_minimization_node( resid ),
			rsd, scminmap, pose );
		/// Now, iterate across all the edges and set them up
		for ( graph::Node::EdgeListIter
				eiter = mingraph.get_node( resid )->edge_list_begin(),
				eiter_end = mingraph.get_node( resid )->edge_list_end();
				eiter != eiter_end; ++eiter ) {
			Size iiresid = (*eiter)->get_other_ind( resid );
			scoring::MinimizationEdge & min_edge( static_cast< scoring::MinimizationEdge & > ( **eiter ));
			if ( resid <= iiresid ) {
				min_edge.reinitialize_active_energy_methods( rsd, *bgres[ iiresid ], pose, true);
				min_edge.setup_for_minimizing( rsd, *bgres[ iiresid ], pose, *scorefxn, scminmap );
			} else {
				min_edge.reinitialize_active_energy_methods( *bgres[ iiresid ], rsd, pose, true);
				min_edge.setup_for_minimizing( *bgres[ iiresid ], rsd, pose, *scorefxn, scminmap );
			}
		}

	}

	/*
	Residue const & iires( ii_atc->residue_atomtree_collection( iiresid ).active_residue() );
	/// Setup the minimization graph for this new restype
	scorefxn->reinitialize_minnode_for_residue(
	* mingraph.get_minimization_node( iiresid ),
	iires, scminmap, pose );
	/// Now, iterate across all the edges and set them up
	for ( graph::Node::EdgeListIter
	eiter = mingraph.get_node( iiresid )->edge_list_begin(),
	eiter_end = mingraph.get_node( iiresid )->edge_list_end();
	eiter != eiter_end; ++eiter ) {
	Size kkresid = (*eiter)->get_other_ind( iiresid );
	MinimizationEdge & min_edge( static_cast< MinimizationEdge & > ( **eiter ));
	if ( iiresid <= kkresid ) {
	min_edge.reinitialize_active_energy_methods( iires, *bgres[ kkresid ], pose, true);
	min_edge.setup_for_minimizing( iires, *bgres[ kkresid ], pose, *scorefxn, scminmap );
	} else {
	min_edge.reinitialize_active_energy_methods( *bgres[ kkresid ], iires, pose, true);
	min_edge.setup_for_minimizing( *bgres[ kkresid ], iires, pose, *scorefxn, scminmap );
	}
	}
	}

	*/

	void debug_nblist_for_respair(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose,
		scoring::ScoreFunction const & sfxn,
		scoring::etable::TableLookupEtableEnergy const & etab,
		scoring::ResiduePairNeighborList const & nblist
	)
	{
		using namespace scoring;
		using namespace scoring::etable;

		utility::vector1< utility::vector1< bool > > pair_visited( rsd1.natoms() );
		for ( Size ii = 1; ii <= rsd1.natoms(); ++ii ) pair_visited[ ii ].resize( rsd2.natoms(), false );

		utility::vector1< SmallAtNb > const & neighbs( nblist.atom_neighbors() );
		for ( Size ii = 1; ii <= neighbs.size(); ++ii ) {
			pair_visited[ neighbs[ ii ].atomno1() ][ neighbs[ ii ].atomno2() ] = true;
			//conformation::Atom const & atom1( rsd1.atom( neighbs[ ii ].atomno1() ) );
			//conformation::Atom const & atom2( rsd2.atom( neighbs[ ii ].atomno2() ) );
			//atom_pair_energy( atom1, atom2, neighbs[ ii ].weight(), emap, dsq );
		}
		EnergyMap emap;
		count_pair::CountPairFunctionCOP cpfxn = etab.get_count_pair_function( rsd1, rsd2, pose, sfxn );
		bool problem( false );
		for ( Size ii = 1; ii <= rsd1.natoms(); ++ii ) {
			for ( Size jj = 1; jj <= rsd2.natoms(); ++jj ) {
				if ( pair_visited[ ii ][ jj ] ) continue;
				Real weight(1.0); Size path_dist( 0 );
				if ( ! cpfxn->count( ii, jj, weight, path_dist ) ) continue;
				emap.zero();

				Real dsq;
				conformation::Atom const & atom1( rsd1.atom( ii ) );
				conformation::Atom const & atom2( rsd2.atom( jj ) );
				etab.atom_pair_energy( atom1, atom2, weight, emap, dsq );

				Real atpairE = sfxn.weights().dot( emap );
				if ( atpairE != 0.0 ) {
					std::cout << "Atom pair " << rsd1.atom_name( ii ) << " & " << rsd2.atom_name( jj ) << " on residues " <<
						rsd1.seqpos() << " & " << rsd2.seqpos() << " have nonzero energy " << atpairE << " at distance " <<
						std::sqrt( dsq ) << std::endl;
					problem = true;
				}
			}
		}
		if ( problem ) {
			for ( Size ii = 1; ii <= neighbs.size(); ++ii ) {
				std::cout << "   neighb: " << rsd1.atom_name( neighbs[ ii ].atomno1()  ) << " & " << rsd2.atom_name( neighbs[ ii ].atomno2() ) << std::endl;
			}
		}
	}

};
