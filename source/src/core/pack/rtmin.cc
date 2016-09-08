// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rtmin.cc
/// @brief  rotamer trials with minimization module header.  Originally concieved of and implemented in Rosetta++ by Chu Wang.
/// @author Ian W. Davis (ian.w.davis@gmail.com)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com) -- reimplemented 8/2010

// Unit headers
#include <core/pack/rtmin.hh>

// Package headers

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/scmin/SCMinMultifunc.hh>
#include <core/pack/scmin/CartSCMinMultifunc.hh>
#include <core/pack/scmin/SCMinMinimizerMap.hh>
#include <core/pack/scmin/AtomTreeSCMinMinimizerMap.hh>
#include <core/pack/scmin/CartSCMinMinimizerMap.hh>
#include <core/pack/scmin/AtomTreeCollection.hh>

// Project headers
#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethod.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/optimization.OptionKeys.gen.hh>

// STL headers

#include <core/kinematics/Jump.hh>
#include <utility/vector0.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end


namespace core {
namespace pack {

using namespace ObjexxFCL::format;


static THREAD_LOCAL basic::Tracer TR( "core.pack.rtmin" );

//forward dec
utility::vector1< uint >
repackable_residues_dup( task::PackerTask const & the_task );

void reinitialize_mingraph_neighborhood_for_residue(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn,
	utility::vector1< conformation::ResidueCOP > const & bgres,
	pack::scmin::SCMinMinimizerMap const & scminmap,
	conformation::Residue const & rsd,
	scoring::MinimizationGraph & mingraph
);


RTMin::RTMin()
: // minimize_ligand_chis_(true),
	// minimize_ligand_jumps_(false),
	nonideal_(false),
	cartesian_(false)
{}

RTMin::RTMin(
	bool /*minimize_ligand_chis*/,
	bool /*minimize_ligand_jumps*/
) // :
// minimize_ligand_chis_(minimize_ligand_chis),
// minimize_ligand_jumps_(minimize_ligand_jumps)
{}

RTMin::~RTMin()= default;

/// @details Don't look, it's not pretty!
void
RTMin::rtmin(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskOP input_task
) const
{
	using namespace conformation;
	using namespace chemical;
	using namespace pack::rotamer_set;
	using namespace pack::scmin;
	using namespace pose;
	using namespace scoring;
	using namespace scoring::methods;
	using namespace optimization;
	using namespace graph;
	using namespace basic::options;

	/// 1st verify that all energy methods are compatible with rtmin.
	/// No energy method that requires whole-structure context to calculate derivatives
	/// can be used within rtmin, unless it is a whole-structure energy, in which
	/// case, RTMin will not complain, BUT the energy will not be minimized.
	bool bad( false );
	 for ( auto const & iter : scfxn.all_methods() ) {
		/// Allow whole-structure energy methods to be present.  They will not be minized!
		/// When would this be good?  If RG or chainbreak were on.
		if ( iter->method_type() != ws  && iter->minimize_in_whole_structure_context( pose ) ) {
			std::cerr << "Scoring term responsible for score types:";
			for ( Size ii = 1; ii <= iter->score_types().size(); ++ii ) {
				std::cerr << " " << iter->score_types()[ ii ];
			}
			std::cerr << " states that it requires whole-structure context to perform minimization, and thus" <<
				" cannot be used in RTMin." << std::endl;
			bad = true;
		}
	}
	if ( bad ) {
		utility_exit_with_message( "Incompatible scoring terms requested in invocation of RTMin" );
	}

	utility::vector1< Size > inactive_neighbors;
	inactive_neighbors.reserve( pose.size() );
	utility::vector1< bool > residue_is_inactive_neighbor( pose.size(), false );
	utility::vector1< bool > active_residue_has_been_visited( pose.size(), false );

	utility::vector1< Size > active_residues = pack::repackable_residues_dup( *input_task );
	numeric::random::random_permutation( active_residues, numeric::random::rg() );

	utility::vector1< conformation::ResidueCOP > bgres( pose.size() );
	graph::GraphOP packer_neighbor_graph = pack::create_packer_graph( pose, scfxn, input_task );
	scoring::MinimizationGraph mingraph( pose.size() );

	SCMinMinimizerMapOP scminmap;
	if ( cartesian_ ) {
		scminmap = SCMinMinimizerMapOP( new CartSCMinMinimizerMap() );
	} else {
		scminmap = SCMinMinimizerMapOP( new AtomTreeSCMinMinimizerMap() );
	}
	scminmap->set_nonideal( nonideal_ );
	scminmap->set_total_residue( pose.size() );

	EnergyMap emap_dummy;

	// true -- nblist, false -- deriv_check, false -- deriv_verbose
	//optimization::MinimizerOptions min_options( "lbfgs_armijo_nonmonotone", 0.1, true, false, false );
	std::string minimizer = "lbfgs";
	Size max_iter=200;
	if ( cartesian_ || nonideal_ ) {
		if ( !scfxn.ready_for_nonideal_scoring() ) {
			utility_exit_with_message( "scorefunction not set up for nonideal/Cartesian scoring" );
		}
	}

	optimization::MinimizerOptions min_options( minimizer, 0.1, true, false, false );
	min_options.max_iter(max_iter);
	min_options.silent(true);


	for ( Size ii = 1; ii <= input_task->num_to_be_packed(); ++ii ) {
		Size iires = active_residues[ ii ];
		for ( graph::Node::EdgeListConstIter
				eiter = packer_neighbor_graph->get_node( iires )->const_edge_list_begin(),
				eiter_end = packer_neighbor_graph->get_node( iires )->const_edge_list_end();
				eiter != eiter_end; ++eiter ) {
			Size jjres = (*eiter)->get_other_ind( iires );
			if ( ! bgres[ jjres ] && ! input_task->being_packed( jjres ) ) {
				inactive_neighbors.push_back( jjres );
				residue_is_inactive_neighbor[ jjres ] = true;
				bgres[ jjres ] = ResidueOP( new Residue( pose.residue( jjres ) ) );
				scminmap->set_natoms_for_residue( jjres, bgres[ jjres ]->natoms() );
				/// Do setup_for_minimizing for background nodes once and leave them alone for
				/// the rest of the trajectory
				scfxn.setup_for_minimizing_for_node(
					* mingraph.get_minimization_node( jjres ), pose.residue( jjres ),
					*scminmap, pose, false, emap_dummy );
			}
			if ( ! input_task->being_packed( jjres ) || iires < jjres ) {
				mingraph.add_edge( iires, jjres ); // add edges, but don't bother calling setup_for_minimization yet
			}
		}

		// LR2B neighbor graphs may not coincide with the
		// packer_neighbor_graph. Make the background nodes the union of
		// the background nodes from the packer_neighbor_graph and the
		// background nodes for each LR2B neighbor graph
		for ( auto
				iter = scfxn.long_range_energies_begin(),
				iter_end = scfxn.long_range_energies_end();
				iter != iter_end; ++iter ) {

			if ( (*iter)->minimize_in_whole_structure_context( pose ) ) continue;

			LREnergyContainerCOP lrec = pose.energies().long_range_container( (*iter)->long_range_type() );
			if ( !lrec || lrec->empty() ) continue;

			EnergyMap dummy_emap;

			// Potentially O(N) operation...
			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_neighbor_iterator_begin( iires ), // traverse both upper and lower neighbors
					rniend = lrec->const_neighbor_iterator_end( iires );
					(*rni) != (*rniend); ++(*rni) ) {
				Size const r1 = rni->lower_neighbor_id();
				Size const r2 = rni->upper_neighbor_id();
				Size const jjres = ( r1 == iires ? r2 : r1 );
				//bool const res_moving_wrt_eachother( true );

				if ( ! bgres[ jjres ] && ! input_task->being_packed( jjres ) ) {
					inactive_neighbors.push_back( jjres );
					residue_is_inactive_neighbor[ jjres ] = true;
					bgres[ jjres ] = ResidueOP( new Residue( pose.residue( jjres ) ) );
					scminmap->set_natoms_for_residue( jjres, bgres[ jjres ]->natoms() );
					// Do setup_for_minimizing for background nodes once and leave them alone for
					// the rest of the trajectory
					scfxn.setup_for_minimizing_for_node(
						* mingraph.get_minimization_node( jjres ), pose.residue( jjres ),
						*scminmap, pose, false, emap_dummy );
				}
				if ( ! input_task->being_packed( jjres ) || iires < jjres ) {
					if ( !mingraph.get_edge_exists(iires, jjres) ) {
						mingraph.add_edge( iires, jjres ); // add edges, but don't bother calling setup_for_minimization yet
					}
				}
			}
		}

	}


	input_task->set_bump_check( false );
	input_task->or_include_current( true );
	input_task->temporarily_fix_everything();

	/// in real rtmin, the active residues will be examined in a random order;
	/// random__shuffle( active_residues );

	Size ndofs = 4;
	if ( nonideal_ ) ndofs=14;
	if ( cartesian_ ) ndofs=75;

	optimization::Multivec chi(ndofs); // guess -- resized smaller

	for ( Size ii = 1; ii <= active_residues.size(); ++ii ) {
		/// Now, build rotamers, prep the nodes and edges of the minimization graph
		/// and build the AtomTreeCollection for this residue;
		Size iiresid = active_residues[ ii ];
		conformation::Residue const & trial_res = pose.residue( iiresid );
		scminmap->activate_residue_dofs( iiresid );

		//pretend this is a repacking and only this residue is being repacked
		//while all other residues are being held fixed.
		input_task->temporarily_set_pack_residue( iiresid, true );

		RotamerSetFactory rsf;
		rotamer_set::RotamerSetOP iirotset = rsf.create_rotamer_set( trial_res );
		iirotset->set_resid( iiresid );
		iirotset->build_rotamers( pose, scfxn, *input_task, packer_neighbor_graph );
		debug_assert( iirotset->id_for_current_rotamer() != 0 );

		AtomTreeCollectionOP ii_atc( new AtomTreeCollection( pose, *iirotset, iiresid ) );
		ii_atc->residue_atomtree_collection( iiresid ).set_active_restype_index( 1 ); // start at the beginning.
		ii_atc->residue_atomtree_collection( iiresid ).set_rescoords( * iirotset->rotamer( 1 ) );
		ii_atc->residue_atomtree_collection( iiresid ).update_atom_tree();
		scminmap->setup( ii_atc );

		{/// SCOPE -- make sure the minimization graph is ready to optimize this residue
			Residue const & iirsd( ii_atc->residue_atomtree_collection( iiresid ).active_residue() );
			if ( ! bgres[ iiresid ] ) {
				// we have not ever done setup for scoring for this residue
				scfxn.setup_for_minimizing_for_node(
					* mingraph.get_minimization_node( iiresid ), iirsd,
					*scminmap, pose, false, emap_dummy );
			} else {
				scfxn.reinitialize_minnode_for_residue(
					* mingraph.get_minimization_node( iiresid ), iirsd,
					*scminmap, pose );
			}
			for ( graph::Node::EdgeListIter
					eiter = mingraph.get_node( iiresid )->edge_list_begin(),
					eiter_end = mingraph.get_node( iiresid )->edge_list_end();
					eiter != eiter_end; ++eiter ) {

				Size jjresid = (*eiter)->get_other_ind( iiresid );
				if ( ! bgres[ jjresid ] ) {
					/// we have an active residue which we have not yet visited in the rtmin traversal
					bgres[ jjresid ] = ResidueOP( new Residue( pose.residue( jjresid ) ) );
					scfxn.setup_for_minimizing_for_node(
						* mingraph.get_minimization_node( jjresid ),
						* bgres[ jjresid ],
						*scminmap, pose, false, emap_dummy );
					scminmap->set_natoms_for_residue( jjresid, bgres[ jjresid ]->natoms() );
				}
				Residue const & jjrsd( * bgres[ jjresid ] );
				MinimizationEdge & min_edge( static_cast< MinimizationEdge & > ( **eiter ));
				//std::cout << "Minedge " << iiresid << " " << jjresid << std::endl;
				if ( jjresid < iiresid ) {
					if ( residue_is_inactive_neighbor[ jjresid ] || ! active_residue_has_been_visited[ jjresid ] ) {
						scfxn.setup_for_minimizing_sr2b_enmeths_for_minedge(
							jjrsd, iirsd, min_edge, *scminmap, pose, true, false, ( EnergyEdge * ) nullptr, emap_dummy );
					} else {
						min_edge.reinitialize_active_energy_methods( iirsd, jjrsd, pose, true);
					}
					/// hold off on the setup_for_minimizing until we know we're done with this edge (no long-range additions)
					///min_edge.setup_for_minimizing( jjrsd, iirsd, pose, scfxn, scminmap );

				} else {
					if ( residue_is_inactive_neighbor[ jjresid ]  || ! active_residue_has_been_visited[ jjresid ] ) {
						scfxn.setup_for_minimizing_sr2b_enmeths_for_minedge(
							iirsd, jjrsd, min_edge, *scminmap, pose, true, false, ( EnergyEdge * ) nullptr, emap_dummy );
					} else {
						min_edge.reinitialize_active_energy_methods( jjrsd, iirsd, pose, true);
					}
					/// hold off on the setup_for_minimizing until we know we're done with this edge (no long-range additions)
					///min_edge.setup_for_minimizing( iirsd, jjrsd, pose, scfxn, scminmap );
				}
			}
			//// LONG RANGE SETUP
			for ( auto
					iter = scfxn.long_range_energies_begin(),
					iter_end = scfxn.long_range_energies_end();
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
					scfxn.setup_for_lr2benmeth_minimization_for_respair(
						lower_res, upper_res, *iter, mingraph, *scminmap, pose,
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
		scminmap->set_natoms_for_residue( iiresid, ii_atc->residue_atomtree_collection( iiresid ).active_residue().natoms()  );
		for ( Size jj = 1, jj_end = iirotset->num_rotamers(); jj <= jj_end; ++jj ) {
			if ( next_restype_index <= iirotset->get_n_residue_types() &&
					iirotset->get_residue_type_begin( next_restype_index ) == jj ) {
				ii_atc->residue_atomtree_collection( iiresid ).set_active_restype_index( next_restype_index );
				++next_restype_index;

				ii_atc->residue_atomtree_collection( iiresid ).set_rescoords( * iirotset->rotamer( jj ));
				ii_atc->residue_atomtree_collection( iiresid ).update_atom_tree();
				scminmap->set_natoms_for_residue( iiresid, iirotset->rotamer( jj )->natoms() );

				scminmap->setup( ii_atc ); // traverse the atom tree and identify dofs
			}

			ii_atc->residue_atomtree_collection( iiresid ).set_rescoords( * iirotset->rotamer( jj ));
			ii_atc->residue_atomtree_collection( iiresid ).update_atom_tree();
			//chi = iirotset->rotamer( jj )->chi();
			scminmap->starting_dofs( chi );

			reinitialize_mingraph_neighborhood_for_residue( pose, scfxn, bgres, *scminmap, scminmap->residue( iiresid ), mingraph );

#ifdef APL_FULL_DEBUG
			pose.replace_residue( iiresid, ii_atc->residue_atomtree_collection( iiresid ).active_residue(), false );
			Real const real_start_score( scfxn( pose ) );
#endif
			//pose.dump_pdb( "rtmin_before_" + utility::to_string( iiresid ) + "_" + utility::to_string( jj ) + ".pdb" );
			/// OK: Minimization graph is initialized.  Now setup the SCMinMultifunc
			//SCMinMultifunc scmin_multifunc( pose, bgres, scfxn, mingraph, *scminmap );
			MultifuncOP scmin_multifunc = scminmap->make_multifunc( pose, bgres, scfxn, mingraph );
			//Real const start_score( scmin_multifunc( chi ) );

			//std::cout << "Starting comparison: " << iiresid << " " << start_score  << " " << iirotset->rotamer( jj )->name() << std::endl;
#ifdef APL_FULL_DEBUG
			deriv_check_for_residue( iiresid, jj, *scmin_multifunc, chi );
			compare_mingraph_and_energy_graph( iiresid, pose, scfxn, mingraph );
#endif

			Minimizer minimizer( *scmin_multifunc, min_options );
			//Real const start_func = (*scmin_multifunc)( chi );
			//Real const end_func =
			minimizer.run( chi );
			/// Note: our neighborlist may have gone out-of-date.  Update now to make sure the best rotamer is placed in the pose
			reinitialize_mingraph_neighborhood_for_residue( pose, scfxn, bgres, *scminmap, scminmap->residue( iiresid ), mingraph );
			Real const end_score = (*scmin_multifunc)( chi );
			//for ( Size kk = 1; kk <= chi.size(); ++kk ) {
			// std::cout << "chi " << kk << " " << chi[ kk ] << " vs "
			//  << ii_atc->residue_atomtree_collection( iiresid ).active_residue().chi()[ kk ]
			//  << " ";
			//}
			//std::cout << std::endl;

#ifdef APL_FULL_DEBUG
			pose.replace_residue( iiresid, ii_atc->residue_atomtree_collection( iiresid ).active_residue(), false );
			Real const real_end_score( scfxn( pose ) );
			//std::cout << "Ending comparison: " << iiresid  << " " << end_score << " " << real_end_score << " " << iirotset->rotamer( jj )->name() << std::endl;
			deriv_check_for_residue( iiresid, jj, *scmin_multifunc, chi );
			compare_mingraph_and_energy_graph( iiresid, pose, scfxn, mingraph );
			//if ( iiresid == 14 && jj == 7 ) {
			//	atom_tree_multifunc_dump( pose, scfxn, chi, ii );
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
		bgres[ iiresid ] = ResidueOP( new Residue( ii_atc->residue_atomtree_collection( iiresid ).active_residue() ) );

		/// NOW, we must call setup_for_scoring_for_residue for this residue we've just replaced, and
		/// for the edges adjacent to this residue and to other non-background residues so that the guarantee
		/// that setup_for_scoring_for_residue has been called on a residue before the next time its score is
		/// evaluated as a two-body energy

		//scfxn.reinitialize_minnode_for_residue( * mingraph.get_minimization_node( iiresid ),
		// *bgres[ iiresid ], scminmap, pose );
		reinitialize_mingraph_neighborhood_for_residue( pose, scfxn, bgres, *scminmap, *bgres[ iiresid ], mingraph );

		/*for ( graph::Graph::EdgeListIter
		edgeit = mingraph.get_node( iiresid )->edge_list_begin(),
		edgeit_end = mingraph.get_node( iiresid )->edge_list_end();
		edgeit != edgeit_end; ++edgeit ) {
		Size const jjresid = (*edgeit)->get_other_ind( iiresid );
		if ( residue_is_inactive_neighbor[ jjresid ] ) continue;

		MinimizationEdge & min_edge = static_cast< MinimizationEdge & > ( (**edgeit) );
		if ( iiresid < jjresid ) {
		min_edge.reinitialize_active_energy_methods( *bgres[ iiresid ], *bgres[ jjresid ], pose, true);
		min_edge.setup_for_minimizing( *bgres[ iiresid ], *bgres[ jjresid ], pose, scfxn, scminmap );
		} else {
		min_edge.reinitialize_active_energy_methods( *bgres[ jjresid ], *bgres[ iiresid ], pose, true);
		min_edge.setup_for_minimizing( *bgres[ jjresid ], *bgres[ iiresid ], pose, scfxn, scminmap );
		}

		}*/
		active_residue_has_been_visited[ iiresid ] = true;
		scminmap->clear_active_dofs();
		pose.replace_residue( iiresid, *bgres[ iiresid ], false );

#ifdef APL_FULL_DEBUG
		for ( Size jj = 1; jj <= bgres[ iiresid ]->natoms(); ++jj ) {
		debug_assert( bgres[ iiresid ]->xyz( jj ).distance( pose.residue( iiresid ).xyz( jj ) ) < 1e-5 );
		}
		Real const ii_final_score( scfxn( pose ) );
	debug_assert( std::abs( best_real_score - ii_final_score ) < 1e-13 );
#endif
		//pose.dump_pdb( "rtmin_selected_" + utility::to_string( iiresid ) + ".pdb" );
		//std::cout << "Round " << ii << " final score: " << scfxn(pose) << std::endl;
	}

}

/*{
using namespace numeric::random;
using namespace core::optimization;

PROF_START( basic::ROTAMER_TRIALS );
pack_scorefxn_pose_handshake( pose, scfxn);
pose.update_residue_neighbors();

utility::vector1< uint > residues_for_trials( repackable_residues_dup( *input_task ));
random_permutation( residues_for_trials, numeric::random::rg() );

task::PackerTaskOP rottrial_task( input_task->clone() );
rottrial_task->set_bump_check( false );
rottrial_task->or_include_current( true );
rottrial_task->temporarily_fix_everything();

// this will call setup fxns for each scoring method, eg HBondEnergy will
// compute backbone hbonds to prepare for hbchecking,
// PairEnergy will update actcoords...
scfxn.setup_for_packing( pose, rottrial_task->repacking_residues(),rottrial_task->designing_residues()  );

rotamer_set::RotamerSetFactory rsf;
graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scfxn, input_task );

Size const num_in_trials = residues_for_trials.size();
for (Size ii = 1; ii <= num_in_trials; ++ii)
{
pose.update_residue_neighbors(); // will return if uptodate

int const resid = residues_for_trials[ ii ];
conformation::Residue const & trial_res = pose.residue( resid );

//pretend this is a repacking and only this residue is being repacked
//while all other residues are being held fixed.
rottrial_task->temporarily_set_pack_residue( resid, true );

rotamer_set::RotamerSetOP rotset = rsf.create_rotamer_set( trial_res );
rotset->set_resid( resid );
rotset->build_rotamers( pose, scfxn, *rottrial_task, packer_neighbor_graph );
scfxn.prepare_rotamers_for_packing( pose, *rotset );
TR.Debug << "working on " << resid << " with " << rotset->num_rotamers() << " rotamers" << std::endl;

// All DOF start false (frozen)
kinematics::MoveMap movemap;
if( !pose.residue_type( resid ).is_ligand() ) movemap.set_chi(resid, true);
else{
if( minimize_ligand_chis_ ) movemap.set_chi(resid, true);
if( minimize_ligand_jumps_ ){
movemap.set_jump( pose.fold_tree().get_jump_that_builds_residue( resid ), true );
}
}
optimization::MinimizerOptions min_options("lbfgs_armijo_nonmonotone", 0.1, true , false, false);

Size best_jj = 0;
Real best_score = 1e99;
conformation::ResidueOP best_rsd;
for ( Size jj = 1; jj <= rotset->num_rotamers(); ++jj ) {
conformation::ResidueOP newresidue( rotset->rotamer( jj )->clone() );

// Assume that protein residues' conformation is fully specified by chi angles.
// This is NOT true for e.g. ligands, which may have e.g. various ring puckers.
if( newresidue->is_protein() && newresidue->type().name() == pose.residue_type(resid).name() ) {
//TR << "Setting chi angles..." << std::endl;
for( Size kk = 1; kk <= newresidue->nchi(); ++kk ) {
pose.set_chi(kk, resid, newresidue->chi(kk));
}
} else {
pose.replace_residue( resid, *newresidue, false );
scfxn.update_residue_for_packing( pose, resid );
}

// Code copied from AtomTreeMinimizer::run()
// This has to be repeated for each residue because the ResidueType may change if we're doing design.
// Even if not, we get a fatal error if we try to do it outside the loop,
// which I think is related to replace_residue() modifying the structure of the atom tree.
// It's important that the structure be scored prior to nblist setup -- why?
// A:  required for graph state == GOOD;  triggers assert in debug mode.
//Real const start_score = scfxn( pose );
// Actually, this appears to be sufficient, and is much cheaper (no twobody energy calc)
pose.scoring_begin( scfxn );
pose.scoring_end( scfxn );
// setup the map of the degrees of freedom
MinimizerMap min_map;
min_map.setup( pose, movemap );
// if we are using the nblist, set it up
if ( min_options.use_nblist() ) {
// setup a mask of the moving dofs
pose.energies().set_use_nblist( pose, min_map.domain_map(), min_options.nblist_auto_update() );
}
scfxn.setup_for_minimizing( pose, min_map );
// setup the function that we will pass to the low-level minimizer
//AtomTreeMultifunc f( pose, min_map, scfxn, min_options.deriv_check(), min_options.deriv_check_verbose() );
SingleResidueMultifunc f( pose, resid, min_map, scfxn, packer_neighbor_graph, min_options.deriv_check(), min_options.deriv_check_verbose() );
// starting position -- "dofs" = Degrees Of Freedom
Multivec dofs( min_map.nangles() );

// Code copied from AtomTreeMinimizer::run()
min_map.copy_dofs_from_pose( pose, dofs );
//Real const start_func = f( dofs );

// This actually caches the hbonds, etc.
for ( scoring::ScoreFunction::AllMethodsIterator it=scfxn.all_energies_begin(),
it_end = scfxn.all_energies_end(); it != it_end; ++it ) {
(*it)->setup_for_scoring( pose, scfxn );
}

// now do the optimization with the low-level minimizer function
Minimizer minimizer( f, min_options );
Real const score = minimizer.run( dofs );
//Real const end_func = f( dofs );
TR.Trace << "Rotamer " << jj << " " << newresidue->name3() <<
" nangles= " << min_map.nangles() <<
//" start_score: " << F(12,3,start_score) <<
//" start_func: " << F(12,3,start_func) <<
" score: "      << F(12,3,score     ) <<
//" end_func: "   << F(12,3,end_func  ) <<
std::endl;

if ( min_options.use_nblist() ) pose.energies().reset_nblist();

if(score < best_score) {
best_jj = jj;
best_score = score;
best_rsd = pose.residue(resid).clone();
}
}

if ( best_jj > 0 ) {
pose.replace_residue ( resid, *best_rsd, false );
scfxn.update_residue_for_packing( pose, resid );
}

rottrial_task->temporarily_set_pack_residue( resid, false );
}
PROF_STOP ( basic::ROTAMER_TRIALS );
}*/

void reinitialize_mingraph_neighborhood_for_residue(
	pose::Pose & pose,
	scoring::ScoreFunction const & scorefxn,
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
	scorefxn.reinitialize_minnode_for_residue(
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
			min_edge.setup_for_minimizing( rsd, *bgres[ iiresid ], pose, scorefxn, scminmap );
		} else {
			min_edge.reinitialize_active_energy_methods( *bgres[ iiresid ], rsd, pose, true);
			min_edge.setup_for_minimizing( *bgres[ iiresid ], rsd, pose, scorefxn, scminmap );
		}
	}

}


utility::vector1< uint >
repackable_residues_dup( task::PackerTask const & the_task )
{
	utility::vector1< int > to_be_packed( the_task.num_to_be_packed() );
	uint count = 0;
	for ( uint ii = 1; ii <= the_task.total_residue(); ++ii ) {
		if ( the_task.pack_residue( ii ) ) {
			++count;
			to_be_packed[ count ] = ii;
		}
	}
	return to_be_packed;
}


} //end namespace core
} //end namespace pack
