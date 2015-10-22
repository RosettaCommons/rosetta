// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/min_pack.cc
/// @brief  pack rotamers with minimization module
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/min_pack.hh>

// Package Headers
#include <core/pack/rtmin.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/scmin/AtomTreeCollection.hh>
#include <core/pack/scmin/SCMinMinimizerMap.hh>
#include <core/pack/scmin/AtomTreeSCMinMinimizerMap.hh>
#include <core/pack/scmin/CartSCMinMinimizerMap.hh>
#include <core/pack/scmin/SidechainStateAssignment.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/ContinuousRotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>

//#include <core/conformation/symmetry/SymmetryInfo.hh>
//#include <core/pose/symmetry/util.hh>


// Project Headers
#include <core/conformation/Residue.hh>
#include <core/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>


// util
#include <basic/Tracer.hh>

/// ObjexxFCL

// Numeric
#include <numeric/random/random.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <utility/vector1.hh>

// APL TEEEEMMMPPP
// #include <ctime>

namespace core {
namespace pack {

//#define APL_FULL_DEBUG

void compare_mingraph_and_energy_graph(
	Size resid,
	pose::Pose const & pose,
	scoring::ScoreFunction const & sfxn,
	scoring::MinimizationGraph const & mingraph
);

static THREAD_LOCAL basic::Tracer TR( "core.pack.min_pack", basic::t_info );

scmin::SCMinMinimizerMapOP
create_scmin_minimizer_map(
	pose::Pose const & pose,
	task::PackerTaskCOP task,
	bool cartesian
)
{
	/// create the scminmap and activate all chi for all residues being packed.
	/// This is important for the initialization of the MinimizationGraph; the
	/// domain map that the scminmap provides to the ScoreFunction must state
	/// that each of the molten residues have color 0.
	scmin::SCMinMinimizerMapOP scminmap;
	if ( cartesian ) {
		scminmap = scmin::SCMinMinimizerMapOP( new scmin::CartSCMinMinimizerMap() );
	} else {
		scminmap = scmin::SCMinMinimizerMapOP( new scmin::AtomTreeSCMinMinimizerMap() );
	}

	scminmap->set_total_residue( pose.total_residue() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( task->being_packed( ii ) ) {
			scminmap->activate_residue_dofs( ii );
		} else {
			scminmap->set_natoms_for_residue( ii, pose.residue( ii ).natoms() );
		}
	}
	return scminmap;
}


scoring::MinimizationGraphOP
create_minimization_graph(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTask const & task,
	graph::Graph const & packer_neighbor_graph,
	scmin::SCMinMinimizerMap const & sc_min_map
)
{
	scoring::MinimizationGraphOP mingraph( new scoring::MinimizationGraph( pose.total_residue() ) );

	/// copy all edges that are incident upon the molten residues; background edges may be ignored.
	for ( graph::Graph::EdgeListConstIter
			eiter = packer_neighbor_graph.const_edge_list_begin(),
			eiter_end = packer_neighbor_graph.const_edge_list_end();
			eiter != eiter_end; ++eiter ) {
		if ( task.being_packed( (*eiter)->get_first_node_ind() ) || task.being_packed( (*eiter)->get_second_node_ind() ) ) {
			mingraph->add_edge( (*eiter)->get_first_node_ind(), (*eiter)->get_second_node_ind() );
		}
	}

	/// now, setup the nodes and edges.
	/// Initialize any node being packed, or that has at least one
	/// neighbor being packed (i.e. it has at least one incident
	/// edge in the minimization graph)

	scoring::EnergyMap dummy;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( task.being_packed( ii ) || mingraph->get_node( ii )->num_edges() != 0 ) {
			sfxn.setup_for_minimizing_for_node(
				* mingraph->get_minimization_node( ii ), pose.residue( ii ),
				sc_min_map, pose, false, dummy );
		}
	}
	for ( graph::Graph::EdgeListIter
			eiter = mingraph->edge_list_begin(), eiter_end = mingraph->edge_list_end();
			eiter != eiter_end; ++eiter ) {
		Size const node1 = (*eiter)->get_first_node_ind();
		Size const node2 = (*eiter)->get_second_node_ind();

		scoring::MinimizationEdge & minedge( static_cast< scoring::MinimizationEdge & > (**eiter) );

		sfxn.setup_for_minimizing_sr2b_enmeths_for_minedge(
			pose.residue( node1 ), pose.residue( node2 ),
			minedge, sc_min_map, pose, true, false,
			static_cast< scoring::EnergyEdge const * > ( 0 ), dummy );
	}

	/// Now initialize the long-range edges
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( ! task.being_packed( ii ) ) continue;
		for ( scoring::ScoreFunction::LR_2B_MethodIterator
				iter = sfxn.long_range_energies_begin(),
				iter_end = sfxn.long_range_energies_end();
				iter != iter_end; ++iter ) {

			if ( (*iter)->minimize_in_whole_structure_context( pose ) ) continue;

			scoring::LREnergyContainerCOP lrec = pose.energies().long_range_container( (*iter)->long_range_type() );
			if ( !lrec || lrec->empty() ) continue;

			// Potentially O(N) operation...
			for ( scoring::ResidueNeighborConstIteratorOP
					rni = lrec->const_neighbor_iterator_begin( ii ), // traverse both upper and lower neighbors
					rniend = lrec->const_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				Size const r1 = rni->lower_neighbor_id();
				Size const r2 = rni->upper_neighbor_id();
				Size const jj = ( r1 == ii ? r2 : r1 );
				bool const res_moving_wrt_eachother( true );

				if ( task.being_packed( jj ) && jj < ii ) continue; // only setup each edge once.

				conformation::Residue const & lower_res( r1 == ii ? pose.residue( ii ) : pose.residue( jj ) );
				conformation::Residue const & upper_res( r1 == ii ? pose.residue( jj ) : pose.residue( ii ) );
				sfxn.setup_for_lr2benmeth_minimization_for_respair(
					lower_res, upper_res, *iter, *mingraph, sc_min_map, pose,
					res_moving_wrt_eachother, false, rni, dummy );
			}
		}
	}
	return mingraph;
}


utility::vector1< conformation::ResidueCOP >
setup_bgres_cops(
	pose::Pose const & pose,
	task::PackerTask const & task
)
{
	utility::vector1< conformation::ResidueCOP > bgres( pose.total_residue() );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( ! task.being_packed( ii ) ) bgres[ ii ] = pose.residue( ii ).clone();
	}
	return bgres;
}


utility::vector1< Real >
initialize_temperatures(
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::packing;

	if ( option[ minpack_temp_schedule ].user() ) {
		return option[ minpack_temp_schedule ]();
	}

	/// "kTsched5" in APL's tuning of the min packer
	utility::vector1< Real > temps( 12 );
	temps[ 1  ] = 33;
	temps[ 2  ] = 10;
	temps[ 3  ] = 5;
	temps[ 4  ] = 1;
	temps[ 5  ] = .6;
	temps[ 6  ] = .45;
	temps[ 7  ] = .4;
	temps[ 8  ] = .3;
	temps[ 9  ] = .2;
	temps[ 10 ] = .1;
	temps[ 11 ] = .05;
	temps[ 12 ] = .03;
	return temps;
}

utility::vector1< Real >
initialize_temperatures_off_rotamer_pack(
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::packing;

	if ( option[ minpack_temp_schedule ].user() ) {
		return option[ minpack_temp_schedule ]();
	}

	/// "kTsched5" in APL's tuning of the min packer
	utility::vector1< Real > temps; temps.reserve( 15 );

	temps.push_back( 20   );
	temps.push_back( 10   );
	temps.push_back( 3    );
	temps.push_back( 1    );
	temps.push_back( .6   );
	temps.push_back( .4   );
	temps.push_back( .3   );
	temps.push_back( .2   );
	temps.push_back( .15  );
	temps.push_back( .1   );
	temps.push_back( .075 );
	temps.push_back( .05  );
	temps.push_back( .03  );
	temps.push_back( .02  );
	temps.push_back( .01  );

	return temps;
}


Real
get_residue_current_energy(
	pose::Pose & pose,
	utility::vector1< conformation::ResidueCOP > const & bgres,
	scoring::ScoreFunction const & sfxn,
	scoring::MinimizationGraph & mingraph,
	scmin::SCMinMinimizerMap & scminmap,
	scmin::SidechainStateAssignment const & curr_state,
	scmin::AtomTreeCollectionOP atc,
	rotamer_set::RotamerSets const & rotsets,
	Size resid,
	Size moltenres_id
#ifdef APL_FULL_DEBUG
	, pose::Pose & debug_pose
#endif
)
{
	/// This will be faster after I subclass the minimization graph and store energies on its edges.
	/// For now, recomupute the energies for this residue

#ifdef APL_FULL_DEBUG

	/// the background residue should already contain the correct coordinates
	for ( Size ii = 1; ii <= debug_pose.total_residue(); ++ii ) {
		if ( bgres[ ii ]->natoms() != debug_pose.residue( ii ).natoms() ) {
			std::cout << "Residue " << ii << " natoms discrepancy" << std::endl;
		}
		for ( Size jj = 1; jj <= debug_pose.residue( ii ).natoms(); ++jj ) {
			if ( bgres[ ii ]->xyz( jj ).distance_squared( debug_pose.residue(ii).xyz(jj) ) > 1e-6 ) {
				std::cout << "Coordinate discrepancy between debug pose and background residue " << ii << " at atom " << jj << std::endl;
			}
		}
	}
	for ( Size ii = 1; ii <= debug_pose.residue( resid ).natoms(); ++ii ) {
		if ( curr_state.momento_for_moltenres( moltenres_id ).coord( ii ).distance_squared( debug_pose.residue( resid ).xyz( ii ) ) > 1e-6 ) {
			std::cout << "  Momento and debug pose coordinate discrepancy " << resid << " " << moltenres_id << " " << ii << std::endl;
		}
	}
	std::cout << " get_residue_current_energy: " << resid << std::endl;
#endif


	/// 1.
	scmin::ResidueAtomTreeCollection & ratc = atc->residue_atomtree_collection( resid );
	Size curr_rot = curr_state.orig_rotamer_id_for_moltenres( moltenres_id );
	Size const restype_index_for_curr_rotamer = rotsets.rotamer_set_for_residue( resid )->get_residue_type_index_for_rotamer( curr_rot );
	ratc.set_active_restype_index( restype_index_for_curr_rotamer );
	ratc.update_from_momento( curr_state.momento_for_moltenres( moltenres_id ) );

#ifdef APL_FULL_DEBUG
	for ( Size ii = 1; ii <= ratc.active_residue().natoms(); ++ii ) {
		// AMW cppcheck: misplaced )
		if ( std::abs( ratc.active_residue().xyz( ii ).distance_squared( curr_state.momento_for_moltenres( moltenres_id ).coord( ii ) ) ) > 1e-6 ) {
			std::cout << "  RATC and momento coordinate discrepancy before func eval " << ii << std::endl;
		}
	}
	for ( Size ii = 1; ii <= ratc.active_residue().chi().size(); ++ii ) {
		Real ratc_chi_ii = ratc.active_residue().chi()[ ii ] <= -180 ? ratc.active_residue().chi()[ ii ] + 360 : ratc.active_residue().chi()[ ii ];
		Real dbpose_chi_ii = debug_pose.residue( resid ).chi()[ ii ] <= -180 ? debug_pose.residue( resid ).chi()[ ii ] + 360 : debug_pose.residue( resid ).chi()[ ii ];
		if ( std::abs( ratc_chi_ii - dbpose_chi_ii ) > 1e-6 && std::abs( std::abs( ratc_chi_ii - dbpose_chi_ii) - 360 ) > 1e-3 ) {
			std::cout << "  Chi discrepancy between momento-restored ratc and debug pose " << ratc_chi_ii << " vs " << dbpose_chi_ii << " " << resid << " " << debug_pose.residue( resid ).name() << " chi: " << ii << std::endl;
		}
	}
#endif

	/// 2.
	scminmap.clear_active_dofs();
	scminmap.activate_residue_dofs( resid );
	scminmap.set_natoms_for_residue( resid, atc->residue_atomtree_collection( resid ).active_restype().natoms() );
	scminmap.setup( atc );

	/// 3.
	reinitialize_mingraph_neighborhood_for_residue( pose, sfxn, bgres, scminmap, ratc.active_residue(), mingraph );

	/// 4.
	optimization::MultifuncOP scmin_func = scminmap.make_multifunc( pose, bgres, sfxn, mingraph );
	utility::vector1< Real > chi;
	scminmap.starting_dofs( chi );

	Real funcval = (*scmin_func)( chi );

#ifdef APL_FULL_DEBUG
	bool bad = false;
	for ( Size ii = 1; ii <= ratc.active_residue().natoms(); ++ii ) {
		if ( std::abs( ratc.active_residue().xyz( ii ).distance_squared( curr_state.momento_for_moltenres( moltenres_id ).coord( ii ) ) ) > 1e-6 ){
			std::cout << "  RATC and momento coordinate discrepancy! " << ii << std::endl;
			bad = true;
		}
	}
	if ( bad ) {
		std::cout << " measured chi " << std::endl;
		conformation::Residue const & r( ratc.active_residue() );
		for ( Size ii = 1; ii <= chi.size(); ++ii ) {
			std::cout << " ratc desired: " << chi[ ii ] << " actual: " << numeric::dihedral_degrees( r.xyz( r.chi_atoms(ii)[ 1 ]), r.xyz( r.chi_atoms(ii)[2]), r.xyz( r.chi_atoms(ii)[3]), r.xyz( r.chi_atoms(ii)[4] ) ) << std::endl;
			std::cout << " momento desired: " << chi[ ii ] << " actual: " << numeric::dihedral_degrees( curr_state.momento_for_moltenres( moltenres_id ).coord( r.chi_atoms(ii)[ 1 ]), curr_state.momento_for_moltenres( moltenres_id ).coord( r.chi_atoms(ii)[2]), curr_state.momento_for_moltenres( moltenres_id ).coord( r.chi_atoms(ii)[3]), curr_state.momento_for_moltenres( moltenres_id ).coord( r.chi_atoms(ii)[4] ) ) << std::endl;
		}
	}

	compare_mingraph_and_energy_graph(resid, debug_pose, sfxn, mingraph );

	for ( Size ii = 1; ii <= chi.size(); ++ii ) {
		if ( std::abs( ratc.active_residue().chi()[ ii ] - chi[ ii ] ) > 1e-3 ) {
			std::cout << "chi discrepancy " << ratc.active_residue().chi()[ ii ] << " vs " << chi[ ii ] << std::endl;
		}
	}
	for ( Size ii = 1; ii <= debug_pose.residue( resid ).natoms(); ++ii ) {
		if ( ratc.active_residue().xyz( ii ).distance_squared( debug_pose.residue( resid ).xyz( ii ) ) > 1e-6 ) {
			std::cout << "Coordinate discrepancy between debug pose and RATC for residue " << resid << " at atom " << ii << std::endl;
		}
	}
	for ( Size ii = 1; ii <= debug_pose.total_residue(); ++ii ) {
		if ( bgres[ ii ]->natoms() != debug_pose.residue( ii ).natoms() ) {
			std::cout << "Residue " << ii << " natoms discrepancy" << std::endl;
		}
		for ( Size jj = 1; jj <= debug_pose.residue( ii ).natoms(); ++jj ) {
			if ( bgres[ ii ]->xyz( jj ).distance_squared( debug_pose.residue(ii).xyz(jj) ) > 1e-6 ) {
				std::cout << "Coordinate discrepancy between debug pose and background residue " << ii << " at atom " << jj << std::endl;
			}
		}
	}

#endif

	return funcval;
}

Real
get_total_energy_for_state(
	pose::Pose & pose,
	utility::vector1< conformation::ResidueCOP > const & bgres,
	scoring::ScoreFunction const & sfxn,
	scoring::MinimizationGraph & mingraph,
	scmin::SCMinMinimizerMap & scminmap,
	scmin::SidechainStateAssignment const & curr_state,
	scmin::AtomTreeCollectionOP atc,
	rotamer_set::RotamerSets const & rotsets
)
{
	/// This will be faster after I subclass the minimization graph and store energies on its edges.
	/// For now, recomupute the energies for this residue

	/// 1.
	for ( Size ii = 1; ii <= rotsets.nmoltenres(); ++ii ) {
		Size iiresid = rotsets.moltenres_2_resid( ii );
		scmin::ResidueAtomTreeCollection & iiratc = atc->residue_atomtree_collection( iiresid );
		Size curr_rot = curr_state.orig_rotamer_id_for_moltenres( ii );
		Size const restype_index_for_curr_rotamer = rotsets.rotamer_set_for_residue( iiresid )->get_residue_type_index_for_rotamer( curr_rot );
		iiratc.set_active_restype_index( restype_index_for_curr_rotamer );
		iiratc.update_from_momento( curr_state.momento_for_moltenres( ii ) );
	}

	/// 2.
	scminmap.clear_active_dofs();
	for ( Size ii = 1; ii <= rotsets.nmoltenres(); ++ii ) {
		Size iiresid = rotsets.moltenres_2_resid( ii );
		scminmap.activate_residue_dofs( iiresid );
		scminmap.set_natoms_for_residue( iiresid, atc->residue_atomtree_collection( iiresid ).active_restype().natoms() );
	}
	scminmap.setup( atc );

	/// 3.
	for ( Size ii = 1; ii <= rotsets.nmoltenres(); ++ii ) {
		Size iiresid = rotsets.moltenres_2_resid( ii );
		scmin::ResidueAtomTreeCollection & iiratc = atc->residue_atomtree_collection( iiresid );
		reinitialize_mingraph_neighborhood_for_residue( pose, sfxn, bgres, scminmap, iiratc.active_residue(), mingraph );
	}

	/// 4.
	//scmin::SCMinMultifunc scmin_func( pose, bgres, sfxn, mingraph, scminmap );
	optimization::MultifuncOP scmin_func = scminmap.make_multifunc( pose, bgres, sfxn, mingraph );
	utility::vector1< Real > allchi;
	scminmap.starting_dofs( allchi );


	return (*scmin_func)( allchi );
}

Real
minimize_alt_rotamer(
	pose::Pose & pose,
	utility::vector1< conformation::ResidueCOP > const & bgres,
	scoring::ScoreFunction const & sfxn,
	scoring::MinimizationGraph & mingraph,
	scmin::SCMinMinimizerMap & scminmap,
	optimization::MinimizerOptions const & min_options,
	scmin::SidechainStateAssignment const & curr_state,
	scmin::AtomTreeCollectionOP atc,
	rotamer_set::RotamerSets const & rotsets,
	Size resid,
	Size moltenres_id,
	Size rotamer_state_on_moltenres
#ifdef APL_FULL_DEBUG
	, pose::Pose & debug_pose
#endif

)
{
	/// 0. (to be made more efficient) get the starting energy for this residue
	/// 1. setup the atom tree collection for the (possibly) new residue type
	/// 2. setup the scminmap so that it has space to store derivatives for all the atoms in the molten residue
	/// 3. setup the mingraph
	/// 4. create the scminmultifunc
	/// 5. run the minimizer
	/// 6. re-evaluate the energy of the new minimized residue (made accurate by a nblist reconstruction)
	/// 7. return the delta-energy

	/// 0.
	Real start_energy_for_residue( 0.0 );
	if ( ! curr_state.any_unassigned() ) {
		start_energy_for_residue = get_residue_current_energy(
			pose, bgres, sfxn, mingraph, scminmap,
			curr_state, atc, rotsets, resid, moltenres_id
#ifdef APL_FULL_DEBUG
			, debug_pose
#endif
		);
	}

	/// 1.
	scmin::ResidueAtomTreeCollection & ratc = atc->residue_atomtree_collection( resid );
	Size const restype_index_for_rotamer = rotsets.rotamer_set_for_residue( resid )->get_residue_type_index_for_rotamer( rotamer_state_on_moltenres );
	ratc.set_active_restype_index( restype_index_for_rotamer );
	ratc.set_rescoords( * rotsets.rotamer_for_moltenres( moltenres_id, rotamer_state_on_moltenres ));
	ratc.update_atom_tree();

	/// 2.
	scminmap.set_natoms_for_residue( resid, atc->residue_atomtree_collection( resid ).active_restype().natoms() );

	if ( curr_state.any_unassigned() ) return 0.0; // quit now; can't minimize in the absence of assigned rotamers

	scminmap.clear_active_dofs();
	scminmap.activate_residue_dofs( resid );
	scminmap.setup( atc );

	/// 3.
	reinitialize_mingraph_neighborhood_for_residue( pose, sfxn, bgres, scminmap, ratc.active_residue(), mingraph );

	/// 4.
	//scmin::SCMinMultifunc scmin_func( pose, bgres, sfxn, mingraph, scminmap );
	optimization::MultifuncOP scmin_func = scminmap.make_multifunc( pose, bgres, sfxn, mingraph );
	//utility::vector1< Real > chi = ratc.active_residue().chi();
	utility::vector1< Real > chi;
	scminmap.starting_dofs( chi );

	/// 5.
	optimization::Minimizer minimizer( *scmin_func, min_options );
	minimizer.run( chi );

	/// 6.
	reinitialize_mingraph_neighborhood_for_residue( pose, sfxn, bgres, scminmap, ratc.active_residue(), mingraph );
	Real final_energy_for_residue = (*scmin_func)( chi );

	/// 7.
	//TR << "start energy for rotamer on " << resid << " " << start_energy_for_residue << " final energy: " << final_energy_for_residue << std::endl;
	return final_energy_for_residue - start_energy_for_residue;

}

bool
pass_metropolis(
	Real temperature,
	Real deltaE
)
{
	/// call this every time, for better numerical stability. Otherwise delta_energies ~ 0 can lead to
	/// instability in the number of calls to the random number generator
	core::PackerEnergy const rg_uniform( numeric::random::rg().uniform() );

	if ( deltaE < 0 ) {
		return true;
	} else { //evaluate prob of substitution
		Real lnprob = deltaE / temperature;
		if ( lnprob < 10.0 ) {
			Real probability = std::exp(-lnprob);
			if ( probability > rg_uniform ) return true;
		}
	}
	return false;
}

Size
n_inner_iterations( Real, Size nrotamers )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::packing;

	int scale = 5;
	if ( option[ minpack_inner_iteration_scale ].user() ) {
		scale = option[ minpack_inner_iteration_scale ]();
		if ( scale < 0 ) scale = 0; // negative iterations make no sense.
	}

	return scale * nrotamers;
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

	bool discrepancy( false );
	// AMW: cppcheck flags this as being reducible in scope, but since it's static
	// I think it should stay as is
	static int n_discreps( 0 );

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
			std::cout << std::endl;
			Real delta = emap.dot( pose.energies().weights() );
			if ( delta > 0 ) {
				discrepancy = true;
			}
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
					if ( kkst != hbond_bb_sc ) {
						discrepancy = true; // the hydrogen bond exclusion rule can lead to funniness.
					}
					if ( kkst == fa_atr || kkst == fa_rep || kkst == fa_sol ) {
						etab_discrepancy = true;
					}
				}
			}
		}
		if ( etab_discrepancy ) {
			std::cout << "etable discrepancy" << std::endl;
			///using namespace scoring;
			///using namespace scoring::etable;
			///methods::EnergyMethodOptions options; // default is fine
			//EtableEnergy etab_energy( *( ScoringManager::get_instance()->etable( options.etable_type() )), options );
			//ResiduePairNeighborList const & nblist( static_cast< ResiduePairNeighborList const & > (minedge->res_pair_min_data().get_data_ref( etab_pair_nblist )) );
			//debug_nblist_for_respair( pose.residue(ii), pose.residue(jj), pose, sfxn, etab_energy, nblist );
		}
	}

	if ( discrepancy ) {
		++n_discreps;
		pose.dump_pdb( "discrepancy_pose_" + utility::to_string( n_discreps ) + ".pdb" );
	}

}

void
assign_random_rotamers(
	pose::Pose & pose,
	utility::vector1< conformation::ResidueCOP > & bgres,
	scoring::ScoreFunction const & sfxn,
	scoring::MinimizationGraphOP mingraph,
	scmin::SCMinMinimizerMapOP scminmap,
	optimization::MinimizerOptions const & min_options,
	scmin::SidechainStateAssignment & curr_state,
	scmin::SidechainStateAssignment & best_state,
	scmin::AtomTreeCollectionOP atc,
	rotamer_set::RotamerSetsCOP rotsets
#ifdef APL_FULL_DEBUG
	, pose::Pose & debug_pose
#endif
)
{
	/// Assign random rotamers to start.
	/// This is the initial assignment; from here on out, we're doing regular sim annealing instead
	/// of waiting until each residue has an assigned rotamer.
	for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
		Size const ii_resid = rotsets->moltenres_2_resid( ii );
		Size const ii_rotamer_on_moltenres = static_cast< Size > ( rotsets->nrotamers_for_moltenres( ii ) * numeric::random::rg().uniform() + 1 );

		/// This is an aborted call; this function knows not to attempt to minimize until curr_state has
		/// no unassigned residues.
		minimize_alt_rotamer(
			pose, bgres, sfxn, *mingraph, *scminmap, min_options, curr_state,
			atc, *rotsets, ii_resid, ii, ii_rotamer_on_moltenres
#ifdef APL_FULL_DEBUG
			, debug_pose
#endif
		);
		atc->moltenres_atomtree_collection( ii ).save_momento( curr_state.state_momento( ii ) );
		curr_state.assign_state( ii, ii_rotamer_on_moltenres );
		bgres[ ii_resid ] = atc->moltenres_atomtree_collection( ii ).active_residue_cop();

#ifdef APL_FULL_DEBUG
      debug_pose.replace_residue( ii_resid, atc->moltenres_atomtree_collection( ii ).active_residue(), false );
#endif
	}

	Real totalE = get_total_energy_for_state( pose, bgres, sfxn, *mingraph, *scminmap, curr_state, atc, *rotsets );
	//std::cout << "Initial assignment total energy: " << totalE << std::endl;
#ifdef APL_FULL_DEBUG
	Real debug_pose_start_score( sfxn( debug_pose ) );
	if ( std::abs( debug_pose_start_score - totalE ) > 1e-3 ) {
		std::cout << "Initial energy assignment disagrees with actual initial energy " << totalE << " " << debug_pose_start_score << std::endl;
		debug_pose.energies().total_energies().show_weighted( std::cout, sfxn.weights() );
		std::cout << std::endl;
		for ( Size kk = 1; kk <= debug_pose.total_residue(); ++kk ) {
			compare_mingraph_and_energy_graph( kk, debug_pose, sfxn, *mingraph );
		}
	}
#endif

	curr_state.assign_energy( totalE );
	best_state = curr_state;

}

void
min_pack(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP input_task,
	bool cartesian,
	bool nonideal
)
{
	using namespace interaction_graph;
	using namespace rotamer_set;
	task::PackerTaskOP task = input_task->clone();

	Real start_score( sfxn( pose ) );

	/// 1. build the packer neighbor graph and construct rotamers
	/// 2. build the SCMinMinimizerMap & the MinimizationGraph
	/// 3. build the atom tree sets for each of the molten residues
	/// 4. assign a random rotamer to each molten residue
	/// 5. begin the temperature loop
	///     6. begin the fixed-temperature loop
	///          7. given a considered rotamer substitution, run a quick minimization
	///          8. if the new rotamer passes the metropolis criterion, save a momento of the network state
	/// 9. construct the pose from the momentos of the lowest scoring network state

	rotamer_set::RotamerSetsOP rotsets;
	scmin::SCMinMinimizerMapOP scminmap;
	scoring::MinimizationGraphOP mingraph;
	scmin::AtomTreeCollectionOP atc;
	optimization::MinimizerOptionsOP min_options;

	min_pack_setup(
		pose, sfxn, task, cartesian, nonideal,
		rotsets, scminmap, mingraph,
		atc, min_options );

	scmin::SidechainStateAssignment best_state( rotsets->nmoltenres() );

	min_pack_optimize(
		pose, sfxn, task, rotsets, scminmap,
		mingraph, atc, *min_options, best_state );


	min_pack_place_opt_rotamers_on_pose( pose, sfxn, rotsets, atc, best_state, start_score );
}

void
min_pack_setup(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	task::PackerTaskOP task,
	bool cartesian,
	bool nonideal,
	rotamer_set::RotamerSetsOP & rotsets,
	scmin::SCMinMinimizerMapOP & scminmap,
	scoring::MinimizationGraphOP & mingraph,
	scmin::AtomTreeCollectionOP & atc,
	optimization::MinimizerOptionsOP & min_options
)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::packing;

	//time_t starttime = time( NULL );

	/// 1.
	pack_scorefxn_pose_handshake( pose, sfxn);
	pose.update_residue_neighbors();
	sfxn.setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );
	if ( option[ minpack_disable_bumpcheck ] ) task->set_bump_check( false );
	graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, sfxn, task );
	rotsets = rotamer_set::RotamerSetsOP( new rotamer_set::RotamerSets() );
	rotsets->set_task( task );
	rotsets->build_rotamers( pose, sfxn, packer_neighbor_graph );
	rotsets->prepare_sets_for_packing( pose, sfxn );

	/// 2.
	scminmap = create_scmin_minimizer_map( pose, task, cartesian );
	scminmap->set_nonideal(nonideal);
	mingraph = create_minimization_graph( pose, sfxn, *task, *packer_neighbor_graph, *scminmap );

	/// 3.
	atc = scmin::AtomTreeCollectionOP( new scmin::AtomTreeCollection( pose, *rotsets ) );

	// true -- nblist, false -- deriv_check, false -- deriv_verbose
	//optimization::MinimizerOptions min_options( "dfpmin", 0.1, true, false, false );
	std::string minimizer = "dfpmin";
	Size max_iter=200;
	if ( cartesian || nonideal ) {
		if ( !sfxn.ready_for_nonideal_scoring() ) {
			utility_exit_with_message( "scorefunction not set up for nonideal/Cartesian scoring" );
		}
		minimizer = "lbfgs_armijo_atol";
		//max_iter = 25;                  // PTC - this doesn't give Cartesian enough time to converge
	}
	min_options = optimization::MinimizerOptionsOP( new optimization::MinimizerOptions( minimizer, 0.1, true, false, false ) );
	min_options->max_iter(max_iter);
	min_options->silent(true);

	//time_t stoptime = time( NULL );
	//std::cerr << "min_pack_setup took " << stoptime - starttime << std::endl;
}

void
min_pack_optimize(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	task::PackerTaskOP task,
	rotamer_set::RotamerSetsOP rotsets,
	scmin::SCMinMinimizerMapOP scminmap,
	scoring::MinimizationGraphOP mingraph,
	scmin::AtomTreeCollectionOP atc,
	optimization::MinimizerOptions const & min_options,
	scmin::SidechainStateAssignment & best_state
)
{
	//time_t starttime = time( NULL );

#ifdef APL_FULL_DEBUG
	pose::Pose debug_pose( pose );
	sfxn( debug_pose );
#endif

	/// 4.
	utility::vector1< conformation::ResidueCOP > bgres = setup_bgres_cops( pose, *task );
	scmin::SidechainStateAssignment curr_state( rotsets->nmoltenres() );
	assign_random_rotamers( pose, bgres, sfxn, mingraph,
		scminmap, min_options, curr_state, best_state, atc, rotsets
#ifdef APL_FULL_DEBUG
		, debug_pose
#endif
	);

	/// 5.
	utility::vector1< Real > temps = initialize_temperatures();

	for ( Size ii = 1; ii <= temps.size(); ++ii ) {
		Real ii_temperature = temps[ ii ];
		/// 6.
		for ( Size jj = 1, jj_end = n_inner_iterations( ii_temperature, rotsets->nrotamers() ); jj <= jj_end; ++jj ) {
			Size const jj_ranrotamer = static_cast< Size > ( rotsets->nrotamers() * numeric::random::rg().uniform() + 1 );
			Size const jj_moltenres_id = rotsets->moltenres_for_rotamer( jj_ranrotamer );
			Size const jj_resid = rotsets->moltenres_2_resid(  jj_moltenres_id );
			Size const jj_rotamer_state_on_moltenres = rotsets->rotid_on_moltenresidue( jj_ranrotamer );

#ifdef APL_FULL_DEBUG
			conformation::ResidueCOP before_sub_rotamer( 0 );
			if ( bgres[ jj_resid ] ) {
				before_sub_rotamer = new conformation::Residue( *bgres[ jj_resid ] );
			}
			if ( ! curr_state.any_unassigned() ) {
				for ( Size kk = 1; kk <= rotsets->nmoltenres(); ++kk ) {
					Size kkresid = rotsets->moltenres_2_resid( kk );
					for ( Size ll = 1; ll <= bgres[ kkresid ]->natoms(); ++ll ) {
						if ( std::abs( bgres[ kkresid ]->xyz( ll ).distance_squared( debug_pose.residue( kkresid ).xyz( ll ) )) > 1e-5 ) {
							std::cout << "  DEBUG error: debug pose and bgres disagree: " << kkresid << " " << bgres[ kkresid ]->atom_name( ll ) << std::endl;
						}
					}
				}
			}
#endif

			/// 7.
			Real deltaE = minimize_alt_rotamer(
				pose, bgres, sfxn, *mingraph, *scminmap, min_options, curr_state,
				atc, *rotsets, jj_resid, jj_moltenres_id, jj_rotamer_state_on_moltenres
#ifdef APL_FULL_DEBUG
				, debug_pose
#endif
			);

#ifdef APL_FULL_DEBUG

			Real debug_before = sfxn( debug_pose );
			scoring::EnergyMap total_before = debug_pose.energies().total_energies();
			debug_pose.replace_residue( jj_resid, atc->moltenres_atomtree_collection( jj_moltenres_id ).active_residue(), false );
			Real debug_after = sfxn( debug_pose );
			scoring::EnergyMap total_after = debug_pose.energies().total_energies();
			Real real_deltaE = debug_after - debug_before;
			if ( !curr_state.any_unassigned() && std::abs( real_deltaE - deltaE ) > 1e-3 ) {
				std::cout << "DeltaE mismatch: " << real_deltaE << " vs " << deltaE << std::endl;
				scoring::EnergyMap total_diff = total_after;
				for ( Size ii = 1; ii <= scoring::n_score_types; ++ii ) {
					total_diff[ (scoring::ScoreType) ii ] -= total_before[ (scoring::ScoreType) ii ];
				}
				std::cout << "Difference in energies: before vs after " << std::endl;
				total_diff.show_weighted( std::cout, sfxn.weights() );
				std::cout << std::endl;
				compare_mingraph_and_energy_graph( jj_resid, debug_pose, sfxn, *mingraph );
			} else if ( !curr_state.any_unassigned() ) {
				//std::cout << "Good deltaE" << std::endl;
			}
#endif

			/// 7.
			if ( curr_state.any_unassigned() || pass_metropolis( ii_temperature, deltaE ) ) {
				//std::cout << "Substitution accepted " << jj_resid << std::endl;
				// accept the rotamer substitution
				if ( curr_state.any_unassigned() ) {
					if ( curr_state.orig_rotamer_id_for_moltenres( jj_moltenres_id ) == 0 && curr_state.n_unassigned() == 1 ) {
						/// This is the last un-assigned residue that needs assigning
						atc->moltenres_atomtree_collection( jj_moltenres_id ).save_momento( curr_state.state_momento( jj_moltenres_id ) );
						curr_state.assign_state( jj_moltenres_id, jj_rotamer_state_on_moltenres );
						bgres[ jj_resid ] = atc->moltenres_atomtree_collection( jj_moltenres_id ).active_residue_cop();
						Real totalE = get_total_energy_for_state( pose, bgres, sfxn, *mingraph, *scminmap, curr_state, atc, *rotsets );
						std::cout << "Initial assignment total energy: " << totalE << std::endl;
#ifdef APL_FULL_DEBUG
						if ( std::abs( debug_after - totalE ) > 1e-3 ) {
							std::cout << "Initial energy assignment disagrees with actual initial energy" << std::endl;
							for ( Size kk = 1; kk <= debug_pose.total_residue(); ++kk ) {
								compare_mingraph_and_energy_graph( kk, debug_pose, sfxn, *mingraph );
							}
						}
#endif

						curr_state.assign_energy( totalE );
						best_state = curr_state;
					} else {
						atc->moltenres_atomtree_collection( jj_moltenres_id ).save_momento( curr_state.state_momento( jj_moltenres_id ) );
						curr_state.assign_state( jj_moltenres_id, jj_rotamer_state_on_moltenres );
						bgres[ jj_resid ] = atc->moltenres_atomtree_collection( jj_moltenres_id ).active_residue_cop();
					}
				} else {
					Real alt_totalE = curr_state.energy() + deltaE;
					curr_state.assign_energy( alt_totalE );
					//std::cout << "Accepted substitution total energy: " << alt_totalE << " temp: " << ii_temperature << std::endl;
#ifdef APL_FULL_DEBUG
					if ( std::abs( alt_totalE - debug_after ) > 1e-3 ) {
						std::cout << "total energy mismatch: " << debug_after << " " << alt_totalE << std::endl;
					}
#endif

					atc->moltenres_atomtree_collection( jj_moltenres_id ).save_momento( curr_state.state_momento( jj_moltenres_id ) );
					curr_state.assign_state( jj_moltenres_id, jj_rotamer_state_on_moltenres );
					bgres[ jj_resid ] = atc->moltenres_atomtree_collection( jj_moltenres_id ).active_residue_cop();
					if ( curr_state.energy() < best_state.energy() ) {
						Real totalE = get_total_energy_for_state( pose, bgres, sfxn, *mingraph, *scminmap, curr_state, atc, *rotsets );
#ifdef APL_FULL_DEBUG
						if ( std::abs( totalE - curr_state.energy() ) > 1e-6 ) {
							std::cout << "totalE and curr_state energy mismatch " << totalE << " " << curr_state.energy() << std::endl;
						}
#endif
						curr_state.assign_energy( totalE );
						if ( totalE < best_state.energy() ) {
							best_state = curr_state;
							//std::cout << "Accepted new best " << best_state.energy() << " " << totalE << " discrepancy: " << best_state.energy() - totalE << std::endl;
						}
					}
				}

			} else {
				// reset the bgres state
				//std::cout << "Substitution rejected " << jj_resid << std::endl;
				atc->moltenres_atomtree_collection( jj_moltenres_id ).
					update_from_momento( curr_state.momento_for_moltenres( jj_moltenres_id ) );
				bgres[ jj_resid ] = atc->moltenres_atomtree_collection( jj_moltenres_id ).active_residue_cop();

				reinitialize_mingraph_neighborhood_for_residue( pose, sfxn, bgres, *scminmap, *bgres[ jj_resid ], *mingraph );

#ifdef APL_FULL_DEBUG
				debug_pose.replace_residue( jj_resid, *bgres[ jj_resid ], false );
				sfxn( debug_pose );
				if ( before_sub_rotamer && before_sub_rotamer->natoms() == bgres[ jj_resid ]->natoms() ) {
					for ( Size kk = 1; kk <= before_sub_rotamer->natoms(); ++kk ) {
						if ( before_sub_rotamer->xyz( kk ).distance_squared( bgres[ jj_resid ]->xyz( kk ) ) > 1e-5 ) {
							std::cout << "Failed to properly restore the previous state after rejecting substitution at residue " << jj_resid << " atom " << kk << std::endl;
						}
					}
				} else {
					std::cout << "NATOMS discrepancy in before-sub rotamer with bgrotamer" << std::endl;
				}
#endif

			}
		}
		//std::cout << "Finished temperature " << ii_temperature << " with energy " << curr_state.energy() << " and best energy " << best_state.energy() << std::endl;
	}
	//time_t stoptime = time( NULL );
	//std::cerr << "min_pack_optimize took " << stoptime - starttime << std::endl;

}

void
min_pack_place_opt_rotamers_on_pose(
	core::pose::Pose & pose,
	core::scoring::ScoreFunction const & sfxn,
	rotamer_set::RotamerSetsOP rotsets,
	scmin::AtomTreeCollectionOP atc,
	scmin::SidechainStateAssignment const & best_state,
	Real start_score
)
{
	/// 8.
	for ( Size ii = 1; ii <= rotsets->nmoltenres(); ++ii ) {
		Size iiresid = rotsets->moltenres_2_resid( ii );
		atc->moltenres_atomtree_collection( ii ).update_from_momento( best_state.momento_for_moltenres( ii ) );
		pose.replace_residue( iiresid, atc->moltenres_atomtree_collection( ii ).active_residue(), false );
	}

	Real final_score( sfxn( pose ) );
	TR << "min pack final score: " << final_score << " start_score: " << start_score << std::endl;

}

void compare_simple_inteaction_graph_alt_state_and_energy_graph(
	Size resid,
	pose::Pose const & pose,
	scoring::ScoreFunction const & sfxn,
	interaction_graph::SimpleInteractionGraph const & ig
)
{
	using namespace pose;
	using namespace scoring;
	using namespace graph;
	using namespace interaction_graph;

	bool discrepancy( false );
	// AMW: cppcheck flags this as being reducible in scope, but since it's static
	// I think it should stay as-is
	static int n_discreps( 0 );

	EnergyMap const & one_body_emap( pose.energies().onebody_energies( resid ));
	Real onebody_energy = sfxn.weights().dot( one_body_emap );
	Real sig_onebody_energy = ig.get_simple_node( resid )->proposed_one_body_energy();

	if ( std::abs( onebody_energy - sig_onebody_energy ) > 1e-5 ) {
		std::cout << " onebody energy discrepancy: " << onebody_energy << " " << sig_onebody_energy << " " << onebody_energy - sig_onebody_energy << std::endl;
	}

	/*for ( Size kk = 1; kk <= total_score; ++kk ) {
	ScoreType kkst = ScoreType(kk);
	if ( sfxn.weights()[ kkst ] != 0.0 ) {
	if ( std::abs( one_body_emap[ kkst ] - min_node_1b[ kkst ] ) > 1e-10 ) {
	std::cout << "   one body discrepancy " << kkst << ": " << one_body_emap[ kkst ] << " " << min_node_1b[ kkst ] << std::endl;
	}
	}
	}*/

	EnergyGraph const & eg( pose.energies().energy_graph() );
	for ( Node::EdgeListConstIter iter = eg.get_node(resid)->const_edge_list_begin(),
			iter_end = eg.get_node(resid)->const_edge_list_end();
			iter != iter_end; ++iter ) {
		Size ii( (*iter)->get_first_node_ind() );
		Size jj( (*iter)->get_second_node_ind() );
		if ( ii != resid && jj != resid ) continue; // only compare nodes that are involved in the current optimization
		SimpleEdge const * simple_edge = static_cast< SimpleEdge const * > ( ig.find_edge( ii, jj ) );
		EnergyEdge const * eedge = static_cast< EnergyEdge const * > ( *iter );
		EnergyMap emap = eedge->fill_energy_map();
		Real edge_energy = sfxn.weights().dot( emap );
		if ( ! simple_edge ) {
			std::cout << "Simple edge " << ii << " " << jj << " missing from simple interaction graph" << std::endl;
			emap.show_if_nonzero_weight( std::cout, pose.energies().weights() );
			std::cout << std::endl;
			Real delta = emap.dot( pose.energies().weights() );
			if ( delta > 0 ) {
				discrepancy = true;
			}
			continue;
		} else {
			Real ig_edge_energy = simple_edge->get_proposed_energy();
			if ( std::abs( edge_energy - ig_edge_energy ) > 1e-5 ) {
				std::cout << " twobody energy discrepancy: " << ii << " " << jj << " " << edge_energy << " " << ig_edge_energy << " " << edge_energy - ig_edge_energy << std::endl;
			}
		}
	}

	if ( discrepancy ) {
		++n_discreps;
		pose.dump_pdb( "discrepancy_pose_" + utility::to_string( n_discreps ) + ".pdb" );
	}

}

/// @details Returns the rotamer index (1 for the current rotamer) for a randomly assigned rotamer
Size
assign_random_continuous_rotamer(
	rotamer_set::ContinuousRotamerSet const & rotset,
	scmin::ResidueAtomTreeCollection & resatc,
	Size ranrot_on_moltenres
)
{
	Size const ran_rotblock_for_ranrot( rotset.get_rotblock_index_for_sampling_rotamer( ranrot_on_moltenres ) );
	resatc.set_active_restype_index( ran_rotblock_for_ranrot );
	if ( ranrot_on_moltenres == rotset.sampling_id_for_current_rotamer() ) {
		resatc.set_rescoords( rotset.current_rotamer_coords() );
		//return set_current_rotamer_for_moltres( ran_moltres, rotsets, atc );
		return 1 + rotset.get_n_baserotamers_for_rotblock( ran_rotblock_for_ranrot );
	} else {
		/// resatc.idealize_active_restype(); -- need this if we ever allow the input sidechain!
		if ( rotset.get_n_baserotamers_for_rotblock( ran_rotblock_for_ranrot ) != 0 ) {
			// i.e. not gly or ala.
			Real rand_btw_0_and_1 = numeric::random::rg().uniform();
			Size const ran_baserot_ind = rotset.pick_baserotamer_from_rotblock( ran_rotblock_for_ranrot, rand_btw_0_and_1 );
			dunbrack::DunbrackRotamerSampleData const & rotdata =
				rotset.baserotamer_data( ran_rotblock_for_ranrot, ran_baserot_ind );
			/// OK: now, sample the dunbrack chi
			for ( Size ii = 1; ii <= rotdata.nchi(); ++ii ) {
				/// First attempt strategy: sample uniformly +/- 1 standard deviation
				Real const iichi_randval = numeric::random::rg().uniform();
				Real const iisd = rotdata.chi_sd()[ ii ];
				resatc.set_chi( ii, rotdata.chi_mean()[ii] - iisd + 2 * iisd * iichi_randval );
			}
		}
		chemical::ResidueTypeCOP ran_restype( rotset.restype_for_rotblock( ran_rotblock_for_ranrot ));
		for ( Size ii = 1; ii  <= ran_restype->n_proton_chi(); ++ii ) {
			Size ii_chi = ran_restype->proton_chi_2_chi( ii );
			utility::vector1< Real > const & ii_chi_samples = ran_restype->proton_chi_samples( ii );
			utility::vector1< Real > const & ii_chi_extra_samples = ran_restype->proton_chi_extra_samples( ii );

			Real iichi_sdev;
			if ( ii_chi_extra_samples.size() == 0 ) {
				iichi_sdev = 10;
			} else {
				// assume extra samples has its largest value in the last position
				iichi_sdev = ii_chi_extra_samples[ ii_chi_extra_samples.size() ];
			}

			Real rand_btw_0_and_1 = numeric::random::rg().uniform();
			Size sample = static_cast< Size > ( ii_chi_samples.size() * rand_btw_0_and_1 ) + 1;
			Real sample_dev = numeric::random::rg().uniform();
			resatc.set_chi( ii_chi, ii_chi_samples[ sample ] - iichi_sdev + 2 * sample_dev * iichi_sdev );
		}
		resatc.update_residue();
		return ran_rotblock_for_ranrot;
	}
}


void
off_rotamer_pack_setup(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP task,
	rotamer_set::ContinuousRotamerSetsOP & rotsets,
	scmin::AtomTreeCollectionOP & atc,
	interaction_graph::SimpleInteractionGraphOP & ig
)
{
	// 1 (misc setup)
	pack_scorefxn_pose_handshake( pose, sfxn);
	sfxn.setup_for_packing( pose, task->repacking_residues(), task->designing_residues() );

	// 1a
	atc = scmin::AtomTreeCollectionOP( new scmin::AtomTreeCollection( pose, *task ) );

	// 1b
	rotsets = rotamer_set::ContinuousRotamerSetsOP( new rotamer_set::ContinuousRotamerSets( pose, *task ) );

	// 1c
	ig = interaction_graph::SimpleInteractionGraphOP( new interaction_graph::SimpleInteractionGraph );
	ig->set_scorefunction( sfxn );
	ig->set_pose_no_initialize( pose );
	graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, sfxn, task );
	//for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
	// ig.get_simple_node( ii )->set_current( pose.residue( ii ).clone() );
	//}
	ig->copy_connectivity( *packer_neighbor_graph );
}

void
off_rotamer_pack_optimize(
	rotamer_set::ContinuousRotamerSets const & rotsets,
	scmin::AtomTreeCollectionOP atc,
	interaction_graph::SimpleInteractionGraph & ig,
	scmin::SidechainStateAssignment & best_state
)
{
#ifdef APL_FULL_DEBUG
	pose::Pose debug_pose( pose );
	sfxn( debug_pose );
#endif


	/// 2.
	Size const n_sample_rots = rotsets.n_sample_rotamers();
	utility::vector1< Real > temps = initialize_temperatures_off_rotamer_pack();
	scmin::SidechainStateAssignment curr_state( rotsets.nmoltenres() );

	for ( Size ii = 1; ii <= temps.size(); ++ii ) {
		Real ii_temperature = temps[ ii ];
		/// 6.
		Size naccepts( 0 );
		Real accum_deltaE( 0 );
		Size const jj_end = n_inner_iterations( ii_temperature, n_sample_rots );
		for ( Size jj = 1; jj <= jj_end; ++jj ) {

			Size const random_sample_rot = static_cast< Size > ( n_sample_rots * numeric::random::rg().uniform() + 1 );
			Size const ran_moltres = rotsets.moltenres_for_sample_rot( random_sample_rot );
			Size const ran_res = rotsets.moltenresid_2_resid( ran_moltres );
			Size const ranrot_on_moltenres = rotsets.full_sample_rot_index_2_moltenres_sample_rot_index( random_sample_rot );

			Size const ranrot_baserot_id = assign_random_continuous_rotamer(
				rotsets.rotamer_set_for_moltenres( ran_moltres ),
				atc->moltenres_atomtree_collection( ran_moltres),
				ranrot_on_moltenres );

			Real const neg_deltaE = ig.consider_substitution( ran_res, atc->moltenres_atomtree_collection( ran_moltres ).active_residue_cop() );
			Real const deltaE = -1 * neg_deltaE;

#ifdef APL_FULL_DEBUG
			Real curE = sfxn( debug_pose );
			debug_pose.replace_residue( ran_res, atc->moltenres_atomtree_collection( ran_moltres ).active_residue(), false );
			Real altE = sfxn( debug_pose );
			Real actual_deltaE = altE - curE;
			if ( ! curr_state.any_unassigned() && std::abs( actual_deltaE - deltaE ) > 1e-5 ) {
				std::cout << "DeltaE discrepancy replacing residue " << ran_res << " "
					<< actual_deltaE << " " << deltaE << " " << actual_deltaE - deltaE << std::endl;
				compare_simple_inteaction_graph_alt_state_and_energy_graph( ran_res, debug_pose, sfxn, ig );
			}
#endif


			if ( curr_state.any_unassigned() || pass_metropolis( ii_temperature, deltaE ) ) {
				/// Accept the substitution
				++naccepts;
				accum_deltaE += deltaE;
				ig.commit_change( ran_res );

				/// Before assigning the state, write down: do we have extra bookkeeping we need to worry about?
				bool const any_previously_unassigned = curr_state.any_unassigned();
				bool const last_unassigned_assigned = any_previously_unassigned &&
					curr_state.orig_rotamer_id_for_moltenres( ran_moltres ) == 0 &&
					curr_state.n_unassigned() == 1;

				/// always record the new state
				atc->moltenres_atomtree_collection( ran_moltres ).save_momento( curr_state.state_momento( ran_moltres ) );
				curr_state.assign_state( ran_moltres, ranrot_baserot_id );

				if ( any_previously_unassigned ) {
					if ( last_unassigned_assigned ) {
						/// This is the last un-assigned residue that needs assigning; save the total energy for current & best
						Real totalE = ig.total_energy();;
						curr_state.assign_energy( totalE );
						best_state = curr_state;
					}
				} else {
					Real alt_totalE = curr_state.energy() + deltaE;
					curr_state.assign_energy( alt_totalE );

					if ( curr_state.energy() < best_state.energy() ) {
						Real totalE = ig.total_energy(); //get_total_energy_for_state( pose, bgres, sfxn, *mingraph, *scminmap, curr_state, atc, *rotsets );
						debug_assert( std::abs( totalE - curr_state.energy() ) < 1e-5 ); // drift does accumulate, but it should be small!
						curr_state.assign_energy( totalE );
						if ( totalE < best_state.energy() ) {
							best_state = curr_state;
							//std::cout << "Accepted new best " << best_state.energy() << " " << totalE << " discrepancy: " << best_state.energy() - totalE << std::endl;
						}
					}

				}
				//Real totalE = ig.total_energy(); //get_total_energy_for_state( pose, bgres, sfxn, *mingraph, *scminmap, curr_state, atc, *rotsets );
				//if( ! curr_state.any_unassigned() && std::abs( totalE - curr_state.energy() ) > 1e-5 ) {
				// std::cout << "Disagreement between totalE and curr_state.energy() " << totalE << " " << curr_state.energy() << " diff " << totalE - curr_state.energy() << std::endl;
				//}
			} else {
				// reject; restore the old residue
				atc->moltenres_atomtree_collection( ran_moltres ).update_from_momento( curr_state.momento_for_moltenres( ran_moltres ) );
				ig.reject_change( ran_res );
#ifdef APL_FULL_DEBUG
				debug_pose.replace_residue( ran_res, atc->moltenres_atomtree_collection( ran_moltres ).active_residue(), false );
#endif
			}
		}
		//std::cout << "Finished temperature " << ii_temperature << " with energy " << curr_state.energy() << " and best energy " << best_state.energy()
		// << " accept rate: " << ((double) naccepts )/ jj_end << " avg deltaE: " << accum_deltaE / ( naccepts == 0 ? 1 : naccepts ) << std::endl;

	}
}

void
off_rotamer_pack_update_pose(
	pose::Pose & pose,
	rotamer_set::ContinuousRotamerSets const & rotsets,
	scmin::AtomTreeCollectionOP atc,
	scmin::SidechainStateAssignment const & best_state
)
{
	/// 8.
	for ( Size ii = 1; ii <= rotsets.nmoltenres(); ++ii ) {
		Size iiresid = rotsets.moltenresid_2_resid( ii );
		atc->moltenres_atomtree_collection( ii ).update_from_momento( best_state.momento_for_moltenres( ii ) );
		pose.replace_residue( iiresid, atc->moltenres_atomtree_collection( ii ).active_residue(), false );
	}

}

void
off_rotamer_pack(
	pose::Pose & pose,
	scoring::ScoreFunction const & sfxn,
	task::PackerTaskCOP task
)
{
	/// 1. Create the classes that will be used in the stochastic packer
	////   (misc setup)
	///    a. The AtomTreeCollection
	///    b. The ContinuousRotamerSet
	///    c. The Packer Neighbor Graph / SimpleIG
	/// 2. Begin annealing
	///    a. Assign random sequence structure
	///    b. loop over temperatures
	///    c.    loop over fixed temperature n iterations
	///    d.        select residue to undergo substitution
	///    e.        select residue type on that residue
	///    f.        select random chi for that restype
	///    g.        set coords for that residue
	///    h.        compute new energy -> deltaE
	///    i.        metropolis accept/reject
	///    j.           restore previous state, or save new state

	Real start_score( sfxn( pose ) );

	rotamer_set::ContinuousRotamerSetsOP rotsets;
	scmin::AtomTreeCollectionOP atc;
	interaction_graph::SimpleInteractionGraphOP ig;

	off_rotamer_pack_setup( pose, sfxn, task, rotsets, atc, ig );

	scmin::SidechainStateAssignment best_state( rotsets->nmoltenres() );
	off_rotamer_pack_optimize( *rotsets, atc, *ig, best_state );

	off_rotamer_pack_update_pose( pose, *rotsets, atc, best_state );

	Real final_score( sfxn( pose ) );
	TR << "off-rotamer pack final score: " << final_score << " start_score: " << start_score << std::endl;

}


} // namespace pack
} // namespace core

