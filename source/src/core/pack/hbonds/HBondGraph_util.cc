// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/hbonds/HBondGraph_util.cc
/// @brief A collections of methods that are useful for dealing with HBondGraphs
/// @author Jack Maguire


#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/pack/hbonds/HBondGraph_util.hh>
#include <core/pack/hbonds/HBondGraphInitializerIG.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>

#include <core/scoring/Energies.hh>

#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/graph/HBondGraph.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/thread_manager/RosettaThreadAssignmentInfo.hh>
#include <basic/thread_manager/RosettaThreadManager.hh>

#include <list>

#include <utility/graph/Graph.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

#include <core/pack/task/PackerTask.hh> // AUTO IWYU For PackerTask

namespace core {
namespace pack {
namespace hbonds {

using namespace scoring::hbonds::graph;

scoring::hbonds::graph::HBondGraphOP create_init_and_create_edges_for_hbond_graph(
	rotamer_set::RotamerSetsOP rotamer_sets,
	scoring::ScoreFunction const & sfxn,
	pose::Pose const & pose,
	Real hydrogen_bond_threshold,
	Real clash_threshold,
	Real hbond_energy_threshold_for_satisfaction /*= -0.25f */,
	bool include_backbone_at_atom_level /*= false*/
){
	scoring::hbonds::graph::HBondGraphOP hbond_graph = create_and_init_hbond_graph( * rotamer_sets );

	HBondGraphInitializerIGOP ig = utility::pointer::make_shared< HBondGraphInitializerIG >( hbond_graph, rotamer_sets, hydrogen_bond_threshold, clash_threshold );
	ig->initialize( *rotamer_sets );

	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, sfxn, rotamer_sets->task() );

	//The following lines precompute the twobody interaction energies, using threads if available:
	//An object to store information about the actual threads that got assigned to the work we'll do.  (Actual threads can be less than requested threads.):
	basic::thread_manager::RosettaThreadAssignmentInfo thread_assignment_info( basic::thread_manager::RosettaThreadRequestOriginatingLevel::CORE_PACK );
	utility::vector1< basic::thread_manager::RosettaThreadFunction > work_vector; //Allocate space for the list of work to be done.
	rotamer_sets->append_two_body_energy_computations_to_work_vector( pose, sfxn, packer_neighbor_graph, ig, work_vector, thread_assignment_info ); //Make a list of work to be done.
	basic::thread_manager::RosettaThreadManager::get_instance()->do_work_vector_in_threads( work_vector, rotamer_sets->task()->ig_threads_to_request(), thread_assignment_info ); //Do the work.
	ig->declare_all_edge_energies_final(); //In a single thread, finalize the edges (not threadsafe, but happens in O(Nedge) time).
	ig->finalize_hbond_graph();

	//If you are running with symmetry, you still need to score one-body interactions

	determine_atom_level_edge_info_for_all_edges(
		* hbond_graph,
		* rotamer_sets,
		* scoring::hbonds::HBondDatabase::get_database(),
		pose.energies().tenA_neighbor_graph(),
		pose,
		false,
		hbond_energy_threshold_for_satisfaction
	);

	utility::vector1< bool > include_all_resids( pose.size(), true );
	determine_atom_level_node_info_for_all_nodes(
		* hbond_graph,
		* rotamer_sets,
		include_all_resids,
		false,
		include_backbone_at_atom_level
	);

	return hbond_graph;
}


void init_node_info( HBondGraph & graph, rotamer_set::RotamerSets const & rotamer_sets ){
	Size const nrot = rotamer_sets.nrotamers();

	for ( Size rot = 1; rot <= nrot; ++rot ) {
		HBondNode * node = graph.get_node( rot );
		Size const mres = rotamer_sets.moltenres_for_rotamer( rot );
		debug_assert( mres );
		node->set_moltenres( MResIDSize( mres ) );
		node->set_local_rotamer_id( RotamerIDSize( rotamer_sets.rotid_on_moltenresidue( rot ) ) );
	}//for rot
}


void find_hbonds_in_residue_pair(
	conformation::Residue const & resA,
	conformation::Residue const & resB,
	scoring::hbonds::HBondDatabase const & database,
	utility::graph::Graph const & tenA_neighbor_graph,
	scoring::hbonds::HBondSet & set
){
	unsigned int const residA = resA.seqpos();
	unsigned short int const num_nbrsA =
		tenA_neighbor_graph.get_node( residA )->num_neighbors_counting_self();
	unsigned short int const num_nbrsB =
		tenA_neighbor_graph.get_node( resB.seqpos() )->num_neighbors_counting_self();

	scoring::hbonds::identify_hbonds_1way(
		database,
		resA,
		resB,
		num_nbrsA,
		num_nbrsB,
		false, //bool const evaluate_derivative,
		false, //bool const exclude_don_bb,
		false, //bool const exclude_don_bsc,
		false, //bool const exclude_acc_scb,
		false, //bool const exclude_acc_sc,
		// output
		set //HBondSet & hbond_set,
		//Real ssdep_weight_factor = 1.0,
		//bool bond_near_wat = false
	);

	scoring::hbonds::identify_hbonds_1way(
		database,
		resB,
		resA,
		num_nbrsB,
		num_nbrsA,
		false, //bool const evaluate_derivative,
		false, //bool const exclude_don_bb,
		false, //bool const exclude_don_bsc,
		false, //bool const exclude_acc_scb,
		false, //bool const exclude_acc_sc,
		// output
		set //HBondSet & hbond_set,
		//Real ssdep_weight_factor = 1.0,
		//bool bond_near_wat = false
	);

}

void determine_atom_level_edge_info_for_all_edges(
	scoring::hbonds::graph::HBondGraph & hb_graph,
	rotamer_set::RotamerSets const & rotamer_sets,
	scoring::hbonds::HBondDatabase const & database,
	utility::graph::Graph const & tenA_neighbor_graph,
	pose::Pose const & pose,
	bool skip_edges_with_degree_zero,
	Real hbond_energy_threshold_for_satisfaction,
	conformation::symmetry::SymmetryInfoCOP symm_info
){
	for ( utility::graph::LowMemEdgeListIter it = hb_graph.edge_list_begin();
			it != hb_graph.edge_list_end(); ++it ) {
		HBondEdge & edge = static_cast< HBondEdge & >( * ( * it ) );

		if ( skip_edges_with_degree_zero ) {
			if ( hb_graph.get_node( edge.get_first_node_ind() )->num_edges() == 1 &&
					hb_graph.get_node( edge.get_second_node_ind() )->num_edges() == 1 ) {
				continue;
			}
		}

		determine_atom_level_edge_info(
			edge,
			rotamer_sets,
			database,
			tenA_neighbor_graph,
			pose,
			hbond_energy_threshold_for_satisfaction,
			symm_info
		);
	}
}

Size get_symm_ind_res(
	pose::Pose const & pose,
	Size const resid,
	conformation::symmetry::SymmetryInfoCOP symm_info,
	std::map< char, std::pair< Size, Size > > chain_bounds
) {
	Size resi_ind( resid );
	if ( resid > symm_info->num_independent_residues() ) {
		char resi_chain = pose.chain( resid );
		if ( symm_info->get_num_components() > 1 ) {
			std::map< std::string, std::pair< Size, Size > > const & component_bounds = symm_info->get_component_bounds();
			std::string resi_comp = symm_info->get_component_of_residue( resid );
			resi_ind = resid - chain_bounds[ resi_chain ].first + component_bounds.at( resi_comp ).first;
		} else {
			resi_ind = resid - chain_bounds[ resi_chain ].first + 1;
		}
	}
	return resi_ind;
}


void determine_atom_level_edge_info(
	HBondEdge & hb_edge,
	rotamer_set::RotamerSets const & rotamer_sets,
	scoring::hbonds::HBondDatabase const & database,
	utility::graph::Graph const & tenA_neighbor_graph,
	pose::Pose const & pose,
	Real hbond_energy_threshold_for_satisfaction,
	conformation::symmetry::SymmetryInfoCOP symm_info
){
	scoring::hbonds::HBondSet set;

	conformation::Residue const & resA = * rotamer_sets.rotamer( hb_edge.get_first_node_ind() );
	conformation::Residue const & resB = * rotamer_sets.rotamer( hb_edge.get_second_node_ind() );
	Size const residA = resA.seqpos();
	Size const residB = resB.seqpos();

	if ( symm_info ) {

		//Stolen from HBNet.cc, authored by Scott Boyken
		std::map< char, std::pair< Size, Size > > chain_bounds;
		for ( Size ic = 1; ic <= pose.conformation().num_chains(); ++ic ) {
			Size ic_begin = pose.conformation().chain_begin( ic );
			Size ic_end = pose.conformation().chain_end( ic );
			char chain = pose.chain( ic_begin );
			chain_bounds[ chain ].first = ic_begin;
			chain_bounds[ chain ].second = ic_end;
		}

		Size const num_symm_subunits = symm_info->subunits();
		utility::vector1< Size > resids_for_A_clones( num_symm_subunits );
		utility::vector1< Size > resids_for_B_clones( num_symm_subunits );

		for ( Size subunit = 1; subunit <= num_symm_subunits; ++subunit ) {
			resids_for_A_clones[ subunit ] = symm_info->equivalent_residue_on_subunit( subunit, residA );
			resids_for_B_clones[ subunit ] = symm_info->equivalent_residue_on_subunit( subunit, residB );
		}

		if ( residA != residB ) {
			find_hbonds_in_residue_pair( resA, resB, database, tenA_neighbor_graph, set );
		}

		//We do not need to do an all-to-all comparison. We just need to iterate over one set of symmetric resids. If A is not symmetric, iterate over B
		if ( ! resids_for_A_clones.empty() ) {
			rotamer_set::symmetry::SymmetricRotamerSet_COP sym_set =
				utility::pointer::dynamic_pointer_cast< rotamer_set::symmetry::SymmetricRotamerSet_ const > ( rotamer_sets.rotamer_set_for_residue( residA ) );
			runtime_assert( sym_set );
			for ( Size const sympos : resids_for_A_clones ) {
				conformation::ResidueOP resA_ii = sym_set->orient_rotamer_to_symmetric_partner( pose, resA, sympos );
				find_hbonds_in_residue_pair( * resA_ii, resB, database, tenA_neighbor_graph, set );
			}
		} else if ( ! resids_for_B_clones.empty() ) {
			rotamer_set::symmetry::SymmetricRotamerSet_COP sym_set =
				utility::pointer::dynamic_pointer_cast< rotamer_set::symmetry::SymmetricRotamerSet_ const > ( rotamer_sets.rotamer_set_for_residue( residB ) );
			runtime_assert( sym_set );
			for ( Size sympos : resids_for_B_clones ) {
				conformation::ResidueOP resB_ii = sym_set-> orient_rotamer_to_symmetric_partner( pose, resB, sympos );
				find_hbonds_in_residue_pair( resA, * resB_ii, database, tenA_neighbor_graph, set );
			}
		}

	} else {
		find_hbonds_in_residue_pair( resA, resB, database, tenA_neighbor_graph, set );
	}

	unsigned short int const num_hbonds = set.nhbonds();
	for ( unsigned short int ii = 1; ii <= num_hbonds; ++ii ) {
		scoring::hbonds::HBond const & hbond = set.hbond( ii );
		if ( hbond.energy() > hbond_energy_threshold_for_satisfaction ) continue;

		bool const first_node_is_donor = hbond.don_res() == residA;
		unsigned short int const Hatm = hbond.don_hatm();
		unsigned short int const Datm = ( first_node_is_donor ? resA.atom_base( Hatm ) : resB.atom_base( Hatm ) );
		unsigned short int const Aatm = hbond.acc_atm();

		if ( false ) { //TODO bridging_waters
			//LKHBondEdge * lk_edge = static_cast< LKHBondEdge * >( lk_edge );
			//LKHBondInfo info = ...
			//lk_edge->register_hbond( info, first_node_is_donor, Aatm, Datm, Hatm );
		} else {
			hb_edge.register_hbond( first_node_is_donor, Aatm, Datm, Hatm, hbond.energy() );
		}
	}// for all hbonds
}

void determine_atom_level_node_info_for_all_nodes(
	scoring::hbonds::graph::HBondGraph & hb_graph,
	rotamer_set::RotamerSets const & rotamer_sets,
	utility::vector1< bool > const & include_these_resids,
	bool skip_nodes_with_no_edges,
	bool include_backbone /*= false*/
){
	Size const num_nodes = hb_graph.num_nodes();
	for ( unsigned int ii=1; ii <= num_nodes; ++ii ) {
		auto * al_node = static_cast< HBondNode * >( hb_graph.get_node( ii ) );
		if ( skip_nodes_with_no_edges && al_node->num_edges() == 0 ) {
			continue;
		}

		determine_atom_level_node_info( * al_node, rotamer_sets, include_these_resids, include_backbone );
	}
}

void determine_atom_level_node_info(
	HBondNode & al_node,
	rotamer_set::RotamerSets const & rotamer_sets,
	utility::vector1< bool > const & include_these_resids,
	bool include_backbone /*= false*/
){
	Size const rot = al_node.global_rotamer_id();
	Size const mres = rotamer_sets.moltenres_for_rotamer( rot );

	conformation::ResidueCOP rotamer = rotamer_sets.rotamer( rot );
	debug_assert( rotamer->seqpos() == rotamer_sets.moltenres_2_resid( mres ) );
	if ( ! include_these_resids[ rotamer->seqpos() ] ) return;

	////////
	//HEAVY
	unsigned short int const nheavyatoms = rotamer->nheavyatoms();
	unsigned short int const first_sidechain_atom = rotamer->first_sidechain_atom();
	for ( unsigned short int ii = ( include_backbone ? 1 : first_sidechain_atom );
			ii <= nheavyatoms;
			++ii ) {

		bool const is_don = rotamer->heavyatom_is_an_acceptor( ii );
		bool const is_acc = rotamer->heavyatom_has_polar_hydrogens( ii );
		if ( is_don || is_acc ) {
			al_node.add_polar_atom(
				ii,
				rotamer->xyz( ii ),
				false, //is_hydrogen
				is_don,
				is_acc,
				rotamer->atom_type( ii ).name() == "OH", //is_hydroxyl
				ii < first_sidechain_atom // is_backbone
			);
		}
	}//for heavy atoms

	//////
	//HPOL
	unsigned short int const natoms = rotamer->natoms();
	for ( unsigned short int jj = nheavyatoms + 1; jj <= natoms; ++jj ) {
		if ( rotamer->atom_is_polar_hydrogen( jj ) ) {
			al_node.add_polar_atom(
				jj,
				rotamer->xyz( jj ),
				true, //is_hydrogen
				false,
				false,
				rotamer->atom_type( rotamer->atom_base( jj ) ).name() == "OH", //is_hydroxyl
				rotamer->atom_is_backbone( jj )
			);
		}
	}

	if ( false ) { //TODO bridging waters
		//LKHBondNode * lkal_node = static_cast< LKHBondNode * >( al_node );
	}
}



void find_satisfying_interactions_with_background(
	HBondNode & node,
	rotamer_set::RotamerSets const & rotamer_sets,
	utility::graph::Graph const & packer_neighbor_graph,
	pose::Pose const & poly_ala_pose,
	Real hbond_energy_threshold_for_satisfaction
) {

	conformation::ResidueCOP const rotamer = rotamer_sets.rotamer( node.global_rotamer_id() );
	unsigned int const resid = rotamer->seqpos();
	debug_assert( resid == rotamer_sets.moltenres_2_resid( node.moltenres() ) );

	utility::graph::Node const * const resid_nbr_node = packer_neighbor_graph.get_node( resid );
	scoring::hbonds::HBondSet set;
	auto const & database = * scoring::hbonds::HBondDatabase::get_database( set.hbond_options().params_database_tag() );

	//////////////////
	// FIND ALL HBONDS
	for ( utility::graph::EdgeListConstIterator it = resid_nbr_node->const_edge_list_begin(),
			end = resid_nbr_node->const_edge_list_end();
			it != end;
			++it ) {

		unsigned int const other_resid = (*it)->get_other_ind( resid );
		conformation::Residue const & other_residue = poly_ala_pose.residue( other_resid );
		find_hbonds_in_residue_pair( *rotamer, other_residue, database, poly_ala_pose.energies().tenA_neighbor_graph(), set );

	}//for all neighbors


	///////////////////////////////////////
	// REMOVE ATOM INFO FOR SATISFIED ATOMS
	//
	// For every rotamer in the hbond graph,
	// determine which atoms are hydrogen bonding
	// with the background and remove them from
	// the container of unsatisfied atoms
	unsigned short int const num_hbonds = set.nhbonds();
	for ( unsigned short int ii = 1; ii <= num_hbonds; ++ii ) {
		scoring::hbonds::HBond const & hbond = set.hbond( ii );
		if ( hbond.energy() > hbond_energy_threshold_for_satisfaction ) continue;

		if ( hbond.don_res() == resid ) {
			unsigned short int const Hatm = hbond.don_hatm();
			unsigned short int const Datm = rotamer->atom_base( Hatm );
			node.remove_atom_info( Hatm );
			node.remove_atom_info( Datm );
		} else {
			debug_assert( hbond.acc_res() == resid );
			unsigned short int const Aatm = hbond.acc_atm();
			node.remove_atom_info( Aatm );
		}//Acceptor
	}// for all hbonds
}

void delete_edges_with_degree_zero( scoring::hbonds::graph::HBondGraph & hb_graph ){
	std::list< utility::graph::LowMemEdge * > edges_to_delete;

	for ( utility::graph::LowMemEdgeListIter it = hb_graph.edge_list_begin(), end = hb_graph.edge_list_end();
			it != end; ++it ) {

		utility::graph::LowMemEdge * const edge = * it;
		if ( hb_graph.get_node( edge->get_first_node_ind() )->num_edges() == 1 &&
				hb_graph.get_node( edge->get_second_node_ind() )->num_edges() == 1 ) {
			edges_to_delete.push_back( edge );
		}
	}

	for ( utility::graph::LowMemEdge * doomed_edge : edges_to_delete ) {
		hb_graph.delete_edge( doomed_edge );
	}
}


/// @brief Construct an HBondGraph from a partial rotsets (like you might see during packing)
/// @detail BB-BB hbonds are only included for the first residue. This means that prolines are not
///         handled correctly. If proline is the first resdiue at a position and other residues
///         are being packed at that position, any hbonds to the other Ns will not be reported.
///   If one wishes to have BB-BB hbonds for all pairs, enable all 4 hbond terms for
///   scorefxn_sc and leave scorefxn_bb as a nullptr (or a blank scorefxn)
///   If your pose is symmetric, this internally desymmetrizes it and returns all hbonds
scoring::hbonds::graph::HBondGraphOP
hbond_graph_from_partial_rotsets(
	pose::Pose const & pose_in,
	pack::rotamer_set::RotamerSets const & original_rotsets,
	scoring::ScoreFunctionOP const & scorefxn_sc, // Only hbond_sc_bb and hbond_sc
	scoring::ScoreFunctionOP const & scorefxn_bb, // Only hbond_lr_bb and hbond_sr_bb
	pack::rotamer_set::RotamerSetsOP & complete_rotsets_out,
	utility::vector1<bool> & position_had_rotset,
	float minimum_hb_cut /* =0 */
) {
	bool separate_bb = bool(scorefxn_bb) && scorefxn_bb->get_nonzero_weighted_scoretypes().size() > 0;

	// Grab the symm_info if we have a symmetric pose.
	conformation::symmetry::SymmetryInfoCOP symm_info;
	if ( core::pose::symmetry::is_symmetric( pose_in ) ) {
		conformation::symmetry::SymmetricConformation const & symm_conf =
			dynamic_cast< conformation::symmetry::SymmetricConformation const & >( pose_in.conformation() );

		symm_info = symm_conf.Symmetry_Info();
	}

	// OK, so first we need to make a RotamerSets that has an entry for every single seqpos.
	//     We'll reuse whatever we can from original_rotsets and then start adding the current rotamers.

	// We need to score the pose. So we should probably make a copy so that we don't confuse the packer.

	pose::Pose pose_sc;
	pose::Pose pose_bb;

	if ( symm_info ) {
		pose_sc = pose_in;
		core::pose::symmetry::make_asymmetric_pose( pose_sc );
		pose_bb = pose_sc;
	} else {
		pose_sc = pose_in;
		pose_bb = pose_in;
	}

	scorefxn_sc->score( pose_sc );
	if ( separate_bb ) scorefxn_bb->score( pose_bb );

	// Pretend like we're going to design every position
	utility::vector1<bool> true_vect( pose_sc.size(), true );
	scorefxn_sc->setup_for_packing( pose_sc, true_vect, true_vect );
	if ( separate_bb ) scorefxn_bb->setup_for_packing( pose_bb, true_vect, true_vect );

	// Blank packer task that says we're going to pack every position
	pack::task::PackerTaskOP task = pack::task::TaskFactory::create_packer_task( pose_sc );

	// Our internal rotsets that the HBondGraph needs
	pack::rotamer_set::RotamerSetsOP rotsets( new pack::rotamer_set::RotamerSets() );
	rotsets->set_task( task );

	// Our rotsets for bb interactions
	pack::rotamer_set::RotamerSetsOP bb_rotsets( new pack::rotamer_set::RotamerSets() );
	bb_rotsets->set_task( task );

	position_had_rotset.clear();
	position_had_rotset.resize( pose_sc.size() );

	// Fill in the rotsets with either the RotamerSet from original_rotsets or make a new one from the current res
	for ( Size resnum = 1; resnum <= pose_sc.size(); resnum++ ) {

		// Get the asymmetric residue if this is symmetric.
		//  Remember, the rotsets that was passed is only for the asymetric unit.
		Size from_resnum = resnum;
		if ( symm_info ) {
			from_resnum = symm_info->bb_follows( resnum );
			if ( from_resnum == 0 ) from_resnum = resnum; // In asymetric unit
		}

		// We need to make a new rotset in either case
		//  Either because there was no rotset for this position
		//  Or because we'll confuse the scorefunction machinery if we don't make a new one
		pack::rotamer_set::RotamerSetOP rotset = utility::pointer::make_shared< pack::rotamer_set::RotamerSet_ >();
		rotset->set_resid( resnum );

		// The bb rotset
		pack::rotamer_set::RotamerSetOP bb_rotset = utility::pointer::make_shared< pack::rotamer_set::RotamerSet_ >();
		bb_rotset->set_resid( resnum );

		// This used to be so pretty before symmetry...
		// If we have a rotamer set in the asymetric unit, then use it. Otherwise just use the residue from the pose
		if ( original_rotsets.has_rotamer_set_for_residue( from_resnum ) ) {

			// Get the current rotset
			pack::rotamer_set::RotamerSetCOP original_rotset = original_rotsets.rotamer_set_for_residue( from_resnum );
			pack::rotamer_set::symmetry::SymmetricRotamerSet_COP symm_rotset;
			if ( symm_info ) {
				symm_rotset = utility::pointer::dynamic_pointer_cast< pack::rotamer_set::symmetry::SymmetricRotamerSet_ const >( original_rotset );
				runtime_assert( symm_rotset );
			}

			// Copy the rotamers from the old rotset to the new rotset
			for ( Size irot = 1; irot <= original_rotset->num_rotamers(); irot++ ) {
				core::conformation::ResidueOP rotamer;
				if ( from_resnum != resnum ) {
					rotamer = symm_rotset->orient_rotamer_to_symmetric_partner( pose_in, *original_rotset->rotamer( irot ), resnum, true );
				} else {
					rotamer = original_rotset->rotamer( irot )->clone();
				}
				rotamer->seqpos( resnum );
				rotset->add_rotamer( *rotamer );
			}

			// Only add the first rotamer to the backbone rotset
			core::conformation::ResidueOP rotamer;
			if ( from_resnum != resnum ) {
				rotamer = symm_rotset->orient_rotamer_to_symmetric_partner( pose_in, *original_rotset->rotamer( 1 ), resnum, true );
			} else {
				rotamer = original_rotset->rotamer( 1 )->clone();
			}
			rotamer->seqpos( resnum );
			bb_rotset->add_rotamer( *rotamer );
			position_had_rotset[resnum] = true;
		} else {

			// If we don't have a rotset. Just use the poses's residue
			rotset->add_rotamer( *pose_sc.residue(resnum).create_rotamer() );
			bb_rotset->add_rotamer( *pose_bb.residue(resnum).create_rotamer() );
			position_had_rotset[resnum] = false;
		}

		scorefxn_sc->prepare_rotamers_for_packing( pose_sc, *rotset );
		if ( separate_bb ) scorefxn_bb->prepare_rotamers_for_packing( pose_bb, *bb_rotset );
		rotsets->set_explicit_rotamers( rotsets->resid_2_moltenres( resnum ), rotset );
		bb_rotsets->set_explicit_rotamers( bb_rotsets->resid_2_moltenres( resnum ), bb_rotset );
	}
	rotsets->update_offset_data();
	bb_rotsets->update_offset_data();

	// Use some of the HBNet machinery to load the HBondGraph for interactions involving sc
	scoring::hbonds::graph::HBondGraphOP hb_graph;
	hb_graph = pack::hbonds::create_init_and_create_edges_for_hbond_graph(
		rotsets,    // our temporary rotsets
		* scorefxn_sc, pose_sc,
		minimum_hb_cut, // allow all hbonds during interation graph generation
		1e6,            // allow all clashes
		minimum_hb_cut, // only allow good hbonds to make it into the final graph
		true            // include backbone atoms ( we want bb_sc too )
	);

	if ( separate_bb ) {

		// Use some of the HBNet machinery to load the HBondGraph for bb_bb only
		scoring::hbonds::graph::HBondGraphOP bb_hb_graph;
		bb_hb_graph = pack::hbonds::create_init_and_create_edges_for_hbond_graph(
			bb_rotsets,    // our temporary rotsets
			* scorefxn_bb, pose_bb,
			minimum_hb_cut, // allow all hbonds during interation graph generation
			1e6,            // allow all clashes
			minimum_hb_cut, // only allow good hbonds to make it into the final graph
			true            // include backbone atoms
		);

		hb_graph->merge( *bb_hb_graph, true );
	}

	complete_rotsets_out = rotsets;

	return hb_graph;
}


} //hbonds
} //pack
} //core
