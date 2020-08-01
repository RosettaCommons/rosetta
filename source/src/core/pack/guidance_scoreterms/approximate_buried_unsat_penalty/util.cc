// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/util.cc
/// @brief  Utility functions for ApproximateBuriedUnsatPenalty
/// @author Brian Coventry (bcov@uw.edu)

#include <core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/util.hh>

#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/hbonds/HBondGraph_util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/atomic_depth/AtomicDepth.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/symmetry/util.hh>

#include <basic/Tracer.hh>

#include <utility/graph/Graph.hh>
#include <utility/pointer/memory.hh>
#include <boost/format.hpp>

#include <chrono>
#include <thread>
#include <cstdint>
#include <iostream>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace approximate_buried_unsat_penalty {

static basic::Tracer TR( "core.pack.guidance_scoreterms.approximate_buried_unsat_penalty" );
static basic::Tracer TR2( "approx" );

using basic::datacache::ResRotPair;
using basic::datacache::ResRotPairHasher;


void
prestore_sc_oversat(
	std::unordered_map< ResRotPair, float, ResRotPairHasher > & this_bur_rotamer_oversat_map,
	Size resnum1, Size rotamer1, bool is_bb1,
	Size resnum2, Size rotamer2, bool is_bb2,
	float score ) {

	// If its bb, we don't want it associated with a rotamer
	//  Also, changing these here does not affect should_store_hbond()
	if ( is_bb1 ) {
		rotamer1 = 0;
	}
	if ( is_bb2 ) {
		rotamer2 = 0;
	}

	bool swap = resnum2 < resnum1;
	Size r1  = swap ? resnum2  : resnum1;
	Size ro1 = swap ? rotamer2 : rotamer1;
	// Size at1 = swap ? atom2    : atom1;

	Size r2  = swap ? resnum1  : resnum2;
	Size ro2 = swap ? rotamer1 : rotamer2;
	// Size at2 = swap ? atom1    : atom2;


	ResRotPair res_rot_pair( r1, ro1, r2, ro2 );

	this_bur_rotamer_oversat_map[ res_rot_pair ] += score;
}

void
accumulate_bur_rotamer_oversats(
	std::unordered_map< ResRotPair, float, ResRotPairHasher > & this_bur_rotamer_oversat_map,
	std::unordered_map< ResRotPair, float, ResRotPairHasher > & this_seqpos_oversat_map
);

void
accumulate_seqpos_oversats(
	basic::datacache::CacheableResRotPairFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & complete_rotsets,
	utility::vector1<bool> const & is_asu,
	ReweightData & reweight_data,
	std::unordered_map< ResRotPair, float, ResRotPairHasher > & this_seqpos_oversat_map,
	core::conformation::symmetry::SymmetryInfoCOP const & symm_info /*can be nullptr*/
);


// How does this work?
//
// Consider the following scenarios:
//
// For these scenarios, consider only A to be buried and the
//   approximate_buried_unsat_penalty to be 1.
//     0              1              2             3
//     |     A +1     |     A +1     |     A +1    |     A +1
//     |              |              |   H -1      |   H -1   H -1
//     |              |      H -1    |  D   H -1   | *D   H -1 D* +1
//     |              |      D       |   +1 D      |   +1 D +1
//     |              |              |             |
// score:    +1             +0             +0            +1
//
// 0. There is a buried acceptor +1.
// 1. There is a buried acceptor +1. There is a hydrogen bond -1.
// 2. There is a buried acceptor +1. There are 2 hydrogen bonds -2.
//     There is one oversaturation penalty. +1
// 3. There is a buried acceptor +1. There are 3 hydrogen bonds -3.
//     There are three oversaturations +3
//
// The approximate_buried_unsat_penalty gets its name for a reason however;
//   as the number of oversaturations at a single atom continues to increase,
//   the oversaturation penalty goes up quadratically
//
//            4              5              6             7
// score:    +3             +6             +10           +15
//
// 4. There is a buried acceptor +1. There are 4 hydrogen bonds -4.
//     There are 4 oversaturations +6
// 5. There is a buried acceptor +1. There are 5 hydrogen bonds -5.
//     There are 5 oversaturations +10
// 6. There is a buried acceptor +1. There are 6 hydrogen bonds -6.
//     There are 6 oversaturations +15
// 7. There is a buried acceptor +1. There are 7 hydrogen bonds -7.
//     There are 7 oversaturations +21
//
// overaturation_penalty_score = ( x * ( x - 1 ) ) / 2
//
//
// Here's where the magic is:
//
// It turns out that it's possible to represent the above algorithm
//   as a twobody, a onebody, and a zerobody energy during packing.
//
// A full three-body calculation is performed before packing and then
//   the two-body, one-body, and zero-body scores are stored in the pose
//   to be looked up during scoring.
//
// These rules generate the tables:
//
//   For buried unsats:
//     Assign one-body penalty if on sidechain of rotamer
//     Assign zero-body penalty if on sidechain of non-packable pose or backbone
//
//   For hbonds:
//     Assign two-body bonus if both are rotamers
//     Assign one-body bonus if one is rotamer and other is non-packable pose
//     Assign zero-body bonus if both are non-packable pose
//
//   For penalty: (two donors are satisfying the same acceptor or vis versa)
//     Assign two-body bonus if both are rotamers
//     Assign one-body bonus if one is rotamer and other is non-packable pose
//     Assign zero-body bonus if both are non-packable pose
//
//
//
// Caveats:
//
//   Lysine: Lysine actually wants 3 hbonds and this is an issue with the above
//      scheme. But consider what happens if its buried unsat penalty is 3,
//      the bonus for making a hbond is 2, and the oversaturation penalty is 1
//
//  hbonds:  0     1      2      3      4      5      6      7
//  score:   3     1      0      0      1      3      6      10
//
//
//
//
//  Rotamers: If there is an oversat against an atom that exists in multiple
//       rotamers at one position (e.g. a backbone atom ), then every single
//       rotamer with that atom, a penalty will be assigned to the oversat
//       edge.
//
//       Fix: Instead of accumulating an oversat, store the information about
//                the oversat residues and the target residue. Then take the max
//                for each unique set and accumulate those
//
//
//
//  Oversat penalty modifications:
//  You can change the penalty for oversaturation with a commandline flag.
//      Here are the tables for what happens
//
//
//  hbonds:  0    1    2    3    4    5    6
//  scores:
// setting:
//
//     0     1    0   -1   -2   -3   -4   -5
//    0.5    1    0  -0.5 -0.5   0    1   2.5
//     1     1    0    0    1    3    6   10
//    1.5    1    0   0.5  2.5   6   11   17.5
//     2     1    0    1    4    9   16   25
//
//

basic::datacache::CacheableResRotPairFloatMapOP
three_body_approximate_buried_unsat_calculation(
	pose::Pose const & pose,
	pack::rotamer_set::RotamerSetsOP const & rotsets,
	scoring::ScoreFunctionOP const & scorefxn_sc, // Only hbond_sc_bb and hbond_sc
	scoring::ScoreFunctionOP const & scorefxn_bb, // Only hbond_lr_bb and hbond_sr_bb
	float atomic_depth_cutoff /*= 4.5f*/,
	float atomic_depth_probe_radius /*= 2.3f*/,
	float atomic_depth_resolution /*= 0.5f*/,
	float minimum_hb_cut /* =0 */,
	bool all_atoms_active /* =false */,
	float oversat_penalty_in /* =1 */,
	bool assume_const_backbone /* =true */,
	UnsatCorrectionOptions const & cor_opt /* = UnsatCorrectionOptions() */,
	HBondBonusOptions const & bonus_opt /* = HBondBonusOptions() */
) {

	using namespace std::chrono;

	const float PENALTY = 1.0f;
	basic::datacache::CacheableResRotPairFloatMapOP score_map ( utility::pointer::make_shared< basic::datacache::CacheableResRotPairFloatMap >() );
	score_map->set_shallow_copy( true );

	milliseconds start = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

	conformation::symmetry::SymmetryInfoCOP symm_info(nullptr); //Set to non-null iff pose is symmetric, below.

	// We don't strictly need to do this. But it will keep the memory requirements down for giant symmetries
	utility::vector1<bool> is_asu( pose.size(), true );
	if ( core::pose::symmetry::is_symmetric( pose ) ) {

		conformation::symmetry::SymmetricConformation const & symm_conf =
#ifndef NDEBUG
			// Debug build: do a slow dynamic cast
			dynamic_cast< conformation::symmetry::SymmetricConformation const & >( pose.conformation() );
#else
			// Release mode: do a fast static cast
			static_cast< conformation::symmetry::SymmetricConformation const & >( pose.conformation() );
#endif

		symm_info = symm_conf.Symmetry_Info();
		for ( Size seqpos = 1; seqpos <= pose.size(); seqpos ++ ) {
			is_asu[seqpos] = symm_info->bb_follows( seqpos ) == 0;
		}
	}

	if ( TR.Trace.visible() ) {
		for ( Size irotset = 1; irotset <= rotsets->nmoltenres(); irotset++ ) {
			Size resnum = rotsets->moltenres_2_resid(irotset);
			pack::rotamer_set::RotamerSetCOP const & rotset = rotsets->rotamer_set_for_residue(resnum);

			for ( Size irot = 1; irot <= rotset->num_rotamers(); irot++ ) {
				conformation::ResidueCOP rotamer = rotset->rotamer( irot );

				TR2 << boost::str(boost::format("Res: %i %s Rot: %i Chis: ")%resnum%rotamer->name3()%irot);
				for ( Size chi = 1; chi <= rotamer->nchi(); chi++ ) {
					TR2 << boost::str(boost::format("%7.3f ")%rotamer->chi(chi) );
				}
				TR2 << std::endl;
			}
		}
	}

	// First we need to get the HBondGraph and its associated complete RotamerSets

	utility::vector1<bool> position_had_rotset;
	pack::rotamer_set::RotamerSetsOP complete_rotsets_out = nullptr;

	TR << "Building hbond graph" << std::endl;

	milliseconds hbond_start = duration_cast< milliseconds >( system_clock::now().time_since_epoch());

	scoring::hbonds::graph::HBondGraphOP hb_graph;
	hb_graph = hbonds::hbond_graph_from_partial_rotsets( pose, *rotsets, scorefxn_sc, scorefxn_bb,
		complete_rotsets_out, position_had_rotset, minimum_hb_cut );


	milliseconds hbond_end = duration_cast< milliseconds >( system_clock::now().time_since_epoch());
	if ( TR.Debug.visible() ) TR.Debug << "TIMING: HBOND: " << (hbond_end - hbond_start).count() << std::endl;


	TR << "Hbond graph has: " << hb_graph->num_edges() << " edges requiring: " << hb_graph->getTotalMemoryUsage() << " bytes" << std::endl;

	pack::rotamer_set::RotamerSetsOP complete_rotsets =
		std::dynamic_pointer_cast<pack::rotamer_set::RotamerSets>( complete_rotsets_out );

	runtime_assert( complete_rotsets );
	runtime_assert( hb_graph->num_nodes() == complete_rotsets->nrotamers() );

	// std::cout << "Sleeping after hbond graph" << std::endl;
	// std::this_thread::sleep_for(std::chrono::milliseconds(5000));

	// SYMMETRY HACK
	// We need to know how many rotamers are available at each position when it comes time to score rotamers
	// Store this information into rotamer 0
	for ( Size i = 1; i <= pose.size(); i++ ) {
		score_map->map()[ResRotPair(i, 0, 0, 0)] = complete_rotsets->nrotamers_for_moltenres( complete_rotsets->resid_2_moltenres( i ) );
	}

	milliseconds depth_start = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );

	// Next, we need to create the AtomicDepth object and figure out which nodes are buried
	scoring::atomic_depth::AtomicDepth atomic_depth( pose, atomic_depth_probe_radius, true, atomic_depth_resolution );

	milliseconds depth_end = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );
	if ( TR.Debug.visible() ) TR.Debug << "TIMING: DEPTH: " << (depth_end - depth_start).count() << std::endl;

	chemical::AtomTypeSet const & atom_type_set = pose.residue(1).type().atom_type_set();

	// std::cout << "Sleeping after atomic depth" << std::endl;
	// std::this_thread::sleep_for(std::chrono::milliseconds(5000));

	// Now we begin the algorithm

	ReweightData reweight_data( pose, rotsets->task() );
	std::unordered_map< ResRotPair, float, ResRotPairHasher > this_bur_rotamer_oversat_map;
	std::unordered_map< ResRotPair, float, ResRotPairHasher > this_seqpos_oversat_map;

	//////////////////// Iterate over all polar heavy atoms looking for buried ones //////////////////////////

	Size last_resnum = 0;

	for ( Size ihbnode = 1; ihbnode <= hb_graph->num_nodes(); ihbnode++ ) {

		// These nodes are one-to-one with the complete_rotsets rotamers
		scoring::hbonds::graph::HBondNode * hbnode = hb_graph->get_node( ihbnode );
		runtime_assert( hbnode );
		conformation::ResidueCOP rotamer = complete_rotsets->rotamer( ihbnode );
		Size node_resnum = complete_rotsets->res_for_rotamer( ihbnode );


		// Performance hack ////////////////////////
		// If we accumulate these after each residue. We can reuse the
		//   vectors and lower our memory footprint.
		if ( node_resnum != last_resnum ) {

			accumulate_seqpos_oversats( score_map, complete_rotsets, is_asu, reweight_data, this_seqpos_oversat_map, symm_info );

			last_resnum = node_resnum;
		}


		// This is all the polar atoms on the rotamer
		boost::container::flat_set< scoring::hbonds::graph::AtomInfo > const & hbnode_atoms =
			hbnode->polar_sc_atoms_not_satisfied_by_background();

		// Prepare the vector of atoms for the burial calculation
		utility::vector1< conformation::Atom > atoms;
		for ( scoring::hbonds::graph::AtomInfo const & info : hbnode_atoms ) {
			atoms.push_back( rotamer->atom( info.local_atom_id() ) );
		}

		utility::vector1<Real> depths = atomic_depth.calcdepth( atoms, atom_type_set );

		for ( Size inodeat = 1; inodeat <= hbnode_atoms.size(); inodeat++ ) {

			scoring::hbonds::graph::AtomInfo const & info =
				* std::next( hbnode_atoms.begin(), inodeat - 1 );

			if ( info.is_hydrogen() ) continue;
			if ( depths[inodeat] < atomic_depth_cutoff ) continue;

			// Now we know we found a buried heavy atom

			// Backbone atom gate 1 //////////////////
			// Only allow backbone atoms from rotamer 1
			if ( assume_const_backbone && info.is_backbone() ) {
				if ( complete_rotsets->rotid_on_moltenresidue( ihbnode ) != 1 ) continue;
			}


			//////////////////// Identify all rotamers that hbond to this atom ////////////////////////////////

			// These aren't quite AtomIDs because they are actually the rotamer number not resnum
			utility::vector1< id::AtomID > hbonding_rotamers;

			for ( utility::graph::LowMemEdgeListIter it = hbnode->edge_list_begin( *hb_graph );
					it != hbnode->edge_list_end( *hb_graph );
					++it ) {

				scoring::hbonds::graph::HBondEdge * hb_edge =
					static_cast< scoring::hbonds::graph::HBondEdge * >( *it );
				Size first_node = hb_edge->get_first_node_ind();
				Size second_node = hb_edge->get_second_node_ind();

				bool we_are_the_first_node = first_node == ihbnode;
				Size other_node = we_are_the_first_node ? second_node : first_node;

				conformation::ResidueCOP other_rotamer = complete_rotsets->rotamer( other_node );

				for ( scoring::hbonds::graph::HBondInfo const & hb_info : hb_edge->hbonds() ) {
					//          first is donor   0        1
					//  we are first   0        don      acc
					//                 1        acc      don
					//                                           xor
					bool we_are_the_donor = ! ( we_are_the_first_node ^ hb_info.first_node_is_donor() );

					Size check_atom = we_are_the_donor ? hb_info.local_atom_id_D() : hb_info.local_atom_id_A();
					if ( info.local_atom_id() != check_atom ) continue;


					// Ok, now we know that this hbond actually involves the buried atom we're looking at.
					//  Which implies that the node we are looking at hbonds to our buried atom

					Size other_atom = we_are_the_donor ? hb_info.local_atom_id_A() : hb_info.local_atom_id_D();

					// Backbone atom gate 2 //////////////////
					// Only allow backbone atoms from rotamer 1
					if ( assume_const_backbone && other_rotamer->atom_is_backbone(other_atom) ) {
						if ( complete_rotsets->rotid_on_moltenresidue( other_node ) != 1 ) {
							continue;
						}
					}

					// For hydroxyl_wants_h, don't count h-bonds to hydroxyl O
					if ( cor_opt.hydroxyl_wants_h && info.is_hydroxyl() && ! we_are_the_donor ) continue;

					hbonding_rotamers.emplace_back( other_atom, other_node );
				}
			}

			//////////////////////////// Perform the scoring bookkeeping //////////////////////////////////////

			// Need to account for NH3 which needs special treatment

			Size num_hs = 0;
			if ( rotamer->heavyatom_has_polar_hydrogens( info.local_atom_id() ) ) {
				num_hs = rotamer->attached_H_end( info.local_atom_id() ) - rotamer->attached_H_begin( info.local_atom_id() ) + 1;
			}

			bool is_N = rotamer->atom_type( info.local_atom_id() ).element() == "N";
			bool is_carboxyl = rotamer->atom_type( info.local_atom_id() ).atom_type_name() == "OOC";

			// default values that look like this. 1 0 0 1 3
			float buried_penalty =  PENALTY;
			float hbond_bonus =     -PENALTY;
			float oversat_penalty = oversat_penalty_in;

			if ( num_hs >= 3 ) { // standard LYS handling. 3 1 0 0 1 3
				buried_penalty =  PENALTY * 3;
				hbond_bonus =    -PENALTY * 2;
				oversat_penalty = oversat_penalty_in;

			} else if ( cor_opt.nh2_wants_2 && is_N && num_hs == 2 ) {   // force NH2 to have 2. 2 0.5 0 0.5 2
				buried_penalty =  PENALTY * 2;
				hbond_bonus =    -PENALTY * 1.5;
				oversat_penalty = oversat_penalty_in;

			} else if ( cor_opt.nh1_wants_1 && is_N && num_hs <= 1 ) {   // force NH1 to have 1. 1 0 1 3
				buried_penalty =  PENALTY;
				hbond_bonus =    -PENALTY;
				oversat_penalty = oversat_penalty_in * 2;

			} else if ( cor_opt.hydroxyl_wants_h && info.is_hydroxyl() ) { // extra tight bound on hydroxyls. 1 0 1 3
				buried_penalty =  PENALTY;
				hbond_bonus =    -PENALTY;
				oversat_penalty = oversat_penalty_in * 2;

			} else if ( cor_opt.carboxyl_wants_2 && is_carboxyl ) { // force carboxyls to have 2. 2 0.5 0 0.5 2
				buried_penalty =  PENALTY * 2;
				hbond_bonus =    -PENALTY * 1.5;
				oversat_penalty = oversat_penalty_in;

			}

			// First we store the buried unsat penalty to whomever owns it

			bool buried_should_store = all_atoms_active;
			buried_should_store |= position_had_rotset[ node_resnum ] && ! ( assume_const_backbone && info.is_backbone() );

			if ( buried_should_store ) {
				add_to_onebody( score_map, complete_rotsets, is_asu, node_resnum, complete_rotsets->rotid_on_moltenresidue( ihbnode ), buried_penalty );
				if ( TR.Trace.visible() ) {
					TR2 << boost::str(boost::format("BUR: OneB Res: %i %s Rot: %i Atom: %s Score: %6.3f")
						%node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
						%rotamer->atom_name(info.local_atom_id())%buried_penalty) << std::endl;
				}
			} else {
				if ( TR.Trace.visible() ) {
					TR2 << boost::str(boost::format("BUR: ZeroB Res: %i %s Rot: %i Atom: %s Score: %6.3f")
						%node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
						%rotamer->atom_name(info.local_atom_id())%buried_penalty) << std::endl;
				}
				// zero_body += buried_penalty;
			}

			// Next we store the hydrogen bond bonuses to the hbonding_rotamers

			utility::vector1<bool> should_store( hbonding_rotamers.size() );
			utility::vector1<bool> is_bb( hbonding_rotamers.size() );

			for ( Size iother = 1; iother <= hbonding_rotamers.size(); iother++ ) {
				id::AtomID const & atid = hbonding_rotamers[ iother ];

				Size other_rotamerid = atid.rsd();
				Size other_atom = atid.atomno();
				Size other_resnum = complete_rotsets->res_for_rotamer( other_rotamerid );
				conformation::ResidueCOP other_rotamer = complete_rotsets->rotamer( other_rotamerid );
				bool other_is_backbone = other_rotamer->atom_is_backbone( other_atom );

				bool other_should_store = all_atoms_active;
				other_should_store |= position_had_rotset[ other_resnum ] && ! ( assume_const_backbone && other_is_backbone );

				is_bb[ iother ] = other_is_backbone;
				should_store[ iother ] = other_should_store;

				bool do_store = true;

				// if ( ! do_store ) std::cout << "Skip: ";
				// if ( ! do_store ) continue;

				if ( buried_should_store && other_should_store ) {
					if ( do_store ) {
						add_to_twobody( score_map, complete_rotsets, is_asu, reweight_data,
							node_resnum, complete_rotsets->rotid_on_moltenresidue( ihbnode ),
							other_resnum, complete_rotsets->rotid_on_moltenresidue( other_rotamerid ),
							hbond_bonus, symm_info );
					}
					if ( TR.Trace.visible() ) {
						TR2 << boost::str(boost::format("HYD: TwoB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
							%node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
							%other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_rotamerid )
							%hbond_bonus) << std::endl;
					}
				} else if ( buried_should_store ) {
					if ( do_store ) add_to_onebody( score_map, complete_rotsets, is_asu, node_resnum, complete_rotsets->rotid_on_moltenresidue( ihbnode ), hbond_bonus );
					if ( TR.Trace.visible() ) {
						TR2 << boost::str(boost::format("HYD: OneB *Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
							%node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
							%other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_rotamerid )
							%hbond_bonus) << std::endl;
					}
				} else if ( other_should_store ) {
					if ( do_store ) add_to_onebody( score_map, complete_rotsets, is_asu, other_resnum, complete_rotsets->rotid_on_moltenresidue( other_rotamerid ), hbond_bonus );
					if ( TR.Trace.visible() ) {
						TR2 << boost::str(boost::format("HYD: OneB Res: %i %s Rot: %i *Res: %i %s Rot: %i Score: %6.3f")
							%node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
							%other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_rotamerid )
							%hbond_bonus) << std::endl;
					}
				} else {
					// zero_body += hbond_bonus;
					if ( TR.Trace.visible() ) {
						TR2 << boost::str(boost::format("HYD: ZeroB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
							%node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
							%other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_rotamerid )
							%hbond_bonus) << std::endl;
					}
				}
			}


			// Finally we store the oversaturation penalties to all pairs of hbonding_rotamers

			// Note!! This for loop is incredibly hot in a performance sense. Only highly optimized code
			//        should go here.
			// Upper triangle for loop
			for ( Size ii = 1; ii <= hbonding_rotamers.size(); ii++ ) {
				for ( Size jj = 1; jj <= hbonding_rotamers.size(); jj++ ) {
					if ( ii >= jj ) continue;

					id::AtomID const & ii_atid = hbonding_rotamers[ ii ];
					Size ii_rotamerid = ii_atid.rsd();
					Size ii_resnum = complete_rotsets->res_for_rotamer( ii_rotamerid ); // potentially slow
					bool ii_should_store = should_store[ ii ];

					id::AtomID const & jj_atid = hbonding_rotamers[ jj ];
					Size jj_rotamerid = jj_atid.rsd();
					Size jj_resnum = complete_rotsets->res_for_rotamer( jj_rotamerid ); // potentiall slow
					bool jj_should_store = should_store[ jj ];

					if ( ii_resnum == jj_resnum && ii_rotamerid != jj_rotamerid && !is_bb[ii] && !is_bb[jj] ) continue;

					bool do_store = true;

					if ( ! all_atoms_active ) {

						if ( assume_const_backbone && info.is_backbone() ) {

							// Only one representative for this atom should come through, so we're good

						} else {
							// ERROR FIX 3 /////////////
							// If the oversat atom is a side chain atom and multiple rotamers at that resid contain
							//   a similarly placed atom, then there will be excessive counting of oversats for
							//   each rotamer that has an atom like that.
							//
							// Fix: Dict key: ov1_resnum, ov1_rotamer, ov2_resnum, ov2_rotamer, tar_resnum
							//      Dict value: vector1< float > [ tar_rotamer ]
							//        And then sum oversats with the correct dict key and vector1 offset
							//
							//      At the end, assign as an oversat the max of that vector
							//
							do_store = false;
							prestore_sc_oversat( this_bur_rotamer_oversat_map,
								ii_resnum, complete_rotsets->rotid_on_moltenresidue( ii_rotamerid ), assume_const_backbone && is_bb[ ii ],
								jj_resnum, complete_rotsets->rotid_on_moltenresidue( jj_rotamerid ), assume_const_backbone && is_bb[ jj ],
								oversat_penalty );
						}
					}

					if ( ii_should_store && jj_should_store ) {
						if ( do_store ) {
							add_to_twobody( score_map, complete_rotsets, is_asu, reweight_data,
								ii_resnum, complete_rotsets->rotid_on_moltenresidue( ii_rotamerid ),
								jj_resnum, complete_rotsets->rotid_on_moltenresidue( jj_rotamerid ),
								oversat_penalty, symm_info );
						}
						if ( TR.Trace.visible() ) {
							TR2 << boost::str(boost::format("OVR: TwoB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
								%ii_resnum%complete_rotsets->rotamer( ii_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( ii_rotamerid )
								%jj_resnum%complete_rotsets->rotamer( jj_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( jj_rotamerid )
								%oversat_penalty) << ( do_store ? "" : " skip") << std::endl;
						}
					} else if ( ii_should_store ) {
						if ( do_store ) add_to_onebody( score_map, complete_rotsets, is_asu, ii_resnum, complete_rotsets->rotid_on_moltenresidue( ii_rotamerid ), oversat_penalty );
						if ( TR.Trace.visible() ) {
							TR2 << boost::str(boost::format("OVR: OneB *Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
								%ii_resnum%complete_rotsets->rotamer( ii_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( ii_rotamerid )
								%jj_resnum%complete_rotsets->rotamer( jj_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( jj_rotamerid )
								%oversat_penalty) << ( do_store ? "" : " skip") << std::endl;
						}
					} else if ( jj_should_store ) {
						if ( do_store ) add_to_onebody( score_map, complete_rotsets, is_asu, jj_resnum, complete_rotsets->rotid_on_moltenresidue( jj_rotamerid ), oversat_penalty );
						if ( TR.Trace.visible() ) {
							TR2 << boost::str(boost::format("OVR: OneB Res: %i %s Rot: %i *Res: %i %s Rot: %i Score: %6.3f")
								%ii_resnum%complete_rotsets->rotamer( ii_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( ii_rotamerid )
								%jj_resnum%complete_rotsets->rotamer( jj_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( jj_rotamerid )
								%oversat_penalty) << ( do_store ? "" : " skip") << std::endl;
						}
					} else {
						// zero_body += hbond_bonus;
						if ( TR.Trace.visible() ) {
							TR2 << boost::str(boost::format("OVR: ZeroB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
								%ii_resnum%complete_rotsets->rotamer( ii_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( ii_rotamerid )
								%jj_resnum%complete_rotsets->rotamer( jj_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( jj_rotamerid )
								%oversat_penalty) << ( do_store ? "" : " skip") << std::endl;
						}
					}
				}
			}
		}
		// Determine if this buried rotamer set any new max values in the oversat map
		accumulate_bur_rotamer_oversats( this_bur_rotamer_oversat_map, this_seqpos_oversat_map );

	}

	// Very important!!! This must be here!!! This catches the final residue.
	accumulate_seqpos_oversats( score_map, complete_rotsets, is_asu, reweight_data, this_seqpos_oversat_map, symm_info );


	//////////////////////////////////// HBond Bonus //////////////////////////////////////
	// This is separate from buried unsats. This is just the best place to put it

	// Iterate accross all rotamers
	// Iterate accross all hbonds from rotamer / Only look at hbonds to greater rotamers
	// if ( assume_const_backbone ):
	//     hbond_to_sidechain: -- twobody this -> that
	//     hbond_to_backbone: -- onebody -> this
	//     hbond_from_backbone -- onebody -> that
	//     backbone_to_backbone -- zerobody
	// else:
	//     hbond_to_sidechain: -- twobody this -> that
	//     hbond_to_backbone: -- twobody this -> that
	//     hbond_from_backbone: -- twobody this -> that
	//     backbone_to_backbone -- twobody this -> that

	if ( bonus_opt.any() ) {

		for ( Size ihbnode = 1; ihbnode <= hb_graph->num_nodes(); ihbnode++ ) {

			// These nodes are one-to-one with the complete_rotsets rotamers
			scoring::hbonds::graph::HBondNode * hbnode = hb_graph->get_node( ihbnode );
			runtime_assert( hbnode );
			conformation::ResidueCOP rotamer = complete_rotsets->rotamer( ihbnode );
			Size node_resnum = complete_rotsets->res_for_rotamer( ihbnode );

			// This is all the polar atoms on the rotamer
			boost::container::flat_set< scoring::hbonds::graph::AtomInfo > const & hbnode_atoms =
				hbnode->polar_sc_atoms_not_satisfied_by_background();

			for ( Size inodeat = 1; inodeat <= hbnode_atoms.size(); inodeat++ ) {

				scoring::hbonds::graph::AtomInfo const & info =
					* std::next( hbnode_atoms.begin(), inodeat - 1 );

				if ( info.is_hydrogen() ) continue;

				// Now we know we found a heavy atom

				// Backbone atom gate 1 //////////////////
				// Only allow backbone atoms from rotamer 1
				if ( assume_const_backbone && info.is_backbone() ) {
					if ( complete_rotsets->rotid_on_moltenresidue( ihbnode ) != 1 ) continue;
				}

				bool we_should_store = all_atoms_active;
				we_should_store |= position_had_rotset[ node_resnum ] && ! ( assume_const_backbone && info.is_backbone() );


				for ( utility::graph::LowMemEdgeListIter it = hbnode->edge_list_begin( *hb_graph );
						it != hbnode->edge_list_end( *hb_graph );
						++it ) {

					scoring::hbonds::graph::HBondEdge * hb_edge =
						static_cast< scoring::hbonds::graph::HBondEdge * >( *it );
					Size first_node = hb_edge->get_first_node_ind();
					Size second_node = hb_edge->get_second_node_ind();

					bool we_are_the_first_node = first_node == ihbnode;
					Size other_node = we_are_the_first_node ? second_node : first_node;

					// Only allow hbonds in the downstream direction
					if ( other_node <= ihbnode ) continue;

					conformation::ResidueCOP other_rotamer = complete_rotsets->rotamer( other_node );

					for ( scoring::hbonds::graph::HBondInfo const & hb_info : hb_edge->hbonds() ) {
						//          first is donor   0        1
						//  we are first   0        don      acc
						//                 1        acc      don
						//                                           xor
						bool we_are_the_donor = ! ( we_are_the_first_node ^ hb_info.first_node_is_donor() );

						Size check_atom = we_are_the_donor ? hb_info.local_atom_id_D() : hb_info.local_atom_id_A();
						if ( info.local_atom_id() != check_atom ) continue;


						// Ok, now we know that this hbond actually involves the buried atom we're looking at.
						//  Which implies that the node we are looking at hbonds to our buried atom

						Size other_atom = we_are_the_donor ? hb_info.local_atom_id_A() : hb_info.local_atom_id_D();

						// Backbone atom gate 2 //////////////////
						// Only allow backbone atoms from rotamer 1
						if ( assume_const_backbone && other_rotamer->atom_is_backbone(other_atom) ) {
							if ( complete_rotsets->rotid_on_moltenresidue( other_node ) != 1 ) {
								continue;
							}
						}

						Size other_resnum = complete_rotsets->res_for_rotamer( other_node );

						bool other_should_store = all_atoms_active;
						other_should_store |= position_had_rotset[ other_resnum ] && ! ( assume_const_backbone && other_rotamer->atom_is_backbone(other_atom) );

						//////////////////////////////////// Calculate and apply HBond Bonus //////////////////////////////////////

						conformation::ResidueCOP donor_res = we_are_the_donor ? rotamer : other_rotamer;
						Size donor_resnum = we_are_the_donor ? node_resnum : other_resnum;
						Size donor_atom = hb_info.local_atom_id_D();

						conformation::ResidueCOP acceptor_res = ! we_are_the_donor ? rotamer : other_rotamer;
						Size acceptor_resnum = ! we_are_the_donor ? node_resnum : other_resnum;
						Size acceptor_atom = hb_info.local_atom_id_A();

						// std::cout << "HERE " << bonus_opt.ser_to_helix_bb() << node_resnum << " " << other_resnum << " " << donor_resnum << " "
						//  << donor_res->atom_name(donor_atom) << " " << acceptor_resnum << " " << acceptor_res->atom_name(acceptor_atom) << std::endl;

						float bonus = 0;

						if ( bonus_opt.cross_chain_bonus() != 0 ) {
							if ( pose.chain( donor_resnum ) != pose.chain( acceptor_resnum ) ) {
								bonus += bonus_opt.cross_chain_bonus();
							}
						}

						if ( bonus_opt.ser_to_helix_bb() != 0 ) {

							// i -> i-4 hbonds where i is SER or THR, donor atom is OG or OG1 and acceptor atom is O
							if ( donor_resnum - acceptor_resnum == 4 &&
									pose.chain( node_resnum ) == pose.chain( other_resnum ) &&
									acceptor_res->atom_name( acceptor_atom ) == " O  " ) {

								if ( donor_res->name1() == 'T' && donor_res->atom_name( donor_atom ) == " OG1" ) {
									bonus += bonus_opt.ser_to_helix_bb();
								}
								if ( donor_res->name1() == 'S' && donor_res->atom_name( donor_atom ) == " OG " ) {
									bonus += bonus_opt.ser_to_helix_bb();
								}
							}
						}

						if ( bonus != 0 ) {

							if ( we_should_store && other_should_store ) {
								add_to_twobody( score_map, complete_rotsets, is_asu, reweight_data,
									node_resnum, complete_rotsets->rotid_on_moltenresidue( ihbnode ),
									other_resnum, complete_rotsets->rotid_on_moltenresidue( other_node ), bonus, symm_info );

								if ( TR.Trace.visible() ) {
									TR2 << boost::str(boost::format("BON: TwoB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
										%node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
										%other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_node )
										%bonus) << std::endl;
								}
							} else if ( we_should_store ) {
								add_to_onebody( score_map, complete_rotsets, is_asu, node_resnum, complete_rotsets->rotid_on_moltenresidue( ihbnode ), bonus );
								if ( TR.Trace.visible() ) {
									TR2 << boost::str(boost::format("BON: OneB *Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
										%node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
										%other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_node )
										%bonus) << std::endl;
								}
							} else if ( other_should_store ) {
								add_to_onebody( score_map, complete_rotsets, is_asu, other_resnum, complete_rotsets->rotid_on_moltenresidue( other_node ), bonus );
								if ( TR.Trace.visible() ) {
									TR2 << boost::str(boost::format("BON: OneB Res: %i %s Rot: %i *Res: %i %s Rot: %i Score: %6.3f")
										%node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
										%other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_node )
										%bonus) << std::endl;
								}
							} else {
								// zero_body += bonus;
								if ( TR.Trace.visible() ) {
									TR2 << boost::str(boost::format("BON: ZeroB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
										%node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
										%other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_node )
										%bonus) << std::endl;
								}
							}
						}

					}
				}
			}
		}
	}

	// // for debugging cache persistence
	// if ( ! all_atoms_active ) {
	//  std::cout << "BUILDING GIANT DEBUG TABLE" << std::endl;
	//  Size scaler = Size(std::pow( 600000000, 0.25 ) );

	//  for ( Size i = 0; i <= scaler; i++ ) {
	//  for ( Size j = 0; j <= scaler; j++ ) {
	//  for ( Size k = 0; k <= scaler; k++ ) {
	//  for ( Size l = 0; l <= scaler; l++ ) {
	//   score_map->map()[{i, j, k, l}] = 5;
	//  }
	//  }
	//  }
	//  }
	// }


	milliseconds end = duration_cast< milliseconds >( system_clock::now().time_since_epoch() );
	if ( TR.Debug.visible() ) TR.Debug<< "TIMING: APPROX: " << milliseconds( end - start).count() << std::endl;




	//////////////////////////////////// Debugging info //////////////////////////////////////

	// Calculate memory use

	Size mem_use = 0;
	// http://jsteemann.github.io/blog/2016/06/14/how-much-memory-does-an-stl-container-use/
	mem_use += score_map->map().size() * ( sizeof(ResRotPair) + sizeof(float) + 20 ) + score_map->four_int_indexed_map().size() * (sizeof( std::tuple< platform::Size, platform::Size, platform::Size, platform::Size > ) + sizeof(float) + 20);
	TR << "Rough mem use: " << mem_use << " bytes" << std::endl;

	return score_map;
}

void
accumulate_bur_rotamer_oversats(
	std::unordered_map< ResRotPair, float, ResRotPairHasher > & this_bur_rotamer_oversat_map,
	std::unordered_map< ResRotPair, float, ResRotPairHasher > & this_seqpos_oversat_map
) {
	for ( auto const & pair : this_bur_rotamer_oversat_map ) {
		ResRotPair const & resrot = pair.first;
		float score = pair.second;

		// Take the max of score and whatever is in the map
		auto iter_new = this_seqpos_oversat_map.emplace( resrot, pair.second );
		if ( ! iter_new.second ) {
			if ( score > (*iter_new.first).second ) {
				(*iter_new.first).second = score;
			}
		}
	}
	this_bur_rotamer_oversat_map.clear();
}


void
accumulate_seqpos_oversats(
	basic::datacache::CacheableResRotPairFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & complete_rotsets,
	utility::vector1<bool> const & is_asu,
	ReweightData & reweight_data,
	std::unordered_map< ResRotPair, float, ResRotPairHasher > & this_seqpos_oversat_map,
	core::conformation::symmetry::SymmetryInfoCOP const & symm_info /*can be nullptr*/
) {

	// Loop over all oversat pairs that were identified
	// Each edge should be listed here at most once
	for ( auto const & pair : this_seqpos_oversat_map ) {

		ResRotPair const & oversat = pair.first;
		float max_penalty = pair.second;

		if ( oversat.first_rot != 0 && oversat.second_rot != 0 ) {
			add_to_twobody( score_map, complete_rotsets, is_asu, reweight_data,
				oversat.first_res, oversat.first_rot,
				oversat.second_res, oversat.second_rot,
				max_penalty, symm_info );
			if ( TR.Trace.visible() ) {
				TR2 << boost::str(boost::format("FINOVR: TwoB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
					%oversat.first_res%complete_rotsets->rotamer_set_for_residue( oversat.first_res )->rotamer( oversat.first_rot )->name3()%oversat.first_rot
					%oversat.second_res%complete_rotsets->rotamer_set_for_residue( oversat.second_res )->rotamer( oversat.second_rot )->name3()%oversat.second_rot
					%max_penalty) << std::endl;
			}
		} else if ( oversat.first_rot != 0 ) {
			add_to_onebody( score_map, complete_rotsets, is_asu, oversat.first_res, oversat.first_rot, max_penalty );
			if ( TR.Trace.visible() ) {
				TR2 << boost::str(boost::format("FINOVR: OneB *Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
					%oversat.first_res%complete_rotsets->rotamer_set_for_residue( oversat.first_res )->rotamer( oversat.first_rot )->name3()%oversat.first_rot
					%oversat.second_res%"BB "%oversat.second_rot
					%max_penalty) << std::endl;
			}
		} else if ( oversat.second_rot != 0 ) {
			add_to_onebody( score_map, complete_rotsets, is_asu, oversat.second_res, oversat.second_rot, max_penalty );
			if ( TR.Trace.visible() ) {
				TR2 << boost::str(boost::format("FINOVR: OneB Res: %i %s Rot: %i *Res: %i %s Rot: %i Score: %6.3f")
					%oversat.first_res%"BB "%oversat.first_rot
					%oversat.second_res%complete_rotsets->rotamer_set_for_residue( oversat.second_res )->rotamer( oversat.second_rot )->name3()%oversat.second_rot
					%max_penalty) << std::endl;
			}
		} else {
			// zero_body += max_penalty;
			if ( TR.Trace.visible() ) {
				TR2 << boost::str(boost::format("FINOVR: ZeroB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
					%oversat.first_res%"BB "%oversat.first_rot
					%oversat.second_res%"BB "%oversat.second_rot
					%max_penalty) << std::endl;
			}
		}
	}

	this_seqpos_oversat_map.clear();
}


void
add_to_onebody(
	basic::datacache::CacheableResRotPairFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const &,
	utility::vector1<bool> const & is_asu,
	Size resnum,
	Size rotamer_id,
	float adder
) {
	if ( ! is_asu[resnum] ) return;

	score_map->map()[ResRotPair(resnum, rotamer_id, 0, 0)] += adder;

	// std::cout << mm( rotamer_id, 1 ) << std::endl;

	// std::cout << "Oneb: " << resnum << " " << adder << std::endl;
}

void
add_to_twobody(
	basic::datacache::CacheableResRotPairFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & rotsets,
	utility::vector1<bool> const & is_asu,
	ReweightData & reweight_data,
	Size resnum1,
	Size rotamer_id1,
	Size resnum2,
	Size rotamer_id2,
	float adder,
	core::conformation::symmetry::SymmetryInfoCOP const & symm_info /*can be nullptr*/
) {

	if ( ! ( is_asu[resnum1] || is_asu[resnum2] ) ) return;

	if ( resnum1 == resnum2 ) {
		if ( rotamer_id1 != rotamer_id2 ) return;
		add_to_onebody( score_map, rotsets, is_asu, resnum1, rotamer_id1, adder );
		return;
	}

	bool swap = resnum1 > resnum2;
	Size r1 = swap ? resnum2 : resnum1;
	Size r2 = swap ? resnum1 : resnum2;

	Size ro1 = swap ? rotamer_id2 : rotamer_id1;
	Size ro2 = swap ? rotamer_id1 : rotamer_id2;

	float reweight = 1.0f;
	if ( reweight_data.edge_reweights ) {
		ResRotPair edge_pair( r1, 0, r2, 0 );
		if ( reweight_data.stored_edge_reweights.count( edge_pair ) == 0 ) {
			if ( reweight_data.task ) {
				reweight = reweight_data.edge_reweights->res_res_weight( reweight_data.pose, *reweight_data.task, r1, r2 );
			} else {
				reweight = reweight_data.edge_reweights->res_res_weight( reweight_data.pose, core::pack::task::PackerTask_(), r1, r2 );
			}
			reweight_data.stored_edge_reweights[ edge_pair ] = reweight;
		} else {
			reweight = reweight_data.stored_edge_reweights.at( edge_pair );
		}
	}

	float const val( adder / reweight ); // Divide by reweight!!! We are trying to cancel them

	score_map->map()[ResRotPair(r1, ro1, r2, ro2)] += val;

	// VKM, 4 April 2020:
	// If we have a symmetric pose and this is a rotamer pair involving an asymmetric unit rotamer and one of its mirrors, store
	// the value indexed by residue addresses.
	if ( symm_info != nullptr ) {
		if ( !( is_asu[ r1 ] && is_asu[ r2 ] ) ) { //Skip if they're both asymmetric unit residues.
			if ( ( symm_info->bb_is_independent(r1) && (symm_info->bb_follows(r2) == r1) ) ||
					( symm_info->bb_is_independent(r2) && (symm_info->bb_follows(r1) == r2) ) ) {
				std::tuple< platform::Size, platform::Size, platform::Size, platform::Size > const myindex( r1, ro1, r2, ro2 );
				if ( score_map->four_int_indexed_map().count( myindex ) == 0 ) {
					score_map->four_int_indexed_map()[myindex] = val;
				} else {
					score_map->four_int_indexed_map()[myindex] += val;
				}
			}
		}
	}
}


}
}
}
}
