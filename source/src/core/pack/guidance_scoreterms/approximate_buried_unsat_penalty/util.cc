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
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/hbonds/HBondGraph_util.hh>
#include <core/scoring/atomic_depth/AtomicDepth.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>

#include <utility/graph/Graph.hh>
#include <boost/format.hpp>

#include <chrono>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace approximate_buried_unsat_penalty {

static basic::Tracer TR( "core.pack.guidance_scoreterms.approximate_buried_unsat_penalty" );



typedef std::pair<Size,Size> ScratchVectorLimits;

template< class T >
class ScratchVectors {

public:
	typedef typename std::vector<T>::iterator iterator;

	ScratchVectors( Size initial_size, T default_value ) :
		default_value_( default_value ),
		next_free_( 0 ),
		data_( initial_size, default_value )
	{}

	ScratchVectorLimits new_vector( Size size ) {
		while ( next_free_ + size > data_.size() ) data_.resize( data_.size() * 10, default_value_ );
		std::pair<Size,Size> offsets( next_free_, next_free_ + size );
		next_free_ += size;
		return offsets;
	}

	iterator iter( Size offset ) {
		return data_.begin() + offset;
	}

	void clear_vectors() {
		next_free_ = 0;
	}

	bool is_data_default() {
		for ( T const & dat : data_ ) {
			if ( dat != default_value_ ) return false;
		}
		return true;
	}

	T default_value_;
	Size next_free_;
	std::vector<T> data_;
};

struct OversatToSidechain {

	OversatToSidechain(
		Size _sc1_resid, Size _sc1_rotid,
		Size _sc2_resid, Size _sc2_rotid,
		Size _bb_resid
	) :
		sc1_resid( _sc1_resid ),
		sc1_rotid( _sc1_rotid ),
		sc2_resid( _sc2_resid ),
		sc2_rotid( _sc2_rotid ),
		bb_resid( _bb_resid )
	{}

	Size sc1_resid;
	Size sc1_rotid;
	Size sc2_resid;
	Size sc2_rotid;
	Size bb_resid;

	bool
	operator<( OversatToSidechain const & ot ) const {
		if ( sc1_resid != ot.sc1_resid ) return sc1_resid < ot.sc1_resid;
		if ( sc1_rotid != ot.sc1_rotid ) return sc1_rotid < ot.sc1_rotid;
		if ( sc2_resid != ot.sc2_resid ) return sc2_resid < ot.sc2_resid;
		if ( sc2_rotid != ot.sc2_rotid ) return sc2_rotid < ot.sc2_rotid;
		if ( bb_resid  != ot.bb_resid  ) return bb_resid  < ot.bb_resid;
		return false;
	}

	bool
	operator==( OversatToSidechain const & ot ) const {
		return  ( sc1_resid == ot.sc1_resid ) &&
			( sc1_rotid == ot.sc1_rotid ) &&
			( sc2_resid == ot.sc2_resid ) &&
			( sc2_rotid == ot.sc2_rotid ) &&
			( bb_resid  == ot.bb_resid  );
	}
};

struct OversatToSidechainHasher {
	std::size_t operator()( OversatToSidechain const & k ) const {
		// Bits allocated to each field:
		// sc1_resid  13
		// sc1_rotid  13
		// sc2_resid  13
		// sc2_rotid  13
		// bb_resid   12 // sorry bb_resid

		return uint64_t( k.sc1_resid  ) << ( 13 + 13 + 13 + 12 )
			^ uint64_t( k.sc1_rotid  ) << (      13 + 13 + 12 )
			^ uint64_t( k.sc2_resid  ) << (           13 + 12 )
			^ uint64_t( k.sc2_rotid  ) << (                12 )
			^ uint64_t( k.bb_resid   );
	}
};


// This function is really complicated. There is no other way.
//
void
prestore_sc_oversat(
	std::unordered_map< OversatToSidechain, ScratchVectorLimits, OversatToSidechainHasher > & oversat_map,
	ScratchVectors<float> & oversat_scratch,
	Size resnum1, Size rotamer1, bool is_bb1,
	Size resnum2, Size rotamer2, bool is_bb2,
	Size sc_resnum, Size sc_rotamer, Size sc_num_rotamers, float score ) {

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


	OversatToSidechain sc_oversat( r1, ro1, r2, ro2, sc_resnum );


	// This is a little complicated because runtime is critical here
	// First try to emplace a fake ScratchVectorLimits into the unordered_map
	// If this succeeds, then that means we have a new key and we need to
	//   allocate a new ScratchVector and store that to the map_iterator
	// Then we use the map_iterator to get our scratch vector and store the score

	typedef std::unordered_map< OversatToSidechain, ScratchVectorLimits, OversatToSidechainHasher >::iterator map_iterator;
	typedef typename ScratchVectors<float>::iterator scratch_iterator;

	std::pair<map_iterator, bool> map_iter_new_key = oversat_map.emplace(std::piecewise_construct,
		std::forward_as_tuple( sc_oversat ),
		std::forward_as_tuple( 0,0 ));

	map_iterator & map_iter = map_iter_new_key.first;
	const bool     new_key  = map_iter_new_key.second;
	if ( new_key ) {
		ScratchVectorLimits new_limits = oversat_scratch.new_vector( sc_num_rotamers + 1 ); // +1 because backbone is 0
		(*map_iter).second = new_limits;
	}

	// Get an iterator to the start of the ScratchVector
	ScratchVectorLimits const & limits = (*map_iter).second;
	scratch_iterator scratch_vector_begin = oversat_scratch.iter( limits.first );

	// Make sure that we aren't going past the end of the allocation
	debug_assert( (int)sc_rotamer < std::distance( scratch_vector_begin, oversat_scratch.iter( limits.second ) ) );

	// Kinda like scratch_vector[ sc_rotamer ] += score;
	*( scratch_vector_begin + sc_rotamer ) += score;
}

void
accumulate_oversats(
	basic::datacache::CacheableUint64MathMatrixFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & complete_rotsets,
	ReweightData & reweight_data,
	std::unordered_map< OversatToSidechain, ScratchVectorLimits, OversatToSidechainHasher > & oversat_map,
	ScratchVectors<float> & oversat_scratch
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

basic::datacache::CacheableUint64MathMatrixFloatMapOP
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
	bool assume_const_backbone /* =true */
) {
	const float PENALTY = 1.0f;
	basic::datacache::CacheableUint64MathMatrixFloatMapOP score_map ( new basic::datacache::CacheableUint64MathMatrixFloatMap() );

	// for ( Size irotset = 1; irotset <= rotsets->nmoltenres(); irotset++ ) {
	//  Size resnum = rotsets->moltenres_2_resid(irotset);
	//  pack::rotamer_set::RotamerSetCOP const & rotset = rotsets->rotamer_set_for_residue(resnum);

	//  for ( Size irot = 1; irot <= rotset->num_rotamers(); irot++ ) {
	//   conformation::ResidueCOP rotamer = rotset->rotamer( irot );

	//   std::cout << boost::str(boost::format("Res: %i %s Rot: %i Chis: ")%resnum%rotamer->name3()%irot);
	//   for ( Size chi = 1; chi <= rotamer->nchi(); chi++ ) {
	//    std::cout << boost::str(boost::format("%7.3f ")%rotamer->chi(chi) );
	//   }
	//   std::cout << std::endl;
	//  }
	// }

	// First we need to get the AtomLevelHBondGraph and its associated complete RotamerSets

	utility::vector1<bool> position_had_rotset;
	pack::rotamer_set::RotamerSetsOP complete_rotsets_out = nullptr;

	TR << "Building hbond graph" << std::endl;

	scoring::hbonds::graph::AtomLevelHBondGraphOP hb_graph;
	hb_graph = hbonds::hbond_graph_from_partial_rotsets( pose, *rotsets, scorefxn_sc, scorefxn_bb,
		complete_rotsets_out, position_had_rotset, minimum_hb_cut );

	TR << "Hbond graph has: " << hb_graph->num_edges() << " edges requiring: " << hb_graph->getTotalMemoryUsage() << " bytes" << std::endl;

	pack::rotamer_set::RotamerSetsOP complete_rotsets =
		std::dynamic_pointer_cast<pack::rotamer_set::RotamerSets>( complete_rotsets_out );

	runtime_assert( complete_rotsets );
	runtime_assert( hb_graph->num_nodes() == complete_rotsets->nrotamers() );


	// Next, we need to create the AtomicDepth object and figure out which nodes are buried
	scoring::atomic_depth::AtomicDepth atomic_depth( pose, atomic_depth_probe_radius, true, atomic_depth_resolution );

	chemical::AtomTypeSet const & atom_type_set = pose.residue(1).type().atom_type_set();

	// Now we begin the algorithm

	ReweightData reweight_data( pose, rotsets->task() );
	std::unordered_map< OversatToSidechain, ScratchVectorLimits, OversatToSidechainHasher > oversat_map;
	// This seems like a good ballpark start
	ScratchVectors<float> oversat_scratch( rotsets->nrotamers(), 0.0f );

	//////////////////// Iterate over all polar heavy atoms looking for buried ones //////////////////////////

	Size last_resnum = 0;

	for ( Size ihbnode = 1; ihbnode <= hb_graph->num_nodes(); ihbnode++ ) {

		// These nodes are one-to-one with the complete_rotsets rotamers
		scoring::hbonds::graph::AtomLevelHBondNode * hbnode = hb_graph->get_node( ihbnode );
		runtime_assert( hbnode );
		conformation::ResidueCOP rotamer = complete_rotsets->rotamer( ihbnode );
		Size node_resnum = complete_rotsets->res_for_rotamer( ihbnode );


		// Performance hack ////////////////////////
		// If we accumulate these after each residue. We can reuse the
		//   vectors and lower our memory footprint.
		if ( node_resnum != last_resnum ) {

			accumulate_oversats( score_map, complete_rotsets, reweight_data, oversat_map, oversat_scratch );

			last_resnum = node_resnum;
		}


		// This is all the polar atoms on the rotamer
		utility::vector1< scoring::hbonds::graph::AtomInfo > const & hbnode_atoms =
			hbnode->polar_sc_atoms_not_satisfied_by_background();

		// Prepare the vector of atoms for the burial calculation
		utility::vector1< conformation::Atom > atoms;
		for ( scoring::hbonds::graph::AtomInfo const & info : hbnode_atoms ) {
			atoms.push_back( rotamer->atom( info.local_atom_id() ) );
		}

		utility::vector1<Real> depths = atomic_depth.calcdepth( atoms, atom_type_set );

		for ( Size inodeat = 1; inodeat <= hbnode_atoms.size(); inodeat++ ) {

			scoring::hbonds::graph::AtomInfo const & info = hbnode_atoms[inodeat];

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

				scoring::hbonds::graph::AtomLevelHBondEdge * hb_edge =
					static_cast< scoring::hbonds::graph::AtomLevelHBondEdge * >( *it );
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

					hbonding_rotamers.emplace_back( other_atom, other_node );
				}
			}

			//////////////////////////// Perform the scoring bookkeeping //////////////////////////////////////

			// Need to account for NH3 which needs special treatment

			Size num_hs = 0;
			if ( rotamer->heavyatom_has_polar_hydrogens( info.local_atom_id() ) ) {
				num_hs = rotamer->attached_H_end( info.local_atom_id() ) - rotamer->attached_H_begin( info.local_atom_id() ) + 1;
			}

			const float buried_penalty =  num_hs >= 3 ?  PENALTY * 3 :  PENALTY;
			const float hbond_bonus =     num_hs >= 3 ? -PENALTY * 2 : -PENALTY;
			const float oversat_penalty = num_hs >= 3 ?  PENALTY     : oversat_penalty_in;

			// First we store the buried unsat penalty to whomever owns it

			bool buried_should_store = all_atoms_active;
			buried_should_store |= position_had_rotset[ node_resnum ] && ! ( assume_const_backbone && info.is_backbone() );

			if ( buried_should_store ) {
				add_to_onebody( score_map, complete_rotsets, node_resnum, complete_rotsets->rotid_on_moltenresidue( ihbnode ), buried_penalty );
				// std::cout << boost::str(boost::format("BUR: OneB Res: %i %s Rot: %i Atom: %s Score: %6.3f")
				//  %node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
				//  %rotamer->atom_name(info.local_atom_id())%buried_penalty) << std::endl;
			} else {
				// std::cout << boost::str(boost::format("BUR: ZeroB Res: %i %s Rot: %i Atom: %s Score: %6.3f")
				//  %node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
				//  %rotamer->atom_name(info.local_atom_id())%buried_penalty) << std::endl;
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
						add_to_twobody( score_map, complete_rotsets, reweight_data,
							node_resnum, complete_rotsets->rotid_on_moltenresidue( ihbnode ),
							other_resnum, complete_rotsets->rotid_on_moltenresidue( other_rotamerid ), hbond_bonus );
					}
					// std::cout << boost::str(boost::format("HYD: TwoB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
					//  %node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
					//  %other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_rotamerid )
					//  %hbond_bonus) << std::endl;
				} else if ( buried_should_store ) {
					if ( do_store ) add_to_onebody( score_map, complete_rotsets, node_resnum, complete_rotsets->rotid_on_moltenresidue( ihbnode ), hbond_bonus );
					// std::cout << boost::str(boost::format("HYD: OneB *Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
					//  %node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
					//  %other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_rotamerid )
					//  %hbond_bonus) << std::endl;
				} else if ( other_should_store ) {
					if ( do_store ) add_to_onebody( score_map, complete_rotsets, other_resnum, complete_rotsets->rotid_on_moltenresidue( other_rotamerid ), hbond_bonus );
					// std::cout << boost::str(boost::format("HYD: OneB Res: %i %s Rot: %i *Res: %i %s Rot: %i Score: %6.3f")
					//  %node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
					//  %other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_rotamerid )
					//  %hbond_bonus) << std::endl;
				} else {
					// zero_body += hbond_bonus;
					// std::cout << boost::str(boost::format("HYD: ZeroB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
					//  %node_resnum%rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( ihbnode )
					//  %other_resnum%other_rotamer->name3()%complete_rotsets->rotid_on_moltenresidue( other_rotamerid )
					//  %hbond_bonus) << std::endl;
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

					if ( ii_resnum == jj_resnum && ii_rotamerid != jj_rotamerid ) continue;

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
							prestore_sc_oversat( oversat_map, oversat_scratch,
								ii_resnum, complete_rotsets->rotid_on_moltenresidue( ii_rotamerid ), assume_const_backbone && is_bb[ ii ],
								jj_resnum, complete_rotsets->rotid_on_moltenresidue( jj_rotamerid ), assume_const_backbone && is_bb[ jj ],
								node_resnum, complete_rotsets->rotid_on_moltenresidue( ihbnode ),
								complete_rotsets->nrotamers_for_moltenres( complete_rotsets->resid_2_moltenres( node_resnum ) ),
								oversat_penalty );
						}
					}

					if ( ii_should_store && jj_should_store ) {
						if ( do_store ) {
							add_to_twobody( score_map, complete_rotsets, reweight_data,
								ii_resnum, complete_rotsets->rotid_on_moltenresidue( ii_rotamerid ),
								jj_resnum, complete_rotsets->rotid_on_moltenresidue( jj_rotamerid ), oversat_penalty );
						}
						// std::cout << boost::str(boost::format("OVR: TwoB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
						//  %ii_resnum%complete_rotsets->rotamer( ii_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( ii_rotamerid )
						//  %jj_resnum%complete_rotsets->rotamer( jj_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( jj_rotamerid )
						//  %oversat_penalty) << std::endl;
					} else if ( ii_should_store ) {
						if ( do_store ) add_to_onebody( score_map, complete_rotsets, ii_resnum, complete_rotsets->rotid_on_moltenresidue( ii_rotamerid ), oversat_penalty );
						// std::cout << boost::str(boost::format("OVR: OneB *Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
						//  %ii_resnum%complete_rotsets->rotamer( ii_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( ii_rotamerid )
						//  %jj_resnum%complete_rotsets->rotamer( jj_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( jj_rotamerid )
						//  %oversat_penalty) << std::endl;
					} else if ( jj_should_store ) {
						if ( do_store ) add_to_onebody( score_map, complete_rotsets, jj_resnum, complete_rotsets->rotid_on_moltenresidue( jj_rotamerid ), oversat_penalty );
						// std::cout << boost::str(boost::format("OVR: OneB Res: %i %s Rot: %i *Res: %i %s Rot: %i Score: %6.3f")
						//  %ii_resnum%complete_rotsets->rotamer( ii_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( ii_rotamerid )
						//  %jj_resnum%complete_rotsets->rotamer( jj_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( jj_rotamerid )
						//  %oversat_penalty) << std::endl;
					} else {
						// zero_body += hbond_bonus;
						// std::cout << boost::str(boost::format("OVR: ZeroB Res: %i %s Rot: %i Res: %i %s Rot: %i Score: %6.3f")
						//  %ii_resnum%complete_rotsets->rotamer( ii_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( ii_rotamerid )
						//  %jj_resnum%complete_rotsets->rotamer( jj_rotamerid )->name3()%complete_rotsets->rotid_on_moltenresidue( jj_rotamerid )
						//  %oversat_penalty) << std::endl;
					}
				}
			}
		}
	}

	// Very important!!! This must be here!!! This catches the final residue.
	accumulate_oversats( score_map, complete_rotsets, reweight_data, oversat_map, oversat_scratch );


	//////////////////////////////////// Debugging info //////////////////////////////////////

	// Calculate memory use

	Size mem_use = 0;
	for ( std::pair<uint64_t, numeric::MathMatrix<float>> pair : score_map->map() ) {
		mem_use += sizeof(float) * pair.second.size();
	}

	TR << "Mem use: " << mem_use << " bytes" << std::endl;







	return score_map;
}


void
accumulate_oversats(
	basic::datacache::CacheableUint64MathMatrixFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & complete_rotsets,
	ReweightData & reweight_data,
	std::unordered_map< OversatToSidechain, ScratchVectorLimits, OversatToSidechainHasher > & oversat_map,
	ScratchVectors<float> & oversat_scratch
) {

	// Loop over all oversat pairs that were identified
	// Each edge should be listed here at most once
	for ( auto const & pair : oversat_map ) {

		OversatToSidechain const & oversat = pair.first;
		ScratchVectorLimits const & limits = pair.second;

		float max_penalty = 0;

		// Assembly level optimization right here. Find the max and clear the vector simultaneously
		for ( auto iter = oversat_scratch.iter( limits.first ), end = oversat_scratch.iter( limits.second );
				iter != end;
				++iter ) {
			max_penalty = std::max<float>( max_penalty, *iter );
			*iter = 0;
		}

		if ( oversat.sc1_rotid != 0 && oversat.sc2_rotid != 0 ) {
			add_to_twobody( score_map, complete_rotsets, reweight_data,
				oversat.sc1_resid, oversat.sc1_rotid,
				oversat.sc2_resid, oversat.sc2_rotid,
				max_penalty );
			// std::cout << boost::str(boost::format("FINOVR: TwoB Res: %i %s Rot: %i Res: %i %s Rot: %i BBRes: %i Score: %6.3f")
			//  %oversat.sc1_resid%complete_rotsets->rotamer_set_for_residue( oversat.sc1_resid )->rotamer( oversat.sc1_rotid )->name3()%oversat.sc1_rotid
			//  %oversat.sc2_resid%complete_rotsets->rotamer_set_for_residue( oversat.sc2_resid )->rotamer( oversat.sc2_rotid )->name3()%oversat.sc2_rotid
			//  %oversat.bb_resid
			//  %max_penalty) << std::endl;
		} else if ( oversat.sc1_rotid != 0 ) {
			add_to_onebody( score_map, complete_rotsets, oversat.sc1_resid, oversat.sc1_rotid, max_penalty );
			// std::cout << boost::str(boost::format("FINOVR: OneB *Res: %i %s Rot: %i Res: %i %s Rot: %i BBRes: %i Score: %6.3f")
			//  %oversat.sc1_resid%complete_rotsets->rotamer_set_for_residue( oversat.sc1_resid )->rotamer( oversat.sc1_rotid )->name3()%oversat.sc1_rotid
			//  %oversat.sc2_resid%"BB "%oversat.sc2_rotid
			//  %oversat.bb_resid
			//  %max_penalty) << std::endl;
		} else if ( oversat.sc2_rotid != 0 ) {
			add_to_onebody( score_map, complete_rotsets, oversat.sc2_resid, oversat.sc2_rotid, max_penalty );
			// std::cout << boost::str(boost::format("FINOVR: OneB Res: %i %s Rot: %i *Res: %i %s Rot: %i BBRes: %i Score: %6.3f")
			//  %oversat.sc1_resid%"BB "%oversat.sc1_rotid
			//  %oversat.sc2_resid%complete_rotsets->rotamer_set_for_residue( oversat.sc2_resid )->rotamer( oversat.sc2_rotid )->name3()%oversat.sc2_rotid
			//  %oversat.bb_resid
			//  %max_penalty) << std::endl;
		} else {
			// zero_body += max_penalty;
			// std::cout << boost::str(boost::format("FINOVR: ZeroB Res: %i %s Rot: %i Res: %i %s Rot: %i BBRes: %i Score: %6.3f")
			//  %oversat.sc1_resid%"BB "%oversat.sc1_rotid
			//  %oversat.sc2_resid%"BB "%oversat.sc2_rotid
			//  %oversat.bb_resid
			//  %max_penalty) << std::endl;
		}
	}

	oversat_scratch.clear_vectors();
	debug_assert( oversat_scratch.is_data_default() );
	oversat_map.clear();

}


void
add_to_onebody(
	basic::datacache::CacheableUint64MathMatrixFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & rotsets,
	Size resnum,
	Size rotamer_id,
	float adder
) {
	uint64_t key = map_key_oneb( resnum );
	if ( score_map->map().count( key ) == 0 ) {

		score_map->map().emplace( std::piecewise_construct, std::forward_as_tuple( key ), std::forward_as_tuple(
			// Add +1 because this is 0 offset
			rotsets->nrotamers_for_moltenres( rotsets->resid_2_moltenres( resnum ) ) + 1, 1+1, 0.0f ));
	}

	numeric::MathMatrix<float> & mm = score_map->map()[key];
	assert( rotamer_id < mm.get_number_rows() );

	mm( rotamer_id, 1 ) += adder;
	// std::cout << mm( rotamer_id, 1 ) << std::endl;

	// std::cout << "Oneb: " << resnum << " " << adder << std::endl;
}

void
add_to_twobody(
	basic::datacache::CacheableUint64MathMatrixFloatMapOP const & score_map,
	pack::rotamer_set::RotamerSetsOP const & rotsets,
	ReweightData & reweight_data,
	Size resnum1,
	Size rotamer_id1,
	Size resnum2,
	Size rotamer_id2,
	float adder
) {

	if ( resnum1 == resnum2 ) {
		if ( rotamer_id1 != rotamer_id2 ) return;
		add_to_onebody( score_map, rotsets, resnum1, rotamer_id1, adder );
		return;
	}

	bool swap = resnum1 > resnum2;
	Size r1 = swap ? resnum2 : resnum1;
	Size r2 = swap ? resnum1 : resnum2;

	Size ro1 = swap ? rotamer_id2 : rotamer_id1;
	Size ro2 = swap ? rotamer_id1 : rotamer_id2;

	uint64_t key = map_key_twob( r1, r2 );
	if ( score_map->map().count( key ) == 0 ) {

		score_map->map().emplace( std::piecewise_construct, std::forward_as_tuple( key ), std::forward_as_tuple(
			// Add +1 because this is 0 offset
			rotsets->nrotamers_for_moltenres( rotsets->resid_2_moltenres( r1 ) ) + 1,
			rotsets->nrotamers_for_moltenres( rotsets->resid_2_moltenres( r2 ) ) + 1, 0.0f ));
	}

	float reweight = 1.0f;
	if ( reweight_data.edge_reweights ) {
		if ( reweight_data.stored_edge_reweights.count( key ) == 0 ) {
			reweight = reweight_data.edge_reweights->res_res_weight( reweight_data.pose, *reweight_data.task, r1, r2 );
			reweight_data.stored_edge_reweights[ key ] = reweight;
		} else {
			reweight = reweight_data.stored_edge_reweights.at( key );
		}
	}


	numeric::MathMatrix<float> & mm = score_map->map()[key];
	assert( ro1 < mm.get_number_rows() );
	assert( ro2 < mm.get_number_cols() );

	mm( ro1, ro2 ) += adder / reweight;   // Divide by reweight!!! We are trying to cancel them
	// std::cout << mm( ro1, ro2 ) << std::endl;
	// std::cout << "Twob: " << r1 << " " << r2 << " " << adder << std::endl;
}


uint64_t
map_key_twob( Size resnum1, Size resnum2 ) {
	debug_assert( resnum1 <= resnum2 );
	debug_assert( resnum1 != 0 );

	return uint64_t(resnum1) << 32 | uint64_t(resnum2);
}

uint64_t
map_key_oneb( Size resnum1 ) {

	return uint64_t(resnum1);
}


}
}
}
}
