// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/interaction_graph/OTFFlexbbInteractionGraph.cc
/// @brief  Class implementation for OTF flexbb IG
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit headers
#include <protocols/flexpack/interaction_graph/OTFFlexbbInteractionGraph.hh>

/// Package headers
#include <protocols/flexpack/OtherContextScoreFunction.hh>

/// Project headers
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <utility/graph/Graph.hh>
#include <core/pose/Pose.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>

#include <utility/string_util.hh>

#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>


namespace protocols {
namespace flexpack {
namespace interaction_graph {

OTFFlexbbNode::OTFFlexbbNode(
	OTFFlexbbInteractionGraph * owner,
	int node_id,
	int num_states
) :
	parent( owner, node_id, num_states ),
	rotamers_( num_states ),
	rotamer_is_proline_( num_states, (unsigned char) 0 ),
	rotamer_is_glycine_( num_states, (unsigned char) 0 )
{}

OTFFlexbbNode::~OTFFlexbbNode() = default;

void
OTFFlexbbNode::print() const
{
	parent::print();
}

unsigned int
OTFFlexbbNode::count_dynamic_memory() const
{
	unsigned int total = sizeof( ResidueCOP ) * rotamers_.size();
	return parent::count_dynamic_memory() + total;
}

void
OTFFlexbbNode::set_rotamer( int state, ResidueCOP rotamer )
{
	rotamers_[ state ] = rotamer;
	rotamer_is_proline_[ state ] = rotamer->aa() == core::chemical::aa_pro;
	rotamer_is_glycine_[ state ] = rotamer->aa() == core::chemical::aa_gly;
}

/// @details compute bounding radii for the rotamers that actually exist -- maybe they're tighter than
/// the nbr_radius for the residue type... maybe I could find the minimal bounding sphere for each.
/// maybe I could keep two nbr_radii, one for all atoms and one for just the sidechain atoms?
void
OTFFlexbbNode::declare_all_rotamers_initialized()
{
	bounding_volumes_for_bb_for_aa_.resize( get_num_distinct_backbones() );
	for ( int ii = 1; ii <= get_num_distinct_backbones(); ++ii ) {
		bounding_volumes_for_bb_for_aa_[ ii ].resize( num_aa_types() );
		std::fill( bounding_volumes_for_bb_for_aa_[ ii ].begin(), bounding_volumes_for_bb_for_aa_[ ii ].end(), 0.0 );
	}

	for ( int ii = 1; ii <= get_num_distinct_backbones(); ++ii ) {
		//std::cout << get_node_index() << " ii: " << ii << std::endl;
		for ( int jj = 1; jj <= num_aa_types(); ++jj ) {
			//std::cout << get_node_index() << " ii: " << ii << " jj: " << jj <<  std::endl;
			Vector nbatm;
			Real bounding_square_radius( 0.0 );
			for ( int kk = 1; kk <= num_states_for_aa_for_bb()( jj, ii ); ++kk ) {
				//std::cout << get_node_index() << " ii: " << ii << " jj: " << jj << " kk: " << kk <<  std::endl;
				Size kkrotid = state_offsets_for_aa_for_bb()( jj, ii ) + kk;
				Residue const & kkrot = *rotamers_[ kkrotid ];
				if ( kk == 1 ) { nbatm = kkrot.xyz( kkrot.nbr_atom() ); }
				for ( Size ll = 1; ll <= kkrot.nheavyatoms(); ++ll ) {
					//std::cout << get_node_index() << " ii: " << ii << " jj: " << jj << " kk: " << kk << " ll: " << ll << std::endl;
					DistanceSquared lldis2 = nbatm.distance_squared( kkrot.xyz( ll ) );
					if ( bounding_square_radius < lldis2 ) {
						bounding_square_radius = lldis2;
					}
				}
			}
			bounding_volumes_for_bb_for_aa_[ ii ][ jj ] = bounding_square_radius;
		}
	}

	/// convert square radii to radii.
	for ( int ii = 1; ii <= get_num_distinct_backbones(); ++ii ) {
		for ( int jj = 1; jj <= num_aa_types(); ++jj ) {
			bounding_volumes_for_bb_for_aa_[ ii ][ jj ] = std::sqrt( bounding_volumes_for_bb_for_aa_[ ii ][ jj ]);
		}
	}

}

OTFFlexbbNode::Real
OTFFlexbbNode::bounding_radius_for_rotamers( int aatype, int bb ) const
{
	return bounding_volumes_for_bb_for_aa_[ bb ][ aatype ];
}


/// EDGE

OTFFlexbbEdge::OTFFlexbbEdge( OTFFlexbbInteractionGraph * owner, int node1, int node2 ) :
	parent( owner, node1, node2 ),
	// compute_bbbb_and_scbb_otf_( false ),
	lr_energies_exist_( false ),
	pose_( get_otfflexbbig_owner()->get_pose() ),
	sfxn_( get_otfflexbbig_owner()->get_scorefxn() )
{
	all_vs_bb_energy_curr_conf_[ 0 ] = all_vs_bb_energy_curr_conf_[ 1 ] = 0.0;
	procorr_curr_conf_[ 0 ] = procorr_curr_conf_[ 1 ] = 0.0;
	glycorr_curr_conf_[ 0 ] = glycorr_curr_conf_[ 1 ] = 0.0;
	scsc_energy_curr_conf_ = 0.0;

	all_vs_bb_energy_alt_conf_[ 0 ] = all_vs_bb_energy_alt_conf_[ 1 ] = 0.0;
	procorr_alt_conf_[ 0 ] = procorr_alt_conf_[ 1 ] = 0.0;
	glycorr_alt_conf_[ 0 ] = glycorr_alt_conf_[ 1 ] = 0.0;
	scsc_energy_alt_conf_ = 0.0;

	if ( nodes_part_of_same_flexseg() ) {
		sr_aa_neighbors_.dimension(
			get_otfflexbbig_owner()->get_num_aa_types(),
			get_otfflexbbig_owner()->get_num_aa_types(),
			num_bb( 0 ), 1); // num compact bb == 1
		all_vs_bb_energies_[ 0 ].dimension( get_num_states_for_node( 0 ), 1 ); // num compact bb == 1
		all_vs_bb_energies_[ 1 ].dimension( get_num_states_for_node( 1 ), 1 ); // num compact bb == 1
		procorr_energies_[ 0 ].dimension( get_num_states_for_node( 0 ), 1 );   // num compact bb == 1
		procorr_energies_[ 1 ].dimension( get_num_states_for_node( 1 ), 1 );   // num compact bb == 1
		glycorr_energies_[ 0 ].dimension( get_num_states_for_node( 0 ), 1 );   // num compact bb == 1
		glycorr_energies_[ 1 ].dimension( get_num_states_for_node( 1 ), 1 );   // num compact bb == 1

	} else {
		sr_aa_neighbors_.dimension(
			get_otfflexbbig_owner()->get_num_aa_types(),
			get_otfflexbbig_owner()->get_num_aa_types(),
			num_bb( 0 ),
			num_bb( 1 ),
			(unsigned char) 0
		);
		all_vs_bb_energies_[ 0 ].dimension( get_num_states_for_node( 0 ), num_bb( 1 ) );
		all_vs_bb_energies_[ 1 ].dimension( get_num_states_for_node( 1 ), num_bb( 0 ) );
		procorr_energies_[ 0 ].dimension( get_num_states_for_node( 0 ), num_bb( 1 ) );
		procorr_energies_[ 1 ].dimension( get_num_states_for_node( 1 ), num_bb( 0 ) );
		glycorr_energies_[ 0 ].dimension( get_num_states_for_node( 0 ), num_bb( 1 ) );
		glycorr_energies_[ 1 ].dimension( get_num_states_for_node( 1 ), num_bb( 0 ) );
	}

	sr_aa_neighbors_ = 0;
	all_vs_bb_energies_[ 0 ] = 0;
	all_vs_bb_energies_[ 1 ] = 0;
	procorr_energies_[ 0 ] = 0;
	procorr_energies_[ 1 ] = 0;
	glycorr_energies_[ 0 ] = 0;
	glycorr_energies_[ 1 ] = 0;

}

OTFFlexbbEdge::~OTFFlexbbEdge() = default;


/// @details compute sc/sc energy while node 0 (first node) is
/// considering a fixed-backbone state substitution.
/// Prerequisit: alternate state must have already been set.
OTFFlexbbEdge::PackerEnergy
OTFFlexbbEdge::compute_samebbconf_alternate_state_energy_first_node()
{
	using namespace core::scoring;

	FlexbbSparseMatrixIndex const & altinfo0( nodes_alt_info(0));
	FlexbbSparseMatrixIndex const & altinfo1( nodes_alt_info(1));

	debug_assert( nodes_cur_state(1) == nodes_alt_state(1) );
	debug_assert( nodes_cur_state(1) == get_flexbb_node(1)->current_state() );
	debug_assert( altinfo0.get_bb() == nodes_cur_info(0).get_bb() || nodes_cur_info(0).get_bb() == 0 );
	debug_assert( nodes_alt_state(0) == get_flexbb_node(0)->alternate_state() );
	debug_assert( !nodes_part_of_same_flexseg() || altinfo0.get_bb() == altinfo1.get_bb() || nodes_cur_state(0) == 0 || nodes_alt_state(1) == 0 );


	if ( nodes_cur_state(1) == 0 || nodes_alt_state( 0 ) == 0 ) {
		set_alt_energy( 0.0 );
		return 0.0;
	}

	int const bb0 = compact_bbindex( altinfo0.get_bb() );
	int const bb1 = compact_bbindex( altinfo1.get_bb() );

	all_vs_bb_energy_alt_conf_[ 0 ] = all_vs_bb_energies_[ 0 ]( nodes_alt_state(0), bb1 );
	procorr_alt_conf_[ 0 ] = state_is_proline( 1, nodes_alt_state(1) ) ?
		procorr_energies_[ 0 ]( nodes_alt_state(0), bb1 ) : 0.0;
	glycorr_alt_conf_[ 0 ] = state_is_glycine( 1, nodes_alt_state(1) ) ?
		glycorr_energies_[ 0 ]( nodes_alt_state(0), bb1 ) : 0.0;

	all_vs_bb_energy_alt_conf_[ 1 ] = all_vs_bb_energy_curr_conf_[ 1 ];
	procorr_alt_conf_[ 1 ] = state_is_proline( 0, nodes_alt_state(0) ) ?
		procorr_energies_[ 1 ]( nodes_alt_state(1), bb0 ) : 0.0;
	glycorr_alt_conf_[ 1 ] = state_is_glycine( 0, nodes_alt_state(0) ) ?
		glycorr_energies_[ 1 ]( nodes_alt_state(1), bb0 ) : 0.0;


	if ( sr_aa_neighbors_( altinfo0.get_aa_type(), altinfo1.get_aa_type(), altinfo0.get_bb(), bb1 ) ) {
		EnergyMap tbemap;
		sfxn_->eval_ci_2b_sc_sc( alt_rot(0), alt_rot(1), *pose_, tbemap );
		sfxn_->eval_cd_2b_sc_sc( alt_rot(0), alt_rot(1), *pose_, tbemap );
		scsc_energy_alt_conf_ = sfxn_->weights().dot( tbemap );
	} else {
		scsc_energy_alt_conf_ = 0.0;
	}

	if ( lr_energies_exist_ ) {
		EnergyMap emap;
		for ( auto iter = sfxn_->long_range_energies_begin(),
				iter_end = sfxn_->long_range_energies_end();
				iter != iter_end; ++iter ) {
			(*iter)->residue_pair_energy( alt_rot(0), alt_rot(1), *pose_, *sfxn_, emap );
		}
		scsc_energy_alt_conf_ += static_cast< PackerEnergy > ( sfxn_->weights().dot( emap ) );
	}

	set_alt_energy( edge_weight() * (scsc_energy_alt_conf_ +
		all_vs_bb_energy_alt_conf_[ 0 ] + all_vs_bb_energy_alt_conf_[ 1 ] +
		procorr_alt_conf_[ 0 ] + procorr_alt_conf_[ 1 ] +
		glycorr_alt_conf_[ 0 ] + glycorr_alt_conf_[ 1 ]));

#ifdef DEBUG_OTF_FLEXBB_ENERGIES_VERBOSE
	std::cout << "node 0 alt is gly?" << state_is_glycine( 0, nodes_alt_state(0) ) << std::endl;
	std::cout << "node 1 alt is gly?" << state_is_glycine( 1, nodes_alt_state(1) ) << std::endl;
	std::cout << "node 0 alt is pro?" << state_is_proline( 0, nodes_alt_state(0) ) << std::endl;
	std::cout << "node 1 alt is pro?" << state_is_proline( 1, nodes_alt_state(1) ) << std::endl;

	std::cout << "compute_samebbconf_alternate_state_energy_first_node: " << get_node_index(0) << " " << get_node_index(1) << " " << alt_energy() << " " << nodes_alt_state(0) << " " << nodes_alt_state(1) << std::endl;
	std::cout << alt_rot(0).aa() << " " << alt_rot(1).aa() << std::endl;
	std::cout << "scsc_energy_alt_conf_ " << scsc_energy_alt_conf_ <<
		" all_vs_bb_energy_alt_conf_[ 0 ] " << all_vs_bb_energy_alt_conf_[ 0 ] <<
		" all_vs_bb_energy_alt_conf_[ 1 ] " << all_vs_bb_energy_alt_conf_[ 1 ] <<
		" procorr_alt_conf_[ 0 ] " << procorr_alt_conf_[ 0 ] <<
		" procorr_alt_conf_[ 1 ] " <<  procorr_alt_conf_[ 1 ] <<
		" glycorr_alt_conf_[ 0 ] " << glycorr_alt_conf_[ 0 ] <<
		" glycorr_alt_conf_[ 1 ] " << glycorr_alt_conf_[ 1 ] << std::endl;

	std::cout << "scsc_energy_curr_conf_ " << scsc_energy_curr_conf_ <<
		" all_vs_bb_energy_curr_conf_[ 0 ] " << all_vs_bb_energy_curr_conf_[ 0 ] <<
		" all_vs_bb_energy_curr_conf_[ 1 ] " << all_vs_bb_energy_curr_conf_[ 1 ] <<
		" procorr_curr_conf_[ 0 ] " << procorr_curr_conf_[ 0 ] <<
		" procorr_curr_conf_[ 1 ] " <<  procorr_curr_conf_[ 1 ] <<
		" glycorr_curr_conf_[ 0 ] " << glycorr_curr_conf_[ 0 ] <<
		" glycorr_curr_conf_[ 1 ] " << glycorr_curr_conf_[ 1 ] << std::endl;
	std::cout << " coord: " << alt_rot(1).xyz( alt_rot(1).nheavyatoms() ).x() <<
						" " << alt_rot(1).xyz( alt_rot(1).nheavyatoms() ).y() <<
						" " << alt_rot(1).xyz( alt_rot(1).nheavyatoms() ).z() << std::endl;
#endif

	// Note: scale the interaction energy by the edge_weight held on this edge.
	return alt_energy();
}

/// @details Prerequisit: alternate state must have already been set.
OTFFlexbbEdge::PackerEnergy
OTFFlexbbEdge::compute_samebbconf_alternate_state_energy_second_node()
{
	using namespace core::scoring;

	FlexbbSparseMatrixIndex const & altinfo0( nodes_alt_info(0));
	FlexbbSparseMatrixIndex const & altinfo1( nodes_alt_info(1));

	debug_assert( nodes_cur_state(0) == nodes_alt_state(0) );
	debug_assert( nodes_cur_state(0) == get_flexbb_node(0)->current_state() );
	debug_assert( altinfo1.get_bb() == nodes_cur_info(1).get_bb() || nodes_cur_info(1).get_bb() == 0 );
	debug_assert( nodes_alt_state(1) == get_flexbb_node(1)->alternate_state() );
	debug_assert( !nodes_part_of_same_flexseg() || altinfo0.get_bb() == altinfo1.get_bb() || nodes_cur_state(0) == 0 || nodes_alt_state(1) == 0 );

	if ( nodes_cur_state(0) == 0 || nodes_alt_state(1) == 0 ) {
		set_alt_energy( 0.0 );
		return 0.0;
	}


	int const bb0 = compact_bbindex( altinfo0.get_bb() );
	int const bb1 = compact_bbindex( altinfo1.get_bb() );

	all_vs_bb_energy_alt_conf_[ 0 ] = all_vs_bb_energy_curr_conf_[ 0 ];
	procorr_alt_conf_[ 0 ] = state_is_proline( 1, nodes_alt_state(1) ) ?
		procorr_energies_[ 0 ]( nodes_alt_state(0), bb1 ) : 0.0;
	glycorr_alt_conf_[ 0 ] = state_is_glycine( 1, nodes_alt_state(1) ) ?
		glycorr_energies_[ 0 ]( nodes_alt_state(0), bb1 ) : 0.0;

	all_vs_bb_energy_alt_conf_[ 1 ] = all_vs_bb_energies_[ 1 ]( nodes_alt_state(1), bb0 );
	procorr_alt_conf_[ 1 ] = state_is_proline( 0, nodes_alt_state(0) ) ?
		procorr_energies_[ 1 ]( nodes_alt_state(1), bb0 ) : 0.0;
	glycorr_alt_conf_[ 1 ] = state_is_glycine( 0, nodes_alt_state(0) ) ?
		glycorr_energies_[ 1 ]( nodes_alt_state(1), bb0 ) : 0.0;

	if ( sr_aa_neighbors_( altinfo0.get_aa_type(), altinfo1.get_aa_type(), altinfo0.get_bb(), bb1 ) ) {
		core::scoring::EnergyMap tbemap;
		sfxn_->eval_ci_2b_sc_sc( alt_rot(0), alt_rot(1), *pose_, tbemap );
		sfxn_->eval_cd_2b_sc_sc( alt_rot(0), alt_rot(1), *pose_, tbemap );
		scsc_energy_alt_conf_ = sfxn_->weights().dot( tbemap );
	} else {
		scsc_energy_alt_conf_ = 0.0;
	}

	if ( lr_energies_exist_ ) {
		EnergyMap emap;
		for ( auto iter = sfxn_->long_range_energies_begin(),
				iter_end = sfxn_->long_range_energies_end();
				iter != iter_end; ++iter ) {
			(*iter)->residue_pair_energy( alt_rot(0), alt_rot(1), *pose_, *sfxn_, emap );
		}
		scsc_energy_alt_conf_ += static_cast< PackerEnergy > ( sfxn_->weights().dot( emap ) );
	}

	set_alt_energy( edge_weight() * (scsc_energy_alt_conf_ +
		all_vs_bb_energy_alt_conf_[ 0 ] + all_vs_bb_energy_alt_conf_[ 1 ] +
		procorr_alt_conf_[ 0 ] + procorr_alt_conf_[ 1 ] +
		glycorr_alt_conf_[ 0 ] + glycorr_alt_conf_[ 1 ]));

#ifdef DEBUG_OTF_FLEXBB_ENERGIES_VERBOSE
	std::cout << "node 0 alt is gly?" << state_is_glycine( 0, nodes_alt_state(0) ) << std::endl;
	std::cout << "node 1 alt is gly?" << state_is_glycine( 1, nodes_alt_state(1) ) << std::endl;
	std::cout << "node 0 alt is pro?" << state_is_proline( 0, nodes_alt_state(0) ) << std::endl;
	std::cout << "node 1 alt is pro?" << state_is_proline( 1, nodes_alt_state(1) ) << std::endl;

	std::cout << "compute_samebbconf_alternate_state_energy_second_node: " << get_node_index(0) << " " << get_node_index(1) << " " << alt_energy() << " " << nodes_alt_state(0) << " " << nodes_alt_state(1) << std::endl;
	std::cout << alt_rot(0).aa() << " " << alt_rot(1).aa() << std::endl;
	std::cout << "scsc_energy_alt_conf_ " << scsc_energy_alt_conf_ <<
		" all_vs_bb_energy_alt_conf_[ 0 ] " << all_vs_bb_energy_alt_conf_[ 0 ] <<
		" all_vs_bb_energy_alt_conf_[ 1 ] " << all_vs_bb_energy_alt_conf_[ 1 ] <<
		" procorr_alt_conf_[ 0 ] " << procorr_alt_conf_[ 0 ] <<
		" procorr_alt_conf_[ 1 ] " <<  procorr_alt_conf_[ 1 ] <<
		" glycorr_alt_conf_[ 0 ] " << glycorr_alt_conf_[ 0 ] <<
		" glycorr_alt_conf_[ 1 ] " << glycorr_alt_conf_[ 1 ] << std::endl;

	std::cout << "scsc_energy_curr_conf_ " << scsc_energy_curr_conf_ <<
		" all_vs_bb_energy_curr_conf_[ 0 ] " << all_vs_bb_energy_curr_conf_[ 0 ] <<
		" all_vs_bb_energy_curr_conf_[ 1 ] " << all_vs_bb_energy_curr_conf_[ 1 ] <<
		" procorr_curr_conf_[ 0 ] " << procorr_curr_conf_[ 0 ] <<
		" procorr_curr_conf_[ 1 ] " <<  procorr_curr_conf_[ 1 ] <<
		" glycorr_curr_conf_[ 0 ] " << glycorr_curr_conf_[ 0 ] <<
		" glycorr_curr_conf_[ 1 ] " << glycorr_curr_conf_[ 1 ] << std::endl;
	std::cout << " coord: " << alt_rot(0).xyz( alt_rot(0).nheavyatoms() ).x() <<
		" " << alt_rot(0).xyz( alt_rot(0).nheavyatoms() ).y() <<
		" " << alt_rot(0).xyz( alt_rot(0).nheavyatoms() ).z() << std::endl;
#endif


	return alt_energy();
}

/// @details Prerequisit: alternate state(s) must have already been set.
OTFFlexbbEdge::PackerEnergy
OTFFlexbbEdge::compute_altbbconf_alternate_state_energy()
{
	using namespace core::scoring;

	FlexbbSparseMatrixIndex const & altinfo0( nodes_alt_info(0));
	FlexbbSparseMatrixIndex const & altinfo1( nodes_alt_info(1));

	debug_assert( nodes_alt_state(0) == get_flexbb_node(0)->alternate_state() );
	debug_assert( nodes_alt_state(1) == get_flexbb_node(1)->alternate_state() );


	if ( nodes_alt_state( 0 ) == 0 || nodes_alt_state( 1 ) == 0 ) {
		set_alt_energy( 0.0 );
		return 0.0;
	}

	int const bb0 = compact_bbindex( altinfo0.get_bb() );
	int const bb1 = compact_bbindex( altinfo1.get_bb() );

	all_vs_bb_energy_alt_conf_[ 0 ] = all_vs_bb_energies_[ 0 ]( nodes_alt_state(0), bb1 );
	procorr_alt_conf_[ 0 ] = state_is_proline( 1, nodes_alt_state(1) ) ?
		procorr_energies_[ 0 ]( nodes_alt_state(0), bb1 ) : 0.0;
	glycorr_alt_conf_[ 0 ] = state_is_glycine( 1, nodes_alt_state(1) ) ?
		glycorr_energies_[ 0 ]( nodes_alt_state(0), bb1 ) : 0.0;

	all_vs_bb_energy_alt_conf_[ 1 ] = all_vs_bb_energies_[ 1 ]( nodes_alt_state(1), bb0 );
	procorr_alt_conf_[ 1 ] = state_is_proline( 0, nodes_alt_state(0) ) ?
		procorr_energies_[ 1 ]( nodes_alt_state(1), bb0 ) : 0.0;
	glycorr_alt_conf_[ 1 ] = state_is_glycine( 0, nodes_alt_state(0) ) ?
		glycorr_energies_[ 1 ]( nodes_alt_state(1), bb0 ) : 0.0;

	if ( sr_aa_neighbors_( altinfo0.get_aa_type(), altinfo1.get_aa_type(), altinfo0.get_bb(), bb1 ) ) {
		core::scoring::EnergyMap tbemap;
		sfxn_->eval_ci_2b_sc_sc( alt_rot(0), alt_rot(1), *pose_, tbemap );
		sfxn_->eval_cd_2b_sc_sc( alt_rot(0), alt_rot(1), *pose_, tbemap );
		scsc_energy_alt_conf_ = static_cast< PackerEnergy > ( sfxn_->weights().dot( tbemap ));
	} else {
		scsc_energy_alt_conf_ = 0.0;
	}

	if ( lr_energies_exist_ ) {
		EnergyMap emap;
		for ( auto iter = sfxn_->long_range_energies_begin(),
				iter_end = sfxn_->long_range_energies_end();
				iter != iter_end; ++iter ) {
			(*iter)->residue_pair_energy( alt_rot(0), alt_rot(1), *pose_, *sfxn_, emap );
		}
		scsc_energy_alt_conf_ += static_cast< PackerEnergy > ( sfxn_->weights().dot( emap ) );
	}

	set_alt_energy( edge_weight() * (scsc_energy_alt_conf_ +
		all_vs_bb_energy_alt_conf_[ 0 ] + all_vs_bb_energy_alt_conf_[ 1 ] +
		procorr_alt_conf_[ 0 ] + procorr_alt_conf_[ 1 ] +
		glycorr_alt_conf_[ 0 ] + glycorr_alt_conf_[ 1 ]));

#ifdef DEBUG_OTF_FLEXBB_ENERGIES_VERBOSE
	std::cout << "compute_altbbconf_alternate_state_energy: " << get_node_index(0) << " " << get_node_index(1) << " " << bb0 << " " << bb1 << " " << alt_energy() << " " << nodes_alt_state(0) << " " << nodes_alt_state(1) << std::endl;
	std::cout << alt_rot(0).aa() << " " << alt_rot(1).aa() << std::endl;
	std::cout << "scsc_energy_alt_conf_ " << scsc_energy_alt_conf_ <<
		" all_vs_bb_energy_alt_conf_[ 0 ] " << all_vs_bb_energy_alt_conf_[ 0 ] <<
		" all_vs_bb_energy_alt_conf_[ 1 ] " << all_vs_bb_energy_alt_conf_[ 1 ] <<
		" procorr_alt_conf_[ 0 ] " << procorr_alt_conf_[ 1 ] <<
		" procorr_alt_conf_[ 1 ] " <<  procorr_alt_conf_[ 1 ] <<
		" glycorr_alt_conf_[ 0 ] " << glycorr_alt_conf_[ 1 ] <<
		" glycorr_alt_conf_[ 1 ] " << glycorr_alt_conf_[ 1 ] << std::endl;

	std::cout << "scsc_energy_curr_conf_ " << scsc_energy_curr_conf_ <<
		" all_vs_bb_energy_curr_conf_[ 0 ] " << all_vs_bb_energy_curr_conf_[ 0 ] <<
		" all_vs_bb_energy_curr_conf_[ 1 ] " << all_vs_bb_energy_curr_conf_[ 1 ] <<
		" procorr_curr_conf_[ 0 ] " << procorr_curr_conf_[ 1 ] <<
		" procorr_curr_conf_[ 1 ] " <<  procorr_curr_conf_[ 1 ] <<
		" glycorr_curr_conf_[ 0 ] " << glycorr_curr_conf_[ 1 ] <<
		" glycorr_curr_conf_[ 1 ] " << glycorr_curr_conf_[ 1 ] << std::endl;
#endif

	return  alt_energy();

}

void
OTFFlexbbEdge::otfedge_note_substitution_accepted()
{
	//std::cout << "otfedge_note_substitution_accepted " << get_node_index(0) << " " << get_node_index(1) << std::endl;
	scsc_energy_curr_conf_           = scsc_energy_alt_conf_;
	all_vs_bb_energy_curr_conf_[ 0 ] = all_vs_bb_energy_alt_conf_[ 0 ];
	all_vs_bb_energy_curr_conf_[ 1 ] = all_vs_bb_energy_alt_conf_[ 1 ];
	procorr_curr_conf_[ 0 ]          = procorr_alt_conf_[ 0 ];
	procorr_curr_conf_[ 1 ]          = procorr_alt_conf_[ 1 ];
	glycorr_curr_conf_[ 0 ]          = glycorr_alt_conf_[ 0 ];
	glycorr_curr_conf_[ 1 ]          = glycorr_alt_conf_[ 1 ];
	/*
	std::cout << "scsc_energy_alt_conf_ " << scsc_energy_alt_conf_ <<
	" all_vs_bb_energy_alt_conf_[ 0 ] " << all_vs_bb_energy_alt_conf_[ 0 ] <<
	" all_vs_bb_energy_alt_conf_[ 1 ] " << all_vs_bb_energy_alt_conf_[ 1 ] <<
	" procorr_alt_conf_[ 0 ] " << procorr_alt_conf_[ 1 ] <<
	" procorr_alt_conf_[ 1 ] " <<  procorr_alt_conf_[ 1 ] <<
	" glycorr_alt_conf_[ 0 ] " << glycorr_alt_conf_[ 1 ] <<
	" glycorr_alt_conf_[ 1 ] " << glycorr_alt_conf_[ 1 ] << std::endl;

	std::cout << "scsc_energy_curr_conf_ " << scsc_energy_curr_conf_ <<
	" all_vs_bb_energy_curr_conf_[ 0 ] " << all_vs_bb_energy_curr_conf_[ 0 ] <<
	" all_vs_bb_energy_curr_conf_[ 1 ] " << all_vs_bb_energy_curr_conf_[ 1 ] <<
	" procorr_curr_conf_[ 0 ] " << procorr_curr_conf_[ 1 ] <<
	" procorr_curr_conf_[ 1 ] " <<  procorr_curr_conf_[ 1 ] <<
	" glycorr_curr_conf_[ 0 ] " << glycorr_curr_conf_[ 1 ] <<
	" glycorr_curr_conf_[ 1 ] " << glycorr_curr_conf_[ 1 ] << std::endl;
	*/
}

unsigned int
OTFFlexbbEdge::count_dynamic_memory() const
{

	unsigned int total = parent::count_dynamic_memory();
	total += all_vs_bb_energies_[ 0 ].size() * sizeof( PackerEnergy );
	total += all_vs_bb_energies_[ 1 ].size() * sizeof( PackerEnergy );
	total += procorr_energies_[ 0 ].size() * sizeof( PackerEnergy );
	total += procorr_energies_[ 1 ].size() * sizeof( PackerEnergy );
	total += glycorr_energies_[ 0 ].size() * sizeof( PackerEnergy );
	total += glycorr_energies_[ 1 ].size() * sizeof( PackerEnergy );
	total += sr_aa_neighbors_.size() * sizeof( unsigned char );
	return total;
}

void
OTFFlexbbEdge::set_ProCorrection_values(
	int node_not_necessarily_proline,
	int state,
	int other_bb,
	PackerEnergy bb_nonprobb_E,
	PackerEnergy bb_probb_E,
	PackerEnergy sc_nonprobb_E,
	PackerEnergy sc_probb_E
)
{
	int const node = which_node( node_not_necessarily_proline );

	all_vs_bb_energies_[ node ]( state, compact_bbindex( other_bb ) ) = sc_nonprobb_E + 0.5 * bb_nonprobb_E;

	procorr_energies_[ node ]( state, compact_bbindex( other_bb ) ) =
		sc_probb_E + 0.5 * bb_probb_E -
		(sc_nonprobb_E + 0.5 * bb_nonprobb_E);

	/*if ( get_node_index(0) == 27 && get_node_index(1) == 28 ) {
	std::cout << "SetProCorr: " << get_node_index(0) << " " << get_node_index(1) << " " << node << " " << state << " " << compact_bbindex( other_bb )
	<< " allvbb " << all_vs_bb_energies_[ node ]( state, compact_bbindex( other_bb ) ) << " "
	<< " procorr " << procorr_energies_[ node ]( state, compact_bbindex( other_bb ) ) << std::endl;
	}*/


}

void
OTFFlexbbEdge::set_GlyCorrection_values(
	int node_not_necessarily_glycine,
	int state,
	int other_bb,
	PackerEnergy bb_nonglybb_E,
	PackerEnergy bb_glybb_E,
	PackerEnergy sc_nonglybb_E,
	PackerEnergy sc_glybb_E
)
{
	int const node = which_node( node_not_necessarily_glycine );

	/// Assume this has already been done when setting the glycine correction?
	/// No.  Proline might not be present.  If both glycine and glycine are present
	/// then this assignment is done twice; that's ok.
	all_vs_bb_energies_[ node ]( state, compact_bbindex( other_bb ) ) = sc_nonglybb_E + 0.5 * bb_nonglybb_E;

	glycorr_energies_[ node ]( state, compact_bbindex( other_bb ) ) =
		sc_glybb_E + 0.5 * bb_glybb_E -
		(sc_nonglybb_E + 0.5 * bb_nonglybb_E);

	/*if ( get_node_index(0) == 27 && get_node_index(1) == 28 ) {
	std::cout << "SetGlyCorr: " << get_node_index(0) << " " << get_node_index(1) << " " << node << " " << state << " " << compact_bbindex( other_bb )
	<< " allvbb " << all_vs_bb_energies_[ node ]( state, compact_bbindex( other_bb ) ) << " "
	<< " glycorr " << glycorr_energies_[ node ]( state, compact_bbindex( other_bb ) ) << std::endl;
	}*/

}


void
OTFFlexbbEdge::prepare_for_simulated_annealing()
{
	//std::cout << "OTFFlexbbEdge::prepare_for_simulated_annealing" << std::endl;

	/// compute sr_aa_neighbors_ table
	using namespace core;

	Real const sfxn_reach = sfxn_->info()->max_atomic_interaction_distance();

	Size const naatypes = sr_aa_neighbors_.size1();
	for ( int ii = 1; ii <= num_bb( 0 ); ++ii ) {
		for ( int jj = 1, jje = compact_bbindex( num_bb( 1 ) ); jj <= jje; ++jj ) {
			int const jjbb = get_nodes_from_same_flexseg() ? ii : jj;
			for ( Size kk = 1; kk <= naatypes; ++kk ) {

				if ( get_otfflexbb_node( 0 )->num_states_for_aa_for_bb()( kk, ii ) == 0 ) continue;
				int kkoffset = get_otfflexbb_node( 0 )->state_offsets_for_aa_for_bb()( kk, ii );
				int kkrotno =  1 + kkoffset;
				Residue const & kkrot = get_otfflexbb_node( 0 )->rotamer( kkrotno );
				Vector const kkpos = kkrot.xyz( kkrot.nbr_atom() );
				Real const kk_reach = get_otfflexbb_node( 0 )->bounding_radius_for_rotamers( kk, ii );

				for ( Size ll = 1; ll <= naatypes; ++ll ) {
					if ( get_otfflexbb_node( 1 )->num_states_for_aa_for_bb()( ll, jjbb ) == 0 ) continue;

					int lloffset = get_otfflexbb_node( 1 )->state_offsets_for_aa_for_bb()( ll, jjbb );
					int llrotno =  1 + lloffset;
					Residue const & llrot = get_otfflexbb_node( 1 )->rotamer( llrotno );
					Vector const llpos = llrot.xyz( llrot.nbr_atom() );
					Real const ll_reach = get_otfflexbb_node( 1 )->bounding_radius_for_rotamers( ll, jjbb );

					// cast boolean d2 comparison to unsigned char -- also should test boolean lookup time
					// as well as integer lookup time... reading boolean FArrays is usually very slow.
					sr_aa_neighbors_( kk, ll, ii, jj ) = (unsigned char)
						(kkpos.distance_squared( llpos ) < std::pow( kk_reach + sfxn_reach + ll_reach, 2 ));
					//std::cout << "EDGE: " << get_node_index(0) << " " << get_node_index(0) << " " << ii << " ";
					//std::cout << jj << " " << kkrot.aa() << " " << llrot.aa() << " " << (bool) sr_aa_neighbors_( kk, ll, ii, jj ) << std::endl;
				}
			}
		}
	}
}

void
OTFFlexbbEdge::note_long_range_interactions_exist()
{
	lr_energies_exist_ = true;
}

inline
bool
OTFFlexbbEdge::state_is_proline( int which_node, int state ) const {
	return get_otfflexbb_node( which_node )->rotamer_is_proline( state );
}

inline
bool
OTFFlexbbEdge::state_is_glycine( int which_node, int state ) const {
	return get_otfflexbb_node( which_node )->rotamer_is_glycine( state );
}

void
OTFFlexbbEdge::zero_state_on_node( int which_node )
{

	parent::set_node_state_to_zero( which_node );
	all_vs_bb_energy_curr_conf_[ 0 ] = all_vs_bb_energy_curr_conf_[ 1 ] = 0.0;
	procorr_curr_conf_[ 0 ] = /*procorr_curr_conf_[ 0 ] =*/ 0.0;  // Was procorr_curr_conf_[ 1 ] meant? ~Labonte
	scsc_energy_curr_conf_ = 0;

	all_vs_bb_energy_alt_conf_[ 0 ] = all_vs_bb_energy_alt_conf_[ 1 ] = 0.0;
	procorr_alt_conf_[ 0 ] = /*procorr_alt_conf_[ 0 ] =*/ 0.0;  // Was procorr_alt_conf_[ 1 ] meant? ~Labonte
	scsc_energy_alt_conf_ = 0;

}

void
OTFFlexbbEdge::print_alt_energies() const
{
	std::cout << "OTFFlexbbEdge::print_alt_energies() " << get_node_index(0) << " " << get_node_index(1) << std::endl;
	std::cout << "scsc_energy_alt_conf_ " << scsc_energy_alt_conf_ <<
		" all_vs_bb_energy_alt_conf_[ 0 ] " << all_vs_bb_energy_alt_conf_[ 0 ] <<
		" all_vs_bb_energy_alt_conf_[ 1 ] " << all_vs_bb_energy_alt_conf_[ 1 ] <<
		" procorr_alt_conf_[ 0 ] " << procorr_alt_conf_[ 1 ] <<
		" procorr_alt_conf_[ 1 ] " <<  procorr_alt_conf_[ 1 ] <<
		" glycorr_alt_conf_[ 0 ] " << glycorr_alt_conf_[ 1 ] <<
		" glycorr_alt_conf_[ 1 ] " << glycorr_alt_conf_[ 1 ] << std::endl;

	std::cout << "scsc_energy_curr_conf_ " << scsc_energy_curr_conf_ <<
		" all_vs_bb_energy_curr_conf_[ 0 ] " << all_vs_bb_energy_curr_conf_[ 0 ] <<
		" all_vs_bb_energy_curr_conf_[ 1 ] " << all_vs_bb_energy_curr_conf_[ 1 ] <<
		" procorr_curr_conf_[ 0 ] " << procorr_curr_conf_[ 1 ] <<
		" procorr_curr_conf_[ 1 ] " <<  procorr_curr_conf_[ 1 ] <<
		" glycorr_curr_conf_[ 0 ] " << glycorr_curr_conf_[ 1 ] <<
		" glycorr_curr_conf_[ 1 ] " << glycorr_curr_conf_[ 1 ] << std::endl;
}


/// GRAPH

OTFFlexbbInteractionGraph::OTFFlexbbInteractionGraph( int num_nodes ) :
	parent( num_nodes ),
	current_pose_energy_( 0.0 ),
	alternate_pose_energy_( 0.0 )
{}

OTFFlexbbInteractionGraph::~OTFFlexbbInteractionGraph() = default;

void
OTFFlexbbInteractionGraph::initialize(
	core::pack_basic::RotamerSetsBase const & rot_sets
)
{
	parent::initialize( rot_sets );

	debug_assert( get_num_nodes() == (int) rot_sets.nmoltenres() );

	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		int ii_rotoffset = rot_sets.nrotamer_offset_for_moltenres( ii );
		for ( int jj = 1, jje = rot_sets.nrotamers_for_moltenres(ii);
				jj <= jje; ++jj ) {
			get_otfflexbb_node(ii)->set_rotamer( jj, rot_sets.rotamer( jj + ii_rotoffset ) );
		}
		get_otfflexbb_node(ii)->declare_all_rotamers_initialized();
	}

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	for ( int ii = 1; ii <= get_num_nodes(); ++ii ) {
		moltenres_2_resid_[ ii ] = rot_sets.moltenres_2_resid( ii );
		resid_2_moltenres_[ moltenres_2_resid_[ ii ] ] = ii;
	}
#endif

}

void
OTFFlexbbInteractionGraph::set_ProCorrection_values_for_edge(
	int node1,
	int node2,
	int node_not_necessarily_proline,
	int state,
	int other_bb,
	PackerEnergy bb_nonprobb_E,
	PackerEnergy bb_probb_E,
	PackerEnergy sc_nonprobb_E,
	PackerEnergy sc_probb_E
)
{
	find_otfflexbb_edge( node1, node2 )->set_ProCorrection_values(
		node_not_necessarily_proline, state, other_bb,
		bb_nonprobb_E, bb_probb_E, sc_nonprobb_E, sc_probb_E );
}

void
OTFFlexbbInteractionGraph::set_GlyCorrection_values_for_edge(
	int node1,
	int node2,
	int node_not_necessarily_glycine,
	int state,
	int other_bb,
	PackerEnergy bb_nonglybb_E,
	PackerEnergy bb_glybb_E,
	PackerEnergy sc_nonglybb_E,
	PackerEnergy sc_glybb_E
)
{
	find_otfflexbb_edge( node1, node2 )->set_GlyCorrection_values(
		node_not_necessarily_glycine, state, other_bb,
		bb_nonglybb_E, bb_glybb_E, sc_nonglybb_E, sc_glybb_E );
}


unsigned int
OTFFlexbbInteractionGraph::count_dynamic_memory() const
{
	/// Note: no memory accounting code for Pose or ScoreFunction!
	return parent::count_dynamic_memory();
}

void
OTFFlexbbInteractionGraph::set_pose( Pose const & pose )
{
	debug_assert( get_num_edges() == 0 );
	pose_ = PoseOP( new Pose( pose ) );

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	current_pose_ = new Pose( pose );
	alternate_pose_ = new Pose( pose );
	resid_2_moltenres_.resize( pose.size() );
	moltenres_2_resid_.resize( get_num_nodes() );
#endif

}


void
OTFFlexbbInteractionGraph::set_scorefxn( ScoreFunction const & sfxn )
{
	debug_assert( get_num_edges() == 0 );
	sfxn_ = sfxn.clone();

#ifdef DEBUG_OTF_FLEXBB_ENERGIES
	using namespace core::scoring;
	debug_assert( pose_ ); // context pose must be set before the score function.

	(*sfxn_)( *pose_ );

	oc_sfxn_ = new OtherContextScoreFunction( *pose_ );

	methods::EnergyMethodOptions opts = sfxn_->energy_method_options();
	opts.decompose_bb_hb_into_pair_energies( true );
	oc_sfxn_->set_energy_method_options( opts );

	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		ScoreType iist = (ScoreType) ii;
		if ( sfxn.weights()[ iist ] != 0.0 ) {
			oc_sfxn_->set_weight( iist, sfxn.weights()[ iist ] );
		}
	}

#endif


}

OTFFlexbbInteractionGraph::PoseCOP
OTFFlexbbInteractionGraph::get_pose() const
{
	debug_assert( pose_ );
	return pose_;
}

OTFFlexbbInteractionGraph::ScoreFunctionCOP
OTFFlexbbInteractionGraph::get_scorefxn() const
{
	debug_assert( sfxn_ );
	return sfxn_;
}

void
OTFFlexbbInteractionGraph::note_long_range_interactions_exist_for_edge(
	int node1,
	int node2
)
{
	OTFFlexbbEdge * edge = find_otfflexbb_edge( node1, node2 );
	if ( edge ) {
		edge->note_long_range_interactions_exist();
	}
}


void
OTFFlexbbInteractionGraph::debug_note_considered_substitution( core::conformation::Residue const & alt_rotamer, int index )
{
#ifndef DEBUG_OTF_FLEXBB_ENERGIES
	utility_exit_with_message( "Do not call OTFFlexbbInteractionGraph::debug_note_considered_substitution unless debugging" );
#endif

	//std::cout << "Replacing residue " << alt_rotamer.seqpos() << "  moltenres: " << resid_2_moltenres_[alt_rotamer.seqpos() ] << std::endl;
	alternate_pose_->replace_residue( alt_rotamer.seqpos(), alt_rotamer, false ); // PHIL PLEASE MAKE THIS O(1)!
	changing_seqpos_.push_back( alt_rotamer.seqpos() );
	alt_rots_.push_back( alt_rotamer.clone() );
	alt_rot_inds_.push_back( index );
}

void
OTFFlexbbInteractionGraph::debug_note_projected_deltaE_of_considered_substitution(
	PackerEnergy deltaE,
	PackerEnergy node_totalE,
	bool require_match
)
{
#ifndef DEBUG_OTF_FLEXBB_ENERGIES
	utility_exit_with_message( "Do not call OTFFlexbbInteractionGraph::debug_note_projected_deltaE_of_considered_substitution unless debugging" );
#endif
	Real alt_E_real = 0;// = (*sfxn_)(*alternate_pose_);
	alternate_pose_energy_ = (*oc_sfxn_)( *alternate_pose_ );

	static int n_correct_since_last_problem = 0;

	Real real_deltaE = alternate_pose_energy_ - current_pose_energy_;
	Real delta_delta = std::abs( deltaE - real_deltaE );
	bool large = std::abs(delta_delta) > 1e-4;
	bool significant = large && std::abs( delta_delta / node_totalE ) > 1e-5 && std::abs( delta_delta / (deltaE + node_totalE )) > 1e-5 && std::abs( delta_delta / real_deltaE ) > 1e-5;

	if ( require_match && ! any_vertex_state_unassigned() && significant && large ) {

		std::cout << "Delta E: predicted -- " << deltaE << " real: "
			<< real_deltaE << " delta_delta: " << delta_delta << " sig: " << std::abs( delta_delta / node_totalE ) << " " <<  std::abs( delta_delta / (deltaE + node_totalE ))
			<< " (since last problem: " << n_correct_since_last_problem << " )"
			<< " ocE: " << alternate_pose_energy_ << " realE " << alt_E_real << " delta: " << alternate_pose_energy_ - alt_E_real  <<  std::endl;
		n_correct_since_last_problem = 0;
		std::cout << "Changing nodes:";
		for ( Size ii = 1; ii <= changing_seqpos_.size(); ++ii ) {
			std::cout << " ( " << changing_seqpos_[ ii ] << " " << alt_rot_inds_[ ii ] << ")";
		}
		std::cout << std::endl;
		for ( Size ii = 1; ii <= changing_seqpos_.size(); ++ii ) {
			Size ii_resid = changing_seqpos_[ ii ];
			Size ii_moltenresid = resid_2_moltenres_[ ii_resid ];
			//std::cout << "Real energies: " << ii_moltenresid;
			Real ii_total( 0.0 );
			for ( utility::graph::Graph::EdgeListIter
					ir  = alternate_pose_->energies().energy_graph().get_node(ii_resid)->edge_list_begin(),
					ire = alternate_pose_->energies().energy_graph().get_node(ii_resid)->edge_list_end();
					ir != ire; ++ir ) {
				Size jj_resid = (*ir)->get_other_ind( ii_resid );
				Size jj_moltenresid = resid_2_moltenres_[ jj_resid ];
				if ( jj_moltenresid == 0 ) continue;

				Real ii_jj_energy = (static_cast< core::scoring::EnergyEdge const *> (*ir))->dot( oc_sfxn_->weights() );
				ii_total += ii_jj_energy;

				//std::cout << "Two-Body Energy Real: " << ii_moltenresid << " " << jj_moltenresid << " " << ii_jj_energy << std::endl;

				OTFFlexbbEdge const * ii_jj_edge = find_otfflexbb_edge( ii_moltenresid, jj_moltenresid );
				if ( !ii_jj_edge  ) {
					if ( ii_jj_energy != 0.0 ) std::cout << "EDGE MISSING FROM INTERACTION GRAPH: ( " << ii_moltenresid << " " << jj_moltenresid << " ) missing energy of " << ii_jj_energy << std::endl;
					continue;
				}

				Real delta = std::abs( ii_jj_edge->alt_energy() - ii_jj_energy );
				if ( delta > 1e-5 ) {

					std::cout << "Two-Body Energy Error: " << ii_moltenresid << " " << jj_moltenresid << " real " << ii_jj_energy ;
					std::cout << " predicted: " << ii_jj_edge->alt_energy() << " ";
					std::cout << alternate_pose_->residue( ii_resid ).aa() << " " << alternate_pose_->residue( jj_resid ).aa() << std::endl;
					std::cout << " bbbb: " << core::pack::rotamer_set::RotamerSets::get_bb_bbE(
						*alternate_pose_, *oc_sfxn_, alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) );
					std::cout << " sc1bb2: " << core::pack::rotamer_set::RotamerSets::get_sc_bbE(
						*alternate_pose_, *oc_sfxn_, alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) );
					std::cout << " bb1sc2: " << core::pack::rotamer_set::RotamerSets::get_sc_bbE(
						*alternate_pose_, *oc_sfxn_, alternate_pose_->residue( jj_resid ), alternate_pose_->residue( ii_resid ) );

					core::scoring::EnergyMap tbemap;
					oc_sfxn_->eval_ci_2b_sc_sc( alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) , *alternate_pose_, tbemap );
					oc_sfxn_->eval_cd_2b_sc_sc( alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) , *alternate_pose_, tbemap );
					Real scsc_energy_alt_conf = oc_sfxn_->weights().dot( tbemap );
					std::cout << " scsc: " << scsc_energy_alt_conf;

					std::cout << std::endl;
					ii_jj_edge->print_alt_energies();

				}

			}
			Real ii_one_body_energy = alternate_pose_->energies().onebody_energies( ii_resid ).dot( oc_sfxn_->weights() );
			Real ii_one_body_and_background = ii_one_body_energy;
			//std::cout << "Real one body energies: " <<  << std::endl;
			//std::cout << "Real two body energy total: " << ii_total << std::endl;
			//std::cout << std::endl;


			//std::cout << "Real energies with background: " << ii_moltenresid;
			ii_total = 0.0;
			for ( utility::graph::Graph::EdgeListIter
					ir  = alternate_pose_->energies().energy_graph().get_node(ii_resid)->edge_list_begin(),
					ire = alternate_pose_->energies().energy_graph().get_node(ii_resid)->edge_list_end();
					ir != ire; ++ir ) {
				Size jj_resid = (*ir)->get_other_ind( ii_resid );
				Size jj_moltenresid = resid_2_moltenres_[ jj_resid ];
				if ( jj_moltenresid != 0 ) continue;
				Real ii_jj_energy = (static_cast< core::scoring::EnergyEdge const *> (*ir))->dot(oc_sfxn_->weights());
				ii_total += ii_jj_energy;
				ii_one_body_and_background += ii_jj_energy;
				//std::cout << " " << jj_resid << " " << ii_jj_energy ;

				//std::cout << " bbbb: " << core::pack::rotamer_set::RotamerSets::get_bb_bbE(
				// *alternate_pose_, *oc_sfxn_, alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) );
				//std::cout << " sc1bb2: " << core::pack::rotamer_set::RotamerSets::get_sc_bbE(
				// *alternate_pose_, *oc_sfxn_, alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) );
				//std::cout << " bb1sc2: " << core::pack::rotamer_set::RotamerSets::get_sc_bbE(
				// *alternate_pose_, *oc_sfxn_, alternate_pose_->residue( jj_resid ), alternate_pose_->residue( ii_resid ) );

				//core::scoring::EnergyMap tbemap;
				//oc_sfxn_->eval_ci_2b_sc_sc( alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) , *alternate_pose_, tbemap );
				//oc_sfxn_->eval_cd_2b_sc_sc( alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) , *alternate_pose_, tbemap );
				//Real scsc_energy_alt_conf = oc_sfxn_->weights().dot( tbemap );
				//std::cout << " scsc: " << scsc_energy_alt_conf;

				//std::cout << std::endl;


			}

			//std::cout << "Real one body energies: " << alternate_pose_->energies().onebody_energies( ii_resid ).dot( oc_sfxn_->weights() ) << std::endl;
			//std::cout << "Real two body background energy total: " << ii_total << std::endl;
			//std::cout << "Real total energy: " << alternate_pose_->energies().residue_total_energies( ii_resid ).dot( oc_sfxn_->weights() ) << std::endl;
			Real one_body_delta = std::abs( ii_one_body_and_background - get_otfflexbb_node( ii_moltenresid )->alternate_state_one_body_energy() );
			if ( one_body_delta > 1e-5 && one_body_delta / ii_one_body_and_background > 1e-5 ) {
				std::cout << "One body energy in error: molt: " << ii_moltenresid << " resid: " << ii_resid
					<< " rot: " << alt_rot_inds_[ ii ] << " " << alternate_pose_->residue( ii_resid ).aa()
					<< " real " << ii_one_body_and_background << " pred "
					<< get_otfflexbb_node( ii_moltenresid )->alternate_state_one_body_energy() << std::endl;
				std::cout << "One body internal " << ii_one_body_energy << " ";
				alternate_pose_->energies().onebody_energies( ii_resid ).show_weighted( std::cout, oc_sfxn_->weights() );
				std::cout << std::endl;
				/*for ( Size jj = 1; jj <= alt_rots_[ ii ]->mainchain_torsions().size(); ++jj ) {
				std::cout << "bb angles: " << jj << " " <<
				alternate_pose_->residue( ii_resid ).mainchain_torsions()[ jj ] << " vs " <<
				alt_rots_[ ii ]->mainchain_torsions()[ jj ] << std::endl;
				}
				for ( Size jj = 1; jj <= alt_rots_[ ii ]->natoms(); ++jj ) {
				std::cout << "coords: " << jj << " (" << alternate_pose_->residue( ii_resid ).xyz( jj ).x()
				<< " " << alternate_pose_->residue( ii_resid ).xyz( jj ).y()
				<< " " << alternate_pose_->residue( ii_resid ).xyz( jj ).z() << ") vs ("
				<< " " << alt_rots_[ ii ]->xyz( jj ).x()
				<< " " << alt_rots_[ ii ]->xyz( jj ).y()
				<< " " << alt_rots_[ ii ]->xyz( jj ).z()
				<< ")" <<  std::endl;
				}*/
				for ( utility::graph::Graph::EdgeListIter
						ir  = alternate_pose_->energies().energy_graph().get_node(ii_resid)->edge_list_begin(),
						ire = alternate_pose_->energies().energy_graph().get_node(ii_resid)->edge_list_end();
						ir != ire; ++ir ) {
					Size jj_resid = (*ir)->get_other_ind( ii_resid );
					Size jj_moltenresid = resid_2_moltenres_[ jj_resid ];
					if ( jj_moltenresid != 0 ) continue;
					Real ii_jj_energy = (static_cast< core::scoring::EnergyEdge const *> (*ir))->dot( oc_sfxn_->weights() );
					std::cout << ii_resid << " " << jj_resid << " " << ii_jj_energy << std::endl;;

					//std::cout << " bbbb: " << core::pack::rotamer_set::RotamerSets::get_bb_bbE(
					// *alternate_pose_, *oc_sfxn_, alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) );
					//std::cout << " sc1bb2: " << core::pack::rotamer_set::RotamerSets::get_sc_bbE(
					// *alternate_pose_, *oc_sfxn_, alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) );
					//std::cout << " bb1sc2: " << core::pack::rotamer_set::RotamerSets::get_sc_bbE(
					// *alternate_pose_, *oc_sfxn_, alternate_pose_->residue( jj_resid ), alternate_pose_->residue( ii_resid ) );

					//core::scoring::EnergyMap tbemap;
					//oc_sfxn_->eval_ci_2b_sc_sc( alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) , *alternate_pose_, tbemap );
					//oc_sfxn_->eval_cd_2b_sc_sc( alternate_pose_->residue( ii_resid ), alternate_pose_->residue( jj_resid ) , *alternate_pose_, tbemap );
					//Real scsc_energy_alt_conf = oc_sfxn_->weights().dot( tbemap );
					//std::cout << " scsc: " << scsc_energy_alt_conf;

					//std::cout << std::endl;


				}

			}

		}


		//for ( Size ii = 1; ii <= moltenres_2_resid_.size(); ++ii ) {
		//std::cout << "moltenres: " << ii << " " << moltenres_2_resid_[ ii ] << std::endl;
		//}
		static int count_bad( 0 );
		alternate_pose_->dump_pdb( "BadPredDeltaE_" + utility::to_string( ++count_bad ) + "_alternate.pdb" );
		current_pose_->dump_pdb( "BadPredDeltaE_" + utility::to_string( count_bad ) + "_current.pdb" );
		debug_assert( ! require_match ||  any_vertex_state_unassigned() || ! ( significant && large) );
		utility_exit_with_message( "Bad predicted deltaE" );
	} else {
		++n_correct_since_last_problem;
	}
}

void
OTFFlexbbInteractionGraph::debug_note_accepted_substitution()
{
#ifndef DEBUG_OTF_FLEXBB_ENERGIES
	utility_exit_with_message( "Do not call OTFFlexbbInteractionGraph::debug_note_accepted_substitution unless debugging" );
#endif
	(*current_pose_) = (*alternate_pose_);
	current_pose_energy_ = alternate_pose_energy_;
	changing_seqpos_.clear();
	alt_rots_.clear();
	alt_rot_inds_.clear();
	//std::cout << "ACCEPTED SUBSTITUTION" << std::endl;
}

void
OTFFlexbbInteractionGraph::debug_note_rejected_substitution()
{
#ifndef DEBUG_OTF_FLEXBB_ENERGIES
	utility_exit_with_message( "Do not call OTFFlexbbInteractionGraph::debug_note_rejected_substitution unless debugging" );
#endif
	(*alternate_pose_) = (*current_pose_);
	changing_seqpos_.clear();
	alt_rots_.clear();
	alt_rot_inds_.clear();
	//std::cout << "REJECTED SUBSTITUTION" << std::endl;
}


}
}
}


