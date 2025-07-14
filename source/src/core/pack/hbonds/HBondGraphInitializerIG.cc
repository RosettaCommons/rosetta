// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/hbonds/HBondGraphInitializerIG.cc
/// @brief Dervied class of PDInteractionGraph that does not save twobody energy calculations but rather passes them directly to a HBondGraph
/// @details This is a HBondGraph creator that is wearing an InteractionGraph disguise so that monte carlo HBNet can collect energy information without having to create custom interfaces in many other classes. This class should not be used as an InteractionGraph because it does not store all of the information that InteractionGraphs need to store. There are a few utility_exit_with_message() calls sprinkled within this class to make sure it is not being misused, but there really is not any need to use it for anything other than HBondGraph creation.
/// @author Jack Maguire

#include <core/conformation/Conformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/hbonds/HBondGraphInitializerIG.hh>
#include <core/pack/interaction_graph/AminoAcidNeighborSparseMatrix.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

#include <utility>

#include <core/scoring/hbonds/graph/HBondGraph.hh> // AUTO IWYU For HBondGraph, HBondEdge, HBondNode

namespace core {
namespace pack {
namespace hbonds {

static basic::Tracer TR( "core.pack.hbonds.HBondGraphInitializerIG" );

BareMinimumPDEdge::BareMinimumPDEdge( interaction_graph::InteractionGraphBase* owner, int first_node_ind, int second_node_ind ) :
	interaction_graph::PDEdge( owner, first_node_ind, second_node_ind )
{
	//Deallocate all the edge table memory
	//Not the best way to do this, but this seems to be the best option that does not alter PDEdge code.
	two_body_energies().drop_all_submatrices();
}


BareMinimumPDEdge::~BareMinimumPDEdge()= default;

void BareMinimumPDEdge::add_to_two_body_energy(int const rot1, int const rot2, PackerEnergy const twobody ){
	auto * owner = static_cast< HBondGraphInitializerIG * > ( get_owner() );
	unsigned int const offset1 = owner->rotamer_sets()->nrotamer_offset_for_moltenres( get_first_node_ind() );
	unsigned int const offset2 = owner->rotamer_sets()->nrotamer_offset_for_moltenres( get_second_node_ind() );
	owner->eval_rot_pair( offset1 + rot1, offset2 + rot2, twobody );
}

void BareMinimumPDEdge::add_to_two_body_energies( ObjexxFCL::FArray2< PackerEnergy > const & res_res_energy_array ){
	auto * owner = static_cast< HBondGraphInitializerIG * > ( get_owner() );
	unsigned int const offset1 = owner->rotamer_sets()->nrotamer_offset_for_moltenres( get_first_node_ind() );
	unsigned int const offset2 = owner->rotamer_sets()->nrotamer_offset_for_moltenres( get_second_node_ind() );
	unsigned int const nrot1 = owner->rotamer_sets()->nrotamers_for_moltenres( get_first_node_ind() );
	unsigned int const nrot2 = owner->rotamer_sets()->nrotamers_for_moltenres( get_second_node_ind() );

	for ( unsigned int rot1 = 1; rot1 <= nrot1; ++rot1 ) {
		for ( unsigned int rot2 = 1; rot2 <= nrot2; ++rot2 ) {
			owner->eval_rot_pair( offset1 + rot1, offset2 + rot2, res_res_energy_array( rot2, rot1 ) );
		}
	}
}

//Constructor
HBondGraphInitializerIG::HBondGraphInitializerIG(
	scoring::hbonds::graph::HBondGraphOP hbond_graph,
	rotamer_set::RotamerSetsCOP rotamer_sets,
	float hbond_threshold,
	float clash_threshold
) :
	interaction_graph::PDInteractionGraph( rotamer_sets->nmoltenres() ),
	hbond_graph_( std::move( hbond_graph ) ),
	rotamer_sets_( rotamer_sets ),
	hbond_threshold_( hbond_threshold ),
	clash_threshold_( clash_threshold )
{}

//Destructor
HBondGraphInitializerIG::~HBondGraphInitializerIG() = default;

void
HBondGraphInitializerIG::find_symmetric_hbonds(
	conformation::symmetry::SymmetryInfo const & symm_info,
	pose::Pose const & pose,
	Real hb_threshold
){
	std::map< char, std::pair< Size, Size > > chain_bounds;
	for ( Size ic = 1; ic <= pose.conformation().num_chains(); ++ic ) {
		Size const ic_begin = pose.conformation().chain_begin( ic );
		Size const ic_end = pose.conformation().chain_end( ic );
		char const chain = pose.chain( ic_begin );
		chain_bounds[ chain ].first = ic_begin;
		chain_bounds[ chain ].second = ic_end;
	}

	Size const scmult_1b ( symm_info.score_multiply_factor() );

	for ( Size mres = 1; mres <= rotamer_sets_->nmoltenres(); ++mres ) {
		Size const resid_raw = rotamer_sets_->moltenres_2_resid( mres );
		Size const resid = get_ind_res( pose, resid_raw, symm_info, chain_bounds );
		if ( symm_info.is_asymmetric_seqpos( resid ) ) continue;

		Size const offset_for_mres = rotamer_sets_->nrotamer_offset_for_moltenres( mres );

		for ( Size rot = 1; rot <= rotamer_sets_->nrotamers_for_moltenres( mres ); ++rot ) {
			Real const one_body_1 = get_one_body_energy_for_node_state( mres, rot ) / scmult_1b;
			if ( one_body_1 < hb_threshold ) {
				scoring::hbonds::graph::HBondEdge * new_edge =
					hbond_graph_->add_edge( offset_for_mres + rot, offset_for_mres + rot );
				new_edge->set_energy( one_body_1 );
			}
		}//rot
	}//mres

}

//Straight-up stolen from HBNet.cc
Size
HBondGraphInitializerIG::get_ind_res(
	pose::Pose const & pose,
	Size const res_i,
	conformation::symmetry::SymmetryInfo const & symm_info,
	std::map< char, std::pair< Size, Size > > & chain_bounds
) const {

	Size resi_ind( res_i );
	if ( res_i > symm_info.num_independent_residues() ) {
		char resi_chain = pose.chain( res_i );
		if ( symm_info.get_num_components() > 1 ) {
			std::map< std::string, std::pair< Size, Size > > const & component_bounds = symm_info.get_component_bounds();
			std::string resi_comp = symm_info.get_component_of_residue( res_i );
			resi_ind = res_i - chain_bounds[ resi_chain ].first + component_bounds.at( resi_comp ).first;
		} else {
			resi_ind = res_i - chain_bounds[ resi_chain ].first + 1;
		}
	}
	return resi_ind;
}


void HBondGraphInitializerIG::eval_rot_pair(
	unsigned int const global_rot1,
	unsigned int const global_rot2,
	PackerEnergy const two_body_energy
){
	using namespace scoring::hbonds::graph;

	if ( two_body_energy > clash_threshold_ ) {
		HBondNode * const node1 =
			static_cast< HBondNode * >( hbond_graph_->get_node( global_rot1 ) );
		node1->register_clash( NodeIDSize( global_rot2 ) );

		HBondNode * const node2 =
			static_cast< HBondNode * >( hbond_graph_->get_node( global_rot2 ) );
		node2->register_clash( NodeIDSize( global_rot1 ) );
		return;
	}
	if ( two_body_energy <= hbond_threshold_ ) {

		bool swap = global_rot1 > global_rot2;
		uint32_t rot1 = swap ? global_rot2 : global_rot1;
		uint32_t rot2 = swap ? global_rot1 : global_rot2;

		std::pair<uint32_t, uint32_t> key( rot1, rot2 );
#ifdef MULTI_THREADED
		utility::thread::WriteLockGuard writeguard( future_edges_mutex_ );
#endif
		if ( future_edges_.count(key) == 0 ) {
			future_edges_[ key ] = two_body_energy;
		} else {
			future_edges_[ key ] += two_body_energy;
		}
	}
}

void HBondGraphInitializerIG::finalize_hbond_graph() {

	using namespace scoring::hbonds::graph;

#ifdef MULTI_THREADED
	utility::thread::WriteLockGuard writeguard( future_edges_mutex_ );
#endif

	for ( auto const & key_val : future_edges_ ) {

		uint32_t rot1 = key_val.first.first;
		uint32_t rot2 = key_val.first.second;

		float val = key_val.second;

		HBondEdge & new_edge =
			static_cast< HBondEdge & >( * hbond_graph_->add_edge( rot1, rot2 ) );
		new_edge.set_energy( val );
	}

	future_edges_.clear();
}



} //hbonds
} //pack
} //core
