// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/hbonds/MCHBNetInteractionGraph.cc
/// @brief Dervied class of PDInteractionGraph that does not save twobody energy calculations but rather passes them directly to a AbstractHBondGraph
/// @details This is a AbstractHBondGraph creator that is wearing an InteractionGraph disguise so that monte carlo HBNet can collect energy information without having to create custom interfaces in many other classes. This class should not be used as an InteractionGraph because it does not store all of the information that InteractionGraphs need to store. There are a few utility_exit_with_message() calls sprinkled within this class to make sure it is not being misused, but there really is not any need to use it for anything other than AbstractHBondGraph creation.
/// @author Jack Maguire, jack@med.unc.edu

#include <core/pack/hbonds/MCHBNetInteractionGraph.hh>
#include <core/scoring/hbonds/graph/HBondGraph.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/interaction_graph/AminoAcidNeighborSparseMatrix.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Conformation.hh>
#include <utility>

namespace core {
namespace pack {
namespace hbonds {

BareMinimumPDEdge::BareMinimumPDEdge( interaction_graph::InteractionGraphBase* owner, int first_node_ind, int second_node_ind ) :
	interaction_graph::PDEdge( owner, first_node_ind, second_node_ind )
{
	//Deallocate all the edge table memory
	//Not the best way to do this, but this seems to be the best option that does not alter PDEdge code.
	two_body_energies().drop_all_submatrices();
}


BareMinimumPDEdge::~BareMinimumPDEdge()= default;

void BareMinimumPDEdge::add_to_two_body_energy(int const rot1, int const rot2, PackerEnergy const twobody ){
	auto * owner = static_cast< MCHBNetInteractionGraph * > ( get_owner() );
	unsigned int const offset1 = owner->rotamer_sets()->nrotamer_offset_for_moltenres( get_first_node_ind() );
	unsigned int const offset2 = owner->rotamer_sets()->nrotamer_offset_for_moltenres( get_second_node_ind() );
	owner->eval_rot_pair( offset1 + rot1, offset2 + rot2, twobody );
}

void BareMinimumPDEdge::add_to_two_body_energies( ObjexxFCL::FArray2< PackerEnergy > const & res_res_energy_array ){
	auto * owner = static_cast< MCHBNetInteractionGraph * > ( get_owner() );
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
MCHBNetInteractionGraph::MCHBNetInteractionGraph( scoring::hbonds::graph::AbstractHBondGraphOP hbond_graph, rotamer_set::RotamerSetsCOP rotamer_sets, float hbond_threshold, float clash_threshold ) :
	interaction_graph::PDInteractionGraph( rotamer_sets->nmoltenres() ),
	hbond_graph_(std::move( hbond_graph )),
	rotamer_sets_( rotamer_sets ),
	hbond_threshold_( hbond_threshold ),
	clash_threshold_( clash_threshold )
{}

//Destructor
MCHBNetInteractionGraph::~MCHBNetInteractionGraph() = default;

void
MCHBNetInteractionGraph::find_symmetric_hbonds(
	core::conformation::symmetry::SymmetryInfo const & symm_info,
	core::pose::Pose const & pose,
	core::Real hb_threshold
){

	std::map< char, std::pair< core::Size, core::Size > > chain_bounds;
	for ( Size ic = 1; ic <= pose.conformation().num_chains(); ++ic ) {
		Size const ic_begin = pose.conformation().chain_begin( ic );
		Size const ic_end = pose.conformation().chain_end( ic );
		char const chain = pose.chain( ic_begin );
		chain_bounds[ chain ].first = ic_begin;
		chain_bounds[ chain ].second = ic_end;
	}

	core::Size const scmult_1b ( symm_info.score_multiply_factor() );

	for ( core::Size mres = 1; mres <= rotamer_sets_->nmoltenres(); ++mres ) {
		core::Size const resid_raw = rotamer_sets_->moltenres_2_resid( mres );
		core::Size const resid = get_ind_res( pose, resid_raw, symm_info, chain_bounds );
		if ( symm_info.is_asymmetric_seqpos( resid ) ) continue;

		core::Size const offset_for_mres = rotamer_sets_->nrotamer_offset_for_moltenres( mres );

		for ( core::Size rot = 1; rot <= rotamer_sets_->nrotamers_for_moltenres( mres ); ++rot ) {
			core::Real const one_body_1 = get_one_body_energy_for_node_state( mres, rot ) / scmult_1b;
			if ( one_body_1 < hb_threshold ) {
				auto * new_edge = static_cast< scoring::hbonds::graph::HBondEdge * >( hbond_graph_->add_edge( offset_for_mres + rot, offset_for_mres + rot ) );
				new_edge->set_energy( one_body_1 );
			}
		}//rot
	}//mres

}

//Straight-up stolen from HBNet.cc
Size
MCHBNetInteractionGraph::get_ind_res(
	pose::Pose const & pose,
	Size const res_i,
	core::conformation::symmetry::SymmetryInfo const & symm_info,
	std::map< char, std::pair< core::Size, core::Size > > & chain_bounds
) const {

	Size resi_ind( res_i );
	if ( res_i > symm_info.num_independent_residues() ) {
		char resi_chain = pose.chain( res_i );
		if ( symm_info.get_num_components() > 1 ) {
			std::map< char, std::pair< Size, Size > > const & component_bounds = symm_info.get_component_bounds();
			char resi_comp = symm_info.get_component_of_residue( res_i );
			resi_ind = res_i - chain_bounds[ resi_chain ].first + component_bounds.at( resi_comp ).first;
		} else {
			resi_ind = res_i - chain_bounds[ resi_chain ].first + 1;
		}
	}
	return resi_ind;
}


} //hbonds
} //pack
} //core
