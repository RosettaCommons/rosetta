// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/hbnet/MCHBNetInteractionGraph.cc
/// @brief Dervied class of PDInteractionGraph that does not save twobody energy calculations but rather passes them directly to a HBondGraph
/// @details This is a HBondGraph creator that is wearing an InteractionGraph disguise so that monte carlo HBNet can collect energy information without having to create custom interfaces in many other classes. This class should not be used as an InteractionGraph because it does not store all of the information that InteractionGraphs need to store. There are a few utility_exit_with_message() calls sprinkled within this class to make sure it is not being misused, but there really is not any need to use it for anything other than HBondGraph creation.
/// @author Jack Maguire, jack@med.unc.edu

#include <protocols/hbnet/MCHBNetInteractionGraph.hh>
#include <protocols/hbnet/HBondGraph.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/interaction_graph/AminoAcidNeighborSparseMatrix.hh>

namespace protocols {
namespace hbnet {

BareMinimumPDEdge::BareMinimumPDEdge( core::pack::interaction_graph::InteractionGraphBase* owner, int first_node_ind, int second_node_ind ) :
	core::pack::interaction_graph::PDEdge( owner, first_node_ind, second_node_ind )
{
	//Deallocate all the edge table memory
	//Not the best way to do this, but this seems to be the best option that does not alter PDEdge code.
	two_body_energies().drop_all_submatrices();
}


BareMinimumPDEdge::~BareMinimumPDEdge(){}

void BareMinimumPDEdge::add_to_two_body_energy(int const rot1, int const rot2, core::PackerEnergy const twobody ){
	MCHBNetInteractionGraph * owner = static_cast< MCHBNetInteractionGraph * > ( get_owner() );
	unsigned int const offset1 = owner->rotamer_sets()->nrotamer_offset_for_moltenres( get_first_node_ind() );
	unsigned int const offset2 = owner->rotamer_sets()->nrotamer_offset_for_moltenres( get_second_node_ind() );
	owner->eval_rot_pair( offset1 + rot1, offset2 + rot2, twobody );
}

void BareMinimumPDEdge::add_to_two_body_energies( ObjexxFCL::FArray2< core::PackerEnergy > const & res_res_energy_array ){
	MCHBNetInteractionGraph * owner = static_cast< MCHBNetInteractionGraph * > ( get_owner() );
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
MCHBNetInteractionGraph::MCHBNetInteractionGraph( HBondGraphOP hbond_graph, core::pack::rotamer_set::RotamerSetsOP rotamer_sets, float hbond_threshold, float clash_threshold ) :
	core::pack::interaction_graph::PDInteractionGraph( rotamer_sets->nmoltenres() ),
	hbond_graph_( hbond_graph ),
	rotamer_sets_( rotamer_sets ),
	hbond_threshold_( hbond_threshold ),
	clash_threshold_( clash_threshold )
{}

//Destructor
MCHBNetInteractionGraph::~MCHBNetInteractionGraph()
{}

} //hbnet
} //protocols
