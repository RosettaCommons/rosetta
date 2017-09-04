// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/hbnet/MCHBNetInteractionGraph.hh
/// @brief Dervied class of PDInteractionGraph that does not save twobody energy calculations but rather passes them directly to a HBondGraph
/// @details This is a HBondGraph creator that is wearing an InteractionGraph disguise so that monte carlo HBNet can collect energy information without having to create custom interfaces in many other classes. This class should not be used as an InteractionGraph because it does not store all of the information that InteractionGraphs need to store. There are a few utility_exit_with_message() calls sprinkled within this class to make sure it is not being misused, but there really is not any need to use it for anything other than HBondGraph creation.
/// @author Jack Maguire, jack@med.unc.edu

#ifndef INCLUDED_protocols_hbnet_MCHBNetInteractionGraph_HH
#define INCLUDED_protocols_hbnet_MCHBNetInteractionGraph_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/hbnet/MCHBNetInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <protocols/hbnet/HBondGraph.hh>

namespace protocols {
namespace hbnet {

class BareMinimumPDEdge : public core::pack::interaction_graph::PDEdge {

public:
	BareMinimumPDEdge( core::pack::interaction_graph::InteractionGraphBase* owner, int first_node_ind, int second_node_ind );
	~BareMinimumPDEdge();

public:

	///@brief look for hbonds and clashes
	void add_to_two_body_energy(int const, int const, core::PackerEnergy const) override;

	///@brief look for hbonds and clashes
	void add_to_two_body_energies( ObjexxFCL::FArray2< core::PackerEnergy > const & res_res_energy_array ) override;

	///@brief identical to add_to_two_body_energies()
	inline void set_two_body_energy(int const rot1, int const rot2, core::PackerEnergy const twobody ) override{
		add_to_two_body_energy( rot1, rot2, twobody );
	}

	///@brief Does nothing in this implementation
	void clear_two_body_energy(int const, int const) override{}

	///@brief Exits if this is called.
	core::PackerEnergy get_two_body_energy( int const, int const ) const override {
		utility_exit_with_message( "MCHBNetInteractionGraph should not be treated like an actual interaction graph" );
	}

	///@brief Does nothing in this implementation
	void declare_energies_final() override {}

	///@brief Does nothing in this implementation
	void prepare_for_simulated_annealing() override {}
};

class MCHBNetInteractionGraph : public core::pack::interaction_graph::PDInteractionGraph {

public:

	//constructor
	MCHBNetInteractionGraph( HBondGraphOP hbond_graph, core::pack::rotamer_set::RotamerSetsOP rotamer_sets, float hbond_threshold, float clash_threshold );

	//destructor
	~MCHBNetInteractionGraph();

public:

	///@brief look for hbonds and clashes
	inline void eval_rot_pair( unsigned int const global_rot1, unsigned int const global_rot2, core::PackerEnergy const two_body_energy );

	core::pack::interaction_graph::EdgeBase * create_new_edge( int index1, int index2) override{
		return new BareMinimumPDEdge( this, index1, index2 );
	}

	inline core::pack::rotamer_set::RotamerSetsOP rotamer_sets(){ return rotamer_sets_; }

private:

	HBondGraphOP hbond_graph_;
	core::pack::rotamer_set::RotamerSetsOP rotamer_sets_;

	float hbond_threshold_;
	float clash_threshold_;

public://methods that we never want to be called

	///@brief Exits if this is called.
	core::PackerEnergy get_two_body_energy_for_edge ( int, int, int, int ) const override {
		utility_exit_with_message( "MCHBNetInteractionGraph should not be treated like an actual interaction graph" );
	}

};

inline void MCHBNetInteractionGraph::eval_rot_pair( unsigned int const global_rot1, unsigned int const global_rot2, core::PackerEnergy const two_body_energy ){
	if ( two_body_energy > clash_threshold_ ) {
		HBondNode * const node1 = static_cast< HBondNode * >( hbond_graph_->get_node( global_rot1 ) );
		node1->register_clash( global_rot2 );

		HBondNode * const node2 = static_cast< HBondNode * >( hbond_graph_->get_node( global_rot2 ) );
		node2->register_clash( global_rot1 );
		return;
	}
	if ( two_body_energy <= hbond_threshold_ ) {
		HBondEdge * new_edge = static_cast< HBondEdge * >( hbond_graph_->add_edge( global_rot1, global_rot2 ) );
		new_edge->set_energy( two_body_energy );
		return;
	}
}

} //hbnet
} //protocols

#endif
