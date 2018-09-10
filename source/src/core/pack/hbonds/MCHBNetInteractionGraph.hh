// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/hbonds/MCHBNetInteractionGraph.hh
/// @brief Dervied class of PDInteractionGraph that does not save twobody energy calculations but rather passes them directly to a AtomLevelHBondGraph
/// @details This is a AtomLevelHBondGraph creator that is wearing an InteractionGraph disguise so that monte carlo HBNet can collect energy information without having to create custom interfaces in many other classes. This class should not be used as an InteractionGraph because it does not store all of the information that InteractionGraphs need to store. There are a few utility_exit_with_message() calls sprinkled within this class to make sure it is not being misused, but there really is not any need to use it for anything other than AtomLevelHBondGraph creation.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_pack_hbonds_MCHBNetInteractionGraph_HH
#define INCLUDED_core_pack_hbonds_MCHBNetInteractionGraph_HH

#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>
#include <core/pack/hbonds/MCHBNetInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.fwd.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>
#include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <boost/functional/hash.hpp>

#include <unordered_map>

namespace core {
namespace pack {
namespace hbonds {

class BareMinimumPDEdge : public interaction_graph::PDEdge {

public:
	BareMinimumPDEdge( interaction_graph::InteractionGraphBase* owner, int first_node_ind, int second_node_ind );
	~BareMinimumPDEdge();

public:

	///@brief look for hbonds and clashes
	void add_to_two_body_energy(int const, int const, PackerEnergy const) override;

	///@brief look for hbonds and clashes
	void add_to_two_body_energies( ObjexxFCL::FArray2< PackerEnergy > const & res_res_energy_array ) override;

	///@brief identical to add_to_two_body_energies()
	void set_two_body_energy(int const rot1, int const rot2, PackerEnergy const twobody ) override{
		add_to_two_body_energy( rot1, rot2, twobody );
	}

	///@brief Does nothing in this implementation
	void clear_two_body_energy(int const, int const) override{}

	///@brief Exits if this is called.
	PackerEnergy get_two_body_energy( int const, int const ) const override {
		utility_exit_with_message( "MCHBNetInteractionGraph should not be treated like an actual interaction graph" );
	}

	///@brief Does nothing in this implementation
	void declare_energies_final() override {}

	///@brief Does nothing in this implementation
	void prepare_for_simulated_annealing() override {}
};

class MCHBNetInteractionGraph : public interaction_graph::PDInteractionGraph {

public:

	//constructor
	MCHBNetInteractionGraph(
		scoring::hbonds::graph::AtomLevelHBondGraphOP hbond_graph,
		rotamer_set::RotamerSetsCOP rotamer_sets,
		float hbond_threshold,
		float clash_threshold );

	//destructor
	~MCHBNetInteractionGraph();

public:

	///@brief look for hbonds and clashes
	void eval_rot_pair(
		unsigned int const global_rot1,
		unsigned int const global_rot2,
		PackerEnergy const two_body_energy );

	interaction_graph::EdgeBase * create_new_edge( int index1, int index2) override{
		return new BareMinimumPDEdge( this, index1, index2 );
	}

	rotamer_set::RotamerSetsCOP rotamer_sets(){
		return rotamer_sets_;
	}

	void finalize_hbond_graph();

private:

	std::unordered_map< std::pair<uint32_t, uint32_t>, float, boost::hash< std::pair<uint32_t, uint32_t>>> future_edges_;
	scoring::hbonds::graph::AtomLevelHBondGraphOP hbond_graph_;
	rotamer_set::RotamerSetsCOP rotamer_sets_;

	float hbond_threshold_;
	float clash_threshold_;

public://methods that we never want to be called

	///@brief Exits if this is called.
	PackerEnergy get_two_body_energy_for_edge ( int, int, int, int ) const override {
		utility_exit_with_message( "MCHBNetInteractionGraph should not be treated like an actual interaction graph" );
	}

	///@brief only twobody energies will get turned into edges. If you are running with symmetry, you need to call this to find edges between resids and their symmetric twins. Must be called after scoring.
	void
	find_symmetric_hbonds(
		conformation::symmetry::SymmetryInfo const & symm_info,
		pose::Pose const & pose,
		Real hb_threshold
	);


private:

	///@brief utility function stolen from HBNet.cc. Calculates independent resid for symmetric resid.
	Size get_ind_res(
		pose::Pose const & pose,
		Size const res_i,
		conformation::symmetry::SymmetryInfo const & symm_info,
		std::map< char, std::pair< Size, Size > > & chain_bounds
	) const;

};

} //hbonds
} //pack
} //core

#endif
