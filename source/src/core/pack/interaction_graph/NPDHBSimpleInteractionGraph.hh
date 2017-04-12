// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/NPDHBSimpleInteractionGraph.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_pack_interaction_graph_NPDHBSimpleInteractionGraph_HH
#define INCLUDED_core_pack_interaction_graph_NPDHBSimpleInteractionGraph_HH

// Unit headers
#include <core/pack/interaction_graph/NPDHBSimpleInteractionGraph.fwd.hh>

// Package headers
#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>
#include <core/pack/interaction_graph/NPDHBondInteractionGraph.fwd.hh>

// Project headers
#include <core/types.hh>
#include <utility/graph/Graph.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/hbonds/NPDHBondSet.fwd.hh>
#include <core/scoring/hbonds/HBondOptions.fwd.hh>
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>
#include <core/scoring/hbonds/types.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <list>

namespace core {
namespace pack {
namespace interaction_graph {

class NPDHBSimpleNode : public SimpleNode
{
public:
	typedef SimpleNode parent;

public:

	//  SimpleNode::SimpleNode();
	NPDHBSimpleNode( utility::graph::Graph * owner, Size resnum );
	virtual ~NPDHBSimpleNode();

	void copy_from( utility::graph::Node const * ) override {}

	void print() const override {}

	platform::Size count_static_memory() const override {return 1;}
	platform::Size count_dynamic_memory() const override {return 1;}

	/// @brief Copy the alternate_residue_ pointer to the current_residue_ pointer;
	/// copy the alternate energies to the current energies for this node and its
	/// incident edges.
	void
	commit_change_npd();

	/// @brief Copy the alternate energies to the current energies for this node
	/// and its incident edges.
	void
	commit_change_no_res_pointer_update();

	/// @brief Reset state on this node so that its incident edges will no longer
	/// think it is considering an alternate conformation.
	void
	reject_change_npd( conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache );

	/// @brief Compute the total energy for the current state assignment
	Real get_curr_upper_hbond_energies();

	/// @brief Compute the change in NPDHBond energy for the given (previously supplied!)
	/// rotamer substitution
	Real get_npdhb_deltaE_for_substitution();

	void prepare_for_neighbors_substitution( NPDHBSimpleNode * changing_nbr );
	void find_hbs_for_nbrs_alt_state_step1( NPDHBSimpleNode * changing_nbr );
	Real find_hbs_for_nbrs_alt_state_step2( NPDHBSimpleNode * changing_nbr );
	void acknowledge_neighbors_substitution();

	void compute_alt_weights_for_hbonds( bool curr_state );

	void reset_hbs();
	void find_curr_hbonds_to_upper_neighbors();
	void calc_curr_hbond_weights();

protected:
	NPDHBSimpleNode * npdhb_cast( utility::graph::Node * ) const;
	NPDHBSimpleNode const * npdhb_cast( utility::graph::Node const * ) const;
	NPDHBSimpleEdge * npdhb_cast( utility::graph::Edge * ) const;
	NPDHBSimpleEdge const * npdhb_cast( utility::graph::Edge const * ) const;
	NPDHBSimpleInteractionGraph const * npdhb_owner() const;
	NPDHBSimpleInteractionGraph * npdhb_owner();

private:

	void
	initialize();

private:
	Size seqpos_;

	utility::vector1< utility::vector1< NPDHBondOP > > curr_atom_hbonds_;
	utility::vector1< utility::vector1< NPDHBondOP > > alt_atom_hbonds_;

	utility::vector1< NPDHBondOP > curr_hbonds_;
	utility::vector1< NPDHBondOP > alt_hbonds_;

	utility::vector1< Real > tmp_energies_;
	utility::vector1< Real > tmp_weights_;


}; //NPDHBSimpleNode


class NPDHBSimpleEdge : public SimpleEdge {
public:
	typedef SimpleEdge parent;

public:

	//SimpleEdge( utility::graph::Graph* owner );
	NPDHBSimpleEdge( utility::graph::Graph* owner, Size res1, Size res2 );
	virtual ~NPDHBSimpleEdge();

	void copy_from( utility::graph::Edge const * ) override {}

	platform::Size count_static_memory() const override {return 1;}
	platform::Size count_dynamic_memory() const override {return 1;}

	//void commit_change() override;

private:

	NPDHBSimpleInteractionGraph *
	get_npdhb_simple_ig_owner();

	NPDHBSimpleInteractionGraph const *
	get_npdhb_simple_ig_owner() const;


private:

}; //SimpleEdge


/// @brief A simple graph class for calculating pairwise decomposable
/// energies as sidechains are moving on a fixed backbone.  This class
/// is responsible for calculating energy changes, but is passive about
/// how the sidechains are changing.  There are two main ways to drive
/// the graph: one where the graph ACTIVELY takes charge of updating pointers
/// to the sidechains, where, each external change of one pointer
/// triggers an update to the energies; and a second, where the graph
/// is PASSIVE wrt the pointers, and they must be maintained by
/// an external driver.
class NPDHBSimpleInteractionGraph : public SimpleInteractionGraph {
public:
	typedef SimpleInteractionGraph parent;

public:

	NPDHBSimpleInteractionGraph();
	virtual ~NPDHBSimpleInteractionGraph();

	void set_scorefunction( scoring::ScoreFunction const & sfxn ) override;

	/// @brief Initialization where the graph adds its own edges
	void initialize( pose::Pose const & pose) override;

	void set_pose_no_initialize( pose::Pose const & pose ) override;

	void setup_after_edge_addition() override;

	void commit_change( Size node_id ) override;
	void reject_change( Size node_id, conformation::ResidueCOP res, basic::datacache::BasicDataCache & residue_data_cache ) override;

	scoring::hbonds::HBondDatabase const & hbond_database() const;
	scoring::hbonds::HBondOptions const & hbond_options() const;
	scoring::hbonds::NPDHBondSet const & npd_hbond_set() const;

	//returns delta-energy
	Real consider_substitution(
		Size node_id,
		conformation::ResidueCOP new_state,
		basic::datacache::BasicDataCache & residue_data_cache
	) override;

	Real total_energy() override;

	//Real npd_hbond_weight() const { return npd_hbond_weight_; }

	utility::vector1< char > & hbonding_to_res();

	NPDHBondOP unused_hbond();
	void return_hbond_to_queue( NPDHBondOP const & hbond );

	Real npd_hb_weight(
		scoring::hbonds::HBEvalType eval_type,
		bool intra_res
	);

	//required functions to override

	Size count_static_memory() const override {return 0;}
	Size count_dynamic_memory() const override {return 0;}

	void delete_edge( utility::graph::Edge * ) override;

	NPDHBSimpleNode *
	get_npdhb_simple_node( Size ind ) {
		return static_cast< NPDHBSimpleNode * > ( get_node( ind ) );
	}

	NPDHBSimpleNode const *
	get_npdhb_simple_node( Size ind ) const {
		return static_cast< NPDHBSimpleNode const * > ( get_node( ind ) );
	}

protected:

	utility::graph::Node* create_new_node( Size node_index ) override;
	utility::graph::Edge* create_new_edge( Size index1, platform::Size index2 ) override;

	NPDHBSimpleNode * npdhb_cast( utility::graph::Node * );
	NPDHBSimpleNode const * npdhb_cast( utility::graph::Node const * ) const;

private:

	scoring::hbonds::HBondDatabaseCOP hbond_database_;
	scoring::hbonds::HBondOptionsCOP hbond_options_;
	scoring::hbonds::NPDHBondSetCOP  npd_hbond_set_;

	Real npd_hbond_weight_;

	utility::vector1< char > hbonding_to_res_;
	std::list< NPDHBondOP > hbonds_queue_;
};


} //namespace interaction_graph
} //namespace pack
} //namespace core

#endif
