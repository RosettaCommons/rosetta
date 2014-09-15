// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/SimpleInteractionGraph.hh
/// @brief
/// @author Liz Kellogg (ekellogg@u.washington.edu)


#ifndef INCLUDED_core_pack_interaction_graph_SimpleInteractionGraph_HH
#define INCLUDED_core_pack_interaction_graph_SimpleInteractionGraph_HH

// Unit headers
#include <core/pack/interaction_graph/SimpleInteractionGraph.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/graph/Graph.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack{
namespace interaction_graph{

class SimpleNode : public graph::Node
{
public:

	//  SimpleNode::SimpleNode();
	SimpleNode( graph::Graph * owner, Size resnum );
	virtual ~SimpleNode();

	virtual  void copy_from( graph::Node const * ) {}

	virtual  void print() const {}

	virtual  platform::Size count_static_memory() const {return 1;}
	virtual  platform::Size count_dynamic_memory() const {return 1;}

	Real
	one_body_energy() const;

	Real
	proposed_one_body_energy() const;

	Real
	current_one_body_energy() const;

	/// @brief Is this node considering a state substitution?
	bool
	moved() const;

	/// @brief Copy the alternate_residue_ pointer to the current_residue_ pointer;
	/// copy the alternate energies to the current energies for this node and its
	/// incident edges.
	void
	commit_change();

	/// @brief Copy the alternate energies to the current energies for this node
	/// and its incident edges.
	void
	commit_change_no_res_pointer_update();

	/// @brief Reset state on this node so that its incident edges will no longer
	/// think it is considering an alternate conformation.
	void
	reject_change();

	/// @brief Set the current residue COP, and follow by computing the energy
	/// for this residue with its neighbors and storing those computed energies
	/// on this node's edges as their "current" energies.
	void set_current( conformation::ResidueCOP res);

	/// @brief Set the alternate residue COP and follow by computing the energy
	/// for this residue with its neighbors and storing those computed energies
	/// on this node's edges as their "proposed" energies
	void set_alternate( conformation::ResidueCOP res);

	/// @brief Passive mode behavior: set the current residue pointer without updating
	/// the current one body or two body energies.
	void set_current_no_E_update( conformation::ResidueCOP res );

	/// @brief Passive mode behavior: set the current residue pointer without updating
	/// the alternate one body or proposed two body energies.
	void set_alternate_no_E_update( conformation::ResidueCOP res );

	void update_energies_after_passive_change();

	/// @brief return the pointer to the current state (might be 0)
	conformation::ResidueCOP
	get_current() const;

	/// @brief return the pointer to the alternate state (might be 0)
	conformation::ResidueCOP
	get_alternate() const;

	Vector const &
	bb_centroid() const;

	Real
	bb_radius() const;

	Vector const &
	curr_sc_centroid() const;

	Real
	curr_sc_radius() const;

	Vector const &
	alt_sc_centroid() const;

	Real
	alt_sc_radius() const;

private:

	void update_current_one_body_energy();

	void update_alternate_one_body_energy();

	void
	initialize();

	/*Vector
	calc_sc_centroid( conformation::Residue const & res ) const;

	Real
	calc_sc_radius( conformation::Residue const & res, Vector const & centroid );*/


private:
	bool moved_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// Size resnum_;
	Real current_one_body_energy_;
	Real alternate_one_body_energy_;

	conformation::ResidueCOP current_residue_;
	conformation::ResidueCOP alternate_residue_;

	Vector bb_centroid_; // should be same for curr and alt
	Real   bb_radius_; // should be same for curr and alt

	Vector curr_sc_centroid_;
	Real   curr_sc_radius_;

	Vector alt_sc_centroid_;
	Real   alt_sc_radius_;

}; //SimpleNode


class SimpleEdge : public graph::Edge {

public:

	//SimpleEdge( graph::Graph* owner );
	SimpleEdge( graph::Graph* owner, Size res1, Size res2 );
	virtual ~SimpleEdge();

	virtual  void copy_from( graph::Edge const * ) {}

	virtual  platform::Size count_static_memory() const {return 1;}
	virtual  platform::Size count_dynamic_memory() const {return 1;}

	void compute_energy( bool use_current_node1, bool use_current_node2 );

	Real get_current_energy() const;

	Real get_proposed_energy() const;

	void update_current_energy();

	void
	update_proposed_energy();

	void commit_change();

	//required functions to override

	/**
	Real curr_scsc_E_;
	Real curr_bbsc_E_;
	Real curr_scbb_E_;

	Real alt_scsc_E_;
	Real alt_bbsc_E_;
	Real alt_scbb_E_;
	**/

private:

	inline
	SimpleInteractionGraph *
	get_simple_ig_owner();

	inline
	SimpleInteractionGraph const *
	get_simple_ig_owner() const;

	Real
	get_bb_E( conformation::Residue const & r1, conformation::Residue const & r2 );

	/// These functions are to carefully gate access to the c-style, index-from-zero arrays.
	Size
	get_bb_index( conformation::Residue const & r ) const;

	bool
	bb_bb_boundaries( Size ind1, Size ind2 ) const;

	bool
	bb_bbE_calced( Size ind1, Size ind2 ) const;

	void
	set_bb_bbE_calced( Size ind1, Size ind2 );

	void
	set_bb_bbE( Size ind1, Size ind2, Real val );

	Real
	bb_bbE( Size ind1, Size ind2 ) const;


private:
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool short_range_energies_exist_; // only evaluate short-range energies for edges between nearby residues
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool long_range_energies_exist_; // only evaluate long range energies for edges that have them

	Real current_energy_;
	Real proposed_energy_;

	bool bb_bbE_calced_[ 3 ][ 3 ]; // indexed 0 for non-pro&non-gly, 1 for pro, 2 for gly
	Real bb_bbE_[ 3 ][ 3 ];        // indexed 0 for non-pro&non-gly, 1 for pro, 2 for gly

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
class SimpleInteractionGraph : public graph::Graph {

public:

	SimpleInteractionGraph();
	virtual ~SimpleInteractionGraph();

	void set_scorefunction( scoring::ScoreFunction const & sfxn );

	scoring::ScoreFunction const &
	scorefunction() const {
	  return *sfxn_;
	}

	pose::Pose const &
	pose() const {
		return *pose_;
	}

	/// @brief Initialization where the graph adds its own edges
	void initialize( pose::Pose const & pose);

	void set_pose_no_initialize( pose::Pose const & pose );

	void commit_change( Size node_id );
	void reject_change( Size node_id );

	//returns delta-energy
	Real consider_substitution( Size node_id, conformation::ResidueCOP new_state );

	Real total_energy();

	//required functions to override

  virtual  Size count_static_memory() const {return 0;}
  virtual  Size count_dynamic_memory() const {return 0;}

  virtual void delete_edge( graph::Edge * );

	SimpleNode *
	get_simple_node( Size ind ) {
		return static_cast< SimpleNode * > ( get_node( ind ) );
	}

	SimpleNode const *
	get_simple_node( Size ind ) const {
		return static_cast< SimpleNode const * > ( get_node( ind ) );
	}

protected:

  virtual graph::Node* create_new_node( Size node_index );
  virtual graph::Edge* create_new_edge( Size index1, platform::Size index2 );

  //virtual Edge* create_new_edge( Edge * example_edge );

private:

	scoring::ScoreFunctionOP sfxn_;
	pose::PoseCOP pose_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// Real accumulated_ediff_; //since data stored in energies
  //object, need to subtract out diffs

}; //SimpleInteractionGraph


inline
SimpleInteractionGraph *
SimpleEdge::get_simple_ig_owner() {
	return static_cast< SimpleInteractionGraph * > ( get_owner() );
}

inline
SimpleInteractionGraph const *
SimpleEdge::get_simple_ig_owner() const {
	return static_cast< SimpleInteractionGraph const * > ( get_owner() );
}


} //namespace interaction_graph
} //namespace pack
} //namespace core

#endif
