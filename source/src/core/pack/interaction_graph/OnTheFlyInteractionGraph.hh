// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/OnTheFlyInteractionGraph.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_OnTheFlyInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_OnTheFlyInteractionGraph_hh

// Unit headers
#include <core/pack/interaction_graph/OnTheFlyInteractionGraph.fwd.hh>

// Package headers
// AUTO-REMOVED #include <core/pack/interaction_graph/AminoAcidNeighborSparseMatrix.hh>
#include <core/pack/interaction_graph/SparseMatrixIndex.hh>
#include <core/pack/interaction_graph/FixedBBInteractionGraph.hh>

#include <core/pack/rotamer_set/RotamerSet.fwd.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#ifdef WIN32
#include <core/conformation/Residue.hh> // WIN32 INCLUDE
#endif

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

/// C++ headers
// AUTO-REMOVED #include <utility> // std::pair

// ObjexxFCL headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>


// Utility headers


namespace core {
namespace pack {
namespace interaction_graph {

enum ResiduePairEvalType {
	whole_whole = 0,
	whole_sc,
	sc_whole,
	sc_sc
};

class OnTheFlyNode : public FixedBBNode
{
public:
	typedef std::pair< Vector, Real > BoundingSphere;

public:
	OnTheFlyNode(
		InteractionGraphBase * owner,
		int node_id,
		int num_states);

	virtual ~OnTheFlyNode();

	void set_rotamers(
		rotamer_set::RotamerSetCOP rotamers
	);

	virtual void zero_one_body_energies();
	virtual void add_to_one_body_energies( ObjexxFCL::FArray1< core::PackerEnergy > & energy1b );
	virtual void update_one_body_energy( int state, core::PackerEnergy energy);
	virtual void set_one_body_energy( int state, core::PackerEnergy energy );
	virtual void add_to_one_body_energy( int state, core::PackerEnergy energy );
	virtual void zero_one_body_energy( int state );

	inline
	int
	get_num_aa_types() const
	{
		return num_aa_types_;
	}

	inline
	utility::vector1< int > &
	get_num_states_for_aa_types()
	{
		return num_states_for_aatype_;
	}

	inline
	utility::vector1< int > const &
	get_num_states_for_aa_types() const
	{
		return num_states_for_aatype_;
	}

	inline
	int
	get_num_states_for_aa_type( int aa_type )
	{
		return num_states_for_aatype_[ aa_type ];
	}

	inline
	SparseMatrixIndex const &
	get_sparse_mat_info_for_state( int state ) const
	{
		assert( state > 0 && state <= get_num_states() );
		return sparse_mat_info_for_state_[ state ];
	}

	inline
	int
	get_state_offset_for_aatype( int aatype ) const {
		return state_offset_for_aatype_[ aatype ];
	}

	inline
	core::PackerEnergy
	get_one_body_energy( int state ) const
	{
		return one_body_energies_[ state ];
	}

	bool
	distinguish_backbone_and_sidechain() const {
		return distinguish_backbone_and_sidechain_;
	}

	void
	distinguish_backbone_and_sidechain( bool setting );

	virtual unsigned int count_static_memory() const = 0;
	virtual unsigned int count_dynamic_memory() const;

	core::PackerEnergy
	compute_rotamer_pair_energy(
		int edge_making_energy_request,
		int state_this,
		int state_other
	) const;

protected:

	inline
	OnTheFlyEdge *
	get_incident_otf_edge( int edge );

	inline
	OnTheFlyEdge const *
	get_incident_otf_edge( int edge ) const;

	inline
	OnTheFlyNode *
	get_adjacent_otf_node( int index );

	inline
	OnTheFlyNode const *
	get_adjacent_otf_node( int index ) const;

	inline
	OnTheFlyInteractionGraph *
	get_on_the_fly_owner();

	inline
	OnTheFlyInteractionGraph const *
	get_on_the_fly_owner() const;

	inline
	conformation::Residue const &
	get_rotamer( int state ) const
	{
		return *rotamers_[ state ];
	}

	inline
	BoundingSphere const &
	sc_bounding_sphere( int state ) const
	{
		return sc_bounding_spheres_[ state ];
	}

	inline
	BoundingSphere const & bb_bounding_sphere() const { return bb_bounding_sphere_; }

private:
	rotamer_set::RotamerSetCOP rotamer_set_;
	utility::vector1< conformation::ResidueCOP > rotamers_;
	utility::vector1< BoundingSphere > sc_bounding_spheres_;
	BoundingSphere bb_bounding_sphere_;

	int num_aa_types_; // rename num_aa_groups_
	utility::vector1< int > num_states_for_aatype_; // rename
	utility::vector1< int > state_offset_for_aatype_;
	utility::vector1< SparseMatrixIndex > sparse_mat_info_for_state_;
	utility::vector1< core::PackerEnergy > one_body_energies_;
	bool distinguish_backbone_and_sidechain_;
};

class OnTheFlyEdge : public FixedBBEdge
{
public:
	virtual ~OnTheFlyEdge();

	OnTheFlyEdge(
		InteractionGraphBase * owner,
		int first_node_ind,
		int second_node_ind
	);

	void
	set_ProCorrection_values(
		int node_not_necessarily_proline,
		int state,
		core::PackerEnergy bb_nonprobb_E,
		core::PackerEnergy bb_probb_E,
		core::PackerEnergy sc_nonprobb_E,
		core::PackerEnergy sc_probb_E
	);

	inline
	core::PackerEnergy
	get_proline_correction_for_node(
		int node_ind,
		int state
	) const
	{
		int which_node = node_ind == get_node_index( 0 ) ? 0 : 1;
		return get_proline_correction( which_node, state );
	}


	virtual unsigned int count_static_memory() const = 0;
	virtual unsigned int count_dynamic_memory() const;

	bool long_range_interactions_exist() const { return long_range_interactions_exist_; }
	bool short_range_interactions_exist() const { return short_range_interactions_exist_; }

	void note_long_range_interactions_exist() { long_range_interactions_exist_ = true; }
	void note_short_range_interactions_exist() { short_range_interactions_exist_ = true; }

	ResiduePairEvalType
	eval_type( int node_index ) const {
		return eval_types_[ which_node( node_index ) ];
	}

protected:

	inline
	core::PackerEnergy
	get_proline_correction(
		int which_node,
		int state
	) const
	{
		return proline_corrections_[ which_node ][ state ];
	}

	inline
	OnTheFlyNode const *
	get_otf_node( int which_node ) const
	{
		return static_cast< OnTheFlyNode const * > (get_node( which_node ));
	}

	inline
	OnTheFlyNode *
	get_otf_node( int which_node )
	{
		return static_cast< OnTheFlyNode * > (get_node( which_node ));
	}

private:

	utility::vector1< core::PackerEnergy > proline_corrections_[ 2 ];
	ResiduePairEvalType eval_types_[ 2 ];
	bool long_range_interactions_exist_;
	bool short_range_interactions_exist_;
};

class OnTheFlyInteractionGraph : public FixedBBInteractionGraph
{
public:
	typedef pose::Pose               Pose;
	typedef pose::PoseOP             PoseOP;
	typedef scoring::ScoreFunction   ScoreFunction;
	typedef scoring::ScoreFunctionOP ScoreFunctionOP;

public:
	OnTheFlyInteractionGraph( int num_nodes );
	~OnTheFlyInteractionGraph();

	virtual void initialize( rotamer_set::RotamerSetsBase const & rot_sets );

	inline
	int get_num_aatypes() const
	{
		return num_aa_types_;
	}

	bool
	distinguish_backbone_and_sidechain_for_node( int node ) const;

	void
	distinguish_backbone_and_sidechain_for_node( int node, bool setting );


	void
	set_score_function( ScoreFunction const & );

	void
	set_pose( Pose const & );

	inline
	Pose const &
	pose() const
	{
		return *pose_;
	}

	/// @brief debugging only -- modify the pose during simulated annealing, if you're so inclined
	inline
	Pose &
	non_const_pose()
	{
		return *pose_;
	}


	inline
	ScoreFunction const &
	score_function() const
	{
		assert( score_function_ );
		return *score_function_;
	}

	/*
	// for using ResidueWeightMap
	inline void set_residue_weight_map(PackerTaskResidueWeightMap const & residue_weight_map_in) {
		residue_weight_map_ = residue_weight_map_in;
	}
	inline float get_residue_weights(int seqpos1, int aa1, int seqpos2, int aa2 ) const {
		{ return residue_weight_map_.get_weight(seqpos1, aa1, seqpos2, aa2); };
	}
	inline bool check_empty_weight_map() { return residue_weight_map_.check_empty_map(); };
	*/


	void
	zero_one_body_energy_for_node_state(
		int node_ind,
		int state
	);

	void
	add_to_one_body_energy_for_node_state(
		int node_ind,
		int state,
		core::PackerEnergy one_body_energy
	);

	void
	set_one_body_energy_for_node_state(
		int node_ind,
		int state,
		core::PackerEnergy one_body_energy
	);

	virtual
	core::PackerEnergy
	get_one_body_energy_for_node_state( int node, int state);

	void
	set_sparse_aa_info_for_edge(
		int node1,
		int node2,
		ObjexxFCL::FArray2_bool const & sparse_conn_info
	);


	void
	reset_rpe_calculations_count();


	Size
	get_num_rpe_calculations_count() const;

	/// @brief to be called by owned OTF node only
	void
	note_rpe_calculated() const {
		++num_rpe_calcs_;
	}

	//virtual bool build_sc_only_rotamer() const = 0;

	void
	set_ProCorrection_values_for_edge(
		int node1,
		int node2,
		int node_not_neccessarily_proline,
		int state,
		core::PackerEnergy bb_nonprobb_E,
		core::PackerEnergy bb_probb_E,
		core::PackerEnergy sc_nonprobb_E,
		core::PackerEnergy sc_probb_E
	);


	void
	note_short_range_interactions_exist_for_edge(
		int node1,
		int node2
	);

	void
	note_long_range_interactions_exist_for_edge(
		int node1,
		int node2
	);

	virtual unsigned int count_dynamic_memory() const;

protected:

	inline
	OnTheFlyNode *
	get_on_the_fly_node( int node_index )
	{
		return ( static_cast< OnTheFlyNode * > (get_node( node_index )) );
	}

	inline
	OnTheFlyNode const *
	get_on_the_fly_node( int node_index ) const
	{
		return ( static_cast< OnTheFlyNode const * > (get_node( node_index )) );
	}

private:
	int num_aa_types_;
	mutable Size num_rpe_calcs_;

	scoring::ScoreFunctionOP score_function_;
	pose::PoseOP pose_;

	// Additional per-residue (or per-residue-per-aa) weights
	// Note: currently computed OUTSIDE the IG, except in the case of OnTheFly IGs
	//PackerTaskResidueWeightMap residue_weight_map_;

};


inline
OnTheFlyEdge *
OnTheFlyNode::get_incident_otf_edge( int edge )
{
	assert( dynamic_cast< OnTheFlyEdge * >  (get_incident_edge( edge )) );
	return static_cast< OnTheFlyEdge * >  (get_incident_edge( edge ));
}

inline
OnTheFlyEdge const *
OnTheFlyNode::get_incident_otf_edge( int edge ) const
{
	assert( dynamic_cast< OnTheFlyEdge const * >  (get_incident_edge( edge )) );
	return static_cast< OnTheFlyEdge const * >  (get_incident_edge( edge ));
}

inline
OnTheFlyNode *
OnTheFlyNode::get_adjacent_otf_node( int index )
{
	assert( dynamic_cast< OnTheFlyNode * > ( get_adjacent_node( index ) ));
	return static_cast< OnTheFlyNode * > ( get_adjacent_node( index ) );
}

inline
OnTheFlyNode const *
OnTheFlyNode::get_adjacent_otf_node( int index ) const
{
	assert( dynamic_cast< OnTheFlyNode const * > ( get_adjacent_node( index ) ));
	return static_cast< OnTheFlyNode const * > ( get_adjacent_node( index ) );
}


inline
OnTheFlyInteractionGraph *
OnTheFlyNode::get_on_the_fly_owner()
{
	assert( dynamic_cast< OnTheFlyInteractionGraph * > ( get_owner() ) );
	return static_cast< OnTheFlyInteractionGraph * > ( get_owner() );
}

inline
OnTheFlyInteractionGraph const *
OnTheFlyNode::get_on_the_fly_owner() const
{
	assert( dynamic_cast< OnTheFlyInteractionGraph const * > ( get_owner() ) );
	return static_cast< OnTheFlyInteractionGraph const * > ( get_owner() );
}

}
}
}

#endif //INCLUDED_core_pack_interaction_graph_OnTheFlyInteractionGraph_HH
