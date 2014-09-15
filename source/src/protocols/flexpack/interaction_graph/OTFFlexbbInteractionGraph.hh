// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/flexpack/interaction_graph/OTFFlexbbIteractionGraph.hh
/// @brief  Declaration for on-the-fly flexible-backbone-packing interaction graph interface & base classes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_flexpack_interaction_graph_OTFFlexbbInteractionGraph_hh
#define INCLUDED_protocols_flexpack_interaction_graph_OTFFlexbbInteractionGraph_hh

//#ifndef DEBUG_OTF_FLEXBB_ENERGIES
//#define DEBUG_OTF_FLEXBB_ENERGIES // disable this line for speed
//#endif


/// Unit headers
#include <protocols/flexpack/interaction_graph/OTFFlexbbInteractionGraph.fwd.hh>

/// Package headers
#include <protocols/flexpack/interaction_graph/FlexbbInteractionGraph.hh>
// AUTO-REMOVED #include <protocols/flexpack/rotamer_set/FlexbbRotamerSet.fwd.hh>
#include <protocols/flexpack/OtherContextScoreFunction.hh> /// CREATE THE .FWD.HH FILE

/// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#ifdef WIN32
#include <core/conformation/Residue.hh>
#endif
/// ObjexxFCL headers
#include <ObjexxFCL/FArray4D.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace flexpack {
namespace interaction_graph {

class OTFFlexbbNode : public FlexbbNode
{
public:
	typedef FlexbbNode parent;
	typedef core::conformation::Residue    Residue;
	typedef core::conformation::ResidueCOP ResidueCOP;
	typedef core::Size                     Size;
	typedef core::PackerEnergy       PackerEnergy;
	typedef core::Real                     Real;
	typedef core::Vector                   Vector;
	typedef core::DistanceSquared          DistanceSquared;

public:

	OTFFlexbbNode( OTFFlexbbInteractionGraph *, int node_id, int num_states );
	virtual ~OTFFlexbbNode();
	virtual void print() const;

	virtual unsigned int count_dynamic_memory() const;

	Residue const &
	rotamer( int index ) const {
		return *rotamers_[ index ];
	}

	void set_rotamer( int state, ResidueCOP rotamer );
	void declare_all_rotamers_initialized();

	Real bounding_radius_for_rotamers( int aatype, int bb ) const;

	bool rotamer_is_proline( int state ) const { return rotamer_is_proline_[ state ]; }
	bool rotamer_is_glycine( int state ) const { return rotamer_is_glycine_[ state ]; }

protected:
	/// Downcast pointers to incident edges, adjacent nodes, and the owning graph
	inline
	OTFFlexbbEdge const * get_incident_otfflexbb_edge( int index ) const;

	inline
	OTFFlexbbEdge * get_incident_otfflexbb_edge( int index );

	inline
	OTFFlexbbNode const * get_adjacent_otfflexbb_node( int index ) const;

	inline
	OTFFlexbbNode * get_adjacent_otfflexbb_node( int index );

	inline
	OTFFlexbbInteractionGraph const * get_otfflexbbig_owner() const;

	inline
	OTFFlexbbInteractionGraph * get_otfflexbbig_owner();

private:
	utility::vector1< ResidueCOP > rotamers_;
	utility::vector1< utility::vector1< Real > > bounding_volumes_for_bb_for_aa_;
	utility::vector1< unsigned char > rotamer_is_proline_;
	utility::vector1< unsigned char > rotamer_is_glycine_;

};

class OTFFlexbbEdge : public FlexbbEdge
{
public:
	typedef FlexbbEdge parent;
	typedef core::conformation::Residue Residue;
	typedef ObjexxFCL::FArray2D< PackerEnergy > FArray2D_PackerEnergy;

public:
	OTFFlexbbEdge( OTFFlexbbInteractionGraph * owner, int node1, int node2 );
	virtual ~OTFFlexbbEdge();

	PackerEnergy
	compute_samebbconf_alternate_state_energy_first_node();

	PackerEnergy
	compute_samebbconf_alternate_state_energy_second_node();

	PackerEnergy
	compute_altbbconf_alternate_state_energy();

	void
	otfedge_note_substitution_accepted();

	virtual unsigned int count_dynamic_memory() const;

	void
	set_ProCorrection_values(
		int node_not_necessarily_proline,
		int state,
		int other_bb,
		PackerEnergy bb_nonprobb_E,
		PackerEnergy bb_probb_E,
		PackerEnergy sc_nonprobb_E,
		PackerEnergy sc_probb_E
	);

	void
	set_GlyCorrection_values(
		int node_not_necessarily_glycine,
		int state,
		int other_bb,
		PackerEnergy bb_nonglybb_E,
		PackerEnergy bb_glybb_E,
		PackerEnergy sc_nonglybb_E,
		PackerEnergy sc_glybb_E
	);

	virtual void prepare_for_simulated_annealing();

	void note_long_range_interactions_exist();

	void print_alt_energies() const;

protected:

	/// Downcasts
	inline
	OTFFlexbbNode const * get_otfflexbb_node( int index ) const;

	inline
	OTFFlexbbNode * get_otfflexbb_node( int index );

	inline
	OTFFlexbbInteractionGraph const * get_otfflexbbig_owner() const;

	inline
	OTFFlexbbInteractionGraph * get_otfflexbbig_owner();

protected:

	Residue const &
	alt_rot( int which_node ) const {
		return get_otfflexbb_node( which_node )->rotamer( nodes_alt_state( which_node ) );
	}

	void
	zero_state_on_node( int which_node );

private:

	//// @brief Is a particular state proline?  Used as a signal for proline correction terms.
	inline
	bool state_is_proline( int which_node, int state ) const;

	inline
	bool
	state_is_glycine( int which_node, int state ) const;

	int compact_bbindex( int index ) const { return nodes_part_of_same_flexseg() ? 1 : index; }

	/// If this flag is "true" then the two tables below are not allocated, and instead
	/// sc/bb and bb/bb energies are calculated alongside sc/sc energies in the OTF
	/// calculations.  Slower, but uses less memory.
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool compute_bbbb_and_scbb_otf_;

	PackerEnergy all_vs_bb_energy_curr_conf_[ 2 ];
	PackerEnergy procorr_curr_conf_[ 2 ];
	PackerEnergy glycorr_curr_conf_[ 2 ];
	PackerEnergy scsc_energy_curr_conf_;

	PackerEnergy all_vs_bb_energy_alt_conf_[ 2 ];
	PackerEnergy procorr_alt_conf_[ 2 ];
	PackerEnergy glycorr_alt_conf_[ 2 ];
	PackerEnergy scsc_energy_alt_conf_;

	// if both nodes are flexible, they store their "1 body" sc/bb and bb/bb energies here.
	// This is relatively expensive.  I'm considering an alternative that computes
	// sc/bb and bb/bb energies on the fly as well...
	// Farray indexed ( state_this, compact bb_other )
	// compact == if the edge connects two nodes that are part of the same flexible segment, then
	// the backbone index is always "1" and the dimension for anything that's indexed by the backbone
	// is 1 (thereby occupying less memory).
	FArray2D_PackerEnergy all_vs_bb_energies_[ 2 ];

	// if either node becomes proline, add in the proline correction for the state on the other node.
	/// Farray indexed ( state_this, compact bb_other )
	FArray2D_PackerEnergy procorr_energies_[ 2 ];


	// if either node becomes glycine, add in the glycine correction for the state on the other node.
	/// Farray indexed ( state_this, compact bb_other )
	FArray2D_PackerEnergy glycorr_energies_[ 2 ];

	ObjexxFCL::FArray4D< unsigned char > sr_aa_neighbors_; // indexed ( aa1, aa2, non-compact bb1, compact bb2 )
	bool lr_energies_exist_;

	core::pose::PoseCOP pose_;
	core::scoring::ScoreFunctionCOP sfxn_;
};

class OTFFlexbbInteractionGraph : public FlexbbInteractionGraph
{
public:
	typedef FlexbbInteractionGraph parent;
	typedef core::Real             Real;
	typedef core::pose::PoseOP     PoseOP;
	typedef core::pose::PoseCOP    PoseCOP;
	typedef core::pose::Pose       Pose;
	typedef core::scoring::ScoreFunctionOP  ScoreFunctionOP;
	typedef core::scoring::ScoreFunctionCOP ScoreFunctionCOP;
	typedef core::scoring::ScoreFunction    ScoreFunction;


public:

	OTFFlexbbInteractionGraph( int num_nodes );
	virtual ~OTFFlexbbInteractionGraph();

	virtual void initialize( core::pack::rotamer_set::RotamerSetsBase const & );

	void
	set_ProCorrection_values_for_edge(
		int node1,
		int node2,
		int node_not_necessarily_proline,
		int state,
		int other_bb,
		PackerEnergy bb_nonprobb_E,
		PackerEnergy bb_probb_E,
		PackerEnergy sc_nonprobb_E,
		PackerEnergy sc_probb_E
	);

	void
	set_GlyCorrection_values_for_edge(
		int node1,
		int node2,
		int node_not_necessarily_glycine,
		int state,
		int other_bb,
		PackerEnergy bb_nonglybb_E,
		PackerEnergy bb_glybb_E,
		PackerEnergy sc_nonglybb_E,
		PackerEnergy sc_glybb_E
	);

	virtual unsigned int count_dynamic_memory() const;

	/// @brief Pose must be set before any edges are added to the graph.
	virtual void set_pose( Pose const & pose );
	/// @brief Score function must be set before any edges are added to the graph.
	virtual void set_scorefxn( ScoreFunction const & sfxn );

	/// @brief Edges request the pose and the score function at the time of their creation.
	PoseCOP          get_pose() const;
	ScoreFunctionCOP get_scorefxn() const;

	/// @brief Informs the edge connecting nodes 1 and 2 that they require long range interactions.
	/// Note -- the edge must already exist.
	void note_long_range_interactions_exist_for_edge( int node1, int node2 );

	void debug_note_considered_substitution( core::conformation::Residue const & alt_rotamer, int index );
	void debug_note_projected_deltaE_of_considered_substitution( PackerEnergy deltaE, PackerEnergy node_alt_total, bool require_match = true );
	void debug_note_accepted_substitution();
	void debug_note_rejected_substitution();

protected:
	/// Downcasts
	OTFFlexbbNode const * get_otfflexbb_node( int index ) const
	{ return static_cast< OTFFlexbbNode const * > (get_node( index )); }

	OTFFlexbbNode * get_otfflexbb_node( int index )
	{ return static_cast< OTFFlexbbNode * > (get_node( index )); }

	OTFFlexbbEdge const * find_otfflexbb_edge( int node1, int node2 ) const
	{
		core::pack::interaction_graph::EdgeBase const * edge = find_edge( node1, node2 );
		if ( edge ) return static_cast< OTFFlexbbEdge const * > ( edge );
		else return 0;
	}

	OTFFlexbbEdge * find_otfflexbb_edge( int node1, int node2 )
	{
		core::pack::interaction_graph::EdgeBase * edge = find_edge( node1, node2 );
		if ( edge ) return static_cast< OTFFlexbbEdge * > ( edge );
		else return 0;
	}

	OTFFlexbbEdge const * cast_otfflexbb_edge( EdgeBase const * edge ) const
	{ assert( mine( edge ) ); return static_cast< OTFFlexbbEdge const * > ( edge ); }

	OTFFlexbbEdge * cast_otfflexbb_edge( EdgeBase * edge )
	{ assert( mine( edge ) ); return static_cast< OTFFlexbbEdge * > ( edge ); }


private:

	PoseOP pose_;
	ScoreFunctionOP sfxn_;
	OtherContextScoreFunctionOP oc_sfxn_;

	/// For debugging purposes
	PoseOP current_pose_;
	Real current_pose_energy_;
	PoseOP alternate_pose_;
	Real alternate_pose_energy_;

	utility::vector1< Size > changing_seqpos_;
	utility::vector1< core::conformation::ResidueOP > alt_rots_;
	utility::vector1< int > alt_rot_inds_;

	utility::vector1< Size > resid_2_moltenres_;
	utility::vector1< Size > moltenres_2_resid_;
};

inline
OTFFlexbbEdge const * OTFFlexbbNode::get_incident_otfflexbb_edge( int index ) const
{	return static_cast< OTFFlexbbEdge const * > ( get_incident_edge( index )); }

inline
OTFFlexbbEdge * OTFFlexbbNode::get_incident_otfflexbb_edge( int index )
{	return static_cast< OTFFlexbbEdge * > ( get_incident_edge( index )); }

inline
OTFFlexbbNode const * OTFFlexbbNode::get_adjacent_otfflexbb_node( int index ) const
{	return static_cast< OTFFlexbbNode const * > ( get_adjacent_node( index )); }

inline
OTFFlexbbNode * OTFFlexbbNode::get_adjacent_otfflexbb_node( int index )
{	return static_cast< OTFFlexbbNode * > ( get_adjacent_node( index )); }

inline
OTFFlexbbInteractionGraph const * OTFFlexbbNode::get_otfflexbbig_owner() const
{	return static_cast< OTFFlexbbInteractionGraph const * > (get_owner()); }

inline
OTFFlexbbInteractionGraph * OTFFlexbbNode::get_otfflexbbig_owner()
{	return static_cast< OTFFlexbbInteractionGraph * > (get_owner()); }


/// Downcasts
inline
OTFFlexbbNode const * OTFFlexbbEdge::get_otfflexbb_node( int index ) const
{	return static_cast< OTFFlexbbNode const * > (get_node( index )); }

inline
OTFFlexbbNode * OTFFlexbbEdge::get_otfflexbb_node( int index )
{	return static_cast< OTFFlexbbNode * > (get_node( index )); }

inline
OTFFlexbbInteractionGraph const * OTFFlexbbEdge::get_otfflexbbig_owner() const
{	return static_cast< OTFFlexbbInteractionGraph const * > (get_owner()); }

inline
OTFFlexbbInteractionGraph * OTFFlexbbEdge::get_otfflexbbig_owner()
{	return static_cast< OTFFlexbbInteractionGraph * > (get_owner()); }


}
}
}

#endif
