// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/NPDHBondInteractionGraph.hh
/// @brief  Interaction graph which implements a non-PD score that optimizes against surface hydrophobic patches.
/// @reference Computational Protein Design with Explicit Consideration of Surface Hydrophobic Patches. R. Jacak, A. Leaver-Fay, and B. Kuhlman. Proteins. 2012 Mar;80(3):825-38.
/// @author Ron Jacak (ron.jacak@gmail.com)
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_interaction_graph_NPDHBondInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_NPDHBondInteractionGraph_hh

// Unit headers
#include <core/pack/interaction_graph/NPDHBondInteractionGraph.fwd.hh>

// Package headers
#include <core/pack/interaction_graph/AdditionalBackgroundNodesInteractionGraph.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>

// Project Headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/conformation/Residue.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/NPDHBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/constants.hh>

// Basic headers
#include <basic/Tracer.hh>

//Utility Headers
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>  // needed to get arg_max - DO NOT AUTO-REMOVE!
#include <utility/exit.hh>
#include <utility/string_util.hh> // needed to get trim - DO NOT AUTOREMOVE!
#include <utility/graph/Graph.hh>

//ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1.io.hh> // needed to stream operator of FArray1Dint, line 4738 - DO NOT AUTOREMOVE!
#include <ObjexxFCL/format.hh> // needed for I() - DO NOT AUTOREMOVE!

//C++ Headers
#include <vector>


namespace core {
namespace pack {
namespace interaction_graph {

static basic::Tracer TR( "core.pack.npd_hbond_ig" );

template < typename V, typename E, typename G > class NPDHBondNode;
template < typename V, typename E, typename G > class NPDHBondBackgroundNode;
template < typename V, typename E, typename G > class NPDHBondEdge;
template < typename V, typename E, typename G > class NPDHBondBackgroundEdge;
template < typename V, typename E, typename G > class NPDHBondInteractionGraph;

struct NPDHBond
{
	Real energy_;
	Real sfxn_wt_;
	Real don_wt_;
	Real acc_wt_;
	Real don_wt_alt_;
	Real acc_wt_alt_;
	Size don_rsd_;
	Size acc_rsd_;
	Size don_atm_;
	Size acc_atm_;
};

typedef utility::pointer::shared_ptr< NPDHBond > NPDHBondOP;

/// @brief Compute the don_wt_alt_s and the acc_wt_alt_s from the perspective of the
/// input residue given the (complete) arrays of NPDHBondOPs for each atom.
/// @remarks Used by both the NPDHBondNode and the NPDHBondBackgroundNode, and so defined
/// in the .cc file.
void
compute_alt_weights_for_npd_hbonds(
	conformation::Residue const & res,
	utility::vector1< utility::vector1< NPDHBondOP > > const & atom_hbonds, // the array is constant, the pointed-at-hbonds are going to be modified
	utility::vector1< Real > & tmp_energies, // the sfxn-weights-weighted energies for the hbonds
	utility::vector1< Real > & temp_weights
);


//----------------------------------------------------------------------------//
//---------------------------- NPDHBond Node Class ---------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// Defines a FirstClass node which will keep track of the hydrogen bonds coming in to the
/// residue and the hydrogen bond score.
/// FirstClassNode is defined and implemented in AdditionalBackgroundNodesInteractionGraph.
///
/// @remarks
/// No public default constructor makes this class uncopyable.
///
template < typename V, typename E, typename G >
class NPDHBondNode : public FirstClassNode< V, E, G > {

public:
	typedef FirstClassNode< V, E, G > parent;
	typedef V grandparent;

public:
	NPDHBondNode( G* owner, int node_index, int num_states );
	virtual ~NPDHBondNode();

	// setter for the rotamers object. called at the very beginning of the NPDHBIG::initialize() method
	void set_rotamers( rotamer_set::RotamerSetCOP rotamers );

	conformation::Residue const & get_rotamer( int state ) const;
	conformation::ResidueCOP get_rotamer_op( int state ) const;
	conformation::Residue const & curr_state_rotamer() const;
	conformation::ResidueOP curr_state_rotamer_op() const;
	conformation::Residue const & alt_state_rotamer() const;
	conformation::ResidueCOP alt_state_rotamer_op() const;

	Size seqpos() const { return seqpos_; }

	virtual void prepare_for_simulated_annealing();

	// TO DO: possible optimization -- precompute hbonds to background
	// TO DO: possible optimization -- precompute bb/bb hbonds if you can prove they won't change

	// for NPDHBIG entry point blanket_assign_state_0()
	virtual void assign_zero_state();
	void acknowledge_neighbors_substitution();
	void acknowledge_neighbors_state_zeroed( Size seqpos );

	core::PackerEnergy calculate_PD_deltaE_for_substitution( int alternate_state, core::PackerEnergy & prev_PDenergies_for_node );
	core::PackerEnergy get_pd_energy_delta();

	/// @brief Return the change in energy induced by substituting an alternate rotamer
	/// at this position.
	Real consider_alternate_state();

	void prepare_for_neighbors_substitution( Size nbrs_seqpos );

	void compute_alt_weights_for_hbonds(
		bool curr_state
	);

	utility::vector1< NPDHBondOP > & current_hbs();
	utility::vector1< NPDHBondOP > & alternate_hbs();
	utility::vector1< utility::vector1< NPDHBondOP > > & alternate_hbs_for_atoms();
	utility::vector1< NPDHBondOP > & alternate_hbs_for_atom( Size atom_index );

	Real get_upper_npd_hbond_energy_totals() const;

	void commit_considered_substitution();
	void reset_after_rejected_substitution();

	virtual unsigned int getMemoryUsageInBytes() const;
	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	virtual void print() const;

protected:
	inline
	NPDHBondEdge< V, E, G > const *
	get_incident_npd_hbond_edge( int index ) const {
		return static_cast< NPDHBondEdge< V, E, G > const * > (parent::get_incident_edge( index ));
	}

	inline
	NPDHBondEdge< V, E, G > *
	get_incident_npd_hbond_edge( int index ) {
		return static_cast< NPDHBondEdge< V, E, G > * > (parent::get_incident_edge( index ));
	}

	inline
	NPDHBondNode< V, E, G > const *
	get_adjacent_npd_hbond_node( int index ) const {
		return static_cast< NPDHBondNode< V, E, G > const * > ( parent::get_adjacent_node( index ));
	}

	inline
	NPDHBondNode< V, E, G > *
	get_adjacent_npd_hbond_node( int index ) {
		return static_cast< NPDHBondNode< V, E, G > * > ( parent::get_adjacent_node( index ));
	}

	inline
	NPDHBondBackgroundEdge< V, E, G > const *
	get_edge_to_npd_hbond_bg_node( int index ) const {
		return static_cast< NPDHBondBackgroundEdge< V, E, G > const * > ( parent::get_edge_to_bg_node( index ));
	}

	inline
	NPDHBondBackgroundEdge< V, E, G > *
	get_edge_to_npd_hbond_bg_node( int index ) {
		return static_cast< NPDHBondBackgroundEdge< V, E, G > * > ( parent::get_edge_to_bg_node( index ));
	}

	inline
	NPDHBondBackgroundNode< V, E, G > const *
	get_adjacent_npd_hbond_bg_node( int index ) const {
		return static_cast< NPDHBondBackgroundNode< V, E, G > const * > ( parent::get_adjacent_background_node( index ));
	}

	inline
	NPDHBondBackgroundNode< V, E, G > *
	get_adjacent_npd_hbond_bg_node( int index ) {
		return static_cast< NPDHBondBackgroundNode< V, E, G > * > ( parent::get_adjacent_background_node( index ));
	}

	inline
	NPDHBondInteractionGraph< V, E, G > const * get_npd_hbond_owner() const {
		return static_cast< NPDHBondInteractionGraph< V, E, G > const *> (parent::get_owner());
	}

	inline
	NPDHBondInteractionGraph< V, E, G > * get_npd_hbond_owner() {
		return static_cast< NPDHBondInteractionGraph< V, E, G > *> ( parent::get_owner() );
	}

private:

	// no default constructor, uncopyable
	NPDHBondNode();
	NPDHBondNode( NPDHBondNode< V, E, G > const & );
	NPDHBondNode< V, E, G > & operator = ( NPDHBondNode< V, E, G > const & );

private:
	core::Size seqpos_;

	utility::vector1< conformation::ResidueCOP > rotamers_vector_;
	utility::vector1< Size > restype_group_for_rotamers_;
	Size max_natoms_;

	utility::vector1< utility::vector1< NPDHBondOP > > curr_atom_hbonds_;
	utility::vector1< utility::vector1< NPDHBondOP > > alt_atom_hbonds_;

	utility::vector1< NPDHBondOP > curr_hbonds_;
	utility::vector1< NPDHBondOP > alt_hbonds_;

	utility::vector1< Real > tmp_energies_;
	utility::vector1< Real > tmp_weights_;

};


//----------------------------------------------------------------------------//
//------------------------- NPDHBond Background Node Class -----------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// Defines a Background Node which will contribute to changes in SASA/hpatchE due to state changes on neighboring nodes,
/// and not because of state changes to it.
/// No default constructor makes this class uncopyable
///
template < typename V, typename E, typename G >
class NPDHBondBackgroundNode : public BackgroundNode< V, E, G > {

public:
	typedef BackgroundNode< V, E, G > parent;

public:
	NPDHBondBackgroundNode( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int node_index );
	virtual ~NPDHBondBackgroundNode();

	void set_rotamer( conformation::ResidueOP const & rotamer );
	conformation::Residue const & get_rotamer() const;
	conformation::ResidueCOP get_rotamer_op() const;
	Size seqpos() const;

	virtual void prepare_for_simulated_annealing();

	void prepare_for_neighbors_substitution( Size nbrs_seqpos );
	void compute_alt_weights_for_hbonds();

	utility::vector1< NPDHBondOP > & current_hbs();
	utility::vector1< NPDHBondOP > & alternate_hbs();
	utility::vector1< utility::vector1< NPDHBondOP > > & alternate_hbs_for_atoms();
	utility::vector1< NPDHBondOP > & alternate_hbs_for_atom( Size atom_index );

	Real get_upper_npd_hbond_energy_totals() const;

	void acknowledge_substitution();
	void acknowledge_neighbors_state_zeroed( Size neighbors_seqpos );


	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	//void write_dot_kinemage( std::ofstream & output_kin );
	virtual void print() const;

protected:
	inline
	NPDHBondBackgroundEdge< V, E, G > * get_npd_hbond_bg_edge( int index ) {
		return (NPDHBondBackgroundEdge< V, E, G > *) parent::get_incident_edge( index );
	}

	inline
	NPDHBondInteractionGraph< V, E, G >* get_npd_hbond_owner() const {
		return (NPDHBondInteractionGraph< V, E, G > *) parent::get_owner();
	}

private:
	conformation::ResidueOP rotamer_;
	Size seqpos_;

	utility::vector1< utility::vector1< NPDHBondOP > > curr_atom_hbonds_;
	utility::vector1< utility::vector1< NPDHBondOP > > alt_atom_hbonds_;

	utility::vector1< NPDHBondOP > curr_hbonds_;
	utility::vector1< NPDHBondOP > alt_hbonds_;

	utility::vector1< Real > tmp_energies_;
	utility::vector1< Real > tmp_weights_;

	bool prepared_for_simA_;

	NPDHBondBackgroundNode();
	NPDHBondBackgroundNode( NPDHBondBackgroundNode< V, E, G > const & );
	NPDHBondBackgroundNode< V, E, G > & operator= ( NPDHBondBackgroundNode< V, E, G > const & );

};


//----------------------------------------------------------------------------//
//------------------------------ NPDHBond Edge Class -----------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// Defines a NPDHBond Edge which connects two first-class NPDHBond Nodes. Edges have to keep some state so that updates
/// to SASA and the hpatch score can be done fast.
///
template < typename V, typename E, typename G >
class  NPDHBondEdge : public FirstClassEdge< V, E, G > {

public:
	typedef  FirstClassEdge< V, E, G >  parent;

public:
	NPDHBondEdge( G * owner, int node1, int node2 );
	virtual ~NPDHBondEdge();

	virtual void prepare_for_simulated_annealing();

	void acknowledge_state_zeroed( int node_index, Size node_seqpos );

	void acknowledge_substitution( bool update_hbonds );

	// Virtual methods from EdgeBase
	virtual void declare_energies_final();

	virtual unsigned int getMemoryUsageInBytes() const;
	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	//Real get_current_two_body_energy() const;

	void consider_alternate_state_step1(
		int node_index,
		int state_index,
		conformation::Residue const & alt_state,
		utility::vector1< NPDHBondOP > & res_hbonds,
		utility::vector1< utility::vector1< NPDHBondOP > > & atom_hbonds,
		utility::vector1< char > & hbonding_to_res
	);

	Real consider_alternate_state_step2(
		utility::vector1< char > const & hbonding_to_res
	);

protected:

	inline
	NPDHBondNode< V, E, G > const * get_npd_hbond_node( int index ) const {
		return static_cast< NPDHBondNode< V, E, G > const * > ( E::get_node( index ));
	}

	inline
	NPDHBondNode< V, E, G > * get_npd_hbond_node( int index ) {
		return static_cast< NPDHBondNode< V, E, G > * > (E::get_node( index ));
	}

	inline
	NPDHBondInteractionGraph< V, E, G > const * get_npd_hbond_owner() const {
		return static_cast< NPDHBondInteractionGraph< V, E, G > const * > ( E::get_owner() );
	}

	inline
	NPDHBondInteractionGraph< V, E, G > * get_npd_hbond_owner() {
		return static_cast< NPDHBondInteractionGraph< V, E, G > * > ( E::get_owner() );
	}

private:
	void inform_non_changing_node_of_neighbors_change();

	//no default constructor, uncopyable
	NPDHBondEdge();
	NPDHBondEdge( NPDHBondEdge< V, E, G > const & );
	NPDHBondEdge< V, E, G > & operator = ( NPDHBondEdge< V, E, G > const & );

private:
	int node_changing_;
	int node_not_changing_;

	int  nodes_curr_states_[2];
	int  nodes_alt_states_[2];

};


//----------------------------------------------------------------------------//
//------------------- NPDHBond Background Edge Class -----------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// Defines an edge between a FirstClass (NPDHBondNode) and a background node (NPDHBondBackgroundNode)
///
/// @details
/// In addition to implementing the virtual base class methods, this class additionally defines methods
/// relating to keeping track of data relating to SASA/hpatch.
///
template < typename V, typename E, typename G >
class NPDHBondBackgroundEdge : public BackgroundToFirstClassEdge< V, E, G > {

public:
	typedef  BackgroundToFirstClassEdge< V, E, G >  parent;

public:
	NPDHBondBackgroundEdge( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int first_class_node_index, int background_node_index );
	virtual ~NPDHBondBackgroundEdge();

	void prepare_for_simulated_annealing();

	//void acknowledge_state_change( int new_state );

	void acknowledge_state_zeroed( Size node_seqpos );
	void acknowledge_substitution( bool update_hbonds );

	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	void consider_alternate_state_step1(
		int state_index,
		conformation::Residue const & alt_state,
		utility::vector1< NPDHBondOP > & res_hbonds,
		utility::vector1< utility::vector1< NPDHBondOP > > & atom_hbonds,
		utility::vector1< char > & hbonding_to_res
	);

	Real consider_alternate_state_step2(
		utility::vector1< char > const & hbonding_to_res
	);

protected:
	inline
	NPDHBondNode< V, E, G > const * get_npd_hbond_node() const {
		return static_cast< NPDHBondNode< V, E, G > const * > ( parent::get_first_class_node() );
	}

	inline
	NPDHBondNode< V, E, G > * get_npd_hbond_node() {
		return static_cast< NPDHBondNode< V, E, G > * > ( parent::get_first_class_node() );
	}

	inline
	NPDHBondBackgroundNode< V, E, G > const * get_npd_hbond_bg_node() const {
		return static_cast< NPDHBondBackgroundNode< V, E, G > const * > ( parent::get_background_node() );
	}

	inline
	NPDHBondBackgroundNode< V, E, G > * get_npd_hbond_bg_node() {
		return static_cast< NPDHBondBackgroundNode< V, E, G > * > ( parent::get_background_node() );
	}

	inline
	NPDHBondInteractionGraph< V, E, G > const * get_npd_hbond_owner() const {
		return static_cast< NPDHBondInteractionGraph< V, E, G > const * > ( parent::get_owner() );
	}

	inline
	NPDHBondInteractionGraph< V, E, G > * get_npd_hbond_owner() {
		return static_cast< NPDHBondInteractionGraph< V, E, G > * > ( parent::get_owner() );
	}

private:
	//no default constructor, uncopyable
	NPDHBondBackgroundEdge();
	NPDHBondBackgroundEdge( NPDHBondBackgroundEdge< V, E, G > const & );
	NPDHBondBackgroundEdge< V, E, G > & operator= ( NPDHBondBackgroundEdge< V, E, G > const & );

private:
	bool prepared_for_simA_;
	Size bg_res_num_atoms_;

	int nodes_curr_state_;
	int nodes_alt_state_;

};


//----------------------------------------------------------------------------//
//--------------------- NPDHBond Interaction Graph -------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// Defines the interaction graph that will keep track of changes to the hpatch score.
///
/// @details
/// In addition to implementing the virtual base class methods, this class additionally defines methods
/// relating to keeping track of data relating to hpatch.
///
template < typename V, typename E, typename G >
class NPDHBondInteractionGraph : public AdditionalBackgroundNodesInteractionGraph< V, E, G > {

public:
	typedef  AdditionalBackgroundNodesInteractionGraph< V, E, G >  parent;

public:
	NPDHBondInteractionGraph( int num_nodes );
	virtual ~NPDHBondInteractionGraph();

	pose::Pose const & pose() const { return *pose_; }
	void set_pose( pose::Pose const & pose );

	scoring::ScoreFunction const & scorefxn() const { return *scorefxn_; }
	void set_score_function( scoring::ScoreFunction const & scorefxn );

	utility::graph::Graph const & packer_neighbor_graph() const { return *neighbor_graph_; }
	void set_packer_neighbor_graph( utility::graph::Graph const & neighbor_graph );

	task::PackerTask const & packer_task() const { return *packer_task_; }
	void set_packer_task( task::PackerTask const & task );

	rotamer_set::RotamerSets const & rotamer_sets() const { return *rotamer_sets_; }
	void set_rotamer_sets( rotamer_set::RotamerSets const & rotsets );

	// void set_score_weight( Real weight ) { npd_hbond_score_weight_ = weight; }

	void initialize( rotamer_set::RotamerSetsBase const & rot_sets );

	void set_num_residues_in_protein( Size num_res );
	void set_num_background_residues( Size num_background_residues );
	void set_residue_as_background_residue( int residue );

	virtual void prepare_for_simulated_annealing();

	virtual void blanket_assign_state_0();
	Real get_npd_hbond_score();
	virtual void set_errorfull_deltaE_threshold( core::PackerEnergy deltaE );

	static void print_npd_hbond_avoidance_stats();
	static void reset_npd_hbond_avoidance_stats();

	virtual void consider_substitution(
		int node_ind,
		int new_state,
		core::PackerEnergy & delta_energy,
		core::PackerEnergy & prev_energy_for_node
	);

	core::PackerEnergy calculate_npd_hbond_deltaE();

	Real calculate_alt_state_npd_hbond_score();
	//void blanket_reset_alt_state_dots();

	virtual core::PackerEnergy commit_considered_substitution();

	virtual core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int & node_states );

	virtual core::PackerEnergy get_energy_current_state_assignment();
	virtual int get_edge_memory_usage() const;
	virtual unsigned int count_static_memory() const;
	virtual unsigned int count_dynamic_memory() const;

	using parent::get_energy_sum_for_vertex_group;

	//virtual void print_current_state_assignment() const;
	virtual Real get_energy_sum_for_vertex_group( Size group_id );

	void print_internal_energies_for_current_state_assignment();
	//void write_dot_kinemage( std::ofstream & output_kin );

	void print() const;

	// method used only by unit tests
	int bg_node_2_resid( Size node_index );
	//Size resid_2_bg_node( Size resid );
	std::vector<int> get_network_state() const;
	void set_observed_sufficient_boolean_true();

	scoring::hbonds::HBondDatabase const & hbond_database() const;
	scoring::hbonds::HBondOptions const & hbond_options() const;
	scoring::hbonds::NPDHBondSet const & npd_hbond_set() const;

	utility::vector1< char > & hbonding_to_res_vector();

	NPDHBondOP unused_hbond();
	void return_hbond_to_queue( NPDHBondOP const & hbond );

	Real npd_hb_weight( scoring::hbonds::HBEvalType eval_type, bool intra_res );

public:

	// need to make these method public for the unit tests
	inline
	NPDHBondNode< V, E, G > const * get_npd_hbond_node( int index ) const {
		return static_cast< NPDHBondNode< V, E, G > const * > (G::get_node( index ));
	}

	inline
	NPDHBondNode< V, E, G > * get_npd_hbond_node( int index ) {
		return static_cast< NPDHBondNode< V, E, G > * > (G::get_node( index ));
	}

	inline
	NPDHBondBackgroundNode< V, E, G > const * get_npd_hbond_bg_node( int index ) const {
		return static_cast< NPDHBondBackgroundNode< V, E, G > const * > (parent::get_background_node( index ));
	}

	inline
	NPDHBondBackgroundNode< V, E, G > * get_npd_hbond_bg_node( int index ) {
		return static_cast< NPDHBondBackgroundNode< V, E, G > * > (parent::get_background_node( index ));
	}

protected:
	virtual NodeBase* create_new_node( int node_index, int num_states);
	virtual EdgeBase* create_new_edge( int index1, int index2);
	virtual BackgroundNode< V, E, G >* create_background_node( int node_index );
	virtual BackgroundToFirstClassEdge< V, E, G >* create_background_edge( int fc_node_index, int bg_node_index);

	void track_npd_hbond_E_min();

	void update_internal_energy_totals_npd_hbond();

private:
	void reset_from_previous_delta_npd_hbond_comp();
	bool decide_procrastinate_npd_hbond_computations( Real const pd_deltaE, Real const threshold ) const;

	NPDHBondInteractionGraph() = delete;
	NPDHBondInteractionGraph( NPDHBondInteractionGraph< V, E, G > const & ) = delete;
	NPDHBondInteractionGraph< V, E, G > & operator = ( NPDHBondInteractionGraph< V, E, G > const & ) = delete;

private:
	pose::PoseOP pose_;
	scoring::ScoreFunctionOP scorefxn_;
	utility::graph::GraphOP neighbor_graph_;
	task::PackerTaskOP packer_task_;
	rotamer_set::RotamerSetsOP rotamer_sets_;

	scoring::EnergyMap npd_hbond_score_weights_;

	Size num_total_residues_ = 0 ;
	Size num_residues_assigned_as_background_ = 0;
	utility::vector1< Size > resid_2_bgenumeration_;
	utility::vector1< Size > bgenumeration_2_resid_;

	bool prepared_for_simulated_annealing_ = false;

	bool observed_sufficient_npd_hbond_E_to_predict_min_ = false;

	static Size num_state_substitutions_considered_;
	static Size num_npd_hbond_comps_procrastinated_;
	static Size num_npd_hbond_comps_later_made_;

	Real npd_hbond_score_min_last_100_ = 0;
	Real npd_hbond_score_min_recent_ = 0;
	Size num_substitutions_since_npd_hbond_min_update_ = 0;

	bool calculated_npd_hbond_deltaE_ = false;
	core::PackerEnergy deltaE_for_substitution_ = 0.0;

	bool last_considered_substitution_accepted_ = false;
	Size node_considering_alt_state_ = 0;
	Size alt_state_being_considered_ = 0;
	Real total_energy_current_state_assignment_ = 0;
	Real total_energy_alternate_state_assignment_ = 0;

	Real npd_hbond_energy_current_state_assignment_ = 0;
	Real npd_hbond_energy_alternate_state_assignment_ = 0;

	int num_commits_since_last_update_ = 0;
	float deltaE_threshold_for_avoiding_npd_hbond_calcs_ = 0.0;

	scoring::hbonds::HBondDatabaseCOP hbond_database_;
	scoring::hbonds::HBondOptionsCOP hbond_options_;
	scoring::hbonds::NPDHBondSetCOP  npd_hbond_set_;

	utility::vector1< char > hbonding_to_res_;

	std::list< NPDHBondOP > hbonds_queue_;

	static const int COMMIT_LIMIT_BETWEEN_UPDATES = 4096;

};


//----------------------------------------------------------------------------//
//-------------------------------- NPDHBond Node Class ---------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// NPDHBondNode constructor
///
template < typename V, typename E, typename G >
NPDHBondNode< V, E, G >::NPDHBondNode( G* owner, int node_index, int num_states ) :
	FirstClassNode< V, E, G > ( owner, node_index, num_states ),
	seqpos_( 0 ),
	max_natoms_( 0 )
{
	rotamers_vector_.resize( num_states );
	restype_group_for_rotamers_.resize( num_states, 0 );
}

template < typename V, typename E, typename G >
NPDHBondNode< V, E, G >::~NPDHBondNode() {}

///
/// @details
/// Need to save a reference to the rotamer_set so that we can determine what a particular state change will do to the score
///
template < typename V, typename E, typename G >
void NPDHBondNode< V, E, G >::set_rotamers( rotamer_set::RotamerSetCOP rotamers ) {

	// get_num_states should call the parent graph's method?
	if ( rotamers->num_rotamers() != (Size) parent::get_num_states() ) {
		utility_exit_with_message( "Number of rotamers is not equal to parents number of states. Quitting.");
	}

	debug_assert( rotamers->num_rotamers() > 0 );
	seqpos_ = rotamers->rotamer( 1 )->seqpos();

	for ( Size ii = 1; ii <= rotamers->num_rotamers(); ++ii ) {
		rotamers_vector_[ ii ] = rotamers->rotamer( ii );
		restype_group_for_rotamers_[ ii ] = rotamers->get_residue_group_index_for_rotamer( ii );
		if ( max_natoms_ < rotamers_vector_[ ii ]->natoms() ) {
			max_natoms_ = rotamers_vector_[ ii ]->natoms();
		}
	}
	curr_atom_hbonds_.resize( max_natoms_ );
	alt_atom_hbonds_.resize( max_natoms_ );
	for ( Size ii = 1; ii <= max_natoms_; ++ii ) {
		curr_atom_hbonds_[ ii ].reserve( 6 );
		alt_atom_hbonds_[ ii ].reserve( 6 );
	}
	curr_hbonds_.reserve( 10 );
	alt_hbonds_.reserve( 10 );

}

template < typename V, typename E, typename G >
conformation::Residue const &
NPDHBondNode< V, E, G >::get_rotamer( int state ) const {
	return *rotamers_vector_[ state ];
}

template < typename V, typename E, typename G >
conformation::ResidueCOP
NPDHBondNode< V, E, G >::get_rotamer_op( int state ) const {
	return rotamers_vector_[ state ];
}

template < typename V, typename E, typename G >
conformation::Residue const &
NPDHBondNode< V, E, G >::curr_state_rotamer() const {
	return *rotamers_vector_[ parent::get_current_state() ];
}

template < typename V, typename E, typename G >
conformation::ResidueOP
NPDHBondNode< V, E, G >::curr_state_rotamer_op() const {
	return rotamers_vector_[ parent::get_current_state() ];
}

template < typename V, typename E, typename G >
conformation::Residue const &
NPDHBondNode< V, E, G >::alt_state_rotamer() const {
	return * rotamers_vector_[ parent::get_alternate_state() ];
}

template < typename V, typename E, typename G >
conformation::ResidueCOP
NPDHBondNode< V, E, G >::alt_state_rotamer_op() const {
	return rotamers_vector_[ parent::get_alternate_state() ];
}

/// @brief invokes V's prep_for_simA method
template < typename V, typename E, typename G >
void
NPDHBondNode< V, E, G >::prepare_for_simulated_annealing() {

	V::prepare_for_simulated_annealing();

	if ( ! parent::get_bg_edge_vector_up_to_date_() ) {
		parent::update_bg_edge_vector();
	}
}


template < typename V, typename E, typename G >
void NPDHBondNode< V, E, G >::assign_zero_state() {

	parent::assign_zero_state();

	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		get_incident_npd_hbond_edge( ii )->acknowledge_state_zeroed( parent::get_node_index(), seqpos_ );
	}

	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		get_edge_to_npd_hbond_bg_node( ii )->acknowledge_state_zeroed( seqpos_ );
	}

}

/// @brief
/// bookkeeping to follow a neighbors state substitution.  this method gets called when a NPDHBondNode commits a sub
/// and then broadcasts that change to all its neighboring fc nodes via the incident NPDHBondEdges.
template < typename V, typename E, typename G >
void NPDHBondNode< V, E, G >::acknowledge_neighbors_substitution()
{
	curr_hbonds_.swap( alt_hbonds_ );
	curr_atom_hbonds_.swap( alt_atom_hbonds_ );
	for ( auto const & hb : curr_hbonds_ ) {
		hb->don_wt_ = hb->don_wt_alt_;
		hb->acc_wt_ = hb->acc_wt_alt_;
	}
}

template < typename V, typename E, typename G >
void NPDHBondNode< V, E, G >::acknowledge_neighbors_state_zeroed( Size seqpos )
{
	// 1. Purge the hbonds to the neighbor from the alt_* set
	// 2. Compute a new set of weights for the remaining hydrogen bonds
	// 3. Then accept the alternate as the current.

	prepare_for_neighbors_substitution( seqpos );
	compute_alt_weights_for_hbonds( true ); // true -> use current rotamer
	acknowledge_neighbors_substitution();
}

///
/// @brief
/// Returns the change in energy induced by changing a node from its current state into some alternate state for the PD energy terms only.
///
/// @details
/// This function always gets called for every substitution. Only the consider_alt_state() call can get procrastinated.
///
template < typename V, typename E, typename G >
core::PackerEnergy NPDHBondNode< V, E, G >::calculate_PD_deltaE_for_substitution(
	int alternate_state,
	core::PackerEnergy & prev_PDenergies_for_node
) {

	debug_assert( alternate_state > 0  && alternate_state <= parent::get_num_states() );

	prev_PDenergies_for_node = parent::get_curr_pd_energy_total();
	parent::calc_deltaEpd( alternate_state );

	return get_pd_energy_delta();
}

/// @brief
/// Returns the deltaE for just the PD terms. Separate method from the one above because this one can be called from within
/// a commit_sub call that didn't go through consider_sub().
template < typename V, typename E, typename G >
core::PackerEnergy NPDHBondNode< V, E, G >::get_pd_energy_delta() {
	return ( parent::get_alt_pd_energy_total() - parent::get_curr_pd_energy_total() );
}

/// @brief compute the non-pairwise-decomposable portion of the deltaE
/// @details This node must have already been informed what the alternate state is.
template < typename V, typename E, typename G >
Real
NPDHBondNode< V, E, G >::consider_alternate_state()
{
	// OK -- so this will proceed in three passes
	// First, we'll visit all the neighbors and let them know this node will be considering
	// a rotamer substitution
	//
	// Second, we'll visit all neighbors and calculate hbond energies to them,
	// and we'll have the neighbors compute new weights. Along the way, we'll
	// figure out which neighbors we're either gaining or losing hydrogen bonds
	// to; only these neighbors need to be iterated across in the second stage.

	// Third, we'll visit all the neighbors we either gained or lost a hydrogen bond
	// to, and we'll ask them to compute their delta( npd hbond energy ), excluding the
	// hydrogen bonds formed with this residue, and the hydrogen bonds formed to
	// "lower neighbors" to other non-changing nodes which either gained or lost a hydrogen
	// bond to this residue: we only want to calculate the new weighted energies for these
	// hydrogen bonds a single time.

	// Finally this node will calculate its own delta( npd hbond energy ) for the hydrogen
	// bonds it is forming with its neighbors

	// Step 0: initialization.
	utility::vector1< char > & hbonding_to_res = get_npd_hbond_owner()->hbonding_to_res_vector();

	// mark all of the residues that are hbonding to this residue already as having
	// possibly-changing energies as a result of this rotamer substitution.
	std::fill( hbonding_to_res.begin(), hbonding_to_res.end(), 0 );
	for ( Size ii = 1; ii <= curr_hbonds_.size(); ++ii ) {
		debug_assert( curr_hbonds_[ ii ]->don_rsd_ == seqpos_ || curr_hbonds_[ ii ]->acc_rsd_ == seqpos_ );
		hbonding_to_res[ curr_hbonds_[ ii ]->don_rsd_ == seqpos_ ? curr_hbonds_[ ii ]->acc_rsd_ : curr_hbonds_[ ii ]->don_rsd_ ] = 1;
	}

	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		get_adjacent_npd_hbond_bg_node( ii )->prepare_for_neighbors_substitution( seqpos_ );
	}
	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		get_adjacent_npd_hbond_node( ii )->prepare_for_neighbors_substitution( seqpos_ );
	}

	alt_hbonds_.resize( 0 );
	for ( Size ii = 1; ii <= max_natoms_; ++ii ) {
		alt_atom_hbonds_[ ii ].resize( 0 );
	}

	// Step 1: compute new hbonds and their weights
	//      1a: for the edges to background nodes
	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		get_edge_to_npd_hbond_bg_node( ii )->consider_alternate_state_step1(
			V::get_alternate_state(),
			alt_state_rotamer(),
			alt_hbonds_,
			alt_atom_hbonds_,
			hbonding_to_res
		);
	}

	//      1a: for the edges to first-class nodes
	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		get_incident_npd_hbond_edge( ii )->consider_alternate_state_step1(
			parent::get_node_index(),
			V::get_alternate_state(),
			alt_state_rotamer(),
			alt_hbonds_,
			alt_atom_hbonds_,
			hbonding_to_res );
	}

	compute_alt_weights_for_hbonds( false /* curr_state? no, for the alt_state*/ );

	// Step 2: sum the weighted deltaE
	//      2a: for the edges to background nodes
	Real delta_npd_hbond_energy( 0 );
	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		if ( hbonding_to_res[ get_adjacent_npd_hbond_bg_node( ii )->seqpos() ] ) {
			delta_npd_hbond_energy +=
				get_edge_to_npd_hbond_bg_node( ii )->consider_alternate_state_step2( hbonding_to_res );
		}
	}

	//      2b: for the edges to first-class nodes
	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		if ( hbonding_to_res[ get_adjacent_npd_hbond_node( ii )->seqpos() ] ) {
			delta_npd_hbond_energy +=
				get_incident_npd_hbond_edge( ii )->consider_alternate_state_step2( hbonding_to_res );
		}
	}

	// Finally, calculate the change for hbonds forming to this residue
	Real curr_hb_E( 0 ), alt_hb_E( 0 );
	for ( auto const & hb : curr_hbonds_ ) {
		//std::cout << " curr hb: " << hb->don_rsd_ << " " << hb->acc_rsd_ << " e: " << hb->energy_ << " wtdE: " << hb->energy_ * hb->sfxn_wt_ * hb->don_wt_ * hb->acc_wt_ << std::endl;
		curr_hb_E += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_ * hb->acc_wt_;
	}
	for ( auto const & hb : alt_hbonds_ ) {
		//std::cout << " alt hb: " << hb->don_rsd_ << " " << hb->acc_rsd_ << " e: " << hb->energy_ << " wtdE: " << hb->energy_ * hb->sfxn_wt_ * hb->don_wt_alt_ * hb->acc_wt_alt_ << std::endl;
		alt_hb_E += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_alt_ * hb->acc_wt_alt_;
	}
	delta_npd_hbond_energy += alt_hb_E - curr_hb_E;
	return delta_npd_hbond_energy;
}

template < typename V, typename E, typename G >
void NPDHBondNode< V, E, G >::prepare_for_neighbors_substitution( Size nbr_seqpos )
{
	alt_hbonds_.resize( 0 );
	for ( Size ii = 1; ii <= max_natoms_; ++ii ) {
		alt_atom_hbonds_[ ii ].resize( 0 );
	}

	for ( Size ii = 1; ii <= curr_hbonds_.size(); ++ii ) {
		NPDHBondOP const & iihb( curr_hbonds_[ ii ] );
		debug_assert( seqpos_ == iihb->don_rsd_ || seqpos_ == iihb->acc_rsd_ );

		// Do not add the hbonds that are formed to the changing neighbor
		if ( nbr_seqpos == iihb->don_rsd_ || nbr_seqpos == iihb->acc_rsd_ ) continue;

		iihb->don_wt_alt_ = iihb->don_wt_;
		iihb->acc_wt_alt_ = iihb->acc_wt_;

		alt_atom_hbonds_[ seqpos_ == iihb->don_rsd_ ? iihb->don_atm_ : iihb->acc_atm_ ].push_back( iihb );
		alt_hbonds_.push_back( iihb );

	}
}

template < typename V, typename E, typename G >
void NPDHBondNode< V, E, G >::compute_alt_weights_for_hbonds(
	bool curr_state
)
{
	if ( ( curr_state && parent::get_current_state() == 0) || ( ! curr_state && parent::get_alternate_state() == 0) ) return;
	//scoring::hbonds::NPDHBondSet const & hbset( get_npd_hbond_owner()->npd_hbond_set() );
	conformation::Residue const & res( curr_state ? curr_state_rotamer() : alt_state_rotamer() );
	compute_alt_weights_for_npd_hbonds( res, alt_atom_hbonds_, tmp_energies_, tmp_weights_ );
}

template < typename V, typename E, typename G >
utility::vector1< NPDHBondOP > &
NPDHBondNode< V, E, G >::current_hbs()
{
	return curr_hbonds_;
}

template < typename V, typename E, typename G >
utility::vector1< NPDHBondOP > &
NPDHBondNode< V, E, G >::alternate_hbs()
{
	return alt_hbonds_;
}

template < typename V, typename E, typename G >
utility::vector1< utility::vector1< NPDHBondOP > > &
NPDHBondNode< V, E, G >::alternate_hbs_for_atoms()
{
	return alt_atom_hbonds_;
}

template < typename V, typename E, typename G >
utility::vector1< NPDHBondOP > &
NPDHBondNode< V, E, G >::alternate_hbs_for_atom( Size atom_index )
{
	return alt_atom_hbonds_[ atom_index ];
}


template < typename V, typename E, typename G >
Real NPDHBondNode< V, E, G >::get_upper_npd_hbond_energy_totals() const
{
	Real etot = 0;
	for ( auto const & hb : curr_hbonds_ ) {
		debug_assert( hb->don_rsd_ == seqpos_ || hb->acc_rsd_ == seqpos_ );
		Size lower_seqpos = std::min( hb->don_rsd_, hb->acc_rsd_ );
		if ( lower_seqpos == seqpos_ ) {
			// i.e., The other residue in this upper residue in this hydrogen bond.
			etot += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_ * hb->acc_wt_;
		}
	}
	return etot;
}

/// @brief
/// Sets the current state to the alternate state this node was asked to
/// consider.  Copies appropriate score information.  Notifies all of its
/// neighbors that it is going through with the state substitution it had been
/// considering.
template < typename V, typename E, typename G >
void NPDHBondNode< V, E, G >::commit_considered_substitution() {

	debug_assert( parent::considering_alternate_state() );

	if ( parent::get_alternate_state() == 0 ) {
		assign_zero_state();
		return;
	}

	//call base class method
	parent::commit_considered_substitution();

	// TO DO!
	curr_hbonds_.swap( alt_hbonds_ );
	curr_atom_hbonds_.swap( alt_atom_hbonds_ );
	for ( auto const & hb : alt_hbonds_ ) {
		get_npd_hbond_owner()->return_hbond_to_queue( hb );
	}

	utility::vector1< char > & hbonding_to_res = get_npd_hbond_owner()->hbonding_to_res_vector();
	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		get_incident_npd_hbond_edge(ii)->acknowledge_substitution( hbonding_to_res[ get_adjacent_npd_hbond_node( ii )->seqpos() ]);

	}
	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		get_edge_to_npd_hbond_bg_node( ii )->acknowledge_substitution( hbonding_to_res[ get_adjacent_npd_hbond_bg_node( ii )->seqpos() ] );
	}

}


template < typename V, typename E, typename G >
void NPDHBondNode< V, E, G >::reset_after_rejected_substitution()
{
	// any new hydrogen bonds that NPDHBondInteractionGraph created for this
	// rotamer can now be recycled for reuse later.
	for ( auto const & hb : alt_hbonds_ ) {
		get_npd_hbond_owner()->return_hbond_to_queue( hb );
	}
	alt_hbonds_.resize( 0 );
}

///
/// @brief
/// Not implemented, but needs to be!
///
template < typename V, typename E, typename G >
unsigned int NPDHBondNode< V, E, G >::getMemoryUsageInBytes() const {
	return 0;
}

///
/// @brief
/// Returns the amount of static memory used by this Node object
///
template < typename V, typename E, typename G >
unsigned int NPDHBondNode< V, E, G >::count_static_memory() const {
	return sizeof ( NPDHBondNode< V, E, G > );
}

///
/// @brief
/// Returns the amount of dynamic memory used by this Node object
///
template < typename V, typename E, typename G >
unsigned int NPDHBondNode< V, E, G >::count_dynamic_memory() const {
	unsigned int total_memory = parent::count_dynamic_memory();

	// TO DO

	return total_memory;
}

///
/// @brief
/// useful for debugging
///
//template < typename V, typename E, typename G >
//void NPDHBondNode< V, E, G >::write_dot_kinemage( std::ofstream & output_kin ) {
// current_state_rotamer_dots_.write_dot_kinemage( output_kin );
//}

///
/// @brief
/// useful for debugging - writes information about a node to the tracer
///
template < typename V, typename E, typename G >
void NPDHBondNode< V, E, G >::print() const {
	parent::print();
	TR << "node " << parent::get_node_index() << ", current_state: " << parent::get_current_state()
		<< ", one body energy: " << parent::get_one_body_energy_current_state() << std::endl;
}


//----------------------------------------------------------------------------//
//----------------------- NPDHBond Background Residue Node Class -----------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// main constructor.  No default constructor, copy constructor or assignment operator
///
template < typename V, typename E, typename G >
NPDHBondBackgroundNode< V, E, G >::NPDHBondBackgroundNode( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int node_index ) :
	BackgroundNode< V, E, G > ( owner, node_index ),
	prepared_for_simA_( false )
{}

///
template < typename V, typename E, typename G >
NPDHBondBackgroundNode< V, E, G >::~NPDHBondBackgroundNode() {
}

///
/// @brief
/// inits the RotamerOP held by a background node. called in the NPDHBondIG::initialize method.
///
template < typename V, typename E, typename G >
void NPDHBondBackgroundNode< V, E, G >::set_rotamer( conformation::ResidueOP const & rotamer ) {
	rotamer_ = rotamer;
	seqpos_ = rotamer_->seqpos();
	curr_atom_hbonds_.resize( rotamer->natoms() );
	alt_atom_hbonds_.resize( rotamer->natoms() );
	curr_hbonds_.reserve( 10 );
	alt_hbonds_.reserve( 10 );
	for ( Size ii = 1; ii <= rotamer->natoms(); ++ii ) {
		curr_atom_hbonds_[ ii ].reserve( 6 );
		alt_atom_hbonds_[ ii ].reserve( 6 );
	}
}

template < typename V, typename E, typename G >
conformation::Residue const &
NPDHBondBackgroundNode< V, E, G >::get_rotamer() const {
	return *rotamer_;
}

template < typename V, typename E, typename G >
conformation::ResidueCOP NPDHBondBackgroundNode< V, E, G >::get_rotamer_op() const {
	return rotamer_;
}

template < typename V, typename E, typename G >
Size NPDHBondBackgroundNode< V, E, G >::seqpos() const {
	return seqpos_;
}


template < typename V, typename E, typename G >
void NPDHBondBackgroundNode< V, E, G >::prepare_for_simulated_annealing() {

	if ( ! prepared_for_simA_ ) {
		if ( ! parent::get_edge_vector_up_to_date() ) {
			parent::update_edge_vector();
		}

		prepared_for_simA_ = true;
	}

}


template < typename V, typename E, typename G >
void NPDHBondBackgroundNode< V, E, G >::prepare_for_neighbors_substitution( Size nbr_seqpos )
{
	alt_hbonds_.resize( 0 );
	for ( Size ii = 1; ii <= rotamer_->natoms(); ++ii ) {
		alt_atom_hbonds_[ ii ].resize( 0 );
	}

	for ( Size ii = 1; ii <= curr_hbonds_.size(); ++ii ) {
		NPDHBondOP const & iihb( curr_hbonds_[ ii ] );
		debug_assert( seqpos_ == iihb->don_rsd_ || seqpos_ == iihb->acc_rsd_ );

		// Do not add the hbonds that are formed to the changing neighbor
		if ( nbr_seqpos == iihb->don_rsd_ || nbr_seqpos == iihb->acc_rsd_ ) continue;

		iihb->don_wt_alt_ = iihb->don_wt_;
		iihb->acc_wt_alt_ = iihb->acc_wt_;

		alt_atom_hbonds_[ seqpos_ == iihb->don_rsd_ ? iihb->don_atm_ : iihb->acc_atm_ ].push_back( iihb );
		alt_hbonds_.push_back( iihb );

	}
}

template < typename V, typename E, typename G >
void NPDHBondBackgroundNode< V, E, G >::compute_alt_weights_for_hbonds()
{
	//scoring::hbonds::NPDHBondSet const & hbset( get_npd_hbond_owner()->npd_hbond_set() );
	compute_alt_weights_for_npd_hbonds( *rotamer_, alt_atom_hbonds_, tmp_energies_, tmp_weights_ );
}

template < typename V, typename E, typename G >
utility::vector1< NPDHBondOP > &
NPDHBondBackgroundNode< V, E, G >::current_hbs()
{
	return curr_hbonds_;
}

template < typename V, typename E, typename G >
utility::vector1< NPDHBondOP > &
NPDHBondBackgroundNode< V, E, G >::alternate_hbs()
{
	return alt_hbonds_;
}

template < typename V, typename E, typename G >
utility::vector1< utility::vector1< NPDHBondOP > > &
NPDHBondBackgroundNode< V, E, G >::alternate_hbs_for_atoms()
{
	return alt_atom_hbonds_;
}

template < typename V, typename E, typename G >
utility::vector1< NPDHBondOP > &
NPDHBondBackgroundNode< V, E, G >::alternate_hbs_for_atom( Size atom_index )
{
	return alt_atom_hbonds_[ atom_index ];
}

template < typename V, typename E, typename G >
Real NPDHBondBackgroundNode< V, E, G >::get_upper_npd_hbond_energy_totals() const
{
	Real etot = 0;
	for ( auto const & hb : curr_hbonds_ ) {
		debug_assert( hb->don_rsd_ == seqpos_ || hb->acc_rsd_ == seqpos_ );
		Size lower_seqpos = std::min( hb->don_rsd_, hb->acc_rsd_ );
		if ( lower_seqpos == seqpos_ ) {
			// i.e., The other residue in this upper residue in this hydrogen bond.
			etot += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_ * hb->acc_wt_;
		}
	}
	return etot;
}



///
/// @brief
/// bookkeeping to reflect a neighboring NPDHBondNode's state substitution
///
template < typename V, typename E, typename G >
void NPDHBondBackgroundNode< V, E, G >::acknowledge_substitution() {
	curr_hbonds_.swap( alt_hbonds_ );
	curr_atom_hbonds_.swap( alt_atom_hbonds_ );
	for ( auto const & hb : curr_hbonds_ ) {
		hb->don_wt_ = hb->don_wt_alt_;
		hb->acc_wt_ = hb->acc_wt_alt_;
	}
}

template < typename V, typename E, typename G >
void NPDHBondBackgroundNode< V, E, G >::acknowledge_neighbors_state_zeroed( Size seqpos )
{
	// 1. Purge the hbonds to the neighbor from the alt_* set
	// 2. Compute a new set of weights for the remaining hydrogen bonds
	// 3. Then accept the alternate as the current.

	prepare_for_neighbors_substitution( seqpos );
	compute_alt_weights_for_hbonds();
	acknowledge_substitution();
}

template < typename V, typename E, typename G >
unsigned int NPDHBondBackgroundNode< V, E, G >::count_static_memory() const {
	return sizeof ( NPDHBondBackgroundNode< V, E, G > );
}

///
template < typename V, typename E, typename G >
unsigned int NPDHBondBackgroundNode< V, E, G >::count_dynamic_memory() const {
	unsigned int total_memory = parent::count_dynamic_memory();

	return total_memory;
}

///
/// @brief
/// used only for debugging
template < typename V, typename E, typename G >
void NPDHBondBackgroundNode< V, E, G >::print() const {
	TR << "bgnode " << parent::get_node_index();
}


//----------------------------------------------------------------------------//
//------------------------------ NPDHBond Edge Class -----------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// main constructor.  No default, or copy constructors, no assignment operator
///
/// @param
/// owner - [in] - the owning interaction graph object
/// node1 - [in] - the index of the lower-indexed NPDHBondNode
/// node2 - [in] - the index of the higher-indexed NPDHBondNode
///
template < typename V, typename E, typename G >
NPDHBondEdge< V, E, G >::NPDHBondEdge( G* owner, int node1, int node2 ) :
	FirstClassEdge< V, E, G > ( owner, node1, node2 ),
	node_changing_( -1 ),
	node_not_changing_( -1 )
{
	nodes_curr_states_[ 0 ] = nodes_curr_states_[ 1 ] = 0;
	nodes_alt_states_[ 0 ] = nodes_alt_states_[ 1 ] = 0;
}

///
template < typename V, typename E, typename G >
NPDHBondEdge< V, E, G >::~NPDHBondEdge() {}

///
/// @brief
/// drops zero submatrices of the AminoAcidNeighborSparseMatrix and if the two_body_energies_ member then holds nothing,
template < typename V, typename E, typename G >
void NPDHBondEdge< V, E, G >::prepare_for_simulated_annealing() {

	parent::prepare_for_simulated_annealing_no_deletion(); // AddtlBGNodesIG doesn't have an Edge method for prep for simA, but PDIG does

	if ( parent::pd_edge_table_all_zeros() ) {
		// TO DO: decide when to delete this edge
		// delete this;
	}
}

///
/// @brief
/// respond to when one of its vertices enters the "unassigned" state.
///
/// @details
/// called during the NPDHBIG::blanket_assign_state_0 -> NPDHBondNode::assign_zero_state() cascade of calls.
///
template < typename V, typename E, typename G >
void NPDHBondEdge< V, E, G >::acknowledge_state_zeroed( int node_that_changed, Size node_seqpos )
{
	// each Edge contains a node_changing_ and node_not_changing_ int (which takes on values of 0 or 1)
	// node_changing_ is the node of this edge that is changing
	node_changing_ = ( node_that_changed == parent::get_node_index(0) ? 0 : 1 );
	node_not_changing_ = ! node_changing_;

	nodes_curr_states_[ node_changing_ ] = 0;

	get_npd_hbond_node( node_not_changing_ )->acknowledge_neighbors_state_zeroed( node_seqpos );
}

///
/// @brief
/// bookkeeping following the decision to substitute a nodes current state with the alternate it was asked to consider.
///
template < typename V, typename E, typename G >
void NPDHBondEdge< V, E, G >::acknowledge_substitution( bool update_hbonds ) {
	if ( update_hbonds ) {
		get_npd_hbond_node( node_not_changing_ )->acknowledge_neighbors_substitution();
	}
	nodes_curr_states_[ node_changing_ ] = nodes_alt_states_[ node_changing_ ];
}


///
/// @brief
/// Reduces memory usage in the two body energy table after the energy
/// calculating function declares that the energies will not change thereafter
///
/// @remarks (all by apl)
/// In the PDEdge's version of this method, after invoking two_body_energies_.drop_zero_submatrices_where_possible();
/// the PDEdge checks if the two body energy table it now holds is empty.  If the table is empty, the edge deletes itself.
template < typename V, typename E, typename G >
void NPDHBondEdge< V, E, G >::declare_energies_final() {
	parent::declare_energies_final_no_deletion();
}

///
/// @remarks
/// Not implemented.
///
template < typename V, typename E, typename G >
unsigned int NPDHBondEdge< V, E, G >::getMemoryUsageInBytes() const {
	return 0;
}

///
template < typename V, typename E, typename G >
unsigned int NPDHBondEdge< V, E, G >::count_static_memory() const {
	return sizeof ( NPDHBondEdge< V, E, G > );
}

///
template < typename V, typename E, typename G >
unsigned int NPDHBondEdge< V, E, G >::count_dynamic_memory() const {

	unsigned int total_memory = parent::count_dynamic_memory();
	// TO DO!

	return total_memory;
}

// template < typename V, typename E, typename G >
// Real NPDHBondEdge< V, E, G >::get_current_two_body_energy() const
// {
//  return current_two_body_energy_;
// }

template < typename V, typename E, typename G >
void
create_hbonds_one_way(
	scoring::hbonds::HBondDatabase const & database,
	scoring::hbonds::HBondOptions const & hbondoptions,
	scoring::hbonds::HBondSet const & hbset,
	NPDHBondInteractionGraph< V, E, G > & ig,
	utility::vector1< char > & hbonding_to_res,
	conformation::Residue const & acc_res,
	utility::vector1< NPDHBondOP > & acc_hbonds,
	utility::vector1< utility::vector1< NPDHBondOP > > & acc_atom_hbonds,
	conformation::Residue const & don_res,
	utility::vector1< NPDHBondOP > & don_hbonds,
	utility::vector1< utility::vector1< NPDHBondOP > > & don_atom_hbonds
)
{
	for ( Size ii_at : acc_res.accpt_pos() ) {
		//Size ii_at = acc_res.accpt_pos()[ ii ];
		Vector const & iixyz = acc_res.xyz( ii_at );
		for ( Size jj_at : don_res.Hpos_polar() ) {
			//Size jj_at = don_res.Hpos_polar()[ jj ];
			Vector const & jjxyz = don_res.xyz( jj_at );
			if ( iixyz.distance_squared( jjxyz ) > core::scoring::hbonds::MAX_R2 ) continue;
			Real e = core::scoring::hbonds::hb_energy( database, hbondoptions, hbset, acc_res, ii_at, don_res, jj_at );
			if ( e >= core::scoring::hbonds::MAX_HB_ENERGY ) continue;
			NPDHBondOP hb = ig.unused_hbond();
			hb->energy_ = e;
			hb->don_wt_ = -1234;
			hb->acc_wt_ = -1234;
			hb->don_wt_alt_ = -1234;
			hb->acc_wt_alt_ = -1234;
			hb->don_rsd_ = don_res.seqpos();
			hb->acc_rsd_ = acc_res.seqpos();
			hb->don_atm_ = jj_at;
			hb->acc_atm_ = ii_at;

			// score function weight should be included in the donor and acceptor weight calculation
			scoring::hbonds::HBEvalTuple hb_eval_tup( don_res.atom_base( jj_at ), don_res, ii_at, acc_res);
			hb->sfxn_wt_ = ig.npd_hb_weight( hb_eval_tup.eval_type(), don_res.seqpos() == acc_res.seqpos());

			acc_hbonds.push_back( hb );
			don_hbonds.push_back( hb );
			acc_atom_hbonds[ ii_at ].push_back( hb );
			don_atom_hbonds[ jj_at ].push_back( hb );
			hbonding_to_res[ don_res.seqpos() ] = 1;
			hbonding_to_res[ acc_res.seqpos() ] = 1;
		}
	}
}


template < typename V, typename E, typename G >
void NPDHBondEdge< V, E, G >::consider_alternate_state_step1(
	int node_index,
	int alternate_state_index,
	conformation::Residue const & alt_state,
	utility::vector1< NPDHBondOP > & res_hbonds,
	utility::vector1< utility::vector1< NPDHBondOP > > & atom_hbonds,
	utility::vector1< char > & hbonding_to_res
)
{
	node_changing_ = parent::which_node( node_index );
	nodes_alt_states_[ node_changing_ ] = alternate_state_index;

	node_not_changing_ =  ! node_changing_;

	// If the other node is assigned state zero, quit
	if ( nodes_curr_states_[ node_not_changing_ ] == 0 ) return;

	conformation::Residue const & unchanging_curr_state( get_npd_hbond_node( node_not_changing_ )->curr_state_rotamer() );
	NPDHBondNode< V, E, G > * unchanging_node = get_npd_hbond_node( node_not_changing_ );
	utility::vector1< utility::vector1< NPDHBondOP > > & unchanging_atom_hbonds =
		unchanging_node->alternate_hbs_for_atoms();


	NPDHBondInteractionGraph< V, E, G > & ig( *get_npd_hbond_owner() );
	scoring::hbonds::HBondDatabase const & database( ig.hbond_database() );
	scoring::hbonds::HBondOptions const & hbondoptions( ig.hbond_options() );
	scoring::hbonds::NPDHBondSet const & hbset( ig.npd_hbond_set() );

	create_hbonds_one_way( database, hbondoptions, hbset, ig, hbonding_to_res,
		alt_state, res_hbonds, atom_hbonds,
		unchanging_curr_state, unchanging_node->alternate_hbs(), unchanging_atom_hbonds );

	create_hbonds_one_way( database, hbondoptions, hbset, ig, hbonding_to_res,
		unchanging_curr_state, unchanging_node->alternate_hbs(), unchanging_atom_hbonds,
		alt_state, res_hbonds, atom_hbonds );

	if ( ! hbonding_to_res[ unchanging_curr_state.seqpos() ] ) return;

	// Now we derive weights for the hbonds on the unchanging node's side.
	get_npd_hbond_node( node_not_changing_ )->compute_alt_weights_for_hbonds( true /* curr_state */ );
}

template < typename V, typename E, typename G >
Real NPDHBondEdge< V, E, G >::consider_alternate_state_step2(
	utility::vector1< char > const & hbonding_to_res
)
{
	Size seqpos_changing = get_npd_hbond_node( node_changing_ )->seqpos();
	Size seqpos_unchanging = get_npd_hbond_node( node_not_changing_ )->seqpos();

	Real orig_hbE( 0 ), alt_hbE( 0 );
	for ( auto const & hb : get_npd_hbond_node( node_not_changing_ )->current_hbs() ) {
		debug_assert( hb->don_rsd_ == seqpos_unchanging || hb->acc_rsd_ == seqpos_unchanging );
		if ( hb->don_rsd_ == seqpos_changing || hb->acc_rsd_ == seqpos_changing   ) { continue; }

		if (
				( hb->don_rsd_ == seqpos_unchanging && ( ! hbonding_to_res[ hb->acc_rsd_ ] || seqpos_unchanging <= hb->acc_rsd_ )) ||
				( hb->acc_rsd_ == seqpos_unchanging && ( ! hbonding_to_res[ hb->don_rsd_ ] || seqpos_unchanging <= hb->don_rsd_ )) ) {
			// count upper hydrogen bonds or hbonds to the unaffected background to ensure that hbonds
			// are counted once and only once
			//std::cout << " orig hb: " << hb->don_rsd_ << " " << hb->acc_rsd_ << " e: " << hb->energy_ << " wtdE: " << hb->energy_ * hb->sfxn_wt_ * hb->don_wt_ * hb->acc_wt_ << std::endl;
			orig_hbE += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_ * hb->acc_wt_;
		}
	}

	for ( auto const & hb : get_npd_hbond_node( node_not_changing_ )->alternate_hbs() ) {
		debug_assert( hb->don_rsd_ == seqpos_unchanging || hb->acc_rsd_ == seqpos_unchanging );
		if ( hb->don_rsd_ == seqpos_changing || hb->acc_rsd_ == seqpos_changing   ) { continue; }

		if (
				( hb->don_rsd_ == seqpos_unchanging && ( ! hbonding_to_res[ hb->acc_rsd_ ] || seqpos_unchanging <= hb->acc_rsd_ )) ||
				( hb->acc_rsd_ == seqpos_unchanging && ( ! hbonding_to_res[ hb->don_rsd_ ] || seqpos_unchanging <= hb->don_rsd_ )) ) {
			// count upper hydrogen bonds or hbonds to the unaffected background to ensure that hbonds
			// are counted once and only once
			// std::cout << " alt hb: " << hb->don_rsd_ << " " << hb->acc_rsd_ << " e: " << hb->energy_ << " wtdE: " << hb->energy_ * hb->sfxn_wt_ * hb->don_wt_alt_ * hb->acc_wt_alt_ << std::endl;
			alt_hbE += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_alt_ * hb->acc_wt_alt_;
		}
	}

	return alt_hbE - orig_hbE;
}

//----------------------------------------------------------------------------//
//----------------------- NPDHBondBackgroundEdge Class ---------------------------//
//----------------------------------------------------------------------------//

///
/// @brief
/// main constructor
///
template < typename V, typename E, typename G >
NPDHBondBackgroundEdge< V, E, G >::NPDHBondBackgroundEdge(
	AdditionalBackgroundNodesInteractionGraph < V, E, G >* owner,
	int first_class_node_index,
	int background_node_index
) :
	BackgroundToFirstClassEdge< V, E, G >( owner, first_class_node_index, background_node_index ),
	prepared_for_simA_( false ),
	nodes_curr_state_( 0 ),
	nodes_alt_state_( 0 )
{}

///
template < typename V, typename E, typename G >
NPDHBondBackgroundEdge< V, E, G >::~NPDHBondBackgroundEdge() {}

template < typename V, typename E, typename G >
void NPDHBondBackgroundEdge< V, E, G >::prepare_for_simulated_annealing() {}


///////
/////// @brief
/////// bookkeeping in response to a NPDHBondNode switching states (without having gone through the usual
/////// consider-substitution/commit-substitution pattern).
////template < typename V, typename E, typename G >
////void NPDHBondBackgroundEdge< V, E, G >::acknowledge_state_change( int new_state ) {
////
//// if ( new_state == nodes_curr_state_ ) { // in the case of the 0-state, just return - don't calculate deltaE
////  return;
//// }
////
//// update_state_at_neighbor( new_state );
//// acknowledge_substitution();
////
////}

///
/// @brief
/// respond to when one of its vertices enters the "unassigned" state.
///
/// @details
/// called during the NPDHBIG::blanket_assign_state_0 -> NPDHBondNode::assign_zero_state() cascade of calls.
///
template < typename V, typename E, typename G >
void NPDHBondBackgroundEdge< V, E, G >::acknowledge_state_zeroed( Size node_seqpos )
{
	get_npd_hbond_bg_node()->acknowledge_neighbors_state_zeroed( node_seqpos );
}

///
/// @brief
/// bookkeeping in response to the NPDHBondNode committing the considered substitution
///
template < typename V, typename E, typename G >
void NPDHBondBackgroundEdge< V, E, G >::acknowledge_substitution( bool update_hbonds ) {
	if ( update_hbonds ) {
		get_npd_hbond_bg_node()->acknowledge_substitution();
	}
	nodes_curr_state_ = nodes_alt_state_;
}

///
template < typename V, typename E, typename G >
unsigned int NPDHBondBackgroundEdge< V, E, G >::count_static_memory() const {
	return sizeof ( NPDHBondBackgroundEdge< V, E, G > );
}

///
template < typename V, typename E, typename G >
unsigned int NPDHBondBackgroundEdge< V, E, G >::count_dynamic_memory() const {

	unsigned int total_memory = parent::count_dynamic_memory();
	// TO DO
	return total_memory;
}

template < typename V, typename E, typename G >
void NPDHBondBackgroundEdge< V, E, G >::consider_alternate_state_step1(
	int state_index,
	conformation::Residue const & alt_state,
	utility::vector1< NPDHBondOP > & res_hbonds,
	utility::vector1< utility::vector1< NPDHBondOP > > & atom_hbonds,
	utility::vector1< char > & hbonding_to_res
)
{
	nodes_alt_state_ = state_index;
	NPDHBondBackgroundNode< V, E, G > * unchanging_node = get_npd_hbond_bg_node();
	conformation::Residue const & unchanging_curr_state( unchanging_node->get_rotamer() );
	utility::vector1< utility::vector1< NPDHBondOP > > & unchanging_atom_hbonds =
		unchanging_node->alternate_hbs_for_atoms();


	NPDHBondInteractionGraph< V, E, G > & ig( *get_npd_hbond_owner() );
	scoring::hbonds::HBondDatabase const & database( ig.hbond_database() );
	scoring::hbonds::HBondOptions const & hbondoptions( ig.hbond_options() );
	scoring::hbonds::NPDHBondSet const & hbset( ig.npd_hbond_set() );

	create_hbonds_one_way( database, hbondoptions, hbset, ig, hbonding_to_res,
		alt_state, res_hbonds, atom_hbonds,
		unchanging_curr_state, unchanging_node->alternate_hbs(), unchanging_atom_hbonds );

	create_hbonds_one_way( database, hbondoptions, hbset, ig, hbonding_to_res,
		unchanging_curr_state, unchanging_node->alternate_hbs(), unchanging_atom_hbonds,
		alt_state, res_hbonds, atom_hbonds );

	if ( ! hbonding_to_res[ unchanging_node->seqpos() ] ) return;

	// Now we derive weights for the hbonds on the unchanging node's side.
	unchanging_node->compute_alt_weights_for_hbonds();

}

template < typename V, typename E, typename G >
Real NPDHBondBackgroundEdge< V, E, G >::consider_alternate_state_step2(
	utility::vector1< char > const & hbonding_to_res
)
{

	Size seqpos_changing = get_npd_hbond_node()->seqpos();
	Size seqpos_unchanging = get_npd_hbond_bg_node()->seqpos();

	Real orig_hbE( 0 ), alt_hbE( 0 );
	for ( auto const & hb : get_npd_hbond_bg_node()->current_hbs() ) {
		debug_assert( hb->don_rsd_ == seqpos_unchanging || hb->acc_rsd_ == seqpos_unchanging );
		if ( hb->don_rsd_ == seqpos_changing || hb->acc_rsd_ == seqpos_changing   ) { continue; }

		if (
				( hb->don_rsd_ == seqpos_unchanging && ( ! hbonding_to_res[ hb->acc_rsd_ ] || seqpos_unchanging <= hb->acc_rsd_ )) ||
				( hb->acc_rsd_ == seqpos_unchanging && ( ! hbonding_to_res[ hb->don_rsd_ ] || seqpos_unchanging <= hb->don_rsd_ )) ) {
			// count upper hydrogen bonds or hbonds to the unaffected background to ensure that hbonds
			// are counted once and only once
			// std::cout << " orig bg hb: " << hb->don_rsd_ << " " << hb->acc_rsd_ << " e: " << hb->energy_ << " wtdE: " << hb->energy_ * hb->sfxn_wt_ * hb->don_wt_ * hb->acc_wt_ << std::endl;
			orig_hbE += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_ * hb->acc_wt_;
		}
	}

	for ( auto const & hb : get_npd_hbond_bg_node()->alternate_hbs() ) {
		debug_assert( hb->don_rsd_ == seqpos_unchanging || hb->acc_rsd_ == seqpos_unchanging );
		if ( hb->don_rsd_ == seqpos_changing || hb->acc_rsd_ == seqpos_changing   ) { continue; }
		if (
				( hb->don_rsd_ == seqpos_unchanging && ( ! hbonding_to_res[ hb->acc_rsd_ ] || seqpos_unchanging <= hb->acc_rsd_ )) ||
				( hb->acc_rsd_ == seqpos_unchanging && ( ! hbonding_to_res[ hb->don_rsd_ ] || seqpos_unchanging <= hb->don_rsd_ )) ) {
			// count upper hydrogen bonds or hbonds to the unaffected background to ensure that hbonds
			// are counted once and only once
			// std::cout << " alt bg hb: " << hb->don_rsd_ << " " << hb->acc_rsd_ << " e: " << hb->energy_ << " wtdE: " << hb->energy_ * hb->sfxn_wt_ * hb->don_wt_alt_ * hb->acc_wt_alt_ << std::endl;
			alt_hbE += hb->energy_ * hb->sfxn_wt_ * hb->don_wt_alt_ * hb->acc_wt_alt_;
		}
	}

	return alt_hbE - orig_hbE;
}

//----------------------------------------------------------------------------//
//--------------------------- NPDHBond Interaction Graph -------------------------//
//----------------------------------------------------------------------------//

template < typename V, typename E, typename G >
Size NPDHBondInteractionGraph< V, E, G >::num_state_substitutions_considered_ = 0;

template < typename V, typename E, typename G >
Size NPDHBondInteractionGraph< V, E, G >::num_npd_hbond_comps_procrastinated_ = 0;

template < typename V, typename E, typename G >
Size NPDHBondInteractionGraph< V, E, G >::num_npd_hbond_comps_later_made_ = 0;

///
/// @details
/// Main constructor. Initializes all member variables to 0 and false.
///
template < typename V, typename E, typename G >
NPDHBondInteractionGraph< V, E, G >::NPDHBondInteractionGraph( int num_nodes ) :
	AdditionalBackgroundNodesInteractionGraph< V, E, G > ( num_nodes )
{}

///
template < typename V, typename E, typename G >
NPDHBondInteractionGraph< V, E, G >::~NPDHBondInteractionGraph() {}

///
/// @details
/// All throughout this class, I refer back to the original pose sequence. To be able to do that, I need to have a
/// handle to the pose in this class.  That's what this method provides.  In IGFactory.cc, this method gets called with
/// the pose object that's being packed/designed.
///
template < typename V, typename E, typename G >
void
NPDHBondInteractionGraph<V, E, G>::set_pose( pose::Pose const & pose ) {

	// call the set_pose function in the LinMemIG class, because it, too, uses the Pose to do its thing
	if ( typeid(G) == typeid( pack::interaction_graph::LinearMemoryInteractionGraph ) ) {
		dynamic_cast<pack::interaction_graph::LinearMemoryInteractionGraph*>(this)->set_pose( pose );
	}

	pose_ = pose::PoseOP( new pose::Pose( pose ) );
	hbonding_to_res_.resize( pose_->total_residue() );
}

///
/// @details
/// All throughout this class, I refer back to the original pose sequence. To be able to do that, I need to have a
/// handle to the pose in this class.  That's what this method provides.  In IGFactory.cc, this method gets called with
/// the pose object that's being packed/designed.
///
template < typename V, typename E, typename G >
void
NPDHBondInteractionGraph<V, E, G>::set_score_function( scoring::ScoreFunction const & sfxn ) {
	debug_assert( pose_ );
	using scoring::EnergiesCacheableDataType::NPD_HBOND_SET;
	using namespace scoring::hbonds;

	// call the set_pose function in the LinMemIG class, because it, too, uses the ScoreFunction to do its thing
	if ( typeid(G) == typeid( pack::interaction_graph::LinearMemoryInteractionGraph ) ) {
		(dynamic_cast< LinearMemoryInteractionGraph * > ( this ))->set_score_function( sfxn );
	}

	scorefxn_ = sfxn.clone();

	npd_hbond_set_ = utility::pointer::static_pointer_cast< NPDHBondSet > ( pose_->energies().data().get_ptr( NPD_HBOND_SET ) );
	if ( ! npd_hbond_set_ ) {
		npd_hbond_set_ = NPDHBondSetOP( new NPDHBondSet );
	}

	hbond_options_ = HBondOptionsOP( new HBondOptions( npd_hbond_set_->hbond_options() ));
	hbond_database_ = HBondDatabase::get_database(hbond_options_->params_database_tag());


}

template < typename V, typename E, typename G >
void
NPDHBondInteractionGraph<V, E, G>::set_packer_neighbor_graph( utility::graph::Graph const & neighbor_graph )
{
	using utility::graph::Graph;
	using utility::graph::GraphOP;
	neighbor_graph_ = GraphOP( new Graph( neighbor_graph ) );
}

/// @brief
/// We need a copy of the packer task to figure out which residues are being packed and/or designed. We have to figure
/// the packing options because it determines whether a residue becomes a FirstClass (NPDHBondNode) node or a background node.
/// This method gets called in IGSupport.cc.
template < typename V, typename E, typename G >
void
NPDHBondInteractionGraph<V, E, G>::set_packer_task( task::PackerTask const & the_task ) {
	packer_task_ = the_task.clone();
}

template < typename V, typename E, typename G >
void
NPDHBondInteractionGraph<V, E, G>::set_rotamer_sets( rotamer_set::RotamerSets const & rotsets ) {
	rotamer_sets_ = rotamer_set::RotamerSetsOP( new rotamer_set::RotamerSets( rotsets ) );
}

///
/// @details
/// This function is the 1st major entry point (well, after the constructor) into the NPDHBIG. It needs to set residues that
/// are not changing as background residues
///
/// Oldest comments:
/// In ++, there's a file InteractionGraphSupport.cc which is an analog of the InteractionGraphFactory in mini.  In ++,
/// the InteractionGraphSupport file instantiates and initializes, depending on the command line switches, the right
/// interaction graph.  For the NPDHBondInteractionGraph, it first initializes the PDInteractionGraph (which is the base)
/// and then calls the NPDHBondIG initialize method.
///
/// The thing is that this initialize method can't be called when the graph is constructed in the InteractionGraphFactory.
/// The reason is that the PDInteractionGraph base initialize() method is NOT called until later in the pack rotamers
/// process.  (Actually it's called from within rotsets->compute_energies().)  Only after the rotsets->compute energies
/// returns can I stick in an initialize() method call for the NPDHBondInteractionGraph (NPDHBIG).  But then that's too late because
/// the rotsets object has already computed some energies. Perhaps that's ok though. The rotsets object only calculates
/// the PD energy terms - it doesn't do anything with non-PD terms.
///
/// If a NPDHBIG init method is called during construction of the NPDHBIG, then the init method that's called by the rotsets object
/// causes all the node info that was set in the NPDHBIG init method to be lost because the rotsets init method recreates all
/// the nodes in the interaction graph when it runs. (That was a fun behaviour to figure out.)
///
/// So the solution I'm going for is to call this init method in the prepare for simulated annealing method of this class.
/// That gets called just before SA starts, so it will do the task above then.  It doesn't really matter when the
/// task gets done as long as it happens before SA starts.  This also means that the NPDHBIG will now have to keep a reference
/// to the Pose, the Task, and the RotamerSets objects since it needs all of these things to do tasks 1) and 2).
/// (For the port of this NPDHBondIG, we might not need the task and rotamer sets objects.)
///
/// prepare_for_simulated_annealing gets called by the FixbbSA::run() method.  Before this method, the
/// rotamersets object has called compute_energies() (the whole process being started in pack_rotamers)
/// which calls initialize() on the IG.  I need to place the NPDHBIG init method directly after the IG init
/// method that the RS object calls.
///
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::initialize( rotamer_set::RotamerSetsBase const & rot_sets_base ) {

	/// TEMP!!!!
	// APL wants to create a parent/base class for RotamerSets called RotamerSetsBase which will hold variables and functions
	// that are common to RotamerSets and FlexbbRotamerSets. Each interaction graph will be passed a RotamerSetsBase object
	// to its initialize method, and functions in this class should only access common base-class functions.  However,
	// this would require changing quite a bit of code and making sure that only base-class RotamerSetsBase methods are used
	// here which I'm not interested in doing right now. So we can cast the RotamerSetsBase object to a RotamerSets object
	// and use it the way it was previously.
	rotamer_set::RotamerSets const & rot_sets( static_cast< rotamer_set::RotamerSets const & > (rot_sets_base) );

	G::initialize( rot_sets );

	// save references to the rotamer_set. these get used later to print information about considered subs, for example.
	for ( Size ii = 1; ii <= rotamer_sets().nmoltenres(); ++ii ) {
		get_npd_hbond_node( ii )->set_rotamers( rotamer_sets().rotamer_set_for_moltenresidue( ii ) );
	}

	// initializes some local variables that translate between the ig enumeration and resid
	set_num_residues_in_protein( pose().size() );

	int nbackground = pose().size() - rot_sets.nmoltenres();
	set_num_background_residues( nbackground );

	Size count_bg( 0 );
	for ( Size ii = 1; ii <= pose().size(); ++ii ) {

		// in our case, first class residues are residues that are designable or packable (see notes in function comment).
		if ( packer_task().being_packed(ii) || packer_task().being_designed(ii) ) {
			// it's probably that being_packed() includes being_designed() so that the conditional needs only to say
			// being_packed(), but it'll short circuit the OR if being_packed() is true so it doesn't matter too much.
			continue;
		} else {
			using core::conformation::Residue;
			using core::conformation::ResidueOP;
			++count_bg;
			set_residue_as_background_residue( ii );
			get_npd_hbond_bg_node( count_bg )->set_rotamer( ResidueOP( new Residue( pose().residue( ii ) )));
		}
	}

	// Now let's add edges to bg residues
	for ( Size ii = 1; ii <= pose().size(); ++ii ) {
		if ( packer_task().being_packed(ii) || packer_task().being_designed(ii) ) {
			// iterate across all of the neighbors, and if they are background residues, then
			// add an edge
			Size ii_node_id = rotamer_sets().resid_2_moltenres( ii );
			for ( auto edge_iter = neighbor_graph_->get_node( ii )->const_edge_list_begin();
					edge_iter != neighbor_graph_->get_node( ii )->const_edge_list_end(); ++edge_iter ) {
				Size jj = (*edge_iter)->get_other_ind( ii );
				if ( ! packer_task().being_packed( jj ) && ! packer_task().being_designed( jj ) ) {
					Size jj_bgnode_id = resid_2_bgenumeration_[ jj ];
					parent::add_background_edge( ii_node_id, jj_bgnode_id );
				}
			}
		}
	}

	for ( Size ii = 1; ii <= rotamer_sets().nmoltenres(); ++ii ) {
		// when G::initialize is called, we recreate all of the FC Nodes so that index ii in RotamerSets maps to FCNode ii.
		rotamer_set::RotamerSetCOP rsop = rotamer_sets().rotamer_set_for_moltenresidue( ii );
		// TO DO: give the FC node the rotamers?
	}

}

/// @brief
/// tells the graph how many residues there are total in the protein
///
/// @details
/// The graph maintains its own enumeration for the background residues; but asks that anyone wanting to refer to them
/// use their original resid. The graph has to switch back and forth between enumeration schemes and must know how many
/// residues there are total to do that efficiently.
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::set_num_residues_in_protein( Size num_res ) {

	num_total_residues_ = num_res;
	resid_2_bgenumeration_.resize( num_total_residues_ );

	for ( Size ii = 1; ii <= num_total_residues_; ++ii ) {
		resid_2_bgenumeration_[ii] = 0;
	}

}

/// @brief
/// tells the graph how many residues there are as part of the protein that are not part of the combinatorial
/// optimization process -- they are part of the background
///
/// @details
/// The other half of switching between enumeration schemes for the background residues is knowing how many background residues there are.
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::set_num_background_residues( Size num_background_residues ) {

	parent::set_num_background_nodes( num_background_residues );
	if ( parent::get_num_background_nodes() == 0 ) {
		return;
	}

	bgenumeration_2_resid_.resize( parent::get_num_background_nodes() );
	for ( Size ii = 1; ii <= (Size)parent::get_num_background_nodes(); ++ii ) {
		bgenumeration_2_resid_[ii] = 0;
	}
}

/// @brief
/// informs the graph that a particular residue is part of the background
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::set_residue_as_background_residue( int residue ) {

	debug_assert( resid_2_bgenumeration_[ residue ] == 0 );

	++num_residues_assigned_as_background_;
	resid_2_bgenumeration_[ residue ] = num_residues_assigned_as_background_;
	bgenumeration_2_resid_[ num_residues_assigned_as_background_ ] = residue;

}


/// @brief
/// Prepares the graph to begin simulated annealing.
///
/// @details
/// Invokes both base-class prepare_for_simulated_annealing subroutines: InteractionGraphBase first, to prepare the
/// NPDHBondNodes and NPDHBondEdges. Then the AdditionalBackgroundNodesInteractionGraph, to prepare the NPDHBondBackgroundNodes,
/// the NPDHBondBackgroundEdges, and to do a little more preparing of the NPDHBondNodes.
/// Also computes background/background overlap.
/// This is the 2nd major entry point into the NPDHBIG.
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::prepare_for_simulated_annealing() {

	if ( prepared_for_simulated_annealing_ ) return;

	// G::prepare() calls InteractionGraphBase::prepare_for_simulated_annealing() - LinmemIG implements one but it also
	// calls the IGBase method.  The IGBase method, in turn, calls prep_for_simA() on all the FC Edges, and then all FC nodes.
	G::prepare_for_simulated_annealing();

	// parent::prepare() calls prep_for_simA() on all the BGNodes
	parent::prepare_for_simulated_annealing();

	prepared_for_simulated_annealing_ = true;

}


/// @brief
/// assigns state 0 -- the unassigned state -- to all (first class) vertices in the graph
///
/// @details
/// This is the 3rd entry point into the NPDHBIG.  It is called by the Annealer just before simulated annealing and rotamer
/// substitutions begin to init the graph to unassigned values everywhere.
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::blanket_assign_state_0() {

	for ( int ii = 1; ii <= parent::get_num_nodes(); ++ii ) {
		get_npd_hbond_node( ii )->assign_zero_state();
	}

	total_energy_current_state_assignment_ = total_energy_alternate_state_assignment_ = 0.0;
	npd_hbond_energy_current_state_assignment_ = npd_hbond_energy_alternate_state_assignment_ = 0.0;

}

///
/// @brief
/// After every 2^10 commits, the graph traverses its nodes and edges and
/// re-tallies the total energy of the current state assignment.  This update
/// prevents the accumulation of numerical drift, increasing accuracy.
///
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::update_internal_energy_totals_npd_hbond() {

	//NPDHBondInteractionGraph< V, E, G >::print_npd_hbond_avoidance_stats();

	parent::update_internal_energy_totals();

	// iterate across all foreground-  and background nodes and sum their current hbond energies
	npd_hbond_energy_current_state_assignment_ = 0;
	for ( int ii = 1; ii <= parent::get_num_nodes(); ++ii ) {
		npd_hbond_energy_current_state_assignment_ += get_npd_hbond_node( ii )->get_upper_npd_hbond_energy_totals();
	}

	for ( int ii = 1; ii <= parent::get_num_background_nodes(); ++ii ) {
		npd_hbond_energy_current_state_assignment_ += get_npd_hbond_bg_node( ii )->get_upper_npd_hbond_energy_totals();
	}

	// npd_hbond_energy_current_state_assignment_ *= npd_hbond_score_weight_;

	total_energy_current_state_assignment_ = parent::get_energy_PD_current_state_assignment() + npd_hbond_energy_current_state_assignment_;

	num_commits_since_last_update_ = 0;

}

///
/// @brief
/// Allows the sim-annealer to specify a deltaE threshold above which, it is no longer necessary to be very accurate.
///
/// @details
/// When the annealer asks the graph to consider a state substitution that produces a large collision, the graph may
/// approximate the hpatch deltaE instead of performing expensive sphere overlap computations.  The deltaE returned by
/// consider_substitution() will be inaccurate, but if the annealer is unlikely to accept the substitution, then time
/// can be saved. The graph guarantees that if the annealer does commit that substitution that it will go back and
/// perform the hpatch computations and return an accurate total energy for the graph.
///
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::set_errorfull_deltaE_threshold( core::PackerEnergy deltaE ) {

	NPDHBondInteractionGraph< V, E, G >::reset_npd_hbond_avoidance_stats();
	deltaE_threshold_for_avoiding_npd_hbond_calcs_ = deltaE;

}

/// @brief
/// reports on the level of success for hpatch score calculation procrastination
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::print_npd_hbond_avoidance_stats() {

	if ( num_state_substitutions_considered_ == 0 ) {
		return;
	}

	TR << "num state substitutions considered: " << num_state_substitutions_considered_ << ", "
		<< "num npd-hbond calcs procrastinated: " << num_npd_hbond_comps_procrastinated_ << ", "
		<< "num npd-hbond calcs later computed: " << num_npd_hbond_comps_later_made_ << std::endl;
	TR << "Percent Avoided: " << (double) (num_npd_hbond_comps_procrastinated_ - num_npd_hbond_comps_later_made_) / num_state_substitutions_considered_ << ", ";

	if ( num_npd_hbond_comps_procrastinated_ != 0 ) {
		TR << "Worthwhile Procrastination: " << (double) (num_npd_hbond_comps_procrastinated_ - num_npd_hbond_comps_later_made_) / num_npd_hbond_comps_procrastinated_ << std::endl;
	} else {
		TR << "Worthwhile Procrastination: " << "N/A" << std::endl;
	}

}

/// @brief
/// resets static member variables of NPDHBondIG that measure how worthwhile hpatch calculation procrastination is.
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::reset_npd_hbond_avoidance_stats() {
	num_state_substitutions_considered_ = 0;
	num_npd_hbond_comps_procrastinated_ = 0;
	num_npd_hbond_comps_later_made_ = 0;
}

///
/// @brief
/// Returns the (possibly approximate) change in energy induced by switching a particular node from its currently assigned state to some alternate state.
///
/// @details
/// First, queries the NPDHBondNode for the pairwise-decomposable (PD) energy. If the PD difference
/// implies a collision, then the NPDHBondIG pretends as if the state substitution causes the best
/// improvement possible in hpatch score, and returns the PD difference + pretend hpatch difference.
/// It will procrastinate computing the actual hpatch score difference until the guiding SimAnnealer
/// decides to commit the substitution. If the SimAnnealer rejects the substitution, then the work
/// to compute the hpatch score is never done. If it is unclear that the SimAnnealer will reject the
/// substitution based on the PD difference, then the Graph goes ahead and computes the change in hpatch
/// score accurately.
///
/// This function is the 4th major entry point from the Annealer into the NPDHBIG.
///
/// Also returns the sum of the two body energies for the node in its current state; the sim-annealer accepts state
/// substitutions at higher chance if the original state was also at a poor energy.
///
/// @param
/// node_ind - [in] - the index of the (first class) node
/// new_state - [in] - the alternate state that the node should consider
/// delta_energy - [out] - the change in energy induced on the entire graph by substituting a node's current state with the alternate.
///       This energy may be inaccurate if it exceeds a threshold set by the sim-annealer.
/// prev_energy_for_node - [out] - the sum of the pair-wise decomposable portion of the energy function for the node's currently assigned state
///
///
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::consider_substitution(
	int node_ind,
	int new_state,
	core::PackerEnergy & delta_energy,
	core::PackerEnergy & prev_energy_for_node
) {

	reset_from_previous_delta_npd_hbond_comp();

	++num_state_substitutions_considered_;

	node_considering_alt_state_ = node_ind;
	alt_state_being_considered_ = new_state;
	last_considered_substitution_accepted_ = false;

	// the below deltaE may be an estimate of the change in energy and not the actual value
	core::PackerEnergy deltaE = get_npd_hbond_node( node_ind )->calculate_PD_deltaE_for_substitution( new_state, prev_energy_for_node );

	calculated_npd_hbond_deltaE_ = false;

	if ( decide_procrastinate_npd_hbond_computations( deltaE, deltaE_threshold_for_avoiding_npd_hbond_calcs_ ) ) {
		++num_npd_hbond_comps_procrastinated_;
	} else {
		deltaE += calculate_npd_hbond_deltaE();
		calculated_npd_hbond_deltaE_ = true;
		deltaE_for_substitution_ = deltaE;
	}

	delta_energy = deltaE;
	total_energy_alternate_state_assignment_ = deltaE + total_energy_current_state_assignment_;
}

///
/// @details
/// Goes through the entire process of calculating the hpatch deltaE for a substitution.
///
template < typename V, typename E, typename G >
core::PackerEnergy NPDHBondInteractionGraph< V, E, G >::calculate_npd_hbond_deltaE()
{
	Real hbond_deltaE = get_npd_hbond_node( node_considering_alt_state_ )->consider_alternate_state();
	npd_hbond_energy_alternate_state_assignment_ = npd_hbond_energy_current_state_assignment_ + hbond_deltaE;
	return hbond_deltaE;
}


/// @brief
/// Commits the substitution that the sim annealer had previously asked the graph to consider.  Returns the accurate total energy for the graph.
template < typename V, typename E, typename G >
core::PackerEnergy NPDHBondInteractionGraph< V, E, G >::commit_considered_substitution()
{

	last_considered_substitution_accepted_ = true;

	core::PackerEnergy npd_hbond_deltaE = 0.0;
	if ( ! calculated_npd_hbond_deltaE_ ) {
		npd_hbond_deltaE = calculate_npd_hbond_deltaE(); // updates all the Nodes and calculates the hpatch score delta

		// Calls to get_pd_energy_delta() must proceed the call to commit_considered_substitution
		deltaE_for_substitution_ = get_npd_hbond_node( node_considering_alt_state_ )->get_pd_energy_delta() + npd_hbond_deltaE;
		++num_npd_hbond_comps_later_made_;
	}

	get_npd_hbond_node( node_considering_alt_state_ )->commit_considered_substitution();

	total_energy_current_state_assignment_ = total_energy_current_state_assignment_ + deltaE_for_substitution_;
	npd_hbond_energy_current_state_assignment_ = npd_hbond_energy_alternate_state_assignment_;

	node_considering_alt_state_ = -1;
	++num_commits_since_last_update_;

	if ( num_commits_since_last_update_ == COMMIT_LIMIT_BETWEEN_UPDATES ) {
		update_internal_energy_totals_npd_hbond();
	}

	track_npd_hbond_E_min();

	return total_energy_current_state_assignment_;
}

/// @brief
/// Switch the state assignment of every first class node in the graph.
/// Useful, for instance, if you want to switch to the best network state that you've found so far.
///
/// @details
/// This function is the last major entry point from the Annealer into the NPDHBIG.
template < typename V, typename E, typename G >
core::PackerEnergy NPDHBondInteractionGraph< V, E, G >::set_network_state( ObjexxFCL::FArray1_int & node_states ) {

	reset_from_previous_delta_npd_hbond_comp();

	for ( int ii = 1; ii <= parent::get_num_nodes(); ++ii ) {
		core::PackerEnergy deltaE = 0.0;
		core::PackerEnergy previousE = 0.0;
		consider_substitution( ii, node_states( ii ), deltaE, previousE ); // might get procrastinated but that's ok
		commit_considered_substitution(); // any procrastinated calc will get updated correctly here
	}

	update_internal_energy_totals_npd_hbond();
	return total_energy_current_state_assignment_;

}

/// @brief
/// returns the energy of the entire graph under the current network state assignment.  Also sends a bunch of information to standard error.
/// Only seems to be called by the MultiCoolAnnealer.
template < typename V, typename E, typename G >
core::PackerEnergy NPDHBondInteractionGraph< V, E, G >::get_energy_current_state_assignment() {
	return total_energy_current_state_assignment_;
}

///
/// @brief
/// Should return a measurement of the memory used by the interaction graph
/// to store the rotamer pair energies.  Unimplemented.
///
template < typename V, typename E, typename G >
int NPDHBondInteractionGraph< V, E, G >::get_edge_memory_usage() const {
	return 0;
}

///
template < typename V, typename E, typename G >
unsigned int NPDHBondInteractionGraph< V, E, G >::count_static_memory() const {
	return sizeof ( NPDHBondInteractionGraph< V, E, G > );
}

///
template < typename V, typename E, typename G >
unsigned int NPDHBondInteractionGraph< V, E, G >::count_dynamic_memory() const {
	unsigned int total_memory = parent::count_dynamic_memory();

	total_memory += resid_2_bgenumeration_.size() * sizeof( Size );
	total_memory += bgenumeration_2_resid_.size() * sizeof( Size );

	return total_memory;
}

/// @brief
/// returns the sum of the PD energy and the hpatch energy for all members first class members of a user-defined
/// vertex subset.  Unimplemented.
template < typename V, typename E, typename G >
Real NPDHBondInteractionGraph< V, E, G >::get_energy_sum_for_vertex_group( Size ) {
	//apl functionality stubbed out for now
	return 0;
}

template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::print_internal_energies_for_current_state_assignment() {

	// print out the one-body and hpatch energies for all first class nodes
	TR << "internal energies: " << std::endl;
	for ( int ii = 1; ii <= parent::get_num_nodes(); ++ii ) {
		Real one_body = get_npd_hbond_node( ii )->get_curr_state_one_body_energy();
		TR << "node " << ii << " 1b: " << one_body;
		Real hbE = get_npd_hbond_node( ii )->get_upper_npd_hbond_energy_totals();
		TR << ", upper hbE = " << hbE;

		if ( ii % 3 == 0 ) {
			TR << std::endl;
		}
	}
	TR << std::endl;

	// print out the hpatch energies for all background nodes
	for ( int ii = 1; ii <= parent::get_num_background_nodes(); ++ii ) {
		Real hbE = get_npd_hbond_bg_node( ii )->get_upper_npd_hbond_energy_totals();
		TR << "bg res: " << bgenumeration_2_resid_[ ii ] << " upper hbE: " << hbE << std::endl;
	}

	// print out the two-body energies for all edges between first-class nodes only?
	int count_edges = 0;
	for ( auto iter = parent::get_edge_list_begin(); iter != parent::get_edge_list_end(); ++iter ) {
		Real edge_energy = (static_cast< E *> (*iter))->get_curr_pd_energy_total();
		TR << "edge: " << edge_energy << " ";

		if ( count_edges % 5 == 0 ) {
			TR << std::endl;
		}
		++count_edges;
	}
}

/// @brief
/// useful for debugging
template< typename V, typename E, typename G >
void
NPDHBondInteractionGraph<V, E, G>::print() const {

	std::cout << "NPDHBond Interaction Graph state: " << std::endl;
	std::cout << "nodes: " << std::endl;
	for ( int jj = 1; jj <= parent::get_num_nodes(); ++jj ) {
		get_npd_hbond_node( jj )->print();
	}

	std::cout << "bgnodes: " << std::endl;
	for ( int ii = 1; ii <= parent::get_num_background_nodes(); ++ii ) {
		get_npd_hbond_bg_node( ii )->print();
	}
}


///
/// @brief
/// Returns the state on each FCNode, but not necessarily in pose resid order. Only used by the unit tests.
///
template< typename V, typename E, typename G >
std::vector< int > NPDHBondInteractionGraph<V, E, G>::get_network_state() const {

	std::vector< int > networkstate;
	for ( int jj = 1; jj <= parent::get_num_nodes(); ++jj ) {
		networkstate.push_back( get_npd_hbond_node(jj)->get_current_state() );
	}
	return networkstate;
}

/// @brief
/// Sets the observed_sufficient_npd_hbond_E_to_predict_min_ to true. Only used by the unit tests.
template< typename V, typename E, typename G >
void NPDHBondInteractionGraph<V, E, G>::set_observed_sufficient_boolean_true() {
	observed_sufficient_npd_hbond_E_to_predict_min_ = true;
}

/// @brief
/// Provides read access to the bg to resid array. Returns -1 if the index is not in bounds.
template < typename V, typename E, typename G >
int NPDHBondInteractionGraph< V, E, G >::bg_node_2_resid( Size node_index ) {

	if ( node_index > num_residues_assigned_as_background_ ) {
		utility_exit_with_message( "Out of bounds array index passed to bg_node_2_resid. Quitting." );
	}
	return bgenumeration_2_resid_[ node_index ];
}


template < typename V, typename E, typename G >
scoring::hbonds::HBondDatabase const & NPDHBondInteractionGraph< V, E, G >::hbond_database() const
{
	return *hbond_database_;
}

template < typename V, typename E, typename G >
scoring::hbonds::HBondOptions const & NPDHBondInteractionGraph< V, E, G >::hbond_options() const
{
	return *hbond_options_;
}

template < typename V, typename E, typename G >
scoring::hbonds::NPDHBondSet const & NPDHBondInteractionGraph< V, E, G >::npd_hbond_set() const
{
	return *npd_hbond_set_;
}

template < typename V, typename E, typename G >
utility::vector1< char > & NPDHBondInteractionGraph< V, E, G >::hbonding_to_res_vector()
{
	return hbonding_to_res_;
}

template < typename V, typename E, typename G >
NPDHBondOP NPDHBondInteractionGraph< V, E, G >::unused_hbond()
{
	// return NPDHBondOP( new NPDHBond ); // TEMP!

	if ( hbonds_queue_.empty() ) {
		return NPDHBondOP( new NPDHBond );
	}
	NPDHBondOP next = hbonds_queue_.front();
	hbonds_queue_.pop_front();
	return next;
}

template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::return_hbond_to_queue( NPDHBondOP const & hbond )
{
	// mark the hbond as not part of any residue; this is mostly for debugging purposes:
	// if an hbond has been recycled but is still being used in any calculation, then something
	// within the NPDHBondInteractionGraph has gone wrong. There are assertions in the code to
	// make sure that the hbond that's about to be used belongs to the residue that's using it
	// so marking the residue as 0 will trip one of those assertions
	hbond->don_rsd_ = hbond->acc_rsd_ = 0;

	hbonds_queue_.push_back( hbond );
}


template < typename V, typename E, typename G >
Real NPDHBondInteractionGraph< V, E, G >::npd_hb_weight(
	scoring::hbonds::HBEvalType eval_type,
	bool intra_res
)
{
	return scoring::hbonds::npd_hb_eval_type_weight( eval_type, scorefxn_->weights(), intra_res );
}



/// @brief
/// factory method pattern for instantiation of NPDHBondNode objects, used by InteractionGraphBase class.
template < typename V, typename E, typename G >
NodeBase* NPDHBondInteractionGraph< V, E, G >::create_new_node( int node_index, int num_states ) {
	return new NPDHBondNode< V, E, G >( this, node_index, num_states );
}

/// @brief
/// factory method pattern for instantiation of NPDHBondEdge objects, used by InteractionGraphBase class.
template < typename V, typename E, typename G >
EdgeBase* NPDHBondInteractionGraph< V, E, G >::create_new_edge( int index1, int index2 ) {
	return new NPDHBondEdge< V, E, G >( this, index1, index2 );
}

/// @brief
/// factory method pattern for instantiation of NPDHBondBackgroundNode objects, used by AdditionalBackgroundNodesInteractionGraph class.
template < typename V, typename E, typename G >
BackgroundNode< V, E, G >* NPDHBondInteractionGraph< V, E, G >::create_background_node( int node_index ) {
	return new NPDHBondBackgroundNode< V, E, G >( this, node_index );
}

/// @brief
/// factory method pattern for instantiation of NPDHBondBackgroundEdge objects, used by AdditionalBackgroundNodesInteractionGraph class.
template < typename V, typename E, typename G >
BackgroundToFirstClassEdge< V, E, G >* NPDHBondInteractionGraph< V, E, G >::create_background_edge( int fc_node_index, int bg_node_index ) {
	return new NPDHBondBackgroundEdge< V, E, G >( this, fc_node_index, bg_node_index );
}

///
/// @brief
/// Keeps track of the minimum hpatch score seen.  Every 100 substitutions, updates the variable npd_hbond_score_min_last_100.
///
template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::track_npd_hbond_E_min() {

	++num_substitutions_since_npd_hbond_min_update_;

	Real alt_hpatchE = npd_hbond_energy_current_state_assignment_;

	if ( npd_hbond_score_min_recent_ > alt_hpatchE ) {
		npd_hbond_score_min_recent_ = alt_hpatchE;
	}

	if ( num_substitutions_since_npd_hbond_min_update_ == 100 ) { // only update min every 100 calls to track_hpatchE_min (aka commits)
		npd_hbond_score_min_last_100_ = npd_hbond_score_min_recent_;
		if ( npd_hbond_energy_current_state_assignment_ < npd_hbond_score_min_last_100_ ) {
			npd_hbond_score_min_last_100_ = npd_hbond_energy_current_state_assignment_;
		}
		observed_sufficient_npd_hbond_E_to_predict_min_ = true;
		num_substitutions_since_npd_hbond_min_update_ = 0;
	}

}

template < typename V, typename E, typename G >
void NPDHBondInteractionGraph< V, E, G >::reset_from_previous_delta_npd_hbond_comp() {
	// TO DO: reclaim the new hydrogen bonds for the previously-substituted residue
	// if the last considered rotamer substitution was rejected.
	if ( last_considered_substitution_accepted_ ) return;
	if ( ! calculated_npd_hbond_deltaE_ ) return;

	// Reclaim the alternate hydrogen bond objects that will not be used in the future
	// (until we reassign them as different hydrogen bonds) because the rotamer that they
	// correspond to has not been accepted
	get_npd_hbond_node( node_considering_alt_state_ )->reset_after_rejected_substitution();

}


/// @details
/// Makes the decision whether or not to procrastinate calculating the NPD-hbond score. Basically, if the PD energy got better (dE < 0)
/// then return false so we don't procrastinate the calculation (because the alternate state probably will be accepted?). If the best
/// guess for the hpatch deltaE also comes back better (dE < 0), then return false.  Finally, if the difference between the deltaE
/// for the PD terms and the (guessed) npd-hbond deltaE is greater than the threshold, return true so we do procrastinate. So basically
/// if the PD energy gets much worse, procrastinate.  Otherwise, don't.
template < typename V, typename E, typename G >
bool NPDHBondInteractionGraph< V, E, G >::decide_procrastinate_npd_hbond_computations(
	Real const pd_deltaE,
	Real const threshold
) const
{

	return false; // TEMP! Short circuit this calculation for now and always compute the hbond deltaE

	Real npd_hbond_deltaE_max = 0;

	if ( ! observed_sufficient_npd_hbond_E_to_predict_min_ ) {
		return false;
	}

	if ( threshold < 0 || pd_deltaE < 0 ) {
		return false;
	}

	npd_hbond_deltaE_max += npd_hbond_energy_current_state_assignment_ - npd_hbond_score_min_last_100_;

	// pd_deltaE must be positive, npd_hbond_deltaE must also be positive.
	// threshold of 5 means, if the PD got more than 5 worse than the guessed hpatch deltaE (e.g. 10 - 3(?) = 7 > 5), procrastinate
	if ( (pd_deltaE - npd_hbond_deltaE_max) > threshold ) {
		return true;
	}
	return false;

}

} //end namespace
} //end namespace
} //end namespace

#endif
