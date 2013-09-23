// -*- Mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/HPatchInteractionGraph.hh
/// @brief  Interaction graph which implements a non-PD score that optimizes against surface hydrophobic patches.
/// @author Ron Jacak (ron.jacak@gmail.com)
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_interaction_graph_HPatchInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_HPatchInteractionGraph_hh

//Rosetta Headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.hh>

#include <core/graph/DisjointSets.hh>

#include <core/pack/interaction_graph/HPatchInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/AdditionalBackgroundNodesInteractionGraph.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>
#include <core/pack/interaction_graph/RotamerDots.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/sasa.hh>
#include <core/pack/interaction_graph/SurfacePotential.hh>

#include <basic/Tracer.hh>

//Utility Headers
// AUTO-REMOVED #include <numeric/xyzVector.io.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>  // needed to get arg_max - DO NOT AUTO-REMOVE!
#include <utility/exit.hh>
#include <utility/string_util.hh> // needed to get trim - DO NOT AUTOREMOVE!

//ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1.io.hh> // needed to stream operator of FArray1Dint, line 4738 - DO NOT AUTOREMOVE!
#include <ObjexxFCL/format.hh> // needed for I() - DO NOT AUTOREMOVE!

//C++ Headers
#include <vector>

//Auto Headers
//#define DOUBLE_CHECK_COUNTS 1
//#define FILE_DEBUG 1

namespace core {
namespace pack {
namespace interaction_graph {

struct
exposed_hydrophobic_data
{
	exposed_hydrophobic_data( Size fc_bg, Size node_ind, Size exhphobe_ind ) :
		fc_bg_( fc_bg ),
		node_index_( node_ind ),
		exhphobe_index_( exhphobe_ind )
	{}

	Size fc_bg_; // 1 for first-class node; 2 for background node
	Size node_index_; // the node index
	Size exhphobe_index_; // the exposed hydrophobic index on that node
};

/// Tracer instance for this file
static basic::Tracer TR_NODE("core.pack.hpatchig.node");
static basic::Tracer TR_EDGE("core.pack.hpatchig.edge");
static basic::Tracer TR_BGNODE("core.pack.hpatchig.bgnode");
static basic::Tracer TR_BGEDGE("core.pack.hpatchig.bgedge");
static basic::Tracer TR_HIG("core.pack.hpatchig.ig");
static basic::Tracer TR_STATS("core.pack.hpatchig.stats");


template < typename V, typename E, typename G > class HPatchNode;
template < typename V, typename E, typename G > class HPatchBackgroundNode;
template < typename V, typename E, typename G > class HPatchEdge;
template < typename V, typename E, typename G > class HPatchBackgroundEdge;
template < typename V, typename E, typename G > class HPatchInteractionGraph;


//----------------------------------------------------------------------------//
//---------------------------- HPatch Node Class -------------------------------//
//----------------------------------------------------------------------------//

///
/// @begin HPatchNode
///
/// @brief
/// Defines a FirstClass node which will keep track of changes in the SASA and hpatch score.
/// FirstClassNode is defined and implemented in AdditionalBackgroundNodesInteractionGraph.
///
/// @remarks
/// No public default constructor makes this class uncopyable.
///
template < typename V, typename E, typename G >
class HPatchNode : public FirstClassNode< V, E, G > {

	public:
		typedef FirstClassNode< V, E, G > parent;
		typedef V grandparent;

	public:
		HPatchNode( G* owner, int node_index, int num_states );
		virtual ~HPatchNode();

		// setter for the rotamers object. called at the very beginning of the HIG::initialize() method
		void set_rotamers( rotamer_set::RotamerSetCOP rotamers );
		conformation::ResidueCOP get_rotamer( int state ) const;
		conformation::ResidueCOP curr_state_rotamer() const {
			return rotamers_vector_[ parent::get_current_state() ];
		}

		conformation::ResidueCOP alt_state_rotamer() const {
			return rotamers_vector_[ parent::get_alternate_state() ];
		}

		virtual void set_rotamer_dots_for_state( Size state, RotamerDots const & rd );

		bool overlaps( HPatchNode< V, E, G >* neighbor );
		bool detect_any_overlap_with_rotamer( RotamerDots const & rotamer ) const;

		virtual void prepare_for_simulated_annealing();
		void initialize_overlap_with_background(
			RotamerDots const & bg_rotamer,
			std::vector< RotamerDotsCache > & overlap,
			std::vector< utility::vector1< utility::vector1< bool > > > & states_atom_atom_overlap_on_bg_res
		);

		// for HIG entry point blanket_assign_state_0()
		virtual void assign_zero_state();
		void acknowledge_neighbors_substitution();

		core::PackerEnergy calculate_PD_deltaE_for_substitution( int alternate_state, core::PackerEnergy & prev_PDenergies_for_node );
		core::PackerEnergy get_pd_energy_delta();
		Real consider_alternate_state();
		Real get_current_state_sasa() const;
		Real get_current_state_sasa( Size atom_index ) const;
		Real get_alternate_state_sasa( Size atom_index ) const;

		utility::vector1< utility::vector1< bool > > const & get_atom_atom_self_overlaps_for_state( Size state ) const;

		inline Size get_current_state_num_atoms() const { return current_state_rotamer_dots_.get_num_atoms(); }
		inline Size get_alt_state_num_atoms() const { return alt_state_rotamer_dots_.get_num_atoms(); }

		inline
		int wt_seqpos_for_node() const {
			return get_hpatch_owner()->rotamer_sets().moltenres_2_resid( parent::get_node_index() );
		}

		inline
		conformation::Residue const & wt_residue_for_node() const {
			return get_hpatch_owner()->pose().residue( (get_hpatch_owner()->rotamer_sets().moltenres_2_resid( parent::get_node_index() )) );
		}

		Real update_state_for_neighbors_substitution( HPatchNode<V,E,G>* node_considering_substitution,
			RotamerDots & neighbors_alternate_state,
			RotamerDotsCache const & neighbors_curr_state_overlap_with_this,
			RotamerDotsCache & this_overlap_with_neighbors_alternate,
			RotamerDotsCache & neighbors_alternate_overlap_with_this,
			utility::vector1< utility::vector1< bool > > & alt_state_atom_atom_overlaps_cache
		);

		void reset_alt_state_dots();

		void commit_considered_substitution();

		virtual unsigned int getMemoryUsageInBytes() const;
		virtual unsigned int count_static_memory() const;
		virtual unsigned int count_dynamic_memory() const;

		//void write_dot_kinemage( std::ofstream & output_kin );
		virtual void print() const;

		// Extra methods only used only for the unit tests.
		RotamerDots const & get_current_state_rotamer_dots();
		RotamerDots const & get_alt_state_rotamer_dots();

		InvRotamerDots const & curr_state_inv_dots() const {
			assert( curr_state_inv_dots_.rotamer() == current_state_rotamer_dots_.rotamer() );
			return curr_state_inv_dots_;
		}
		InvRotamerDots const & alt_state_inv_dots() const {
			assert( alt_state_inv_dots_.rotamer() == alt_state_rotamer_dots_.rotamer() );
			return alt_state_inv_dots_;
		}

		Size
		max_hphobe_atoms_any_restype() const {
			return max_hphobes_;
		}

		utility::vector1< Size > const &
		curr_state_hphobes() const {
			return hphobe_ats_for_restype_group_[ restype_group_for_rotamers_[ grandparent::get_current_state() ]];
		}

		utility::vector1< Size > const &
		alt_state_hphobes() const {
			return hphobe_ats_for_restype_group_[ restype_group_for_rotamers_[ grandparent::get_alternate_state() ]];
		}

		Size n_curr_state_hphobes() const { return curr_state_hphobes().size(); }
		Size n_alt_state_hphobes() const { return alt_state_hphobes().size(); }

		utility::vector1< Size > const & curr_state_exp_hphobes() const { return curr_state_exp_hphobes_; }
		utility::vector1< Size > const & alt_state_exp_hphobes() const { return alt_state_exp_hphobes_; }
		Size n_curr_state_exp_hphobes() const { return curr_state_exp_hphobes_.size(); }
		Size n_alt_state_exp_hphobes() const { return alt_state_exp_hphobes_.size(); }

	protected:
		inline
		HPatchEdge< V, E, G >* get_incident_hpatch_edge( int index ) const {
			return (HPatchEdge< V, E, G > *) parent::get_incident_edge( index );
		}

		inline
		HPatchBackgroundEdge< V, E, G >* get_edge_to_hpatch_bg_node( int index ) const {
			return (HPatchBackgroundEdge< V, E, G > *) parent::get_edge_to_bg_node( index );
		}

		inline
		HPatchInteractionGraph< V, E, G >* get_hpatch_owner() const {
			return (HPatchInteractionGraph< V, E, G > *) parent::get_owner();
		}

		Real get_sasa_difference() const;

	private:
		void update_alt_state_exphphobes();
		void initialize_self_overlap();
		void initialize_atom_atom_overlap_cache();

		// no default constructor, uncopyable
		HPatchNode();
		HPatchNode( HPatchNode< V, E, G > const & );
		HPatchNode< V, E, G > & operator = ( HPatchNode< V, E, G > const & );

	private:
		utility::vector1< conformation::ResidueCOP > rotamers_vector_;
		utility::vector1< Size > restype_group_for_rotamers_;
		utility::vector1< utility::vector1< Size > > hphobe_ats_for_restype_group_;
		Size max_hphobes_;

		// outer vector is all states possible at this Node, inner two vectors are ii and jj atoms
		utility::vector1< utility::vector1< utility::vector1< bool > > > self_atom_atom_overlaps_;

		// vector1 of RotamerDots objects. for every state on every molten residue there is a RotamerDots object.
		utility::vector1< RotamerDots > self_and_bg_dots_for_states_;

		RotamerDots current_state_rotamer_dots_;
		RotamerDots alt_state_rotamer_dots_;

		InvRotamerDots curr_state_inv_dots_;
		InvRotamerDots alt_state_inv_dots_;

		utility::vector1< Size > curr_state_exp_hphobes_;
		utility::vector1< Size > alt_state_exp_hphobes_;

		bool alt_state_dots_matches_current_state_dots_;

};



//----------------------------------------------------------------------------//
//------------------------- HPatch Background Node Class -----------------------//
//----------------------------------------------------------------------------//

///
/// @begin HPatchBackgroundNode
///
/// @brief
/// Defines a Background Node which will contribute to changes in SASA/hpatchE due to state changes on neighboring nodes,
/// and not because of state changes to it.
///	No default constructor makes this class uncopyable
///
template < typename V, typename E, typename G >
class HPatchBackgroundNode : public BackgroundNode< V, E, G > {

	public:
		typedef BackgroundNode< V, E, G > parent;

	public:
		HPatchBackgroundNode( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int node_index );
		virtual ~HPatchBackgroundNode();

		void set_rotamer( conformation::ResidueOP const & rotamer );
		conformation::ResidueCOP get_rotamer() const;

		void set_rotamer_dots( RotamerDots const & bg_rd );

		bool detect_overlap( HPatchNode< V, E, G >* node ) const;
		virtual void prepare_for_simulated_annealing();

		void initialize_bg_bg_overlap( HPatchBackgroundNode< V, E, G > & other );

		Real update_state_for_substitution(
			HPatchNode< V, E, G >* fc_node_changing,
			RotamerDotsCache const & nodes_curr_overlap_with_bg_res,
			RotamerDotsCache const & nodes_alt_overlap_with_bg_res
		);

		void reset_alt_state_dots();
		void acknowledge_substitution();
		Real get_current_sasa() const;
		Real get_current_sasa( Size atom_index ) const;
		Real get_alternate_sasa() const;
		Real get_alternate_sasa( Size atom_index ) const;

		utility::vector1< utility::vector1< bool > > const & get_atom_atom_self_overlaps() const;
		utility::vector1< Size > const & get_hphobes() const { return hphobe_ats_; }
		Size n_hphobes() const { return n_hphobes_; }

		virtual unsigned int count_static_memory() const;
		virtual unsigned int count_dynamic_memory() const;

		//void write_dot_kinemage( std::ofstream & output_kin );
		virtual void print() const;

		// Only used for the unit tests.
		RotamerDots const & get_current_state_rotamer_dots();
		RotamerDots const & get_alt_state_rotamer_dots();

		InvRotamerDots const & curr_state_inv_dots() const {
			assert( curr_state_inv_dots_.rotamer() == alt_state_rotamer_dots_.rotamer() );
			return curr_state_inv_dots_;
		}

		InvRotamerDots const & alt_state_inv_dots() const {
			assert( alt_state_inv_dots_.rotamer() == alt_state_rotamer_dots_.rotamer() );
			return alt_state_inv_dots_;
		}

		utility::vector1< Size > const & curr_state_exp_hphobes() const { return curr_state_exp_hphobes_; }
		utility::vector1< Size > const & alt_state_exp_hphobes() const { return alt_state_exp_hphobes_; }

		Size n_curr_state_exp_hphobes() const { return curr_state_exp_hphobes_.size(); }
		Size n_alt_state_exp_hphobes() const { return alt_state_exp_hphobes_.size(); }

	protected:
		inline
		HPatchBackgroundEdge< V, E, G > * get_hpatch_bg_edge( int index ) {
			return (HPatchBackgroundEdge< V, E, G > *) parent::get_incident_edge( index );
		}

		inline
		HPatchInteractionGraph< V, E, G >* get_hpatch_owner() const {
			return (HPatchInteractionGraph< V, E, G > *) parent::get_owner();
		}

	private:
		void update_alt_state_exphphobes();
		void initialize_self_overlap();
		void initialize_atom_atom_overlaps();

		conformation::ResidueOP rotamer_;
		utility::vector1< Size > hphobe_ats_;
		Size n_hphobes_;

		bool prepared_for_simA_;

		RotamerDots current_state_rotamer_dots_;
		RotamerDots alt_state_rotamer_dots_;

		InvRotamerDots curr_state_inv_dots_;
		InvRotamerDots alt_state_inv_dots_;

		utility::vector1< Size > curr_state_exp_hphobes_;
		utility::vector1< Size > alt_state_exp_hphobes_;

		bool alt_state_dots_matches_current_state_dots_;

		utility::vector1< utility::vector1< bool > > self_atom_atom_overlaps_;

		HPatchBackgroundNode();
		HPatchBackgroundNode( HPatchBackgroundNode< V, E, G > const & );
		HPatchBackgroundNode< V, E, G > & operator= ( HPatchBackgroundNode< V, E, G > const & );

};


//----------------------------------------------------------------------------//
//------------------------------ HPatch Edge Class -----------------------------//
//----------------------------------------------------------------------------//

///
/// @begin HPatchEdge
///
/// @brief
/// Defines a HPatch Edge which connects two first-class HPatch Nodes. Edges have to keep some state so that updates
/// to SASA and the hpatch score can be done fast.
///
template < typename V, typename E, typename G >
class  HPatchEdge : public FirstClassEdge< V, E, G > {

	public:
		typedef  FirstClassEdge< V, E, G >  parent;

	public:
		HPatchEdge( G * owner, int node1, int node2 );
		virtual ~HPatchEdge();

		virtual void prepare_for_simulated_annealing();

		void acknowledge_state_zeroed( int node_index );

		Real update_state_at_neighbor(
			int node_considering_substitution,
			int alt_state,
			RotamerDots & alt_state_dots
		);

		void acknowledge_substitution();

		utility::vector1< utility::vector1< bool > > const & get_current_state_atom_atom_overlaps() const;
		utility::vector1< utility::vector1< bool > > const & get_alt_state_atom_atom_overlaps() const;

		// Virtual methods from EdgeBase
		virtual void declare_energies_final();

		virtual unsigned int getMemoryUsageInBytes() const;
		virtual unsigned int count_static_memory() const;
		virtual unsigned int count_dynamic_memory() const;

		Real get_current_two_body_energy() const;

	protected:
		inline
		HPatchNode< V, E, G >* get_hpatch_node( int index ) {
			return (HPatchNode< V, E, G >*) E::get_node( index );
		}

	private:
		void inform_non_changing_node_of_neighbors_change();

		//no default constructor, uncopyable
		HPatchEdge();
		HPatchEdge( HPatchEdge< V, E, G > const & );
		HPatchEdge< V, E, G > & operator = ( HPatchEdge< V, E, G > const & );

	private:
		int node_changing_;
		int node_not_changing_;
		int nodes_curr_states_[2];
		int nodes_alt_states_[2];

		RotamerDotsCache nodes_curr_pair_dot_counts_[2];
		RotamerDotsCache nodes_alt_pair_dot_counts_[2];

		utility::vector1< utility::vector1< bool > > current_state_atom_atom_overlaps_;
		utility::vector1< utility::vector1< bool > > alt_state_atom_atom_overlaps_;

		/// pairs of hphobes that have exposed overlap
		utility::vector1< std::pair< Size, Size > > curr_state_exolap_hphobes_;
		utility::vector1< std::pair< Size, Size > > alt_state_exolap_hphobes_;

};



//----------------------------------------------------------------------------//
//------------------- HPatch Background Edge Class -----------------------//
//----------------------------------------------------------------------------//

///
/// @begin HPatchBackgroundEdge
///
/// @brief
/// Defines an edge between a FirstClass (HPatchNode) and a background node (HPatchBackgroundNode)
///
/// @details
/// In addition to implementing the virtual base class methods, this class additionally defines methods
/// relating to keeping track of data relating to SASA/hpatch.
///
template < typename V, typename E, typename G >
class HPatchBackgroundEdge : public BackgroundToFirstClassEdge< V, E, G > {

	public:
		typedef  BackgroundToFirstClassEdge< V, E, G >  parent;

	public:
		HPatchBackgroundEdge( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int first_class_node_index, int background_node_index );
		virtual ~HPatchBackgroundEdge();

		void prepare_for_simulated_annealing();
		void initialize_overlap_cache( RotamerDots const & bg_residue );

		void acknowledge_state_change( int new_state );

		Real update_state_at_neighbor( int alt_state );
		void acknowledge_substitution();

		utility::vector1< utility::vector1< bool > > const & get_atom_atom_overlaps_for_state( Size state ) const;

		virtual unsigned int count_static_memory() const;
		virtual unsigned int count_dynamic_memory() const;

	protected:
		inline
		HPatchNode< V, E, G >* get_hpatch_node() const {
			return ( HPatchNode< V, E, G >* ) parent::get_first_class_node();
		}

		inline
		HPatchBackgroundNode< V, E, G >* get_hpatch_bg_node() const {
			return ( HPatchBackgroundNode< V, E, G >* ) parent::get_background_node();
		}

	private:
		//no default constructor, uncopyable
		HPatchBackgroundEdge();
		HPatchBackgroundEdge( HPatchBackgroundEdge< V, E, G > const & );
		HPatchBackgroundEdge< V, E, G > & operator= ( HPatchBackgroundEdge< V, E, G > const & );

	private:
		bool prepared_for_simA_;
		Size bg_res_num_atoms_;

		//uses index 0; do not change to utility::vector1
		std::vector< RotamerDotsCache > node_states_coverage_of_bg_res_;

		int nodes_curr_state_;
		int nodes_alt_state_;

		RotamerDotsCache curr_dots_cache_;
		RotamerDotsCache alt_dots_cache_;

		// uses index 0; do not change
		std::vector< utility::vector1< utility::vector1< bool > > > node_states_overlap_with_bg_res_;

		/// pairs of hphobes that have exposed overlap
		utility::vector1< std::pair< Size, Size > > curr_state_exolap_hphobes_;
		utility::vector1< std::pair< Size, Size > > alt_state_exolap_hphobes_;

};



//----------------------------------------------------------------------------//
//--------------------- HPatch Interaction Graph -------------------------//
//----------------------------------------------------------------------------//

///
/// @begin HPatchInteractionGraph
///
/// @brief
/// Defines the interaction graph that will keep track of changes to the hpatch score.
///
/// @details
/// In addition to implementing the virtual base class methods, this class additionally defines methods
/// relating to keeping track of data relating to hpatch.
///
template < typename V, typename E, typename G >
class HPatchInteractionGraph : public AdditionalBackgroundNodesInteractionGraph< V, E, G > {

	public:
		typedef  AdditionalBackgroundNodesInteractionGraph< V, E, G >  parent;

	public:
		HPatchInteractionGraph( int num_nodes );
		virtual ~HPatchInteractionGraph();

		pose::Pose const & pose() const { return *pose_; }
		void set_pose( pose::Pose const & pose );

		task::PackerTask const & packer_task() const { return *packer_task_; }
		void set_packer_task( task::PackerTask const & task );

		rotamer_set::RotamerSets const & rotamer_sets() const { return *rotamer_sets_; }
		void set_rotamer_sets( rotamer_set::RotamerSets const & rotsets );

		void set_score_weight( Real weight ) { hpatch_score_weight_ = weight; }

		void initialize( rotamer_set::RotamerSetsBase const & rot_sets );

		void set_num_residues_in_protein( Size num_res );
		void set_num_background_residues( Size num_background_residues );
		void set_residue_as_background_residue( int residue );
		void set_background_residue_rotamer_dots( Size residue, conformation::Residue const & rotamer );
		void set_rotamer_dots_for_node_state( Size node_index, Size state, conformation::Residue const & rotamer );

		virtual void prepare_for_simulated_annealing();

		virtual void blanket_assign_state_0();
		Real get_hpatch_score();
		virtual void set_errorfull_deltaE_threshold( core::PackerEnergy deltaE );

		static void print_hpatch_avoidance_stats();
		static void reset_hpatch_avoidance_stats();

		virtual void consider_substitution( int node_ind, int new_state, core::PackerEnergy & delta_energy, core::PackerEnergy & prev_energy_for_node );
		core::PackerEnergy calculate_hpatch_deltaE();

		void register_fc_node_in_state0();
		void register_fc_node_affected_by_rotsub( int fc_node_ind );
		void register_bg_node_affected_by_rotsub( int bg_node_ind );

		// used only for intra-residue connections
		void update_disjoint_sets_using_cache(
			conformation::Residue const & rsd,
			InvRotamerDots const & invdots,
			utility::vector1< Size > const & exp_hphobes,
			Size residue_djs_offset,
			utility::vector1< utility::vector1< bool > > const & atom_atom_self_overlap,
			graph::DisjointSets & ds
		);

		// used for inter-residue connections
		void update_disjoint_sets_using_cache(
			conformation::Residue const & rsd1,
			InvRotamerDots const & invdots1,
			utility::vector1< Size > const & exp_hphobes1,
			Size djs_offset_1,
			conformation::Residue const & rsd2,
			InvRotamerDots const & invdots2,
			utility::vector1< Size > const & exp_hphobes2,
			Size djs_offset_2,
			utility::vector1< utility::vector1< bool > > const & atom_atom_overlaps,
			graph::DisjointSets & ds
		);

		Real calculate_alt_state_hpatch_score();
		//void blanket_reset_alt_state_dots();

		virtual core::PackerEnergy commit_considered_substitution();

		virtual core::PackerEnergy set_network_state( FArray1_int & node_states );

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
		std::vector<int> get_network_state() const;
		void set_observed_sufficient_boolean_true();
		void get_all_sasas( utility::vector1< Real > & node_sasas, utility::vector1< Real > & bgnode_sasas );

		/*utility::vector1< Size > const & get_fc_nodes_near_rotsub();
		utility::vector1< bool > const & get_fc_nodes_near_rotsub_bool();
		utility::vector1< Size > const & get_bg_nodes_near_rotsub();
		utility::vector1< bool > const & get_bg_nodes_near_rotsub_bool();

		utility::vector1< Size > const & get_fc_n_exp_hphobes();
		utility::vector1< Size > const & get_bg_n_exp_hphobes();

		utility::vector1< Size > const & get_fc_exp_hphobe_djs_offsets();
		utility::vector1< Size > const & get_bg_exp_hphobe_djs_offsets();*/

		// need to make these method public for the unit tests
		inline
		HPatchNode< V, E, G >* get_hpatch_node( int index ) const {
			return (HPatchNode< V, E, G > *) G::get_node( index );
		}

		inline
		HPatchBackgroundNode< V, E, G >* get_hpatch_bg_node( int index ) const {
			return (HPatchBackgroundNode< V, E, G > *) parent::get_background_node( index );
		}

	protected:
		virtual NodeBase* create_new_node( int node_index, int num_states);
		virtual EdgeBase* create_new_edge( int index1, int index2);
		virtual BackgroundNode< V, E, G >* create_background_node( int node_index );
		virtual BackgroundToFirstClassEdge< V, E, G >* create_background_edge( int fc_node_index, int bg_node_index);

		void track_hpatch_E_min();

		void update_internal_energy_totals_hpatch();

	private:
		void reset_from_previous_deltaHpatch_comp();
		bool decide_procrastinate_hpatch_computations( Real const pd_deltaE, Real const threshold ) const;

		void detect_background_residue_and_first_class_residue_overlap();
		void initialize_bg_bg_overlaps();
		void initialize_bg_bg_atom_atom_overlaps();
		utility::vector1< utility::vector1< bool > > const & get_bg_bg_atom_atom_overlaps( Size node1_index, Size node2_index );

		void init_SASA_radii_from_database();

#ifdef DOUBLE_CHECK_COUNTS
		void verify_sasas_correct();
#endif

		HPatchInteractionGraph();
		HPatchInteractionGraph( HPatchInteractionGraph< V, E, G > const & );
		HPatchInteractionGraph< V, E, G > & operator = ( HPatchInteractionGraph< V, E, G > const & );

	private:
		pose::PoseOP pose_;
		task::PackerTaskOP packer_task_;
		rotamer_set::RotamerSetsOP rotamer_sets_;

		core::Real hpatch_score_weight_;

		Size num_total_residues_;
		Size num_residues_assigned_as_background_;
		utility::vector1< Size > resid_2_bgenumeration_;
		utility::vector1< Size > bgenumeration_2_resid_;

		utility::vector1< utility::vector1< Size > > bg_bg_respairs_w_hphobe_olap_;
		utility::vector1< utility::vector1< utility::vector1< utility::vector1< bool > > > > bg_bg_atom_atom_overlaps_;
		utility::vector1< utility::vector1< utility::vector1< utility::vector1< bool > > > > curr_bg_bg_exhpobeolap_;
		utility::vector1< utility::vector1< utility::vector1< utility::vector1< bool > > > > alt_bg_bg_exhpobeolap_;

		// keep track of which residues are effected by a considered rotamer substitution
		bool some_node_in_state_0_;
		utility::vector1< Size > fc_nodes_near_rotsub_;
		utility::vector1< bool > fc_nodes_near_rotsub_bool_;
		utility::vector1< Size > bg_nodes_near_rotsub_;
		utility::vector1< bool > bg_nodes_near_rotsub_bool_;

		/// "djs" for DisJoint Sets.  An abuse of the idea of an acronym.
		utility::vector1< Size > fc_exp_hphobe_djs_offsets_;
		utility::vector1< Size > fc_n_exp_hphobes_;
		utility::vector1< Size > bg_exp_hphobe_djs_offsets_;
		utility::vector1< Size > bg_n_exp_hphobes_;

		utility::vector1< exposed_hydrophobic_data > djs_id_2_hphobe_index_;
		//utility::vector1< Real > sasa_for_djs_node_;
		utility::vector1< Real > ep_sasa_for_djs_node_;

		/// extract the disjoint sets connected component information by hand (faster)
		utility::vector1< Size > reps_for_nonzero_rank_ccs_;
		utility::vector1< Size > djs_rep_node_index_2_cc_index_;
		utility::vector1< Real > sasa_for_cc_;

		bool prepared_for_simulated_annealing_;

		bool observed_sufficient_hpatch_E_to_predict_min_;

		static Size num_state_substitutions_considered_;
		static Size num_hpatch_comps_procrastinated_;
		static Size num_hpatch_comps_later_made_;

		Real hpatch_score_min_last_100_;
		Real hpatch_score_min_recent_;
		Size num_substitutions_since_hpatch_min_update_;

		bool calculated_hpatch_deltaE_;
		core::PackerEnergy deltaE_for_substitution_;

		Size node_considering_alt_state_;
		Size alt_state_being_considered_;
		Real total_energy_current_state_assignment_;
		Real total_energy_alternate_state_assignment_;

		Real hpatch_energy_current_state_assignment_;
		Real hpatch_energy_alternate_state_assignment_;

		int num_commits_since_last_update_;
		float deltaE_threshold_for_avoiding_hpatch_calcs_;

		static utility::vector1< Real > radii_;
		static bool initialized_SASA_radii;
		static const int COMMIT_LIMIT_BETWEEN_UPDATES = 4096;

		static std::string carbon_atom;
		static std::string sulfur_atom;
};



//----------------------------------------------------------------------------//
//-------------------------------- HPatch Node Class ---------------------------//
//----------------------------------------------------------------------------//

///
/// @begin HPatchNode< V, E, G >::HPatchNode
///
/// @brief
/// HPatchNode constructor
///
template < typename V, typename E, typename G >
HPatchNode< V, E, G >::HPatchNode( G* owner, int node_index, int num_states ) :
	FirstClassNode< V, E, G > ( owner, node_index, num_states ),
	current_state_rotamer_dots_(), // does this call the RD constructor? Yes!
	alt_state_rotamer_dots_(),
	alt_state_dots_matches_current_state_dots_( true )
{
	rotamers_vector_.resize( num_states );
	restype_group_for_rotamers_.resize( num_states, 0 );
	self_and_bg_dots_for_states_.resize( parent::get_num_states() );
}


///
/// @begin HPatchNode< V, E, G >::~HPatchNode
///
/// @brief
/// destructor -- no dynamically allocated data, does nothing
///
template < typename V, typename E, typename G >
HPatchNode< V, E, G >::~HPatchNode() {
	//TR_NODE << "called destructor" << std::endl;
}


///
/// @begin HPatchNode::set_rotamers
///
/// @details
/// Need to save a reference to the rotamer_set so that we can determine what a particular state change will do to the score
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::set_rotamers( rotamer_set::RotamerSetCOP rotamers ) {

	// get_num_states should call the parent graph's method?
	if ( rotamers->num_rotamers() != (Size) parent::get_num_states() ) {
		utility_exit_with_message( "Number of rotamers is not equal to parents number of states. Quitting.");
	}

#ifdef FILE_DEBUG
	TR_NODE << "set_rotamers: adding " << rotamers->num_rotamers() << " to local rotamers vector." << std::endl;
#endif

	for ( Size ii = 1; ii <= rotamers->num_rotamers(); ++ii ) {
		rotamers_vector_[ ii ] = rotamers->rotamer( ii );
		restype_group_for_rotamers_[ ii ] = rotamers->get_residue_group_index_for_rotamer( ii );
	}

	{ // scope -- I should make this a function
	std::string carbon( "C" ), sulfur( "S" );
	max_hphobes_ = 0;

	// here we set every index of hphobe_ats_for_restype_group_ with the atom ids of the hydrophobic atoms in that restype.
	// this is done by first looking at that residue type and counting up the number of C and S atoms to reserve that number
	// of elements in the inner vector and then going back over the list of heavy atoms and setting the index of that atom id.
	// this can't be done in one loop because we also want to set the max_hphobes_ variable and this function only runs once
	// at initializations so it's ok.

	hphobe_ats_for_restype_group_.resize( restype_group_for_rotamers_[ rotamers->num_rotamers() ] );
	for ( Size ii = 1; ii <= hphobe_ats_for_restype_group_.size(); ++ii ) {
		chemical::ResidueType const & ii_restype = rotamers_vector_[ rotamers->get_residue_type_begin( ii ) ]->type();
		Size count_cs = 0;
		for ( Size jj = 1; jj <= ii_restype.nheavyatoms(); ++jj ) {
			if ( ii_restype.atom_type( jj ).element() == carbon || ii_restype.atom_type( jj ).element() == sulfur ) {
				++count_cs;
			}
		}
		if ( max_hphobes_ < count_cs ) {
			max_hphobes_ = count_cs;
		}
		hphobe_ats_for_restype_group_[ ii ].reserve( count_cs );
		for ( Size jj = 1; jj <= ii_restype.nheavyatoms(); ++jj ) {
			if ( ii_restype.atom_type( jj ).element() == carbon || ii_restype.atom_type( jj ).element() == sulfur ) {
				hphobe_ats_for_restype_group_[ ii ].push_back( jj );
			}
		}
	}
	} // end scope

	curr_state_exp_hphobes_.reserve( max_hphobes_ );
	alt_state_exp_hphobes_.reserve( max_hphobes_ );

#ifdef FILE_DEBUG
	/*TR_NODE << "rotamers_vector_: [ ";
	for ( Size ii = 1; ii <= rotamers_vector_.size(); ++ii ) {
		TR_NODE << ii << ":" << rotamers_vector_[ ii ]->name1() << "-" << rotamers_vector_[ ii ]->seqpos() << ", ";
	}
	TR_NODE << "]" << std::endl;*/
#endif

}

///
/// @begin HPatchNode::get_rotamer
///
/// @brief
/// Returns a constant OP to the rotamer/residue object for the given state.
///
/// @details
/// Need to save a reference to the rotamer_set so that we can determine what a particular state change will do to the score
///
template < typename V, typename E, typename G >
conformation::ResidueCOP
HPatchNode< V, E, G >::get_rotamer( int state ) const {
	return rotamers_vector_[ state ];
}


///
/// @begin HPatchNode::set_rotamer_dots_for_state
///
/// @brief
/// stores the coordinates for a state.
///
/// @detailed
/// Currently RotamerCoords stores the atoms in the order they are created for the pose. In the future, this might be
/// changes to store the atoms in special trie order to reduce the number of sphere calculations necessary.
///
/// self_and_bg_dots_for_states_ is a vector1 of RotamerDots objects. the size of the vector is set during Node
/// construction to parent::get_num_states().
///
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::set_rotamer_dots_for_state( Size state, RotamerDots const & rd ) {
	assert( state > 0 && state <= (Size)parent::get_num_states() );
	self_and_bg_dots_for_states_[ state ] = rd;
}


///
/// @begin HPatchNode::overlaps
///
/// @brief
/// returns true if any sphere for any atom of any state on this vertex overlaps with any atom on any sphere on a neighboring vertex.
///
/// @detailed
/// neighbor - [in] - the vertex that neighbors this vertex
///
/// @remarks
/// This method could be faster if I built a rotamer trie for each vertex, so I've got code in here to measure just how
/// much time I spend on this task. In a complete redesign of UBQ, I spend 8 to 10 seconds in this method.
/// In comparison, I spend 15 minutes in simulated annealing. (apl)
///
template < typename V, typename E, typename G >
bool HPatchNode< V, E, G >::overlaps( HPatchNode< V, E, G >* neighbor ) {

	//static double time_spent_here( 0.0 );
	bool found_overlap = false;

	//TR_NODE << "overlaps(): num_states: " << parent::get_num_states() << ", neighbors num_states: " << neighbor->get_num_states() << std::endl;

	//std::clock_t starttime = clock();
	for ( Size ii = 1; ii <= (Size)parent::get_num_states(); ++ii ) {
		//for (Size jj = 1; jj <= (Size)parent::get_num_states(); ++jj ) { // WRONG!
		for (Size jj = 1; jj <= (Size)neighbor->get_num_states(); ++jj ) {
			if ( self_and_bg_dots_for_states_[ ii ].overlaps( neighbor->self_and_bg_dots_for_states_[ jj ] ) ) {
				found_overlap = true;
				break;
			}
		}
		if (found_overlap)
			break;
	}
	//std::clock_t stoptime = clock();
	//time_spent_here += ((double)(stoptime - starttime)) / CLOCKS_PER_SEC;
	//TR_NODE << "time spent computing overlaps: " << time_spent_here << ", found overlap? " << found_overlap << std::endl;

	return found_overlap;
}


///
/// @begin HPatchNode::detect_any_overlap_with_rotamer
///
/// @brief
/// determines if any atom for any rotamer of this vertex overlaps with any atom from some background residue.
/// called by BGNodes for the detect_background_residue_and_first_class_residue_overlap phase of the prep for simA
/// call in the HPatchIG.
///
template < typename V, typename E, typename G >
bool HPatchNode< V, E, G >::detect_any_overlap_with_rotamer( RotamerDots const & bg_dots ) const {

	// for every rotamer possible at this molten res position/node...
	for ( Size ii = 1; ii <= (Size)parent::get_num_states(); ++ii ) {
		// RotamerDots objects know how to determine if they overlap with another RotamerDots object
		if ( self_and_bg_dots_for_states_[ ii ].overlaps( bg_dots ) ) {
			//TR_NODE << "state " << ii << " of node " << parent::get_node_index() << " overlaps with bg node." << std::endl;
			return true; // as soon as one overlap is found, this method immediately returns
		}
	}
	return false;
}


///
/// @begin HPatchNode::prepare_for_simulated_annealing
///
/// @brief
/// invokes V's prep_for_simA method. Also computes the "self-overlap" that each of its rotamers has.
///
template < typename V, typename E, typename G >
void
HPatchNode< V, E, G >::prepare_for_simulated_annealing() {

#ifdef FILE_DEBUG
	TR_NODE << "prepare_for_simulated_annealing(): calling V::prep for simA for node " << parent::get_node_index() << std::endl;
#endif

	// parent is AddtlBGNodesIG; V is either PDIG or LinmemIG; here we want to call the templated classes prep_for_simA
	// method so that it can get ready for all the pairwise-decomposable energy term stuff
	// PDIG::prep_for_simA() calls update_edge_vectors which eventually leads to a seq fault when it requests the edge
	// table pointer.
	V::prepare_for_simulated_annealing();

	if ( ! parent::get_bg_edge_vector_up_to_date_() ) {
#ifdef FILE_DEBUG
		TR_NODE << "prepare_for_simulated_annealing(): bg edge vector not up to date, calling parent::update_bg_edge_vector()" << std::endl;
#endif
		parent::update_bg_edge_vector();
	}

#ifdef FILE_DEBUG
	TR_NODE << "prepare_for_simulated_annealing(): calling initialize_self_overlap()" << std::endl;
#endif
	initialize_self_overlap();
#ifdef FILE_DEBUG
	TR_NODE << "prepare_for_simulated_annealing(): ... done incrementing self overlap for all states" << std::endl;
#endif

	initialize_atom_atom_overlap_cache();
}


///
/// @begin HPatchNode::initialize_self_overlap
///
/// @brief
/// Initializes rotamer self overlap; called by HPatchInteractionGraph;:prepare_for_sim_annealing right before simulated
/// annealing begins.
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::initialize_self_overlap() {

	for ( Size ii=1; ii <= (Size)parent::get_num_states(); ++ii ) {
		//TR_NODE << "initialize_self_overlap(): self_and_bg_dots_for_states_[ " << ii << " ]..." << std::endl;
		self_and_bg_dots_for_states_[ ii ].increment_self_overlap();
	}

}


///
/// @begin HPatchNode::initialize_atom_atom_overlap_cache
///
/// @brief
/// Initializes the atom_atom_overlaps vector.
///
/// @detailed
/// The atom_atom_overlaps vector stores a boolean for the intra-residue atom-atom overlap for every state possible at
/// this Node. During simulated annealing, the IG has to determine the connected components after every sub. To do this,
/// it has to know which atoms are exposed as well as overlapping. Instead of recomputing atom-atom overlaps after every
/// sub, do it once and cache the values on the Node.
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::initialize_atom_atom_overlap_cache() {

	// resize the outer-most vector to the number of states possible
	self_atom_atom_overlaps_.resize( (Size)parent::get_num_states() );

	for ( Size state = 1; state <= (Size)parent::get_num_states(); ++state ) {

		//get_rotamer( state ) returns a ResidueCOP; we can use that to get residue info
		conformation::ResidueCOP rotamer_ = get_rotamer( state );

		self_atom_atom_overlaps_[ state ].resize( rotamer_->nheavyatoms(), utility::vector1< bool >( rotamer_->nheavyatoms(), false ) );

		std::string carbon_atom = "C";
		std::string sulfur_atom = "S";
		Real const probe_radius = 1.4;

		RotamerDots const & rd = self_and_bg_dots_for_states_[ state ]; // the same RD object can be used for both res ii and jj

		for ( Size iia=1; iia <= rotamer_->nheavyatoms(); ++iia ) {
			// immediately continue if not a hydrophobic atom; no point in computing overlaps for polar atoms
			if ( rotamer_->atom_type( iia ).element() != carbon_atom && rotamer_->atom_type( iia ).element() != sulfur_atom )
				continue;

			Real const iia_atom_radius = rd.get_atom_radius( iia ) + probe_radius;

			// for intra-residue we only have iterate over the greater-indexed heavyatoms
			for ( Size jja = iia + 1; jja <= rotamer_->nheavyatoms(); ++jja ) {
				if ( rotamer_->atom_type( jja ).element() != carbon_atom && rotamer_->atom_type( jja ).element() != sulfur_atom )
					continue;

				Real const jja_atom_radius = rd.get_atom_radius( jja ) + probe_radius;

				// check if the two atoms overlap, and if that degree of overlap exceeds the threshold
				Vector const & iia_atom_xyz = rotamer_->atom( iia ).xyz();
				Vector const & jja_atom_xyz = rotamer_->atom( jja ).xyz();

				// exit if large probe radii do not touch
				Real const distance_squared = iia_atom_xyz.distance_squared( jja_atom_xyz );
				//if ( distance_squared <= (iia_atom_radius + iia_atom_radius) * (jja_atom_radius + jja_atom_radius) ) {
				if ( distance_squared <= (iia_atom_radius + jja_atom_radius) * (iia_atom_radius + jja_atom_radius) ) {

					Real const distance_ijxyz = std::sqrt( distance_squared );
					int degree_of_overlap;
					core::scoring::get_overlap( iia_atom_radius, jja_atom_radius, distance_ijxyz, degree_of_overlap );
					if ( degree_of_overlap >= 15 ) {

						self_atom_atom_overlaps_[ state ][ iia ][ jja ] = true;
#ifdef FILE_DEBUG
						//TR_NODE << "initialize_atom_atom_overlap_cache(): overlapping intra-residue atom pair: "
						//	<< rotamer_->seqpos() << "/" << utility::trim( rotamer_->atom_name( iia ) ) << " - " << rotamer_->seqpos() << "/" << utility::trim( rotamer_->atom_name( jja ) )
						//	<< ", degree of overlap: " << degree_of_overlap << std::endl;
#endif

					}
				} // end if distance
			} // for loop over all jj heavyatoms
		} // for loop over all ii heavyatoms
	} // for loop over all states


/*
#ifdef FILE_DEBUG
	TR_NODE << "node " << parent::get_node_index() << std::endl;
	for ( Size state = 1; state <= (Size)parent::get_num_states(); ++state ) {
		TR_NODE << "state " << state << " self_atom_atom_overlaps_: [ " << std::endl;
		for ( Size aa=1; aa <= self_atom_atom_overlaps_[ state ].size(); ++aa ) {
			TR_NODE << "self_atom_atom_overlaps_[ " << aa << " ]: [ ";
			for ( Size bb=1; bb <= self_atom_atom_overlaps_[ state ][ aa ].size(); ++bb ) {
				TR_NODE << self_atom_atom_overlaps_[ state ][ aa ][ bb ] << ", ";
			}
			TR_NODE << "]" << std::endl;
		}
		TR_NODE << "]" << std::endl;
	}
	TR_NODE << std::endl;
#endif
*/

}


///
/// @begin HPatchNode< V, E, G >::initialize_overlap_with_background
///
/// @brief
/// This method computes and caches each set of overlaps this node's states have with a background residue.
///
/// @detailed
/// The node stores these overlaps in two places: 1) in the input vector from the HPatchBackgroundEdge object, and 2) in
/// its own array of dot coverage counts for the self-and-background sphere overlaps.
///
/// @param
/// bg_rotamer - [in] - the RotamerDots object for the background residue
/// node_states_coverage_of_bg_res_ -[out] - the array of RotamerDotsCache objects that the HPatchBackgroundEdge stores for the overlap of each state in this residue on the background residue.
/// states_atom_atom_overlap_on_bg_res -[out] - an array that holds the atom-atom overlap each state on this Node has on the background residue
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::initialize_overlap_with_background(
	RotamerDots const & bg_rotamer_dots,
	std::vector< RotamerDotsCache > & node_states_coverage_of_bg_res_,
	std::vector< utility::vector1< utility::vector1< bool > > > & states_atom_atom_overlap_on_bg_res
) {

	std::vector< RotamerDotsCache* > v_dots_for_bg_res( (Size)parent::get_num_states() + 1, static_cast< RotamerDotsCache* >(0) );

	// since there's no way to default init the atom-atom overlap vector, we have to make sure we set every index, even 0
	utility::vector1< utility::vector1< bool > > zero_vector( bg_rotamer_dots.get_num_atoms(), utility::vector1< bool >( bg_rotamer_dots.get_num_atoms(), false ) );
	states_atom_atom_overlap_on_bg_res[ 0 ] = zero_vector;

	// too much output during initialization
	//TR_NODE << "initialize_overlap_with_background(): states on FC node " << parent::get_node_index() << " which overlap bg residue include: ";
	for ( Size ii = 1; ii <= (Size)parent::get_num_states(); ++ii ) {

		// now resize the atom-atom overlap vector to the right sizes. the values will get filled by the Dots method increment_this_and_cache().
		Size ii_state_num_atoms = self_and_bg_dots_for_states_[ ii ].get_num_atoms();
		states_atom_atom_overlap_on_bg_res[ ii ].resize( ii_state_num_atoms, utility::vector1< bool >( bg_rotamer_dots.get_num_atoms(), false ) );

		// do any of the states possible at this FCNode overlap with the bg_rotamer?
		if ( self_and_bg_dots_for_states_[ ii ].overlaps( bg_rotamer_dots ) ) {
			v_dots_for_bg_res[ ii ] = new RotamerDotsCache( bg_rotamer_dots.get_num_atoms() );
			//self_and_bg_dots_for_states_[ ii ].increment_this_and_cache( bg_rotamer_dots, *v_dots_for_bg_res[ ii ], states_atom_atom_overlap_on_bg_res[ ii ] );
			// can't just make this increment_this_and_cache call because then the atom_atom_overlap array is in the wrong "order"

			utility::vector1< utility::vector1< bool > > temp_atat_overlap( bg_rotamer_dots.get_num_atoms(), utility::vector1< bool >( ii_state_num_atoms, false ) );
			self_and_bg_dots_for_states_[ ii ].increment_this_and_cache( bg_rotamer_dots, *v_dots_for_bg_res[ ii ], temp_atat_overlap );

			/// now transpose the table for future use.
			for ( Size jj = 1; jj <= bg_rotamer_dots.get_num_atoms(); ++jj ) {
				for ( Size kk = 1; kk <= ii_state_num_atoms; ++kk ) {
					states_atom_atom_overlap_on_bg_res[ ii ][ kk ][ jj ] = temp_atat_overlap[ jj ][ kk ];
				}
			}

			// now save this overlap to the passed in (by reference) vector - BUT ALLOCATE SPACE FOR IT FIRST
			node_states_coverage_of_bg_res_[ ii ] = *(v_dots_for_bg_res[ ii ]);
			//TR_NODE << ii << ", ";

			// 'delete' RotamerDotsCache objects; frees the memory allocated to RotamerDotsCache objects created with 'new' above
			delete v_dots_for_bg_res[ ii ];
		}

	}
	//TR_NODE << std::endl;

}


///
/// @begin HPatchNode::assign_zero_state
///
/// @brief
/// Assign the node to state 0 -- the "unassigned" state.
///
/// @detailed
/// A node in state 0 produces no hpatch score.  Its neighbors have to adjust their scores appropriately. This method
/// iterates over all the edges emanating from this node and tells them to acknowledge that they've been zeroed out.
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::assign_zero_state() {

	// depending on the type of interaction graph (linmem, standard pd) used for the simulation, this will call
	// the appropriate assign_zero_state() method and initialize the one- and two-body energy term variables to zero
	// parent refers to the AddtlBGNodesIG; G refers to either PDIG or LinmemIG.  AddtlBGNodesIG doesn't define
	// a assign_zero_state() method! But it extends from the other two classes so it will eventually find the right method.
	parent::assign_zero_state();  // was this a bug previously?  Shouldn't it be V? apparently not...
	//V::assign_zero_state();

#ifdef FILE_DEBUG
	TR_NODE << "assign_zero_state(): calling acknowledge_state_zeroed() on all incident edges." << std::endl;
#endif
	for ( Size ii = 1; ii <= (Size)parent::get_num_incident_edges(); ++ii ) {
		get_incident_hpatch_edge( ii )->acknowledge_state_zeroed( parent::get_node_index() );
	}

	parent::update_bg_edge_vector();

#ifdef FILE_DEBUG
	TR_NODE << "assign_zero_state(): calling acknowledge_state_change() on all incident bg edges." << std::endl;
#endif
	for ( Size ii = 1; ii <= (Size)parent::get_num_edges_to_background_nodes(); ++ii ) {
		get_edge_to_hpatch_bg_node( ii )->acknowledge_state_change( 0 );
	}

	RotamerDots rd; // create an empty, uninit'd RotamerDots object
	current_state_rotamer_dots_ = rd;
	alt_state_rotamer_dots_ = rd;
	alt_state_dots_matches_current_state_dots_ = true;
	// is setting alt_state dots necessary? Yes! Because when other FC nodes enter the unassigned state, they will
	// trigger calls to acknowledge_neighbors_sub which sets current state to alt state. in theory, I could even put
	// just alt_state = emtpy RD and the current state dots would eventually get overwritten.

	get_hpatch_owner()->register_fc_node_in_state0();

}


///
/// @begin HPatchNode::acknowledge_neighbors_substitution
///
/// @brief
/// bookkeeping to follow a neighbors state substitution.  this method gets called when a HPatchNode commits a sub
/// and then broadcasts that change to all its neighboring fc nodes via the incident HPatchEdges. basically we need
/// to set current state equal to alt state here. (Hopefully alt state is still correct!!)  Since there's no way for
/// a HPatchNode to know what other HPatchNodes are connected to it except via HPatchEdges, the calls seem a bit
/// complicated.  A HPatchNode has to call acknowledge_state on each Edge.  The Edges have to figure out which Node
/// is changing/not changing and then they call an inform_non_changing node of change method.  That method then makes
/// the call to this method on the correct HPatchNode. The inform_non_changing method can not be removed, because
/// it's used during the substitution evaluations as well.
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::acknowledge_neighbors_substitution() {

#ifdef FILE_DEBUG
	// too much output...
	//if ( current_state_rotamer_dots_ != alt_state_rotamer_dots_ ) {
	//	TR_NODE << "acknowledge_neighbors_substitution - node " << parent::get_node_index() << " changing curr_state_dots from "
	//		<< current_state_rotamer_dots_ << " to " << alt_state_rotamer_dots_ << std::endl;
	//}
#endif

	current_state_rotamer_dots_ = alt_state_rotamer_dots_; // what happens here in the 0-state case?
	curr_state_inv_dots_ = alt_state_inv_dots_;
	curr_state_exp_hphobes_ = alt_state_exp_hphobes_;
	alt_state_dots_matches_current_state_dots_ = true;
}


///
/// @begin HPatchNode::get_current_state_sasa()
///
/// @brief
/// returns the amount of sasa this node has in its current state assignment
///
/// @detailed
/// RotamerDots objects know when they are in the unassigned state, so nothing special needs to be done to handle
/// the 0 state. But, we could also query parent using get_current_state() to check if it's nonzero as an alternative.
///
template < typename V, typename E, typename G >
Real HPatchNode< V, E, G >::get_current_state_sasa() const {
	return current_state_rotamer_dots_.get_sasa();
}


///
/// @begin HPatchNode::get_current_state_sasa()
///
/// @brief
/// Returns the current state SASA for the passed in atom index.
///
template < typename V, typename E, typename G >
Real HPatchNode< V, E, G >::get_current_state_sasa( Size atom_index ) const {
	return current_state_rotamer_dots_.get_atom_sasa( atom_index );
}


///
/// @begin HPatchNode::get_alternate_state_sasa()
///
/// @brief
/// Returns the alternate state SASA for the passed in atom index.
///
template < typename V, typename E, typename G >
Real HPatchNode< V, E, G >::get_alternate_state_sasa( Size atom_index ) const {
	return alt_state_rotamer_dots_.get_atom_sasa( atom_index );
}


///
/// @begin HPatchNode< V, E, G >::get_atom_atom_self_overlaps_for_state
///
/// @brief
/// Returns a const reference to the atom-x-atom-pair vector-of-vectors of bools that specifies which atoms are overlapping, in the given state.
///
template < typename V, typename E, typename G >
utility::vector1< utility::vector1 < bool > > const &
HPatchNode< V, E, G >::get_atom_atom_self_overlaps_for_state( Size state ) const {
	return self_atom_atom_overlaps_[ state ];
}


///
/// @begin HPatchNode::calculate_PD_deltaE_for_substitution
///
/// @brief
/// Returns the change in energy induced by changing a node from its current state into some alternate state for the PD energy terms only.
///
/// @detailed
/// This function always gets called for every substitution. Only the consider_alt_state() call can get procrastinated.
///
template < typename V, typename E, typename G >
core::PackerEnergy HPatchNode< V, E, G >::calculate_PD_deltaE_for_substitution( int alternate_state, core::PackerEnergy & prev_PDenergies_for_node ) {

	assert( alternate_state > 0  && alternate_state <= parent::get_num_states() );

	prev_PDenergies_for_node = parent::get_curr_pd_energy_total();

	// have base class perform deltaE computations for PD portion of energy function
	// parent refers to AddtlBGNodesIG; G refers to either PDIG or LinmemIG.  Which do I want here? AddtlBGNodesIG doesn't
	// implement the calc_deltaEpd method but it extends either PDIG or LinmemIG so the runtime will eventually find the
	// right function to call.
	parent::calc_deltaEpd( alternate_state );

	return get_pd_energy_delta();

}


///
/// @begin HPatchNode< V, E, G >::get_pd_energy_delta()
///
/// @brief
/// Returns the deltaE for just the PD terms. Separate method from the one above because this one can be called from within
/// a commit_sub call that didn't go through consider_sub().
///
template < typename V, typename E, typename G >
core::PackerEnergy HPatchNode< V, E, G >::get_pd_energy_delta() {
	return ( parent::get_alt_pd_energy_total() - parent::get_curr_pd_energy_total() );
}


///
/// @begin HPatchNode< V, E, G >::consider_alternate_state()
///
/// @brief
/// Instructs the Node to update the alt state information held by it and its neighbors in response to switching from the
/// current state to an alternate state.
///
template < typename V, typename E, typename G >
Real HPatchNode< V, E, G >::consider_alternate_state() {

	get_hpatch_owner()->register_fc_node_affected_by_rotsub( parent::get_node_index() );

#ifdef FILE_DEBUG
	TR_NODE << "consider_alternate_state(): current_state:" << parent::get_current_state();
	if ( parent::get_current_state() == 0 ) {
		TR_NODE << " (" << wt_residue_for_node().name3() << "-";
		if ( wt_residue_for_node().is_polar() ) { TR_NODE << "P)"; } else { TR_NODE << "HP)"; }
	} else {
		TR_NODE << " (" << get_rotamer( parent::get_current_state() )->name() << "-";
		if ( get_rotamer( parent::get_current_state() )->is_polar() ) { TR_NODE << "P)"; } else { TR_NODE << "HP)"; }
	}
	TR_NODE << ", alternate_state: " << parent::get_alternate_state() << " (" << get_rotamer( parent::get_alternate_state() )->name() << "-";
	if ( get_rotamer( parent::get_alternate_state() )->is_polar() ) { TR_NODE << "P)"; } else { TR_NODE << "HP)"; }
	TR_NODE << std::endl;
#endif

	// self_and_bg_dots_for_states has the RotamerDots object for a state, with the overlap counts set for all BG nodes
	// that overlap with that state
	//runtime_assert_msg( alt_state_dots_matches_current_state_dots_, "alt_state_dots and current_state_dots do not match for FC node " + ObjexxFCL::format::I( 3, parent::get_node_index() ) ); // causes crash in pmut_scan run
	assert( alt_state_dots_matches_current_state_dots_ );
	alt_state_dots_matches_current_state_dots_ = false;

	alt_state_rotamer_dots_ = self_and_bg_dots_for_states_[ parent::get_alternate_state() ];

	Real delta_sasa = 0;
	for ( Size ii = 1; ii <= (Size)parent::get_num_incident_edges(); ++ii ) {
		delta_sasa += get_incident_hpatch_edge(ii)->update_state_at_neighbor( parent::get_node_index(), parent::get_alternate_state(), alt_state_rotamer_dots_ );
	}
	for ( Size ii = 1; ii <= (Size)parent::get_num_edges_to_background_nodes(); ++ii ) {
		delta_sasa += get_edge_to_hpatch_bg_node( ii )->update_state_at_neighbor( parent::get_alternate_state() );
	}

	update_alt_state_exphphobes();

	//alt_state_inv_dots_.setup_from_rotamer_dots( alt_state_rotamer_dots_ );
	alt_state_inv_dots_.setup_from_rotamer_dots( alt_state_rotamer_dots_, alt_state_exp_hphobes_ );

#ifdef FILE_DEBUG
	if ( get_sasa_difference() != 0 ) {
		TR_NODE << "consider_alternate_state(): curr state dots: " << current_state_rotamer_dots_.get_sasa()
			<< ", alt state dots: " << alt_state_rotamer_dots_.get_sasa()
			<< ", sasa difference: " << get_sasa_difference() << std::endl;
	}
#endif
	delta_sasa += get_sasa_difference();

	return delta_sasa;
}


///
/// @begin HPatchNode< V, E, G >::update_alt_state_exphphobes
///
/// @brief
/// Updates the vector alt_state_exp_hphobes_ by checking the sasa of every atom in the residue currently on this Node.
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::update_alt_state_exphphobes() {
	utility::vector1< Size > const & alt_hphobes( alt_state_hphobes() ); // get all the hydrophobic atoms in this residue first
	alt_state_exp_hphobes_.clear();
	for ( Size ii = 1; ii <= alt_hphobes.size(); ++ii ) { // then figure out if they're exposed by checking their SASA
		Size const ii_atom = alt_hphobes[ ii ];
		if ( alt_state_rotamer_dots_.get_atom_sasa( ii_atom ) > 0.0 ) {
			alt_state_exp_hphobes_.push_back( ii_atom );
		}
	}
}


///
/// @begin HPatchNode< V, E, G >::update_state_for_neighbors_substitution
///
/// @brief
/// returns the change in sasa for this node induced by a state substitution at a neighboring node.  The node
/// increments the dot coverage count for the RotamerDots object representing the alternate state dots for that neighbor.
///
/// @detailed
/// The procedure is simple at heart -- the caching makes it complex.
///
/// alt_state_rotamer_dots_ = current_state_rotamer_dots_;      //copy dot coverage counts
/// alt_state_rotamer_dots_.decrement( neighbors_curr_state_overlap_with_this );
/// alt_state_rotamer_dots_.increment_both( neighbors_alternate_state );
/// return ( alt_state_rotamer_dots_.get_score() - current_state_rotamer_dots_.get_score() );
///
/// Extensive caching techniques save time. Each HPatchEdge stores the dot coverage for each pair of HPatchNodes in their
/// current state. These are stored in the RotamerDotsCache objects. Let's name the vertices: this vertex is vertex B.
/// Vertex B is projecting its hpatch deltaE while vertex A is considering a substitution from one state to another.
/// The edge connecting A and B provides node B with the set of masks for the overlap by the atoms of A's alternate
/// state on all atoms of vertex B.
///
/// @param
/// neighbors_alternate_state - [in/out] - the RotamerDots object representing the alternate state vertex A (the neighbor) is considering.
/// neighbors_curr_state_overlap_with_this - [in] - the RotamerDotsCache object held by the HPatchEdge connecting A and B that represents the coverage of A's current state of B's atoms.
/// this_overlap_with_neighbors_alternate - [out] - the dot cache to store the results of computing the overlap of B's current state on A's alt state
/// neighbors_alternate_overlap_with_this - [out] - the dot cache to store the results of computing the overlap of A's alt state on B's current state
/// mask_this_covered_by_other - [in/out] - starting set of overlaps
/// mask_other_covered_by_this - [in/out] - starting set of overlaps
/// num_atoms_same - [in] - the number of atoms shared between A's current and alternate states.
///
/// What do we do when the current state is unassigned. Then that results in the alt_state starting from being unassigned
/// and we get problems because the RotamerDots object doesn't have any memory assigned to it if it's in the unassigned state.
///
/// it's probable that a sub that was considered before wasn't committed, so the alt state count at this node needs to go
/// back to what the current state count is.  the problem that may crop up with the reset here is that on a previous consider
/// call, a set of nodes will update their counts.  if commit is not called, and a second consider is called, then this node
/// will have its alt state count reset but what about all the other nodes in the previous set.  how will their counts get
/// reset to the current state count?  perhaps a check in the commit method can be added.  nope, a new method has been
/// added to Nodes and BGNodes to handle this case.
///	alt_state_total_hASA_ = curr_state_total_hASA_;
///
///
template < typename V, typename E, typename G >
Real HPatchNode< V, E, G >::update_state_for_neighbors_substitution (
	HPatchNode<V,E,G>* node_considering_substitution,
	RotamerDots & neighbors_alternate_state,
	RotamerDotsCache const & neighbors_curr_state_overlap_with_this,
	RotamerDotsCache & this_overlap_with_neighbors_alternate,
	RotamerDotsCache & neighbors_alternate_overlap_with_this,
	utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache
)
{

	get_hpatch_owner()->register_fc_node_affected_by_rotsub( parent::get_node_index() );
	parent::set_alternate_state( parent::get_current_state() );

	// don't bother doing sasa updates if this node is in the unassigned state. things really only start to get interesting
	// once all nodes have something for their current state.
	if ( current_state_rotamer_dots_.state_unassigned() ) {
#ifdef FILE_DEBUG
		TR_NODE << "update_state_for_neighbors_substitution(): node " << parent::get_node_index() << " in unassigned state. skipping sasa calculations." << std::endl;
#endif
		return 0.0;
	}

#ifdef FILE_DEBUG
	//TR_NODE << "update_state_for_neighbors_substitution(): current_state_rotamer_dots_" << std::endl;
	//current_state_rotamer_dots_.print( std::cout );
#endif

	// everything has to start from the current state. then we'll update the alt state to be correct.
	//runtime_assert_msg( alt_state_dots_matches_current_state_dots_, "alt_state_dots and current_state_dots do not match for FC node " + ObjexxFCL::format::I( 3, parent::get_node_index() ) ); // causes crash in pmut_scan run
	assert( alt_state_dots_matches_current_state_dots_ );
	/// APL -- this should be unnecessary
	/// APL TEMP alt_state_rotamer_dots_ = current_state_rotamer_dots_;
	alt_state_dots_matches_current_state_dots_ = false;

#ifdef FILE_DEBUG
	//TR_NODE << "update_state_for_neighbors_substitution(): set alt_state dots to curr_state dots: " << std::endl;
	//alt_state_rotamer_dots_.print( std::cout );
#endif

	// at the beginning of annealing, nodes are in the unassigned state.
	// there will be times when the overlap cached on the edge between this node and the one undergoing substitution will not be
	// initialized. those times will be when the node undergoing substitution is going from state 0 to some other state.
	// in that case, we have to skip this call because the neighbors curr_state overlap is uninit'd and will results in bus errors.

	if ( node_considering_substitution->get_current_state() != 0 ) {
		// subtract out the cached overlap that the current state (on the other node) had with this node
		alt_state_rotamer_dots_.decrement_from_cached( neighbors_curr_state_overlap_with_this );
#ifdef FILE_DEBUG
		//TR_NODE << "update_state_for_neighbors_substitution(): decrementing neighbors_curr_state_overlap_with_this" << std::endl;
		//neighbors_curr_state_overlap_with_this.print( std::cout );
#endif
	}

	// then add in the overlap that the alt state (on the other node) would have on this nodes state and cache the overlap.
	// the outer vector of alt_state_atom_atom_overlaps_cache has the changing nodes atoms, the inner vector has this nodes
	// (the non-changing nodes) atoms.
	alt_state_rotamer_dots_.increment_both_and_cache( neighbors_alternate_state, this_overlap_with_neighbors_alternate,
		neighbors_alternate_overlap_with_this, atom_atom_overlaps_cache );

	update_alt_state_exphphobes();
	alt_state_inv_dots_.setup_from_rotamer_dots( alt_state_rotamer_dots_, alt_state_exp_hphobes_ );

#ifdef FILE_DEBUG
	//TR_NODE << "update_state_for_neighbors_substitution(): neighbors_alternate_overlap_on_this: " << std::endl;
	//neighbors_alternate_overlap_with_this.print( std::cout );
#endif

#ifdef FILE_DEBUG
	// this is almost certaintly going to print, right?
	if ( current_state_rotamer_dots_ != alt_state_rotamer_dots_ ) {
		TR_NODE << "update_state_for_neighbors_substitution(): node " << parent::get_node_index()
			<< " calc. deltaE for changing node " << node_considering_substitution->get_node_index()
			<< ", curr state sasa: " << current_state_rotamer_dots_.get_sasa() << ", alt state sasa: " << alt_state_rotamer_dots_.get_sasa()
				<< ", sasa difference: " << get_sasa_difference() << std::endl;
	}
#endif

	return get_sasa_difference();
}


///
/// @begin HPatchNode::get_sasa_difference
///
/// @brief
/// Returns alt_state_rotamer_dots_.get_sasa() - curr_state_dots.get_sasa() except when either the current state or the
/// alternate state is the unassigned state; state 0.
///
/// @detailed
/// This method requires that the variables current_state_ and alternate_state_  correspond to the rotamers held in
/// current_state_rotamer_dots_ and alt_state_rotamer_dots_.  Usually, alternate_state_ holds meaningful data only if a vertex is considering
/// a state substitution.  This method will be invoked by a vertex as it considers a state substitition; it will also be
/// invoked by each of its neighbors.  That means the statement alternate_state_ = current_state_; must be present in
/// HPatchNode::update_state_for_neighbors_substitution().
///
template < typename V, typename E, typename G >
Real HPatchNode< V, E, G >::get_sasa_difference() const {

	if ( parent::get_current_state() == 0 && parent::get_alternate_state() == 0 )
		return 0.0;
	else if ( parent::get_current_state() != 0 && parent::get_alternate_state() == 0 )
		return -1 * current_state_rotamer_dots_.get_sasa();
	else if ( parent::get_current_state() == 0 && parent::get_alternate_state() != 0 )
		return alt_state_rotamer_dots_.get_sasa();

	// both current and alt are nonzero
	return alt_state_rotamer_dots_.get_sasa() - current_state_rotamer_dots_.get_sasa();
}


///
/// @begin HPatchNode::reset_alt_state_dots
///
/// @brief
/// Sets the alt state rotamer dots to the current state rotamer dots.  See comments in SIG and commit_considered_substitution
/// for more information about why this method exists.
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::reset_alt_state_dots() {

	if ( alt_state_rotamer_dots_.state_unassigned() && current_state_rotamer_dots_.state_unassigned() ) return;
	// since we're only incrementing and decrementing the alt_state hASA, we need a way to reset it each iteration of consider()
	// so we don't get alt_state SASA/overlaps of nonsensical values which screw up the score.
	//
	//if ( alt_state_rotamer_dots_ != current_state_rotamer_dots_ ) { // alot of time is spent doing RotamerDots comparison
	if ( ! alt_state_dots_matches_current_state_dots_ ) {
		alt_state_rotamer_dots_ = current_state_rotamer_dots_;
		alt_state_inv_dots_ = curr_state_inv_dots_;
		alt_state_exp_hphobes_ = curr_state_exp_hphobes_;
		alt_state_dots_matches_current_state_dots_ = true;
#ifdef FILE_DEBUG
		TR_NODE << "reset_alt_state_dots(): node " << parent::get_node_index() << " alt state dots set equal to current state dots." << std::endl;
#endif
	}

}


///
/// @begin HPatchNode::commit_considered_substitution
///
/// @brief
/// Sets the current state to the alternate state this node was asked to
/// consider.  Copies appropriate score information.  Notifies all of its
/// neighbors that it is going through with the state substitution it had been
/// considering.
///
/// @detailed
/// There's a potential situation with considers() and commits() that needs to be checked for here. It's possible that
/// a consider() call is made which causes a set of Nodes to update their alt states correspondingly.  Since the
/// consider() only gets processed by (or actually the call only goes out to) Nodes which are neighbor graph neighbors,
/// only those Nodes will have their counts updated. Assume that the first consider() is really bad and no commit goes out.
/// If we then consider() another sub, the alt state counts at the previous consider()'s set of nodes are incorrect.  If
/// get_deltaE is called by this consider() on some of those nodes, some of them will have their alt states reset.  But
/// when the commit goes out to ALL Nodes that are neighbors (NOT just the neighbor graph neighbors) it's possible that
/// some of the Nodes will save the wrong alt state count.  One way I think this can be avoided to is check at this node,
/// if the alt state count is different from current, whether the node that originally changed is a neighbor graph neighbor
/// of this node.  If it is, that means the counts changed because of that node.  If it's not, then this Node must have
/// been one of the ones that fell out of sync.
///
/// Oooooh, I just thought of another way.  When the SIG consider() method is called, I could have a reset alt state
/// counts method that will deal with Nodes that are out of sync.  That's more elegant than yet another if statement here!
///
/// I'm not sure the above is really a problem in the case of this IG.  It seems like even if a non-committed consider()
/// call occurs that alters the alt state counts/dots at some set of nodes, when a following sub does get commit'd(), then
/// the nodes that will get the commit message should have all been updated. Some nodes will still have alt_state counts
/// that are weird, but when a consider() call comes back around to them, it should reset the alt_state before doing anything.
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::commit_considered_substitution() {

	assert( parent::considering_alternate_state() );

	if ( parent::get_alternate_state() == 0 ) {
		assign_zero_state();
		return;
	}

	//call base class method
	parent::commit_considered_substitution();

	current_state_rotamer_dots_ = alt_state_rotamer_dots_;
	curr_state_inv_dots_ = alt_state_inv_dots_;
	curr_state_exp_hphobes_ = alt_state_exp_hphobes_;
	alt_state_dots_matches_current_state_dots_ = true;
#ifdef FILE_DEBUG
	TR_NODE << "Committed substitution node " << parent::get_node_index() << std::endl;
	//current_state_rotamer_dots_.print( std::cout );
#endif

	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		get_incident_hpatch_edge(ii)->acknowledge_substitution();
	}
	for (int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii) {
		get_edge_to_hpatch_bg_node( ii )->acknowledge_substitution();
	}

}


///
/// @begin HPatchNode::getMemoryUsageInBytes
///
/// @brief
/// Not implemented, but needs to be!
///
template < typename V, typename E, typename G >
unsigned int HPatchNode< V, E, G >::getMemoryUsageInBytes() const {
	return 0;
}

///
/// @begin HPatchNode::count_static_memory
///
/// @brief
/// Returns the amount of static memory used by this Node object
///
template < typename V, typename E, typename G >
unsigned int HPatchNode< V, E, G >::count_static_memory() const {
	return sizeof ( HPatchNode< V, E, G > );
}

///
/// @begin HPatchNode::count_dynamic_memory
///
/// @brief
/// Returns the amount of dynamic memory used by this Node object
///
template < typename V, typename E, typename G >
unsigned int HPatchNode< V, E, G >::count_dynamic_memory() const {
	unsigned int total_memory = parent::count_dynamic_memory();

	total_memory += rotamers_vector_.size() * sizeof( conformation::ResidueCOP );
	total_memory += self_and_bg_dots_for_states_.size() * sizeof( RotamerDots );

	return total_memory;
}


///
/// @begin HPatchNode::write_dot_kinemage
///
/// @brief
/// useful for debugging
///
//template < typename V, typename E, typename G >
//void HPatchNode< V, E, G >::write_dot_kinemage( std::ofstream & output_kin ) {
//	current_state_rotamer_dots_.write_dot_kinemage( output_kin );
//}


///
/// @begin HPatchNode::print
///
/// @brief
/// useful for debugging - writes information about a node to the tracer
///
template < typename V, typename E, typename G >
void HPatchNode< V, E, G >::print() const {
	TR_NODE << "node " << parent::get_node_index() << ", current_state: " << parent::get_current_state()
			<< ", one body energy: " << parent::get_one_body_energy_current_state() << std::endl;
	current_state_rotamer_dots_.print( std::cout );
	alt_state_rotamer_dots_.print( std::cout );
}


//----------------------------------------------------------------------------//
//----------------------- HPatch Background Residue Node Class -----------------//
//----------------------------------------------------------------------------//

///
/// @begin HPatchBackgroundNode::HPatchBackgroundNode
///
/// @brief
/// main constructor.  No default constructor, copy constructor or assignment operator
///
template < typename V, typename E, typename G >
HPatchBackgroundNode< V, E, G >::HPatchBackgroundNode( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int node_index ) :
	BackgroundNode< V, E, G > ( owner, node_index ),
	prepared_for_simA_( false ),
	current_state_rotamer_dots_(),
	alt_state_rotamer_dots_(),
	alt_state_dots_matches_current_state_dots_( true )
{}


///
/// @begin HPatchBackgroundNode::~HPatchBackgroundNode
///
template < typename V, typename E, typename G >
HPatchBackgroundNode< V, E, G >::~HPatchBackgroundNode() {
	//TR_BGNODE << "called destructor" << std::endl;
}


///
/// @begin HPatchBackgroundNode::set_rotamer
///
/// @brief
/// inits the RotamerOP held by a background node. called in the HPatchIG::initialize method.
///
template < typename V, typename E, typename G >
void HPatchBackgroundNode< V, E, G >::set_rotamer( conformation::ResidueOP const & rotamer ) {
	rotamer_ = rotamer;
	{ // scope -- I should make this a function
	std::string carbon( "C" ), sulfur( "S" );
	n_hphobes_ = 0;
	chemical::ResidueType const & restype = rotamer_->type();
	for ( Size ii = 1; ii <= restype.nheavyatoms(); ++ii ) {
		if ( restype.atom_type( ii ).element() == carbon || restype.atom_type( ii ).element() == sulfur ) {
			++n_hphobes_;
		}
	}
	hphobe_ats_.reserve( n_hphobes_ );
	for ( Size ii = 1; ii <= restype.nheavyatoms(); ++ii ) {
		if ( restype.atom_type( ii ).element() == carbon || restype.atom_type( ii ).element() == sulfur ) {
			hphobe_ats_.push_back( ii );
		}
	}
	} // end scope
	curr_state_exp_hphobes_.reserve( n_hphobes_ );
	alt_state_exp_hphobes_.reserve( n_hphobes_ );


}


///
/// @begin HPatchBackgroundNode::get_rotamer
///
/// @brief
/// returns a const reference to the ResidueOP object held by this background Node. needed by the initialize_bg_bg_overlap
/// call in this class so that one BGNode can get at the residue/rotamer information on another BGNode instance.
///
template < typename V, typename E, typename G >
conformation::ResidueCOP HPatchBackgroundNode< V, E, G >::get_rotamer() const {
	return rotamer_;
}


///
/// @begin HPatchBackgroundNode::set_rotamer_dots
///
/// @brief
/// inits the RotamerDots object held by a background node. called in the HPatchIG::initialize method.
///
template < typename V, typename E, typename G >
void HPatchBackgroundNode< V, E, G >::set_rotamer_dots( RotamerDots const & bg_rd ) {
	current_state_rotamer_dots_ = bg_rd;
	alt_state_dots_matches_current_state_dots_ = false;
}


///
/// @begin HPatchBackgroundNode::detect_overlap
///
/// @brief
/// returns true if this background residue overlaps with any atom on any rotamer of a HPatchNode.
///
/// @detailed
/// Uses the HPatchNode function detect_any_overlap_with_rotamer. It calls the function on the passed in HPatchNode with the
/// RotamerDots object this BGNode keeps in current_state_rotamer_dots_.
///
template < typename V, typename E, typename G >
bool HPatchBackgroundNode< V, E, G >::detect_overlap( HPatchNode< V, E, G >* node ) const {
	return node->detect_any_overlap_with_rotamer( current_state_rotamer_dots_ );
}


///
/// @begin HPatchBackgroundNode::prepare_for_simulated_annealing
///
/// @brief
/// detects self overlap and asks its incident HPatchBackgroundEdges to compute and cache the overlaps of this residue
/// with all states on the neighboring HPatchNode.
///
/// @detailed
/// if we've already gone through this method once, then don't update the parent classes edge vectors. but we do need to
/// reinitialize the self overlap on the BGNode as well as the overlap caused by all other bg nodes.
///
template < typename V, typename E, typename G >
void HPatchBackgroundNode< V, E, G >::prepare_for_simulated_annealing() {

	if ( ! prepared_for_simA_ ) {
		if ( ! parent::get_edge_vector_up_to_date() )
			parent::update_edge_vector();

		prepared_for_simA_ = true;
	}

	// this call is not necessary the first time through packing, but is necessary the 2nd and higher times through because
	// the current state rotamer dots will have overlap from FC nodes mixed in with self and bg node overlap.
 	current_state_rotamer_dots_.zero();
	alt_state_dots_matches_current_state_dots_ = false;

	initialize_self_overlap();
	initialize_atom_atom_overlaps();
//#ifdef FILE_DEBUG
//	TR_BGNODE << "prepare_for_simulated_annealing(): initializing overlap on bg node " << parent::get_node_index() << std::endl;
//#endif
	for ( Size ii = 1; ii <= (Size)parent::get_num_incident_edges(); ++ii ) {
		get_hpatch_bg_edge( ii )->initialize_overlap_cache( current_state_rotamer_dots_ );
	}

}


///
/// @begin HPatchBackgroundNode::initialize_self_overlap
///
/// @brief
/// initializes the self overlap for its RotamerDots object
///
template < typename V, typename E, typename G >
void HPatchBackgroundNode< V, E, G >::initialize_self_overlap() {
//#ifdef FILE_DEBUG
//	if ( parent::get_node_index() == 5 ) {
//		TR_BGNODE << "initialize_self_overlap(): current_state_rotamer_dots_: " << std::endl;
//		current_state_rotamer_dots_.print( std::cout );
//	}
//#endif

	alt_state_dots_matches_current_state_dots_ = false;
	current_state_rotamer_dots_.increment_self_overlap();

//#ifdef FILE_DEBUG
//	if ( parent::get_node_index() == 5 ) {
//		TR_BGNODE << "initialize_self_overlap(): current_state_rotamer_dots_ after increment_self_overlap: " << std::endl;
//		current_state_rotamer_dots_.print( std::cout );
//	}
//#endif

}


///
/// @begin HPatchBackgroundNode::initialize_atom_atom_overlaps
///
/// @brief
/// initializes the atom-atom overlap vector stored by this background node. called once during prep_for_simA.
///
template < typename V, typename E, typename G >
void HPatchBackgroundNode< V, E, G >::initialize_atom_atom_overlaps() {

	self_atom_atom_overlaps_.resize( rotamer_->nheavyatoms(), utility::vector1< bool >( rotamer_->nheavyatoms(), false ) );

	std::string carbon_atom = "C";
	std::string sulfur_atom = "S";
	Real const probe_radius = 1.4;

	for ( Size iia=1; iia <= rotamer_->nheavyatoms(); ++iia ) {
		// immediately continue if not a hydrophobic atom; no point in computing overlaps for polar atoms
		if ( rotamer_->atom_type( iia ).element() != carbon_atom && rotamer_->atom_type( iia ).element() != sulfur_atom )
			continue;

		Real const iia_atom_radius = current_state_rotamer_dots_.get_atom_radius( iia ) + probe_radius;

		// for intra-residue we only have iterate over the greater-indexed heavyatoms
		for ( Size jja = iia + 1; jja <= rotamer_->nheavyatoms(); ++jja ) {
			if ( rotamer_->atom_type( jja ).element() != carbon_atom && rotamer_->atom_type( jja ).element() != sulfur_atom )
				continue;

			Real const jja_atom_radius = current_state_rotamer_dots_.get_atom_radius( jja ) + probe_radius;

			// check if the two atoms overlap, and if that degree of overlap exceeds the threshold
			Vector const & iia_atom_xyz = rotamer_->atom( iia ).xyz();
			Vector const & jja_atom_xyz = rotamer_->atom( jja ).xyz();

			Real const distance_squared = iia_atom_xyz.distance_squared( jja_atom_xyz );
			//if ( distance_squared <= (iia_atom_radius + iia_atom_radius) * (jja_atom_radius + jja_atom_radius) ) {
			if ( distance_squared <= (iia_atom_radius + jja_atom_radius) * (iia_atom_radius + jja_atom_radius) ) {

				Real const distance_ijxyz = std::sqrt( distance_squared );
				int degree_of_overlap;
				core::scoring::get_overlap( iia_atom_radius, jja_atom_radius, distance_ijxyz, degree_of_overlap );
				if ( degree_of_overlap >= 15 ) {
#ifdef FILE_DEBUG
					//TR_BGNODE << "initialize_self_overlap(): overlapping intra-residue atom pair: "
					//	<< rotamer_->seqpos() << "/" << utility::trim( rotamer_->atom_name( iia ) ) << " - " << rotamer_->seqpos() << "/" << utility::trim( rotamer_->atom_name( jja ) )
					//	<< ", degree of overlap: " << degree_of_overlap << std::endl;
#endif
					self_atom_atom_overlaps_[ iia ][ jja ] = true;
				}
			} // end if distance
		} // for loop over all jj heavyatoms
	} // for loop over all ii heavyatoms


/*
#ifdef FILE_DEBUG
	TR_BGNODE << "background node " << parent::get_node_index() << " self_atom_atom_overlaps_: [ " << std::endl;
	for ( Size aa=1; aa <= self_atom_atom_overlaps_.size(); ++aa ) {
		TR_BGNODE << "self_atom_atom_overlaps_[ " << aa << " ]: [ ";
		for ( Size bb=1; bb <= self_atom_atom_overlaps_[ aa ].size(); ++bb ) {
			TR_BGNODE << self_atom_atom_overlaps_[ aa ][ bb ] << ", ";
		}
		TR_BGNODE << "]" << std::endl;
	}
	TR_BGNODE << std::endl;
#endif
*/

}


///
/// @begin HPatchBackgroundNode::initialize_bg_bg_overlap
///
/// @brief
/// stores the sphere overlap for a pair of background nodes
///
template < typename V, typename E, typename G >
void HPatchBackgroundNode< V, E, G >::initialize_bg_bg_overlap( HPatchBackgroundNode< V, E, G > & other ) {
	current_state_rotamer_dots_.increment_both( other.current_state_rotamer_dots_ );
	alt_state_dots_matches_current_state_dots_ = false;

//#ifdef FILE_DEBUG
//	if ( parent::get_node_index() == 5 ) {
//		TR_BGNODE << "initialize_bg_bg_overlap(): current_state_rotamer_dots_ after adding in others overlap: " << std::endl;
//		current_state_rotamer_dots_.print( std::cout );
//	}
//#endif
}


///
/// @begin HPatchBackgroundNode::update_state_for_substitution
///
/// @brief
/// returns the change in sasa induced by a HPatchNode undergoing a state substitution. The overlap between the HPatchNode
/// and this HPatchBackgroundNode has been precomputed and stored in the HPatchBackgroundEdge that connects them.
/// There is little work to do in this subroutine.
///
template < typename V, typename E, typename G >
Real HPatchBackgroundNode< V, E, G >::update_state_for_substitution( HPatchNode< V, E, G >* fc_node_changing,
	RotamerDotsCache const & nodes_curr_overlap_with_bg_res, RotamerDotsCache const & nodes_alt_overlap_with_bg_res ) {

	get_hpatch_owner()->register_bg_node_affected_by_rotsub( parent::get_node_index() );

	fc_node_changing->wt_seqpos_for_node(); // do something with fc_node_changing or otherwise we get unused var warning

	//alt_state_rotamer_dots_.copy_spheres_not_doubly_covered( &current_state_rotamer_dots_ );
	//alt_state_rotamer_dots_ = current_state_rotamer_dots_;
	if ( ! alt_state_dots_matches_current_state_dots_ ) {
		alt_state_rotamer_dots_ = current_state_rotamer_dots_;
	}
	alt_state_dots_matches_current_state_dots_ = false;

#ifdef FILE_DEBUG
	//TR_BGNODE << "update_state_for_substitution(): alt_state_rotamer_dots_" << std::endl;
	//alt_state_rotamer_dots_.print( std::cout );
	//TR_BGNODE << "update_state_for_substitution(): nodes_curr_overlap_with_bg_res" << std::endl;
	//nodes_curr_overlap_with_bg_res.print( std::cout );
#endif

	alt_state_rotamer_dots_.decrement_from_cached( nodes_curr_overlap_with_bg_res );

#ifdef FILE_DEBUG
	//TR_BGNODE << "update_state_for_substitution(): nodes_alt_overlap_with_bg_res" << std::endl;
	//nodes_alt_overlap_with_bg_res.print( std::cout );
#endif

	alt_state_rotamer_dots_.increment_from_cached( nodes_alt_overlap_with_bg_res );

	update_alt_state_exphphobes();
	alt_state_inv_dots_.setup_from_rotamer_dots( alt_state_rotamer_dots_, alt_state_exp_hphobes_ );
#ifdef FILE_DEBUG
	//TR_BGNODE << "update_state_for_substitution(): alt_state_rotamer_dots_" << std::endl;
	//alt_state_rotamer_dots_.print( std::cout );
#endif

	Real sasa_difference = alt_state_rotamer_dots_.get_sasa() - current_state_rotamer_dots_.get_sasa();

#ifdef FILE_DEBUG
	Size fc_node_resid = fc_node_changing->wt_seqpos_for_node();
	TR_BGNODE << "update_state_for_substitution: bg node " << parent::get_node_index()
		<< " calculated deltaE for changing hpatch resid " << fc_node_resid << "; "
		<< "curr state sasa: " << current_state_rotamer_dots_.get_sasa() << ", alt state sasa: " << alt_state_rotamer_dots_.get_sasa()
		<< ", sasa difference: " << sasa_difference << std::endl;
#endif

	return sasa_difference;
}


///
/// @begin HPatchBackgroundNode::reset_alt_state_dots
///
/// @brief
/// Sets the alt state dots to the current state dots.  See comments in HIG and commit_considered_substitution
/// for more information about why this method exists.
///
template < typename V, typename E, typename G >
void HPatchBackgroundNode< V, E, G >::reset_alt_state_dots() {

	//if ( alt_state_rotamer_dots_ != current_state_rotamer_dots_ ) { // alot of time spent doing comparison. just set them equal irregardless.
	if ( ! alt_state_dots_matches_current_state_dots_ ) {
		alt_state_rotamer_dots_ = current_state_rotamer_dots_;
		alt_state_inv_dots_ = curr_state_inv_dots_;
		alt_state_exp_hphobes_ = curr_state_exp_hphobes_;
		alt_state_dots_matches_current_state_dots_ = true;
#ifdef FILE_DEBUG
		TR_BGNODE << "reset_alt_state_dots(): node " << parent::get_node_index() << " alt state rotamer dots set equal to current state rotamer dots." << std::endl;
#endif
	}

}


///
/// @begin HPatchBackgroundNode::acknowledge_substitution
///
/// @brief
/// bookkeeping to reflect a HPatchNode's state substitution. uses the RotamerDots class method operator=.
///
template < typename V, typename E, typename G >
void HPatchBackgroundNode< V, E, G >::acknowledge_substitution() {
	current_state_rotamer_dots_ = alt_state_rotamer_dots_;
	curr_state_inv_dots_ = alt_state_inv_dots_;
	curr_state_exp_hphobes_ = alt_state_exp_hphobes_;
	alt_state_dots_matches_current_state_dots_ = true;
}


///
/// @begin HPatchBackgroundNode::update_alt_state_exphphobes
///
/// @brief
/// Updates the vector alt_state_exp_hphobes_ with the atom id (??) of the exposed hydrophobic atoms in this residue.
///
template < typename V, typename E, typename G >
void HPatchBackgroundNode< V, E, G >::update_alt_state_exphphobes() {
	alt_state_exp_hphobes_.clear();
	for ( Size ii = 1; ii <= hphobe_ats_.size(); ++ii ) {
		Size const ii_atom = hphobe_ats_[ ii ];
		if ( alt_state_rotamer_dots_.get_atom_sasa( ii_atom ) > 0.0 ) {
			alt_state_exp_hphobes_.push_back( ii_atom );
		}
	}
}


///
/// @begin HPatchBackgroundNode< V, E, G >::get_current_sasa
///
/// @brief
/// returns the total SASA under the current state assignment
///
template < typename V, typename E, typename G >
Real HPatchBackgroundNode< V, E, G >::get_current_sasa() const {
	return current_state_rotamer_dots_.get_sasa();
}


///
/// @begin HPatchBackgroundNode< V, E, G >::get_current_sasa
///
/// @brief
/// Returns the current state SASA for the passed in atom index
///
template < typename V, typename E, typename G >
Real
HPatchBackgroundNode< V, E, G >::get_current_sasa( Size atom_index ) const {
	return current_state_rotamer_dots_.get_atom_sasa( atom_index );
}


///
/// @begin HPatchBackgroundNode< V, E, G >::get_alternate_sasa
///
/// @brief
/// returns the total SASA under the alternate state assignment
///
template < typename V, typename E, typename G >
Real
HPatchBackgroundNode< V, E, G >::get_alternate_sasa() const {
	return alt_state_rotamer_dots_.get_sasa();
}


///
/// @begin HPatchBackgroundNode< V, E, G >::get_alternate_sasa
///
/// @brief
/// Returns the alternate state SASA for the passed in atom index
///
template < typename V, typename E, typename G >
Real
HPatchBackgroundNode< V, E, G >::get_alternate_sasa( Size atom_index ) const {
	return alt_state_rotamer_dots_.get_atom_sasa( atom_index );
}


///
/// @begin HPatchBackgroundNode< V, E, G >::get_atom_atom_self_overlaps
///
/// @brief
/// Returns a const reference to the atom x atom pair vector of vectors of bools that specifies which atoms are overlapping.
///
template < typename V, typename E, typename G >
utility::vector1< utility::vector1 < bool > > const & HPatchBackgroundNode< V, E, G >::get_atom_atom_self_overlaps() const {
	return self_atom_atom_overlaps_;
}


///
/// @begin HPatchBackgroundNode::count_static_memory
///
template < typename V, typename E, typename G >
unsigned int HPatchBackgroundNode< V, E, G >::count_static_memory() const {
	return sizeof ( HPatchBackgroundNode< V, E, G > );
}

///
/// @begin HPatchBackgroundNode::count_dynamic_memory
///
template < typename V, typename E, typename G >
unsigned int HPatchBackgroundNode< V, E, G >::count_dynamic_memory() const {
	unsigned int total_memory = parent::count_dynamic_memory();

	return total_memory;
}

///
/// @begin HPatchBackgroundNode::print
///
/// @brief
/// used only for debugging
//template < typename V, typename E, typename G >
//void HPatchBackgroundNode< V, E, G >::write_dot_kinemage( std::ofstream & output_kin ) {
//	current_state_rotamer_dots_.write_dot_kinemage( output_kin );
//}


///
/// @begin HPatchBackgroundNode::print
///
/// @brief
/// used only for debugging
template < typename V, typename E, typename G >
void HPatchBackgroundNode< V, E, G >::print() const {
	TR_BGNODE << "bgnode " << parent::get_node_index() << ", current state sasa: " << current_state_rotamer_dots_.get_sasa()
		<< ", alt state sasa: " << alt_state_rotamer_dots_.get_sasa() << std::endl;
	current_state_rotamer_dots_.print( std::cout );
}


//----------------------------------------------------------------------------//
//------------------------------ HPatch Edge Class -----------------------------//
//----------------------------------------------------------------------------//


///
/// @begin HPatchEdge::HPatchEdge
///
/// @brief
/// main constructor.  No default, or copy constructors, no assignment operator
///
/// @param
/// owner - [in] - the owning interaction graph object
/// node1 - [in] - the index of the lower-indexed HPatchNode
/// node2 - [in] - the index of the higher-indexed HPatchNode
///
template < typename V, typename E, typename G >
HPatchEdge< V, E, G >::HPatchEdge( G* owner, int node1, int node2 ) :
	FirstClassEdge< V, E, G > ( owner, node1, node2 ),
	node_changing_( -1 ),
	node_not_changing_( -1 )
	//nodes_curr_pair_dot_counts_ // calls default RDC constructor, which is fine
	//nodes_alt_pair_dot_counts_ // calls default RDC constructor, which is fine
	//current_state_atom_atom_overlaps_ // vectors are init'd to nothing
	//alt_state_atom_atom_overlaps_;
{
	nodes_curr_states_[ 0 ] = nodes_curr_states_[ 1 ] = 0;
	nodes_alt_states_[ 0 ] = nodes_alt_states_[ 1 ] = 0;

}


///
/// @begin HPatchEdge::~HPatchEdge
///
template < typename V, typename E, typename G >
HPatchEdge< V, E, G >::~HPatchEdge() {
	//TR_EDGE << "called destructor" << std::endl;
}


///
/// @begin HPatchEdge::prepare_for_simulated_annealing
///
/// @brief
/// drops zero submatrices of the AminoAcidNeighborSparseMatrix and if the two_body_energies_ member then holds nothing,
/// it checks whether or not its incident nodes have any sphere overlaps.  If they don't then the edge deletes itself.
///
template < typename V, typename E, typename G >
void HPatchEdge< V, E, G >::prepare_for_simulated_annealing() {

#ifdef FILE_DEBUG
	TR_EDGE << "prepare_for_simulated_annealing called on edge e(" << parent::get_node_index( 0 ) << "," << parent::get_node_index( 1 ) << ")" << std::endl;
#endif

	parent::prepare_for_simulated_annealing_no_deletion(); // AddtlBGNodesIG doesn't have an Edge method for prep for simA, but PDIG does

	if ( parent::pd_edge_table_all_zeros() ) {
		if ( ! (get_hpatch_node(0)->overlaps( get_hpatch_node(1)) ) ) {
#ifdef FILE_DEBUG
			TR_EDGE << "prepare_for_simulated_annealing - dropping edge e(" << parent::get_node_index( 0 )
				<< "," << parent::get_node_index( 1 ) << ")" << std::endl;
#endif
			delete this;
		}
	}
}


///
/// @begin HPatchEdge::acknowledge_state_zeroed
///
/// @brief
/// respond to when one of its vertices enters the "unassigned" state.
///
/// @detailed
/// called during the HIG::blanket_assign_state_0 -> HPatchNode::assign_zero_state() cascade of calls.
///
template < typename V, typename E, typename G >
void HPatchEdge< V, E, G >::acknowledge_state_zeroed( int node_that_changed ) {

	// each Edge contains a node_changing_ and node_not_changing_ int (which takes on values of 0 or 1)
	// node_changing_ is the node of this edge that is changing
	node_changing_ = ( node_that_changed == parent::get_node_index(0) ? 0 : 1 );
	node_not_changing_ = ! node_changing_;

	nodes_curr_states_[ node_changing_ ] = 0;

	nodes_curr_pair_dot_counts_[ node_changing_ ].resize( 0 );
	nodes_curr_pair_dot_counts_[ node_not_changing_ ].resize( 0 );

	nodes_alt_pair_dot_counts_[ node_changing_ ].resize( 0 );
	nodes_alt_pair_dot_counts_[ node_not_changing_ ].resize( 0 );

	// clear the atom-atom overlap caches; not sure if this is necessary since at construct time these vectors will be
	// empty. the unit tests will call blanket_assign_state_0 on the same graph multiple times though, so clear() calls
	// here will ensure nothing remains.
	// apl -- removing since it breaks monotonicity assumption -- current_state_atom_atom_overlaps_.clear();
	// apl -- removing since it breaks monotonicity assumption -- alt_state_atom_atom_overlaps_.clear();

	inform_non_changing_node_of_neighbors_change();
}


///
/// @begin HPatchEdge::inform_non_changing_node_of_neighbors_change
///
/// @brief
/// tells the node that isn't considering a substitution or changing state that its neighbor who is has changed.
///
template < typename V, typename E, typename G >
inline
void HPatchEdge< V, E, G >::inform_non_changing_node_of_neighbors_change() {
	get_hpatch_node( node_not_changing_ )->acknowledge_neighbors_substitution();
}


///
/// @begin HPatchEdge::update_state_at_neighbor
///
/// @brief
/// returns the change in sasa for the neighbor of a node that is produced by the state substitution it is considering.
///
/// @detailed
/// Very complicated.  See HPatchNode::project_deltaE_for_neighbors_state_sub
/// This edge collects cached sphere overlaps and hands this cached data to
/// the node that is not considering the state substitution.  That node
/// computes more sphere overlaps and returns a delta sasa; this method then
/// passes that delta sasa along.
///
/// See more comments inline.
///
/// @param
/// changing_node_alt_state_dots - [in] - the RotamerDots object for the alternate state at the changing Node
///
template < typename V, typename E, typename G >
Real HPatchEdge< V, E, G >::update_state_at_neighbor( int node_considering_substitution, int alt_state, RotamerDots & changing_node_alt_state_dots ) {

	using namespace utility;

	node_changing_ = ( node_considering_substitution == parent::get_node_index(0) ? 0 : 1 );
	node_not_changing_ = ! node_changing_;

	nodes_alt_states_[ node_changing_ ] = alt_state;
	nodes_alt_states_[ node_not_changing_ ] = nodes_curr_states_[ node_not_changing_ ];

	// zero out the RotamerDotsCache object stored for the non_changing FC node
	nodes_alt_pair_dot_counts_[ node_not_changing_ ].zero();

	// the alternate state RDC object for the changing node should be the curr state dot counts, up to the number of atoms that are the same
	// as the current state?  but since we're not using trie ordering, the number of atoms same as current will be zero.
	nodes_alt_pair_dot_counts_[ node_changing_ ] = nodes_curr_pair_dot_counts_[ node_changing_ ];
	// nodes_alt_pair[ changing ] now contains the RDC object which represents the overlap b/t the two FC Nodes on this Edge in the
	// context of the current state!
	nodes_alt_pair_dot_counts_[ node_changing_ ].zero();

	// nodes_curr_pair_dot_counts_[ node_not_changing_ ] is the overlap the non-changing node experiences because of the changing Node
	// in its current state. (aka neighbors_curr_state_overlap_with_this) - this gets subtracted

	// nodes_alt_pair_dot_counts_[ node_changing_ ] is the overlap the non_changing node produces on the changing Node in its alt state
	// (from perspective of non-changing node: this_overlap_with_neighbors_alternate)
	// this, in turn, gets passed to increment_both_and_cache:
	// other_dots_covered_by_this - the dot coverage cache represnting dots on the surface of other_rotamer that
	// are covered by this rotamer

	// nodes_alt_pair_dot_counts_[ node_not_changing ] is the overlap the changing node in its alt state produces on the non_changing FC
	// node. (from perspective of nonchanging node: neighbors_alternate_overlap_with_this, or this_dots_covered_by_other)

	// have to invalidate/resize the caches stored on the Edge since they will both change
	nodes_alt_pair_dot_counts_[ node_changing_ ].resize( changing_node_alt_state_dots.get_num_atoms() );
	//if ( get_hpatch_node( node_not_changing_ )->get_current_state() != 0 )
		nodes_alt_pair_dot_counts_[ node_not_changing_ ].resize( get_hpatch_node( node_not_changing_ )->get_current_state_num_atoms() );
	// the RDC object for the node_not_changing may get resized to 0, but it doesn't matter because the method update_state_
	// for_neighbors_sub will immediately return if the current state on the non_changing node is 0.

	// also have to update the atom-atom overlap cache for the considered sub.
	// we have to resize the vector of vectors here, but we have to be consistent with which one comes first. let's always make
	// the fc node that's changing the outer vector and the non-changing node the inner vector.
	// the values in this vector will get set by the non-changing node when it's updating the sasa information. no point in
	// doing the sphere overlap calculations twice.
	Size const node_not_changing_num_atoms = get_hpatch_node( node_not_changing_ )->get_current_state_num_atoms();
	Size const node_changing_alt_state_num_atoms = changing_node_alt_state_dots.get_num_atoms();

	//alt_state_atom_atom_overlaps_.resize( node_changing_alt_state_num_atoms, utility::vector1< bool >( node_not_changing_num_atoms, false ) );
	/// APL -- don't resize these arrays; too expensive; just leave it large
	if ( alt_state_atom_atom_overlaps_.size() < node_changing_alt_state_num_atoms ) {
		alt_state_atom_atom_overlaps_.resize( node_changing_alt_state_num_atoms, utility::vector1< bool >( node_not_changing_num_atoms, false ) );
	}

	for ( Size ii = 1; ii <= node_changing_alt_state_num_atoms; ++ii ) {
		if ( alt_state_atom_atom_overlaps_[ii].size() < node_not_changing_num_atoms ) {
			alt_state_atom_atom_overlaps_[ii].resize( node_not_changing_num_atoms, false );
		}
		for ( Size jj = 1; jj <= node_not_changing_num_atoms; ++jj ) {
			alt_state_atom_atom_overlaps_[ ii ][ jj ] = false;
		}
	}

	Real const delta_sasa = get_hpatch_node( node_not_changing_ )->update_state_for_neighbors_substitution(
		get_hpatch_node( node_changing_ ), changing_node_alt_state_dots,
		nodes_curr_pair_dot_counts_[ node_not_changing_ ],
		nodes_alt_pair_dot_counts_[ node_changing_ ],
		nodes_alt_pair_dot_counts_[ node_not_changing_ ],
		alt_state_atom_atom_overlaps_ );

	// in get_hpatch_deltaE_for_nbs_state_sub(), non-changing node calls...
	// alt_state_rotamer_dots_.increment_both_and_cache( neighbors_alternate_state, this_overlap_with_neighbors_alternate,
	//	neighbors_alternate_overlap_with_this, mask_this_covered_by_other, mask_other_covered_by_this );

#ifdef FILE_DEBUG
	// this line just generates waaaay too much output...
	//TR_EDGE << "update_state_at_neighbor() for changing node: " << node_considering_substitution
	//	<< " and node: " << parent::get_node_index( node_not_changing_ ) << " returning delta sasa of " << delta_sasa << std::endl;
#endif

	return delta_sasa;
}


///
/// @begin HPatchEdge::acknowledge_substitution
///
/// @brief
/// bookkeeping following the decision to substitute a nodes current state with the alternate it was asked to consider.
///
template < typename V, typename E, typename G >
void HPatchEdge< V, E, G >::acknowledge_substitution() {

	inform_non_changing_node_of_neighbors_change();

	nodes_curr_states_[ node_changing_ ] = nodes_alt_states_[ node_changing_ ];

	nodes_curr_pair_dot_counts_[0] = nodes_alt_pair_dot_counts_[0];
	nodes_curr_pair_dot_counts_[1] = nodes_alt_pair_dot_counts_[1];

	//current_state_atom_atom_overlaps_ = alt_state_atom_atom_overlaps_;
	if ( node_changing_ == 0 ) {
		if ( current_state_atom_atom_overlaps_.size() < alt_state_atom_atom_overlaps_.size() ) {
			current_state_atom_atom_overlaps_.resize( alt_state_atom_atom_overlaps_.size() );
		}
		for ( Size ii = 1; ii <= alt_state_atom_atom_overlaps_.size(); ++ii ) {
			if ( current_state_atom_atom_overlaps_[ii].size() < alt_state_atom_atom_overlaps_[ii].size() ) {
				current_state_atom_atom_overlaps_[ii].resize( alt_state_atom_atom_overlaps_[ii].size() );
			}
			current_state_atom_atom_overlaps_[ii] = alt_state_atom_atom_overlaps_[ii];
		}
	} else {

		/// Preserve the convention that the current_state_atom_atom_overlaps_ dimensions 1st by node 0
		/// and 2nd by node 1.  This means that if we just considered a substitution at node 1, that the
		/// alt_state_atom_atom_overlaps_ will be dimensioned 1st by node 1 and 2nd by node 0, so we need
		/// to copy curr[ii][jj] from alt[jj][ii];
		if ( alt_state_atom_atom_overlaps_.size() == 0 ) {
			current_state_atom_atom_overlaps_ = alt_state_atom_atom_overlaps_;
			return;
		}

		// debug
		assert( get_hpatch_node( node_changing_ )->get_alt_state_num_atoms() <= alt_state_atom_atom_overlaps_.size() );
		assert( get_hpatch_node( node_not_changing_ )->get_current_state_num_atoms() <= alt_state_atom_atom_overlaps_[1].size() );

		if ( current_state_atom_atom_overlaps_.size() < alt_state_atom_atom_overlaps_[1].size() ) {
			current_state_atom_atom_overlaps_.resize( alt_state_atom_atom_overlaps_[1].size() );
		}

		for ( Size ii=1; ii <= alt_state_atom_atom_overlaps_[1].size(); ++ii ) {
			if ( current_state_atom_atom_overlaps_[ii].size() < alt_state_atom_atom_overlaps_.size() ){
				current_state_atom_atom_overlaps_[ii].resize( alt_state_atom_atom_overlaps_.size() );
			}
		}

		for ( Size ii=1; ii <= alt_state_atom_atom_overlaps_.size(); ++ii ) {
			/// strict monotone growth
			assert( alt_state_atom_atom_overlaps_[ ii ].size() <= alt_state_atom_atom_overlaps_[ 1 ].size() );
			assert( ii > get_hpatch_node( node_changing_ )->get_alt_state_num_atoms()  ||
				get_hpatch_node( node_not_changing_ )->get_current_state_num_atoms() <= alt_state_atom_atom_overlaps_[ii].size() );
			for ( Size jj = 1; jj <= alt_state_atom_atom_overlaps_[ ii ].size(); ++jj ) {
				current_state_atom_atom_overlaps_[ jj ][ ii ] = alt_state_atom_atom_overlaps_[ ii ][ jj ];
			}
		}
	}

	return;
}


///
/// @begin HPatchEdge< V, E, G >::get_current_state_atom_atom_overlaps
///
/// @brief
/// Returns a const reference to the atom-x-atom-pair vector-of-vectors of bools that specifies which atoms are overlapping,
/// assuming the current state assignment.
///
template < typename V, typename E, typename G >
utility::vector1< utility::vector1 < bool > > const &
HPatchEdge< V, E, G >::get_current_state_atom_atom_overlaps() const {
	return current_state_atom_atom_overlaps_;
}


///
/// @begin HPatchEdge< V, E, G >::get_alt_state_atom_atom_overlaps
///
/// @brief
/// Returns a const reference to the atom-x-atom-pair vector-of-vectors of bools that specifies which atoms are overlapping,
/// assuming the alternate state assignment.
///
template < typename V, typename E, typename G >
utility::vector1< utility::vector1 < bool > > const &
HPatchEdge< V, E, G >::get_alt_state_atom_atom_overlaps() const {
	return alt_state_atom_atom_overlaps_;
}


///
/// @begin HPatchEdge::declare_energies_final
///
/// @brief
/// Reduces memory usage in the two body energy table after the energy
/// calculating function declares that the energies will not change thereafter
///
/// @remarks (all by apl)
/// In the PDEdge's version of this method, after invoking two_body_energies_.drop_zero_submatrices_where_possible();
/// the PDEdge checks if the two body energy table it now holds is empty.  If the table is empty, the edge deletes itself.
///
/// A HPatchEdge should not delete itself if the pair energies are all zero since the Minkowski sum of a water and a van
/// der Waal's sphere extends further out from an atoms center than its (lj_atr, lj_rep, lksolv) interaction sphere.
/// However, if a HPatchEdge holds no pair energies, it's a very good candidate for removal  -- it just first needs to check
/// that no (vdw + 1.4 A) spheres overlap between any pair of rotamers on the edges it connects.
///
template < typename V, typename E, typename G >
void HPatchEdge< V, E, G >::declare_energies_final() {
	parent::declare_energies_final_no_deletion();
}


///
/// @begin HPatchEdge::getMemoryUsageInBytes
///
/// @remarks
/// Not implemented.
///
template < typename V, typename E, typename G >
unsigned int HPatchEdge< V, E, G >::getMemoryUsageInBytes() const {
	return 0;
}


///
/// @begin HPatchEdge::count_static_memory
///
template < typename V, typename E, typename G >
unsigned int HPatchEdge< V, E, G >::count_static_memory() const {
	return sizeof ( HPatchEdge< V, E, G > );
}

///
/// @begin HPatchEdge::count_dynamic_memory
///
template < typename V, typename E, typename G >
unsigned int HPatchEdge< V, E, G >::count_dynamic_memory() const {

	unsigned int total_memory = parent::count_dynamic_memory();

	//assert( alt_state_atom_atom_overlaps_.size() == current_state_atom_atom_overlaps_.size() );
	total_memory += sizeof( utility::vector1< bool > ) * alt_state_atom_atom_overlaps_.size();
	for ( Size ii = 1; ii <= alt_state_atom_atom_overlaps_.size(); ++ii ) {
		total_memory += sizeof( bool ) * alt_state_atom_atom_overlaps_[ii].size();
	}
	total_memory += sizeof( utility::vector1< bool > ) * current_state_atom_atom_overlaps_.size();
	for ( Size ii = 1; ii <= current_state_atom_atom_overlaps_.size(); ++ii ) {
		total_memory += sizeof( bool ) * current_state_atom_atom_overlaps_[ii].size();
	}

	return total_memory;
}



//----------------------------------------------------------------------------//
//----------------------- HPatchBackgroundEdge Class ---------------------------//
//----------------------------------------------------------------------------//

///
/// @begin HPatchBackgroundEdge< V, E, G >::HPatchBackgroundEdge
///
/// @brief
/// main constructor
///
template < typename V, typename E, typename G >
HPatchBackgroundEdge< V, E, G >::HPatchBackgroundEdge( AdditionalBackgroundNodesInteractionGraph < V, E, G >* owner, int first_class_node_index, int background_node_index ) :
	BackgroundToFirstClassEdge< V, E, G >( owner, first_class_node_index, background_node_index ),
	prepared_for_simA_( false ),
	node_states_coverage_of_bg_res_(),
	nodes_curr_state_( 0 ),
	nodes_alt_state_( 0 ),
	node_states_overlap_with_bg_res_()
	//curr_dots_cache_( 0 ), // default constructed, which means they need to be sized before use
	//alt_dots_cache_( 0 ) // default constructed
{}


///
/// @begin HPatchBackgroundEdge::~HPatchBackgroundEdge
///
template < typename V, typename E, typename G >
HPatchBackgroundEdge< V, E, G >::~HPatchBackgroundEdge() {
	//TR_EDGE << "called destructor" << std::endl;
	// node_states_coverage_of_bg_res_ is a vector, so when it goes out of scope the memory will be freed
}


///
/// @begin HPatchBackgroundEdge::prepare_for_simulated_annealing
///
/// @brief
/// Invoked by AdditionalBackgroundNodesInteractionGraph::prepare_for_simulated_annealing.
///
/// @remarks
/// The HPatchBackgroundEdge has no responsibilities in this function.  However, when the AdditionalBackgroundNodesInteractionGraph invokes
/// prepare_for_simulated_annealing on the HPatchBackgroundNode that this edge is incident upon, that node will invoke initialize_overlap_cache
/// on this edge.
///
template < typename V, typename E, typename G >
void HPatchBackgroundEdge< V, E, G >::prepare_for_simulated_annealing() {}


///
/// @begin HPatchBackgroundEdge::initialize_overlap_cache
///
/// @brief
/// compute the sphere overlaps of the background node with every state on the first class node. The HPatchBackgroundEdge
/// hands its stl vector of RotamerDotsCache objects (node_states_coverage_of_bg_res) to the HPatchNode
///
/// Called during the prep for simA method in HPatchBGNodes.
/// This method in turn calls a HPatchNode method, init_overlap_with_background to set the vector of RDC object pointers.
///
template < typename V, typename E, typename G >
void HPatchBackgroundEdge< V, E, G >::initialize_overlap_cache( RotamerDots const & bg_residue_dots ) {

	//using namespace utility;

	bg_res_num_atoms_ = bg_residue_dots.get_num_atoms();

	if ( ! prepared_for_simA_ ) {
		// init every element of this std::vector with an instance of a RDC of size bg res num atoms
		// an alternative way to do this would be to make this a vector of pointers, and allocate only the pointers for
		// which there is nonzero overlap. but then we'd have to be extra careful about freeing the memory in the destructor.
		// if I find that too much memory is being used by this, I'll change it later.
		node_states_coverage_of_bg_res_.resize( (Size)(get_hpatch_node()->get_num_states() + 1), RotamerDotsCache( bg_res_num_atoms_ ) );

		// resize the atom-atom overlaps vector. but we can't init the values to anything because we don't know how many atoms
		// are in each state on the first class node. just have to init the values in the node method.
		node_states_overlap_with_bg_res_.resize( (Size)(get_hpatch_node()->get_num_states() + 1) );

		get_hpatch_node()->initialize_overlap_with_background( bg_residue_dots, node_states_coverage_of_bg_res_, node_states_overlap_with_bg_res_ );

		prepared_for_simA_ = true;
	}

	// since we now know how many atoms the bg residue has, resize the two member RDC objects to that number of atoms
	curr_dots_cache_.resize( bg_res_num_atoms_ );
	alt_dots_cache_.resize( bg_res_num_atoms_ );

}


///
/// @begin HPatchBackgroundEdge::acknowledge_state_change
///
/// @brief
/// bookkeeping in response to a HPatchNode switching states (without having gone through the usual
/// consider-substitution/commit-substitution pattern).
///
/// @detailed
/// State "0" is handled by the HPatchBackgroundEdge just like any other state. The dot overlap cache is simply empty:
/// the unassigned state induces no overlap on the background node.  The HPatchBackgroundEdge keeps a stl vector of
/// RotamerDotCaches.  Position 0 in that vector holds a Cache object of nothing but 0's.
///
template < typename V, typename E, typename G >
void HPatchBackgroundEdge< V, E, G >::acknowledge_state_change( int new_state ) {

	if ( new_state == nodes_curr_state_ ) // in the case of the 0-state, just return - don't calculate deltaE
		return;

	update_state_at_neighbor( new_state );
	acknowledge_substitution();

}


///
/// @begin HPatchBackgroundEdge::update_state_at_neighbor
///
/// @brief
/// returns the change in hpatch energy produced by a background node in response to a considered state substitution of
/// the first class node
///
template < typename V, typename E, typename G >
Real HPatchBackgroundEdge< V, E, G >::update_state_at_neighbor( int alt_state ) {

#ifdef FILE_DEBUG
	//TR_BGEDGE << "update_state_at_neighbor() for changing fc node: " << parent::get_first_class_node_index()
	//	<< " and bg node: " << parent::get_background_node_index() << std::endl;
#endif

	// node_states_coverage_of_bg_res is a vector that contains whether or not a given state on a FCNode (that this edge
	// is connected to) overlaps with this bg residue. so if the value at the alt_state state index is nonzero, then there's
	// some overlap with the alt_state on this bg residue. in that case, we need to set the non_zero_overlap array and
	// then do some updating. (in the case of zero state being assigned, that will point at index 0 which has nothing in
	// it) resulting in nothing being done.
	alt_dots_cache_.zero();
	alt_dots_cache_ = node_states_coverage_of_bg_res_[ alt_state ];

	nodes_alt_state_ = alt_state;

	// for when blanket assign zero state gets called; it doesn't really matter what we return since it's an unassigned state
	// ah, but we do because later this 0-state sub gets "committed" and current_state gets set to alt_state. if we don't
	// first "consider" the 0-state sub in update_state_for_sub() then we reinit incorrectly.
	//if ( alt_state == 0 ) return 0.0;

	Real const delta_sasa = get_hpatch_bg_node()->update_state_for_substitution( get_hpatch_node(), curr_dots_cache_, alt_dots_cache_ );

	return delta_sasa;
}


///
/// @begin HPatchBackgroundEdge::acknowledge_substitution
///
/// @brief
/// bookkeeping in response to the HPatchNode committing the considered substitution
///
template < typename V, typename E, typename G >
void HPatchBackgroundEdge< V, E, G >::acknowledge_substitution() {

	get_hpatch_bg_node()->acknowledge_substitution();

	nodes_curr_state_ = nodes_alt_state_;
	curr_dots_cache_ = alt_dots_cache_;

}


///
/// @begin HPatchBackgroundEdge< V, E, G >::get_atom_atom_self_overlaps_for_state
///
/// @brief
/// Returns a const reference to the atom-x-atom-pair vector-of-vectors of bools that specifies which atoms are overlapping,
/// assuming the alternate state assignment.
///
template < typename V, typename E, typename G >
utility::vector1< utility::vector1 < bool > > const &
HPatchBackgroundEdge< V, E, G >::get_atom_atom_overlaps_for_state( Size state ) const {
	assert( state <= node_states_overlap_with_bg_res_.size() );
	return node_states_overlap_with_bg_res_[ state ];
}


///
/// @begin HPatchBackgroundEdge::count_static_memory
///
template < typename V, typename E, typename G >
unsigned int HPatchBackgroundEdge< V, E, G >::count_static_memory() const {
	return sizeof ( HPatchBackgroundEdge< V, E, G > );
}

///
/// @begin HPatchBackgroundEdge::count_dynamic_memory
///
template < typename V, typename E, typename G >
unsigned int HPatchBackgroundEdge< V, E, G >::count_dynamic_memory() const {

	unsigned int total_memory = parent::count_dynamic_memory();
	total_memory += node_states_coverage_of_bg_res_.size() * sizeof ( RotamerDotsCache );
	total_memory += node_states_overlap_with_bg_res_.size() * sizeof( bool ); // underestimate; each element is a vector, too! (ronj)
	return total_memory;
}



//----------------------------------------------------------------------------//
//--------------------------- HPatch Interaction Graph -------------------------//
//----------------------------------------------------------------------------//

template < typename V, typename E, typename G >
Size HPatchInteractionGraph< V, E, G >::num_state_substitutions_considered_ = 0;

template < typename V, typename E, typename G >
Size HPatchInteractionGraph< V, E, G >::num_hpatch_comps_procrastinated_ = 0;

template < typename V, typename E, typename G >
Size HPatchInteractionGraph< V, E, G >::num_hpatch_comps_later_made_ = 0;

template < typename V, typename E, typename G >
bool HPatchInteractionGraph< V, E, G >::initialized_SASA_radii = false;

template < typename V, typename E, typename G >
utility::vector1< Real > HPatchInteractionGraph< V, E, G >::radii_;

template < typename V, typename E, typename G >
std::string HPatchInteractionGraph< V, E, G >::carbon_atom = "C";

template < typename V, typename E, typename G >
std::string HPatchInteractionGraph< V, E, G >::sulfur_atom = "S";

///
/// @begin HPatchInteractionGraph::HPatchInteractionGraph
///
/// @detailed
/// Main constructor. Initializes all member variables to 0 and false.
///
template < typename V, typename E, typename G >
HPatchInteractionGraph< V, E, G >::HPatchInteractionGraph( int num_nodes ) :
	AdditionalBackgroundNodesInteractionGraph< V, E, G > ( num_nodes ),
	hpatch_score_weight_( 1.0 ),
	num_total_residues_( 0 ),
	num_residues_assigned_as_background_( 0 ),
	some_node_in_state_0_( true ),
	fc_nodes_near_rotsub_( num_nodes, 0 ),
	fc_nodes_near_rotsub_bool_( num_nodes, true ),
	fc_exp_hphobe_djs_offsets_( num_nodes, 0 ),
	fc_n_exp_hphobes_( num_nodes, 0 ),
	prepared_for_simulated_annealing_( false ),
	observed_sufficient_hpatch_E_to_predict_min_( false ),
	hpatch_score_min_last_100_( 0 ),
	hpatch_score_min_recent_( 0 ),
	num_substitutions_since_hpatch_min_update_( 0 ),
	calculated_hpatch_deltaE_( false ),
	deltaE_for_substitution_( 0.0f ),
	node_considering_alt_state_( 0 ),
	alt_state_being_considered_( 0 ),
	total_energy_current_state_assignment_( 0 ),
	total_energy_alternate_state_assignment_( 0 ),
	hpatch_energy_current_state_assignment_( 0 ),
	hpatch_energy_alternate_state_assignment_( 0 ),
	num_commits_since_last_update_( 0 ),
	deltaE_threshold_for_avoiding_hpatch_calcs_( -1.0f )
{
	/// set all nodes as participating in a rotamer substitution if any node's state is 0 (unassigned)
	for ( Size ii = 1; ii <= fc_nodes_near_rotsub_.size(); ++ii ) { fc_nodes_near_rotsub_[ ii ] = ii; }
}


///
/// @begin HPatchInteractionGraph::~HPatchInteractionGraph
///
template < typename V, typename E, typename G >
HPatchInteractionGraph< V, E, G >::~HPatchInteractionGraph() {
	//TR_HIG << "called destructor" << std::endl;
}


///
/// @begin HPatchInteractionGraph::set_pose
///
/// @detailed
/// All throughout this class, I refer back to the original pose sequence. To be able to do that, I need to have a
/// handle to the pose in this class.  That's what this method provides.  In IGFactory.cc, this method gets called with
/// the pose object that's being packed/designed.
///
template < typename V, typename E, typename G >
void
HPatchInteractionGraph<V, E, G>::set_pose( pose::Pose const & pose ) {

#ifdef FILE_DEBUG
	TR_HIG << "set_pose() called: typeid() of base class G returned: " << typeid(G).name() << std::endl;
#endif

	// call the set_pose function in the LinMemIG class, because it, too, uses the Pose to do its thing
	if ( typeid(G) == typeid( pack::interaction_graph::LinearMemoryInteractionGraph ) ) {
		dynamic_cast<pack::interaction_graph::LinearMemoryInteractionGraph*>(this)->set_pose( pose );
	}

	pose_ = new pose::Pose( pose );
}


///
/// @begin HPatchInteractionGraph::set_packer_task
///
/// @brief
/// We need a copy of the packer task to figure out which residues are being packed and/or designed. We have to figure
/// the packing options because it determines whether a residue becomes a FirstClass (HPatchNode) node or a background node.
/// This method gets called in IGSupport.cc.
///
template < typename V, typename E, typename G >
void
HPatchInteractionGraph<V, E, G>::set_packer_task( task::PackerTask const & the_task ) {
	packer_task_ = the_task.clone();
}


///
/// @begin HPatchInteractionGraph::set_rotamer_sets
///
/// @detailed
/// It's nice to be able to print out information about the possible rotamer during initialization of the IG.
/// This method gets called in IGSupport.cc.
///
template < typename V, typename E, typename G >
void
HPatchInteractionGraph<V, E, G>::set_rotamer_sets( rotamer_set::RotamerSets const & rotsets ) {
	rotamer_sets_ = new rotamer_set::RotamerSets( rotsets );
}


///
/// @begin HPatchInteractionGraph::initialize()
///
/// @detailed
/// This function is the 1st major entry point (well, after the constructor) into the HIG. It needs to set residues that
/// are not changing as background residues
///
/// Oldest comments:
/// In ++, there's a file InteractionGraphSupport.cc which is an analog of the InteractionGraphFactory in mini.  In ++,
/// the InteractionGraphSupport file instantiates and initializes, depending on the command line switches, the right
/// interaction graph.  For the HPatchInteractionGraph, it first initializes the PDInteractionGraph (which is the base)
/// and then calls the HPatchIG initialize method.
///
/// The thing is that this initialize method can't be called when the graph is constructed in the InteractionGraphFactory.
/// The reason is that the PDInteractionGraph base initialize() method is NOT called until later in the pack rotamers
/// process.  (Actually it's called from within rotsets->compute_energies().)  Only after the rotsets->compute energies
/// returns can I stick in an initialize() method call for the HPatchInteractionGraph (HIG).  But then that's too late because
/// the rotsets object has already computed some energies. Perhaps that's ok though. The rotsets object only calculates
/// the PD energy terms - it doesn't do anything with non-PD terms.
///
/// If a HIG init method is called during construction of the HIG, then the init method that's called by the rotsets object
/// causes all the node info that was set in the HIG init method to be lost because the rotsets init method recreates all
/// the nodes in the interaction graph when it runs. (That was a fun behaviour to figure out.)
///
/// So the solution I'm going for is to call this init method in the prepare for simulated annealing method of this class.
/// That gets called just before SA starts, so it will do the task above then.  It doesn't really matter when the
/// task gets done as long as it happens before SA starts.  This also means that the HIG will now have to keep a reference
/// to the Pose, the Task, and the RotamerSets objects since it needs all of these things to do tasks 1) and 2).
/// (For the port of this HPatchIG, we might not need the task and rotamer sets objects.)
///
///	prepare_for_simulated_annealing gets called by the FixbbSA::run() method.  Before this method, the
/// rotamersets object has called compute_energies() (the whole process being started in pack_rotamers)
/// which calls initialize() on the IG.  I need to place the HIG init method directly after the IG init
/// method that the RS object calls.
///
/// Newest comments:
/// Don't call the parent classes initialize method because that calls create_new_node for all designable positions, making
/// every node in the graph a first class residue (erroneously!). Then after that, some residues are set as background
/// residues so I'm surprised things were even working at all!  Unfortunately, the class hierarchy makes setting designable
/// residues that are surface-exposed be the only FirstClassNodes difficult.  What has to happen is that non-designable
/// positions will be the BGNodes and molten positions will be FCNodes.  Then, each type of node will store a boolean flagging
/// whether it is surface-exposed or not.  Then, when a consider sub call is made, the nodes will immediately check the
/// value of the boolean and return 0 immediately if it's not surface-exposed.
///
/// 07/16/2008- Unfortunately, I can't just replicate all of the functionality of the parent classes initialize method
/// here because it's different for PDIGs and LMIGs.  And some of the their corresponding function calls are not defined
/// in both classes so I get compile errors if I try to do it this way.  There are a few possible solutions I could use
/// for the problem of schizophrenic node (nodes which are getting assigned as BG and FC nodes).  1) I could just use
/// Andrew's LFs designation that any designable or packable node is a FC node and any node not being packed or designed
/// is a BGNode.  This approach is the easiest to implement but since nodes that are only being packed aren't changing in
/// sequence, their surface score shouldn't be changing so it's an inefficiency.  2) I could change the "create_node"
/// call in this class to return either a Node or BGNode depending on what the packer_task has in it.  Then the parent
/// templated IG classes would get the right kind of node.  But then, I would have to potentially change the AddtlBGNodesIG
/// class method "setNumResiduesAsBackground" because in that call is a for loop which creates background nodes.  This
/// would complicate creating new non-PD terms in mini, because it would be specific for this case (surface). 3) Figure
/// out some OO way of handling the distinction between packable and nonpackable residues. The OO-correct way to do this
/// would probably be to subclass AdditionalBackgroundNodesInteractionGraph to something like
/// SASAExposedNodesInteractionGraph and then make the HPatchInteractionGraph inherit that one.  Then
/// "FirstClassNodes" would be nodes that are designable and surface-exposed; all other nodes would be "BackgroundNodes".
///
/// Since it would require alot of work to write another derived class for probably not that significant (or necessary)
/// performance improvement, I'll just make packable residues be FC nodes also.  So continue to call the parent classes
/// initialize() method.  Just set residues which are not packable nor designable to BG nodes.
///
/// Newest Newest comments": The above text was the case for the SurfaceIG. Turns out that we really do want all packable
/// and designable positions to be FC Nodes. So there's no inefficiency here, as was the case for the SurfaceIG.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::initialize( rotamer_set::RotamerSetsBase const & rot_sets_base ) {

	/// TEMP!!!!
	// APL wants to create a parent/base class for RotamerSets called RotamerSetsBase which will hold variables and functions
	// that are common to RotamerSets and FlexbbRotamerSets. Each interaction graph will be passed a RotamerSetsBase object
	// to its initialize method, and functions in this class should only access common base-class functions.  However,
	// this would require changing quite a bit of code and making sure that only base-class RotamerSetsBase methods are used
	// here which I'm not interested in doing right now. So we can cast the RotamerSetsBase object to a RotamerSets object
	// and use it the way it was previously.
	rotamer_set::RotamerSets const & rot_sets( static_cast< rotamer_set::RotamerSets const & > (rot_sets_base) );

#ifdef FILE_DEBUG
	TR_HIG << "initialize() called" << std::endl;
	TR_HIG << "calling set_rotamers on " << rotamer_sets().nmoltenres() << " molten residues." << std::endl;
#endif

	// parent refers to AdditionalBackgroundNodesIG; G refers to the templated IG: PD or LMIG
	// call the parent initialize() method
	G::initialize( rot_sets );

	// save references to the rotamer_set. these get used later to print information about considered subs, for example.
	for ( Size ii = 1; ii <= rotamer_sets().nmoltenres(); ++ii ) {
		get_hpatch_node( ii )->set_rotamers( rotamer_sets().rotamer_set_for_moltenresidue( ii ) );
	}

	// initializes some local variables that translate between the ig enumeration and resid
	set_num_residues_in_protein( pose().total_residue() );

	int nbackground = pose().total_residue() - rot_sets.nmoltenres();
	set_num_background_residues( nbackground );

	for ( Size ii = 1; ii <= pose().total_residue(); ++ii) {

		// in our case, first class residues are residues that are designable or packable (see notes in function comment).
		if ( packer_task().being_packed(ii) || packer_task().being_designed(ii) ) {
			// it's probably that being_packed() includes being_designed() so that the conditional needs only to say
			// being_packed(), but it'll short circuit the OR if being_packed() is true so it doesn't matter too much.
			continue;
		} else {
			set_residue_as_background_residue( ii );
			// init the RotamerDots object that each BGNode keeps for itself with its rotamer info. the line below calls
			// the IG method which then calls the BGNode method.
			set_background_residue_rotamer_dots( ii, pose().residue(ii) );
		}
	}

	for ( Size ii = 1; ii <= rotamer_sets().nmoltenres(); ++ii ) {
		// when G::initialize is called, we recreate all of the FC Nodes so that index ii in RotamerSets maps to FCNode ii.
		rotamer_set::RotamerSetCOP rsop = rotamer_sets().rotamer_set_for_moltenresidue( ii );

#ifdef FILE_DEBUG
		TR_HIG << "initialize: calling set_rotamer_dots_for_node_state for node " << ii << " with rotamers [";
#endif
		for ( Size jj = 1; jj <= rsop->num_rotamers(); ++jj ) {
			set_rotamer_dots_for_node_state( ii, jj, *(rsop->rotamer( jj )) );
		}
#ifdef FILE_DEBUG
		TR_HIG << "] " << std::endl;
#endif
	}

	//j setup the radii array, indexed by the atom type int. atom index for looking up an extra data type stored in the AtomTypes
	//ronj reads the values out of the database file sasa_radii.txt in the extras folder of atom_type_sets and stores the values
	//ronj for each atom type into the radii_ member variable. each index of the radii array corresponds to some atom type.
	init_SASA_radii_from_database();

#ifdef FILE_DEBUG
	TR_HIG << "initialize_hpatch_ig: DONE with initialization.\n---" << std::endl;
#endif

}


///
/// @begin HPatchInteractionGraph::create_new_node
///
/// @brief
/// factory method pattern for instantiation of HPatchNode objects, used by InteractionGraphBase class.
///
template < typename V, typename E, typename G >
NodeBase* HPatchInteractionGraph< V, E, G >::create_new_node( int node_index, int num_states ) {
#ifdef FILE_DEBUG
	TR_HIG << "create_new_node called with node_index " << node_index << " and num_states " << num_states << std::endl;
#endif
	return new HPatchNode< V, E, G >( this, node_index, num_states );
}


///
/// @begin HPatchInteractionGraph::create_new_edge
///
/// @brief
/// factory method pattern for instantiation of HPatchEdge objects, used by InteractionGraphBase class.
///
template < typename V, typename E, typename G >
EdgeBase* HPatchInteractionGraph< V, E, G >::create_new_edge( int index1, int index2 ) {
#ifdef FILE_DEBUG
	TR_HIG << "create_new_edge() called for indices " << index1 << " and " << index2 << std::endl;
#endif
	return new HPatchEdge< V, E, G >( this, index1, index2 );
}


///
/// @begin HPatchInteractionGraph::create_background_node
///
/// @brief
/// factory method pattern for instantiation of HPatchBackgroundNode objects, used by AdditionalBackgroundNodesInteractionGraph class.
///
template < typename V, typename E, typename G >
BackgroundNode< V, E, G >* HPatchInteractionGraph< V, E, G >::create_background_node( int node_index ) {
#ifdef FILE_DEBUG
	TR_HIG << "create_background_node() called for index " << node_index << std::endl;
#endif
	return new HPatchBackgroundNode< V, E, G >( this, node_index );
}


///
/// @begin HPatchInteractionGraph::create_background_edge
///
/// @brief
/// factory method pattern for instantiation of HPatchBackgroundEdge objects, used by AdditionalBackgroundNodesInteractionGraph class.
///
template < typename V, typename E, typename G >
BackgroundToFirstClassEdge< V, E, G >* HPatchInteractionGraph< V, E, G >::create_background_edge( int fc_node_index, int bg_node_index ) {
#ifdef FILE_DEBUG
	TR_HIG << "create_background_edge() called for indices " << fc_node_index << " and " << bg_node_index << std::endl;
#endif
	return new HPatchBackgroundEdge< V, E, G >( this, fc_node_index, bg_node_index );
}


///
/// @begin HPatchInteractionGraph::set_num_residues_in_protein
///
/// @brief
/// tells the graph how many residues there are total in the protein
///
/// @detailed
/// The graph maintains its own enumeration for the background residues; but asks that anyone wanting to refer to them
/// use their original resid. The graph has to switch back and forth between enumeration schemes and must know how many
/// residues there are total to do that efficiently.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::set_num_residues_in_protein( Size num_res ) {

	num_total_residues_ = num_res;
	resid_2_bgenumeration_.resize( num_total_residues_ );

	for ( Size ii = 1; ii <= num_total_residues_; ++ii ) {
		resid_2_bgenumeration_[ii] = 0;
	}

}


///
/// @begin HPatchInteractionGraph::set_num_background_residues
///
/// @brief
/// tells the graph how many residues there are as part of the protein that are not part of the combinatorial
/// optimization process -- they are part of the background
///
/// @detailed
/// The other half of switching between enumeration schemes for the background residues is knowing how many background residues there are.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::set_num_background_residues( Size num_background_residues ) {

#ifdef FILE_DEBUG
	TR_HIG << "set_num_background_residues: setting num background residues to " << num_background_residues << std::endl;
#endif

	parent::set_num_background_nodes( num_background_residues );
	if ( parent::get_num_background_nodes() == 0)
		return;

	bgenumeration_2_resid_.resize( parent::get_num_background_nodes() );
	bg_nodes_near_rotsub_.resize( parent::get_num_background_nodes() );
	for ( Size ii = 1; ii <= bg_nodes_near_rotsub_.size(); ++ii ) { bg_nodes_near_rotsub_[ ii ] = ii; }
	bg_nodes_near_rotsub_bool_.resize( parent::get_num_background_nodes(), true );
	bg_exp_hphobe_djs_offsets_.resize( parent::get_num_background_nodes(), 0 );
	bg_n_exp_hphobes_.resize( parent::get_num_background_nodes(), 0 );
	for ( Size ii = 1; ii <= (Size)parent::get_num_background_nodes(); ++ii ) {
		bgenumeration_2_resid_[ii] = 0;
	}
}


///
/// @begin HPatchInteractionGraph::set_residue_as_background_residue
///
/// @brief
/// informs the graph that a particular residue is part of the background
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::set_residue_as_background_residue( int residue ) {

	assert( resid_2_bgenumeration_[ residue ] == 0 );

	++num_residues_assigned_as_background_;
	resid_2_bgenumeration_[ residue ] = num_residues_assigned_as_background_;
	bgenumeration_2_resid_[ num_residues_assigned_as_background_ ] = residue;

#ifdef FILE_DEBUG
	TR_HIG << "set_residue_as_background_residue: set residue " << pose().residue( residue ).name3() << " "
		<< residue << " as background node " << num_residues_assigned_as_background_ << ". bgenumeration_2_resid_: [ ";
	for ( Size ii=1; ii <= bgenumeration_2_resid_.size(); ++ii ) {
		TR_HIG << bgenumeration_2_resid_[ ii ] << ", ";
	}
	TR_HIG << " ]" << std::endl;
#endif

}


///
/// @begin HPatchInteractionGraph::set_background_residue_rotamer_dots
///
/// @brief
/// Creates and inits a RotamerDots object for a background residue.  the residue must first have been declared to be a background residue.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::set_background_residue_rotamer_dots( Size residue, conformation::Residue const & rotamer ) {

	Size bgid = resid_2_bgenumeration_[ residue ];
	// the use of an OP here seems silly. though, storing a reference in a class is not as easy as it sounds.
	// it seems like a raw pointer here would make the most sense, since I'm not allocating any memory, but oh well.
	conformation::ResidueOP rotamer_op = new conformation::Residue( rotamer ); // calls copy constructor
	get_hpatch_bg_node( bgid )->set_rotamer( rotamer_op );

	bool exclude_hydrogens = true;
	bool use_expanded_polar_atom_radii = true;
	RotamerDots rd( rotamer_op, exclude_hydrogens, use_expanded_polar_atom_radii );
	get_hpatch_bg_node( bgid )->set_rotamer_dots( rd );

}


///
/// @begin HPatchInteractionGraph::set_rotamer_dots_for_node_state
///
/// @brief
/// store the coordinates for a particular rotamer
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::set_rotamer_dots_for_node_state( Size node_index, Size state, conformation::Residue const & rotamer ) {

	conformation::ResidueOP rotamer_op = new conformation::Residue( rotamer ); // calls copy constructor
	bool exclude_hydrogens = true;
	bool use_expanded_polar_atom_radii = true;
	RotamerDots rd( rotamer_op, exclude_hydrogens, use_expanded_polar_atom_radii );
#ifdef FILE_DEBUG
	TR_NODE << state << ":" << rotamer.name1() << "-" << rotamer.seqpos() << ", ";
#endif
	get_hpatch_node( node_index )->set_rotamer_dots_for_state( state, rd );

}


///
/// @begin HPatchInteractionGraph::prepare_for_simulated_annealing
///
/// @brief
/// Prepares the graph to begin simulated annealing.
///
/// @detailed
/// Invokes both base-class prepare_for_simulated_annealing subroutines: InteractionGraphBase first, to prepare the
/// HPatchNodes and HPatchEdges. Then the AdditionalBackgroundNodesInteractionGraph, to prepare the HPatchBackgroundNodes,
/// the HPatchBackgroundEdges, and to do a little more preparing of the HPatchNodes.
/// Also computes background/background overlap.
/// This is the 2nd major entry point into the HIG.
///
/// @remarks
/// As it is written, it should only be executed once; that will change, however if we make HPatchInteractionGraph work
/// with ligands (ligands that stay fixed during any single sim annealing process, but that move between anealings.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::prepare_for_simulated_annealing() {

	if ( prepared_for_simulated_annealing_ ) {

		// the function below just figures out which nodes and bgnodes to place edges between. this doesn't need to be repeated.
		//detect_background_residue_and_first_class_residue_overlap();

		// chains up to the IGBase method which, in turn, calls prep_for_simA() on all the FC Edges, and then all FC nodes.
		//G::prepare_for_simulated_annealing();
		// Edge::prep_for_simA() drops edges if 1) the two-body energies stored on that node are empty and 2) the incident
		// nodes on that edge don't have any sphere overlaps. None of this needs to be redone.
		//
		// Node::prep_for_simA() updates all the edge vectors in the parent classes (either PDIG or LinMemIG) and then
		// initializes the member variable self_and_bg_dots_for_states_. That is a vector1 of RotamerDots objects that contain
		// the self overlap each state causes and the overlap that all of that-Nodes-neighboring-BGNodes cause to a state.
		// this variable doesn't change during a trajectory so no need to reinitialize it.

		// parent::prepare() calls prep_for_simA() on all the BGEdges and BGNodes. the parent classes edge vectors are updated
		// and then initialize self overlap is called. the state on the BGNodes needs to be reinit'd to the point where the
		// self overlap and overlap from other bg nodes is stored. BGNode::prep_for_simA() also calls initialize_overlap_cache
		// to fill the node_states_coverage_of_bg_res_ vector. That part doesn't really need to be repeated, so it only gets
		// done if the BGEdge::prep_for_simA boolean isn't set. The reason we do need to call initialize_overlap_cache() though
		// is that the RotamerDotsCache objects get resized (and hence cleared) in that function.
		parent::prepare_for_simulated_annealing();

		// this method should get repeated to get BGNodes back to the right state
		initialize_bg_bg_overlaps();

		return;
	}

	detect_background_residue_and_first_class_residue_overlap();

	// G::prepare() calls InteractionGraphBase::prepare_for_simulated_annealing() - LinmemIG implements one but it also
	// calls the IGBase method.  The IGBase method, in turn, calls prep_for_simA() on all the FC Edges, and then all FC nodes.
	G::prepare_for_simulated_annealing();

	// parent::prepare() calls prep_for_simA() on all the BGNodes
	parent::prepare_for_simulated_annealing();

#ifdef FILE_DEBUG
	TR_HIG << "prepare_for_simulated_annealing(): initializing background-background overlap" << std::endl;
#endif

	initialize_bg_bg_overlaps();
	initialize_bg_bg_atom_atom_overlaps();

	prepared_for_simulated_annealing_ = true;

#ifdef FILE_DEBUG
	TR_HIG << "prepare_for_simulated_annealing: number edges in graph: " << parent::get_num_edges() << std::endl;
#endif

}


///
/// @begin HPatchInteractionGraph::detect_background_residue_and_first_class_residue_overlap
///
/// @brief
/// iterates across all pairs of first- and second-class nodes to determine which share sphere overlaps.
/// Adds a HPatchBackgroundEdge between any pair that do.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::detect_background_residue_and_first_class_residue_overlap() {

	for ( Size ii = 1; ii <= (Size)parent::get_num_background_nodes(); ++ii ) {

#ifdef FILE_DEBUG
		TR_HIG << "detect_bg_and_fc_residue_neighbors: checking for neighbors of background residue " << pose().residue( bgenumeration_2_resid_[ ii ] ).name3()
			<< " " << pose().residue( bgenumeration_2_resid_[ ii ] ).seqpos() << std::endl;
#endif

		for ( Size jj = 1; jj <= (Size)parent::get_num_nodes(); ++jj ) {

			// ii: background node index, jj: first-class node index
			// the background nodes and first-class nodes should be set at this point due to the initalize method above.
			// if the background node overlaps with the first-class node, we need to add an edge between them.
			if ( get_hpatch_bg_node( ii )->detect_overlap( get_hpatch_node( jj ) ) ) {
#ifdef FILE_DEBUG
				TR_HIG << "detect_bg_residue_and_fc_residue_overlap: --- adding FC/BG edge: fc node id:" << jj
					<< " / bg node: " << ii << ", bg resid:" << bgenumeration_2_resid_[ ii ] << std::endl;
#endif
				parent::add_background_edge( jj, ii );
			}
		}
	}

#ifdef FILE_DEBUG
	TR_HIG << "DONE detecting background and first class overlap.\n---" << std::endl;
#endif

}


///
/// @begin HPatchInteractionGraph::initialize_bg_bg_overlaps
///
/// @brief
/// computes the background overlaps for all background node pairs
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::initialize_bg_bg_overlaps() {

	for ( Size ii = 1; ii <= (Size)parent::get_num_background_nodes(); ++ii ) {
		for ( Size jj = ii + 1; jj <= (Size)parent::get_num_background_nodes(); ++jj ) {
			// get_hpatch_bg_node() returns a pointer, so we must dereference it to pass it
			get_hpatch_bg_node( ii )->initialize_bg_bg_overlap( *get_hpatch_bg_node( jj ) );
		}
	}

}


///
/// @begin HPatchBackgroundNode::get_bg_bg_atom_atom_overlaps
///
//template < typename V, typename E, typename G >
template < typename V, typename E, typename G >
utility::vector1< utility::vector1< bool > > const &
HPatchInteractionGraph< V, E, G >::get_bg_bg_atom_atom_overlaps( Size node1_index, Size node2_index ) {

	// only the positions where node2_index is greater than node1_index will have atom-atom overlap information
	assert( node1_index < node2_index );

	return bg_bg_atom_atom_overlaps_[ node1_index ][ node2_index ];

}

///
/// @begin HPatchBackgroundNode::initialize_bg_bg_atom_atom_overlaps
///
/// @brief
/// initializes the atom-atom overlap vector stored by the IG for all the bg-bg node overlaps.
///
/// @detailed
/// During simulated annealing, the IG has to determine the connected components after every sub. To do this, it has to
/// check a very large number of atom pairs for overlap.  These pairs include intra-Node atom-pairs, intra-BGNode atom-pairs,
/// BGNode-BGNode atom-pairs, BGNode-FCNode atom-pairs, and FCNode-FCNode atom-pairs.  Now intra-Node/BGNode atom-pairs
/// can be stored on the Nodes and BGNodes. Similarly, BG-Node-FCNode and FCNode-FCNode atom-pairs can be held on the
/// BGEdges and FCEdges, respectively.  Unfortunately, edges between BGNodes do not exist in the IG.  So, we have to init
/// a vector that will be held by the IG when we initialize the dot counts and sphere overlaps for all the BGNode pairs.
/// The reason we do this is so that we don't have to go and repeat the calculation of all these sphere overlaps after
/// every sub.  These will remain constant during the course of a simulation.
///
/// The hard part is figuring out what kind of data structure I can use here that will make it easy for the IG to get
/// what it needs when it's calculating a score.
///
/// Making the data structure a 4D vector, i.e. a vector of vectors of vectors of vectors of bools.  The outer two vectors
/// are the two bgnode ids we're checking overlap for.  The innermost two vectors are the atom-atom overlap information
/// for two given bgnode ids. Only the higher indexed nodes have to contain vectors of vectors of bools. For instance,
/// the outer two vectors will have (1,2), (1,3), (1,4), ... (2,3), (2,4), (2,5), ..., (3,4), (3,5), (3,6), etc
/// The inner two vectors is just an all-atoms-in-ii x all-atoms-in-jj.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::initialize_bg_bg_atom_atom_overlaps() {

	Real const probe_radius = 1.4;

	bg_bg_respairs_w_hphobe_olap_.resize( (Size) parent::get_num_background_nodes() );
	bg_bg_atom_atom_overlaps_.resize(     (Size) parent::get_num_background_nodes() );
	curr_bg_bg_exhpobeolap_.resize(       (Size) parent::get_num_background_nodes() );
	alt_bg_bg_exhpobeolap_.resize(        (Size) parent::get_num_background_nodes() );

	for ( Size bgnode_ii = 1; bgnode_ii <= (Size)parent::get_num_background_nodes(); ++bgnode_ii ) {

		conformation::ResidueCOP bgnode_ii_rotamer = get_hpatch_bg_node( bgnode_ii )->get_rotamer();
		Size const ii_natoms = bgnode_ii_rotamer->nheavyatoms();
		utility::vector1< Size > const & ii_hphobes( get_hpatch_bg_node( bgnode_ii )->get_hphobes() );

		for ( Size bgnode_jj = bgnode_ii + 1; bgnode_jj <= (Size)parent::get_num_background_nodes(); ++bgnode_jj ) {

			conformation::ResidueCOP bgnode_jj_rotamer = get_hpatch_bg_node( bgnode_jj )->get_rotamer();
			Size const jj_natoms = bgnode_jj_rotamer->nheavyatoms();
			utility::vector1< Size > const & jj_hphobes( get_hpatch_bg_node( bgnode_jj )->get_hphobes() );

			utility::vector1< utility::vector1< Size > > ii_jj_overlap( ii_natoms, utility::vector1< Size >( jj_natoms, false ) );

			bool any_hphobe_olap = false;
			for ( Size ii=1; ii <= ii_hphobes.size(); ++ii ) {
				Size const iia = ii_hphobes[ ii ];

				conformation::Atom const & iia_atom( bgnode_ii_rotamer->atom( iia ) );
				utility::vector1< Real > const & atom_radii = radii_;
				Real const iia_atom_radius = atom_radii[ iia_atom.type() ] + probe_radius; // radii_ is a class member variable

				for ( Size jj=1; jj <= jj_hphobes.size(); ++jj ) {
					Size const jja = jj_hphobes[ jj ];

					conformation::Atom const & jja_atom( bgnode_jj_rotamer->atom( jja ) );
					Real const jja_atom_radius = atom_radii[ jja_atom.type() ] + probe_radius;

					// check if the two atoms overlap at all; use distance squared over distance to keep things faster.
					Vector const & iia_atom_xyz = bgnode_ii_rotamer->atom( iia ).xyz();
					Vector const & jja_atom_xyz = bgnode_jj_rotamer->atom( jja ).xyz();

					Real const distance_squared = iia_atom_xyz.distance_squared( jja_atom_xyz );

					if ( distance_squared <= (iia_atom_radius + jja_atom_radius) * (iia_atom_radius + jja_atom_radius) ) {

						Real const distance_ijxyz = std::sqrt( distance_squared );
						int degree_of_overlap;
						core::scoring::get_overlap( iia_atom_radius, jja_atom_radius, distance_ijxyz, degree_of_overlap );
						if ( degree_of_overlap >= 15 ) {
#ifdef FILE_DEBUG
							/*TR_HIG << "initialize_bg_bg_atom_atom_overlaps(): overlapping bg-bg atom-atom pair: "
								<< bgnode_ii_rotamer->seqpos() << "/" << utility::trim( bgnode_ii_rotamer->atom_name( iia ) ) << " - "
								<< bgnode_jj_rotamer->seqpos() << "/" << utility::trim( bgnode_jj_rotamer->atom_name( jja ) )
								<< ", degree of overlap: " << degree_of_overlap << std::endl;*/
#endif
							ii_jj_overlap[ iia ][ jja ] = true;
							any_hphobe_olap = true;
						}
					} // end if distance
				} // for loop over all jj hphobes
			} // for loop over all ii hphobes

			if ( any_hphobe_olap ) {
				bg_bg_respairs_w_hphobe_olap_[ bgnode_ii ].push_back( bgnode_jj );
				bg_bg_atom_atom_overlaps_[ bgnode_ii ].push_back( ii_jj_overlap );
				// just get the dimensions right for now, set these values to false later
				curr_bg_bg_exhpobeolap_[ bgnode_ii ].push_back( ii_jj_overlap );
				alt_bg_bg_exhpobeolap_[ bgnode_ii ].push_back( ii_jj_overlap );
			}
		} // node_jj
	} // node_ii



#ifdef FILE_DEBUG
	/*for ( Size bgnode_ii = 1; bgnode_ii <= (Size)parent::get_num_background_nodes(); ++bgnode_ii ) {
		TR_HIG << "background nodes overlapping with background node " << bgnode_ii << ": ";
		for ( Size jj = 1; jj <= bg_bg_respairs_w_hphobe_olap_[ bgnode_ii ].size(); ++jj ) {
			TR_HIG << bg_bg_respairs_w_hphobe_olap_[ bgnode_ii ][ jj ] << ", ";
		}
		TR_HIG << std::endl;

		for ( Size bgnode_jj = 1; bgnode_jj <= bg_bg_respairs_w_hphobe_olap_[ bgnode_ii ].size(); ++bgnode_jj ) {
			TR_HIG << "background node " << bgnode_ii << " x background node " << bg_bg_respairs_w_hphobe_olap_[ bgnode_ii ][ bgnode_jj ] << " overlap: " << std::endl;
			for ( Size ii_atom = 1; ii_atom <= bg_bg_atom_atom_overlaps_[ bgnode_ii ][ bgnode_jj ].size(); ++ii_atom ) {
				TR_HIG << "bgnode " << bgnode_ii << " atom " << ii_atom << ": [ ";
				for ( Size jj_atom = 1; jj_atom <= bg_bg_atom_atom_overlaps_[ bgnode_ii ][ bgnode_jj ][ ii_atom ].size(); ++jj_atom ) {
					TR_HIG << bg_bg_atom_atom_overlaps_[ bgnode_ii ][ bgnode_jj ][ ii_atom ][ jj_atom ] << ", ";
				}
				TR_HIG << "]" << std::endl;
			}
			TR_HIG << std::endl;
		}
		TR_HIG << std::endl;
	}*/
#endif


}


///
/// @begin HPatchInteractionGraph::blanket_assign_state_0
///
/// @brief
/// assigns state 0 -- the unassigned state -- to all (first class) vertices in the graph
///
/// @detailed
/// This is the 3rd entry point into the HIG.  It is called by the Annealer just before simulated annealing and rotamer
/// substitutions begin to init the graph to unassigned values everywhere.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::blanket_assign_state_0() {

	for ( Size ii = 1; ii <= (Size)parent::get_num_nodes(); ++ii ) {
#ifdef FILE_DEBUG
		TR_HIG << "blanket_assign_state_0() calling assign_zero_state() on node " << ii << " of " << (Size)parent::get_num_nodes() << std::endl;
#endif
		get_hpatch_node( ii )->assign_zero_state();
	}

	//update_internal_energy_totals_hpatch();
	// instead of calling update_internal, just reset the cached energy values to zero. the values should already be
	// at zero because they are default init'd by the constructor to be 0, but the unit tests re-use the same IG and
	// just call prep_for_simA() and this method to reinit the IG.
	total_energy_current_state_assignment_ = total_energy_alternate_state_assignment_ = 0.0;
	hpatch_energy_current_state_assignment_ = hpatch_energy_alternate_state_assignment_ = 0.0;

	/// set all nodes as participating in a rotamer substitution if any node's state is 0 (unassigned)
	some_node_in_state_0_ = true;
	fc_nodes_near_rotsub_.resize( (Size)parent::get_num_nodes() );
	bg_nodes_near_rotsub_.resize( (Size)parent::get_num_background_nodes() );
	
	for ( Size ii = 1; ii <= fc_nodes_near_rotsub_.size(); ++ii ) { fc_nodes_near_rotsub_[ ii ] = ii; }
	for ( Size ii = 1; ii <= fc_nodes_near_rotsub_bool_.size(); ++ii ) { fc_nodes_near_rotsub_bool_[ ii ] = true; }

#ifdef FILE_DEBUG
		TR_HIG << "blanket_assign_state_0() reset fc_nodes_near_rotsub_ to [ ";
		for ( Size ii = 1; ii <= fc_nodes_near_rotsub_.size(); ++ii ) { TR_HIG << fc_nodes_near_rotsub_[ ii ] << ", "; }
		TR_HIG << "]" << std::endl;

		TR_HIG << "blanket_assign_state_0() reset fc_nodes_near_rotsub_bool_ to [ ";
		for ( Size ii = 1; ii <= fc_nodes_near_rotsub_bool_.size(); ++ii ) { TR_HIG << fc_nodes_near_rotsub_bool_[ ii ] << ", "; }
		TR_HIG << "]" << std::endl;
#endif

	for ( Size ii = 1; ii <= bg_nodes_near_rotsub_.size(); ++ii ) { bg_nodes_near_rotsub_[ ii ] = ii; }
	for ( Size ii = 1; ii <= bg_nodes_near_rotsub_bool_.size(); ++ii ) { bg_nodes_near_rotsub_bool_[ ii ] = true; }

#ifdef FILE_DEBUG
		TR_HIG << "blanket_assign_state_0() reset bg_nodes_near_rotsub_ to [ ";
		for ( Size ii = 1; ii <= bg_nodes_near_rotsub_.size(); ++ii ) { TR_HIG << bg_nodes_near_rotsub_[ ii ] << ", "; }
		TR_HIG << "]" << std::endl;

		TR_HIG << "blanket_assign_state_0() reset bg_nodes_near_rotsub_bool_ to [ ";
		for ( Size ii = 1; ii <= bg_nodes_near_rotsub_bool_.size(); ++ii ) { TR_HIG << bg_nodes_near_rotsub_bool_[ ii ] << ", "; }
		TR_HIG << "]" << std::endl;
#endif

}


///
/// @begin HPatchInteractionGraph::update_internal_energy_totals_hpatch
///
/// @brief
/// After every 2^10 commits, the graph traverses its nodes and edges and
/// re-tallies the total energy of the current state assignment.  This update
/// prevents the accumulation of numerical drift, increasing accuracy.
///
/// @detailed
/// this function becomes less necessary in this implementation since the graph itself calculates the score. no longer do we
/// have the situation where the nodes/bgnodes keep the score and the IG just applies deltaE's every commit. now the IG calculates
/// the score for every consider/commit. the only "drift" that may accumulate would result from not calculating de novo
/// the PD current state energy for a long time. that's why I'm leaving this function in.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::update_internal_energy_totals_hpatch() {

	HPatchInteractionGraph< V, E, G >::print_hpatch_avoidance_stats();

	parent::update_internal_energy_totals();
	total_energy_current_state_assignment_ = parent::get_energy_PD_current_state_assignment() + hpatch_energy_current_state_assignment_;

	num_commits_since_last_update_ = 0;

#ifdef DOUBLE_CHECK_COUNTS
	verify_sasas_correct();
#endif

}


#ifdef DOUBLE_CHECK_COUNTS
///
/// @begin HPatchInteractionGraph::verify_sasas_correct
///
/// @brief
/// Verifies that the SASAs held on the Nodes and BGNodes are correct. Only runs every 1000 commits and if DOUBLE_CHECK_COUNTS is defined.
///
/// @detailed
/// Constructs a new Pose object based on the current state of the IG and runs the calc_per_atom_sasa function in sasa.cc
/// as an independent way of making sure that the Nodes and BGNodes are correctly keeping track of SASA. Very slow.
/// Should only be used for testing.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::verify_sasas_correct() {

	using namespace ObjexxFCL::format;

	utility::vector1< Real > node_sasas( parent::get_num_nodes(), 0.0 );
	utility::vector1< Real > bgnode_sasas( parent::get_num_background_nodes(), 0.0 );

	//TR_HIG << "verify_sasas_correct(): Node SASAs: [ ";
	for ( Size ii=1; ii <= (Size)parent::get_num_nodes(); ++ii ) {
		node_sasas[ ii ] = (get_hpatch_node(ii)->get_current_state_rotamer_dots()).get_sasa();
		//TR_HIG << node_sasas[ ii ] << ", ";
	}
	//TR_HIG << " ]" << std::endl;

	//TR_HIG << "verify_sasas_correct(): BGNode SASAs: [ ";
	for ( Size ii=1; ii <= (Size)parent::get_num_background_nodes(); ++ii ) {
		bgnode_sasas[ ii ] = (get_hpatch_bg_node(ii)->get_current_state_rotamer_dots()).get_sasa();
		//TR_HIG << bgnode_sasas[ ii ] << ", ";
	}
	//TR_HIG << " ]" << std::endl;

	// create a copy of the passed in pose...
	pose::Pose pose_copy = pose();

	// ...and place the current state rotamers on it
	for ( core::uint ii = 1; ii <= rotamer_sets().nmoltenres(); ++ii ) {
		core::uint iiresid = rotamer_sets().moltenres_2_resid( ii );
		core::uint iicurrrot = get_hpatch_node(ii)->get_current_state();
		conformation::ResidueCOP bestrot( rotamer_sets().rotamer_set_for_moltenresidue( ii )->rotamer( iicurrrot ) );

		conformation::ResidueOP newresidue( bestrot->create_residue() );
		pose_copy.replace_residue( iiresid, *newresidue, false );
	}

	Real total_sasa = 0.0;
	id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa;
	Real probe_radius = 1.4;

	// create an atom_subset mask that calculates sasa for heavyatoms and not hydrogens
	id::AtomID_Map< bool > atom_subset;
	atom_subset.resize( pose_copy.n_residue() );
	// init all heavy atoms to true and all H's to false
	for ( Size ii=1; ii <= pose_copy.n_residue(); ++ii ) {
		atom_subset.resize( ii, pose_copy.residue_type(ii).natoms(), false );
		for ( Size jj = 1; jj <= pose_copy.residue_type(ii).nheavyatoms(); ++jj ) {
			atom_subset[ ii ][ jj ] = true;
		}
	}

	total_sasa = core::scoring::calc_per_atom_sasa( pose_copy, atom_sasa, rsd_sasa, probe_radius, false /* no big polar H */, atom_subset, true /* use_naccess_sasa_radii */ );

	bool incorrect_sasa_found = false;

	//TR_HIG << "verify_sasas_correct(): calc_per_atom_sasa Node SASAs: [ ";
	for ( Size ii=1; ii <= node_sasas.size(); ++ii ) {
		//TR_HIG << rsd_sasa[ rotamer_sets().moltenres_2_resid( ii ) ] << ", ";
		if ( fabs( node_sasas[ ii ] - rsd_sasa[ rotamer_sets().moltenres_2_resid( ii ) ] ) > 0.1 )
			incorrect_sasa_found = true;
	}
	//TR_HIG << " ]" << std::endl;

	//TR_HIG << "verify_sasas_correct(): calc_per_atom_sasa BGNode SASAs: [ ";
	// all the background nodes should have some value for SASA
	for ( Size ii=1; ii <= bgnode_sasas.size(); ++ii ) {
		//TR_HIG << rsd_sasa[ bg_node_2_resid( ii ) ] << ", ";
		if ( fabs( bgnode_sasas[ ii ] - rsd_sasa[ bg_node_2_resid( ii ) ] ) > 0.1 )
			incorrect_sasa_found = true;
	}
	//TR_HIG << " ]" << std::endl;

	if ( incorrect_sasa_found )
		utility_exit_with_message( "SASA values are out of sync. Something is wrong. Quitting." );

	TR_HIG << "verify_sasas_correct() called.  and checked out" << std::endl;

}
#endif

///
/// @begin HPatchInteractionGraph::set_errorfull_deltaE_threshold
///
/// @brief
/// Allows the sim-annealer to specify a deltaE threshold above which, it is no longer necessary to be very accurate.
///
/// @detailed
/// When the annealer asks the graph to consider a state substitution that produces a large collision, the graph may
/// approximate the hpatch deltaE instead of performing expensive sphere overlap computations.  The deltaE returned by
/// consider_substitution() will be inaccurate, but if the annealer is unlikely to accept the substitution, then time
/// can be saved. The graph guarantees that if the annealer does commit that substitution that it will go back and
/// perform the hpatch computations and return an accurate total energy for the graph.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::set_errorfull_deltaE_threshold( core::PackerEnergy deltaE ) {

#ifdef FILE_DEBUG
	TR_HIG << "set_errorfull_deltaE_threshold: setting threshold to " << deltaE << std::endl;
	// leave this inside since it's debugging output
	HPatchInteractionGraph< V, E, G >::print_hpatch_avoidance_stats();
#endif

	HPatchInteractionGraph< V, E, G >::reset_hpatch_avoidance_stats();
	deltaE_threshold_for_avoiding_hpatch_calcs_ = deltaE;

}


///
/// @begin HPatchInteractionGraph::print_hpatch_avoidance_stats
///
/// @brief
/// reports on the level of success for hpatch score calculation procrastination
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::print_hpatch_avoidance_stats() {

	if ( num_state_substitutions_considered_ == 0 )
		return;

	TR_STATS << "num state substitutions considered: " << num_state_substitutions_considered_ << ", "
			<< "num hpatch calcs procrastinated: " << num_hpatch_comps_procrastinated_ << ", "
			<< "num hpatch calcs later computed: " << num_hpatch_comps_later_made_ << std::endl;
	TR_STATS << "Percent Avoided: " << (double) (num_hpatch_comps_procrastinated_ - num_hpatch_comps_later_made_) / num_state_substitutions_considered_ << ", ";

	if ( num_hpatch_comps_procrastinated_ != 0 ) {
		TR_STATS << "Worthwhile Procrastination: " << (double) (num_hpatch_comps_procrastinated_ - num_hpatch_comps_later_made_) / num_hpatch_comps_procrastinated_ << std::endl;
	} else {
		TR_STATS << "Worthwhile Procrastination: " << "N/A" << std::endl;
	}

}

///
/// @begin HPatchInteractionGraph::reset_hpatch_avoidance_stats
///
/// @brief
/// resets static member variables of HPatchIG that measure how worthwhile hpatch calculation procrastination is.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::reset_hpatch_avoidance_stats() {
	num_state_substitutions_considered_ = 0;
	num_hpatch_comps_procrastinated_ = 0;
	num_hpatch_comps_later_made_ = 0;
}


///
/// @begin HPatchInteractionGraph::consider_substitution
///
/// @brief
/// Returns the (possibly approximate) change in energy induced by switching a particular node from its currently assigned state to some alternate state.
///
/// @detailed
/// First, queries the HPatchNode for the pairwise-decomposable (PD) energy. If the PD difference
/// implies a collision, then the HPatchIG pretends as if the state substitution causes the best
/// improvement possible in hpatch score, and returns the PD difference + pretend hpatch difference.
/// It will procrastinate computing the actual hpatch score difference until the guiding SimAnnealer
/// decides to commit the substitution. If the SimAnnealer rejects the substitution, then the work
/// to compute the hpatch score is never done. If it is unclear that the SimAnnealer will reject the
/// substitution based on the PD difference, then the Graph goes ahead and computes the change in hpatch
/// score accurately.
///
/// This function is the 4th major entry point from the Annealer into the HIG.
///
/// Also returns the sum of the two body energies for the node in its current state; the sim-annealer accepts state
/// substitutions at higher chance if the original state was also at a poor energy.
///
/// @param
/// node_ind - [in] - the index of the (first class) node
/// new_state - [in] - the alternate state that the node should consider
/// delta_energy - [out] - the change in energy induced on the entire graph by substituting a node's current state with the alternate.
///							This energy may be inaccurate if it exceeds a threshold set by the sim-annealer.
/// prev_energy_for_node - [out] - the sum of the pair-wise decomposable portion of the energy function for the node's currently assigned state
///
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::consider_substitution( int node_ind, int new_state, core::PackerEnergy & delta_energy, core::PackerEnergy & prev_energy_for_node ) {

#ifdef FILE_DEBUG
	TR_HIG << "---" << std::endl;
	TR_HIG << "consider_substitution(): trying new state " << new_state << " ("
		<< get_hpatch_node( node_ind )->get_rotamer( new_state )->name() << ") on node/molten res " << node_ind << " (wt: "
		<< pose().residue( rotamer_sets().moltenres_2_resid( node_ind ) ).name3() << " " << rotamer_sets().moltenres_2_resid( node_ind ) << ") " << std::endl;
#endif

	reset_from_previous_deltaHpatch_comp();

	++num_state_substitutions_considered_;

	node_considering_alt_state_ = node_ind;
	alt_state_being_considered_ = new_state;

	// the below deltaE may be an estimate of the change in energy and not the actual value
	core::PackerEnergy deltaE = get_hpatch_node( node_ind )->calculate_PD_deltaE_for_substitution( new_state, prev_energy_for_node );

	calculated_hpatch_deltaE_ = false;

	if ( decide_procrastinate_hpatch_computations( deltaE, deltaE_threshold_for_avoiding_hpatch_calcs_ ) ) {
		//TR_HIG << "procrastinated" << std::endl;
		++num_hpatch_comps_procrastinated_;
	} else {
#ifdef FILE_DEBUG
		TR_HIG << "deltaE for PD: " << deltaE << std::endl;
#endif
		deltaE += calculate_hpatch_deltaE();
		calculated_hpatch_deltaE_ = true;;
	}

	delta_energy = deltaE;
	total_energy_alternate_state_assignment_ = deltaE + total_energy_current_state_assignment_;

#ifdef FILE_DEBUG
	TR_HIG << "total deltaE for substitution: " << delta_energy << std::endl;
#endif

}


///
/// @begin HPatchInteractionGraph< V, E, G >::calculate_hpatch_deltaE()
///
/// @detailed
/// Goes through the entire process of calculating the hpatch deltaE for a substitution.
///
template < typename V, typename E, typename G >
core::PackerEnergy HPatchInteractionGraph< V, E, G >::calculate_hpatch_deltaE() {

	// update the sasa information on all the relevant nodes and bg nodes
#ifndef FILE_DEBUG
	get_hpatch_node( node_considering_alt_state_ )->consider_alternate_state(); // don't assign deltaSASA to a variable
#else
	Real delta_sasa = get_hpatch_node( node_considering_alt_state_ )->consider_alternate_state();
	if ( delta_sasa != 0.0 ) {
		TR_HIG << "delta sasa: " << delta_sasa << std::endl;
	}
#endif

#ifdef FILE_DEBUG
	// on the very first consider sub call after all nodes get an assigned state, the some_node_in_state_0_ variable will still be true.
	// that's because it's value doesn't change until after the calculate_alt_state_hpatch_score() function call below. so we'll miss out
	// on which nodes/bgnodes will be affected by the very first substitution.  that's ok.
	if ( ! some_node_in_state_0_ ) {
		TR_HIG << "calculate_hpatch_deltaE(): FC nodes affected by substitution currently being considered include: ";
		for ( Size ii = 1; ii <= fc_nodes_near_rotsub_.size(); ++ii ) { TR_HIG << fc_nodes_near_rotsub_[ ii ] << " "; } TR_HIG << std::endl;

		TR_HIG << "calculate_hpatch_deltaE(): BG nodes affected by substitution currently being considered include: ";
		for ( Size ii = 1; ii <= bg_nodes_near_rotsub_.size(); ++ii ) { TR_HIG << bg_nodes_near_rotsub_[ ii ] << " "; } TR_HIG << std::endl;

		TR_HIG.flush();
	}
#endif

	// determine the new patch score given the updated sasa information
	hpatch_energy_alternate_state_assignment_ = calculate_alt_state_hpatch_score();
	core::PackerEnergy hpatch_deltaE = hpatch_energy_alternate_state_assignment_ - hpatch_energy_current_state_assignment_;
#ifdef FILE_DEBUG
	TR_HIG << "calculate_hpatch_deltaE(): hpatchE current state: " << hpatch_energy_current_state_assignment_
		<< ", alt state: " << hpatch_energy_alternate_state_assignment_ << ", hpatch deltaE: " << hpatch_deltaE << std::endl;
#endif

	return hpatch_deltaE;
}


///
/// @begin HPatchInteractionGraph< V, E, G >:: register_fc_node_in_state0()
/// 
/// @detailed
/// Initialized to true in the constructor, and also set to true after a call to blanket_assign_state0.  After all nodes
/// get assigned a state, then this boolean is set to false.  It is used to reduce the number of operations that are performed
/// when annealing is just beginning. Functionality added by Andrew - best commenting I can do. -RJ
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >:: register_fc_node_in_state0() {
	some_node_in_state_0_ = true;
}


///
/// @begin HPatchInteractionGraph::register_fc_node_affected_by_rotsub()
///
/// @brief
/// Called by HPatchNodes to specify which first-class Nodes are affected by a given rotamer substitution. If a given
/// consider gets rejected, then fc_nodes_near_rotsub_ gets cleared when reset_from_previous_deltaHpatch_comp.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::register_fc_node_affected_by_rotsub( int fc_node_ind ) {

//#ifdef FILE_DEBUG
	//TR_HIG << "register_fc_node_affected_by_rotsub(): node " << fc_node_ind << " registered as being affected by substitution being considered." << std::endl;
//#endif

	// not sure this conditional makes sense. why would we ever have a first class node try to register themselves twice with the IG?  if the IG is
	// keeping track of nodes/bgnodes correctly, a single FC node should only ever get a chance to call this method one time. also, because fc_nodes_near_rotsub_bool_
	// is reset to true for all positions in the blanket_assign_state_0() call, we potentially neglect to add some FC nodes to the fc_nodes_near_rotsub_ vector
	// because the value is already true.  that's bad. -ronj
	// 2/5/2013:  if you don't comment out the conditional below (and the analogous one in the function after this one), the IG begins to fail when multiple
	// packing runs are used (e.g. in pmut_scan_scan protocol).  however, if you do comment them out, some of the unit tests begin to fail.  need to update the 
	// unit tests with a more appropriate test. or figure out why the unit tests fail when they are commented out and fix the IG. 
	
	//if ( ! fc_nodes_near_rotsub_bool_[ fc_node_ind ] ) {
		fc_nodes_near_rotsub_.push_back( fc_node_ind );
		fc_nodes_near_rotsub_bool_[ fc_node_ind ] = true;
	//}

	// we instead could place an if statement above that says if any node is in the unassigned state, then don't bother updating either of these vectors. but 
	// as these functions (this one and the one below) get called millions of times, better not to stick an if statement that only applies during the beginning of
	// annealing smack-dab in the middle of them.  the implication of this, though, is that after the first consider-substitution-that-happens-once-every-node-goes-into
	// the-assigned-state (regardless of whether it gets committed or not), the fc_nodes_near_rotsub_ vectors will have duplicate values in them (lots of duplicates,
	// potentially).  that's because they will be getting values added to them while the annealer slowly gets all the nodes to assigned states. The multi cool annealer 
	// will reach this very quickly as it force inits all nodes to assigned states. the standard fixbb annealer could take a long time to reach assigned states for every
	// node.  after the first substitutions processing is done, the IG will reset after that sub consideration.  in there, it will iterate through these vectors and reset
	// the nodes and bgnodes.  it will visit some nodes and bgnodes more than once, but at the beginning of the reset...() function calls in nodes and bgnodes, there's 
	// a check for whether anything needs to be done.  after the first visit, repeat visits to nodes/bgnodes won't require much runtime.  I'm thinking this expense that
	// occurs only once when the IG become fully assigned is better than having two extra if statements in this function and the one below (which get called millions
	// of times).  the only real way to know would be to check runtimes of a couple design runs, but I don't feel like doing that. so I'm going with taking out the if
	// statements and accepting the one-time cost that occurs once the IG is fully assigned. crossing my fingers. -ronj
	// 

}

template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::register_bg_node_affected_by_rotsub( int bg_node_ind ) {

//#ifdef FILE_DEBUG
	//TR_HIG << "register_bg_node_affected_by_rotsub(): bgnode " << bg_node_ind << " registered as being affected by substitution being considered." << std::endl;
//#endif
	//if ( ! bg_nodes_near_rotsub_bool_[ bg_node_ind ] ) {
		bg_nodes_near_rotsub_.push_back( bg_node_ind );
		bg_nodes_near_rotsub_bool_[ bg_node_ind ] = true;
	//}
}


///
/// @begin HPatchInteractionGraph::init_SASA_radii_from_database()
///
/// @brief
/// Puts the SASA radii values in the database into the passed in radii array. Temporary location. Should only be done once.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::init_SASA_radii_from_database() {

	core::chemical::AtomTypeSet const & atom_type_set = pose().residue(1).atom_type_set();
	radii_.resize( atom_type_set.n_atomtypes(), 0.0 );

	core::Size SASA_RADIUS_INDEX = atom_type_set.extra_parameter_index( "NACCESS_SASA_RADIUS" );
	for ( core::Size ii=1; ii <= radii_.size(); ++ii ) {
		radii_[ii] = atom_type_set[ii].extra_parameter( SASA_RADIUS_INDEX );
	}

	initialized_SASA_radii = true;
}


///
/// @begin HPatchInteractionGraph::update_disjoint_sets_using_cache()
///
/// @brief
/// Helper function for calculating the hpatch score. This function is specific for intra-residue overlaps, the one
/// below is for inter-residue overlaps. The reason they're separate is because the inner loops have different start locations.
/// intra-residue nested loops only have to go over greater-indexed atoms because of commutativity. inter-residue nested
/// loops have to do all ii-res x jj-res atom-atom overlaps.
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::update_disjoint_sets_using_cache(
	conformation::Residue const &, // rsd,
	InvRotamerDots const & invdots,
	utility::vector1< Size > const & exp_hphobes,
	Size residue_djs_offset,
	utility::vector1< utility::vector1< bool > > const & atom_atom_overlaps,
	graph::DisjointSets & ds
) {

	for ( Size ii=1, iiend = exp_hphobes.size(); ii <= iiend; ++ii ) {

		Size const iia = exp_hphobes[ii];
		Size const ii_djs_id = residue_djs_offset + ii;

		// start iterating jja from iia + 1!
		// for intra-residue atom pairs, the entries below and on the diagonal are not necessary
		for ( Size jj = ii + 1; jj <= iiend; ++jj ) {
			Size const jja( exp_hphobes[jj] );
			Size const jj_djs_id = residue_djs_offset + jj;

			if ( atom_atom_overlaps[ iia ][ jja ] ) {

				/*TR_HIG << "update_disjoint_sets_using_cache(): overlapping atom pair: "
					<< rsd.aa() << " " << node_index << "/" << utility::trim( rsd.atom_name( iia ) ) << " - "
					<< rsd.aa() << " " << node_index << "/" << utility::trim( rsd.atom_name( jja ) ) << std::endl;

				Real const ii_rad = RotamerDots::radius_for_attype( rsd.atom(iia).type() ) + 1.4;
				Real const jj_rad = RotamerDots::radius_for_attype( rsd.atom(jja).type() ) + 1.4;
				if ( rsd.xyz( iia ).distance_squared( rsd.xyz( jja ) ) > (ii_rad+jj_rad)*(ii_rad+jj_rad) ) {
					std::cerr << "ERROR discrepancy between atom_atom_overlaps array and actual overlap information" << std::endl;
					std::cerr << "rsd.seqpos() " << rsd.seqpos() << std::endl;
					std::cerr << "Atom ii= " << iia << " rad: " << ii_rad << " Atom jj= " << jja << " rad: " << jj_rad << " sum: " << ii_rad+jj_rad << std::endl;
					std::cerr << "Distance " << rsd.xyz(iia).distance( rsd.xyz(jja) ) << std::endl;
				}*/

				if ( ds.ds_find( ii_djs_id ) == ds.ds_find( jj_djs_id ) ) continue; // fast
				if (  ! invdots.atom_overlap_is_exposed( iia, jja ) ) continue;     // slow
				ds.ds_union( ii_djs_id, jj_djs_id);
			}
		}
	}

	return;
}


///
/// @begin HPatchInteractionGraph::update_disjoint_sets_using_cache()
///
/// @brief
/// Helper function for calculating the hpatch score.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::update_disjoint_sets_using_cache(
	conformation::Residue const & , // rsd1,
	InvRotamerDots const & invdots1,
	utility::vector1< Size > const & exp_hphobes1,
	Size djs_offset_1,
	conformation::Residue const & , // rsd2,
	InvRotamerDots const & invdots2,
	utility::vector1< Size > const & exp_hphobes2,
	Size djs_offset_2,
	utility::vector1< utility::vector1< bool > > const & atom_atom_overlaps,
	graph::DisjointSets & ds
)
{

	for ( Size ii=1, iiend = exp_hphobes1.size(), jjend = exp_hphobes2.size(); ii <= iiend; ++ii ) {
		Size const iia = exp_hphobes1[ ii ];
		Size const ii_djs_index = djs_offset_1 + ii;

		//TR_HIG << "update_disjoint_sets_using_cache: iia: " << iia << std::endl;
		// start iterating jja from 1!  this are inter-residue atom pairs so all of them are important

		for ( Size jj = 1; jj <= jjend; ++jj ) {
			Size const jja = exp_hphobes2[ jj ];
			Size const jj_djs_index = djs_offset_2 + jj;

			//TR_HIG << "\tupdate_disjoint_sets_using_cache: jja: " << jja << std::endl;
			if ( atom_atom_overlaps[ iia ][ jja ] ) {

				/*TR_HIG << "update_disjoint_sets_using_cache(): overlapping atom pair: "
					<< rsd1.aa() << " " << rsd1.seqpos() << "/" << utility::trim( rsd1.atom_name( iia ) ) << " - "
					<< rsd2.aa() << " " << rsd2.seqpos() << "/" << utility::trim( rsd2.atom_name( jja ) ) << std::endl;

				//Real const ii_rad = RotamerDots::radius_for_attype( rsd1.atom(iia).type() ) + 1.4;
				//Real const jj_rad = RotamerDots::radius_for_attype( rsd2.atom(jja).type() ) + 1.4;

				Real const ii_rad = RotamerDots::epradius_for_attype( rsd1.atom(iia).type() ) + 1.4;
				Real const jj_rad = RotamerDots::epradius_for_attype( rsd2.atom(jja).type() ) + 1.4;

				if ( rsd1.xyz( iia ).distance_squared( rsd2.xyz( jja ) ) > (ii_rad + jj_rad) * (ii_rad + jj_rad) ) {
					std::cerr << "ERROR discrepancy between atom_atom_overlaps array and actual overlap information" << std::endl;
					std::cerr << "rsd1.seqpos() " << rsd1.seqpos() << std::endl;
					std::cerr << "rsd2.seqpos() " << rsd2.seqpos() << std::endl;
					std::cerr << "Atom ii= " << iia << " rad: " << ii_rad << " Atom jj= " << jja << " rad: " << jj_rad << " sum: " << ii_rad+jj_rad << std::endl;
					std::cerr << "Distance " << rsd1.xyz(iia).distance( rsd2.xyz(jja) ) << std::endl;
					std::cerr << "Distance swapped? ";
					if ( iia <= rsd2.natoms() && jja <= rsd1.natoms() ) {
						std::cerr << rsd2.xyz(iia).distance( rsd1.xyz(jja) ) << " w/ iirad = " <<
						RotamerDots::radius_for_attype( rsd2.atom(iia).type() ) + 1.4 << " and jjrad = " <<
						RotamerDots::radius_for_attype( rsd1.atom(jja).type() ) + 1.4 << " sum: " <<
						RotamerDots::radius_for_attype( rsd2.atom(iia).type() ) + 1.4 + RotamerDots::radius_for_attype( rsd1.atom(jja).type() ) + 1.4;
					}
					std::cerr << std::endl;

				}*/

				if ( ds.ds_find( ii_djs_index ) == ds.ds_find( jj_djs_index ) )  continue; // fast
				if (  ! invdots1.atom_overlap_is_exposed( iia, invdots2, jja ) ) continue; // slow
				ds.ds_union( ii_djs_index, jj_djs_index );
			}
		}
	}

	return;
}


///
/// @begin HPatchInteractionGraph::calculate_alt_state_hpatch_score()
///
/// @brief
/// Constructs an atom-level graph for the alternate state assignment on the IG and then runs the union-find algorithm on
/// it to obtain the connected components.  Each connected component is a surface hydrophobic patch.  What the total energy
/// of the IG should be hasn't yet been decided.  One possibility will be to assign every CC/patch a score and return the
/// sum as the total score.
///
/// Queries all of the Nodes and BGNodes to determine the connectivity assuming the alternate state.
///
template < typename V, typename E, typename G >
Real HPatchInteractionGraph< V, E, G >::calculate_alt_state_hpatch_score() {

	// any_vertex_state_unassigned() is an O(N) operation -- only perform this check if a few times
	// at the beginning of simA
	if ( some_node_in_state_0_ && parent::any_vertex_state_unassigned() ) {
#ifdef FILE_DEBUG
		if ( some_node_in_state_0_ ) {
			TR_HIG << "calculate_alt_state_hpatch_score(): some_node_in_state_0_ is true. returning 0.0." << std::endl;
		}
#endif
		return 0.0; // don't bother running union-find if any of the nodes are still unassigned
	}

	some_node_in_state_0_ = false;

	if ( !initialized_SASA_radii ) {
		// setup radii array, temporarily located here
		//j setup the radii array, indexed by the atom type int. atom index for looking up an extra data type stored in the AtomTypes
		//ronj reads the values out of the database file sasa_radii.txt in the extras folder of atom_type_sets and stores the values
		//ronj for each atom type into the radii array. each index of the radii array corresponds to some atom type.
		init_SASA_radii_from_database();
	}

	// Figure out how many exposed hydrophobic atoms exist in the structure, and create a mapping between
	// these exp-hphobes and indices for nodes in a disjoint-sets data structure.
	// Determining this value so we can init the DisjointSets object with that number of "nodes". Instead of init'ing it to ALL
	// heavyatoms in a pose, init it to just all exposed hydrophobic atoms.

	Size tot_exp_hphobes( 0 ); // stands for total number of exposed hydrophobic ATOMS
	djs_id_2_hphobe_index_.clear(); // is a vector1 of exposed_hydrophobic_data structs

	for ( Size ii = 1; ii <= (Size) parent::get_num_nodes(); ++ii ) {
		fc_exp_hphobe_djs_offsets_[ ii ] = tot_exp_hphobes; // fc node 1 will have an offset of 0, node 2 will have offset however many exposed hp atoms node 1 has, etc

		// get the number of exposed hp atoms on this fc node. the number depends on whether this node is near the node considering substitution.
		Size ii_n_exhphobes( fc_nodes_near_rotsub_bool_[ ii ] ? get_hpatch_node( ii )->n_alt_state_exp_hphobes() : get_hpatch_node( ii )->n_curr_state_exp_hphobes() );
		fc_n_exp_hphobes_[ ii ] = ii_n_exhphobes;

		tot_exp_hphobes += ii_n_exhphobes;
		for ( Size jj = 1; jj <= ii_n_exhphobes; ++jj ) {
			// exposed_hydrophobic_data is a struct; 1 indicates first-class node, ii and jj are the node id and index into the exp hp atoms vector respectively
			djs_id_2_hphobe_index_.push_back( exposed_hydrophobic_data( 1, ii, jj ) );
		}
	}

	for ( Size ii = 1; ii <= (Size) parent::get_num_background_nodes(); ++ii ) {
		bg_exp_hphobe_djs_offsets_[ ii ] = tot_exp_hphobes;

		// get the number of exposed hp atoms on this bg node.
		Size ii_n_exhphobes( bg_nodes_near_rotsub_bool_[ ii ] ? get_hpatch_bg_node( ii )->n_alt_state_exp_hphobes() : get_hpatch_bg_node( ii )->n_curr_state_exp_hphobes() );
		bg_n_exp_hphobes_[ ii ] = ii_n_exhphobes;

		tot_exp_hphobes += ii_n_exhphobes;
		for ( Size jj = 1; jj <= ii_n_exhphobes; ++jj ) {
			// exposed_hydrophobic_data is a struct; 2 indicates background node, ii and jj are the node id and index into the exp hp atoms vector respectively
			djs_id_2_hphobe_index_.push_back( exposed_hydrophobic_data( 2, ii, jj ) );
		}
	}

	//sasa_for_djs_node_.resize( tot_exp_hphobes );
	ep_sasa_for_djs_node_.resize( tot_exp_hphobes );

	// now fill in the values for the ep_sasa_for_djs_node_ array. instead of storing the SASA/EPSASA for all heavy atoms in
	// the pose, the exposed hp vectors allow us to only store the SASA/EPSASA of exposed hydrophobic atoms. Uses less memory,
	// and is much faster to do operations on.  Remember a djs node is really an exposed hydrophobic atom.
	for ( Size ii = 1; ii <= (Size) parent::get_num_nodes(); ++ii ) {
		Size const ii_offset = fc_exp_hphobe_djs_offsets_[ ii ];
		Size const ii_n_exhphobes( fc_n_exp_hphobes_[ ii ] ); // not strictly necessary but this allows for preventing vector resize operations
		utility::vector1< Size > const & ii_exhphobes( fc_nodes_near_rotsub_bool_[ ii ] ? get_hpatch_node( ii )->alt_state_exp_hphobes() : get_hpatch_node( ii )->curr_state_exp_hphobes() );

		// for the number of exposed hp atoms, get the atom index from the ii_exhphobes vector and convert the atom index
		// to a djs node id by adding the offset for first-class node ii
		for ( Size jj = 1; jj <= ii_n_exhphobes; ++jj ) {
			Size const jj_atom = ii_exhphobes[ jj ];
			Size const jj_djs_node_id = ii_offset + jj;
			//sasa_for_djs_node_[ jj_djs_node_id ] = get_hpatch_node( ii )->get_alt_state_rotamer_dots().get_atom_sasa( jj_atom );
			ep_sasa_for_djs_node_[ jj_djs_node_id ] = get_hpatch_node( ii )->get_alt_state_rotamer_dots().get_atom_sasa( jj_atom );
		}
	}

	for ( Size ii = 1; ii <= (Size) parent::get_num_background_nodes(); ++ii ) {
		Size const ii_offset = bg_exp_hphobe_djs_offsets_[ ii ];
		Size const ii_n_exhphobes( bg_n_exp_hphobes_[ ii ] );
		utility::vector1< Size > const & ii_exhphobes( bg_nodes_near_rotsub_bool_[ ii ] ? get_hpatch_bg_node( ii )->alt_state_exp_hphobes() : get_hpatch_bg_node( ii )->curr_state_exp_hphobes() );

		for ( Size jj = 1; jj <= ii_n_exhphobes; ++jj ) {
			Size const jj_atom = ii_exhphobes[ jj ];
			Size const jj_djs_node_id = ii_offset + jj;
			//sasa_for_djs_node_[ jj_djs_node_id ] = get_hpatch_bg_node( ii )->get_alt_state_rotamer_dots().get_atom_sasa( jj_atom );
			ep_sasa_for_djs_node_[ jj_djs_node_id ] = get_hpatch_bg_node( ii )->get_alt_state_rotamer_dots().get_atom_sasa( jj_atom );
		}
	}

	// the above loops have determined the number of exposed hydrophobic atoms that are present. now init the disjointsets object.
	graph::DisjointSets ds( tot_exp_hphobes );

	// intra-residue connections
	for ( Size ii = 1; ii <= (Size)parent::get_num_nodes(); ++ii ) {
		if ( fc_n_exp_hphobes_[ ii ] == 0 ) continue;

		conformation::ResidueCOP rsd( ii == node_considering_alt_state_ ? get_hpatch_node( ii )->curr_state_rotamer() : get_hpatch_node( ii )->alt_state_rotamer() );
		Size const ii_state = ( ii == node_considering_alt_state_ ? alt_state_being_considered_ : get_hpatch_node( ii )->get_current_state() );
		utility::vector1< utility::vector1< bool > > const & atom_atom_self_overlaps = get_hpatch_node( ii )->get_atom_atom_self_overlaps_for_state( ii_state );
		utility::vector1< Size > const & ii_exhphobes( fc_nodes_near_rotsub_bool_[ ii ] ? get_hpatch_node( ii )->alt_state_exp_hphobes() : get_hpatch_node( ii )->curr_state_exp_hphobes() );

		update_disjoint_sets_using_cache( *rsd, get_hpatch_node( ii )->alt_state_inv_dots(), ii_exhphobes, fc_exp_hphobe_djs_offsets_[ ii ], atom_atom_self_overlaps, ds );

	}
	for ( Size ii = 1; ii <= (Size) parent::get_num_background_nodes(); ++ii ) {
		if ( bg_n_exp_hphobes_[ ii ] == 0 ) continue;

		utility::vector1< Size > const & ii_exp_hphobes( get_hpatch_bg_node( ii )->alt_state_exp_hphobes() );
		conformation::Residue const & ii_rsd( * get_hpatch_bg_node( ii )->get_rotamer() );
		utility::vector1< utility::vector1< bool > > const & atom_atom_self_overlaps = get_hpatch_bg_node( ii )->get_atom_atom_self_overlaps();

		update_disjoint_sets_using_cache( ii_rsd, get_hpatch_bg_node( ii )->alt_state_inv_dots(), ii_exp_hphobes, bg_exp_hphobe_djs_offsets_[ ii ], atom_atom_self_overlaps, ds );
	}


	/// For all of the inter-residue connections, we have a vector of ints that's stored on each of the Nodes and BGNodes.
	/// The vector specifies which atom indexes in that ResidueType are carbon or sulfur atoms *and* are exposed. Then,
	/// when iterating over all atom-atom pairs, we don't have to check each time to see if an atom is hydrophobic and
	/// exposed. We end up never even visiting the polar atoms which is awesome.

	// bgnode-bgnode connections
	// apl -- why not keep a graph?
	for ( Size ii = 1; ii <= (Size) parent::get_num_background_nodes(); ++ii ) {
		if ( bg_n_exp_hphobes_[ ii ] == 0 ) continue;
		if ( bg_bg_respairs_w_hphobe_olap_[ ii ].size() == 0 ) continue;

		Size const ii_djs_offset = bg_exp_hphobe_djs_offsets_[ ii ];

		conformation::Residue const & ii_rsd( * get_hpatch_bg_node( ii )->get_rotamer() );
		utility::vector1< Size > const & ii_exp_hphobes( get_hpatch_bg_node( ii )->alt_state_exp_hphobes() );

		// instead of iterating over all higher-indexed background nodes, only iterate over the ones that bgnode ii has overlap with.
		// this can be a significant time savings if there are alot of bgnodes in a simulation.
		for ( Size jj = 1; jj <= bg_bg_respairs_w_hphobe_olap_[ ii ].size(); ++jj ) {

			// all of the vectors are indexed on node id - not on the value of jj!
			// not true. bg_bg_respairs_w_hphobe_olap_ and bg_bg_atom_atom_overlaps_ are indexed on jj, not the node id
			Size jj_bg_node_ind = bg_bg_respairs_w_hphobe_olap_[ ii ][ jj ];

			// bg_bg_respairs_w_hphobe_olap_ holds the bg node ids. bg_bg_respairs_w_hphobe_olap_[ii][jj] will be the bg node id
			// that has overlap with bgnode ii.
			//if ( bg_n_exp_hphobes_[ jj ] == 0 ) continue;
			if ( bg_n_exp_hphobes_[ jj_bg_node_ind ] == 0 ) continue; // bg_n_exp_hphobes_ is indexed on node id, not jj

			//Size const jj_djs_offset = bg_exp_hphobe_djs_offsets_[ ii ];
			Size const jj_djs_offset = bg_exp_hphobe_djs_offsets_[ jj_bg_node_ind ];

			conformation::Residue const & jj_rsd( * get_hpatch_bg_node( jj_bg_node_ind )->get_rotamer() );
			utility::vector1< Size > const & jj_exp_hphobes( get_hpatch_bg_node( jj_bg_node_ind )->alt_state_exp_hphobes() );

			utility::vector1< utility::vector1< bool > > const & atom_atom_overlaps = bg_bg_atom_atom_overlaps_[ ii ][ jj ];

#ifdef FILE_DEBUG
			/*TR_HIG << "bg_bg_respairs_w_hphobe_olap_: [ ";
			for ( Size aa = 1; aa <= bg_bg_respairs_w_hphobe_olap_[ ii ].size(); ++aa ) {
				TR_HIG << bg_bg_respairs_w_hphobe_olap_[ ii ][ aa ] << ", ";
			}
			TR_HIG << "]" << std::endl;

			TR_HIG << "ii_exp_hphobes: [ ";
			for ( Size aa = 1; aa <= ii_exp_hphobes.size(); ++aa ) {
				TR_HIG << ii_exp_hphobes[ aa ] << ", ";
			}
			TR_HIG << "]" << std::endl;

			TR_HIG << "jj_exp_hphobes: [ ";
			for ( Size aa = 1; aa <= jj_exp_hphobes.size(); ++aa ) {
				TR_HIG << jj_exp_hphobes[ aa ] << ", ";
			}
			TR_HIG << "]" << std::endl;

			TR_HIG << "background node " << ii << " x background node " << jj_bg_node_ind << " overlap: " << std::endl;
			for ( Size ii_atom = 1; ii_atom <= ii_exp_hphobes.size(); ++ii_atom ) {
				TR_HIG << "bgnode " << ii << " atom " << ii_exp_hphobes[ ii_atom ] << ": [ ";
				for ( Size jj_atom = 1; jj_atom <= jj_exp_hphobes.size(); ++jj_atom ) {
					TR_HIG << atom_atom_overlaps[ ii_exp_hphobes[ ii_atom ] ][ jj_exp_hphobes[ jj_atom ] ] << ", ";
				}
				TR_HIG << "]" << std::endl;
			}
			TR_HIG << std::endl;*/
#endif

			update_disjoint_sets_using_cache( ii_rsd, get_hpatch_bg_node( ii )->alt_state_inv_dots(), ii_exp_hphobes, ii_djs_offset,
				jj_rsd, get_hpatch_bg_node( jj_bg_node_ind )->alt_state_inv_dots(), jj_exp_hphobes, jj_djs_offset, atom_atom_overlaps, ds );
		}
	}

#ifdef FILE_DEBUG
	TR_HIG << "calculate_alt_state_hpatch_score(): iterating over first-class edges" << std::endl;
	TR_HIG.flush();
#endif

	//
	// first-class edges
	// now we have to iterate over every edge (and bg edge) in the IG and do pairwise-atom comparisons to see if two atoms should be
	// assigned to the same connected component. some of these edges will be connected to the node considering a substitution. that
	// node should have broadcast the sub being considered to all of its neighboring FC nodes. therefore, these edges will have
	// different values for the overlap between the current state and alt state. have to make sure we get the alt state overlap
	// for these edges.  for all other FC edges, the alt state could be anything. it could be the same as the current state, or
	// it could be some alt_state that was considering a long time ago that wasn't commit'd(). so for these edges, we have to use
	// the overlap values for the current state rotamers at both nodes.
	//

	for ( std::list<EdgeBase*>::iterator iter = parent::get_edge_list_begin(); iter != parent::get_edge_list_end(); ++iter ) {

		Size node0_index = ((HPatchEdge< V, E, G >*)(*iter))->get_first_node_ind();
		Size node1_index = ((HPatchEdge< V, E, G >*)(*iter))->get_second_node_ind();

		// if the number of exposed hydrophobic atoms on either FC node is zero, then no point in doing work on this edge.
		// why? because, if none of the atoms are exposed, then there's no way we can connect nodes in the disjoint sets
		// object for the nodes on this edge.
		if ( fc_n_exp_hphobes_[ node0_index ] == 0 || fc_n_exp_hphobes_[ node1_index ] == 0 ) continue;

		conformation::ResidueCOP ii_rsd( node0_index == node_considering_alt_state_ ?
			get_hpatch_node( node0_index )->get_rotamer( alt_state_being_considered_ ) :
			get_hpatch_node( node0_index )->get_rotamer( get_hpatch_node( node0_index )->get_current_state()  ) );

		utility::vector1< Size > const & ii_exp_hphobes( get_hpatch_node( node0_index )->alt_state_exp_hphobes() );
		Size const ii_djs_offset = fc_exp_hphobe_djs_offsets_[ node0_index ];

#ifdef FILE_DEBUG
		//Size node0_current_state = get_hpatch_node( node0_index )->get_current_state();
		//TR_HIG << "calculate_alt_state_hpatch_score(): E(" << node0_index << "," << node1_index << "), ";
		//TR_HIG << "node " << node0_index << ": " << node0_current_state << " " << ii_rsd->name3();
#endif

		conformation::ResidueCOP jj_rsd( node1_index == node_considering_alt_state_ ?
			get_hpatch_node( node1_index )->get_rotamer( alt_state_being_considered_ ) :
			get_hpatch_node( node1_index )->get_rotamer( get_hpatch_node( node1_index )->get_current_state()  ) );

		utility::vector1< Size > const & jj_exp_hphobes( get_hpatch_node( node1_index )->alt_state_exp_hphobes() );
		Size const jj_djs_offset = fc_exp_hphobe_djs_offsets_[ node1_index ];

#ifdef FILE_DEBUG
		//Size node1_current_state = get_hpatch_node( node1_index )->get_current_state();
		//TR_HIG << ", node " << node1_index << ": " << node1_current_state << " " << jj_rsd->name3();
		//TR_HIG << std::endl;
#endif

		// it's possible we have to reorder the nodes so that the atom-atom overlap information is interpreted correctly.
		// this situation arises only when one of the nodes is the one that's considering a substitution.
		if ( node0_index == node_considering_alt_state_ || node1_index == node_considering_alt_state_ ) {
			utility::vector1< utility::vector1< bool > > const & atom_atom_overlaps = ((HPatchEdge< V, E, G >*)(*iter))->get_alt_state_atom_atom_overlaps();
			if ( node0_index == node_considering_alt_state_ ) {
				/// order the function parameters so that node0 is the outer-loop residue
				update_disjoint_sets_using_cache(
					*ii_rsd, get_hpatch_node( node0_index )->alt_state_inv_dots(), ii_exp_hphobes, ii_djs_offset,
					*jj_rsd, get_hpatch_node( node1_index )->alt_state_inv_dots(), jj_exp_hphobes, jj_djs_offset,
					atom_atom_overlaps, ds );
			} else {
				/// order the function parameters so that node1 is the outer-loop residue
				update_disjoint_sets_using_cache(
					*jj_rsd, get_hpatch_node( node1_index )->alt_state_inv_dots(), jj_exp_hphobes, jj_djs_offset,
					*ii_rsd, get_hpatch_node( node0_index )->alt_state_inv_dots(), ii_exp_hphobes, ii_djs_offset,
					atom_atom_overlaps, ds );
			}
		} else {
			utility::vector1< utility::vector1< bool > > const & atom_atom_overlaps = ((HPatchEdge< V, E, G >*)(*iter))->get_current_state_atom_atom_overlaps();

			/// order the function parameters so that node1 is the inner-loop residue
			update_disjoint_sets_using_cache(
				*ii_rsd, get_hpatch_node( node0_index )->alt_state_inv_dots(), ii_exp_hphobes, ii_djs_offset,
				*jj_rsd, get_hpatch_node( node1_index )->alt_state_inv_dots(), jj_exp_hphobes, jj_djs_offset,
				atom_atom_overlaps, ds );
		}

	} // for loop over all edges

#ifdef FILE_DEBUG
	TR_HIG << "calculate_alt_state_hpatch_score(): iterating over background edges" << std::endl;
	TR_HIG.flush();
#endif

	//
	// background edges
	// now we have to iterate over every bgedge and do the same as for the FC edges. bgedges are susceptible to the same problem
	// as described above for fc edges. it may be that the node of a bgedge is the node that is considering a sub for this
	// consider() call. in that case, we have to use the overlap that the bgnode and the alt_state on the node have for this
	// function.  alternatively, if the node isn't the one considering a sub, then we need to use the current state on the
	// fc node and bgnode overlap.
	//

	typename std::list< core::pack::interaction_graph::BackgroundToFirstClassEdge< V, E, G >* >::const_iterator iter;
	for ( iter = parent::get_bg_edge_list_begin(); iter != parent::get_bg_edge_list_end(); ++iter ) {

		Size fc_node_index = ((HPatchBackgroundEdge< V, E, G >*)(*iter))->get_first_class_node_index();
		Size bg_node_index = ((HPatchBackgroundEdge< V, E, G >*)(*iter))->get_background_node_index();
		if ( fc_n_exp_hphobes_[ fc_node_index ] == 0 || bg_n_exp_hphobes_[ bg_node_index ] == 0 ) continue;

		Size fc_node_current_state = get_hpatch_node( fc_node_index )->get_current_state();

		conformation::ResidueCOP ii_rsd( fc_node_index == node_considering_alt_state_ ?
			get_hpatch_node( fc_node_index )->get_rotamer( alt_state_being_considered_ ) :
			get_hpatch_node( fc_node_index )->get_rotamer( fc_node_current_state ) );

		utility::vector1< Size > const & fc_exp_hphobes( get_hpatch_node( fc_node_index )->alt_state_exp_hphobes() );
		Size const fc_djs_offset = fc_exp_hphobe_djs_offsets_[ fc_node_index ];

#ifdef FILE_DEBUG
		//TR_HIG << "calculate_alt_state_hpatch_score(): bgE(" << fc_node_index << "," << bg_node_index << "), ";
		//TR_HIG << "node " << fc_node_index << ": " << fc_node_current_state << " " << ii_rsd->name3();
#endif
		conformation::ResidueCOP jj_rsd( pose().residue( bgenumeration_2_resid_[ bg_node_index ] ) );
		utility::vector1< Size > const & bg_exp_hphobes( get_hpatch_bg_node( bg_node_index )->alt_state_exp_hphobes() );
		Size const bg_djs_offset = bg_exp_hphobe_djs_offsets_[ bg_node_index ];


#ifdef FILE_DEBUG
		//TR_HIG << ", bgnode " << bg_node_index << ": " << jj_rsd->name3() << std::endl;
#endif

		// get the cached atom-atom overlap in the context of the alt_state if the fc node on this bgedge is the node considering a sub
		utility::vector1< utility::vector1< bool > > const & atom_atom_overlaps(
			fc_node_index == node_considering_alt_state_ ?
			((HPatchBackgroundEdge< V, E, G >*)(*iter))->get_atom_atom_overlaps_for_state( alt_state_being_considered_ ) :
			((HPatchBackgroundEdge< V, E, G >*)(*iter))->get_atom_atom_overlaps_for_state( fc_node_current_state ) );

		update_disjoint_sets_using_cache(
			*ii_rsd, get_hpatch_node( fc_node_index )->alt_state_inv_dots(), fc_exp_hphobes, fc_djs_offset,
			*jj_rsd, get_hpatch_bg_node( bg_node_index )->alt_state_inv_dots(), bg_exp_hphobes, bg_djs_offset,
			atom_atom_overlaps, ds );

	} // for loop over all bg edges


	//
	// finally, reprint the largest patch size atoms for easy access
	//
#ifdef FILE_DEBUG
	utility::vector1< Size > set_sizes = ds.disjoint_set_sizes();
	Size const largest_set_size_index = utility::arg_max( set_sizes ); // vector1.functions.hh
	TR_HIG << "num atoms in largest patch: : " << set_sizes[ largest_set_size_index ] << std::endl;
#endif

	Size n_ccs( 0 ); // number of connected components
	reps_for_nonzero_rank_ccs_.resize( 0 );
	djs_rep_node_index_2_cc_index_.resize( tot_exp_hphobes, 0 ); // use "0" to say "I'm not the representative for any CC"
	sasa_for_cc_.resize( 0 ); // the patch area

	for ( Size ii = 1; ii <= tot_exp_hphobes; ++ii ) {
		Size ii_rep = ds.ds_find( ii );
		if ( ds.node( ii_rep ).rank == 0 ) continue; /// 1-atom connected component.
		Size ii_cc = djs_rep_node_index_2_cc_index_[ ii_rep ];
		if ( ii_cc == 0 ) {
			++n_ccs;
			reps_for_nonzero_rank_ccs_.push_back( ii_rep );
			sasa_for_cc_.push_back( 0.0 );
			djs_rep_node_index_2_cc_index_[ ii_rep ] = n_ccs;
			ii_cc = n_ccs;
		}
		sasa_for_cc_[ ii_cc ] += ep_sasa_for_djs_node_[ ii ];
	}

	core::Real total_alt_state_hpatch_score = 0.0;
	for ( Size ii = 1; ii <= sasa_for_cc_.size(); ++ii ) {

		Real const patch_area = sasa_for_cc_[ ii ];
		Real score = 0.0;
		if ( patch_area > SurfacePotential::MAX_HPATCH_AREA ) {
			score = hpatch_score_weight_ * SurfacePotential::MAX_HPATCH_SCORE;
		} else {
			score = hpatch_score_weight_ * SurfacePotential::get_instance()->hpatch_score( patch_area );
		}

		total_alt_state_hpatch_score += score;

	} // end loop over all patches/connected components

	/// now unwind the reps_for_nonzero_rank_ccs_ data:
	for ( Size ii = 1; ii <= reps_for_nonzero_rank_ccs_.size(); ++ii ) {
		djs_rep_node_index_2_cc_index_[ reps_for_nonzero_rank_ccs_[ ii ]] = 0;
	}


#ifdef FILE_DEBUG
		std::map< Size, utility::vector1< Size > > sets = ds.sets();
		std::map< Size, utility::vector1< Size > >::iterator it;

		core::Real patch_area = 0.0;
		for ( it = sets.begin() ; it != sets.end(); it++ ) {
			if ( (*it).second.size() < 2 ) continue; // don't print 1-atom patches
			patch_area = 0.0; // reset the patch area for each patch
			TR_HIG << "atoms in patch: [ ";

			/*std::cout << "(*it).second: [ ";
			for ( Size ii = 1; ii <= (*it).second.size(); ++ii ) {
				std::cout << (*it).second[ ii ] << ", ";
			}
			std::cout << "]" << std::endl;*/

			// (*it).second has all the djs nodes that are in this connected component. we have to translate the djs node
			// into IG node and residue/atom type. djs_id_2_hphobe_index_ provides a mapping from djs node id to a exposed_hydrophobic_data
			// struct. the first struct field, fc_bg_, is a 1 for first-class node and 2 for bgnode.
			// the second field is node_index_ and the third is the atom: exhphobe_index_

			for ( Size ii = 1; ii <= (*it).second.size(); ++ii ) {
				Size const node_index = djs_id_2_hphobe_index_[ (*it).second[ ii ] ].node_index_;
				Size const exphobe_index = djs_id_2_hphobe_index_[ (*it).second[ ii ] ].exhphobe_index_;
				//std::cout << "node_index: " << node_index << ", exphobe_index: " << exphobe_index;

				if ( djs_id_2_hphobe_index_[ (*it).second[ ii ] ].fc_bg_ == 1 ) { // 1 == FC
					conformation::Residue const & rsd( node_considering_alt_state_ == node_index ?
						*get_hpatch_node( node_index )->get_rotamer( alt_state_being_considered_ ) :
						*get_hpatch_node( node_index )->get_rotamer( get_hpatch_node( node_index )->get_current_state() ) );

					// can't use node_considering_alt_state here; have to use the fc_nodes_near_rotsub_bool to see if a particular
					// IG node is near the position considering a substitution
					//utility::vector1< Size > const & exhphobes( node_considering_alt_state_ == node_index ?
					//	get_hpatch_node( node_index )->alt_state_exp_hphobes() : get_hpatch_node( node_index )->curr_state_exp_hphobes() );
					utility::vector1< Size > const & exhphobes( fc_nodes_near_rotsub_bool_[ node_index ] ?
						get_hpatch_node( node_index )->alt_state_exp_hphobes() : get_hpatch_node( node_index )->curr_state_exp_hphobes() );

					/*std::cout << ", exhphobes: [ ";
					for ( Size aa=1; aa <= exhphobes.size(); ++aa ) { std::cout << exhphobes[ aa ] << ", "; }
					std::cout << "]" << std::endl;*/
					if ( exphobe_index > exhphobes.size() ) {
						std::cout << "node_considering_alt_state_: " << node_considering_alt_state_ << std::endl;
						std::cout << "node_considering_alt_state_ == node_index: " << ( node_considering_alt_state_ == node_index ? "yes" : "no" ) << std::endl;
						std::cout << "fc_nodes_near_rotsub_bool_[ ii ]: " << fc_nodes_near_rotsub_bool_[ node_index ] << std::endl;

						get_hpatch_node( node_index )->print(); // will print both current and alt state dots

						utility::vector1< Size > const & curr_state_exhphobes = get_hpatch_node( node_index )->curr_state_exp_hphobes();
						std::cout << "curr state exp hphobes: [ ";
						for ( Size aa=1; aa <= curr_state_exhphobes.size(); ++aa ) { std::cout << curr_state_exhphobes[ aa ] << ", "; }
						std::cout << "], ";
						utility::vector1< Size > const & alt_state_exhphobes = get_hpatch_node( node_index )->alt_state_exp_hphobes();
						std::cout << "alt state exp hphobes: [ ";
						for ( Size aa=1; aa <= alt_state_exhphobes.size(); ++aa ) { std::cout << alt_state_exhphobes[ aa ] << ", "; }
						std::cout << "]" << std::endl;
					}

					Size atom_index = exhphobes[ exphobe_index ];

					assert( rsd.atom_type( atom_index ).element() == carbon_atom || rsd.atom_type( atom_index ).element() == sulfur_atom );
					TR_HIG << rotamer_sets().moltenres_2_resid( node_index ) << "/" << utility::trim( rsd.atom_name( atom_index ) ) << " + ";
					patch_area += ep_sasa_for_djs_node_[ (*it).second[ ii ] ];

				} else { // 2 == BG
					conformation::Residue const & rsd = pose().residue( bgenumeration_2_resid_[ node_index ] );

					utility::vector1< Size > const & exhphobes( get_hpatch_bg_node( node_index )->alt_state_exp_hphobes() );
					//utility::vector1< Size > const & exhphobes( bg_nodes_near_rotsub_bool_[ node_index ] ?
					//	get_hpatch_bg_node( node_index )->alt_state_exp_hphobes() : get_hpatch_bg_node( node_index )->curr_state_exp_hphobes() );

					/*std::cout << ", exhphobes: [ ";
					for ( Size aa=1; aa <= exhphobes.size(); ++aa ) { std::cout << exhphobes[ aa ] << ", "; }
					std::cout << "]" << std::endl;*/

					Size atom_index = exhphobes[ exphobe_index ];

					assert( rsd.atom_type( atom_index ).element() == carbon_atom || rsd.atom_type( atom_index ).element() == sulfur_atom	);
					TR_HIG << bgenumeration_2_resid_[ node_index ] << "/" << utility::trim( rsd.atom_name( atom_index ) ) << " + ";
					patch_area += ep_sasa_for_djs_node_[ (*it).second[ ii ] ];
				}
			}

			Real score = 0.0;
			if ( patch_area > pack::interaction_graph::SurfacePotential::MAX_HPATCH_AREA ) {
				score = hpatch_score_weight_ * pack::interaction_graph::SurfacePotential::MAX_HPATCH_SCORE;
			} else {
				score = hpatch_score_weight_ * pack::interaction_graph::SurfacePotential::get_instance()->hpatch_score( patch_area );
			}
			TR_HIG << "], patch_area: " << patch_area << ", score: " << score << std::endl;

		}
		TR_HIG << std::endl;
		TR_HIG.flush();
#endif

	return total_alt_state_hpatch_score;

}


///
/// @begin HPatchInteractionGraph< V, E, G >::decide_procrastinate_hpatch_computations()
///
/// @detailed
/// Makes the decision whether or not to procrastinate calculating the hpatch score. Basically, if the PD energy got better (dE < 0)
/// then return false so we don't procrastinate the calculation (because the alternate state probably will be accepted?). If the best
/// guess for the hpatch deltaE also comes back better (dE < 0), then return false.  Finally, if the difference between the deltaE
/// for the PD terms and the (guessed) hpatch deltaE is greater than the threshold, return true so we do procrastinate. So basically
/// if the energy (especially the PD energy) gets worse, procrastinate.  Otherwise, don't.
///
template < typename V, typename E, typename G >
bool HPatchInteractionGraph< V, E, G >::decide_procrastinate_hpatch_computations( Real const pd_deltaE, Real const threshold ) const {

	Real hpatch_deltaE_max = 0;

	if ( ! observed_sufficient_hpatch_E_to_predict_min_ )
		return false;

	if ( threshold < 0 || pd_deltaE < 0 )
		return false;

	hpatch_deltaE_max += hpatch_energy_current_state_assignment_ - hpatch_score_min_last_100_;

	// pd_deltaE must be positive, hpatch_deltaE must also be positive.
	// threshold of 5 means, if the PD got more than 5 worse than the guessed hpatch deltaE (e.g. 10 - 3(?) = 7 > 5), procrastinate
#ifdef FILE_DEBUG
	TR_HIG << "decide_procrastinate_hpatch_computations(): pd_deltaE: " << pd_deltaE << ", hpatch_deltaE_max: " << hpatch_deltaE_max
			<< ", threshold: " << threshold << std::endl;
#endif
	if ( (pd_deltaE - hpatch_deltaE_max) > threshold ) {
		return true;
	}
	return false;

}


///
/// @begin HPatchInteractionGraph::reset_from_previous_deltaHpatch_comp
///
/// @brief
/// Iterates through all FCNodes and BGNodes affected by the last deltaHpatch computation, and resets their
/// state.
///
/// @detailed
/// the Node consider() function is the main entry point from the HIG when the annealer considers a sub. need to make sure that the
/// node's alt_state is the same as the current state before we start incrementing/decrementing things or otherwise you wind
/// up with total SASAs that are wrong. alt state will be the same as current state on the first time through because of
/// prep for simA call, but not necessarily the case on following substitutions coming from annealer. any given sub only sets the alt state
/// SASA on that particular (changing) node.  all of the other nodes which were called on to consider a substitution and
/// whose SASAs changed need to reset their alt state dots also (or, otherwise, it appears as if the substitution was
/// committed).  More about this problem can be read in the comments for commit_considered_substitution()
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::reset_from_previous_deltaHpatch_comp() {

	/// This reset is necessary only if the last hpatch deltaE was calculated
	///if ( some_node_in_state_0_ ) return;

#ifdef FILE_DEBUG
	if ( ! some_node_in_state_0_ ) {
		TR_HIG << "reset_from_previous_deltaHpatch_comp: calling reset_alt_state_dots on all FC nodes." << std::endl;
		TR_HIG << "reset_from_previous_deltaHpatch_comp: FC nodes affected by substitution previously considered include: ";
		for ( Size ii = 1; ii <= fc_nodes_near_rotsub_.size(); ++ii ) { TR_HIG << fc_nodes_near_rotsub_[ ii ] << " "; } TR_HIG << std::endl;
	} else {
		TR_HIG << "reset_from_previous_deltaHpatch_comp: unnecessary because not all nodes are in an assigned state." << std::endl;
	}
#endif

	for ( Size ii = 1; ii <= fc_nodes_near_rotsub_.size(); ++ii ) {
		Size const ii_fc_node = fc_nodes_near_rotsub_[ ii ];
		assert( fc_nodes_near_rotsub_bool_[ ii_fc_node ] );
		get_hpatch_node( ii_fc_node )->reset_alt_state_dots();
	}

	if ( ! some_node_in_state_0_ ) {
		for ( Size ii = 1; ii <= fc_nodes_near_rotsub_.size(); ++ii ) {
			Size const ii_fc_node = fc_nodes_near_rotsub_[ ii ];
			fc_nodes_near_rotsub_bool_[ ii_fc_node ] = false;
		}
		fc_nodes_near_rotsub_.clear();
	}

#ifdef FILE_DEBUG
	if ( ! some_node_in_state_0_ ) {
		TR_HIG << "reset_from_previous_deltaHpatch_comp: calling reset_alt_state_dots on all BG nodes." << std::endl;
		TR_HIG << "reset_from_previous_deltaHpatch_comp: BG nodes affected by substitution previously considered include: ";
		for ( Size ii = 1; ii <= bg_nodes_near_rotsub_.size(); ++ii ) { TR_HIG << bg_nodes_near_rotsub_[ ii ] << " "; } TR_HIG << std::endl;
	} else {
		TR_HIG << "reset_from_previous_deltaHpatch_comp: unnecessary because not all nodes are in an assigned state." << std::endl;
	}
#endif

	for ( Size ii = 1; ii <= bg_nodes_near_rotsub_.size(); ++ii ) {
		Size const ii_bg_node = bg_nodes_near_rotsub_[ ii ];
		assert( bg_nodes_near_rotsub_bool_[ ii_bg_node ] );
		get_hpatch_bg_node( ii_bg_node )->reset_alt_state_dots();
	}

	if ( ! some_node_in_state_0_ ) {
		for ( Size ii = 1; ii <= bg_nodes_near_rotsub_.size(); ++ii ) {
			Size const ii_bg_node = bg_nodes_near_rotsub_[ ii ];
			bg_nodes_near_rotsub_bool_[ ii_bg_node ] = false;
		}
		bg_nodes_near_rotsub_.clear();
	}
}


///
/// @begin HPatchInteractionGraph::commit_considered_substitution
///
/// @brief
/// Commits the substitution that the sim annealer had previously asked the graph to consider.  Returns the accurate total energy for the graph.
///
template < typename V, typename E, typename G >
core::PackerEnergy HPatchInteractionGraph< V, E, G >::commit_considered_substitution() {

	// std::cerr << "Committing substitution at node " << node_considering_alt_state_ << std::endl;

#ifdef FILE_DEBUG
	TR_HIG << "commit_considered_substitution(): committing sub on node " << node_considering_alt_state_ << std::endl;
#endif

	core::PackerEnergy hpatch_deltaE = 0.0;
	if ( ! calculated_hpatch_deltaE_ ) {
		hpatch_deltaE = calculate_hpatch_deltaE(); // updates all the Nodes and calculates the hpatch score delta
		++num_hpatch_comps_later_made_;
	}

	// the call to calculate_PD_deltaE_for_substitution() always happens so parent::get_alt_pd_energy_total()
	// will return the correct value. to get the right deltaE we have to add the hpatch deltaE to the
	// PD terms deltaE. deltaE_for_substitution_ is a class member Real variable.
	deltaE_for_substitution_ = get_hpatch_node( node_considering_alt_state_ )->get_pd_energy_delta() + hpatch_deltaE;

	// have to get the pd deltaE before making this call
	get_hpatch_node( node_considering_alt_state_ )->commit_considered_substitution();

	total_energy_current_state_assignment_ = total_energy_current_state_assignment_ + deltaE_for_substitution_;
	hpatch_energy_current_state_assignment_ = hpatch_energy_alternate_state_assignment_;

	node_considering_alt_state_ = -1;
	++num_commits_since_last_update_;

	if ( num_commits_since_last_update_ == COMMIT_LIMIT_BETWEEN_UPDATES ) {
		update_internal_energy_totals_hpatch();
	}

	track_hpatch_E_min();

	return total_energy_current_state_assignment_;
}


///
/// @begin HPatchInteractionGraph< V, E, G >::track_hpatch_E_min
///
/// @brief
/// Keeps track of the minimum hpatch score seen.  Every 100 substitutions, updates the variable hpatch_score_min_last_100.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::track_hpatch_E_min() {

	++num_substitutions_since_hpatch_min_update_;

	Real alt_hpatchE = hpatch_energy_current_state_assignment_;

	if ( hpatch_score_min_recent_ > alt_hpatchE )
		hpatch_score_min_recent_ = alt_hpatchE;

	if ( num_substitutions_since_hpatch_min_update_ == 100 ) { // only update min every 100 calls to track_hpatchE_min (aka commits)
		hpatch_score_min_last_100_ = hpatch_score_min_recent_;
		if ( hpatch_energy_current_state_assignment_ < hpatch_score_min_last_100_ )
			hpatch_score_min_last_100_ = hpatch_energy_current_state_assignment_;
		observed_sufficient_hpatch_E_to_predict_min_ = true;
		num_substitutions_since_hpatch_min_update_ = 0;
	}

}


///
/// @begin HPatchInteractionGraph::set_network_state
///
/// @brief
/// Switch the state assignment of every first class node in the graph.
/// Useful, for instance, if you want to switch to the best network state that you've found so far.
///
/// @detailed
/// This function is the last major entry point from the Annealer into the HIG.
///
template < typename V, typename E, typename G >
core::PackerEnergy HPatchInteractionGraph< V, E, G >::set_network_state( FArray1_int & node_states ) {

#ifdef FILE_DEBUG
	TR_HIG << "set_network_state() called with states: " << node_states << std::endl;
#endif
	for ( Size ii = 1; ii <= (Size)parent::get_num_nodes(); ++ii ) {
		core::PackerEnergy deltaE = 0.0;
		core::PackerEnergy previousE = 0.0;
		consider_substitution( ii, node_states( ii ), deltaE, previousE ); // might get procrastinated but that's ok
		commit_considered_substitution(); // it will get updated correctly here
	}

	// not necessary because the code above goes through the consider/commit process which updates the internal
	// energy totals.
 	//update_internal_energy_totals_hpatch();

	return total_energy_current_state_assignment_;

}


///
/// @begin HPatchInteractionGraph::get_energy_current_state_assignment
///
/// @brief
/// returns the energy of the entire graph under the current network state assignment.  Also sends a bunch of information to standard error.
/// Only seems to be called by the MultiCoolAnnealer.
///
template < typename V, typename E, typename G >
core::PackerEnergy HPatchInteractionGraph< V, E, G >::get_energy_current_state_assignment() {
	return total_energy_current_state_assignment_;
}



///
/// @begin HPatchInteractionGraph::print_current_state_assignment
///
/// @brief
/// Should write the state assigned to each first class vertex to the screen.
///
/*template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::print_current_state_assignment() const {

	// print out the one-body and hpatch energies for all first class nodes
	TR_HIG << "internal energies: " << std::endl;
	for (int ii = 1; ii <= parent::get_num_nodes(); ++ii) {
		Real one_body = get_hpatch_node( ii )->get_curr_state_one_body_energy();
		TR_HIG << "node " << ii << " 1b: " << one_body;
		Real sasa = get_hpatch_node( ii )->get_current_state_sasa();
		TR_HIG << ", sasa = " << sasa;

		if ( ii % 3 == 0) {
			TR_HIG << std::endl;
		}
	}

	TR_HIG << std::endl;

	// print out the hpatch energies for all background nodes
	//for ( Size ii = 1; ii <= (Size)parent::get_num_background_nodes(); ++ii ) {
	//	Real bg_sasa = get_hpatch_bg_node( ii )->get_current_sasa();
	//	TR_HIG << "bg res: " << bgenumeration_2_resid_[ ii ] << " sasa: " << bg_sasa << std::endl;
	//}

	// print out the two-body energies for all edges between first-class nodes only?
	int count_edges = 0;
	for (std::list< core::pack::interaction_graph::EdgeBase*>::const_iterator iter = parent::get_edge_list_begin(); iter != parent::get_edge_list_end(); ++iter) {
		Real edge_energy = ((HPatchEdge< V, E, G >*) (*iter))->get_current_two_body_energy();
		TR_HIG << "edge: " << edge_energy << " ";

		if ( count_edges % 5 == 0)
			TR_HIG << std::endl;
		++count_edges;
	}

}*/


///
/// @begin HPatchInteractionGraph::get_edge_memory_usage
///
/// @brief
/// Should return a measurement of the memory used by the interaction graph
/// to store the rotamer pair energies.  Unimplemented.
///
template < typename V, typename E, typename G >
int HPatchInteractionGraph< V, E, G >::get_edge_memory_usage() const {
	return 0;
}


///
/// @begin HPatchInteractionGraph::count_static_memory
///
template < typename V, typename E, typename G >
unsigned int HPatchInteractionGraph< V, E, G >::count_static_memory() const {
	return sizeof ( HPatchInteractionGraph< V, E, G > );
}


///
/// @begin HPatchInteractionGraph::count_dynamic_memory
///
template < typename V, typename E, typename G >
unsigned int HPatchInteractionGraph< V, E, G >::count_dynamic_memory() const {
	unsigned int total_memory = parent::count_dynamic_memory();

	total_memory += resid_2_bgenumeration_.size() * sizeof( Size );
	total_memory += bgenumeration_2_resid_.size() * sizeof( Size );

	return total_memory;
}


///
/// @begin HPatchInteractionGraph::get_energy_sum_for_vertex_group
///
/// @brief
/// returns the sum of the PD energy and the hpatch energy for all members first class members of a user-defined
/// vertex subset.  Unimplemented.
///
template < typename V, typename E, typename G >
Real HPatchInteractionGraph< V, E, G >::get_energy_sum_for_vertex_group( Size ) {
	//apl functionality stubbed out for now
	return 0;
}


///
/// @begin HPatchInteractionGraph::print_internal_energies_for_current_state_assignment
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::print_internal_energies_for_current_state_assignment() {

	// print out the one-body and hpatch energies for all first class nodes
	TR_HIG << "internal energies: " << std::endl;
	for ( Size ii = 1; ii <= parent::get_num_nodes(); ++ii ) {
		Real one_body = get_hpatch_node( ii )->get_curr_state_one_body_energy();
		TR_HIG << "node " << ii << " 1b: " << one_body;
		Real sasa = get_hpatch_node( ii )->get_current_state_sasa();
		TR_HIG << ", sasa = " << sasa;

		if ( ii % 3 == 0) {
			TR_HIG << std::endl;
		}
	}
	TR_HIG << std::endl;

	// print out the hpatch energies for all background nodes
	for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii) {
		Real bg_sasa = get_hpatch_bg_node( ii )->get_current_sasa();
		TR_HIG << "bg res: " << bgenumeration_2_resid_[ ii ] << " sasa: " << bg_sasa << std::endl;
	}

	// print out the two-body energies for all edges between first-class nodes only?
	int count_edges = 0;
	for (std::list< core::pack::interaction_graph::EdgeBase*>::const_iterator iter = parent::get_edge_list_begin(); iter != parent::get_edge_list_end(); ++iter) {
		Real edge_energy = ((HPatchEdge< V, E, G >*) (*iter))->get_current_two_body_energy();
		TR_HIG << "edge: " << edge_energy << " ";

		if ( count_edges % 5 == 0)
			TR_HIG << std::endl;
		++count_edges;
	}
}


///
/// @begin HPatchInteractionGraph::write_dot_kinemage
///
/*template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::write_dot_kinemage( std::ofstream & output_kin ) {

	output_kin << "@group {dots} off" << std::endl;
	output_kin << "@subgroup {molten_residues} dominant" << std::endl;

	for (int ii = 1; ii <= parent::get_num_nodes(); ++ii) {
		get_hpatch_node( ii )->write_dot_kinemage( output_kin );
	}

	output_kin << "@subgroup {background_residues} dominant" << std::endl;
	for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii) {
		get_hpatch_bg_node( ii )->write_dot_kinemage( output_kin );
	}

}*/


///
/// @begin HPatchInteractionGraph::print
///
/// @brief
/// useful for debugging
///
template< typename V, typename E, typename G >
void
HPatchInteractionGraph<V, E, G>::print() const {

	std::cout << "HPatch Interaction Graph state: " << std::endl;
	std::cout << "nodes: " << std::endl;
	for (int jj = 1; jj <= parent::get_num_nodes(); ++jj) {
		get_hpatch_node( jj )->print();
	}

	std::cout << "bgnodes: " << std::endl;
	for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii) {
		get_hpatch_bg_node( ii )->print();
	}
}


// The below functions are only used for the unit tests. However, since most developers these days are running the unit
// tests using release mode, I can't #ifdef these functions (to leave them out of release mode builds) or the unit
// tests don't compile. So just compile them regardless of build mode.

/// @begin HPatchNode::get_current_state_rotamer_dots
///
/// @brief
/// Returns current state. Only used by the unit tests.
///
template< typename V, typename E, typename G >
RotamerDots const & HPatchNode<V, E, G>::get_current_state_rotamer_dots() { return current_state_rotamer_dots_; }

/// @begin HPatchNode::get_alt_state_rotamer_dots
///
/// @brief
/// Returns current state. Only used by the unit tests.
///
template< typename V, typename E, typename G >
RotamerDots const & HPatchNode<V, E, G>::get_alt_state_rotamer_dots() { return alt_state_rotamer_dots_; }


/// @begin HPatchBackgroundNode::get_current_state_rotamer_dots
///
/// @brief
/// Returns current state. Only used by the unit tests.
///
template< typename V, typename E, typename G >
RotamerDots const & HPatchBackgroundNode<V, E, G>::get_current_state_rotamer_dots() { return current_state_rotamer_dots_; }


/// @begin HPatchBackgroundNode::get_alt_state_rotamer_dots
///
/// @brief
/// Returns current state. Only used by the unit tests.
///
template< typename V, typename E, typename G >
RotamerDots const & HPatchBackgroundNode<V, E, G>::get_alt_state_rotamer_dots() { return alt_state_rotamer_dots_; }


///
/// @begin HPatchInteractionGraph::get_network_state
///
/// @brief
/// Returns the state on each FCNode, but not necessarily in pose resid order. Only used by the unit tests.
///
template< typename V, typename E, typename G >
std::vector< int > HPatchInteractionGraph<V, E, G>::get_network_state() const {

	std::vector< int > networkstate;
	for ( int jj = 1; jj <= parent::get_num_nodes(); ++jj ) {
		networkstate.push_back( get_hpatch_node(jj)->get_current_state() );
	}
	return networkstate;
}


///
/// @begin HPatchInteractionGraph::set_observed_sufficient_boolean_true
///
/// @brief
/// Sets the observed_sufficient_hpatch_E_to_predict_min_ to true. Only used by the unit tests.
///
template< typename V, typename E, typename G >
void HPatchInteractionGraph<V, E, G>::set_observed_sufficient_boolean_true() {
	observed_sufficient_hpatch_E_to_predict_min_ = true;
}


///
/// @begin HPatchInteractionGraph< V, E, G >::get_all_sasas
///
/// @brief
/// Iterates over all nodes and bgnodes
/// brute-force recounting.
///
template < typename V, typename E, typename G >
void HPatchInteractionGraph< V, E, G >::get_all_sasas( utility::vector1< Real > & node_sasas, utility::vector1< Real > & bgnode_sasas ) {

	for ( Size ii=1; ii <= (Size)parent::get_num_nodes(); ++ii ) {
		node_sasas[ ii ] = ((get_hpatch_node(ii))->get_current_state_rotamer_dots()).get_sasa();
	}

	for ( Size ii=1; ii <= (Size)parent::get_num_background_nodes(); ++ii ) {
		bgnode_sasas[ ii ] = ((get_hpatch_bg_node(ii))->get_current_state_rotamer_dots()).get_sasa();
	}

}


///
/// @begin HPatchInteractionGraph::bg_node_2_resid
///
/// @brief
/// Provides read access to the bg to resid array. Returns -1 if the index is not in bounds.
///
template < typename V, typename E, typename G >
int HPatchInteractionGraph< V, E, G >::bg_node_2_resid( Size node_index ) {

	if ( node_index > num_residues_assigned_as_background_ ) {
		utility_exit_with_message( "Out of bounds array index passed to bg_node_2_resid. Quitting." );
	}
	return bgenumeration_2_resid_[ node_index ];
}


/*///
/// @begin HPatchInteractionGraph::get_fc_nodes_near_rotsub
///
/// @brief
/// Read access to the vector fc_nodes_near_rotsub_.
///
template < typename V, typename E, typename G >
utility::vector1< Size > const & HPatchInteractionGraph< V, E, G >::get_fc_nodes_near_rotsub() {
	return fc_nodes_near_rotsub_;
}

///
/// @begin HPatchInteractionGraph::get_fc_nodes_near_rotsub_bool
///
/// @brief
/// Read access to the vector fc_nodes_near_rotsub_bool_.
///
template < typename V, typename E, typename G >
utility::vector1< bool > const & HPatchInteractionGraph< V, E, G >::get_fc_nodes_near_rotsub_bool() {
	return fc_nodes_near_rotsub_bool_;
}

///
/// @begin HPatchInteractionGraph::get_bg_nodes_near_rotsub
///
/// @brief
/// Read access to the vector bg_nodes_near_rotsub_.
///
template < typename V, typename E, typename G >
utility::vector1< Size > const & HPatchInteractionGraph< V, E, G >::get_bg_nodes_near_rotsub() {
	return bg_nodes_near_rotsub_;
}

///
/// @begin HPatchInteractionGraph::get_bg_nodes_near_rotsub_bool
///
/// @brief
/// Read access to the vector bg_nodes_near_rotsub_bool_.
///
template < typename V, typename E, typename G >
utility::vector1< bool > const & HPatchInteractionGraph< V, E, G >::get_bg_nodes_near_rotsub_bool() {
	return bg_nodes_near_rotsub_bool_;
}


///
/// @begin HPatchInteractionGraph::get_fc_exp_hphobe_djs_offsets
///
/// @brief
/// Read access to the vector fc_exp_hphobe_djs_offsets_.
///
template < typename V, typename E, typename G >
utility::vector1< Size > const & HPatchInteractionGraph< V, E, G >::get_fc_exp_hphobe_djs_offsets() {
	return fc_exp_hphobe_djs_offsets_;
}


///
/// @begin HPatchInteractionGraph::get_bg_exp_hphobe_djs_offsets
///
/// @brief
/// Read access to the vector bg_exp_hphobe_djs_offsets_.
///
template < typename V, typename E, typename G >
utility::vector1< Size > const & HPatchInteractionGraph< V, E, G >::get_bg_exp_hphobe_djs_offsets() {
	return bg_exp_hphobe_djs_offsets_;
}


///
/// @begin HPatchInteractionGraph::get_fc_n_exp_hphobes
///
/// @brief
/// Read access to the vector fc_n_exp_hphobes_.
///
template < typename V, typename E, typename G >
utility::vector1< Size > const & HPatchInteractionGraph< V, E, G >::get_fc_n_exp_hphobes() {
	return fc_n_exp_hphobes_;
}

///
/// @begin HPatchInteractionGraph::get_bg_n_exp_hphobes
///
/// @brief
/// Read access to the vector bg_n_exp_hphobes_.
///
template < typename V, typename E, typename G >
utility::vector1< Size > const & HPatchInteractionGraph< V, E, G >::get_bg_n_exp_hphobes() {
	return bg_n_exp_hphobes_;
}*/



} //end namespace
} //end namespace
} //end namespace

#endif

