// -*- Mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/SurfaceInteractionGraph.hh
/// @brief  Interaction graph which implements a non-PD, environment-dependent score for surface residues
/// @author Ron Jacak (ron.jacak@gmail.com)

#ifndef INCLUDED_core_pack_interaction_graph_SurfaceInteractionGraph_hh
#define INCLUDED_core_pack_interaction_graph_SurfaceInteractionGraph_hh

//Rosetta Headers
#include <core/pack/interaction_graph/AdditionalBackgroundNodesInteractionGraph.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>
#include <core/pack/interaction_graph/SurfaceInteractionGraph.fwd.hh>
#include <core/pack/interaction_graph/SurfacePotential.hh>

#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/TenANeighborGraph.hh>

//Utility Headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>

//ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1.io.hh>

//C++ Headers
#include <vector>
// AUTO-REMOVED #include <typeinfo> //required by GCC 4.3.2

#include <utility/vector1.hh>


/// Tracer instance for this file
static thread_local basic::Tracer TR_NODE( "core.pack.surfaceig.node" );
static thread_local basic::Tracer TR_EDGE( "core.pack.surfaceig.edge" );
static thread_local basic::Tracer TR_BGNODE( "core.pack.surfaceig.bgnode" );
static thread_local basic::Tracer TR_BGEDGE( "core.pack.surfaceig.bgedge" );
static thread_local basic::Tracer TR_SIG( "core.pack.surfaceig.sig" );
static thread_local basic::Tracer TR_STATS( "core.pack.surfaceig.stats" );

//#define DOUBLE_CHECK_COUNTS 1
//#define FILE_DEBUG 1

namespace core {
namespace pack {
namespace interaction_graph {

template < typename V, typename E, typename G > class SurfaceNode;
template < typename V, typename E, typename G > class SurfaceBackgroundNode;
template < typename V, typename E, typename G > class SurfaceEdge;
template < typename V, typename E, typename G > class SurfaceBackgroundEdge;
template < typename V, typename E, typename G > class SurfaceInteractionGraph;


//----------------------------------------------------------------------------//
//---------------------------- Surface Node Class ----------------------------//
//----------------------------------------------------------------------------//

///
/// @begin SurfaceNode
///
/// @brief
/// Defines a FirstClass node which will keep track of changes in the surface energy.
/// FirstClassNode is defined and implemented in AdditionalBackgroundNodesInteractionGraph.
///
/// @remarks
/// No public default constructor makes this class uncopyable.
///
template < typename V, typename E, typename G >
class SurfaceNode : public FirstClassNode< V, E, G > {

	public:
		typedef FirstClassNode< V, E, G > parent;

	public:
		SurfaceNode( G* owner, int node_index, int num_states );
		virtual ~SurfaceNode();

		virtual void assign_zero_state();
		virtual bool state_unassigned() const { return parent::get_current_state() == 0; }

		void assign_state_surface( int state );
		Real get_curr_state_surface_energy() const;
		Real project_deltaE_for_substitution_surface( int alternate_state, core::PackerEnergy & prev_energy_for_node, float deltaE_thresh_for_avoiding_surface_calcs );

		Real get_surface_deltaE_for_neighbors_state_substitution( SurfaceNode<V,E,G>* node_considering_substitution, int changing_nodes_curr_state, int changing_nodes_alt_state );

		Real commit_considered_substitution_surface();
		void acknowledge_neighbors_substitution_surface();

		bool detect_neighborship_with_node( int node_id, bool first_class ) const;

		static void print_surface_avoidance_stats();
		static void reset_surface_avoidance_stats();

		// virtual methods from NodeBase class
		virtual void print() const;

		virtual void prepare_for_simulated_annealing();
		virtual unsigned int getMemoryUsageInBytes() const;
		virtual unsigned int count_static_memory() const;
		virtual unsigned int count_dynamic_memory() const;

		// setter for the rotamers object.
		void set_rotamers( rotamer_set::RotamerSetCOP rotamers );
		conformation::ResidueCOP get_rotamer( int state ) const;

		// hold on to these so we can look up is_hydrophobic on state changes
		rotamer_set::RotamerSetCOP rotamer_set_;
		utility::vector1< conformation::ResidueCOP > rotamers_vector_;

		inline
		int wt_seqpos_for_node() const {
			return get_surface_owner()->rotamer_sets().moltenres_2_resid( parent::get_node_index() );
		}
		inline
		conformation::Residue wt_residue_for_node() const {
			return get_surface_owner()->pose().residue( (get_surface_owner()->rotamer_sets().moltenres_2_resid( parent::get_node_index() )) );
		}


		Real get_surface_score_difference() const;

		bool is_surface_exposed() const;
		void surface_exposed( bool value );
		bool is_below_buried_residue_no_hsasa_cutoff() const;
		void is_below_buried_residue_no_hsasa_cutoff( bool value );

		void reset_alt_state_total_hASA();

		void initialize_num_neighbors_counting_self() const;
		int num_neighbors_counting_self() const;

		void init_hASA_variables();
		Real calculate_amount_total_hydrophobic_ASA();
		Real average_residue_hASA() const;
		Real average_residue_hASA( chemical::AA residue_type, Size num_nbs ) const;
		Real hASA_energy( Real patch_area ) const;

		void verify_patch_areas_correct( int node_id, int previous_state, Real previous_state_hASA );

		// Extra methods only used only for the unit tests.
		void set_observed_sufficient_boolean_true();
		std::vector<Real> get_hASA_for_node_and_nbs();
		std::vector<Real> get_alt_state_hASA_for_node_and_nbs();
		Real get_current_hASA();
		Real get_alt_hASA();

	protected:
		inline
		SurfaceEdge< V, E, G >* get_incident_surface_edge( int index ) const {
			return (SurfaceEdge< V, E, G >*) parent::get_incident_edge( index );
		}

		inline
		SurfaceBackgroundEdge< V, E, G >* get_edge_to_surface_bg_node( int index ) const {
			return (SurfaceBackgroundEdge< V, E, G >*) parent::get_edge_to_bg_node( index );
		}

		inline
		SurfaceInteractionGraph< V, E, G >* get_surface_owner() const {
			return (SurfaceInteractionGraph< V, E, G >*) parent::get_owner();
		}

	private:
		Real project_surface_deltaE();
		void track_surface_E_min();

		bool decide_procrastinate_surface_computations( Real const pd_deltaE, Real const threshold ) const;

		bool calculated_surface_deltaE_;
		Real deltaE_for_substitution_;

		Real curr_state_total_hASA_;
		Real alt_state_total_hASA_;

		bool have_prepared_for_simA_;

		Real surface_score_min_last_100_;
		Real surface_score_min_recent_;
		int num_substitutions_since_surface_min_update_;
		bool observed_sufficient_surface_E_to_predict_min_;

		bool surface_exposed_;
		bool is_below_buried_residue_no_hsasa_cutoff_;
		std::map< std::pair<int, int >, int > fc_neighbor_map;
		std::map< std::pair<int, int >, int > bg_neighbor_map;

		mutable int num_neighbors_counting_self_;

		static Real surface_energy_weight_;
		static const int MAX_SURFACE_ENERGY;
		static const int INTERACTION_RADIUS = 10;
		static const int SURFACE_EXPOSED_CUTOFF = 20;
		static const int BURIED_RESIDUE_NO_HSASA_CUTOFF = 24;
		static const int MAX_PATCH_SURFACE_AREA = 1100;

		static int num_state_substitutions_considered_;
		static int num_surface_comps_procrastinated_;
		static int num_surface_comps_later_made_;

		//no default constructor, uncopyable
		SurfaceNode();
		SurfaceNode( SurfaceNode< V, E, G > const & );
		SurfaceNode< V, E, G > & operator = ( SurfaceNode< V, E, G > const & );


};


//----------------------------------------------------------------------------//
//------------------- Surface Background Node Class -----------------------//
//----------------------------------------------------------------------------//

///
/// @begin SurfaceBackgroundNode
///
/// @brief
/// Defines a BackgroundResidue node which will contribute to changes in surface energy
/// due to state changes on neighboring nodes, and not because of state changes to it.
///
template < typename V, typename E, typename G >
class SurfaceBackgroundNode : public BackgroundNode< V, E, G > {

	public:
		typedef BackgroundNode< V, E, G > parent;

	public:

		SurfaceBackgroundNode( AdditionalBackgroundNodesInteractionGraph< V, E, G > * owner, int node_index );
		virtual ~SurfaceBackgroundNode();

		bool detect_neighborship( SurfaceNode< V, E, G >* node ) const;

		Real project_surface_deltaE_for_substitution( SurfaceNode< V, E, G >* fc_node_changing,
				int changing_nodes_curr_state, int changing_nodes_alt_state );

		void acknowledge_substitution_surface();
		Real get_surface_score() const;

		virtual void prepare_for_simulated_annealing();
		void print() const;

		virtual unsigned int count_static_memory() const;
		virtual unsigned int count_dynamic_memory() const;

		inline
		conformation::Residue const & wt_residue_for_node() const {
			return get_surface_owner()->pose().residue( (get_surface_owner()->bg_node_2_resid(parent::get_node_index())) );
		}

		bool is_surface_exposed() const;
		void surface_exposed( bool value );
		bool is_below_buried_residue_no_hsasa_cutoff() const;
		void is_below_buried_residue_no_hsasa_cutoff( bool value );

		void reset_alt_state_total_hASA();

		void initialize_num_neighbors_counting_self() const;
		int num_neighbors_counting_self() const;

		void init_hASA_variables();
		Real calculate_amount_total_hydrophobic_ASA();
		Real average_residue_hASA() const;
		Real average_residue_hASA( chemical::AA residue_type, Size num_nbs ) const;
		Real hASA_energy( Real patch_area ) const;

		// Only used for the unit tests.
		Real get_current_hASA();
		Real get_alt_hASA();

	protected:
		inline
		SurfaceInteractionGraph< V, E, G >* get_surface_owner() const {
			return (SurfaceInteractionGraph< V, E, G >*) parent::get_owner();
		}

		inline
		SurfaceBackgroundEdge< V, E, G >* get_surface_bg_edge( int index ) {
			return (SurfaceBackgroundEdge< V, E, G >*) parent::get_incident_edge( index );
		}

	private:
		static Real surface_energy_weight_;
		static const int MAX_SURFACE_ENERGY;
		static const int INTERACTION_RADIUS = 10;
		static const int SURFACE_EXPOSED_CUTOFF = 20;
		static const int BURIED_RESIDUE_NO_HSASA_CUTOFF = 24;
		static const int MAX_PATCH_SURFACE_AREA = 1100;

		Real curr_state_total_hASA_;
		Real alt_state_total_hASA_;

		bool have_prepared_for_simA_;

		bool surface_exposed_;
		bool is_below_buried_residue_no_hsasa_cutoff_;
		mutable int num_neighbors_counting_self_;

		//no default constructor, uncopyable
		SurfaceBackgroundNode();
		SurfaceBackgroundNode( SurfaceBackgroundNode< V, E, G > const & );
		SurfaceBackgroundNode< V, E, G > & operator = ( SurfaceBackgroundNode< V, E, G > const & );

};


//----------------------------------------------------------------------------//
//-------------------------- Surface Edge Class ---------------------------//
//----------------------------------------------------------------------------//

///
/// @begin SurfaceEdge
///
/// @brief
/// Defines a Surface edge which will be used in determining surface energy.
///
template < typename V, typename E, typename G >
class  SurfaceEdge : public FirstClassEdge< V, E, G > {

	public:
		typedef  FirstClassEdge< V, E, G >  parent;

	public:
		SurfaceEdge( G* owner, int node1, int node2 );
		virtual ~SurfaceEdge();

		void acknowledge_state_zeroed_surface( int node_index );

		Real get_surface_deltaE_for_neighbor( int node_considering_substitution, int alt_state );

		Real get_current_two_body_energy() const;
		void acknowledge_substitution_surface();

		//Virtual methods from EdgeBase
		virtual void declare_energies_final();
		virtual void prepare_for_simulated_annealing();
		virtual unsigned int getMemoryUsageInBytes() const;

		Real get_max_surface_deltaE_guess( int node_changing ) const;

		virtual unsigned int count_static_memory() const;
		virtual unsigned int count_dynamic_memory() const;

		// this method used only for testing
		inline
		void set_max_surface_deltaE() {
			max_surface_deltaE_last_50_commits_[0] = max_surface_deltaE_last_50_commits_[1] = 0.1;
		}

//	protected:
	public:
		// Need to make this method public so that SurfaceNodes can call methods on the SurfaceBackgroundNode connected to them
		// via this Edge.  The class hierarchy/design does not necessitate this; just some debugging checks in the SurfaceNode that need
		// access to the BackgroundNode object.
		inline
		SurfaceNode< V, E, G >* get_surface_node( int index ) {
			return (SurfaceNode< V, E, G >*) E::get_node( index );
		}

	private:
		inline
		void inform_non_changing_node_of_neighbors_change();

		int node_changing_;
		int node_not_changing_;
		int nodes_curr_states_[2];
		int nodes_alt_states_[2];

		Real max_surface_deltaE_last_50_commits_[2];
		Real max_surface_deltaE_recent_50_commits_[2];
		Real magnitude_last_surface_deltaE_[2];
		int num_surface_deltaE_observations_since_update_[2];

		void track_max_magnitude_surface_deltaE();

		//no default constructor, uncopyable
		SurfaceEdge();
		SurfaceEdge( SurfaceEdge< V, E, G > const & );
		SurfaceEdge< V, E, G > & operator = ( SurfaceEdge< V, E, G > const & );
};



//----------------------------------------------------------------------------//
//------------------- Surface Background Edge Class -----------------------//
//----------------------------------------------------------------------------//

///
/// @begin SurfaceBackgroundEdge
///
/// @brief
/// Defines an edge between a FirstClass (SurfaceNode) and a background node (SurfaceBackgroundNode)
///
/// @detailed
/// In addition to implementing the virtual base class methods, this class additionally defines methods
/// relating to keeping track of data relating to surface.
///
template < typename V, typename E, typename G >
class SurfaceBackgroundEdge : public BackgroundToFirstClassEdge< V, E, G > {

	public:
		typedef  BackgroundToFirstClassEdge< V, E, G > parent;

	public:
		SurfaceBackgroundEdge( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int first_class_node_index, int background_node_index );
		virtual ~SurfaceBackgroundEdge();

		void prepare_for_simulated_annealing();

		Real get_surface_deltaE_for_substitution( int alt_state );
		void acknowledge_substitution_surface();
		void acknowledge_state_change( int new_state );
		Real get_max_surface_deltaE_guess() const;

		virtual unsigned int count_static_memory() const;
		virtual unsigned int count_dynamic_memory() const;

		// this function only used for unit tests
		inline
		void set_max_surface_deltaE() {
			max_surface_deltaE_last_50_commits_ = 0.1;
		}

	public:
		// Making this method public also so that a BGNode can call a method on a Node via this Edge.
		inline
		SurfaceNode< V, E, G >* get_surface_node() const {
			return (SurfaceNode< V, E, G >*) parent::get_first_class_node();
		}

	public:
		// Need to make this method public so that SurfaceNodes can call methods on the SurfaceBackgroundNode connected to them
		// via this Edge.  The class hierarchy/design does not necessitate this; just some debugging checks in the SurfaceNode that need
		// access to the BackgroundNode object.
		inline
		SurfaceBackgroundNode< V, E, G >* get_surface_bg_node() const {
			return (SurfaceBackgroundNode< V, E, G >*) parent::get_background_node();
		}

	private:
		int fc_node_curr_state_;
		int fc_node_alt_state_;

		Real max_surface_deltaE_last_50_commits_;
		Real max_surface_deltaE_recent_50_commits_;
		Real magnitude_last_surface_deltaE_;
		int num_surface_deltaE_observations_since_update_;

		void track_max_magnitude_surface_deltaE();

		//no default constructor, uncopyable
		SurfaceBackgroundEdge();
		SurfaceBackgroundEdge( SurfaceBackgroundEdge< V, E, G > const & );
		SurfaceBackgroundEdge< V, E, G > & operator = ( SurfaceBackgroundEdge< V, E, G > const & );

};




//----------------------------------------------------------------------------//
//--------------------- Surface Interaction Graph -------------------------//
//----------------------------------------------------------------------------//

///
/// @begin SurfaceInteractionGraph
///
/// @brief
/// Defines the interaction graph that will keep track of changes to the surface score.
///
/// @detailed
/// In addition to implementing the virtual base class methods, this class additionally defines methods
/// relating to keeping track of data relating to surface.
///
template < typename V, typename E, typename G >
class SurfaceInteractionGraph : public AdditionalBackgroundNodesInteractionGraph< V, E, G > {

	public:
		typedef  AdditionalBackgroundNodesInteractionGraph< V, E, G >  parent;

	public:
		SurfaceInteractionGraph( int num_nodes );
		virtual ~SurfaceInteractionGraph();

		void initialize( rotamer_set::RotamerSetsBase const & rot_sets );

		// Virtual public methods from InteractionGraphBase
		virtual void prepare_for_simulated_annealing();
		virtual void  blanket_assign_state_0();
		virtual core::PackerEnergy set_state_for_node( int node_ind, int new_state );
		virtual core::PackerEnergy set_network_state( ObjexxFCL::FArray1_int& node_states );

		virtual void consider_substitution( int node_ind, int new_state, core::PackerEnergy & delta_energy, core::PackerEnergy & prev_energy_for_node );
		virtual core::PackerEnergy commit_considered_substitution();
		virtual core::PackerEnergy get_energy_current_state_assignment();

		using parent::set_errorfull_deltaE_threshold;

		virtual void set_errorfull_deltaE_threshold( Real deltaE );

		void set_num_residues_in_protein( int num_res );
		void set_num_background_residues( int num_background_residues );
		void set_residue_as_background_residue( int residue );

		void print_internal_energies_for_current_state_assignment();

		virtual int get_edge_memory_usage() const;
		virtual unsigned int count_static_memory() const;
		virtual unsigned int count_dynamic_memory() const;
		void print() const;

		inline
		pose::Pose const & pose() const { return *pose_; }
		void set_pose( pose::Pose const & pose );

		inline
		task::PackerTask const & packer_task() const { return *packer_task_; }
		void set_packer_task( task::PackerTask const & task );

		void set_surface_score_weight( Real weight ) { surface_score_weight_ = weight; }
		Real surface_score_weight() { return surface_score_weight_; }

		inline
		rotamer_set::RotamerSets const & rotamer_sets() const { return *rotamer_sets_; }
		void set_rotamer_sets( rotamer_set::RotamerSets const & rotsets );

		int bg_node_2_resid( int node_index );

		// methods used only by unit tests
		std::vector<int> get_network_state() const;
		void set_observed_sufficient_boolean_true();
		std::vector<Real> get_hASA_for_node_and_nbs(int index);
		std::vector<Real> get_alt_state_hASA_for_node_and_nbs(int index);

	protected:

		//Factory Methods:
		//From InteractionGraphBase
		virtual core::pack::interaction_graph::NodeBase* create_new_node( int node_index, int num_states );
		virtual core::pack::interaction_graph::EdgeBase* create_new_edge( int index1, int index2);

		//From AdditionalBackgroundNodesInteractionGraph
		virtual BackgroundNode< V, E, G >* create_background_node( int node_index );
		virtual BackgroundToFirstClassEdge< V, E, G >* create_background_edge( int fc_node_index, int bg_node_index);

		inline
		SurfaceNode< V, E, G >* get_surface_node( int index ) const {
			return (SurfaceNode< V, E, G >*) G::get_node( index );
		}

		inline
		SurfaceBackgroundNode< V, E, G >* get_surface_bg_node( int index ) const {
			return (SurfaceBackgroundNode< V, E, G >*) parent::get_background_node( index );
		}

		void update_internal_energy_totals_surface();

	private:
		void detect_background_residue_and_first_class_residue_neighbors();
		void blanket_reset_alt_state_total_hASAs();

		int num_total_residues_;
		int num_residues_assigned_as_background_;

		utility::vector1< int > resid_2_bgenumeration_;
		utility::vector1< int > bgenumeration_2_resid_;

		int num_commits_since_last_update_;
		core::PackerEnergy total_energy_current_state_assignment_;
		core::PackerEnergy total_energy_alternate_state_assignment_;
		int node_considering_alt_state_;
		Real deltaE_threshold_for_avoiding_surface_calcs_;
		bool prepared_for_simulated_annealing_;

		static const int COMMIT_LIMIT_BETWEEN_UPDATES = 1024; // 2^10
		static const int SURFACE_EXPOSED_CUTOFF = 20;
		static const int BURIED_RESIDUE_NO_HSASA_CUTOFF = 24;

		pose::PoseOP pose_;
		task::PackerTaskOP packer_task_;
		rotamer_set::RotamerSetsOP rotamer_sets_;
		core::Real surface_score_weight_;

		//no default constructor, uncopyable
		SurfaceInteractionGraph();
		SurfaceInteractionGraph( SurfaceInteractionGraph< V, E, G > const & );
		SurfaceInteractionGraph< V, E, G > & operator = ( SurfaceInteractionGraph< V, E, G > const & );

};


// Begin implementation: must be contained within this header file due to use of templates!


//----------------------------------------------------------------------------//
//-------------------------- Surface Node Class ---------------------------//
//----------------------------------------------------------------------------//

template < typename V, typename E, typename G >
Real SurfaceNode< V, E, G >::surface_energy_weight_ = 1.0;

template < typename V, typename E, typename G >
const int SurfaceNode< V, E, G >::MAX_SURFACE_ENERGY = 100;

template < typename V, typename E, typename G >
int SurfaceNode< V, E, G >::num_state_substitutions_considered_( 0 );

template < typename V, typename E, typename G >
int SurfaceNode< V, E, G >::num_surface_comps_procrastinated_( 0 );

template < typename V, typename E, typename G >
int SurfaceNode< V, E, G >::num_surface_comps_later_made_( 0 );


///
/// @begin SurfaceNode< V, E, G >::SurfaceNode
///
/// @brief
/// SurfaceNode constructor
///
/// @param
/// owner - [in] - the owning interaction graph
/// node_id - [in] - the index for this node amongst its owners set
/// num_states - [in] - the number of states for this node
///
template < typename V, typename E, typename G >
SurfaceNode< V, E, G >::SurfaceNode( G* owner, int node_index, int num_states ) :
	FirstClassNode< V, E, G > ( owner, node_index, num_states ),
	rotamers_vector_( num_states ),
	calculated_surface_deltaE_( false ),
	deltaE_for_substitution_( 0.0f ),
	curr_state_total_hASA_( 0.0 ),
	alt_state_total_hASA_( 0 ),
	have_prepared_for_simA_( false ),
	surface_score_min_last_100_( 0 ),
	surface_score_min_recent_( 0 ),
	num_substitutions_since_surface_min_update_( 0 ),
	observed_sufficient_surface_E_to_predict_min_( false ),
	surface_exposed_( false ),
	is_below_buried_residue_no_hsasa_cutoff_( false ),
	num_neighbors_counting_self_(-1)
{
	// the weight to apply to the surface score is stored in the packer task. this allows users to change the
	// weight with a command line option
	surface_energy_weight_ = get_surface_owner()->surface_score_weight();
#ifdef FILE_DEBUG
	TR_NODE << "Setting surface_energy_weight to " << surface_energy_weight_ << std::endl;
#endif
}


///
/// @begin SurfaceNode< V, E, G >::~SurfaceNode
///
/// @brief
/// destructor -- no dynamically allocated data, does nothing
///
template < typename V, typename E, typename G >
SurfaceNode< V, E, G >::~SurfaceNode() {}


///
/// @begin SurfaceNode::prepare_for_simulated_annealing
///
/// @brief
/// invokes V's prep_for_simA method
/// Also populates a map defining which nodes are neighbors according to the tenA neighbor graph.  This map will be
/// checked whenever iterating over all edges defined by the energy graph to see if a Node really needs to be updated
/// for the surface score.
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::prepare_for_simulated_annealing() {

	// parent is AddtlBGNodesIG; V is either PDIG or LinmemIG; here we want to call the templated classes prep_for_simA
	// method so that it can get ready for all the pairwise-decomposable energy term stuff
	V::prepare_for_simulated_annealing();

	if ( ! parent::get_bg_edge_vector_up_to_date_() ) {
		parent::update_bg_edge_vector();
	}

	have_prepared_for_simA_ = true;

	// populate the neighbor map for fast lookups of "are they neighbors" later
	// Important: this map is also used in the calculate_amount_total_hydrophobic_ASA method, so make sure to init the map before
	// calling this method.
	int this_nodes_index = parent::get_node_index();

	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		int other_nodes_index = get_incident_surface_edge(ii)->get_other_ind( this_nodes_index );
		// set i,j and j,i in case when we do the lookup one key pair is not set
		if ( this->detect_neighborship_with_node( other_nodes_index, true /* this node is a FC node */ ) ) {
			fc_neighbor_map[ std::pair<int,int>(this_nodes_index, other_nodes_index) ] = 1;
			fc_neighbor_map[ std::pair<int,int>(other_nodes_index, this_nodes_index) ] = 1;
		}
	}

	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		int bg_nodes_index = get_edge_to_surface_bg_node( ii )->get_other_ind( this );
		if ( this->detect_neighborship_with_node( bg_nodes_index, false /* not a FC node */ ) ) {
			bg_neighbor_map[ std::pair<int,int>(this_nodes_index, bg_nodes_index) ] = 1;
			bg_neighbor_map[ std::pair<int,int>(bg_nodes_index, this_nodes_index) ] = 1;
		}
	}

	init_hASA_variables();

}


///
/// @begin SurfaceNode::detect_neighborship_with_node
///
/// @brief
/// Determines if this vertex neighbors the calling node
/// Called by all other nodes, and if true is returned an Edge of some sort is created between the nodes.
///
/// @detailed
/// To determine whether the passed in fc node is neighbors with this one, we need to iterate through the
/// list of neighbors in the tenA neighbor graph (a ContextGraph) that's stored in the Pose Energies object. Assumes
/// that the Pose has already been scored by some ScoreFunction.
/// The tenA neighbor graph is cool if we want to look at an interaction distance of 10A. But what if we want to only
/// find patches of exposed hydrophobics within an 8A sphere.  Then we can't use the tenA neighbor.  Instead use an
/// expensive distance measurement between the centroids of the wild-type residues to determine "are they neighbors".
///
/// The only problem with making this method generic for fc nodes or bg nodes is that when an index is passed in, it's
/// not possible to really tell if it's a fc node index or bg node index.  Rather than making two separate methods for
/// detecting neighborship, add a second parameter to this function which specificies if the node being checked is fc
/// or not.  That way we can get the pose residue object correct.
///
/// This method basically makes the decision for whether to create a SurfaceEdge or SurfaceBackgroundEdge between
/// two Nodes.  It really is what defines how much of the surface we design for.
///
template < typename V, typename E, typename G >
bool SurfaceNode< V, E, G >::detect_neighborship_with_node( int node_index, bool first_class ) const {

	// owner is the graph which has a reference to the Pose
	pose::Pose poseRef = get_surface_owner()->pose();
	// core::scoring::TenANeighborGraph const & tenA_neighbor_graph( poseRef.energies().tenA_neighbor_graph() );

	int fc_node_index = parent::get_node_index();

	// since this method is generic for detecting neighborship with either bg nodes or fc nodes, we need to figure out
	// whether the passed in node_index is for a background node or fc node.  Only way I can think of to do that here
	// is to check the bg node array in the SIG and see if it returns a non-negative value. If it return -1, then that
	// node_index was not found and it must be a fc node. Right?
	// Changing this method to have a second parameter. Use the second parameter boolean to tell what type of node it is
	int other_node_resid;
	if ( first_class ) {
		other_node_resid = get_surface_owner()->rotamer_sets().moltenres_2_resid( node_index );
	} else {
		other_node_resid = get_surface_owner()->bg_node_2_resid( node_index );
	}

	// we need to determine whether the centroid of the side chain (is that the same as the action coordinate?)
	// of this residue has a distance less than the centroid of the res1 sidechain.  This method is a brute
	// force way of determining "neighbor".
	conformation::Residue const & rsd1 = poseRef.residue( get_surface_owner()->rotamer_sets().moltenres_2_resid( fc_node_index ) );
	conformation::Residue const & rsd2 = poseRef.residue( other_node_resid );

	Real distanceBetweenAtoms =  rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) );

#ifdef FILE_DEBUG
	// toooo much output
	//TR_NODE << "detect_neighborship_with_node: fc node " << fc_node_index << " (resid: " <<
	//	get_surface_owner()->rotamer_sets().moltenres_2_resid( fc_node_index ) << ") checking if ";
	//if ( first_class ) { TR_NODE << "fc node "; } else { TR_NODE << "bg node "; }
	//TR_NODE	<< node_index << " (resid: " << other_node_resid << ") is a neighbor; distance between atoms " << rsd1.name3() << " " << rsd1.seqpos()
	//	<< " " << rsd1.atom_name( rsd1.nbr_atom() ) << " - " << rsd2.atom_name( rsd2.nbr_atom() )
	//	<< rsd2.name3() << " " << rsd2.seqpos() << ": " << distanceBetweenAtoms << std::endl;
#endif

	if ( distanceBetweenAtoms <= INTERACTION_RADIUS ) {
		return true;
	}

	// for every Edge in the neighbor graph, figure out if that residue is surface exposed *and* hydrophobic
	/* for ( core::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node( fc_node_index )->const_edge_list_begin(),
		eli_end = tenA_neighbor_graph.get_node( fc_node_index )->const_edge_list_end(); eli != eli_end; ++eli ) {
		if ( (*eli)->get_other_ind( fc_node_index ) == node_index ) {
			return true;
		}
	} */

	// if not within the interaction radius, not neighbors
	return false;

}



///
/// @begin SurfaceNode::assign_zero_state
///
/// @brief
/// Assign the node to state 0 -- the "unassigned" state.
///
/// @detailed
/// A node in state 0 produces no surface score.  Its neighbors have to adjust their scores appropriately. This method
/// iterates over all the edges emanating from this node and tells them to acknowledge that they've been zeroed out.
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::assign_zero_state() {

#ifdef FILE_DEBUG
	//TR_NODE << "assign_zero_state - node " << parent::get_node_index() << " going to state 0." << std::endl;
#endif

	// depending on the type of interaction graph (linmem, standard pd) used for the simulation, this will call
	// the appropriate assign_zero_state() method and initialize the one- and two-body energy term variables to zero
	// parent refers to the AddtlBGNodesIG; G refers to either PDIG or LinmemIG.  AddtlBGNodesIG doesn't define
	// a assign_zero_state() method! But it extends from the other two classes so it will eventually find the right method.
	parent::assign_zero_state();  // was this a bug previously?  apparently not...
	//V::assign_zero_state();

	// get_incident_surface edge() is an inlined protected method which just calls the parents get_incident_edge() method
	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		get_incident_surface_edge( ii )->acknowledge_state_zeroed_surface( parent::get_node_index() );
	}

	parent::update_bg_edge_vector();

	for (int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		get_edge_to_surface_bg_node( ii )->acknowledge_state_change( 0 );
	}

}


///
/// @begin SurfaceNode::acknowledge_neighbors_substitution_surface
///
/// @brief
/// bookkeeping to follow a neighbors state substitution.  this method gets called when a SurfaceNode commits a sub
/// and then broadcasts that change to all its neighboring fc nodes via the incident SurfaceEdges. basically we need
/// to set current state equal to alt state here. (Hopefully alt state is still correct!!)  Since there's no way for
/// a SurfaceNode to know what other SurfaceNodes are connected to it except via SurfaceEdges, the calls seem a bit
/// complicated.  A SurfaceNode has to call acknowledge_state on each Edge.  The Edges have to figure out which Node
/// is changing/not changing and then they call an inform_non_changing node of change method.  That method then makes
/// the call to this method on the correct SurfaceNode. The inform_non_changing method can not be removed, because
/// it's used during the substitution evaluations as well.
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::acknowledge_neighbors_substitution_surface() {

#ifdef FILE_DEBUG
	// too much output...
	//if ( std::fabs( curr_state_total_hASA_ - alt_state_total_hASA_ ) > 0.01 ) {
	//	TR_NODE << "acknowledge_neighbors_substitution_surface - node " << parent::get_node_index() << " changing curr_state_total_hASA_ from "
	//		<< curr_state_total_hASA_ << " to " << alt_state_total_hASA_ << std::endl;
	//}
#ifdef DOUBLE_CHECK_COUNTS
	Real previous_state_hASA = curr_state_total_hASA_;
#endif
#endif

	curr_state_total_hASA_ = alt_state_total_hASA_;

#ifdef DOUBLE_CHECK_COUNTS
	if ( std::fabs( curr_state_total_hASA_ - previous_state_hASA ) > 0.01 ) {
		// only double-check the counts for nodes that will be assigned a surface score, i.e. surface_exposed bool set to true
		if ( parent::get_current_state() != 0 && is_surface_exposed() ) {
			verify_patch_areas_correct( parent::get_node_index(), parent::get_current_state(), previous_state_hASA );
		}
	}
#endif

}


///
/// @begin SurfaceNode::assign_state_surface
///
/// @brief
/// Assigns the node to one of the states in its state space or to the unassigned state.
///
/// @detailed
/// Piggy backs on the standard project_deltaE/commit_considered_substitution protocol. Slightly inefficient; however,
/// this code is executed once per simulated annealing run.  Not worth optimizing. (apl)
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::assign_state_surface( int state ) {

#ifdef FILE_DEBUG
	TR_NODE << "assign_state_surface(): node " << parent::get_node_index() << " being assigned state " << state << std::endl;
#endif

	if ( state == 0 ) {
		assign_zero_state();
		return;
	}

	// the inefficiency lies in the fact that every time assign_state_surface is called with a nonzero value,
	// the above conditional is evaluated before project deltaE is called.
	float temp(0.0f);
	project_deltaE_for_substitution_surface( state, temp, -1.0f /* threshold for avoiding calculations; basically, don't procrastinate */ );
	commit_considered_substitution_surface();
}


///
/// @begin SurfaceNode::get_curr_state_surface_energy()
///
/// @brief
/// returns the surface energy for the node in its current state assignment
///
template < typename V, typename E, typename G >
Real SurfaceNode< V, E, G >::get_curr_state_surface_energy() const {

	if ( !(this->is_surface_exposed()) ) {
		return 0.0;
	}

	// the surface energy should be zero in the unassigned state?  yup!
	if ( parent::get_current_state() == 0 ) {
		return 0.0;
	}

	if ( curr_state_total_hASA_ > MAX_PATCH_SURFACE_AREA ) {
		return surface_energy_weight_ * MAX_SURFACE_ENERGY;
	}

	return surface_energy_weight_ * hASA_energy( curr_state_total_hASA_ );
}


///
/// @begin SurfaceNode::project_deltaE_for_substitution_surface
///
/// @brief
/// Returns the (possibly approximate) change in energy induced by changing a node from its current state into some alternate state.
///
/// @detailed
/// Iterates across the edges first to determine the sum of the pairwise
/// decomposable (PD) energy.  If the PD difference implies a collision, then
/// the SurfaceNode pretends as if the state substitution causes the best
/// improvement possible in surface score for it and its neighbors, returns
/// the PD difference + pretend surface difference.  It will procrastinate
/// computing the actual surface score difference until the guiding SimAnnealer
/// decides to commit the substitution.  If the SimAnnealer rejects the
/// substitution, then the work to compute the surface score is never done.
/// If it is unclear that the SimAnnealer will reject the substitution
/// based on the PD difference, then the vertex computes the change in
/// surface score that it induces on its neighbors.
/// The surface score calculation really isn't that intensive but it doesn't
/// make sense to calculate it needlessly.  (apl)
///
/// Bug source: Incident edges to a node are edges which have some energy on them based on the score function used (meaning they connect
/// to nodes that interact with this node).  This set of nodes is larger than the set of nodes that the tenA neighbor graph in the Pose
/// considers neighbors for this node.  Thus, when calculating the true number of se hp neighbors, some nodes are updated when
/// they shouldn't be leading to incorrect counts.  Possible solutions include 1) keeping some type of data structure (e.g. a map)
/// which specifies which residues are neighbors (in the tenA neighbor graph sort of way) and whenever iterating over incident edges
/// checking to make sure nodes are really neighbors before doing any updating of the counts; 2) forget about which nodes are truly
/// neighbors (as determined by the tenA neighbor graph) and just use whatever the energy function thinks neighbors should be. I don't
/// like this approach though, because I want to have the freedom to adjust what kind of radius I use for residue neighbors for future
/// versions of the surface energy term.  Not just assume that the energy function picks neighors correctly.
///
/// Where should the neighbor cache be kept?  If I keep it on the Node, then the lookups would have to occur in the Node methods that
/// change the counts.  If I keep it on the Edges though, I can check for neighborship when iterating over edges and not ever even
/// go into the Node methods. Keep it on the Nodes.  That's where it makes sense to keep it.
///
///
/// @param
/// alternate_state - [in] - the alternate state that the node must consider
/// prev_PDenergies_for_node - [out] - the sum of the PD energies (one-and two-body) for the node in its current state
/// deltaE_thresh_for_avoiding_surface_calcs - [in] - the minimum deltaE for which the SimAnnealer will accept slopyness in the deltaE projections
///
template < typename V, typename E, typename G >
Real SurfaceNode< V, E, G >::project_deltaE_for_substitution_surface(
	int alternate_state, core::PackerEnergy & prev_PDenergies_for_node, float deltaE_thresh_for_avoiding_surface_calcs ) {

	runtime_assert( alternate_state > 0  && alternate_state <= parent::get_num_states() );

	++num_state_substitutions_considered_;
	calculated_surface_deltaE_ = false;

	prev_PDenergies_for_node = parent::get_curr_pd_energy_total();

	// have base class perform deltaE computations for PD portion of energy function
	// parent refers to AddtlBGNodesIG; G refers to either PDIG or LinmemIG.  Which do I want here? AddtlBGNodesIG doesn't
	// implement the calc_deltaEpd method but it extends either PDIG or LinmemIG so the runtime will eventually find the
	// right function to call.
	parent::calc_deltaEpd( alternate_state );

	deltaE_for_substitution_ = parent::get_alt_pd_energy_total() - parent::get_curr_pd_energy_total();

	if ( decide_procrastinate_surface_computations( deltaE_for_substitution_, deltaE_thresh_for_avoiding_surface_calcs ) ) {
		++num_surface_comps_procrastinated_;
	} else {
#ifdef FILE_DEBUG
		TR_NODE << "deltaE for PD: " << deltaE_for_substitution_ << std::endl;
#endif
		deltaE_for_substitution_ += project_surface_deltaE();
	}

	return deltaE_for_substitution_;
}


///
/// @begin SurfaceNode< V, E, G >::project_surface_deltaE()
///
/// @brief
/// returns the change in surface score for this node and all of its neighbors produced by switching from its current
/// state to an alternate state.
///
template < typename V, typename E, typename G >
Real SurfaceNode< V, E, G >::project_surface_deltaE() {

#ifdef FILE_DEBUG
	TR_NODE << "project_surface_deltaE(): current_state:" << parent::get_current_state();
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


	// Next is where we update the total hASA of *all* se nbs by incrementing or decrementing the cached amount (instead of
	// reiterating over all edges and recalculating).  The cached amount is just a float of the amount of hASA total this
	// node has.  We should only change that amount if *this* node is SE and changing states. One bug I previously had
	// was that the SE check was not being done and amounts were getting incremented too high!

	// only change the total hASA (at this node) if this residue is surface exposed! but don't quit here because the change
	// at this node could influence the total hASA at a neighboring node!
	if ( this->is_surface_exposed() ) {

		// Previously, hASA calculations were only done for certain state changes. But in this version of the SIG, we're
		// going to consider all neighbors NOT just the hydrophobic ones. So every substitution that causes a state change
		// will lead to a change in the hASA

		if ( parent::get_current_state() == 0 ) { // runs only once per sim annealing run

			chemical::AA curr_AA = wt_residue_for_node().aa();
			chemical::AA alt_AA = get_rotamer( parent::get_alternate_state() )->aa();
			Size nbs = num_neighbors_counting_self();

			alt_state_total_hASA_ = curr_state_total_hASA_ - average_residue_hASA( curr_AA, nbs ) + average_residue_hASA( alt_AA, nbs );

			if ( alt_state_total_hASA_ < -0.001 ) {
				this->print();
				utility_exit_with_message( "Nonsensical alt_state reached (alt state total hASA < 0). Exiting. " );
			}

		} else { // node wasn't in the unassigned state - 99% of calls will follow this branch

			chemical::AA curr_AA = get_rotamer( parent::get_current_state() )->aa();
			chemical::AA alt_AA = get_rotamer( parent::get_alternate_state() )->aa();
			Size nbs = num_neighbors_counting_self();

			alt_state_total_hASA_ = curr_state_total_hASA_ - average_residue_hASA( curr_AA, nbs ) + average_residue_hASA( alt_AA, nbs );

		}
	}


	// iterate over all the edges coming from this node and get the deltaE for surface from the neighboring nodes
	Real surface_deltaE = 0;

	// the incident edges of the "parent" are the edges of the PDNode or LinearMemNode
#ifdef FILE_DEBUG
	TR_NODE << "project_surface_deltaE(): calculating deltaE for " << parent::get_num_incident_edges() << " first class and "
			<< parent::get_num_edges_to_background_nodes() << " background neighbors." << std::endl;
#endif
	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		int other_nodes_index = get_incident_surface_edge(ii)->get_other_ind( parent::get_node_index() );
		// "search" through the map to make sure this edge is in the map (if not found, find returns map.end() which is what I check for)
		if ( fc_neighbor_map.find( std::pair<int,int>( parent::get_node_index(), other_nodes_index) ) != fc_neighbor_map.end() ) {
			surface_deltaE += get_incident_surface_edge(ii)->get_surface_deltaE_for_neighbor( parent::get_node_index(), parent::get_alternate_state() );
		}
	}

	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		int other_nodes_index = get_edge_to_surface_bg_node(ii)->get_other_ind( this );
		// "search" through the map to make sure this edge is in the map (if not found, find returns map.end() which is what I check for)
		if ( bg_neighbor_map.find( std::pair<int,int>( parent::get_node_index(), other_nodes_index) ) != bg_neighbor_map.end() ) {
			surface_deltaE += get_edge_to_surface_bg_node( ii )->get_surface_deltaE_for_substitution( parent::get_alternate_state() );
		}
	}

	calculated_surface_deltaE_ = true;

	surface_deltaE += get_surface_score_difference();

#ifdef FILE_DEBUG
	if ( get_surface_score_difference() != 0 ) {
		TR_NODE << "project_surface_deltaE(): curr state hASA: " << curr_state_total_hASA_ << ", alt state hASA: " << alt_state_total_hASA_
			<< ", score difference: " << get_surface_score_difference() << std::endl;
	}
#endif

#ifdef FILE_DEBUG
	if ( surface_deltaE > 100000 || surface_deltaE < -100000 ) {
		std::cerr << "project_surface_deltaE(): deltaE for PD: " << deltaE_for_substitution_ << std::endl;
		std::cerr << "project_surface_deltaE(): deltaE for surfaceE: " << surface_deltaE << std::endl;

		std::cerr << "project_surface_deltaE(): current_state:" << parent::get_current_state();
		if ( parent::get_current_state() == 0 ) {
			std::cerr << " (" << wt_residue_for_node().name3() << "-";
			if ( wt_residue_for_node().is_polar() ) { std::cerr << "P)"; } else { std::cerr << "HP)"; }
		} else {
			std::cerr << " (" << get_rotamer( parent::get_current_state() )->name() << "-";
			if ( get_rotamer( parent::get_current_state() )->is_polar() ) { std::cerr << "P)"; } else { std::cerr << "HP)"; }
		}
		std::cerr << ", alternate_state: " << parent::get_alternate_state() << " (" << get_rotamer( parent::get_alternate_state() )->name() << "-";
		if ( get_rotamer( parent::get_alternate_state() )->is_polar() ) { std::cerr << "P)"; } else { std::cerr << "HP)"; }
		std::cerr << std::endl;

		get_surface_owner()->print();
		utility_exit_with_message( "surface deltaE went crazy. Terminating." );
	}
#endif

#ifdef FILE_DEBUG
	if ( surface_deltaE != 0 ) {
		TR_NODE << "deltaE for surfaceE: " << surface_deltaE << std::endl;
	}
#endif

	return surface_deltaE;

}


///
/// @begin SurfaceNode< V, E, G >::decide_procrastinate_surface_computations()
///
/// @detailed
/// Makes the decision whether or not to procrastinate calculating the surface score. Basically, if the PD energy got better (dE < 0)
/// then return false so we don't procrastinate the calculation (because the alternate state probably will be accepted?). If the best
/// guess for the surface deltaE also comes back better (dE < 0), then return false.  Finally, if the difference between the deltaE
/// for the PD terms and the (guessed) surface deltaE is greater than the threshold, return true so we do procrastinate. So basically
/// if the energy (especially the PD energy) gets worse, procrastinate.  Otherwise, don't.
///
template < typename V, typename E, typename G >
bool SurfaceNode< V, E, G >::decide_procrastinate_surface_computations( Real const pd_deltaE, Real const threshold ) const {

	Real surface_deltaE_max = 0;

	if ( ! observed_sufficient_surface_E_to_predict_min_ )
		return false;

	if ( threshold < 0 || pd_deltaE < 0 )
		return false;

	for( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii) {
		Real mag_deltaE = get_incident_surface_edge(ii)->get_max_surface_deltaE_guess( parent::get_node_index() );
		if ( mag_deltaE < 0.0f ) {
			return false;
		}
		surface_deltaE_max += mag_deltaE;
	}

	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		Real mag_deltaE = get_edge_to_surface_bg_node( ii )->get_max_surface_deltaE_guess();
		if ( mag_deltaE < 0.0f ) {
			return false;
		}
		surface_deltaE_max += mag_deltaE;
	}

	surface_deltaE_max += (surface_energy_weight_ * hASA_energy( curr_state_total_hASA_ ) ) - surface_score_min_last_100_;

	// pd_deltaE must be positive, surface deltaE must also be positive.
	// threshold of 5 means, if the PD got more than 5 worse than the guessed surface deltaE (e.g. 10 - 3(?) = 7 > 5), procrastinate
#ifdef FILE_DEBUG
	TR_NODE << "decide_procrastinate_surface_computations(): pd_deltaE: " << pd_deltaE << ", surface_deltaE_max: " << surface_deltaE_max
			<< ", threshold: " << threshold << std::endl;
#endif
	if ( (pd_deltaE - surface_deltaE_max) > threshold ) {
		return true;
	}
	return false;

}



///
/// @begin SurfaceNode< V, E, G >::get_surface_deltaE_for_neighbors_state_substitution
///
/// @brief
/// returns the change in surface score for this node induced by a state substitution at a neighboring node.
///
/// @detailed
/// So a neighboring first class node that is changing is invoking (via the SurfaceEdge that connects them) this method
/// so that this node can update its counts and then return a change in surface score for this node. This method needs
/// to check if the alt_state that the connected node is considering is going to change the counts here and then perform
/// the lookup for the change in score.
///
/// is the rotamer id (state) the same for every node? no. need to either do the work in the edge class method
/// or pass a reference to the node.  oooooh, we can use the SIG to get the node changing, too!!!  NO!  Can't. The
/// SIG get_surface_node( index ) method is protected.  Lame.
/// Since I need a way to check the rotamer vector for the changing node (to see what the changing states are) I have to
/// make this method take a raw pointer to the node changing.  The alternative would be to do all the state checking in the
/// edge class method, but I feel like this sort of state checking should be done by a Node, not an Edge.  Further, this
/// frees the Edge class from having to keep any of the data that the Nodes contains (i.e. the counts of neighbors).
/// The edges just pass along the message (from the changing node) to the other Node, not doing any kind of Node-specific
/// behaviour.
///
/// the current state of this node may be unassigned (0). if so, we need to initialize the hASA properly.
/// as of 7/22, the current amount of hydrophobic surface area should be init'd after prep for simA runs. alt state also.
/// if ( parent::get_current_state() == 0 )
///	curr_state_total_hASA_ = get_surface_owner()->calculate_amount_total_hydrophobic_ASA;
///
/// 07/24/08 problem situation
/// it's probable that a sub that was considered before wasn't committed, so the alt state count at this node needs to go
/// back to what the current state count is.  the problem that may crop up with the reset here is that on a previous consider
/// call, a set of nodes will update their counts.  if commit is not called, and a second consider is called, then this node
/// will have its alt state count reset but what about all the other nodes in the previous set.  how will their counts get
/// reset to the current state count?  perhaps a check in the commit method can be added.  nope, a new method has been
/// added to Nodes and BGNodes to handle this case.
///	alt_state_total_hASA_ = curr_state_total_hASA_;
///
/// 02/20/09 Changed all variable names above to reflect hASA instead of counts.
///
template < typename V, typename E, typename G >
Real SurfaceNode< V, E, G >::get_surface_deltaE_for_neighbors_state_substitution(
			SurfaceNode<V,E,G>* node_considering_substitution, int changing_nodes_curr_state, int changing_nodes_alt_state ) {

	// If this node is not surface exposed, then it doesn't have a surface score. Even if another neighboring node is
	// changing state, this node should be unaffected because it's not on the surface.  Only residues which are on the
	// surface have a surface score associated with them.  So don't bother updating the hASA if we're just going to return
	// 0.0 anyway.  Can't think of any reason for updating the hASA, thus I'm moving this check up here.
	// Is this the correct usage of is_se()? Yes. This node won't have a surface deltaE because it has fewer nbs than the
	// cutoff we're using. But this node *could* still add hASA to other nodes that will get surface scores.
	if ( !(this->is_surface_exposed()) ) {
		return 0.0;
	}

	// only change the count(hASA) here if the node that is changing is surface exposed!!!!  otherwise, leave it at what it was
	// initialized in the hASA function to be.  this nodes hASA only needs to change if the changing node is both surface
	// exposed and is changing state.
	// WRONG!  The neighboring Node might not be surface exposed, but if it's below the buried residue cutoff it could still
	// contribute some hASA to this Node.
	//if ( node_considering_substitution->is_surface_exposed() ) {
	if ( node_considering_substitution->is_below_buried_residue_no_hsasa_cutoff() ) {

		// it may also be that the node undergoing substitution is in the unassigned state, right?
		// if that's the case, we need to look up the wild-type residue info for that node (NOT FOR THIS NODE!)
		if ( changing_nodes_curr_state == 0 ) {

			chemical::AA curr_AA = node_considering_substitution->wt_residue_for_node().aa();
			chemical::AA alt_AA = node_considering_substitution->get_rotamer( changing_nodes_alt_state )->aa();
			Size nbs = node_considering_substitution->num_neighbors_counting_self();

			alt_state_total_hASA_ = curr_state_total_hASA_ - average_residue_hASA( curr_AA, nbs ) + average_residue_hASA( alt_AA, nbs );

			if ( alt_state_total_hASA_ < -0.001 ) {
				this->print();
				utility_exit_with_message( "Nonsensical alt_state reached (total hASA < 0). Exiting. " );
			}

		} else { // node wasn't in the unassigned state - 99% of calls will follow this branch

			chemical::AA curr_AA = node_considering_substitution->get_rotamer( changing_nodes_curr_state )->aa();
			chemical::AA alt_AA = node_considering_substitution->get_rotamer( changing_nodes_alt_state )->aa();
			Size nbs = node_considering_substitution->num_neighbors_counting_self();

			alt_state_total_hASA_ = curr_state_total_hASA_ - average_residue_hASA( curr_AA, nbs ) + average_residue_hASA( alt_AA, nbs );
		}
	}

#ifdef FILE_DEBUG
	if ( curr_state_total_hASA_ != alt_state_total_hASA_ ) {
		TR_NODE << "get_surface_deltaE_for_neighbors_state_substitution: node " << parent::get_node_index()
			<< " calc. deltaE for changing node " << node_considering_substitution->get_node_index()
			<< ", curr hASA: " << curr_state_total_hASA_ << ", alt hASA: " << alt_state_total_hASA_
				<< ", score_difference: " << get_surface_score_difference() << std::endl;
	}
#endif

	return get_surface_score_difference();
}


///
/// @begin SurfaceNode::get_surface_score_difference
///
/// @brief
/// Returns the difference alt state score - current state score
///
/// @detailed
/// This method requires that the variables current_state_ and alternate_state_ hold the correct values. An assert guard
/// checks to make sure an out-of-bounds exception won't occur.
///
///
template < typename V, typename E, typename G >
Real SurfaceNode< V, E, G >::get_surface_score_difference() const {

	if ( curr_state_total_hASA_ == alt_state_total_hASA_ ) { return 0.0; }

	if ( curr_state_total_hASA_ < -0.001 ) {
		this->print();
		utility_exit_with_message( "Nonsensical curr_state reached (total hASA < 0). Exiting. " );
	}
	if ( alt_state_total_hASA_ < -0.001 ) {
		this->print();
		utility_exit_with_message( "Nonsensical alt_state reached (total hASA < 0). Exiting. " );
	}

	if ( ( curr_state_total_hASA_ > MAX_PATCH_SURFACE_AREA ) && ( alt_state_total_hASA_ > MAX_PATCH_SURFACE_AREA ) ) {
		// this sub could be good or bad, but let the contributions of the other nodes decide what happens
		return 0.0;
	}

	// if the alt state is greater than 800A^2, return a really bad energy!
	if ( alt_state_total_hASA_ > MAX_PATCH_SURFACE_AREA ) {
		return surface_energy_weight_ * MAX_SURFACE_ENERGY;
	}
	// alt state must be 800 or less, which is in bounds; moving on

	if ( curr_state_total_hASA_ > MAX_PATCH_SURFACE_AREA ) {
		// then, the current state has more hASA than the alt state, which is a favorable change
		return surface_energy_weight_ * -1 * MAX_SURFACE_ENERGY;  // yeah, it's a *really* favorable change, but oh well
	}

	return surface_energy_weight_ * ( hASA_energy(alt_state_total_hASA_) - hASA_energy(curr_state_total_hASA_) );

}


///
/// @begin SurfaceNode::commit_considered_substitution_surface
///
/// @brief
/// Sets the current state to the alternate state this node was asked to
/// consider.  Copies appropriate score information.  Notifies all of its
/// neighbors that it is going through with the state substitution it had been
/// considering.
///
/// @detailed
/// If the node had procrastinated the surface deltaE computation because it
/// thought the sim-annealer unlikely to accept the substitition, then the
/// node must now do the work it procrastinated.
///
/// There's a potential situation with considers() and commits() that needs to be checked for here. It's possible that
/// a consider() call is made which causes a set of Nodes to update their alt state counts correspondingly.  Since the
/// consider() only gets processed by (or actually the call only goes out to) Nodes which are neighbor graph neighbors,
/// only those Nodes will have their counts updated. Assume that the first consider() is really bad and no commit goes out.
/// If we then consider() another sub, the alt state counts at the previous consider()'s set of nodes are incorrect.  If
/// get_deltaE is called by this consider() on some of those nodes, some of them will have their alt states reset.  But
/// when the commit goes out to ALL Nodes that are neighbors (NOT just the neighbor graph neighbors) it's possible that
/// some of the Nodes will save the wrong alt state count.  One way I think this can be avoided to is check at this node,
/// if the alt state count is different from current, whether the node that originally changed is a neighbor graph neighbor
/// of this node.  If it is, that means the counts changed because of that node.  If it's not, then this Node must have
/// been one of the ones that fell out of sync.
/// Oooooh, I just thought of another way.  When the SIG consider() method is called, I could have a reset alt state
/// counts method that will deal with Nodes that are out of sync.  That's more elegant than yet another if statement here!
///
template < typename V, typename E, typename G >
Real SurfaceNode< V, E, G >::commit_considered_substitution_surface() {

	assert( parent::considering_alternate_state() );

#ifdef DOUBLE_CHECK_COUNTS
	int previous_state = parent::get_current_state();
	Real previous_state_hASA = curr_state_total_hASA_;
#endif

	if ( ! calculated_surface_deltaE_ ) {
		++num_surface_comps_later_made_;
		Real temp = project_surface_deltaE();
		deltaE_for_substitution_ = parent::get_alt_pd_energy_total() - parent::get_curr_pd_energy_total() + temp;
	}

	// call base class method
	parent::commit_considered_substitution();

	curr_state_total_hASA_ = alt_state_total_hASA_;


#ifdef DOUBLE_CHECK_COUNTS
	if ( ( parent::get_current_state() != 0 ) && is_surface_exposed() ) {
		verify_patch_areas_correct( parent::get_node_index(), previous_state, previous_state_hASA );
	}
#endif

#ifdef FILE_DEBUG
	TR_NODE << "Committed substitution node " << parent::get_node_index() << std::endl;
#endif

	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		get_incident_surface_edge(ii)->acknowledge_substitution_surface();
	}
	for (int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii) {
		get_edge_to_surface_bg_node(ii)->acknowledge_substitution_surface();
	}

	track_surface_E_min();

	return deltaE_for_substitution_;
}



///
/// @begin SurfaceNode< V, E, G >::verify_patch_areas_correct
///
/// @brief
/// Checks to make sure that the current state num se hp nbs is correct, by going through all residues in the Pose and
/// brute-force recounting.
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::verify_patch_areas_correct( int node_id, int previous_state, Real previous_state_hASA ) {

	// do a little double checking here to make sure things aren't getting off track somehow
	// keep all this wrapped within a #ifdef because it's very, very slow and inefficient
	// the only nodes that are going to have counts that are correct, though, will be the ones that are surface exposed;
	// nodes that aren't surface exposed don't update their counts but instead return immediately to save computational time

	Real total_hASA = 0.0;
	if ( this->is_surface_exposed() ) {
		//total_hASA += average_residue_hASA( wt_residue_for_node().aa(), num_neighbors_counting_self() );
		total_hASA += average_residue_hASA(); // not sure why this was wt_residue before
		TR_NODE << "verify_patch_areas_correct(): adding " << average_residue_hASA() << " to total_hASA for self; total_hASA: " << total_hASA << std::endl;
	}

	// have to iterate through the edges emanating from this node in the tenA neighbor graph stored in the Pose.
	// Using the edges that are connected to this SurfaceNode (via SurfaceEdges) is not correct because some
	// of the EnergyMethods produce energies between Nodes that are not considered neighbors in the tenA neighor graph
	// and hence persist as SurfaceEdges.  But just because these two Nodes have an Edge between them in the PDGraph
	// doesn't mean we should count them as neighbors from the Pose neighbor standpoint.
	//
	// the tenA neighbor graph gives us the index of the other node, but how do we get that Node as an object so we
	// can call is_se() on it?  do the really slow operation of iterating over all the SurfaceEdges until we find the
	// to the Node we're interested in.
	//
	// as of 7/25, changing the above to iterate over all residues that fall within the radius of interaction instead of
	// whatever the tenA neighbor graph has in it.  instead of performing the expensive distance calculation though, use
	// the distance neighbor map that was created earlier to determine which nodes are close by.  No!  Can't do that
	// because then we don't have two independent methods of checking the true count.  Iterate over all residues here
	// and redo the expensive calculation.  The map is used to speed things up during the actual simulation.
	//

	//core::scoring::TenANeighborGraph const & tenA_neighbor_graph( get_surface_owner()->pose().energies().tenA_neighbor_graph() );

	//for ( core::graph::EdgeListConstIterator eli = tenA_neighbor_graph.get_node( parent::get_node_index() )->const_edge_list_begin(),
	//		eli_end = tenA_neighbor_graph.get_node( parent::get_node_index() )->const_edge_list_end(); eli != eli_end; ++eli ) {

	pose::Pose poseRef = get_surface_owner()->pose();

	conformation::Residue const & rsd1 = poseRef.residue( get_surface_owner()->rotamer_sets().moltenres_2_resid( parent::get_node_index() ) );
	Real distanceBetweenAtoms = 0.0;

	for ( Size rr = 1; rr <= get_surface_owner()->rotamer_sets().nmoltenres(); ++rr ) {

		if ( (int)rr == parent::get_node_index() ) { continue; } // skip if we're on the current residue

		// moltenres_2_resid will go from the node/mr_id to the pose residue id
		conformation::Residue const & rsd2 = poseRef.residue( get_surface_owner()->rotamer_sets().moltenres_2_resid( rr ) );

		distanceBetweenAtoms = rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) );
		if ( distanceBetweenAtoms > INTERACTION_RADIUS ) {
			continue; // that residue is not a neighbor
		}
		// if less than the cutoff, then we need to figure out which Edge to signal. the index 'rr' is just the moltenres_id which
		// should be the same as the node id.

		// save the value to simplify code ahead
		// int other_node_index = (*eli)->get_other_ind( parent::get_node_index() );
		int other_node_index = rr;
		for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {

			if ( get_incident_surface_edge(ii)->get_other_ind( parent::get_node_index() ) == other_node_index ) {
				// got the right SurfaceEdge, now call the right methods on the Node object - but how do we get access
				// to the Node.  We need to determine which index (0 or 1) the changing node is.
				// An inefficient way to do it is to get the node with index 0 and call get_node_index on it, and check if that
				// node has the index that we have here.  If not, get the node with index 1.  There has to be a better way, but
				// this way will work.  This is just debugging code after all.
				int edge_index = ( other_node_index == get_incident_surface_edge(ii)->get_surface_node(0)->get_node_index() ? 0 : 1 );

				//if ( get_incident_surface_edge(ii)->get_surface_node( edge_index )->is_surface_exposed() ) {

					// don't bother checking the patch areas if there's still unassigned nodes
					if ( (get_incident_surface_edge(ii)->get_surface_node( edge_index ))->get_current_state() == 0 ) {
						return;
					}
					total_hASA += (get_incident_surface_edge(ii)->get_surface_node( edge_index ))->average_residue_hASA();
					TR_NODE << "verify_patch_areas_correct(): adding " << (get_incident_surface_edge(ii)->get_surface_node( edge_index ))->average_residue_hASA()
						<< " to total_hASA for incident edge " << ii << ", total_hASA: " << total_hASA << std::endl;
				//}
			}
		}
	}

	// now do all of this over again for the background nodes!  the total number of background nodes is just the total number
	// of residue - number of moltenresidues.  need to add a method which allows access to the bgenumeration_2_resid array
	// kept in the SIG.
	int nbackground = poseRef.total_residue() - get_surface_owner()->rotamer_sets().nmoltenres();
	for ( int id = 1; id <= nbackground; ++id ) {

		conformation::Residue const & rsd2 = poseRef.residue( get_surface_owner()->bg_node_2_resid( id ) );

		distanceBetweenAtoms =  rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) );
		if ( distanceBetweenAtoms > INTERACTION_RADIUS ) {
			continue; // that residue is not a neighbor
		}
		// if less than the cutoff, then we need to figure out which BGEdge to signal. we can iterate over all the BGEdges
		// and look for the one that has a node id the same as the index 'id'.

		int other_node_index = id;
		for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
			if ( get_edge_to_surface_bg_node(ii)->get_other_ind( this ) == other_node_index ) {

				// got the right SurfaceBackgroundEdge, now call the right methods on the BackgroundNode object
				//if ( get_edge_to_surface_bg_node(ii)->get_surface_bg_node()->is_surface_exposed() ) {
					total_hASA += (get_edge_to_surface_bg_node(ii)->get_surface_bg_node())->average_residue_hASA();
					TR_NODE << "verify_patch_areas_correct(): adding " << (get_edge_to_surface_bg_node(ii)->get_surface_bg_node())->average_residue_hASA()
						<< " to total_hASA for incident bg edge " << ii << ", total_hASA: " << total_hASA << std::endl;
				//}
			}
		}

	}

	if ( std::fabs( total_hASA - curr_state_total_hASA_ ) > 0.01 ) {

		TR_NODE << "verify_patch_areas_correct(): calling node " << node_id << "; this node current_state: " << parent::get_current_state();
		TR_NODE << " (" << get_rotamer( parent::get_current_state() )->name() << "-";
		if ( get_rotamer( parent::get_current_state() )->is_polar() ) { TR_NODE << "P)"; } else { TR_NODE << "HP)"; }

		TR_NODE << ", previous_state: " << previous_state;
		if ( previous_state != 0 ) {
			TR_NODE << " (" << get_rotamer( previous_state )->name() << "-";
			if ( get_rotamer( previous_state )->is_polar() ) { TR_NODE << "P)"; } else { TR_NODE << "HP)"; }
		}
		TR_NODE << ", previous_state hASA: " << previous_state_hASA << std::endl;

		std::cout << "This node:" << std::endl;
		this->print();

		std::cout << "Neighboring nodes:" << std::endl;
		for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
			int other_nodes_index = get_incident_surface_edge(ii)->get_other_ind( parent::get_node_index() );
			int edge_index = ( other_nodes_index == get_incident_surface_edge(ii)->get_surface_node(0)->get_node_index() ? 0 : 1 );
			// "search" through the map to make sure this edge is in the map (if not found, find returns map.end() which is what I check for)
			if ( fc_neighbor_map.find( std::pair<int,int>( parent::get_node_index(), other_nodes_index) ) != fc_neighbor_map.end() ) {

				if ( get_incident_surface_edge(ii)->get_surface_node( edge_index )->is_surface_exposed() ) {
					std::cout << "**";
				}
				get_incident_surface_edge(ii)->get_surface_node( edge_index )->print();
			}
		}

		for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
			int other_nodes_index = get_edge_to_surface_bg_node(ii)->get_other_ind( this );
			// "search" through the map to make sure this edge is in the map (if not found, find returns map.end() which is what I check for)
			if ( bg_neighbor_map.find( std::pair<int,int>( parent::get_node_index(), other_nodes_index) ) != bg_neighbor_map.end() ) {
				if ( get_edge_to_surface_bg_node(ii)->get_surface_bg_node()->is_surface_exposed() ) {
					std::cout << "**";
				}
				get_edge_to_surface_bg_node(ii)->get_surface_bg_node()->print();
			}
		}

		std::cout << "Nodes counted iterating over all residues:" << std::endl;

		{ // complete duplication of code above, but this is all debugging code!
				for ( Size rr = 1; rr <= get_surface_owner()->rotamer_sets().nmoltenres(); ++rr ) {

				if ( (int)rr == parent::get_node_index() ) { continue; } // skip if we're on the current residue

				// moltenres_2_resid will go from the node/mr_id to the pose residue id
				conformation::Residue const & rsd2 = poseRef.residue( get_surface_owner()->rotamer_sets().moltenres_2_resid( rr ) );

				distanceBetweenAtoms = rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) );
				if ( distanceBetweenAtoms > INTERACTION_RADIUS ) {
					continue;
				}
				// if less than the cutoff, then we need to figure out which Edge to signal. the index 'rr' is just the moltenres_id which
				// should be the same as the node id.

				// save the value to simplify code ahead
				// int other_node_index = (*eli)->get_other_ind( parent::get_node_index() );
				int other_node_index = rr;
				for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {

					if ( get_incident_surface_edge(ii)->get_other_ind( parent::get_node_index() ) == other_node_index ) {
						// got the right SurfaceEdge, now call the right methods on the Node object - but how do we get access
						// to the Node.  We need to determine which index (0 or 1) the changing node is.
						// An inefficient way to do it is to get the node with index 0 and call get_node_index on it, and check if that
						// node has the index that we have here.  If not, get the node with index 1.  There has to be a better way, but
						// this way will work.  This is just debugging code after all.
						int edge_index = ( other_node_index == get_incident_surface_edge(ii)->get_surface_node(0)->get_node_index() ? 0 : 1 );
						if ( get_incident_surface_edge(ii)->get_surface_node( edge_index )->is_surface_exposed() ) {
							std::cout << "**";
						}
						get_incident_surface_edge(ii)->get_surface_node( edge_index )->print();
					}
				}
			}

			// now do all of this over again for the background nodes!  the total number of background nodes is just the total number
			// of residue - number of moltenresidues.  need to add a method which allows access to the bgenumeration_2_resid array
			// kept in the SIG.
			int nbackground = poseRef.total_residue() - get_surface_owner()->rotamer_sets().nmoltenres();
			for ( int id = 1; id <= nbackground; ++id ) {

				conformation::Residue const & rsd2 = poseRef.residue( get_surface_owner()->bg_node_2_resid( id ) );

				distanceBetweenAtoms =  rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) );
				if ( distanceBetweenAtoms > INTERACTION_RADIUS ) {
					continue; // that residue is not a neighbor
				}
				// if less than the cutoff, then we need to figure out which BGEdge to signal. we can iterate over all the BGEdges
				// and look for the one that has a node id the same as the index 'id'.

				int other_node_index = id;
				for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
					if ( get_edge_to_surface_bg_node(ii)->get_other_ind( this ) == other_node_index ) {
						if ( get_edge_to_surface_bg_node(ii)->get_surface_bg_node()->is_surface_exposed() ) {
							std::cout << "**";
						}
						get_edge_to_surface_bg_node(ii)->get_surface_bg_node()->print();
					}
				}

			}
		}

		TR_NODE << "total_hASA: " << total_hASA << ", curr_state_total_hASA_: " << curr_state_total_hASA_ << std::endl;
		utility_exit_with_message( "total hASA values are out of sync. Something is wrong. Quitting." );
	}
}


///
/// @begin SurfaceNode< V, E, G >::track_surface_E_min
///
/// @brief
/// Keeps track of the minimum surface score seen.  Every 100 substitutions, updates the variable surface_score_min_last_100.
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::track_surface_E_min() {

	++num_substitutions_since_surface_min_update_;

	Real alt_surfaceE = surface_energy_weight_ * hASA_energy( alt_state_total_hASA_ );

	if ( surface_score_min_recent_ > alt_surfaceE )
		surface_score_min_recent_ = alt_surfaceE;

	// every 100 substitutions, update the variable surface_score_min_last_100 with the value held in
	// surface_score_min_recent_.  That explains the difference between last_100 and recent.
	if ( num_substitutions_since_surface_min_update_ == 100 ) {
		surface_score_min_last_100_ = surface_score_min_recent_;
		Real curr_surfaceE = surface_energy_weight_ * hASA_energy( curr_state_total_hASA_ );
		if (curr_surfaceE < surface_score_min_last_100_ ) {
			surface_score_min_last_100_ = curr_surfaceE;
		}
		observed_sufficient_surface_E_to_predict_min_ = true;
		num_substitutions_since_surface_min_update_ = 0;
	}

}


///
/// @begin SurfaceNode::print
///
/// @brief
/// useful for debugging - writes information about a node to the tracer
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::print() const {

	std::cout << "(node id:" << parent::get_node_index() << ", current_state: " << parent::get_current_state();

	if ( parent::get_current_state() == 0 ) {
		std::cout << " (" << wt_residue_for_node().name3() << "-";
		if ( wt_residue_for_node().is_polar() ) { std::cout << "P)"; } else { std::cout << "HP)"; }
	} else {
		std::cout << " (" << get_rotamer( parent::get_current_state() )->name() << "-";
		if ( get_rotamer( parent::get_current_state() )->is_polar() ) { std::cout << "P)"; } else { std::cout << "HP)"; }
	}

	if ( parent::get_current_state() == 0 ) {
		std::cout << ", hASA: " << average_residue_hASA( wt_residue_for_node().aa(), num_neighbors_counting_self() );
	} else {
		std::cout << ", hASA: " << average_residue_hASA();
	}

	std::cout << ", one body energy: " << parent::get_one_body_energy_current_state()
		<< ", curr: " << curr_state_total_hASA_ << ", alt: " << alt_state_total_hASA_ << ", surface energy: " << get_curr_state_surface_energy()
		<< ", se: " << surface_exposed_ << std::endl;

	//for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
	//	std::cout << "e:" << get_incident_surface_edge(ii)->get_first_node_ind() << "-" << get_incident_surface_edge(ii)->get_second_node_ind() << ", ";
	//}
	//std::cout << std::endl << std::endl;
}


///
/// @begin SurfaceNode::count_static_memory
///
/// @brief
/// Returns the amount of static memory used by this Node object
///
template < typename V, typename E, typename G >
unsigned int SurfaceNode< V, E, G >::count_static_memory() const {
	return sizeof( SurfaceNode< V, E, G > );
}

///
/// @begin SurfaceNode::count_dynamic_memory
///
/// @brief
/// Returns the amount of dynamic memory used by this Node object
///
template < typename V, typename E, typename G >
unsigned int SurfaceNode< V, E, G >::count_dynamic_memory() const {

	unsigned int total_memory = parent::count_dynamic_memory();

	// additional memory used?

	return total_memory;
}


///
/// @begin SurfaceNode::getMemoryUsageInBytes
///
/// @brief
/// Not implemented, but needs to be!
///
template < typename V, typename E, typename G >
unsigned int SurfaceNode< V, E, G >::getMemoryUsageInBytes() const {
	return 0;
}


///
/// @begin SurfaceNode::print_surface_avoidance_stats
///
/// @brief
/// reports on the level of success for surface score calculation procrastination
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::print_surface_avoidance_stats() {

	if (num_state_substitutions_considered_ == 0)
		return;

	TR_STATS << "SurfaceE Calculation Avoidance Statistics:" << std::endl;
	TR_STATS << "num state substitutions considered: " << num_state_substitutions_considered_ << ", "
			<< "num surface calcs procrastinated: " << num_surface_comps_procrastinated_ << ", "
			<< "num surface calcs later computed: " << num_surface_comps_later_made_ << std::endl;
	TR_STATS << "Percent Avoided: " << (double)(num_surface_comps_procrastinated_ - num_surface_comps_later_made_) / num_state_substitutions_considered_ << ", ";

	if ( num_surface_comps_procrastinated_ != 0 ) {
		TR_STATS << "Percent Procrastinated worthwhile: " <<
			(double)(num_surface_comps_procrastinated_ - num_surface_comps_later_made_) / num_surface_comps_procrastinated_ << std::endl;
	} else {
		TR_STATS << "Percent Procrastinated worthwhile: " << "N/A" << std::endl;
	}

}


///
/// @begin SurfaceNode::reset_surface_avoidance_stats
///
/// @brief
/// resets static member variables of SurfaceNode that measure how worthwhile
/// surface score calculation procrastination is.
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::reset_surface_avoidance_stats() {
	num_state_substitutions_considered_ = 0;
	num_surface_comps_procrastinated_ = 0;
	num_surface_comps_later_made_ = 0;
}


///
/// @begin SurfaceNode::set_rotamers
///
/// @details
/// Need to save a reference to the rotamer_set so that we can determine what a particular state change will do to the score
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::set_rotamers( rotamer_set::RotamerSetCOP rotamers ) {

	// get_num_states should call the parent graph's method?
	if ( rotamers->num_rotamers() != (Size) parent::get_num_states() ) {
		utility_exit_with_message( "Number of rotamers is not equal to parents number of states. Quitting.");
	}

#ifdef FILE_DEBUG
	TR_NODE << "set_rotamers: adding " << rotamers->num_rotamers() << " to local rotamers vector." << std::endl;
#endif

	for ( Size ii = 1; ii <= rotamers->num_rotamers(); ++ii ) {
		rotamers_vector_[ ii ] = rotamers->rotamer( ii );
	}

#ifdef FILE_DEBUG
	TR_NODE << "rotamers_vector_: [ ";
	for ( Size ii = 1; ii <= rotamers_vector_.size(); ++ii ) {
		TR_NODE << ii << ":" << rotamers_vector_[ ii ]->name1() << "-" << rotamers_vector_[ ii ]->seqpos() << ", ";
	}
	TR_NODE << "]" << std::endl;
#endif

}

///
/// @begin SurfaceNode::get_rotamer
///
/// @details
/// Need to save a reference to the rotamer_set so that we can determine what a particular state change will do to the score
///
template < typename V, typename E, typename G >
conformation::ResidueCOP SurfaceNode< V, E, G >::get_rotamer( int state ) const {
	return rotamers_vector_[ state ];
}




///
/// @begin SurfaceNode::is_surface_exposed
///
/// @brief
/// Returns the value of surface_exposed_ which gets set during the SIG initialize() method.
///
template < typename V, typename E, typename G >
bool SurfaceNode< V, E, G >::is_surface_exposed() const {
	return surface_exposed_;
}

///
/// @begin SurfaceNode::surface_exposed
///
/// @brief
/// setter for the surface_exposed_ bool
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::surface_exposed( bool value ) {
	surface_exposed_ = value;
}


///
/// @begin SurfaceNode::is_below_buried_residue_no_hsasa_cutoff
///
/// @brief
/// Returns the value of is_below_buried_residue_no_hsasa_cutoff_ which gets set during the SIG initialize() method.
///
template < typename V, typename E, typename G >
bool SurfaceNode< V, E, G >::is_below_buried_residue_no_hsasa_cutoff() const {
	return is_below_buried_residue_no_hsasa_cutoff_;
}

///
/// @begin SurfaceNode::is_below_buried_residue_no_hsasa_cutoff
///
/// @brief
/// setter for the is_below_buried_residue_no_hsasa_cutoff_ bool
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::is_below_buried_residue_no_hsasa_cutoff( bool value ) {
	is_below_buried_residue_no_hsasa_cutoff_ = value;
}


///
/// @begin SurfaceNode::reset_alt_state_total_hASA
///
/// @brief
/// Sets the alt state total hASA to the current state total hASA.  See comments in SIG and commit_considered_substitution_surface
/// for more information about why this method exists.
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::reset_alt_state_total_hASA() {

	// since we're only incrementing and decrementing the alt_state hASA, we need a way to reset it each iteration of consider()
	// so we don't get alt_state total hASA of nonsenical values which screw up the score.
	//
	if ( alt_state_total_hASA_ != curr_state_total_hASA_ ) {
		alt_state_total_hASA_ = curr_state_total_hASA_;
#ifdef FILE_DEBUG
		TR_NODE << "reset_alt_state_total_hASA(): alt state total hASA was not equal to current state total hASA, but is now." << std::endl;
#endif
	}

}


///
/// @begin SurfaceNode::initialize_num_neighbors_counting_self
///
/// @brief
/// Returns the number of neighbors within the interaction threshold, counting self.

/// @detailed
/// A member variable stores the value so that this method only has to run once and lookups can be performed thereafter.
/// The reason we need to know the number of neighbors is that determining whether a residue is surface exposed depends
/// on how many neighbors it has. For most applications, we have been using the approximation that residues with fewer
/// than 16 neighbors in a 10A radius are surface-exposed.
/// The original method of doing this was done with the tenA nb graph.  Then I switched to using cB-cB distances and
/// using the INTERACTION_RADIUS to define neighbors.  I'm reverting back to using the tenA neighbor graph.  The number
/// of neighbors variable is only used to determine whether a residue is surface exposed or not.  It's not like the detect
/// neighborship between Node objects method that effectively defines how much of the surface gets looked at per substitution.
/// Thus, using the tenA nb graph here is ok.
///
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::initialize_num_neighbors_counting_self() const {

	/* pose::Pose poseRef = get_surface_owner()->pose();

	conformation::Residue const & rsd1 = poseRef.residue( get_surface_owner()->rotamer_sets().moltenres_2_resid( parent::get_node_index() ) );
	Real distanceBetweenAtoms = 0.0;

	for ( Size ii=1; ii < poseRef.total_residue(); ++ii ) {

		if ( ii == get_surface_owner()->rotamer_sets().moltenres_2_resid( parent::get_node_index() ) ) { continue; }

		conformation::Residue const & rsd2 = poseRef.residue( ii );

		distanceBetweenAtoms = rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) );
		if ( distanceBetweenAtoms <= INTERACTION_RADIUS ) {
			num_neighbors_counting_self_++;
		}

	}

	// finally, count ourselves...
	num_neighbors_counting_self_++;
	*/

	core::scoring::TenANeighborGraph const & tenA_neighbor_graph( get_surface_owner()->pose().energies().tenA_neighbor_graph() );

	num_neighbors_counting_self_ =
		tenA_neighbor_graph.get_node( get_surface_owner()->rotamer_sets().moltenres_2_resid( parent::get_node_index() ) )->num_neighbors_counting_self();


#ifdef FILE_DEBUG
	TR_NODE << "initialize_num_neighbors_counting_self(): fc node index " << parent::get_node_index() << " (PDB: "
			<< wt_residue_for_node().seqpos() << "): " << num_neighbors_counting_self_ << " nbs" << std::endl;
#endif

}


///
/// @begin SurfaceNode::num_neighbors_counting_self
///
/// @brief
/// Returns the value stored in the member variable.
///
template < typename V, typename E, typename G >
int SurfaceNode< V, E, G >::num_neighbors_counting_self() const {

	if ( num_neighbors_counting_self_ == -1 ) {
		initialize_num_neighbors_counting_self();
	}

	return num_neighbors_counting_self_;
}


///
/// @begin SurfaceNode::init_hASA_variables
///
template < typename V, typename E, typename G >
void SurfaceNode< V, E, G >::init_hASA_variables() {
	curr_state_total_hASA_ = calculate_amount_total_hydrophobic_ASA();
	alt_state_total_hASA_ = curr_state_total_hASA_;

}

///
/// @begin SurfaceNode::calculate_amount_total_hydrophobic_ASA
///
/// @detailed
/// Iterates over all the edges emanating from a residues context graph and determines how many surface-exposed
/// total neighbors that residue has including itself.  Uses the context graph in the Pose
/// object because this function is only used for initialization.  After being calculated once, incrementing and
/// decrementing can be used instead of reiterating over all edges.
///
/// Bug source.  This function is not just used once during initialization.
///
/// As of 8/1/08, changing this function to iterate over all residues in the pose to determine the number of se
/// nbs instead of using the tenA neighbor graph.  This change is so that different interaction thresholds can be
/// used.  For example, say I want to find all patches within 8A, not 10A. Can't use the 10A to initialize with
/// or certain correctness checks fail in the Node methods.
///
/// Big change: Making this method be part of Node/BGNodes.  It's just too difficult to keep this method as part of the
/// SIG.  And I can't immediately think of a good reason to keep it in the SIG.  It's more a functionality which the
/// nodes should handle.
///
/// 02/20/09, Changing the method to return the total amount of hydrophobic accessibile surface area instead of the
/// number of hydrophobic neighbors. But most of the code remains the same.  We want to sum up the average residue hASA
/// for all of the surface-exposed neighbors within 10A.
///
template < typename V, typename E, typename G >
Real SurfaceNode<V, E, G>::calculate_amount_total_hydrophobic_ASA() {

	Real total_hASA = 0.0;

#ifdef FILE_DEBUG
	TR_NODE << "calculate_amount_total_hydrophobic_ASA(): ";
#endif

	if ( num_neighbors_counting_self() <= SURFACE_EXPOSED_CUTOFF ) {
		total_hASA += average_residue_hASA( wt_residue_for_node().aa(), num_neighbors_counting_self() );
	}

	// now iterate through all the other nodes.  have to use the nodes/bgnodes because they have se defined.
	// neighboring nodes should be reachable via edge objects that were init'd before.
	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {

		int other_node_index = get_incident_surface_edge(ii)->get_other_ind( parent::get_node_index() );

		// "search" through the map to make sure this edge is in the map (if not found, find returns map.end() which is what I check for)
		if ( fc_neighbor_map.find( std::pair<int,int>( parent::get_node_index(), other_node_index ) ) != fc_neighbor_map.end() ) {

#ifdef FILE_DEBUG
			TR_NODE << other_node_index << "-";
#endif
			int edge_index = ( other_node_index == get_incident_surface_edge(ii)->get_surface_node(0)->get_node_index() ? 0 : 1 );
			chemical::AA other_nodes_type = (get_incident_surface_edge(ii)->get_surface_node( edge_index )->wt_residue_for_node()).aa();
			Size nbs = get_incident_surface_edge(ii)->get_surface_node( edge_index )->num_neighbors_counting_self();
			Real hASA = average_residue_hASA( other_nodes_type, nbs );

#ifdef FILE_DEBUG
			TR_NODE << ":" << hASA;
#endif
			total_hASA += hASA;
#ifdef FILE_DEBUG
			TR_NODE << ", ";
#endif
		}
	}

	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		int other_node_index = get_edge_to_surface_bg_node(ii)->get_other_ind( this );
		// "search" through the map to make sure this edge is in the map (if not found, find returns map.end() which is what I check for)
		if ( bg_neighbor_map.find( std::pair<int,int>( parent::get_node_index(), other_node_index) ) != bg_neighbor_map.end() ) {

#ifdef FILE_DEBUG
			TR_NODE << "b" << other_node_index << "-";
#endif
			chemical::AA bg_nodes_type = (get_edge_to_surface_bg_node(ii)->get_surface_bg_node()->wt_residue_for_node()).aa();
			Size nbs = get_edge_to_surface_bg_node(ii)->get_surface_bg_node()->num_neighbors_counting_self();
			Real hASA = average_residue_hASA( bg_nodes_type, nbs );
#ifdef FILE_DEBUG
			TR_NODE << ":" << hASA;
#endif

			total_hASA += hASA;
#ifdef FILE_DEBUG
		TR_NODE << ", ";
#endif
		}
	}

	//if ( total_hASA > 1200 ) {
	//	utility_exit_with_message( "Nonsensical amount of hydrophobic accessible surface area found for residue. Quitting." );
	//}

#ifdef FILE_DEBUG
	TR_NODE << "returning patch area of " << total_hASA << " for node index " << parent::get_node_index()
			<< " (PDB: " << wt_residue_for_node().seqpos() << ")" << std::endl;
#endif
	return total_hASA;

}

///
/// @begin SurfaceNode::average_residue_hASA
///
/// @brief
/// Returns the amount of hydrophobic ASA this node/residue adds to a patch.
///
/// @detailed
/// This method is oblivious to whether a given residue should contribute to the patch area.  This function just does
/// the lookup for the hASA (using the SurfacePotential class).  Calls to this function should be surrounded by is-this-
/// a-surface-exposed-hydrophobic checks.
///
template < typename V, typename E, typename G >
Real SurfaceNode<V, E, G>::average_residue_hASA() const {
	return average_residue_hASA( get_rotamer( parent::get_current_state() )->aa(), num_neighbors_counting_self() /* this nodes num nbs */ );
}

///
/// @begin SurfaceNode::average_residue_hASA
///
/// @brief
/// Same function as above, but this one allows you to specify what type and num_nbs. This function is meant to be used
/// only when initializing all of the Nodes/BGNodes because at that point we have to use the wild-type residue type, not
/// the new rotamers residue type (which is what the above function will do).
///
template < typename V, typename E, typename G >
Real SurfaceNode<V, E, G>::average_residue_hASA( chemical::AA residue_type, Size num_nbs ) const {
	if ( num_nbs > (Size)BURIED_RESIDUE_NO_HSASA_CUTOFF ) return 0.00;
	return SurfacePotential::get_instance()->average_residue_hASA( residue_type, num_nbs );
}


///
/// @begin SurfaceNode::hASA_energy()
///
/// @brief
/// Calls a SurfacePotential class method to get the energy for a given patch area.
/// Only nodes with number of nbs less than or equal to SURFACE_EXPOSED_CUTOFF will get a score. All other nodes will
/// have a score of 0.0.
///
template < typename V, typename E, typename G >
Real SurfaceNode< V, E, G >::hASA_energy( Real patch_area ) const {
	if ( patch_area > MAX_PATCH_SURFACE_AREA ) {
		return surface_energy_weight_ * MAX_SURFACE_ENERGY;
	}
	if ( num_neighbors_counting_self() > SURFACE_EXPOSED_CUTOFF ) {
		return 0.0;
	}
	return SurfacePotential::get_instance()->hASA_patch_energy( patch_area, (Size)num_neighbors_counting_self() );
}



//----------------------------------------------------------------------------//
//----------------- Surface Background Node Class -------------------------//
//----------------------------------------------------------------------------//


template < typename V, typename E, typename G >
Real SurfaceBackgroundNode< V, E, G >::surface_energy_weight_ = 1.0;

template < typename V, typename E, typename G >
const int SurfaceBackgroundNode< V, E, G >::MAX_SURFACE_ENERGY = 100;


///
/// @begin SurfaceBackgroundNode::SurfaceBackgroundNode
///
/// @brief
/// main constructor.  No default constructor, copy constructor or assignment operator
///
/// @param
/// owner - [in] - the owning interaction graph
/// node_id - [in] - the index for this node amongst its owners set
///
template < typename V, typename E, typename G >
SurfaceBackgroundNode< V, E, G >::SurfaceBackgroundNode( AdditionalBackgroundNodesInteractionGraph< V, E, G >* owner, int node_index ) :
	BackgroundNode< V, E, G > ( owner, node_index ),
	curr_state_total_hASA_( 0.0 ),
	alt_state_total_hASA_( 0.0 ),
	have_prepared_for_simA_( false ),
	surface_exposed_( false ),
	is_below_buried_residue_no_hsasa_cutoff_( false ),
	num_neighbors_counting_self_(-1)
{
	surface_energy_weight_ = get_surface_owner()->surface_score_weight();
#ifdef FILE_DEBUG
	TR_BGNODE << "Setting surface_energy_weight to " << surface_energy_weight_ << std::endl;
#endif
}


///
/// @begin SurfaceBackgroundNode::~SurfaceBackgroundNode
///
/// @brief
/// destructor - no dynamicly allocated memory to manage -> empty
///
template < typename V, typename E, typename G >
SurfaceBackgroundNode< V, E, G >::~SurfaceBackgroundNode() {}


///
/// @begin SurfaceBackgroundNode::prepare_for_simulated_annealing
///
/// @brief
/// Sets the hASA of surface exposed hydrophobic neighbors using the pose's neighbor graph. There's
/// no need for BGNodes to keep track of which Nodes are really neighbors (energy graph vs neighbor graph) because they
/// never broadcast changes to their neighboring Nodes.
///
template < typename V, typename E, typename G >
void SurfaceBackgroundNode< V, E, G >::prepare_for_simulated_annealing() {

	if ( ! parent::get_edge_vector_up_to_date() )
		parent::update_edge_vector();

	// initialize the areas for this background node properly
	// the calculate ASA method assumes that BGNodes have edges to other FC nodes by this point.
	init_hASA_variables(); // make this a separate function so that we can redo this initialization if prep_for_simA gets called again

	have_prepared_for_simA_ = true;

}


///
/// @begin SurfaceBackgroundNode::detect_neighborship
///
/// @brief
/// returns true if this background residue neighbors with any rotamer of a SurfaceNode
///
/// @detailed
/// If this node is in the list of neighbors for the passed in SurfaceNode, return true.  Leave the responsibility
/// of figuring out neighborship to the first class node.
///
template < typename V, typename E, typename G >
bool SurfaceBackgroundNode< V, E, G >::detect_neighborship( SurfaceNode< V, E, G >* node ) const {
	return node->detect_neighborship_with_node( parent::get_node_index(), false /* i.e. this node is not a FC node */ );
}


///
/// @begin SurfaceBackgroundNode::project_surface_deltaE_for_substitution
///
/// @brief
/// returns the change in surface induced by a SurfaceNode undergoing a state substitution.
///
/// @detailed
/// Is it possible for a background node to have a state of zero (or unassigned state).  Not really - there's no state
/// for background nodes since they're not changing.
/// What happens if the first class node undergoing sub is in the unassigned state?  We have to figure out what the wild-type
/// residue at the changing position is.  Should I pass in the node index to this method?  I have to pass in a pointer
/// to the node because that's how to get the rotamer information stored there.  But get_node_index() is a protected
/// method so I can't call that using the pointer.
///
/// Bug fix: Turns out that if you do multiple loops of pack_rotamers, then the method blanket assign state 0 gets called.
/// Have to handle the case where the current state is some nonzero value and the alt_state is going to the zero state.
/// If not, you get seg faults.  But, if we're going into the unassigned state, does it make sense in actually calculating
/// a "correct" surface deltaE.  Doesn't seem like it.  Going to just add a instant return of a bad value!
///
/// There's a potential problem with the setting of the alt state count that can creep up on multiple rounds of consider()
/// without concomitant commit(). A separate method has been added to this class to handle these cases.  It basically
/// executes the following line each round of consider().
/// alt_state_num_se_hp_nbs_ = curr_state_num_se_hp_nbs_;
///
/// 02/20/2009/ Updating the function to do everything with hASA rather than hp counts.
///
template < typename V, typename E, typename G >
Real SurfaceBackgroundNode< V, E, G >::project_surface_deltaE_for_substitution(
	SurfaceNode< V, E, G >* fc_node_changing, int changing_nodes_curr_state, int changing_nodes_alt_state ) {

	// keep the value of the alt state in case this substitution is committed later. is this even used?
	// parent::set_alternate_state( parent::get_current_state() );
	// activating this as of 02/23/09 because parent::get_alternate_state() is called in at least one place
	// deactivating becase the method does not exist (at least in parent)

	// it may also be that the node undergoing substitution is in the unassigned state, right?
	// if that's the case, we need to look up the wild-type residue info
	//
	// the node undergoing sub is a first class node. so we need to get that nodes index so we can look it up
	// using the pose object.  get_index_of_adjacent_first_class_node is a method in BackgroundNode which looks in a vector
	// that holds all the adjacent first class node.  Only problem is that this particular bg node may have a number of adjacent
	// first class nodes.  How do we determine the index of the one that is changing?!?
	// get_node_index() is a protected method of Nodes!
	// Have to pass the index in to this method. Hopefully the BGEdge can get the index!

	parent::update_edge_vector();

	// for when blanket assign zero state gets called; it doesn't really matter what we return since it's an unassigned state
	if ( changing_nodes_alt_state == 0 ) {
		return MAX_SURFACE_ENERGY;
	}

	// If it's not a surface exposed residue then don't bother updating any of the hASA and stuff, just return 0.0.
	if ( !(this->is_surface_exposed()) ) {
		return 0.0;
	}

	// this is where we update the hASA for this node.
	// rather than recalculating the sum by going over all neighboring nodes, we can just take away the average contribution
	// to the patch of the current residue and add in the alt state residues contribution

	if ( changing_nodes_curr_state == 0 ) {
		// curr_state_total_hASA_ should have been init'd during the call to 'prep for simA'
		// likewise, alt_state = curr_state above

		// the hASA here should only be changed if the fc node that's changing is surface exposed
		// WRONG. the hASA can still change if the nb node isn't surface exposed. it might have more nbs than the surface
		// exposed cutoff but still not be completely buried.
		//if ( fc_node_changing->is_surface_exposed() ) {
		if ( fc_node_changing->is_below_buried_residue_no_hsasa_cutoff() ) {

			chemical::AA curr_AA = fc_node_changing->wt_residue_for_node().aa();
			chemical::AA alt_AA = fc_node_changing->get_rotamer( changing_nodes_alt_state )->aa();
			Size nbs = fc_node_changing->num_neighbors_counting_self();

			alt_state_total_hASA_ = curr_state_total_hASA_ - average_residue_hASA( curr_AA, nbs ) + average_residue_hASA( alt_AA, nbs );

			if ( alt_state_total_hASA_ < -0.001 ) {
				TR_BGNODE << "curr_state_total_hASA_: " << curr_state_total_hASA_ << ", alt_state_total_hASA_: " << alt_state_total_hASA_
					<< ", curr_AA: " << curr_AA << ", alt_AA: " << alt_AA << std::endl;
				this->print();
				utility_exit_with_message("Nonsensical alt_state reached (alt_state total hASA < 0). Exiting.");
			}
		}

	} else { // node wasn't in the unassigned state - 99% of calls will follow this branch

		// the hASA here should only be changed if the fc node that's changing is surface exposed
		//if ( fc_node_changing->is_surface_exposed() ) { // WRONG!
		if ( fc_node_changing->is_below_buried_residue_no_hsasa_cutoff() ) {

			chemical::AA curr_AA = fc_node_changing->get_rotamer( changing_nodes_curr_state )->aa();
			chemical::AA alt_AA = fc_node_changing->get_rotamer( changing_nodes_alt_state )->aa();
			Size nbs = fc_node_changing->num_neighbors_counting_self();

			alt_state_total_hASA_ = curr_state_total_hASA_ - average_residue_hASA( curr_AA, nbs ) + average_residue_hASA( alt_AA, nbs );
		}
	}

	Real difference = 0.0;

	if ( curr_state_total_hASA_ == alt_state_total_hASA_ ) { return 0.0; }  // not sure if this is possible anymore

	if ( curr_state_total_hASA_ < -0.001 ) { this->print(); utility_exit_with_message( "Nonsensical curr_state reached (curr_state total hASA < 0). Exiting." ); }
	if ( alt_state_total_hASA_ < -0.001 ) { this->print(); utility_exit_with_message( "Nonsensical alt_state reached (alt_state total hASA < 0). Exiting." ); }

	if ( ( curr_state_total_hASA_ > MAX_PATCH_SURFACE_AREA ) && ( alt_state_total_hASA_ > MAX_PATCH_SURFACE_AREA ) ) {
		// this sub could be good or bad, but let the contributions of the other nodes decide what happens
		return 0.0;
	}

	// if the alt state is greater than 800A2,
	if ( alt_state_total_hASA_ > MAX_PATCH_SURFACE_AREA ) {  // cutoff the distribution at 800A^2
		difference = surface_energy_weight_ * MAX_SURFACE_ENERGY;
#ifdef FILE_DEBUG
		TR_BGNODE << "project_surface_deltaE_for_substitution: bg node " << parent::get_node_index()
			<< " calculated deltaE for changing surface resid " << fc_node_changing->wt_seqpos_for_node() << "; "
			<< "curr state hASA: " << curr_state_total_hASA_ << ", alt hASA: " << alt_state_total_hASA_
			<< ", score difference: " << difference << std::endl;
#endif
		return difference;
	}
	// alt state must be 800A2 or less, which is in bounds

	if ( curr_state_total_hASA_ > MAX_PATCH_SURFACE_AREA ) {
		// then, the current state has less hASA than the alt state, so this would be a favorable change
		difference = surface_energy_weight_ * -1 * MAX_SURFACE_ENERGY;
#ifdef FILE_DEBUG
	TR_BGNODE << "project_surface_deltaE_for_substitution: bg node " << parent::get_node_index()
		<< " calculated deltaE for changing surface resid " << fc_node_changing->wt_seqpos_for_node() << "; "
		<< "curr state hASA: " << curr_state_total_hASA_ << ", alt hASA: " << alt_state_total_hASA_
		<< ", score difference: " << difference << std::endl;
#endif
		return difference;
	}

	difference = surface_energy_weight_ * ( hASA_energy( alt_state_total_hASA_ ) - hASA_energy( curr_state_total_hASA_ ) );

#ifdef FILE_DEBUG
	TR_BGNODE << "project_surface_deltaE_for_substitution: bg node " << parent::get_node_index()
		<< " calculated deltaE for changing surface resid " << fc_node_changing->wt_seqpos_for_node() << "; "
		<< "curr state hASA: " << curr_state_total_hASA_ << ", alt hASA: " << alt_state_total_hASA_
		<< ", score difference: " << difference << std::endl;
#endif

	return difference;

}


///
/// @begin SurfaceBackgroundNode::acknowledge_substitution
///
/// @brief
/// bookkeeping to reflect a SurfaceNode's state substitution.
///
template < typename V, typename E, typename G >
void SurfaceBackgroundNode< V, E, G >::acknowledge_substitution_surface() {
	curr_state_total_hASA_ = alt_state_total_hASA_;
}


///
/// @begin SurfaceBackgroundNode< V, E, G >::get_surface_score
///
/// @brief
/// returns the surface score under the current state assignment
///
template < typename V, typename E, typename G >
Real SurfaceBackgroundNode< V, E, G >::get_surface_score() const {

	if ( !(this->is_surface_exposed()) ) {
		return 0.0;
	}

	if ( curr_state_total_hASA_ > MAX_PATCH_SURFACE_AREA ) {
		return surface_energy_weight_ * MAX_SURFACE_ENERGY;
	}

	return surface_energy_weight_ * hASA_energy( curr_state_total_hASA_ );
}


///
/// @begin SurfaceBackgroundNode::print
///
/// @brief
/// used only for debugging
///
template < typename V, typename E, typename G >
void SurfaceBackgroundNode< V, E, G >::print() const {

	std::cout << "(bgnode id:" << parent::get_node_index() << ", hASA: " << average_residue_hASA()
		<< ", curr: " << curr_state_total_hASA_ << ", alt: " << alt_state_total_hASA_
		<< ", score: " << get_surface_score() << ", se: " << surface_exposed_ << std::endl;
}

///
/// @begin SurfaceBackgroundNode::count_static_memory
///
template < typename V, typename E, typename G >
unsigned int SurfaceBackgroundNode< V, E, G >::count_static_memory() const {
	return sizeof( SurfaceBackgroundNode< V, E, G > );
}


///
/// @begin SurfaceBackgroundNode::count_dynamic_memory
///
template < typename V, typename E, typename G >
unsigned int SurfaceBackgroundNode< V, E, G >::count_dynamic_memory() const {
	unsigned int total_memory = parent::count_dynamic_memory();

	// any other dynamic memory that needs to be counted?
	return total_memory;
}


///
/// @begin SurfaceBackgroundNode::is_surface_exposed
///
/// @brief
/// Returns the value of surface_exposed_ which gets set during the SIG initialize() method.
///
template < typename V, typename E, typename G >
bool SurfaceBackgroundNode< V, E, G >::is_surface_exposed() const {
	return surface_exposed_;
}


///
/// @begin SurfaceBackgroundNode::surface_exposed
///
/// @brief
/// using this for debugging
///
template < typename V, typename E, typename G >
void SurfaceBackgroundNode< V, E, G >::surface_exposed( bool value ) {
	surface_exposed_ = value;
}


///
/// @begin SurfaceBackgroundNode::is_below_buried_residue_no_hsasa_cutoff
///
/// @brief
/// Returns the value of is_below_buried_residue_no_hsasa_cutoff_ which gets set during the SIG initialize() method.
///
template < typename V, typename E, typename G >
bool SurfaceBackgroundNode< V, E, G >::is_below_buried_residue_no_hsasa_cutoff() const {
	return is_below_buried_residue_no_hsasa_cutoff_;
}

///
/// @begin SurfaceBackgroundNode::is_below_buried_residue_no_hsasa_cutoff
///
/// @brief
/// setter for the is_below_buried_residue_no_hsasa_cutoff_ bool
///
template < typename V, typename E, typename G >
void SurfaceBackgroundNode< V, E, G >::is_below_buried_residue_no_hsasa_cutoff( bool value ) {
	is_below_buried_residue_no_hsasa_cutoff_ = value;
}


///
/// @begin SurfaceBackgroundNode::reset_alt_state_total_hASA
///
/// @brief
/// Sets the alt state total hASA to the current state total hASA.  See comments in SIG and commit_considered_substitution_surface
/// for more information about why this method exists.
///
template < typename V, typename E, typename G >
void SurfaceBackgroundNode< V, E, G >::reset_alt_state_total_hASA() {

	if ( alt_state_total_hASA_ != curr_state_total_hASA_ ) {
		alt_state_total_hASA_ = curr_state_total_hASA_;
#ifdef FILE_DEBUG
		TR_BGNODE << "reset_alt_state_total_hASA(): alt state total hASA was not equal to current state total hASA, but is now." << std::endl;
#endif
	}

}


///
/// @begin SurfaceBackgroundNode::initialize_num_neighbors_counting_self
///
/// @brief
/// Returns the number of neighbors within the interaction threshold, counting self.  A member variable stores the
/// value so that this method only has to run once and lookups can be performed thereafter.
///
/// @detailed
/// See comments in SurfaceNode::initialize_num_neighbors_counting_self for details.
///
template < typename V, typename E, typename G >
void SurfaceBackgroundNode< V, E, G >::initialize_num_neighbors_counting_self() const {

	/*pose::Pose poseRef = get_surface_owner()->pose();

	conformation::Residue const & rsd1 = poseRef.residue( get_surface_owner()->bg_node_2_resid( parent::get_node_index() ) );
	Real distanceBetweenAtoms = 0.0;

	for ( Size ii=1; ii < poseRef.total_residue(); ++ii ) {

		if ( ii == get_surface_owner()->rotamer_sets().moltenres_2_resid( parent::get_node_index() ) ) { continue; }

		conformation::Residue const & rsd2 = poseRef.residue( ii );

		distanceBetweenAtoms = rsd1.xyz( rsd1.nbr_atom() ).distance( rsd2.xyz( rsd2.nbr_atom() ) );
		if ( distanceBetweenAtoms <= INTERACTION_RADIUS ) {
			num_neighbors_counting_self_++;
		}

	}

	// finally, count ourselves...
	num_neighbors_counting_self_++;
	*/

	core::scoring::TenANeighborGraph const & tenA_neighbor_graph( get_surface_owner()->pose().energies().tenA_neighbor_graph() );

	num_neighbors_counting_self_ = tenA_neighbor_graph.get_node( get_surface_owner()->bg_node_2_resid( parent::get_node_index() ) )->num_neighbors_counting_self();

#ifdef FILE_DEBUG
	TR_BGNODE << "initialize_num_neighbors_counting_self(): bg node index " << parent::get_node_index() << " (PDB: "
			<< get_surface_owner()->bg_node_2_resid( parent::get_node_index() )
			<< "): " << num_neighbors_counting_self_ << " nbs" << std::endl;
#endif

}


///
/// @begin SurfaceBackgroundNode::num_neighbors_counting_self
///
/// @brief
/// Returns the value stored in the member variable.
///
template < typename V, typename E, typename G >
int SurfaceBackgroundNode< V, E, G >::num_neighbors_counting_self() const {

	if ( num_neighbors_counting_self_ == -1 ) {
		initialize_num_neighbors_counting_self();
	}

	return num_neighbors_counting_self_;
}


///
/// @begin SurfaceBackgroundNode::init_hASA_variables
///
template < typename V, typename E, typename G >
void SurfaceBackgroundNode< V, E, G >::init_hASA_variables() {
	curr_state_total_hASA_ = calculate_amount_total_hydrophobic_ASA();
	alt_state_total_hASA_ = curr_state_total_hASA_;

}


///
/// @begin SurfaceBackgroundNode::calculate_amount_total_hydrophobic_ASA
///
/// @detailed
/// See comments in SurfaceNode for more info.
///
template < typename V, typename E, typename G >
Real SurfaceBackgroundNode<V, E, G>::calculate_amount_total_hydrophobic_ASA() {

	Real total_hASA = 0.0;

#ifdef FILE_DEBUG
	TR_BGNODE << "calculate_amount_total_hydrophobic_ASA(): ";
#endif
	if ( num_neighbors_counting_self() <= SURFACE_EXPOSED_CUTOFF ) {
		total_hASA += average_residue_hASA( wt_residue_for_node().aa(), num_neighbors_counting_self() );
	}

	// now iterate through all the other nodes.  have to use the nodes/bgnodes because they have se defined.
	// neighboring nodes should be reachable via edge objects that were init'd before.
	// only need to iterate over fc nodes edging with this node, because bgnode-bgnode edges don't exist
	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {

		// num_incident_edges may contain more edges than just the ones we want here, so call detect neighborship to determine
		// if we're in the cutoff since the neighbor_map doesn't exist in this class
		if ( get_surface_bg_edge(ii)->get_surface_node()->detect_neighborship_with_node( parent::get_node_index(), false ) ) {

#ifdef FILE_DEBUG
			int other_node_index = get_surface_bg_edge(ii)->get_other_ind( this );
			TR_BGNODE << other_node_index << "-";
#endif
			chemical::AA fc_nodes_type = ( get_surface_bg_edge(ii)->get_surface_node()->wt_residue_for_node() ).aa();
			Size nbs = get_surface_bg_edge(ii)->get_surface_node()->num_neighbors_counting_self();
			Real hASA = average_residue_hASA( fc_nodes_type, nbs );
#ifdef FILE_DEBUG
			TR_BGNODE << ":" << hASA;
#endif

			total_hASA += hASA;
#ifdef FILE_DEBUG
		TR_BGNODE << ", ";
#endif
		}
	}

	//if ( total_hASA > 1200 ) {
	//	TR_BGNODE << "calculate_amount_total_hydrophobic_ASA(): calculating patch area for residue " << wt_residue_for_node().name3() << " "
	//			<< get_surface_owner()->bg_node_2_resid( parent::get_node_index() ) << std::endl;
	//	print();
	//}

#ifdef FILE_DEBUG
	TR_BGNODE << "returning patch area of " << total_hASA << " for bg node index " << parent::get_node_index()
			<< " (PDB " << get_surface_owner()->bg_node_2_resid( parent::get_node_index() ) << ")" << std::endl;
#endif
	return total_hASA;

}


///
/// @begin SurfaceBackgroundNode::average_residue_hASA
///
/// @brief
/// Returns the amount of hydrophobic ASA this residue adds to a patch.
///
/// @detailed
/// This method is oblivious to whether a given residue should contribute to the patch area.  This function just does
/// the lookup for the hASA (using the SurfacePotential class).  Calls to this function should be surrounded by is-this-
/// a-surface-exposed-hydrophobic checks.
/// Since background nodes don't actually change type, whatever the wild type residue's aa type is what we use for the score.
///
template < typename V, typename E, typename G >
Real SurfaceBackgroundNode<V, E, G>::average_residue_hASA() const {
	return average_residue_hASA( wt_residue_for_node().aa(), num_neighbors_counting_self() /* this nodes num nbs */ );
}

///
/// @begin SurfaceBackgroundNode::average_residue_hASA
///
/// @brief
/// Same function as above, but this one allows you to specify what type and num_nbs. This function is meant to be used
/// only when initializing all of the Nodes/BGNodes because at that point we have to use the wild-type residue type, not
/// the new rotamers residue type (which is what the above function will do).
///
template < typename V, typename E, typename G >
Real SurfaceBackgroundNode<V, E, G>::average_residue_hASA( chemical::AA residue_type, Size num_nbs ) const {
	if ( num_nbs > (Size)BURIED_RESIDUE_NO_HSASA_CUTOFF ) return 0.00;
	return SurfacePotential::get_instance()->average_residue_hASA( residue_type, num_nbs );
}


///
/// @begin SurfaceBackgroundNode::hASA_energy()
///
/// @brief
/// Calls a SurfacePotential class method to get the energy for a given patch area.
/// Not sure if I should take num neighbors as a parameter to the function, or just use the class method
/// num_neighbors_counting_self().
///
template < typename V, typename E, typename G >
Real SurfaceBackgroundNode< V, E, G >::hASA_energy( Real patch_area ) const {
	if ( patch_area > MAX_PATCH_SURFACE_AREA ) {
		return surface_energy_weight_ * MAX_SURFACE_ENERGY;
	}
	if ( num_neighbors_counting_self() > SURFACE_EXPOSED_CUTOFF ) {
		return 0.0;
	}
	return SurfacePotential::get_instance()->hASA_patch_energy( patch_area, (Size)num_neighbors_counting_self() );
}



//----------------------------------------------------------------------------//
//-------------------------- Surface Edge Class ---------------------------//
//----------------------------------------------------------------------------//

///
/// @begin SurfaceEdge::SurfaceEdge
///
/// @brief
/// main constructor.  No default, or copy constructors, no assignment operator
///
/// @param
/// owner - [in] - the owning interaction graph object
/// node1 - [in] - the index of the lower-indexed SurfaceNode
/// node2 - [in] - the index of the higher-indexed SurfaceNode
///
template < typename V, typename E, typename G >
SurfaceEdge< V, E, G >::SurfaceEdge( G* owner, int node1, int node2 ) :
	FirstClassEdge< V, E, G > ( owner, node1, node2 ),
	node_changing_( -1 ),
	node_not_changing_( -1 )
{

	nodes_curr_states_[ 0 ] = nodes_curr_states_[ 1 ] = 0;
	nodes_alt_states_[ 0 ] = nodes_alt_states_[ 1 ] = 0;

	//variables for tracking the magnitude of the surface_deltaE's
	max_surface_deltaE_last_50_commits_[0] = max_surface_deltaE_last_50_commits_[1] = -1.0f;
	max_surface_deltaE_recent_50_commits_[0] = max_surface_deltaE_recent_50_commits_[1] = -1.0f;
	magnitude_last_surface_deltaE_[0] = magnitude_last_surface_deltaE_[1] = -1.0f;
	num_surface_deltaE_observations_since_update_[0] = num_surface_deltaE_observations_since_update_[1] = 0;

}


///
/// @begin SurfaceEdge::~SurfaceEdge
///
template < typename V, typename E, typename G >
SurfaceEdge< V, E, G >::~SurfaceEdge() {}


///
/// @begin SurfaceEdge::prepare_for_simulated_annealing
///
/// @brief
/// drops zero submatrices of the AminoAcidNeighborSparseMatrix and if the two_body_energies_ member then holds nothing,
/// it checks whether or not its incident nodes have any sphere overlaps.  If they don't then the edge deletes itself.
///
template < typename V, typename E, typename G >
void SurfaceEdge< V, E, G >::prepare_for_simulated_annealing() {

	parent::prepare_for_simulated_annealing_no_deletion();

	if ( parent::pd_edge_table_all_zeros() ) {

		// if the two fc nodes aren't neighbors, then we don't need to keep this edge
		if ( !( get_surface_node(0)->detect_neighborship_with_node( parent::get_node_index(1), true ) ) ) {

		// if the two fc nodes aren't neighbors, then we don't need to keep this edge
//#ifdef FILE_DEBUG
//			TR_EDGE << "prepare_for_simulated_annealing - dropping edge e(" << parent::get_node_index( 0 )
//				<< "," << parent::get_node_index( 1 ) << ")" << std::endl;
//#endif
			delete this;
		}
	}
}


///
/// @begin SurfaceEdge::acknowledge_state_zeroed_surface
///
/// @brief
/// respond to when one of its vertices enters the "unassigned" state
///
/// @param
/// node_that_changed - [in] - the index of the node that changed
///
template < typename V, typename E, typename G >
void SurfaceEdge< V, E, G >::acknowledge_state_zeroed_surface( int node_that_changed ) {

	// node_changing_ is the node of this edge that is changing
	node_changing_ = ( node_that_changed == parent::get_node_index(0) ? 0 : 1 );
	node_not_changing_ = !node_changing_;

	nodes_curr_states_[ node_changing_ ] = 0;
	inform_non_changing_node_of_neighbors_change();
}


///
/// @begin SurfaceEdge::get_surface_deltaE_for_neighbor
///
/// @brief
/// returns the change in surface score for the neighbor of a node that is
/// produced by the state substitution it is considering.
///
/// @detailed
/// Need to tell the node that's not changing which node is considering the change and what the alternate state it
/// might change to is.  Then that node can recompute it's change in se hp nb count and return a change in score.
///
///
template < typename V, typename E, typename G >
Real SurfaceEdge< V, E, G >::get_surface_deltaE_for_neighbor( int node_considering_substitution, int alt_state ) {

	node_changing_ = ( node_considering_substitution == parent::get_node_index(0) ? 0 : 1 );
	node_not_changing_ = ! node_changing_;

	// save the index to the alt_state in a local array called nodes_alt_states_
	// the nodes_curr_states_ and nodes_alt_state_ are int arrays of size 2 that just hold the rotamer id (state) for
	// convenient access
	nodes_alt_states_[ node_changing_ ] = alt_state;
	nodes_alt_states_[ node_not_changing_ ] = nodes_curr_states_[ node_not_changing_ ];

	// nodes_curr_states_[ node_changing_ ] is the current state at the node considering substitution

	Real surface_deltaE = get_surface_node( node_not_changing_ )->get_surface_deltaE_for_neighbors_state_substitution(
			get_surface_node( node_changing_ ), nodes_curr_states_[ node_changing_ ], alt_state );

	magnitude_last_surface_deltaE_[ node_changing_ ] = std::abs( surface_deltaE );

#ifdef FILE_DEBUG
	// this line just generates waaaay too much output
	//TR_EDGE << "get_surface_deltaE_for_neighbor() for changing node: " << node_considering_substitution
	//	<< " and node: " << parent::get_node_index( node_not_changing_ ) << " returning deltaE of " << surface_deltaE << std::endl;
#endif

	return surface_deltaE;

}


///
/// @begin SurfaceEdge::acknowledge_substitution_surface
///
/// @brief
/// bookkeeping following the decision to substitute a nodes current state with the alternate it was asked to consider.
///
template < typename V, typename E, typename G >
void SurfaceEdge< V, E, G >::acknowledge_substitution_surface() {

	inform_non_changing_node_of_neighbors_change();

	nodes_curr_states_[ node_changing_ ] = nodes_alt_states_[ node_changing_ ];

	track_max_magnitude_surface_deltaE();

	return;
}

///
/// @begin SurfaceEdge::inform_non_changing_node_of_neighbors_change
///
/// @brief
/// tells the node that isn't considering a substitution or changing state
/// that its neighbor who is has changed.
///
template < typename V, typename E, typename G >
inline
void SurfaceEdge< V, E, G >::inform_non_changing_node_of_neighbors_change() {
	get_surface_node( node_not_changing_ )->acknowledge_neighbors_substitution_surface();
}



///
/// @begin SurfaceEdge::track_max_magnitude_surface_deltaE
///
/// @brief
/// Keeps track of the maximum surface deltaE seen in the last 50 commits
///
template < typename V, typename E, typename G >
void SurfaceEdge< V, E, G >::track_max_magnitude_surface_deltaE() {

	++num_surface_deltaE_observations_since_update_[ node_changing_ ];

	//TR_EDGE << "(e: " << get_node_index( 0 ) << " " << get_node_index( 1 );
	//TR_EDGE << "-- "<< get_node_index( node_changing_ ) << ", " <<  num_surface_deltaE_observations_since_update_[ node_changing_ ] << ") ";

	if ( magnitude_last_surface_deltaE_[ node_changing_ ] > max_surface_deltaE_recent_50_commits_[ node_changing_ ] ) {
		max_surface_deltaE_recent_50_commits_[ node_changing_ ] = magnitude_last_surface_deltaE_[ node_changing_ ];
	}

	if ( num_surface_deltaE_observations_since_update_[ node_changing_] == 50 ) {
		//TR_EDGE << "Edge: " << get_node_index( 0 ) << " " << get_node_index( 1 );
		//TR_EDGE << " max_surface_delatE for node "<< get_node_index( node_changing_ ) << ": " << max_surface_deltaE_recent_50_commits_[ node_changing_ ] << std::endl;
		max_surface_deltaE_last_50_commits_[ node_changing_ ] = max_surface_deltaE_recent_50_commits_[ node_changing_ ];
		num_surface_deltaE_observations_since_update_[ node_changing_] = 0;
		max_surface_deltaE_recent_50_commits_[ node_changing_ ] = -1.0f;
	}
}


///
/// @begin SurfaceEdge::declare_energies_final
///
/// @brief
/// Reduces memory usage in the two body energy table after the energy
/// calculating function declares that the energies will not change thereafter
///
/// @remarks (all by apl)
/// In the PDEdge's version of this method, after invoking two_body_energies_.drop_zero_submatrices_where_possible();
/// the PDEdge checks if the two body energy table it now holds is empty.  If the table is empty, the edge deletes itself.
///
/// A SASAEdge should not delete itself if the pair energies are all zero since the Minkowski sum of a water and a van
/// der Waal's sphere extends further out from an atoms center than its (lj_atr, lj_rep, lksolv) interaction sphere.
/// However, if a SASAEdge holds no pair energies, it's a very good candidate for removal  -- it just first needs to check
/// that no (vdw + 1.4 A) spheres overlap between any pair of rotamers on the edges it connects.
///
template < typename V, typename E, typename G >
void SurfaceEdge< V, E, G >::declare_energies_final() {
	parent::declare_energies_final_no_deletion();
}


///
/// @begin SurfaceEdge::get_max_surface_deltaE_guess
///
template < typename V, typename E, typename G >
Real SurfaceEdge< V, E, G >::get_max_surface_deltaE_guess( int /*node_index*/ ) const {
	return max_surface_deltaE_last_50_commits_[ node_changing_ ];
}


///
/// @begin SurfaceEdge::getMemoryUsageInBytes
///
/// @remarks
/// Not implemented.
///
template < typename V, typename E, typename G >
unsigned int SurfaceEdge< V, E, G >::getMemoryUsageInBytes() const {
	return 0;
}


///
/// @begin SurfaceEdge::count_static_memory
///
template < typename V, typename E, typename G >
unsigned int SurfaceEdge< V, E, G >::count_static_memory() const {
	return sizeof( SurfaceEdge< V, E, G > );
}


///
/// @begin SurfaceEdge::count_dynamic_memory
///
template < typename V, typename E, typename G >
unsigned int SurfaceEdge< V, E, G >::count_dynamic_memory() const {

	unsigned int total_memory = parent::count_dynamic_memory();

	// any additional memory that needs to be counted?
	return total_memory;
}



//----------------------------------------------------------------------------//
//-------------------- SurfaceBackgroundEdge Class ------------------------//
//----------------------------------------------------------------------------//

///
/// @begin SurfaceBackgroundEdge< V, E, G >::SurfaceBackgroundEdge
///
/// @brief
/// main constructor
///
/// @param
/// owner - [in] - the owning graph
/// first_class_node_index - [in] - the index of the first class node upon which this new edge is incident
/// second_class_node_index - [in] - the index of the second class node upon which this new edge is incident
///
template < typename V, typename E, typename G >
SurfaceBackgroundEdge< V, E, G >::SurfaceBackgroundEdge
	( AdditionalBackgroundNodesInteractionGraph < V, E, G >* owner, int first_class_node_index, int background_node_index ):
	BackgroundToFirstClassEdge< V, E, G >( owner, first_class_node_index, background_node_index ),
	fc_node_curr_state_( 0 ),
	fc_node_alt_state_( 0 )
{
	max_surface_deltaE_last_50_commits_ = -1.0f;
	max_surface_deltaE_recent_50_commits_ = -1.0f;
	magnitude_last_surface_deltaE_ = -1.0f;
	num_surface_deltaE_observations_since_update_ = 0;
}


///
/// @begin SurfaceBackgroundEdge::~SurfaceBackgroundEdge
///
/// @brief
///
template < typename V, typename E, typename G >
SurfaceBackgroundEdge< V, E, G >::~SurfaceBackgroundEdge() {}


///
/// @begin SurfaceBackgroundEdge::prepare_for_simulated_annealing
///
/// @brief
/// Invoked by AdditionalBackgroundNodesInteractionGraph::prepare_for_simulated_annealing.
///
/// @remarks
/// The SurfaceBackgroundEdge has no responsibilities in this function.  However,
/// when the AdditionalBackgroundNodesInteractionGraph invokes
/// prepare_for_simulated_annealing on the SurfaceBackgroundNode that this edge
/// is incident upon, that node will invoke initialize_overlap_cache on
/// this edge
///
template < typename V, typename E, typename G >
void SurfaceBackgroundEdge< V, E, G >::prepare_for_simulated_annealing() {}


///
/// @begin SurfaceBackgroundEdge::acknowledge_state_change
///
/// @brief
/// bookkeeping in response to a SurfaceNode switching states (without having
/// gone through the usual consider-substitution/commit-substitution pattern).
///
/// @param
/// new_state - [in] - the state the SurfaceNode switched to
///
template < typename V, typename E, typename G >
void SurfaceBackgroundEdge< V, E, G >::acknowledge_state_change( int new_state ) {

	if ( new_state == fc_node_curr_state_ )
		return;

	get_surface_deltaE_for_substitution( new_state );
	acknowledge_substitution_surface();
}


///
/// @begin SurfaceBackgroundEdge::get_surface_deltaE_for_substitution
///
/// @brief
/// returns the change in surface energy produced by a background node in response to a considered state
/// substitution of the first class node
///
/// @detailed
/// Need to tell the bg node that's not changing that the fc node is considering a change and what the alternate state it
/// might change to is.  Then that node can recompute it's change in se hp nb total hASA and return a change in score.
/// Called once by acknowledge_state_change method in this class when blanket setting the state of all the nodes
/// to zero during initialization.
///
/// @param
/// alt_state - [in] - the alternate state (i.e. rotamer id) that the first class node is considering
///
template < typename V, typename E, typename G >
Real SurfaceBackgroundEdge< V, E, G >::get_surface_deltaE_for_substitution( int alt_state ) {

	// keep the alt_state around in case the fc node commits this sub being considered
	fc_node_alt_state_ = alt_state;

	// so what we need to do here is update the alt state hASA based on the first class nodes alt state.
	// do we want the edge to handle this update or leave it to the node? for the surface (first class) edge class
	// method comparable to this one, we passed a reference to the first class node and current state of the node that
	// was changing. the current state (and alt_state) of the changing node is kept as a member variable in this class.

	Real const surface_deltaE = get_surface_bg_node()->project_surface_deltaE_for_substitution(
			get_surface_node(), fc_node_curr_state_, fc_node_alt_state_ );

	magnitude_last_surface_deltaE_ = std::abs( surface_deltaE );

	return surface_deltaE;

}


///
/// @begin SurfaceBackgroundEdge::acknowledge_substitution_surface
///
/// @brief
/// bookkeeping in response to the SurfaceNode committing the considered substitution
///
template < typename V, typename E, typename G >
void SurfaceBackgroundEdge< V, E, G >::acknowledge_substitution_surface() {

	fc_node_curr_state_ = fc_node_alt_state_;

	get_surface_bg_node()->acknowledge_substitution_surface();

	track_max_magnitude_surface_deltaE();
}


///
/// @begin SurfaceBackgroundEdge::track_max_magnitude_surface_deltaE
///
/// @brief
/// Keeps track of the maximum surface deltaE seen in the last 50 commits. Used possibly for procrastinating
/// recalculations of the surface score.
///
template < typename V, typename E, typename G >
void SurfaceBackgroundEdge< V, E, G >::track_max_magnitude_surface_deltaE() {

	++num_surface_deltaE_observations_since_update_;

	if ( magnitude_last_surface_deltaE_ > max_surface_deltaE_recent_50_commits_) {
		max_surface_deltaE_recent_50_commits_ = magnitude_last_surface_deltaE_;
	}

	if ( num_surface_deltaE_observations_since_update_ == 1000 ) {
		max_surface_deltaE_last_50_commits_ = max_surface_deltaE_recent_50_commits_;
		num_surface_deltaE_observations_since_update_= 0;
		max_surface_deltaE_recent_50_commits_ = -1.0f;
	}
}


///
/// @begin SurfaceBackgroundEdge::get_max_surface_deltaE_guess
///
/// @brief
/// Returns the value of max_surface_deltaE_last_50_commits_
///
template < typename V, typename E, typename G >
Real SurfaceBackgroundEdge< V, E, G >::get_max_surface_deltaE_guess() const {
	return max_surface_deltaE_last_50_commits_;
}


///
/// @begin SurfaceBackgroundEdge::count_static_memory
///
template < typename V, typename E, typename G >
unsigned int SurfaceBackgroundEdge< V, E, G >::count_static_memory() const {
	return sizeof( SurfaceBackgroundEdge< V, E, G > );
}


///
/// @begin SurfaceBackgroundEdge::count_dynamic_memory
///
template < typename V, typename E, typename G >
unsigned int SurfaceBackgroundEdge< V, E, G >::count_dynamic_memory() const {

	unsigned int total_memory = parent::count_dynamic_memory();
	return total_memory;
}



//----------------------------------------------------------------------------//
//--------------------- Surface Interaction Graph -------------------------//
//----------------------------------------------------------------------------//


///
/// @begin SurfaceInteractionGraph::SurfaceInteractionGraph
///
/// @brief
/// Main constructor
///
/// @param
/// num_nodes - [in] - the number of first class nodes in the graph.  The number
///   of second class nodes is reported to the graph much later.
///
template < typename V, typename E, typename G >
SurfaceInteractionGraph< V, E, G >::SurfaceInteractionGraph( int num_nodes ) :
	AdditionalBackgroundNodesInteractionGraph< V, E, G > ( num_nodes ),
	num_total_residues_( 0 ),
	num_residues_assigned_as_background_( 0 ),
	num_commits_since_last_update_( 0 ),
	total_energy_current_state_assignment_( 0 ),
	total_energy_alternate_state_assignment_( 0 ),
	node_considering_alt_state_( 0 ),
	deltaE_threshold_for_avoiding_surface_calcs_( -1.0f ),
	prepared_for_simulated_annealing_( false ),
	surface_score_weight_( 1 )
{}


///
/// @begin SurfaceInteractionGraph::~SurfaceInteractionGraph
///
template < typename V, typename E, typename G >
SurfaceInteractionGraph< V, E, G >::~SurfaceInteractionGraph() {}


///
/// @begin SurfaceInteractionGraph::initialize()
///
/// @detailed
/// This initialize() method needs to do two things: 1) it needs to set residues that are not changing as background
/// residues and 2) it needs to save references to the rotamer set or sets so that nodes can later figure out whether
/// a particular rotamer change is a hydrophobic or not.
///
/// In ++, there's a file InteractionGraphSupport.cc which is an analog of the InteractionGraphFactory in mini.  In ++,
/// the InteractionGraphSupport file instantiates and initializes, depending on the command line switches, the right
/// interaction graph.  For the SASAInteractionGraph, it first initializes the PDInteractionGraph (which is the base)
/// and then calls the SASAIG initialize method.
///
/// The thing is that this initialize method can't be called when the graph is constructed in the InteractionGraphFactory.
/// The reason is that the PDInteractionGraph base initialize() method is NOT called until later in the pack rotamers
/// process.  (Actually it's called from within rotsets->compute_energies().)  Only after the rotsets->compute energies
/// returns can I stick in an initialize() method call for the SurfaceInteractionGraph (SIG).  But then that's too late because
/// the rotsets object has already computed some energies. Perhaps that's ok though. The rotsets object only calculates
/// the PD energy terms - it doesn't do anything with non-PD terms.
///
/// If a SIG init method is called during construction of the SIG, then the init method that's called by the rotsets object
/// causes all the node info that was set in the SIG init method to be lost because the rotsets init method recreates all
/// the nodes in the interaction graph when it runs. (That was a fun behaviour to figure out.)
///
/// So the solution I'm going for is to call this init method in the prepare for simulated annealing method of this class.
/// That gets called just before SA starts, so it will do tasks 1) and 2) above then.  It doesn't really matter when the
/// tasks get done as long as it happens before SA starts.  This also means that the SIG will now have to keep a reference
/// to the Pose, the Task, and the RotamerSets objects since it needs all of these things to do tasks 1) and 2).
///
///
///	prepare_for_simulated_annealing gets called by the FixbbSA::run() method.  Before this method, the
/// rotamersets object has called compute_energies() (the whole process being started in pack_rotamers)
/// which calls initialize() on the IG.  I need to place the SIG init method directly after the IG init
/// method that the RS object calls.
///
/// Older comments:
/// This interaction graph is templated so that it can be used with the original PD or the OTF graphs.  If it's an OTF
/// graph, then a reference to the pose object is available with the pose() method.  However, I can't rely on having that
/// method for the case when an OTF graph is not used. Thus, take the pose reference in as a parameter and refactor later
/// if necessary.
///
/// Newest comments:
/// Don't call the parent classes initialize method because that calls create_new_node for all designable position, making
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
/// SurfaceExposedNodesInteractionGraph and then make the SurfaceInteractionGraph inherit that one.  Then
/// "FirstClassNodes" would be nodes that are designable and surface-exposed; all other nodes would be "BackgroundNodes".
///
/// Since it would require alot of work to write another derived class for probably not that significant (or necessary)
/// performance improvement, I'll just make packable residues be FC nodes also.  So continue to call the parent classes
/// initialize() method.  Just set residues which are not packable nor designable to BG nodes.
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::initialize( rotamer_set::RotamerSetsBase const & rot_sets_base ) {

	/// TEMP!!!!
	// APL wants to create a parent/base class for RotamerSets called RotamerSetsBase which will hold variables and functions
	// that are common to RotamerSets and FlexbbRotamerSets. Each interaction graph will be passed a RotamerSetsBase object
	// to its initialize method, and functions in this class should use only access common base-class functions.  However,
	// this would require changing quite a bit of code and making sure that only base-class RotamerSetsBase methods are used
	// here which I'm not interested in doing right now. So we can cast the RotamerSetsBase object to a RotamerSets object
	// and use it the way it was previously.
	rotamer_set::RotamerSets const & rot_sets( static_cast< rotamer_set::RotamerSets const & > (rot_sets_base) );

#ifdef FILE_DEBUG
	TR_SIG << "initialize() called" << std::endl;
	TR_SIG << "calling set_rotamers on " << rotamer_sets().nmoltenres() << " molten residues." << std::endl;
#endif

	// parent refers to AdditionalBackgroundNodesIG; G refers to the templated IG: PD or LMIG
	// call the parent initialize() method
	G::initialize( rot_sets );

	// save references to the rotamer_set because this will be needed later to figure out whether a new rotamer/state is hydrophobic or not
	for ( Size ii = 1; ii <= rotamer_sets().nmoltenres(); ++ii ) {
		get_surface_node( ii )->set_rotamers( rotamer_sets().rotamer_set_for_moltenresidue( ii ) );
	}

	// initializes some local variables that translate between the ig enumeration and resid
	set_num_residues_in_protein( pose().total_residue() );

	// unfortunately can't just use nmolten res here because moltenres includes residues which are being just repacked
	// and those shouldn't be considered first class nodes (would be very inefficient to do so!)
	//int num_designed = 0;
	//for (Size ii = 1; ii <= pose().total_residue(); ++ii) {
	//	if ( packer_task().being_designed(ii) ) {
	//		num_designed++;
	//	}
	//}
	//int nbackground = pose().total_residue() - num_designed;

	// 7/15/08 can just use nmolten res because that includes NATAA (packable) residues
	int nbackground = pose().total_residue() - rot_sets.nmoltenres();
	set_num_background_residues( nbackground );

	for (Size ii = 1; ii <= pose().total_residue(); ++ii) {

		// IS THE NEIGHBOR GRAPH even going to be populated at this point?  Hopefully, yes.  It seems that at the
		// start of pack_rotamers, the scorefunction is evaulated on the pose which should populate the graphs.

		// in our case, first class residues are residues that are designable or packable (see notes in function comment).
		// This decision leads to the possibility of designable/packable residues that are not surface-exposed but are
		// still considered FC nodes.  That's ok because in most of the functions, there are checks for being surface-exposed
		// and if those fail, the nodes immediately return 0.0.
		if ( packer_task().being_packed(ii) || packer_task().being_designed(ii) ) {
			// it's probably that being_packed() includes being_designed() so that the conditional needs only to say
			// being_packed(), but it'll short circuit the OR if being_packed() is true so it doesn't matter too much.
			continue;
		} else {
			set_residue_as_background_residue( ii );
		}
	}


	//core::scoring::TenANeighborGraph const & tenA_neighbor_graph( pose().energies().tenA_neighbor_graph() );

	// set the surface_exposed_ and is_below_buried_residue_no_hsasa booleans for all FC nodes, then BG nodes
	Size node_num_nbs = 0;
	for (int ii = 1; ii <= parent::get_num_nodes(); ++ii) {
		node_num_nbs = get_surface_node( ii )->num_neighbors_counting_self();

		// Note: we have to set two booleans here: surface_exposed_ and is_below_buried_residue_no_hsasa_cutoff_.
		// a Node may have more than the SURFACE_EXPOSED_CUTOFF number of neighbors, but still contribute hSASA to other nodes.

		if ( node_num_nbs <= (Size)SURFACE_EXPOSED_CUTOFF ) {
		// old way commented below
		// if ( tenA_neighbor_graph.get_node(ii)->num_neighbors_counting_self() <= SURFACE_EXPOSED_CUTOFF ) {
			get_surface_node( ii )->surface_exposed( true );
			get_surface_node( ii )->is_below_buried_residue_no_hsasa_cutoff( true );
#ifdef FILE_DEBUG
	TR_SIG << "initialize_surface_ig: setting surface_exposed and is_below_buried_residue_no_hsasa_cutoff bools to true for FC node " << ii << std::endl;
#endif
		} else if ( node_num_nbs > (Size)SURFACE_EXPOSED_CUTOFF && node_num_nbs <= (Size)BURIED_RESIDUE_NO_HSASA_CUTOFF ) {
			get_surface_node( ii )->is_below_buried_residue_no_hsasa_cutoff( true );
#ifdef FILE_DEBUG
	TR_SIG << "initialize_surface_ig: setting is_below_buried_residue_no_hsasa_cutoff bool to true for FC node " << ii << std::endl;
#endif
		}

	}

	for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii) {
		node_num_nbs = get_surface_bg_node( ii )->num_neighbors_counting_self();
		if ( node_num_nbs <= (Size)SURFACE_EXPOSED_CUTOFF ) {
		// if ( tenA_neighbor_graph.get_node(ii)->num_neighbors_counting_self() <= SURFACE_EXPOSED_CUTOFF ) {
			get_surface_bg_node( ii )->surface_exposed( true );
			get_surface_bg_node( ii )->is_below_buried_residue_no_hsasa_cutoff( true );
#ifdef FILE_DEBUG
	TR_SIG << "initialize_surface_ig: setting surface_exposed and is_below_buried_residue_no_hsasa_cutoff bools to true for BG node " << ii << std::endl;
#endif
		} else if ( node_num_nbs > (Size)SURFACE_EXPOSED_CUTOFF && node_num_nbs <= (Size)BURIED_RESIDUE_NO_HSASA_CUTOFF ) {
			get_surface_node( ii )->is_below_buried_residue_no_hsasa_cutoff( true );
#ifdef FILE_DEBUG
	TR_SIG << "initialize_surface_ig: setting is_below_buried_residue_no_hsasa_cutoff bool to true for BG node " << ii << std::endl;
#endif
		}
	}


#ifdef FILE_DEBUG
	TR_SIG << "initialize_surface_ig: DONE with initialization.\n---" << std::endl;
#endif

}


///
/// @begin SurfaceInteractionGraph::set_num_residues_in_protein
///
/// @brief
/// tells the graph how many residues there are total in the protein
///
/// @detailed
/// The graph maintains its own enumeration for the background residues;
/// but asks that anyone wanting to refer to them use their original resid.
/// The graph has to switch back and forth between enumeration schemes and
/// must know how many residues there are total to do that efficiently.
///
/// How does this method handle cofactors, ligands that get included in pose().total_residue()?
/// It definitely allocates space in the array for them, and those residues will apparently be made into first-class nodes.
/// But when it comes time to do the packing, they don't matter because they're not hydrophobic - they're not
/// even amino acids - so they don't get updated.
///
/// @param
/// num_res - [in] - the total number of residues
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::set_num_residues_in_protein( int num_res ) {

	num_total_residues_ = num_res;
	resid_2_bgenumeration_.resize( num_total_residues_ );
	for (int ii = 1; ii <= num_total_residues_; ++ii) {
		resid_2_bgenumeration_[ii] = 0;
	}
}

///
/// @begin SurfaceInteractionGraph::set_num_background_residues
///
/// @brief
/// tells the graph how many residues there are as part of the protein that are
/// not part of the combinatorial optimization process -- they are part of the
/// background
///
/// @detailed
/// The other half of switching between enumeration schemes for the background
/// residues is knowing how many background residues there are.
///
/// @param
/// num_background_residues - [in] - the number of background residues in the protein.
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::set_num_background_residues( int num_background_residues ) {

#ifdef FILE_DEBUG
	TR_SIG << "set_num_background_residues: setting num background residues to " << num_background_residues << std::endl;
#endif

	parent::set_num_background_nodes( num_background_residues );
	if ( parent::get_num_background_nodes() == 0 )
		return;

	bgenumeration_2_resid_.resize( parent::get_num_background_nodes() );
	for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii) {
		bgenumeration_2_resid_[ii] = 0;
	}
}


///
/// @begin SurfaceInteractionGraph::set_residue_as_background_residue
///
/// @brief
/// informs the graph that a particular residue is part of the background.
/// first class residues, in terms of surface, will be residues that are surface exposed and are varying. background
/// residues will be ones that are surface exposed but not varying since they will contribute to the surface score
/// of the first class residues and whose own surface score will change depending on the first class residues.
///
/// @param
/// residue - [in] - the residue identifier
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::set_residue_as_background_residue( int residue ) {

	assert( resid_2_bgenumeration_[ residue ] == 0 );

#ifdef FILE_DEBUG
	TR_SIG << "set_residue_as_background_residue: setting residue " << pose().residue( residue ).name3() << " " << residue << " as background node." << std::endl;
#endif

	++num_residues_assigned_as_background_;
	resid_2_bgenumeration_[ residue ] = num_residues_assigned_as_background_;
	bgenumeration_2_resid_[ num_residues_assigned_as_background_ ] = residue;

#ifdef FILE_DEBUG
	//TR_SIG << bgenumeration_2_resid_ << std::endl; // does not compile
#endif
}


///
/// @begin SurfaceInteractionGraph::bg_node_2_resid
///
/// @brief
/// Provides read access to the bg to resid array. Returns -1 if the index is not in bounds.
///
template < typename V, typename E, typename G >
int SurfaceInteractionGraph< V, E, G >::bg_node_2_resid( int node_index ) {

	if ( node_index > num_residues_assigned_as_background_ ) {
		utility_exit_with_message( "Out of bounds array index passed to bg_node_2_resid. Quitting." );
	}
	return bgenumeration_2_resid_[ node_index ];
}


///
/// @begin SurfaceInteractionGraph::prepare_for_simulated_annealing
///
/// @brief
/// Prepares the graph to begin simulated annealing.
///
/// @detailed
/// Invokes both base-class prepare_for_simulated_annealing subroutines:
/// InteractionGraphBase first, to prepare the SurfaceNodes and SurfaceEdges.
/// Then the AdditionalBackgroundNodesInteractionGraph, to prepare the
/// SurfaceBackgroundNodes, the SurfaceBackgroundEdges, and to do a little more
/// preparing of the SurfaceNodes.
///
/// @remarks
/// As it is written, it should only be executed once; that will change, however
/// if we make SurfaceInteractionGraph work with ligands (ligands that stay fixed
/// during any single sim annealing process, but that move between anealings.
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::prepare_for_simulated_annealing() {

	if ( prepared_for_simulated_annealing_ ) {

		// for the unit tests and potentially other use cases, prepare_for_simulated_annealing may be called more than once.
		// the graph will go into an unassigned state when blanket_assign_state_0 is called, but the starting hASA values
		// could be wrong if any consider_sub() calls were called after the first prep_for_simA() call but before the second
		// prep_for_simA() call.  But we don't need to go through and recalculate neighbors and creating Nodes/Edges and all
		// that.  All we really need to do is reset/re-initialize the starting hASA variables.  So, if we've already gone
		// through this method once, just tell all Nodes/BGNodes to reinit themselves.
		for (int ii = 1; ii <= parent::get_num_nodes(); ++ii)
			get_surface_node( ii )->init_hASA_variables();
		for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii)
			get_surface_bg_node( ii )->init_hASA_variables();

		return;
	}

	// since the count of se nbs gets determined in the prep_for_simA method, this call needs to occur before the
	// Nodes and BGNodes get their prep_for_simA() methods called.  Otherwise, there are no edges to traverse out to
	// when trying to determine how many se hp nbs a given Node has.
	detect_background_residue_and_first_class_residue_neighbors();

	// G::prepare() calls InteractionGraphBase::prepare_for_simulated_annealing() - LinmemIG implements one but it also
	// calls the IGBase method.  The IGBase method, in turn, calls prep_for_simA() on all the FC nodes.
	G::prepare_for_simulated_annealing();

	// parent::prepare() calls prep_for_simA() on all the BGNodes
	parent::prepare_for_simulated_annealing();

	prepared_for_simulated_annealing_ = true;

#ifdef FILE_DEBUG
	TR_SIG << "prepare_for_simulated_annealing: number edges in graph: " << parent::get_num_edges() << std::endl;
#endif

}


///
/// @begin SurfaceInteractionGraph::detect_background_residue_and_first_class_residue_neighbors
///
/// @brief
/// iterates across all pairs of first- and second-class nodes to determine which should be considered neighbors.
/// Adds a SurfaceBackgroundEdge between any pair that are.
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::detect_background_residue_and_first_class_residue_neighbors() {

	for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii) {

#ifdef FILE_DEBUG
		TR_SIG << "detect_bg_and_fc_residue_neighbors: checking for neighbors of background residue " << pose().residue( bgenumeration_2_resid_[ ii ] ).name3()
			<< " " << pose().residue( bgenumeration_2_resid_[ ii ] ).seqpos() << std::endl;
#endif

		for (int jj = 1; jj <= parent::get_num_nodes(); ++jj) {

			// ii: background node index, jj: first-class node index
			// the background nodes and first-class nodes should be set at this point due to the initalize method above.
			// if the background node is a "neighbor" of the first-class node, we need to add an edge between them.
			if ( get_surface_bg_node( ii )->detect_neighborship( get_surface_node( jj ) ) ) {
#ifdef FILE_DEBUG
				TR_SIG << "detect_bg_and_fc_residue_neighbors: --- adding FC/BG edge: fc node id:" << pose().residue(jj).seqpos()
						<< " / bg resid:" << bgenumeration_2_resid_[ ii ] << std::endl;
#endif
				parent::add_background_edge( jj, ii );
			}
		}
	}

#ifdef FILE_DEBUG
	TR_SIG << "DONE detecting background and first class neighbors.\n---" << std::endl;
#endif

}


///
/// @begin SurfaceInteractionGraph::blanket_assign_state_0
///
/// @brief
/// assigns state 0 -- the unassigned state -- to all (first class) vertices
/// in the graph
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::blanket_assign_state_0() {

#ifdef FILE_DEBUG
	TR_SIG << "blanket_assign_state_0() called" << std::endl;
#endif
	for (int ii = 1; ii <= parent::get_num_nodes(); ++ii) {
		get_surface_node( ii )->assign_zero_state();
	}
	update_internal_energy_totals_surface();
}


///
/// @begin SurfaceInteractionGraph::update_internal_energy_totals_surface
///
/// @brief
/// After every 2^10 commits, the graph traverses its nodes and edges and
/// re-tallies the total energy of the current state assignment.  This update
/// prevents the accumulation of numerical drift, increasing accuracy.
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::update_internal_energy_totals_surface() {

#ifdef FILE_DEBUG
	TR_SIG << "update_internal_energy_totals_surface" << std::endl;
	// leave the call to the print routine inside the #ifdef b/c this output is debugging related
	SurfaceNode< V, E, G >::print_surface_avoidance_stats();
#endif

	parent::update_internal_energy_totals();
	total_energy_current_state_assignment_ = parent::get_energy_PD_current_state_assignment();

	for (int ii = 1; ii <= parent::get_num_nodes(); ++ii) {
		Real fc_surface = get_surface_node( ii )->get_curr_state_surface_energy();
		total_energy_current_state_assignment_ += fc_surface;
	}

	for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii) {
		Real bg_surface = get_surface_bg_node( ii )->get_surface_score();
		total_energy_current_state_assignment_ += bg_surface;
	}

	num_commits_since_last_update_ = 0;
}


///
/// @begin SurfaceInteractionGraph::set_errorfull_deltaE_threshold
///
/// @brief
/// Allows the sim-annealer to specify a deltaE threshold above which, it is
/// no longer necessary to be very accurate.
///
/// @detailed
/// When the annealer asks the graph to consider a state substitution that
/// produces a large collision, the graph may approximate the surface deltaE
/// instead of recalculating it.  The deltaE returned by consider_substitution() will be
/// inaccurate, but if the annealer is unlikely to accept the substitution, then time can be saved.
/// The graph guarantees that if the annealer does commit that substitution that it will go back and
/// perform the surface computations and return an accurate total energy for the graph.
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::set_errorfull_deltaE_threshold( Real deltaE ) {

#ifdef FILE_DEBUG
	TR_SIG << "set_errorfull_deltaE_threshold: setting threshold to " << deltaE << std::endl;
	// leave this inside since it's debugging output
	SurfaceNode< V, E, G >::print_surface_avoidance_stats();
#endif

	SurfaceNode< V, E, G >::reset_surface_avoidance_stats();
	deltaE_threshold_for_avoiding_surface_calcs_ = deltaE;
}


///
/// @begin SurfaceInteractionGraph::consider_substitution
///
/// @brief
/// returns the change in energy induced by switching a particular node from its
/// currently assigned state to some alternate state.
///
/// @detailed
/// Also returns the sum of the two body energies for the node in its current
/// state; the sim-annealer accepts state substitutions at higher chance if the
/// original state was also at a poor energy.
///
/// @param
/// node_ind - [in] - the index of the (first class) node
/// new_state - [in] - the alternate state that the node should consider
/// delta_energy - [out] - the change in energy induced on the entire graph
///   by substituting a node's current state with the alternate.  This energy
///   may be inaccurate if it exceeds a threshold set by the sim-annealer.
/// prev_energy_for_node - [out] - the sum of the pair-wise decomposable portion
///   of the energy function for the node's currently assigned state
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::consider_substitution( int node_ind, int new_state, core::PackerEnergy & delta_energy, core::PackerEnergy & prev_energy_for_node ) {

	blanket_reset_alt_state_total_hASAs();

#ifdef FILE_DEBUG
	TR_SIG << "---" << std::endl;
	TR_SIG << "consider_substitution(): trying new state " << new_state << " ("
		<< get_surface_node( node_ind )->get_rotamer( new_state )->name() << ") on node/molten res " << node_ind << " (wt: "
		<< pose().residue( rotamer_sets().moltenres_2_resid( node_ind ) ).name3() << " " << rotamer_sets().moltenres_2_resid( node_ind ) << ") " << std::endl;
#endif

	node_considering_alt_state_ = node_ind;

	delta_energy = get_surface_node( node_ind )->project_deltaE_for_substitution_surface(
			new_state, prev_energy_for_node, deltaE_threshold_for_avoiding_surface_calcs_ );

#ifdef FILE_DEBUG
	TR_SIG << "consider_substitution: delta_energy: " << delta_energy << std::endl;
#endif
	total_energy_alternate_state_assignment_ = delta_energy + total_energy_current_state_assignment_;

}


///
/// @begin SurfaceInteractionGraph::blanket_reset_alt_state_total_hASAs
///
/// @brief
/// Iterates through all Nodes and BGNodes and resets their alt state total hASA to the current state total hASA.
///
/// @detailed
/// the Node consider() function is the main entry point from the SIG when the annealer considers a sub. need to make sure that the
/// node's alt_state hASA is the same as the current state before we start incrementing/decrementhing things or otherwise I wind
/// up with total hASAs that are wrong. alt state will be the same as current state on the first time through because of
/// prep for simA call, but not necessarily the case on next substitution coming from annealer.  this only sets the alt state
/// hASA on this particular (changing) node.  all of the other nodes which were called on to consider a substitution and
/// whose hASAs changed need to reset their alt state hASA also (or, otherwise, it appears as if the substitution was
/// committed).  More about this problem can be read in the comments for commit_considered_substitution_surface()
///
/// 02/23/2009 Updated this function and the functions it calls to do everything in terms of hASA instead of se hp counts.
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::blanket_reset_alt_state_total_hASAs() {

#ifdef FILE_DEBUG
	TR_SIG << "blanket_reset_alt_state_total_hASAs: calling reset_alt_state_total_hASA on all FC nodes." << std::endl;
#endif
	for (int ii = 1; ii <= parent::get_num_nodes(); ++ii) {
		get_surface_node( ii )->reset_alt_state_total_hASA();
	}

#ifdef FILE_DEBUG
	TR_SIG << "blanket_reset_alt_state_total_hASAs: calling reset_alt_state_total_hASA on all BG nodes." << std::endl;
#endif
	for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii) {
		get_surface_bg_node( ii )->reset_alt_state_total_hASA();
	}
}


///
/// @begin SurfaceInteractionGraph::commit_considered_substitution
///
/// @brief
/// Commits the substitution that the sim annealer had previously asked the graph to consider.  Returns the accurate total
/// energy for the graph.
///
template < typename V, typename E, typename G >
core::PackerEnergy SurfaceInteractionGraph< V, E, G >::commit_considered_substitution() {

#ifdef FILE_DEBUG
	TR_SIG << "commit_considered_substitution(): committing sub on node " << node_considering_alt_state_ << std::endl;
#endif
	core::PackerEnergy deltaE = get_surface_node( node_considering_alt_state_ )->commit_considered_substitution_surface();

	total_energy_current_state_assignment_ = total_energy_current_state_assignment_ + deltaE;

	node_considering_alt_state_ = -1;
	++num_commits_since_last_update_;
	if ( num_commits_since_last_update_ == COMMIT_LIMIT_BETWEEN_UPDATES ) {
		update_internal_energy_totals_surface();
	}

	return total_energy_current_state_assignment_;
}


///
/// @begin SurfaceInteractionGraph::set_network_state
///
/// @brief
/// Switch the state assignment of every first class node in the graph.
/// Useful, for instance, if you want to switch to the best network state
/// that you've found so far.  Like after simulated annealing is done with!
///
/// @remarks
/// Considered removing the use of FArrays here, but unfortunately it's what the Annealers use and I'm
/// not qualified to be refactoring the annealers.
///
/// @param
/// node_states - [in] - the array of states, one for each vertex in the graph
///
template < typename V, typename E, typename G >
core::PackerEnergy SurfaceInteractionGraph< V, E, G >::set_network_state( ObjexxFCL::FArray1_int & node_states) {

#ifdef FILE_DEBUG
	TR_SIG << "set_network_state() called with states: " << node_states << std::endl;
#endif

	blanket_reset_alt_state_total_hASAs();  // why is this needed here again?

	for (int ii = 1; ii <= parent::get_num_nodes(); ++ii) {
		get_surface_node( ii )->assign_state_surface( node_states(ii) );
	}
	update_internal_energy_totals_surface();

	return total_energy_current_state_assignment_;
}


///
/// @begin SurfaceInteractionGraph::set_state_for_node
///
/// @brief
/// Assign new_state to the first class node with index node_ind
///
/// @param
/// node_ind - [in] - the index of the (first class) node
/// new_state - [in] - the new state assigned to the node
///
template < typename V, typename E, typename G >
core::PackerEnergy SurfaceInteractionGraph< V, E, G >::set_state_for_node( int node_ind, int new_state ) {

	get_surface_node( node_ind )->assign_state_surface( new_state );

	update_internal_energy_totals_surface();
	return total_energy_current_state_assignment_;
}


///
/// @begin SurfaceInteractionGraph::get_energy_current_state_assignment
///
/// @brief
/// returns the energy of the entire graph under the current network state
/// assignment.
///
template < typename V, typename E, typename G >
core::PackerEnergy SurfaceInteractionGraph< V, E, G >::get_energy_current_state_assignment() {
	return total_energy_current_state_assignment_;
}


///
/// @begin SurfaceInteractionGraph::create_new_node
///
/// @brief
/// factory method pattern for instantiation of SurfaceNode objects, used by
/// InteractionGraphBase class.
///
/// @detailed
/// The IGBase class calls this method for all residues that are designable.
///
/// @param
/// node_index - [in] - the index of the (first class) node
/// num_states - [in] - the number of states for that node
///
template < typename V, typename E, typename G >
core::pack::interaction_graph::NodeBase*
SurfaceInteractionGraph< V, E, G >::create_new_node( int node_index, int num_states) {
#ifdef FILE_DEBUG
	TR_SIG << "create_new_node called with node_index " << node_index << " and num_states " << num_states << std::endl;
#endif
	return new SurfaceNode< V, E, G >( this, node_index, num_states );
}


///
/// @begin SurfaceInteractionGraph::create_new_edge
///
/// @brief
/// factory method pattern for instantiation of SurfaceEdge objects, used by
/// InteractionGraphBase class.
///
/// @param
/// index1 - [in] - the index of the lower-indexed (first class) node
/// index2 - [in] - the index of the higher-indexed (first class) node
///
template < typename V, typename E, typename G >
core::pack::interaction_graph::EdgeBase*
SurfaceInteractionGraph< V, E, G >::create_new_edge( int index1, int index2) {
#ifdef FILE_DEBUG
	TR_SIG << "create_new_edge() called for indices " << index1 << " and " << index2 << std::endl;
#endif
	return new SurfaceEdge< V, E, G >( this, index1, index2 );
}


///
/// @begin SurfaceInteractionGraph::create_background_node
///
/// @brief
/// factory method pattern for instantiation of SurfaceBackgroundResidueNode
/// objects, used by AdditionalBackgroundNodesInteractionGraph class.
///
/// @param
/// node_index - [in] - the index of the (second class) node
///
template < typename V, typename E, typename G >
BackgroundNode< V, E, G >* SurfaceInteractionGraph< V, E, G >::create_background_node( int node_index ) {
#ifdef FILE_DEBUG
	TR_SIG << "create_background_node() called for index " << node_index << std::endl;
#endif
	return new SurfaceBackgroundNode< V, E, G >( this, node_index );
}


///
/// @begin SurfaceInteractionGraph::create_background_edge
///
/// @brief
/// factory method pattern for instantiation of SurfaceBackgroundEdge
/// objects, used by AdditionalBackgroundNodesInteractionGraph class.
///
/// @param
/// fc_node_index - [in] - the index of the first class node
/// sc_node_index - [in] - the index of the second class node
///
/// @remarks
/// Returns a pointer to the created SurfaceBackgroundEdge
///
template < typename V, typename E, typename G >
BackgroundToFirstClassEdge< V, E, G >* SurfaceInteractionGraph< V, E, G >::create_background_edge( int fc_node_index, int bg_node_index ) {
	return new SurfaceBackgroundEdge< V, E, G >( this, fc_node_index, bg_node_index );
}


///
/// @begin SurfaceInteractionGraph::print_internal_energies_for_current_state_assignment
///
/// @brief
/// Prints out the one and two-body energies for every node and bgnode in the graph.
///
template < typename V, typename E, typename G >
void SurfaceInteractionGraph< V, E, G >::print_internal_energies_for_current_state_assignment() {

	// print out the one-body and surface energies for all first class nodes
	TR_SIG << "internal energies: " << std::endl;
	for (int ii = 1; ii <= parent::get_num_nodes(); ++ii) {
		Real one_body = get_surface_node( ii )->get_curr_state_one_body_energy();
		TR_SIG << "node " << ii << " 1b: " << one_body;

		Real surface = get_surface_node( ii )->get_curr_state_surface_energy();
		TR_SIG << ", surface = " << surface;

		if ( ii % 3 == 0) {
			TR_SIG << std::endl;
		}
	}

	TR_SIG << std::endl;

	// print out the surface energies for all background nodes
	for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii) {

		Real bg_surface = get_surface_bg_node( ii )->get_surface_score();
		TR_SIG << "bg res: " << bgenumeration_2_resid_[ ii ] << " surface: " << bg_surface << std::endl;
	}

	// print out the two-body energies for all edges between first-class nodes only?
	int count_edges = 0;
	for (std::list< core::pack::interaction_graph::EdgeBase*>::const_iterator iter = parent::get_edge_list_begin(); iter != parent::get_edge_list_end(); ++iter) {
		Real edge_energy = ((SurfaceEdge< V, E, G >*) (*iter))->get_current_two_body_energy();
		TR_SIG << "edge: " << edge_energy << " ";

		if ( count_edges % 5 == 0)
			TR_SIG << std::endl;
		++count_edges;
	}
}


///
/// @begin SurfaceInteractionGraph::get_edge_memory_usage
///
/// @brief
/// Should return a measurement of the memory used by the interaction graph
/// to store the rotamer pair energies.  Unimplemented.
///
template < typename V, typename E, typename G >
int SurfaceInteractionGraph< V, E, G >::get_edge_memory_usage() const {
	return 0;
}


///
/// @begin SurfaceInteractionGraph::count_static_memory
///
template < typename V, typename E, typename G >
unsigned int SurfaceInteractionGraph< V, E, G >::count_static_memory() const {
	return sizeof( SurfaceInteractionGraph< V, E, G > );
}


///
/// @begin SurfaceInteractionGraph::count_dynamic_memory
///
template < typename V, typename E, typename G >
unsigned int SurfaceInteractionGraph< V, E, G >::count_dynamic_memory() const {

	unsigned int total_memory = parent::count_dynamic_memory();
	total_memory += resid_2_bgenumeration_.size() * sizeof ( int );
	total_memory += bgenumeration_2_resid_.size() * sizeof ( int );

	return total_memory;
}


///
/// @begin SurfaceInteractionGraph::set_pose
///
template < typename V, typename E, typename G >
void
SurfaceInteractionGraph<V, E, G>::set_pose( pose::Pose const & pose ) {

#ifdef FILE_DEBUG
	TR_SIG << "set_pose() called: typeid() of base class G returned: " << typeid(G).name() << std::endl;
#endif

	// call the set_pose function in the LinMemIG class, because it, too, uses the Pose to do its thing
	if ( typeid(G) == typeid( pack::interaction_graph::LinearMemoryInteractionGraph ) ) {
		dynamic_cast<pack::interaction_graph::LinearMemoryInteractionGraph*>(this)->set_pose( pose );
	}

	pose_ = pose::PoseOP( new pose::Pose( pose ) );
}


///
/// @begin SurfaceInteractionGraph::set_packer_task
///
/// @brief
/// We need a copy of the packer task for 2 reasons: 1) to figure out what (if any) the weight on the surface score
/// is; and 2) to determine which residues are being packed and/or designed. We have to figure the packing options
/// because it determines whether a residue becomes a FirstClass (SurfaceNode) node or a background node.
///
template < typename V, typename E, typename G >
void
SurfaceInteractionGraph<V, E, G>::set_packer_task( task::PackerTask const & the_task ) {
	packer_task_ = the_task.clone();
}


///
/// @begin SurfaceInteractionGraph::set_rotamer_sets
///
template < typename V, typename E, typename G >
void
SurfaceInteractionGraph<V, E, G>::set_rotamer_sets( rotamer_set::RotamerSets const & rotsets ) {
	rotamer_sets_ = rotamer_set::RotamerSetsOP( new rotamer_set::RotamerSets( rotsets ) );
}


///
/// @begin SurfaceInteractionGraph::print
///
/// @brief
/// useful for debugging
///
template< typename V, typename E, typename G >
void
SurfaceInteractionGraph<V, E, G>::print() const {

	std::cout << "Surface Interaction Graph state: " << std::endl;
	std::cout << "nodes: " << std::endl;
	for (int jj = 1; jj <= parent::get_num_nodes(); ++jj) {
		get_surface_node( jj )->print();
	}

	std::cout << "bgnodes: " << std::endl;
	for (int ii = 1; ii <= parent::get_num_background_nodes(); ++ii) {
		get_surface_bg_node( ii )->print();
	}
}



// The below functions are only used for the unit tests. However, since most developers these days are running the unit
// tests using release mode, I can't #ifdef these functions (to leave them out of release mode builds) or the unit
// tests don't compile. So just compile them regardless of build mode.

///
/// @begin SurfaceNode::set_observed_sufficient_boolean_true
///
/// @brief
/// Sets the observed_sufficient_surface_E_to_predict_min_ to true. Only used by the unit tests.
///
template< typename V, typename E, typename G >
void SurfaceNode<V, E, G>::set_observed_sufficient_boolean_true() {
	observed_sufficient_surface_E_to_predict_min_ = true;

	for( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii)
		get_incident_surface_edge(ii)->set_max_surface_deltaE(); // defined in header, sets value to 0.1

	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii )
		get_edge_to_surface_bg_node(ii)->set_max_surface_deltaE(); // defined in header, sets value to 0.1

}

///
/// @begin SurfaceNode::get_hASA_for_node_and_nbs
///
/// @brief
/// Returns a vector of Reals containing the hASA of se hp nbs for this node in the first index, and then
/// the hASA for all neighboring nodes and bgnodes. Note: Order is not guaranteed here.
/// Only used by the unit tests.  Should not exist in release mode.
///
template< typename V, typename E, typename G >
std::vector< Real > SurfaceNode<V, E, G>::get_hASA_for_node_and_nbs() {

	std::vector<Real> hASA_vector;
	hASA_vector.push_back( curr_state_total_hASA_ );

	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		int other_nodes_index = get_incident_surface_edge(ii)->get_other_ind( parent::get_node_index() );
		int edge_index = ( other_nodes_index == get_incident_surface_edge(ii)->get_surface_node(0)->get_node_index() ? 0 : 1 );
		// "search" through the map to make sure this edge is in the map (if not found, find returns map.end() which is what I check for)
		if ( fc_neighbor_map.find( std::pair<int,int>( parent::get_node_index(), other_nodes_index) ) != fc_neighbor_map.end() ) {

			hASA_vector.push_back( get_incident_surface_edge(ii)->get_surface_node( edge_index )->get_current_hASA() );
		}
	}

	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		int other_nodes_index = get_edge_to_surface_bg_node(ii)->get_other_ind( this );
		// "search" through the map to make sure this edge is in the map (if not found, find returns map.end() which is what I check for)
		if ( bg_neighbor_map.find( std::pair<int,int>( parent::get_node_index(), other_nodes_index) ) != bg_neighbor_map.end() ) {

			hASA_vector.push_back( get_edge_to_surface_bg_node(ii)->get_surface_bg_node()->get_current_hASA() );
		}
	}

	return hASA_vector;
}

/// @begin SurfaceNode::get_current_hASA
///
/// @brief
/// Returns current hASA. Only used by the unit tests.  Should not exist in release mode.
///
template< typename V, typename E, typename G >
Real SurfaceNode<V, E, G>::get_current_hASA() { return curr_state_total_hASA_; }

/// @begin SurfaceBackgroundNode::get_current_count
///
/// @brief
/// Returns current hASA. Only used by the unit tests.  Should not exist in release mode.
///
template< typename V, typename E, typename G >
Real SurfaceBackgroundNode<V, E, G>::get_current_hASA() { return curr_state_total_hASA_; }



///
/// @begin SurfaceNode::get_alt_state_hASA_for_node_and_nbs
///
/// @brief
/// Returns a vector of Reals containing the hASA of se hp nbs for this node in the first index, and then
/// the hASA for all neighboring nodes and bgnodes. Note: Order is not guaranteed here.
/// Only used by the unit tests.  Should not exist in release mode.
///
template< typename V, typename E, typename G >
std::vector< Real > SurfaceNode<V, E, G>::get_alt_state_hASA_for_node_and_nbs() {

	std::vector< Real > alt_state_hASA_vector;
	alt_state_hASA_vector.push_back( alt_state_total_hASA_ );

	for ( int ii = 1; ii <= parent::get_num_incident_edges(); ++ii ) {
		int other_nodes_index = get_incident_surface_edge(ii)->get_other_ind( parent::get_node_index() );
		int edge_index = ( other_nodes_index == get_incident_surface_edge(ii)->get_surface_node(0)->get_node_index() ? 0 : 1 );
		// "search" through the map to make sure this edge is in the map (if not found, find returns map.end() which is what I check for)
		if ( fc_neighbor_map.find( std::pair<int,int>( parent::get_node_index(), other_nodes_index) ) != fc_neighbor_map.end() ) {

			alt_state_hASA_vector.push_back( get_incident_surface_edge(ii)->get_surface_node( edge_index )->get_alt_hASA() );
		}
	}

	for ( int ii = 1; ii <= parent::get_num_edges_to_background_nodes(); ++ii ) {
		int other_nodes_index = get_edge_to_surface_bg_node(ii)->get_other_ind( this );
		// "search" through the map to make sure this edge is in the map (if not found, find returns map.end() which is what I check for)
		if ( bg_neighbor_map.find( std::pair<int,int>( parent::get_node_index(), other_nodes_index) ) != bg_neighbor_map.end() ) {

			alt_state_hASA_vector.push_back( get_edge_to_surface_bg_node(ii)->get_surface_bg_node()->get_alt_hASA() );
		}
	}

	return alt_state_hASA_vector;
}


/// @begin SurfaceNode::get_alt_hASA
///
/// @brief
/// Returns current hASA. Only used by the unit tests.  Should not exist in release mode.
///
template< typename V, typename E, typename G >
Real SurfaceNode<V, E, G>::get_alt_hASA() { return alt_state_total_hASA_; }

/// @begin SurfaceBackgroundNode::get_alt_hASA
///
/// @brief
/// Returns current hASA. Only used by the unit tests.  Should not exist in release mode.
///
template< typename V, typename E, typename G >
Real SurfaceBackgroundNode<V, E, G>::get_alt_hASA() { return alt_state_total_hASA_; }


///
/// @begin SurfaceInteractionGraph::get_network_state
///
/// @brief
/// Returns the state on each FCNode, but not necessarily in pose resid order. Only used by the unit tests.
///
template< typename V, typename E, typename G >
std::vector< int > SurfaceInteractionGraph<V, E, G>::get_network_state() const {

	std::vector< int > networkstate;
	for (int jj = 1; jj <= parent::get_num_nodes(); ++jj) {
		networkstate.push_back( get_surface_node(jj)->get_current_state() );
	}
	return networkstate;
}

///
/// @begin SurfaceInteractionGraph::get_hASA_for_node_and_nbs
///
/// @brief
/// Sets the observed_sufficient_surface_E_to_predict_min_ to true. Only used by the unit tests.
///
template< typename V, typename E, typename G >
std::vector< Real > SurfaceInteractionGraph<V, E, G>::get_hASA_for_node_and_nbs( int index ) {
	return get_surface_node( index )->get_hASA_for_node_and_nbs();
}


///
/// @begin SurfaceInteractionGraph::get_alt_state_hASA_for_node_and_nbs
///
/// @brief
/// Sets the observed_sufficient_surface_E_to_predict_min_ to true. Only used by the unit tests.
///
template< typename V, typename E, typename G >
std::vector< Real > SurfaceInteractionGraph<V, E, G>::get_alt_state_hASA_for_node_and_nbs( int index ) {
	return get_surface_node( index )->get_alt_state_hASA_for_node_and_nbs();
}


///
/// @begin SurfaceInteractionGraph::set_observed_sufficient_boolean_true
///
/// @brief
/// Sets the observed_sufficient_surface_E_to_predict_min_ to true. Only used by the unit tests.
///
template< typename V, typename E, typename G >
void SurfaceInteractionGraph<V, E, G>::set_observed_sufficient_boolean_true() {
	for (int jj = 1; jj <= parent::get_num_nodes(); ++jj) {
		get_surface_node(jj)->set_observed_sufficient_boolean_true();
	}
}


}
}
} //end all namespaces


#endif
