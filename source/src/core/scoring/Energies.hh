// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Energies.hh
/// @brief  Energies class
/// @author Phil Bradley


#ifndef INCLUDED_core_scoring_Energies_hh
#define INCLUDED_core_scoring_Energies_hh

// Package Headers
#include <core/scoring/Energies.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/scoring/MinimizationGraph.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>

#ifdef WIN32 //VC++ needs full class declaration
#include <core/scoring/ContextGraph.hh> // WIN32 INCLUDE
#include <core/scoring/LREnergyContainer.hh> // WIN32 INCLUDE
#endif

// Project headers
#include <core/graph/Graph.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <core/conformation/PointGraph.fwd.hh>

#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>

// Utility Headers
#include <utility/py/PyAssert.hh>
#include <utility/pointer/ReferenceCount.hh>

// STL Headers
#include <map>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.hh>

//Auto Headers
#include <core/scoring/EnergyGraph.fwd.hh>
namespace core {
namespace scoring {


/// A cached energies object

/**
@li  Stores total, residue, and residue-pair energies, as well as
residue neighbor information.

@li  Meant to replace fullatom_energies:: namespace.

@li  Stores residue neighbor information as well as cached residue
pair energies in O(N) space using a graph.

@li  Also stores a DomainMap object which is used during scoring
to know which rsd pairs have changed relative orientation and
which residues have changed internally

@note  We distinguish between two kinds of per-residue (1D) energy:
onebody residue energies and twobody residue energies. Onebody
residue energies are things like dunbrack, intrares, Paa, which
depend only on the state of the residue in question. Twobody
residue energies (like the residue atr energy) are summations
of twobody interactions involving a single residue. The onebody
energies can be reused at positions whose internal phi/psi/chi
conformation hasn't changed. Twobody residue energies, on the other
hand must be invalidated if the structure has changed at all.
**/

class Energies : public utility::pointer::ReferenceCount
{

public:
	typedef scoring::EnergyGraph EnergyGraph;
	typedef scoring::NeighborList NeighborList;
	typedef scoring::ScoreType ScoreType;
	typedef scoring::ScoreFunctionInfo ScoreFunctionInfo;
	typedef basic::datacache::BasicDataCache BasicDataCache;

	//typedef scoring::NeighborEnergies NeighborEnergies;

	typedef id::AtomID AtomID;
	typedef id::AtomID_Mask AtomID_Mask;
	typedef kinematics::DomainMap DomainMap;

	typedef graph::Graph Graph;

	typedef scoring::EnergyMap EnergyMap;

public:
	/// ctor -- ensure correct initial state
	Energies();

	// Explicit copy ctor to avoid #include .hh's
	Energies( Energies const & src );

	// Explicit assignmnet operator to avoid #include .hh's
	virtual Energies const & operator = ( Energies const & rhs );

	/// @details determine whether my type is the same as another Conformation's
	virtual
	bool
	same_type_as_me( Energies const & other, bool recurse = true ) const;

	/// dtor
	virtual
	~Energies();

	virtual
	EnergiesOP
	clone() const;


	/// @brief Pose must claim its Energies object; this should happen once,
	/// at the time the Pose is allocated.  Future copying of the Energies object
	/// will not change ownership.  The purpose of ownership is to allow lazy context-graph
	/// creation.  When context graphs are requested that have not been created, they
	/// must be created and their edges must be updated to represent the current conformation
	/// of the pose.
	void set_owner( pose::Pose * owner );

	//////////////////////////
	/// energy access
	//////////////////////////

	/// @brief Returns the total score
	///
	/// example(s):
	///     pose.energies().total_energy()
	/// See also:
	///     Energies
	///     Energies.residue_total_energy
	///     Energies.residue_total_energies
	///     Pose
	///     ScoreFunction
	///     ScoreFunction.show
	///     create_score_function
	Real
	total_energy() const
	{
		return total_energy_;
	}

	Real &
	total_energy()
	{
		return total_energy_;
	}

	/// @brief Returns the total_energies EnergyMap after first computing the
	/// component energies if they are not up-to-date
	EnergyMap const &
	total_energies() const;

	/// @brief Returns a non-const reference to the total_energies EnergyMap
	/// so that external sources may append additional information to the Energies
	/// object.  This is primarily useful for outputting score data with structures
	/// when those terms are not part of the energy function.
	/// This function will update the component energies if they are not up-to-date.
	EnergyMap &
	total_energies();

	/// @brief Read access to the components of the one-body energies.
	EnergyMap const &
	onebody_energies( int const seqpos ) const
	{
		return onebody_energies_[ seqpos ];
	}

	/// @brief Write access to the components of the one-body energies.
	/// This access is intended only for the ScoreFunction.
	EnergyMap &
	onebody_energies( int const seqpos )
	{
		return onebody_energies_[ seqpos ];
	}

	/// @brief Read access to the components of the "finalized" energies;
	/// These will include any score component calculated in the finalize
	/// phase of score function evaluation.  These energies are copied
	/// between Energies objects, and are not evaluated during the component-
	/// energy update.
	EnergyMap const &
	finalized_energies() const {
		return finalized_energies_;
	}

	/// @brief Write access to the components of the "finalized" energies.
	/// This access is intended only for the ScoreFunction.
	EnergyMap &
	finalized_energies() {
		return finalized_energies_;
	}

	/// @brief Returns the unweighted total_energies EnergyMap for
	/// Residue  <seqpos>
	/// @note: only evaluated when requested (lazy!), inaccessible during
	/// minimization, EnergyMap is an EMapVector
	///
	/// example(s):
	///     r3 = pose.energies().residue_total_energies(3)
	///     r3[fa_sol]
	/// See also:
	///     Energies
	///     Energies.residue_total_energy
	///     Pose
	///     ScoreFunction
	///     ScoreFunction.show
	///     create_score_function
	///     EMapVector
	EnergyMap const &
	residue_total_energies( int const seqpos ) const
	{
		debug_assert( !use_nblist() && energies_updated() );
		//  PyAssert( (!use_nblist()) && (energies_updated()), "Energies::residue_total_energies(): the Energies object isn't ready! Has it been scored?" );
		PyAssert( (seqpos>0) && (seqpos<=int(size())), "Energies::residue_total_energies( int const seqpos ): variable seqpos is out of range!" );
		if ( ! residue_total_energies_uptodate_ ) accumulate_residue_total_energies();
		return residue_total_energies_[ seqpos ];
	}

	/// @brief Returns the weighted total energy of residue  <seqpos>
	///
	/// example(s):
	///     pose.energies().residue_total_energy(3)
	/// See also:
	///     Energies
	///     Energies.residue_total_energies
	///     Pose
	///     ScoreFunction
	///     ScoreFunction.show
	///     create_score_function
	Real
	residue_total_energy( int const seqpos ) const
	{
		debug_assert( !use_nblist() && energies_updated() );
		//  PyAssert( (!use_nblist()) && (energies_updated()), "Energies::residue_total_energy(): the Energies object isn't ready! Has it been scored?" );
		PyAssert( (seqpos>0) && (seqpos<=int(size())), "Energies::residue_total_energy( int const seqpos ): variable seqpos is out of range!" );
		if ( ! residue_total_energy_uptodate_ ) accumulate_residue_total_energy();
		return residue_total_energy_[ seqpos ];
	}

	/// @brief Read access to the EnergyGraph.
	EnergyGraph const &
	energy_graph() const;

	/// @brief Write access to the EnergyGraph.
	EnergyGraph &
	energy_graph();

	/// @brief get the graph encoding # neighbors within 10 Angstroms
	/// If the graph has not been requested up until this point, then it will
	/// be instantiated and filled.  If the pose has changed size since the last
	/// score function evaluation (or if the pose has never been scored) this
	/// function will exit.
	TenANeighborGraph const &
	tenA_neighbor_graph() const;

	/// @brief Write access to the graph encoding # neighbors within 10 Angstroms
	/// If the graph has not been requested up until this point, then it will
	/// be instantiated and filled.  If the pose has changed size since the last
	/// score function evaluation (or if the pose has never been scored) this
	/// function will exit.
	TenANeighborGraph &
	tenA_neighbor_graph();

	/// @brief get the graph encoding # neighbors within 12 Angstroms
	scoring::TwelveANeighborGraph const &
	twelveA_neighbor_graph() const;

	scoring::TwelveANeighborGraph &
	twelveA_neighbor_graph();


	scoring::ContextGraphCOP
	context_graph( scoring::ContextGraphType type ) const;

	scoring::ContextGraphOP
	context_graph( scoring::ContextGraphType type );

	/// @brief Allows non-scorefunction components of Rosetta to impose requirements on
	/// the context graphs that this object maintains.
	void
	require_context_graph( scoring::ContextGraphType ) const;

	/// @brief Returns an EnergyMap of the ScoreFunction weights from the last
	/// scoring
	///
	/// example(s):
	///     we = pose.energies().weights()
	///     we[fa_atr]
	/// See also:
	///     Energies
	///     Energies.residue_total_energies
	///     Energies.residue_total_energy
	///     Pose
	///     ScoreFunction
	///     ScoreFunction.show
	///     create_score_function
	///     EMapVector
	EnergyMap weights() const;

	// tex 10/31/08 - added this method for access by the silent-file
	// classes, so that we can initialize weights from silent-files.
	/// @brief Setter for the weights in this Energies object.
	void weights( EnergyMap new_weights );

	/// @brief update the residue neighbors
	void
	update_residue_neighbors(
		DomainMap const & domain_map_in,
		pose::Pose const & pose
	);


	bool
	residue_neighbors_updated() const
	{
		return ( graph_state_ == GOOD );
	}

	/// @brief Returns true if the score is up-to-date
	///
	/// example(s):
	///     pose.energies().energies_updated()
	/// See also:
	///     Energies
	///     Energies.residue_total_energy
	///     Energies.residue_total_energies
	///     Pose
	///     ScoreFunction
	///     ScoreFunction.show
	///     create_score_function
	bool
	energies_updated() const
	{
		return ( energy_state_ == GOOD );
	}


	/// @brief check if rsd has changed internal conformation, necessitating,  recomputation of 1d energies like dun,intra,prob,etc
	//
	bool
	res_moved( int const seqpos ) const;

	void
	reset_res_moved( int const seqpos );

	/// @brief for debugging -- forget all stored energies, does not change size
	void
	clear_energies();


	/// @brief Returns the number of held residue energies
	///
	/// example(s):
	///     r3 = pose.energies().residue_total_energies(3)
	///     r3[fa_sol]
	/// See also:
	///     Energies
	///     Energies.residue_total_energies
	///     Energies.residue_total_energy
	///     Pose
	///     ScoreFunction
	///     ScoreFunction.show
	///     create_score_function
	inline
	Size
	size() const
	{
		return size_;
	}

	/// are we in the midst of a scoring calculation?
	bool
	scoring() const
	{
		return scoring_;
	}


	void
	show( std::ostream & out ) const;


	void
	show( std::ostream & out, Size res ) const;

	//wrapper function of energies.show() for Pyrosetta
	void
	show() const {show(std::cout);};

	/// @brief Shows the energy information of residue  <seqpos>
	/// @note: wrapper function of energies.show(Size) for Pyrosetta
	///
	/// example(s):
	///     pose.energies().show(3)
	/// See also:
	///     Energies
	///     Energies.residue_total_energies
	///     Energies.residue_total_energy
	///     Pose
	///     ScoreFunction
	///     ScoreFunction.show
	///     create_score_function
	void
	show(Size res) const {show(std::cout, res);};

	void
	show_totals( std::ostream & out ) const;

	void
	show_total_headers( std::ostream & out ) const;

	friend std::ostream & operator<<(std::ostream& out, const Energies& e );

	/// @brief called (eg by pose) to notify us of a change to the structure
	/**
	Triggers clearing of the total energies and the twobody rsd energies
	PHIL -- should also mark the neighbor links as invalid somehow...
	Called by pose when someone tries to access us, if the Conformation
	indicates that the structure has moved since the last score evaluation

	const b/c called inside const access methods
	**/
	void
	structure_has_moved( Size const nres ) const;

	/// @brief  Notification of the start of a scoring calculation.

	void
	scoring_begin(
		scoring::ScoreFunction const & sfxn,
		pose::Pose const & pose // for the nbr calculation
	);

	/// @brief signal from the scorefxn that scoring is over
	void
	scoring_end(
		scoring::ScoreFunction const & scorefxn
	);

	/// @brief get scorefxn info
	scoring::ScoreFunctionInfo const &
	get_scorefxn_info() const
	{
		return *scorefxn_info_;
	}

	/// @brief kill everything (that nobody forgot about)
	void
	clear();


	bool
	use_nblist() const
	{
		return use_nblist_;
	}

	bool
	use_nblist_auto_update() const {
		return use_nblist_auto_update_;
	}

	MinimizationGraphOP
	minimization_graph();

	MinimizationGraphCOP
	minimization_graph() const;

	void
	set_minimization_graph( MinimizationGraphOP );


	scoring::NeighborList const &
	nblist( EnergiesCacheableDataType::Enum const & type ) const;


	void
	set_nblist(
		EnergiesCacheableDataType::Enum const & type,
		scoring::NeighborListOP nblist_in
	);


	void
	set_use_nblist(
		pose::Pose const & pose,
		DomainMap const & domain_map_in,
		bool const use_nblist_auto_update
	);


	void
	reset_nblist();

	/// @brief BasicDataCache indexed by enum in core/scoring/EnergiesCacheableDataType.hh
	BasicDataCache const &
	data() const
	{
		return data_cache_;
	}

	/// @brief BasicDataCache indexed by enum in core/scoring/EnergiesCacheableDataType.hh
	BasicDataCache &
	data()
	{
		return data_cache_;
	}

	/// @brief instructs Pose whether the domain map info in the Conformation object should be discarded
	bool discard_conformation_domain_map() const;

	/// @brief Return the color assigned to a particular residue (index = pos) as held in the
	/// domain map describing how this residue has moved with respect to the other residues in the
	/// pose.
	///
	/// CAUTION new behavior: domain_map may not return 0 for residues that have undergone internal
	/// degree of freedom changes since the last scoring.  Ask the res_moved() method for that information
	int
	domain_map( int const pos ) const
	{
		require_scoring();
		return domain_map_(pos);
	}

	/// @brief Read access for the domain map.
	DomainMap const &
	domain_map() const
	{
		require_scoring();
		return domain_map_;
	}

	void
	set_long_range_container( methods::LongRangeEnergyType, LREnergyContainerOP );

	LREnergyContainerOP
	nonconst_long_range_container( methods::LongRangeEnergyType );

	LREnergyContainerCOP
	long_range_container( methods::LongRangeEnergyType ) const;

protected:

	/// @brief get access to the point graph. For derived classes
	conformation::PointGraphOP
	point_graph();

	/// @brief Write access to the EnergyGraph.
	EnergyGraph &
	energy_graph_no_state_check();

	utility::vector1< ContextGraphOP > &
	context_graphs() const;

	utility::vector1< bool > &
	required_context_graphs() const;

	Real
	max_context_neighbor_cutoff() const;

	void
	set_max_context_neighbor_cutoff( Real val ) const;

	pose::Pose *
	owner() const;

	/// @brief Read access for the domain map. There is only one
	/// difference to the public interface to the domain_map_:
	/// we don't reequire that scoring is performed. This function
	/// is used during minimization by derived classes of the minimizer
	int
	domain_map_during_minimization( int const pos ) const
	{
		return domain_map_(pos);
	}


	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// private methods
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////

private:

	/// @brief should push this into energygraph once it stabilizes
	enum EnergyState {
		BAD, MOD, GOOD
	};

	/// @brief Internal method to resize data that is dependent on the number
	/// of residues in the Pose.
	void
	set_size( Size const size_in );


	/// @brief will die if not in the middle of an energy evaluation
	inline
	void
	require_scoring() const;

	/// @brief Lazy component evaluation for residue pair energies;
	void
	update_component_energies() const;

	/// @brief Lazy update of residue total component energies; only compute them when they are requested.
	/// Forces an update of the energy_graph_.
	void
	accumulate_residue_total_energies() const;

	/// @brief Lazy update of the per-residue total weighted energy.  Only computed when requested.
	/// Does not force an update of the energy_graph_
	void
	accumulate_residue_total_energy() const;


	/// @brief sum the residue energies to this type to get a total energy
	void
	accumulate_residue_energies( ScoreType const & type ) const;

	/// @brief update the context graphs and the energy graph according to a new
	/// score function type
	void
	set_scorefxn_info( scoring::ScoreFunctionInfoOP info );

	/// @brief Save state information from the domain map and wipe dirty energies
	void
	internalize_new_domain_map();

	/// @brief Delete edges from energy_graph and tenA_neighbor_graph according to
	/// the movement history described by the domain map.
	void prepare_neighbor_graphs();

	//// @brief Delete edges for a particular graph using the domain map data
	void delete_graph_edges_using_domain_map( graph::Graph & g);

	/// @brief Reset the "already computed" status for pairs of residues represented
	/// in a particular long-range energy container using the domain map.
	void update_domainmap_for_lr_energy_container( LREnergyContainerOP lrec );

	/// @brief Detect the new set of neighbors given the structure of the Pose
	/// (find_neighbors()) and add new edges to the neighbor graphs so that the
	/// status of the neighbor graphs is current wrt the current conformation.
	virtual
	void update_neighbor_links( pose::Pose const & pose );

	/// @brief Create a point graph representing the xyz coordinates of the "neighbor
	/// atoms."  (Each residue_type indicates one of its atoms to be used for
	/// neighbor detection -- e.g. CB for most amino acids, CA for GLY.)  The point
	/// graph will then be used in a call to find_neighbors to add upper-edges to
	/// neighboring residues.
	virtual
	void fill_point_graph( pose::Pose const & pose, conformation::PointGraphOP pg ) const;

	/// @brief During Energies copy-ctor and assignment operator, copy over the
	/// neighbor lists objects (clone them).
	void copy_nblists( Energies const & other );

	/// @brief During Energies copy-ctor and assignment operator, copy over the
	/// context graphs and the historical information about who required them
	/// (a score term, e.g. the HBondEnergy, or some other function, e.g. pack_rotamers)
	void copy_context_graphs( Energies const & other );

	//// @brief During Energies copy-ctor and assignment operator, copy over
	/// the long-range energy containers from the source Energies object.
	void copy_lr_energy_containers( Energies const & other );

	/// @brief Internal method that handles the bookkeeping for intializing a
	/// new context graph, either requested by a scoring term (external = false)
	/// or some other function (external = true ).
	virtual
	void
	require_context_graph_( scoring::ContextGraphType type, bool external ) const;


private:

	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////

	/// our internal nres
	/**
	Used for dimensioning new Energy1D's, eg
	and for checking the validity of the NeighborEnergies by comparing
	against their size
	**/
	Size size_;

	/// In order to do lazy context-graph creation, the Energies object must be able to access
	/// coordinates stored in the pose. This now means an Energies object cannot
	/// live independently of a Pose!
	// Should eventually be PoseAP once Pose's c'tor is protected and an instantiator
	// method is implemented (w/ copy via clone).
	pose::Pose * owner_;

	mutable EnergyGraphOP energy_graph_;

	/// @brief The collection of context graphs used by the context-dependent energy components
	/// which the Energies object is responsible for maintaining.  Context graphs
	/// are allocated as requested, (e.g. in a call to tenA_neighbor_graph() ), if
	/// they have not yet been allocated.  If a portion of the code requires a context
	/// graph (e.g the packer), it may call require_context_graph( int cgtype ), which
	/// will create the context graph immediately
	mutable utility::vector1< ContextGraphOP > context_graphs_;
	/// those required by non-score function entities (e.g. the packer).
	mutable utility::vector1< bool > externally_required_context_graphs_;
	/// OR of the sfxn required context graphs and the externally required ones.
	mutable utility::vector1< bool > required_context_graphs_;
	/// The maximum neighbor cutoff for all required context graphs;
	mutable Real max_context_neighbor_cutoff_;

	utility::vector1< LREnergyContainerOP > long_range_energy_containers_;

	/// atom-atom neighborlists
	std::map< EnergiesCacheableDataType::Enum, scoring::NeighborListOP > nblist_;
	bool use_nblist_;
	bool use_nblist_auto_update_;
	MinimizationGraphOP minimization_graph_;

	/// cached onebody energies -- expensive -- shortly to be replaced by an FArray2D which can be smarly indexed into
	/// to access only the active one-body energy terms.
	utility::vector1< EnergyMap > onebody_energies_;

	/// cached energy component totals -- summed twobody and one body interactions for each residue
	mutable bool residue_total_energies_uptodate_;
	mutable utility::vector1< EnergyMap > residue_total_energies_;
	/// cached energy totals -- only the weighted sum, no components.
	mutable bool residue_total_energy_uptodate_;
	mutable utility::vector1< Real > residue_total_energy_; // The weighted sum of all interaction energies for each residue

	/// cached total energies
	mutable EnergyMap total_energies_;
	mutable Real total_energy_;

	/// Energies computed during the finalize() stage of scoring.
	EnergyMap finalized_energies_;

	/// info about last score evaluation
	scoring::ScoreFunctionInfoOP scorefxn_info_;

	/// ScoreFunction weights from the last scoring call
	EnergyMap scorefxn_weights_;

	/// Domain map, stores information about the rigid-bodies whose internal conformation is unchanged since
	/// the last score calc'n
	/**
	If domain_map_(i) == 0 then residue i has changed internal conformation. -- no longer true
	If domain_map_(i) > 0 and domain_map_(i) == domain_map_(j),  then residues
	i and j are unchanged wrt one another
	**/
	DomainMap domain_map_;

	/// are we within a scoring evaluation?
	bool scoring_;

	mutable EnergyState energy_state_;
	mutable EnergyState graph_state_;

	/// @brief BasicDataCache indexed by enum in core/scoring/EnergiesCacheableDataType.hh
	/// @warning DataCache must always be initialized with the number of cacheable
	///  data types -- see the last enum entry.
	BasicDataCache data_cache_;
	//std::map< CacheableDataType, CacheableDataOP > cached_data_;

	/// Keep this guy between score function evaluations to
	/// avoid the expense of recreating it each time; this data does not
	/// need to be copied in either the copy-ctor or the assignment operator
	/// Its purpose is solely to improve performance and the data is used
	/// only inside the neighbor calculation function call.
	conformation::PointGraphOP point_graph_;
};


inline
void
Energies::require_scoring() const
{
	if ( !scoring_ ) {
		utility_exit_with_message(
			"Energies:: operation only permitted during scoring." );
	}
}


} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_Energies_HH
