// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/Energies.cc
/// @brief  Energies class to store cached energies and track the residue
/// neighbor relationships
/// @author Phil Bradley
/// @author Andrew Leaver-Fay

// Unit Headers
#include <core/scoring/Energies.hh>

#include <basic/Tracer.hh>

// Package Headers
#include <core/scoring/ContextGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/ContextGraphFactory.hh>
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/scoring/NeighborList.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>

#include <core/pose/symmetry/util.hh>


// Project Headers
#include <core/conformation/PointGraph.hh>
#include <core/conformation/find_neighbors.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/id/AtomID.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// Numeric headers
#include <numeric/numeric.functions.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>

//Auto Headers
#include <utility/vector1.hh>
//Auto Headers
#include <core/conformation/PointGraphData.hh>
#include <core/graph/ArrayPool.hh>
#include <core/graph/UpperEdgeGraph.hh>
#include <core/scoring/EnergyGraph.hh>


static THREAD_LOCAL basic::Tracer tr( "core.scoring.Energies" );

using namespace ObjexxFCL::format;

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>
#include <utility/serialization/ObjexxFCL/FArray1D.srlz.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {

Energies::Energies()
: utility::pointer::ReferenceCount(),
	size_(0),
	owner_( 0 ),
	energy_graph_( EnergyGraphOP( new EnergyGraph ) ),
	context_graphs_( scoring::num_context_graph_types, 0 ),
	externally_required_context_graphs_( scoring::num_context_graph_types, false ),
	required_context_graphs_( scoring::num_context_graph_types, false ),
	max_context_neighbor_cutoff_( 0.0 ),
	long_range_energy_containers_( scoring::methods::n_long_range_types, 0 ),
	use_nblist_(false),
	use_nblist_auto_update_(false),
	minimization_graph_( /* 0 */ ),
	residue_total_energies_uptodate_( false ),
	residue_total_energy_uptodate_( false ),
	total_energy_( 0.0 ),
	scorefxn_info_( scoring::ScoreFunctionInfoOP( new ScoreFunctionInfo ) ),
	scorefxn_weights_(),
	scoring_(false),
	energy_state_( BAD ),
	graph_state_( BAD ),
	data_cache_( EnergiesCacheableDataType::num_cacheable_data_types ),
	point_graph_( /* 0 */ )
{}


/// copy ctor -- deep copy
Energies::Energies( Energies const & other )
: utility::pointer::ReferenceCount(),
	size_( other.size_ ),
	owner_( 0 ),
	energy_graph_( EnergyGraphOP( new EnergyGraph( *other.energy_graph_ ) ) ),
	context_graphs_( scoring::num_context_graph_types, 0 ),
	externally_required_context_graphs_( other.externally_required_context_graphs_ ),
	required_context_graphs_( other.required_context_graphs_ ),
	max_context_neighbor_cutoff_( other.max_context_neighbor_cutoff_ ),
	long_range_energy_containers_( scoring::methods::n_long_range_types, 0 ),
	use_nblist_( other.use_nblist_ ),
	use_nblist_auto_update_( other.use_nblist_auto_update_ ),
	minimization_graph_( other.minimization_graph_ ? new MinimizationGraph( * other.minimization_graph_ ) : 0 ),
	onebody_energies_( other.onebody_energies_ ),
	residue_total_energies_uptodate_( false ),
	residue_total_energies_( size_ ),
	residue_total_energy_uptodate_( false ),
	residue_total_energy_( size_, 0.0 ),
	total_energies_( other.total_energies_ ),
	total_energy_( other.total_energy_ ),
	finalized_energies_( other.finalized_energies_ ),
	scorefxn_info_( scoring::ScoreFunctionInfoOP( new ScoreFunctionInfo( *(other.scorefxn_info_) ) )),
	scorefxn_weights_( other.scorefxn_weights_ ),
	domain_map_( other.domain_map_ ),
	scoring_( other.scoring_ ),
	energy_state_( other.energy_state_ ),
	graph_state_( other.graph_state_ ),
	data_cache_( other.data_cache_ ),
	point_graph_( /* 0 */ )
{
	copy_nblists( other );
	copy_context_graphs( other );
	copy_lr_energy_containers( other );
}

/// assignment operator -- deep copy
Energies &
Energies::operator = ( Energies const & rhs )
{
	if ( this == &rhs ) return *this;

	size_ =  rhs.size_;
	(*energy_graph_) = (*rhs.energy_graph_);
	context_graphs_.resize( scoring::num_context_graph_types);
	use_nblist_ =  rhs.use_nblist_;
	use_nblist_auto_update_ = rhs.use_nblist_auto_update_;
	minimization_graph_ = MinimizationGraphOP( rhs.minimization_graph_ ? new MinimizationGraph( * rhs.minimization_graph_ ) : 0 );
	onebody_energies_ =  rhs.onebody_energies_;
	residue_total_energies_uptodate_ = false;
	//std::fill( residue_total_energies_.begin(), residue_total_energies_.end(), EnergyMap() ); // unncessary
	residue_total_energy_uptodate_ = rhs.residue_total_energy_uptodate_;
	residue_total_energy_ = rhs.residue_total_energy_;
	total_energies_ = rhs.total_energies_;
	total_energy_ = rhs.total_energy_;
	finalized_energies_ = rhs.finalized_energies_;
	if ( *scorefxn_info_ == *rhs.scorefxn_info_ ) {
		// noop; sfxn info matches, so no need to duplicate the score function.
	} else {
		scorefxn_info_ = scoring::ScoreFunctionInfoOP( new ScoreFunctionInfo( *(rhs.scorefxn_info_)) );
	}
	scorefxn_weights_ = rhs.scorefxn_weights_;
	domain_map_ =  rhs.domain_map_;
	scoring_ =  rhs.scoring_;
	energy_state_ =  rhs.energy_state_;
	graph_state_ =  rhs.graph_state_;
	data_cache_ = rhs.data_cache_;

	copy_nblists( rhs );
	copy_context_graphs( rhs );
	copy_lr_energy_containers( rhs );

	/// NOTE: point_graph_ is intentionally not copied here ////

	return *this;
}

/// @details If recurse is true, then this is the first call to same_type_as_me;
// determine if the other object is also an Energies object.  If recurse is false
// then the other object is also of type Energies, so return true.
bool
Energies::same_type_as_me( Energies const & other, bool recurse /* = true */ ) const
{
	if ( recurse ) {
		return other.same_type_as_me( *this, false );
	} else {
		return true;
	}
}


Energies::~Energies()
{
	//std::cout << "energies dstor" << std::endl;
}

/// @details make a copy of this Energies( allocate actual memory for it )
EnergiesOP
Energies::clone() const
{
	return EnergiesOP( new Energies( *this ) );
}


void
Energies::set_owner( pose::Pose * owner ) {
	if ( owner == 0 ) {
		owner_ = 0;
	} else if ( owner_ != 0 ) {
		utility_exit_with_message( "Attempted to set the owner twice, once with " + utility::to_string( owner_ ) + " and now with " + utility::to_string( owner ) );
	} else {
		owner_ = owner;
	}
}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// accessors
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


EnergyMap const &
Energies::total_energies() const
{
	return total_energies_;
}

EnergyMap &
Energies::total_energies()
{
	return total_energies_;
}

///////////////////////////////////////////////////////////////////////////////
/// @brief get the energy graph (const)
///
/// @details assert that the graph state reflects the current set of neighbors; if the
/// structure has moved, the neighbors may have changed.
/// update_residue_neighbors should have been called first.  If nb_list_ is true
/// then update_residue_neighbors will set the graph state to GOOD without
/// updating the neighbors.
Energies::EnergyGraph const &
Energies::energy_graph() const
{
	debug_assert( graph_state_ == GOOD );
	return *energy_graph_;
}

/// @brief get the energy graph - see comments for const version of this method
Energies::EnergyGraph &
Energies::energy_graph()
{
	debug_assert( graph_state_ == GOOD );
	return *energy_graph_;
}

/// @details IA: Note derived classes may need access to graph before it has been initalized!
/// In particular a derived class might overload update_residue_neighbors.
Energies::EnergyGraph &
Energies::energy_graph_no_state_check()
{
	return *energy_graph_;
}

/// @details convenience function -- this function violates the idea that the
/// energies object doesn't know the kind of context information its
/// maintaining.  There is no reason to add a method like this for future
/// context graphs.
scoring::TenANeighborGraph const &
Energies::tenA_neighbor_graph() const
{
	using namespace scoring;

	debug_assert( graph_state_ == GOOD );

	if ( ! context_graphs_[ ten_A_neighbor_graph ] ) {
		require_context_graph_( ten_A_neighbor_graph, true );
	}

	debug_assert( dynamic_cast< TenANeighborGraph const * > (context_graphs_[ ten_A_neighbor_graph ].get()) );

	return static_cast< TenANeighborGraph const & > ( *context_graphs_[ ten_A_neighbor_graph ] );
}

/// @brief get the graph encoding # neighbors within 10 Angstroms -- see comments for const version of this method
scoring::TenANeighborGraph &
Energies::tenA_neighbor_graph()
{
	using namespace scoring;

	debug_assert( graph_state_ == GOOD );

	if ( ! context_graphs_[ ten_A_neighbor_graph ] ) {
		require_context_graph_( ten_A_neighbor_graph, true );
	}

	debug_assert( dynamic_cast< TenANeighborGraph const * > (context_graphs_[ ten_A_neighbor_graph ].get()) );

	return static_cast< TenANeighborGraph & > ( *context_graphs_[ ten_A_neighbor_graph ] );
}

/// @details convenience function -- this function violates the idea that the
/// energies object doesn't know the kind of context information its
/// maintaining.  There is no reason to add a method like this for future
/// context graphs.
scoring::TwelveANeighborGraph const &
Energies::twelveA_neighbor_graph() const
{
	using namespace scoring;
	debug_assert( graph_state_ == GOOD );

	if ( ! context_graphs_[ twelve_A_neighbor_graph ] ) {
		require_context_graph_( twelve_A_neighbor_graph, true );
	}

	debug_assert( dynamic_cast< TwelveANeighborGraph const * > (context_graphs_[ twelve_A_neighbor_graph ].get()) );

	return static_cast< TwelveANeighborGraph const & > ( *context_graphs_[ twelve_A_neighbor_graph ] );
}

/// @brief get the graph encoding # CB (or Gly CA) neighbors within 12 Angstroms
/// of the side-chain "centroid" atom for each residue.
scoring::TwelveANeighborGraph &
Energies::twelveA_neighbor_graph()
{
	using namespace scoring;

	debug_assert( graph_state_ == GOOD );

	if ( ! context_graphs_[ twelve_A_neighbor_graph ] ) {
		require_context_graph_( twelve_A_neighbor_graph, true );
	}

	debug_assert( dynamic_cast< TwelveANeighborGraph const * > (context_graphs_[ twelve_A_neighbor_graph ].get()) );

	return static_cast< TwelveANeighborGraph & > ( *context_graphs_[ twelve_A_neighbor_graph ] );
}


/// @brief non-const access to a particular context graph
scoring::ContextGraphOP
Energies::context_graph( scoring::ContextGraphType type )
{
	debug_assert( graph_state_ == GOOD );
	if ( context_graphs_[ type ] == 0 ) require_context_graph_( type, true );
	return context_graphs_[ type ];
}

/// @brief const access to a particular context graph
scoring::ContextGraphCOP
Energies::context_graph( scoring::ContextGraphType type ) const
{
	debug_assert( graph_state_ == GOOD );
	if ( context_graphs_[ type ] == 0 ) require_context_graph_( type, true );
	return context_graphs_[ type ];
}

void
Energies::require_context_graph( scoring::ContextGraphType type ) const
{
	require_context_graph_( type, true );
}


/// @details May not be called during scoring, to discourage access by energy methods.
EnergyMap
Energies::weights() const
{
	if ( scoring_ ) {
		utility_exit_with_message(
			"Energies:: operation NOT permitted during scoring." );
	}

	return scorefxn_weights_;
}


void
Energies::weights( EnergyMap new_weights ) {
	scorefxn_weights_ = new_weights;
}

/// @brief get access to the point graph. For derived classes
conformation::PointGraphOP
Energies::point_graph()
{
	return point_graph_;
}

utility::vector1< ContextGraphOP > &
Energies::context_graphs() const
{
	return context_graphs_;
}

utility::vector1< bool > &
Energies::required_context_graphs() const
{
	return required_context_graphs_;
}

Real
Energies::max_context_neighbor_cutoff() const
{
	return max_context_neighbor_cutoff_;
}

void
Energies::set_max_context_neighbor_cutoff( Real val ) const
{
	max_context_neighbor_cutoff_ = val;
}

pose::Pose *
Energies::owner() const
{
	return owner_;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// private
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
void
Energies::clear_energies()
{
	domain_map_ = 0;
	total_energies_.clear();
	onebody_energies_.clear();
	residue_total_energies_.clear();
	residue_total_energies_uptodate_ = false;
	residue_total_energy_.clear();
	residue_total_energy_uptodate_ = false;
	if ( size_ ) {
		onebody_energies_.resize( size_ );
		residue_total_energies_.resize( size_ );
		residue_total_energy_.resize( size_ );
		std::fill( residue_total_energy_.begin(), residue_total_energy_.end(), (Real) 0.0 );
	}
	energy_graph_->drop_all_edges();
	for ( uint ii = 1; ii <= context_graphs_.size(); ++ii ) {
		if ( context_graphs_[ ii ] ) context_graphs_[ ii ]->drop_all_edges();
	}
	// is this really what we want to do here?
	for ( Size ii = 1; ii <= long_range_energy_containers_.size(); ++ii ) {
		long_range_energy_containers_[ ii ] = 0;
	}
	graph_state_ = BAD;
	energy_state_ = BAD;
}


///////////////////////////////////////////////////////////////////////////////
void
Energies::prepare_neighbor_graphs()
{
	// The graph state should either be GOOD (no work requried) or MOD (correct, mod the domain map)
	debug_assert( graph_state_ != BAD );
	delete_graph_edges_using_domain_map( *energy_graph_ );

	// Later, the absence of edges from the energy_graph_ is taken as a signal
	// that all residues have moved with respect to each other, and therefore
	// all found neighbor pairs will have edges added in the context graphs.
	// The context graphs must be empty to avoid an edge duplication event.
	if ( energy_graph_->num_edges() != 0 ) {
		for ( uint ii = 1; ii <= context_graphs_.size(); ++ii ) {
			if ( context_graphs_[ ii ] ) delete_graph_edges_using_domain_map( *context_graphs_[ ii ] );
		}
	} else {
		for ( uint ii = 1; ii <= context_graphs_.size(); ++ii ) {
			if ( context_graphs_[ ii ] ) context_graphs_[ ii ]->drop_all_edges();
		}
	}

	// Tell LR graphs about the new domain-map
	for ( Size ii = 1; ii <= long_range_energy_containers_.size(); ++ii ) {
		if ( long_range_energy_containers_[ ii ] && ! long_range_energy_containers_[ ii ]->empty() ) {
			update_domainmap_for_lr_energy_container( long_range_energy_containers_[ ii ] );
		}
	}
}


///////////////////////////////////////////////////////////////////////////////
/// @brief removes all edges in the graph that represent connections that have changed;
///
/// @details O(N) as it iterates across the edges that already exist, instead over over all pairs
/// that moved with respect to each other (which would be O(N^2)) and deleting any
/// obsolete edges
void Energies::delete_graph_edges_using_domain_map( Graph & g )
{
	using namespace graph;
	for ( Graph::EdgeListIter iter = g.edge_list_begin(),
			iter_end = g.edge_list_end(); iter != iter_end; /* no increment statement*/ ) {
		Graph::EdgeListIter iter_next = iter;
		++iter_next;

		int const n1( (*iter)->get_first_node_ind() ), n2( (*iter)->get_second_node_ind() );
		if ( domain_map_( n1 ) == 0 || domain_map_( n2 ) == 0 || domain_map_( n1 ) != domain_map_( n2 ) ) {
			g.delete_edge(*iter); //drop the edge from the graph.
		}
		iter = iter_next;
	}
}

void Energies::update_domainmap_for_lr_energy_container(
	LREnergyContainerOP lrec
)
{
	// Potentially O(N^2) operation...
	for ( Size ii = 1; ii <= size_; ++ii ) {
		int iimap = domain_map_( ii );
		for ( ResidueNeighborIteratorOP
				rni = lrec->upper_neighbor_iterator_begin( ii ),
				rniend = lrec->upper_neighbor_iterator_end( ii );
				(*rni) != (*rniend); ++(*rni) ) {
			if ( iimap == 0 || iimap != domain_map_( rni->upper_neighbor_id() ) ) {
				rni->mark_energy_uncomputed();
			}
		}
	}
}

MinimizationGraphOP
Energies::minimization_graph()
{
	return minimization_graph_;
}

MinimizationGraphCOP
Energies::minimization_graph() const
{
	return minimization_graph_;
}

void
Energies::set_minimization_graph( MinimizationGraphOP mingraph )
{
	minimization_graph_ = mingraph;
}


scoring::NeighborList const &
Energies::nblist( EnergiesCacheableDataType::Enum const & type ) const
{
	debug_assert( use_nblist_ && nblist_.find( type ) != nblist_.end() );
	return *( nblist_.find( type )->second );
}


void
Energies::set_nblist(
	EnergiesCacheableDataType::Enum const & type,
	scoring::NeighborListOP nblist_in
)
{
	debug_assert( use_nblist_ && nblist_.find( type ) == nblist_.end() );
	nblist_[ type ] = nblist_in;
}

///////////////////////////////////////////////////////////////////////////////
void
Energies::set_use_nblist(
	pose::Pose const & ,
	DomainMap const & domain_map_in,
	bool const use_nblist_auto_update
)
{
	debug_assert( !use_nblist_ && !scoring_ );
	debug_assert( energy_state_ == GOOD ); // we should have just scored

	use_nblist_ = true;
	use_nblist_auto_update_ = use_nblist_auto_update;

	domain_map_ = domain_map_in;
	internalize_new_domain_map();

	// setup the energy graph and context graphs
	// this deletes pair-moved edges
	// APL MOD 6.29.2010 -- no longer delete edges holding CD scores -- prepare_neighbor_graphs();

	// this adds appropriate pair-moved edges
	/// APL MOD 6.29.2010 -- do not try to add new edges either -- update_neighbor_links( pose );

	// so at this point, the neighbor links are correct and complete
	// but the neighbor-links for moving pairs are empty (no cached energies)
	// (here moving pairs are defined wrt the passed-in set of moving
	// dofs)
	//
	// this is what we want, since inside atom_pair_energy we want to
	// recover cached energies for non-moving rsd pairs while
	// calculating interactions for moving rsd pairs using the atom-atom nblist
	//
}


/// @details  Show the per-residue and total energies, only for scores with a non-zero weight in the last
/// energy evaluation

void
Energies::show( std::ostream & out ) const
{
	debug_assert( energy_state_ == GOOD );

	if ( ! residue_total_energies_uptodate_ ) accumulate_residue_total_energies();

	// show the header (names of scoretypes)
	out << "E       "; // "E" plus spaces to account for "(i)" plus the four character alignment of the residue number
	for ( EnergyMap::const_iterator it=total_energies_.begin(),
			it_end = total_energies_.end(); it != it_end; ++it ) {
		ScoreType const scoretype = ScoreType( it - total_energies_.begin() + 1 ); // hacky
		if ( scorefxn_weights_[ scoretype ] != 0 ) {
			out << A(14, name_from_score_type(scoretype));
		}
	}
	out << std::endl;

	// show the onebody energies
	for ( Size i=1; i<= size(); ++i ) {
		out << "E(i)" << I(4,i);
		EnergyMap const & emap( residue_total_energies_[i] );
		for ( EnergyMap::const_iterator it = emap.begin(), it_end = emap.end(); it != it_end; ++it ) {
			ScoreType const scoretype = ScoreType( it - emap.begin() + 1 ); // hacky
			if ( scorefxn_weights_[ scoretype ] != 0 ) {
				out << F(14,2,*it);
			}
		}
		out << std::endl;
	}

	// show the total energies
	for ( EnergyMap::const_iterator it=total_energies_.begin(),
			it_end = total_energies_.end(); it != it_end; ++it ) {
		ScoreType const scoretype = ScoreType( it - total_energies_.begin() + 1 ); // hacky
		if ( scorefxn_weights_[ scoretype ] != 0 ) {
			//out << "total_energy " << scoretype << ' ' << F(12,3,*it) << std::endl;
			out << "total_energy" << A(14, name_from_score_type(scoretype)) << F(12,3,*it) << std::endl;
		}
	}
}

/// @details  Show the residue energies, only for scores with a non-zero weight in the last
/// energy evaluation
void
Energies::show( std::ostream & out, Size res ) const
{
	debug_assert( energy_state_ == GOOD );

	if ( ! residue_total_energies_uptodate_ ) accumulate_residue_total_energies();
	// show the header (names of scoretypes)
	out << "E       "; // "E" plus spaces to account for "(i)" plus the four character alignment of the residue number
	for ( EnergyMap::const_iterator it=total_energies_.begin(),
			it_end = total_energies_.end(); it != it_end; ++it ) {
		ScoreType const scoretype = ScoreType( it - total_energies_.begin() + 1 ); // hacky
		if ( scorefxn_weights_[ scoretype ] != 0 ) {
			out << A(14, name_from_score_type(scoretype));
		}
	}
	out << std::endl;

	out << "E(i)" << I(4,res);
	EnergyMap const & emap( residue_total_energies_[res] );
	for ( EnergyMap::const_iterator it = emap.begin(), it_end = emap.end(); it != it_end; ++it ) {
		ScoreType const scoretype = ScoreType( it - emap.begin() + 1 ); // hacky
		if ( scorefxn_weights_[ scoretype ] != 0 ) {
			out << F(14,2,*it);
		}
	}
	out << std::endl;
}

std::ostream & operator<<(std::ostream & out, const Energies& e )
{
	if ( e.size() == 0 ) {
		out << "The pose must be scored first to generate an energy table." << std::endl;
	} else {
		// per-residue energies
		out << A( 4, "res" );
		//EnergyMap const & temp_emap( e.residue_total_energies( 1 ) );
		EnergyMap const & temp_emap( e.get_scorefxn_info().scores_present() );
		for ( EnergyMap::const_iterator it = temp_emap.begin(),
				it_end = temp_emap.end(); it != it_end; ++it ) {
			if ( *it ) {
				out << A( 15, name_from_score_type( ScoreType( it - temp_emap.begin() + 1 ) ) );
			}
		}
		out << std::endl;

		for ( Size i=1; i <= e.size(); ++i ) {
			out << I( 4, i );
			EnergyMap const & emap( e.residue_total_energies( i ) );
			for ( EnergyMap::const_iterator it = emap.begin(),
					it_end = emap.end(); it != it_end; ++it ) {
				if ( temp_emap[ ScoreType( it - emap.begin() + 1) ] ) out << F( 15, 3, *it );
			}
			out << std::endl;
		}

		// total energies
		out << A( 4, "tot" );
		for ( EnergyMap::const_iterator it = e.total_energies().begin(),
				it_end = e.total_energies().end(); it != it_end; ++it ) {
			if ( temp_emap[ ScoreType( it - e.total_energies().begin() + 1 ) ] ) {
				out << F(15,3,*it);
			}
			//out << "total_energy ";
			//  << ScoreType( it - e.total_energies().begin() + 1 ) << ' '
		}
		out << std::endl;
	}

	return out;
}


void
Energies::structure_has_moved( Size const nres )
const {
	if ( scoring_ ) {
		utility_exit_with_message("No structure mods allowed during scoring!");
	}

	residue_total_energies_uptodate_ = false;
	residue_total_energy_uptodate_ = false;

	if ( energy_state_ == GOOD ) energy_state_ = MOD;
	if ( graph_state_ == GOOD ) graph_state_ = MOD;

	total_energy_ = 0.0;

	if ( nres != size_ ) {
		energy_state_ = BAD;
		graph_state_ = BAD;
	}
}

void
Energies::show_total_headers( std::ostream & out ) const
{
	// total energies
	EnergyMap const & temp_emap( get_scorefxn_info().scores_present() );
	for ( EnergyMap::const_iterator it = total_energies().begin(),
			it_end = total_energies().end(); it != it_end; ++it ) {
		if ( temp_emap[ ScoreType( it - total_energies().begin() + 1 ) ] ) {
			// out << A( 15, ScoreType( it - total_energies().begin() + 1 ) );
			out << ScoreType( it - total_energies().begin() + 1 ) << ' ';
		}
	}
}

void
Energies::show_totals( std::ostream & out ) const
{
	// total energies

	EnergyMap const & temp_emap( get_scorefxn_info().scores_present() );
	for ( EnergyMap::const_iterator it = total_energies().begin(),
			it_end = total_energies().end(); it != it_end; ++it ) {
		if ( temp_emap[ ScoreType( it - total_energies().begin() + 1 ) ] ) {
			// out << ScoreType( it - total_energies().begin() + 1 ) << ' '
			//    << F(15,3,*it);
			out << F( 15,3, *it );
		}

	}
}

void
Energies::clear()
{
	clear_energies();
	domain_map_.clear();
	size_ = 0;

	set_scorefxn_info( scoring::ScoreFunctionInfoOP( new ScoreFunctionInfo ) );
	scorefxn_weights_.zero();
	return;
}

///////////////////////////////////////////////////////////////////////////////
void
Energies::reset_nblist()
{
	debug_assert( use_nblist_ && !scoring_ );
	use_nblist_ = false;
	nblist_.clear();

	// trigger complete recalculation of nbr links
	// since we have not been updating these during minimization
	//
	//fpd  Is setting graphstate to BAD necessary? It seems to cause recalculation
	//fpd  of every energy term, even when a very small portion of the graph moved during
	//fpd  minimization.  Changing this to MOD
	graph_state_ = MOD;

	/// APL destroy the minimization graph
	minimization_graph_.reset();

	structure_has_moved( size_ );
}

void
Energies::reset_res_moved( int const seqpos )
{
	energy_graph_->get_energy_node( seqpos )->moved( false );
}

///////////////////////////////////////////////////////////////////////////////
void
Energies::set_size( Size const new_size )
{
	if ( size_ == new_size ) {
		debug_assert( domain_map_.size() == size_ &&
			Size(energy_graph_->num_nodes()) == size_ );
		return;
	}
	size_ = new_size;
	energy_graph_->set_num_nodes( size_ );
	for ( uint ii = 1; ii <= context_graphs_.size(); ++ii ) {
		if ( context_graphs_[ ii ] ) context_graphs_[ ii ]->set_num_nodes( size_ );
	}
	for ( uint ii = 1; ii <= long_range_energy_containers_.size(); ++ii ) {
		if ( long_range_energy_containers_[ ii ] ) long_range_energy_containers_[ ii ]->set_num_nodes( size_ );
	}
	domain_map_.dimension( size_ );
	onebody_energies_.clear();
	onebody_energies_.resize( size_ );
	total_energies_.clear();
	residue_total_energies_uptodate_ = false;
	residue_total_energies_.clear();
	residue_total_energies_.resize( size_ );
	residue_total_energy_uptodate_ = false;
	residue_total_energy_.clear();
	residue_total_energy_.resize( size_ );
	std::fill( residue_total_energy_.begin(), residue_total_energy_.end(), (Real) 0.0 );
}

///////////////////////////////////////////////////////////////////////////////
void
Energies::update_residue_neighbors(
	DomainMap const & domain_map_in,
	pose::Pose const & pose
)
{

	if ( graph_state_ == GOOD ) return;

	if ( size_ != pose.total_residue() ) {
		set_size( pose.total_residue() );
		graph_state_ = BAD;
	}

	if ( graph_state_ == BAD ) {
		graph_state_ = MOD;
		domain_map_ = 0;
	} else {
		domain_map_ = domain_map_in;
	}

	internalize_new_domain_map();

	// no updates to the neighbors during minimization
	// the neighbors should have been calculated inside set_use_nblist
	if ( use_nblist_ ) {
		debug_assert( size_ == pose.total_residue() && graph_state_ != BAD );
		graph_state_ = GOOD;  // pretend
		return;
	}

	prepare_neighbor_graphs();
	update_neighbor_links( pose );
	graph_state_ = GOOD;

}


/// @details because the Energies object does not update the EnergyGraph during minimization, the Pose should
/// not direct the conformation object to discard its movement data.  Movement data should be discarded
/// immediately after a call to update_neighbors where use_nblist_ has been false.
bool
Energies::discard_conformation_domain_map() const
{
	return ! use_nblist_;
}

void
Energies::set_long_range_container( methods::LongRangeEnergyType lrtype, LREnergyContainerOP lrec)
{
	long_range_energy_containers_[ lrtype ] = lrec;
}

LREnergyContainerOP
Energies::nonconst_long_range_container( methods::LongRangeEnergyType lrtype )
{
	return long_range_energy_containers_[ lrtype ];
}

LREnergyContainerCOP
Energies::long_range_container( methods::LongRangeEnergyType lrtype ) const
{
	return long_range_energy_containers_[ lrtype ];
}


void
Energies::internalize_new_domain_map()
{

	for ( uint ii = 1, ii_end = size_; ii <= ii_end; ++ii ) {

		if ( domain_map_(ii) == 0 ) {
			//energy_graph_->get_energy_node( ii )->moved( true );
			energy_graph_->get_energy_node( ii )->moved( true );

			/// moved from scoring_begin()
			// onebody residue energies are still valid unless bb/chi changed
			// for that sequence position
			//std::cout << "Energies::scoring_begin() res_moved: " << i << std::endl;
			onebody_energies_[ii].clear();
		}
	}
}

void
Energies::accumulate_residue_total_energies() const
{
	if ( residue_total_energies_uptodate_ ) return;

	if ( residue_total_energies_.size() != size_ ) {
		residue_total_energies_.resize( size_ );
	}

	if ( energy_state_ != GOOD ) {
		for ( Size ii = 1; ii <= residue_total_energies_.size(); ++ii ) {
			residue_total_energies_[ ii ].zero();
		}
		/// Consider the residue_total_energies up to date in the sense
		/// that they reflect the current total_energies emap; not in the
		/// sense that they reflect the preserved energies held in the
		/// energy graph.
		residue_total_energies_uptodate_ = true;

		return;
	}

	ScoreTypes twobody_score_types;
	ScoreTypes all_score_types;
	for ( Size ii = 1; ii <= n_score_types; ++ii ) {
		/// ODDITY: The former code accumulated all energies, even those with weights of 0;
		/// I'll duplicate that behavior here, though it does not seem to be the correct meaning
		/// of residue total energies.
		/// uncomment the if block to instead accumulate only those residue energies with a non-zero weight.
		//if ( scorefxn_weights_[ (ScoreType) ii ] != 0.0 ) {
		if ( ii <= n_shortranged_2b_score_types ) twobody_score_types.push_back( (ScoreType) ii );
		all_score_types.push_back( (ScoreType) ii );
		//}
	}

	// debug
	debug_assert( !use_nblist() && energies_updated() );

	// start with the one-body energies
	for ( Size i=1, i_end = residue_total_energies_.size(); i<= i_end; ++i ) {
		residue_total_energies_[i] = onebody_energies_[i];
	}

	for ( Size i=1, i_end = residue_total_energies_.size(); i<= i_end; ++i ) {
		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph_->get_node(i)->const_upper_edge_list_begin(),
				irue = energy_graph_->get_node(i)->const_upper_edge_list_end();
				iru != irue; ++iru ) {

			EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );

			Size const j( edge.get_second_node_ind() );

			// accumulate energies
			for ( Size ii = 1; ii <= twobody_score_types.size(); ++ii ) {
				residue_total_energies_[ i ][ twobody_score_types[ ii ]] += 0.5 * edge[ twobody_score_types[ ii ]];
				residue_total_energies_[ j ][ twobody_score_types[ ii ]] += 0.5 * edge[ twobody_score_types[ ii ]];
			}

		} // nbrs of i
	} // i=1,nres

	for ( Size ii = 1; ii <= long_range_energy_containers_.size(); ++ii ) {
		LREnergyContainerCOP lrec = long_range_energy_containers_[ ii ];
		if ( lrec == 0 ) continue;
		if ( lrec->empty() ) continue;

		// Potentially O(N^2) operation...
		for ( Size i = 1; i <= residue_total_energies_.size(); ++i ) {
			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( i ),
					rniend = lrec->const_upper_neighbor_iterator_end( i );
					(*rni) != (*rniend); ++(*rni) ) {
				Size j = rni->upper_neighbor_id();
				EnergyMap emap;
				rni->retrieve_energy( emap );

				residue_total_energies_[ i ].accumulate( emap, all_score_types, 0.5 );
				residue_total_energies_[ j ].accumulate( emap, all_score_types, 0.5 );
			}
		}
	}

	for ( Size i=1, i_end = residue_total_energies_.size(); i<= i_end; ++i ) {
		residue_total_energies_[i][ total_score ] = residue_total_energies_[i].dot( scorefxn_weights_ );
	}

	residue_total_energies_uptodate_ = true;
}

void
Energies::accumulate_residue_total_energy() const
{
	if ( residue_total_energy_uptodate_ ) return;

	ScoreTypes twobody_score_types;
	for ( Size ii = 1; ii <= n_shortranged_2b_score_types; ++ii ) {
		if ( scorefxn_weights_[ (ScoreType) ii ] != 0.0 ) {
			twobody_score_types.push_back( (ScoreType) ii );
		}
	}

	for ( Size ii = 1; ii <= size_; ++ii ) {
		residue_total_energy_[ ii ] = scorefxn_weights_.dot( onebody_energies_[ ii ] );
	}

	for ( Size ii = 1; ii <= size_; ++ii ) {
		//conformation::Residue const & resl( pose.residue( i ) );

		for ( graph::Graph::EdgeListConstIter
				iru  = energy_graph_->get_node(ii)->const_upper_edge_list_begin(),
				irue = energy_graph_->get_node(ii)->const_upper_edge_list_end();
				iru != irue; ++iru ) {
			EnergyEdge const & edge( static_cast< EnergyEdge const & > (**iru) );
			Size const jj( edge.get_second_node_ind() );

			Real half = 0.5 * edge.dot( scorefxn_weights_ );

			residue_total_energy_[ ii ] += half;
			residue_total_energy_[ jj ] += half;
		} // nbr jj of ii
	} // ii=1,nres

	for ( Size ii = 1; ii <= long_range_energy_containers_.size(); ++ii ) {
		LREnergyContainerCOP lrec = long_range_energy_containers_[ ii ];
		if ( lrec == 0 ) continue;
		if ( lrec->empty() ) continue;

		// Potentially O(N^2) operation...
		for ( Size ii = 1; ii <= size_; ++ii ) {
			for ( ResidueNeighborConstIteratorOP
					rni = lrec->const_upper_neighbor_iterator_begin( ii ),
					rniend = lrec->const_upper_neighbor_iterator_end( ii );
					(*rni) != (*rniend); ++(*rni) ) {
				Size jj = rni->upper_neighbor_id();
				EnergyMap emap;
				rni->retrieve_energy( emap );
				Real half = 0.5 * scorefxn_weights_.dot( emap );
				residue_total_energy_[ ii ] += half;
				residue_total_energy_[ jj ] += half;
			}
		}
	}


	residue_total_energy_uptodate_ = true;
}


/// @brief update the context graphs and the energy graph according to a new
/// score function type
///
void
Energies::set_scorefxn_info( scoring::ScoreFunctionInfoOP info )
{
	using namespace scoring;
	bool deleted_a_graph( false ); // update the max_context_neighbor_cutoff_ if deleting a graph
	for ( int ii = 1; ii <= num_context_graph_types; ++ii ) {

		if ( scorefxn_info_->requires_context_graph( ContextGraphType( ii )) != required_context_graphs_[ ii ] ) {
			if ( required_context_graphs_[ ii ] ) {
				if ( ! externally_required_context_graphs_[ ii ] ) {
					debug_assert( context_graphs_[ ii ] );
					/// The new score function does not require that the context graph be maintained, nor
					/// did any external (non-energy-method ) peice of code.  It is safe to delete this
					/// graph and to stop maintaining it during score evaluations.
					required_context_graphs_[ ii ] = false;
					context_graphs_[ ii ] = 0;
					deleted_a_graph = true;
				}
			} else {
				/// Required by the score function, but not externally required, nor required by the previous
				/// score function.
				debug_assert( context_graphs_[ ii ] == 0 );
				debug_assert( externally_required_context_graphs_[ ii ] );
				require_context_graph_( ContextGraphType (ii), false );
			}
		}
		// else -- noop, score function and I already agree on context graph ii
	}

	if ( deleted_a_graph ) {
		max_context_neighbor_cutoff_ = 0;
		for ( Size ii = 1; ii <= num_context_graph_types; ++ii ) {
			if ( context_graphs_[ ii ] ) {
				if ( max_context_neighbor_cutoff_ < context_graphs_[ ii ]->neighbor_cutoff() ) {
					max_context_neighbor_cutoff_ = context_graphs_[ ii ]->neighbor_cutoff() ;
				}
			}
		}
	}

	/// Inform the EnergyGraph of the currently-active terms.  Active means that the term has a
	/// non-zero weight.
	ScoreTypes active;
	active.reserve( n_shortranged_2b_score_types );
	for ( Size ii = 1; ii <= n_shortranged_2b_score_types; ++ii ) {
		if ( info->scores_present()[ (ScoreType) ii ] != 0.0 ) {
			active.push_back( (ScoreType) ii );
		}
	}
	energy_graph_->active_score_types( active );
	scorefxn_info_ = info;
}


/// @brief Add edges to the energy_graph and the context graphs according to domain map
///
/// @details Precondition: if the graph contains any edges, then all neighbor relationships between
/// pairs of residues that have not moved relative to each other are represented by the
/// presence or absence of edges in the energy graph.  If there are no edges in the
/// energy graph, then all pair inforamtion must be calculated
/// Precondition: if two residues have changed relative to one another according to
/// the domain map, then there are no edges between them.
void
Energies::update_neighbor_links(
	pose::Pose const & pose
)
{
	using namespace graph;
	using namespace scoring;

	runtime_assert( !core::pose::symmetry::is_symmetric( pose ) );

	//std::cout << "update_neighbor_links: interaction dist: " << scorefxn_info_->max_atomic_interaction_distance() <<
	// std::endl;

	if ( point_graph_ == 0 ) {
		point_graph_ = conformation::PointGraphOP( new core::conformation::PointGraph );
	}
	fill_point_graph( pose, point_graph_ );

	// According to the domain map, add some of the edges detected by the octree to
	// the energy graph and to the context graphs

	bool all_moved( energy_graph_->num_edges() == 0 );

	utility::vector1< ContextGraphOP > context_graphs_present;
	for ( uint ii = 1, ii_end = context_graphs_.size(); ii <= ii_end; ++ii ) {
		if ( context_graphs_[ ii ] ) context_graphs_present.push_back( context_graphs_[ ii ] );
	}

	for ( uint ii = 1, ii_end = pose.total_residue(); ii <= ii_end; ++ii ) {

		int const ii_map( domain_map_(ii) );
		bool const ii_moved( ii_map == 0 || all_moved );

		Distance const iiradius( pose.residue_type( ii ).nbr_radius() );
		Distance const ii_intxn_radius( iiradius +
			scorefxn_info_->max_atomic_interaction_distance() );

		for ( core::conformation::PointGraph::UpperEdgeListConstIter
				ii_iter = point_graph_->get_vertex( ii ).upper_edge_list_begin(),
				ii_end_iter = point_graph_->get_vertex( ii ).upper_edge_list_end();
				ii_iter != ii_end_iter; ++ii_iter ) {
			uint const jj = ii_iter->upper_vertex();
			if ( ( domain_map_(jj) != ii_map ) || ii_moved ) {

				Distance const jjradius( pose.residue_type( jj ).nbr_radius() );
				DistanceSquared const square_distance( ii_iter->data().dsq() );

				// How about we simply make sure the radii sum is positive instead of paying for a sqrt
				// if ( std::sqrt( square_distance ) < ( ii_intxn_radius + jj_res.nbr_radius() ) ) {
				if ( ii_intxn_radius + jjradius > 0 ) {
					if ( square_distance < (ii_intxn_radius + jjradius )*(ii_intxn_radius + jjradius ) ) {
						energy_graph_->add_energy_edge( ii, jj, square_distance );
					}
					for ( uint kk = 1; kk <= context_graphs_present.size(); ++kk ) {
						context_graphs_present[ kk ]->conditionally_add_edge( ii, jj, square_distance );
					}
				}
			}
		}
	}
}

/// @brief determine distance cutoff threshold based on scorefxn_info_ and
/// then add edges to the PointGraph class
void
Energies::fill_point_graph( pose::Pose const & pose, conformation::PointGraphOP pg ) const {

	core::conformation::residue_point_graph_from_conformation( pose.conformation(), *pg );

	Distance const max_pair_radius = pose::pose_max_nbr_radius( pose );
	Distance const energy_neighbor_cutoff = 2 * max_pair_radius + scorefxn_info_->max_atomic_interaction_distance();

	Distance const context_cutoff = max_context_neighbor_cutoff_;

	Distance const neighbor_cutoff = numeric::max( energy_neighbor_cutoff, context_cutoff );

	// Stuarts O( n log n ) octree algorithm
	core::conformation::find_neighbors<core::conformation::PointGraphVertexData,core::conformation::PointGraphEdgeData>( pg, neighbor_cutoff );
}

void Energies::copy_nblists( Energies const & other )
{
	for ( std::map< EnergiesCacheableDataType::Enum, scoring::NeighborListOP >::const_iterator
			other_nblist_iter = other.nblist_.begin(),
			other_nblist_end = other.nblist_.end();
			other_nblist_iter != other_nblist_end; ++other_nblist_iter ) {
		nblist_[ other_nblist_iter->first ] = other_nblist_iter->second->clone();
	}
}

void Energies::copy_context_graphs( Energies const & other )
{
	debug_assert( other.context_graphs_.size() == context_graphs_.size() );
	externally_required_context_graphs_ = other.externally_required_context_graphs_;
	required_context_graphs_ = other.required_context_graphs_;
	max_context_neighbor_cutoff_ = other.max_context_neighbor_cutoff_;
	for ( Size ii = 1; ii <= other.context_graphs_.size(); ++ii ) {
		if ( other.context_graphs_[ ii ] ) {
			context_graphs_[ ii ] = other.context_graphs_[ ii ]->clone();
		} else {
			context_graphs_[ ii ] = 0;
		}
	}
}

void Energies::copy_lr_energy_containers( Energies const & other )
{
	for ( Size ii = 1; ii <= other.long_range_energy_containers_.size(); ++ii ) {
		if ( other.long_range_energy_containers_[ ii ] ) {
			long_range_energy_containers_[ ii ] = other.long_range_energy_containers_[ ii ]->clone();
		} else {
			long_range_energy_containers_[ ii ] = 0;
		}
	}
}

/// @brief Create a context graph.  If the requirement is external, someone other than a ScoreFunction
/// has declared that the context graph is needed (possibly by asking for it right now) so the graph
/// must also be made up-to-date.
void
Energies::require_context_graph_( scoring::ContextGraphType type, bool external ) const
{
	debug_assert( context_graphs_[ type ] == 0 );
	required_context_graphs_[ type ] = true;
	context_graphs_[ type ] = ContextGraphFactory::create_context_graph( type );
	if ( context_graphs_[ type ] == 0 ) {
		utility_exit_with_message( "Error: Null returned from ContextGraphFactory::create_context_graph( " + utility::to_string( type ) + ")" );
	}
	context_graphs_[ type ]->set_num_nodes( size_ );

	if ( max_context_neighbor_cutoff_ < context_graphs_[ type ]->neighbor_cutoff() ) {
		max_context_neighbor_cutoff_ = context_graphs_[ type ]->neighbor_cutoff() ;
	}

	if ( external ) {
		externally_required_context_graphs_[ type ] = true;

		using namespace graph;

		core::conformation::PointGraphOP point_graph( new core::conformation::PointGraph );
		fill_point_graph( *owner_, point_graph );
		for ( uint ii = 1, ii_end = size_; ii <= ii_end; ++ii ) {
			Distance const iiradius( owner_->residue_type( ii ).nbr_radius() );
			Distance const ii_intxn_radius( iiradius +
				scorefxn_info_->max_atomic_interaction_distance() );
			for ( core::conformation::PointGraph::UpperEdgeListConstIter
					ii_iter = point_graph->get_vertex( ii ).upper_edge_list_begin(),
					ii_end_iter = point_graph->get_vertex( ii ).upper_edge_list_end();
					ii_iter != ii_end_iter; ++ii_iter ) {
				uint const jj = ii_iter->upper_vertex();
				Distance const jjradius( owner_->residue_type( jj ).nbr_radius() );

				// Make sure to add edges only if the intereaction radius
				if ( ii_intxn_radius + jjradius > 0 ) {
					DistanceSquared const square_distance( ii_iter->data().dsq() );
					context_graphs_[ type ]->conditionally_add_edge( ii, jj, square_distance );
				}
			}
		}
	}


}


///////////////////////////////////////////////////////////////////////////////
/// @brief  This is called by ScoreFunction at the start of an energy calculation
///
/// @details
/// the responsibilities are:
/// 1. to check for size/scorefxn mismatches,
/// 2. reset some of the cached energies
/// Expect that Energies::structure_has_moved( ... ) has been called
/// already if the stucture has changed, ie if domain_map != 1


void
Energies::scoring_begin(
	scoring::ScoreFunction const & sfxn,
	pose::Pose const & pose // for the neighbor calculation
)
{
	// set our scoring flag to true
	scoring_ = true;

	total_energy_ = 0.0;
	finalized_energies_.zero();

	// check for info mismatch and save scorefxn_info_ before (possibly) calling set_size()
	// context graphs and energy graph invalidated with change of scoring function
	bool scorefunction_changed( false );
	if ( *(sfxn.info()) != (*scorefxn_info_) ) {
		//std::cout << "ScoreFunction changed!" << std::endl;
		energy_state_ = BAD;
		graph_state_ = BAD;
		set_scorefxn_info( sfxn.info() );
		scorefunction_changed = true;
	}

	/*if ( ! scorefunction_changed ) {
	/// Force a rescore if weights have changed.
	for ( Size ii = 1; ii <= n_shortranged_2b_score_types; ++ii ) {
	if ( scorefxn_weights_[ (ScoreType) ii ] != sfxn.weights()[ (ScoreType) ii ] ) {
	energy_state_ = BAD;
	graph_state_ = BAD;
	scorefunction_changed = true;
	break;
	}
	}
	}*/

	// check size
	if ( size_ != pose.total_residue() ) {
		energy_state_ = BAD;
		graph_state_ = BAD;
		set_size( pose.total_residue() );
	}

	// set our internal domain_map if the size has changed
	if ( energy_state_ == BAD ) {
		domain_map_ = 0;
		energy_state_ = MOD; // since everything in domain_map has moved, we're
		graph_state_ = MOD;  // completely correct modulo the domain_map
		internalize_new_domain_map();
		if ( scorefunction_changed ) {
			prepare_neighbor_graphs();
			update_neighbor_links( pose );
			// update connectivity in EnergyGraph, but do not mark the graph state as good.
		}
	}

	// reset total_energies_ -- should probably have been done
	total_energies_.clear();

	residue_total_energies_uptodate_ = false;
	residue_total_energy_uptodate_ = false;

	scorefxn_weights_ = sfxn.weights();

}


///////////////////////////////////////////////////////////////////////////////
void
Energies::scoring_end( scoring::ScoreFunction const &  )
{
	scoring_ = false;
	energy_state_ = GOOD;
}

bool
Energies::res_moved( int const seqpos ) const
{
	require_scoring();
	return energy_graph_->get_energy_node( seqpos )->moved();
}


} // scoring
} // core

// ///////////////////////////////////////////////////////////////////////////////
// // this is called by ScoreFunction at the start of an energy calculation
// //
// //
// void
// Energies::scoring_begin(
//  DomainMap const & domain_map_in,
//  ScoreFunctionInfo const & info,
//  bool update_neighbors,
//  pose::Pose const & pose // for the neighbor calculation
// )
// {

//  // check size
//  if ( size_ != pose.total_residue() ) {
//   energy_state_ = BAD;
//   graph_state_ = BAD;
//   set_size( pose.total_residue() );
//  }

//  // check for info mismatch
//  if ( info != scorefxn_info_ ) {
//   energy_state_ = BAD;
//   scorefxn_info_ = info;
//  }

//  // set our scoring flag to true
//  scoring_ = true;

//  // update the domain_map
//  if ( ( update_neighbors && graph_state_ == BAD ) ||
//     ( energy_state_ == BAD ) ) {
//   domain_map_ = 0;
//   energy_state_ = MOD;
//   graph_state_ = MOD;
//  } else {
//   domain_map_ = domain_map_in;
//  }


//  // reset energy/neighbor links
//  if ( !use_nblist_ ) {
//   // deletes the moving links
//   prepare_neighbor_graphs();
//  } else {
//   // pass -- preserve the nbr graphs during minimization
//  }


//  // total and twobody energies have to be reset if anything has changed
//  for ( Size i=1; i<= size_; ++i ) {
//   if ( domain_map_(i) == 0 || domain_map_(i) != domain_map_(1) ) {
//    total_energies_.clear();
//    break;
//   }
//  }

//  // onebody residue energies are still valid unless bb/chi changed
//  // for that sequence position
//  for ( Size i=1; i<= size_; ++i ) {
//   if ( res_moved(i) ) {
//    onebody_energies_[i].clear();
//   }
//  }


//  ////////////////////////////////////////////
//  // update the neighbor links if desired
//  if ( use_nblist_ ) {
//   // don't modify the graphs at all
//   // we pretend that the graph_state is good
//   if ( update_neighbors ) graph_state_ = GOOD;
//  } else if ( update_neighbors ) {
//   scoring::update_neighbor_energy_links
//    ( pose, energy_graph_, tenA_neighbor_graph_, domain_map_ );
//   graph_state_ = GOOD;
//  } else {
//   //apl do we really want to do this?
//   //if graph state was good going into this, why would it no longer be good?
//   //  rotamer trials wants to use the neighbor graph that already exists
//   //  (and it has to call this function once for each residue it modifies)
//   //
//   //graph_state_ = BAD;
//  }


// }

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::Energies::save( Archive & arc ) const {
	arc( CEREAL_NVP( size_ ) ); // Size
	// The owner_ pointer will be set by the Pose that deserializes this Energies object
	// arc( CEREAL_NVP( owner_ ) ); // pose::Pose *; raw pointer: pose::Pose *
	// EXEMPT owner_

	arc( CEREAL_NVP( energy_graph_ ) ); // EnergyGraphOP
	arc( CEREAL_NVP( context_graphs_ ) ); // utility::vector1<ContextGraphOP>
	arc( CEREAL_NVP( externally_required_context_graphs_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( required_context_graphs_ ) ); // utility::vector1<_Bool>
	arc( CEREAL_NVP( max_context_neighbor_cutoff_ ) ); // Real
	arc( CEREAL_NVP( long_range_energy_containers_ ) ); // utility::vector1<LREnergyContainerOP>
	arc( CEREAL_NVP( nblist_ ) ); // std::map<EnergiesCacheableDataType::Enum, scoring::NeighborListOP>
	arc( CEREAL_NVP( use_nblist_ ) ); // _Bool
	arc( CEREAL_NVP( use_nblist_auto_update_ ) ); // _Bool
	arc( CEREAL_NVP( minimization_graph_ ) ); // MinimizationGraphOP
	arc( CEREAL_NVP( onebody_energies_ ) ); // utility::vector1<EnergyMap>
	arc( CEREAL_NVP( residue_total_energies_uptodate_ ) ); // _Bool
	arc( CEREAL_NVP( residue_total_energies_ ) ); // utility::vector1<EnergyMap>
	arc( CEREAL_NVP( residue_total_energy_uptodate_ ) ); // _Bool
	arc( CEREAL_NVP( residue_total_energy_ ) ); // utility::vector1<Real>
	arc( CEREAL_NVP( total_energies_ ) ); // EnergyMap
	arc( CEREAL_NVP( total_energy_ ) ); // Real
	arc( CEREAL_NVP( finalized_energies_ ) ); // EnergyMap
	arc( CEREAL_NVP( scorefxn_info_ ) ); // scoring::ScoreFunctionInfoOP
	arc( CEREAL_NVP( scorefxn_weights_ ) ); // EnergyMap
	arc( CEREAL_NVP( domain_map_ ) ); // DomainMap
	arc( CEREAL_NVP( scoring_ ) ); // _Bool
	arc( CEREAL_NVP( energy_state_ ) ); // enum core::scoring::Energies::EnergyState
	arc( CEREAL_NVP( graph_state_ ) ); // enum core::scoring::Energies::EnergyState
	arc( CEREAL_NVP( data_cache_ ) ); // BasicDataCache
	arc( CEREAL_NVP( point_graph_ ) ); // conformation::PointGraphOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::Energies::load( Archive & arc ) {
	arc( size_ ); // Size

	// arc( owner_ ); // pose::Pose *; raw pointer: pose::Pose *
	// EXEMPT owner_

	arc( energy_graph_ ); // EnergyGraphOP
	arc( context_graphs_ ); // utility::vector1<ContextGraphOP>
	arc( externally_required_context_graphs_ ); // utility::vector1<_Bool>
	arc( required_context_graphs_ ); // utility::vector1<_Bool>
	arc( max_context_neighbor_cutoff_ ); // Real
	arc( long_range_energy_containers_ ); // utility::vector1<LREnergyContainerOP>
	arc( nblist_ ); // std::map<EnergiesCacheableDataType::Enum, scoring::NeighborListOP>
	arc( use_nblist_ ); // _Bool
	arc( use_nblist_auto_update_ ); // _Bool
	arc( minimization_graph_ ); // MinimizationGraphOP
	arc( onebody_energies_ ); // utility::vector1<EnergyMap>
	arc( residue_total_energies_uptodate_ ); // _Bool
	arc( residue_total_energies_ ); // utility::vector1<EnergyMap>
	arc( residue_total_energy_uptodate_ ); // _Bool
	arc( residue_total_energy_ ); // utility::vector1<Real>
	arc( total_energies_ ); // EnergyMap
	arc( total_energy_ ); // Real
	arc( finalized_energies_ ); // EnergyMap
	arc( scorefxn_info_ ); // scoring::ScoreFunctionInfoOP
	arc( scorefxn_weights_ ); // EnergyMap
	arc( domain_map_ ); // DomainMap
	arc( scoring_ ); // _Bool
	arc( energy_state_ ); // enum core::scoring::Energies::EnergyState
	arc( graph_state_ ); // enum core::scoring::Energies::EnergyState
	arc( data_cache_ ); // BasicDataCache
	arc( point_graph_ ); // conformation::PointGraphOP
}
SAVE_AND_LOAD_SERIALIZABLE( core::scoring::Energies );
CEREAL_REGISTER_TYPE( core::scoring::Energies )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_Energies )
#endif // SERIALIZATION
