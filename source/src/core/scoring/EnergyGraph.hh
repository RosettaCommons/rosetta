// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/EnergyGraph.hh
/// @brief  Energy graph class declaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_EnergyGraph_hh
#define INCLUDED_core_scoring_EnergyGraph_hh

// Unit Headers
#include <core/scoring/EnergyGraph.fwd.hh>

// Project Headers
#include <core/scoring/EnergyMap.hh>
#include <core/graph/Graph.hh>
#include <core/graph/ArrayPool.hh>
#include <core/scoring/ScoreType.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

/// Class EnergyNode holds the result of a domainmap update from the
/// Conformation object held by a pose; if the internal degrees of freedom
/// for a residue (corresponding to a node in this graph) have changed
/// (and are marked with color "0" in the domainmap), then the EnergyNode
/// object will hold that information for the ScoringFunction to retrieve
class EnergyNode : public graph::Node
{
public:
	typedef graph::Node parent;

public:
	EnergyNode( graph::Graph * owner, Size index );
	virtual ~EnergyNode();
	virtual void copy_from( parent const * source );

	virtual void print() const;
	virtual Size count_static_memory() const;
	virtual Size count_dynamic_memory() const;

	bool moved() const; // getter
	void moved( bool ); // setter

#ifdef SERIALIZATION
	/// @brief Serialization function, but one that is called by the graph so that
	/// the EnergyNode object can save and then later restore its data.
	/// The archive is not expected to construct the EnergyNode iself, so instead
	/// of calling this function "save" as the archive would need if it were to
	/// perform the construction itself, it's named "save_to_archive".
	template < class Archive >
	void save_to_archive( Archive & archive ) const;

	/// @brief Deserialization function, but one that is called on the already-constructed
	/// EnergyNode object instead of being called by the archive directly.
	template < class Archive >
	void load_from_archive( Archive & archive );
#endif // SERIALIZATION

private:
	bool moved_;

};

/// Class EnergyEdge holds scores for pair interactions for short-ranged energy
/// functions.  It also records whether or not it has been scored -- when an
/// edge is added to the graph, it sets "energies_computed_" as false, and the
/// ScoreFunction class marks edges as having their energies computed once it
/// computes them.
class EnergyEdge : public graph::Edge
{
public:
	typedef graph::Edge parent;
	typedef scoring::EnergyMap EnergyMap;
public:
	EnergyEdge( EnergyGraph * owner, Size n1, Size n2 );
	EnergyEdge( EnergyGraph * owner, EnergyEdge const & example_edge );
	virtual ~EnergyEdge();

	/// @brief Copy the data held on the example edge, source.
	/// The source edge must be castable to class EnergyEdge.
	virtual void copy_from( parent const * source );

	/// @brief Store the energies held in the input emap on this edge; only
	/// those ScoreTypes which are active are stored.
	inline void store_active_energies( EnergyMap const & emap );

	/// @brief Store the intersection of the energies held in the input emap on this edge:
	/// The intersection is between the set of active ScoreTypes and the ScoreTypes given
	/// in the input "subset" list.  subset need not be sorted.
	inline void store_active_energies( EnergyMap const & emap, ScoreTypes const & subset );

	/// @brief Load an energy map with the non-zero
	inline EnergyMap fill_energy_map() const;

	/// @brief Add the non-zero elements into the energy map
	inline void add_to_energy_map( EnergyMap & emap ) const;

	/// @brief Add the non-zero elements into the energy map
	inline void add_to_energy_map( EnergyMap & emap, ScoreTypes const & subset ) const;

	/// @brief Read the value stored on this edge for a particular score type.
	inline Real operator[] ( ScoreType st ) const;

	/// @brief Compute the weighted energy for the components stored on this edge
	inline Real dot( EnergyMap const & weights ) const;

	//EnergyMap & energy_map() { return energy_map_;}
	//EnergyMap const & energy_map() const { return energy_map_;}

	DistanceSquared square_distance() const { return dsqr_; }
	void square_distance( DistanceSquared dsqr ) { dsqr_ = dsqr; }

	void mark_energies_computed();
	void mark_energies_uncomputed();

	bool energies_not_yet_computed() const { return energies_not_yet_computed_;}

	virtual Size count_static_memory() const;
	virtual Size count_dynamic_memory() const;

#ifdef    SERIALIZATION
	/// @brief Serialization function, but one that is called by the graph so that
	/// the EnergyEdge object can save and then later restore its data.
	/// The archive is not expected to construct the EnergyEdge iself, so instead
	/// of calling this function "save" as the archive would need if it were to
	/// perform the construction itself, it's named "save_to_archive".
	template < class Archive >
	void save_to_archive( Archive & archive ) const;

	/// @brief Deserialization function, but one that is called on the already-constructed
	/// EnergyEdge object instead of being called by the archive directly.
	template < class Archive >
	void load_from_archive( Archive & archive );
#endif // SERIALIZATION

protected:
	/// Downcasts

	inline
	EnergyGraph const *
	get_energy_owner() const;

	inline
	EnergyGraph *
	get_energy_owner();

private:

	bool energies_not_yet_computed_;
	DistanceSquared dsqr_; // measured between the nbr atoms
	graph::ArrayPoolElement< Real > array_;
};

/// @brief Class to hold the component energies between pairs of residues.
/// Each node represents a residue in its corresponding structure.
/// Each edge in the graph holds a two-body energy map representing the
/// unweighted components of the energy function for those terms with non-zero
/// weight.  The EnergyGraph may be accessed from the pose's Energies object,
/// but at a price of an extra score evaluation.  This second score evaluation
/// may be avoided if you use the ScoreFunction::score_components( pose ) method.
class EnergyGraph : public graph::Graph
{
public:
	typedef graph::Graph parent;

public:
	virtual ~EnergyGraph();

	EnergyGraph( Size num_nodes );
	EnergyGraph();
	EnergyGraph( EnergyGraph const & src );

	EnergyGraph & operator = ( EnergyGraph const & rhs );

	//bool energy_exists( scoring::ScoreType const & type ) const;

	Size n_active_score_types() const { return active_2b_score_types_.size(); }

	/// @brief Set the active score types, and return true if the new score types
	/// are the same as the old score types and the graph state is still good.
	/// Returns false if the score types have changed, indicating that the graph
	/// has dropped all of its edges;
	bool active_score_types( ScoreTypes const & active );

	/// @brief Add an energy edge to the graph and store the square distance
	void
	add_energy_edge( Size index1, Size index2, DistanceSquared dsq );

	/// @brief Add an energy edge to the graph and set the energies
	/// for the non-zero-weighted components of the input emap.
	/// void
	/// add_energy_edge( Size index1, Size index2, EnergyMap const & emap );

	inline
	EnergyNode const *
	get_energy_node( Size index ) const
	{
		return static_cast< EnergyNode const * > ( get_node( index ));
	}

	inline
	EnergyNode *
	get_energy_node( Size index )
	{
		return static_cast< EnergyNode * > ( get_node( index ));
	}


	EnergyEdge * find_energy_edge( Size n1, Size n2);
	EnergyEdge const * find_energy_edge( Size n1, Size n2) const;

	virtual void delete_edge( graph::Edge * edge );

	/// @brief As an edge is deleted from the graph, it must reliquish hold over its
	/// array-pool element so that the element may be reused by new edges.
	void deallocate_arraypoolelement( graph::ArrayPoolElement< Real > & element );

	utility::vector1< int > const &
	score_type_2_active() const {
		return score_type_2_active_;
	}

	ScoreTypes const &
	active_2b_score_types() const {
		return active_2b_score_types_;
	}

	/// @brief Give non-const access to the array pool -- this function should only
	/// be called by class EnergyEdge.  I wish C++ let me declare this function private
	/// and that EnergyEdge could be a "friend" of this function.
	graph::ArrayPool< Real > & array_pool() { return energy_array_pool_; }

#ifdef     SERIALIZATION
	template < class Archive >
	void save( Archive & archive ) const;

	template < class Archive >
	void load( Archive & archive );
#endif // SERIALIZATION

protected:
	virtual Size count_static_memory() const;
	virtual Size count_dynamic_memory() const;

	virtual graph::Node * create_new_node( Size index );
	virtual graph::Edge * create_new_edge( Size index1, Size index2 );
	virtual graph::Edge * create_new_edge( graph::Edge const * example_edge );

private:

	boost::unordered_object_pool< EnergyEdge > * energy_edge_pool_;
	graph::ArrayPool< Real >                     energy_array_pool_;

	ScoreTypes              active_2b_score_types_;
	/// @brief these are flag values; <0 has a special meaning, so they need to be ints
	utility::vector1< int > score_type_2_active_;
};


/// @details Store only the active terms -- active meaning with non-zero weight.
inline
void EnergyEdge::store_active_energies( EnergyMap const & emap ) {
	ScoreTypes const & active( get_energy_owner()->active_2b_score_types() );

	for ( Size ii = 1, iilag = 0; ii <= active.size(); ++ii, ++iilag ) {
		array_[ iilag ] = emap[ active[ ii ] ];
	}

}

/// @details Subset may specify ScoreTypes that are not active, but only the active
/// ScoreTypes have their eneriges stored.
inline void
EnergyEdge::store_active_energies( EnergyMap const & emap, ScoreTypes const & subset )
{
	utility::vector1< int > const & st2active( get_energy_owner()->score_type_2_active() );

	for ( Size ii = 1, iilag = 0; ii <= subset.size(); ++ii, ++iilag ) {
		if ( st2active[ subset[ ii ]] >= 0 ) {
			array_[ st2active[ subset[ ii ]] ] = emap[ subset[ ii ] ];
		}
	}

}


/// @details Only load the non-zero terms -- the zero entries are already
inline
EnergyMap
EnergyEdge::fill_energy_map() const
{
	EnergyMap emap;
	ScoreTypes const & active( get_energy_owner()->active_2b_score_types() );

	for ( Size ii = 1, iilag = 0; ii <= active.size(); ++ii, ++iilag ) {
		emap[ active[ ii ] ] = array_[ iilag ];
	}
	return emap;
}

/// @brief Add the non-zero elements into the energy map
inline
void EnergyEdge::add_to_energy_map( EnergyMap & emap ) const
{
	ScoreTypes const & active( get_energy_owner()->active_2b_score_types() );

	for ( Size ii = 1, iilag = 0; ii <= active.size(); ++ii, ++iilag ) {
		emap[ active[ ii ] ] += array_[ iilag ];
	}
}

/// @brief Add the non-zero elements into the energy map
inline
void EnergyEdge::add_to_energy_map( EnergyMap & emap, ScoreTypes const & subset ) const
{
	utility::vector1< int > const & st2active( get_energy_owner()->score_type_2_active() );

	for ( Size ii = 1, iilag = 0; ii <= subset.size(); ++ii, ++iilag ) {
		if ( st2active[ subset[ ii ]] >= 0 ) {
			emap[ subset[ ii ] ] += array_[ st2active[ subset[ ii ]] ];
		}
	}

}


/// @details The owner stores a map from score types to indices indicating
/// score types not represented in the mapping with an index of -1.
inline
Real
EnergyEdge::operator[ ] ( ScoreType st ) const {
	int aid = get_energy_owner()->score_type_2_active()[ st ];
	if ( aid >= 0 ) {
		return array_[ aid ];
	} else {
		return 0;
	}
}

/// @brief Compute the weighted energy for the components stored on this edge
inline
Real
EnergyEdge::dot( EnergyMap const & weights ) const
{
	Real weighted_sum( 0.0 );
	ScoreTypes const & active( get_energy_owner()->active_2b_score_types() );

	for ( Size ii = 1, iilag = 0; ii <= active.size(); ++ii, ++iilag ) {
		weighted_sum += weights[ active[ ii ] ] * array_[ iilag ];
	}
	return weighted_sum;
}


EnergyGraph const *
EnergyEdge::get_energy_owner() const {
	return static_cast< EnergyGraph const * > (get_owner());
}


EnergyGraph *
EnergyEdge::get_energy_owner() {
	return static_cast< EnergyGraph * > (get_owner());
}


} //namespace scoring
} //namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_EnergyGraph )
#endif

#endif

