// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/ConstraintsEnergyContainer.hh
/// @brief  Constraints Energy Container class declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_constraints_ConstraintEnergyContainer_hh
#define INCLUDED_core_scoring_constraints_ConstraintEnergyContainer_hh

// Unit headers
#include <core/scoring/constraints/ConstraintEnergyContainer.fwd.hh>

// Package headers
#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/constraints/ConstraintGraph.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <core/scoring/EnergyMap.fwd.hh>


// Project headers
#include <core/pose/Pose.fwd.hh>
#include <utility/graph/Graph.hh>

#include <utility/vector1.hh>
#include <iosfwd>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace constraints {

class CstResNeighbIterator : public ResidueNeighborIterator {
	CstResNeighbIterator & operator = (CstResNeighbIterator const & );
public:
	typedef utility::graph::Node::EdgeListIter EdgeListIter;

public:

	CstResNeighbIterator( Size focused_node, EdgeListIter edge_iter);

	~CstResNeighbIterator() override;

	ResidueNeighborIterator & operator = ( ResidueNeighborIterator const & ) override;
	ResidueNeighborIterator const & operator ++ () override;
	bool operator == ( ResidueNeighborIterator const & ) const override;
	bool operator != ( ResidueNeighborIterator const & ) const override;

	Size upper_neighbor_id() const override;
	Size lower_neighbor_id() const override;

	Size residue_iterated_on() const override;
	Size neighbor_id() const override;

	void save_energy( EnergyMap const & ) override;
	void retrieve_energy( EnergyMap & ) const override;
	void accumulate_energy( EnergyMap & ) const override;

	void mark_energy_computed() override;
	void mark_energy_uncomputed() override;

	bool energy_computed() const override;
protected:

	static ConstraintEdge * downcast_cstedge( utility::graph::Edge * );

private:
	Size focused_node_;
	EdgeListIter edge_iter_;

};

class CstResNeighbConstIterator : public ResidueNeighborConstIterator {
	CstResNeighbConstIterator & operator = (CstResNeighbConstIterator const & );
public:
	typedef utility::graph::Node::EdgeListConstIter EdgeListConstIter;

public:
	CstResNeighbConstIterator( Size focused_node, EdgeListConstIter edge_iter);

	~CstResNeighbConstIterator() override;

	ResidueNeighborConstIterator & operator = ( ResidueNeighborConstIterator const & ) override;
	ResidueNeighborConstIterator const & operator ++ () override;
	bool operator == ( ResidueNeighborConstIterator const & ) const override;
	bool operator != ( ResidueNeighborConstIterator const & ) const override;

	Size upper_neighbor_id() const override;
	Size lower_neighbor_id() const override;

	Size residue_iterated_on() const override;
	Size neighbor_id() const override;

	void retrieve_energy( EnergyMap & ) const override;
	void accumulate_energy( EnergyMap & ) const override;

	bool energy_computed() const override;

protected:
	static ConstraintEdge const * downcast_cstedge( utility::graph::Edge const * );

private:
	Size focused_node_;
	utility::graph::Node::EdgeListConstIter edge_iter_; // no need to store a const iterator if this class guarantees that no non-const methods are called

};

class CstEnergyContainer : public LREnergyContainer {

public:
	CstEnergyContainer();

	CstEnergyContainer( pose::Pose const & );

	~CstEnergyContainer() override;

	bool empty() const override;

	LREnergyContainerOP clone() const override;

	void
	set_num_nodes( Size ) override;

	bool
	any_neighbors_for_residue( int resid ) const override;

	bool
	any_upper_neighbors_for_residue( int resid ) const override;

	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const override;

	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const override;

	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const override;

	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_end( int resid ) const override;

	ResidueNeighborIteratorOP
	neighbor_iterator_begin( int resid ) override;

	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid ) override;

	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid ) override;

	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid ) override;

	/// Does the constraint graph that this CEC defines match up with the constraint set stored
	/// in the pose?
	bool
	matches( ConstraintSetCOP cst_set );

private:
	ConstraintGraphOP cst_graph_;

	Size cst_set_revision_id_;
	ConstraintSetCOP constraint_set_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

}
}
}

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_constraints_ConstraintEnergyContainer )
#endif // SERIALIZATION


#endif
