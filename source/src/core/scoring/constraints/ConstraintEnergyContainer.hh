// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
// AUTO-REMOVED #include <core/scoring/ScoreFunction.fwd.hh>

// AUTO-REMOVED #include <core/scoring/ScoreType.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/graph/Graph.hh>

#include <utility/vector1.hh>
#include <iostream>


namespace core {
namespace scoring {
namespace constraints {

class CstResNeighbIterator : public ResidueNeighborIterator {
public:
	typedef graph::Node::EdgeListIter EdgeListIter;

public:

	CstResNeighbIterator( Size focused_node, EdgeListIter edge_iter);

	virtual ~CstResNeighbIterator();

	virtual ResidueNeighborIterator const & operator = ( ResidueNeighborIterator const & );
	virtual ResidueNeighborIterator const & operator ++ ();
	virtual bool operator == ( ResidueNeighborIterator const & ) const;
	virtual bool operator != ( ResidueNeighborIterator const & ) const;

	virtual Size upper_neighbor_id() const;
	virtual Size lower_neighbor_id() const;

	virtual Size residue_iterated_on() const;
	virtual Size neighbor_id() const;

	virtual void save_energy( EnergyMap const & );
	virtual void retrieve_energy( EnergyMap & ) const;
	virtual void accumulate_energy( EnergyMap & ) const;

	virtual void mark_energy_computed();
	virtual void mark_energy_uncomputed();

	virtual bool energy_computed() const;
protected:

	static ConstraintEdge * downcast_cstedge( graph::Edge * );

private:
	Size focused_node_;
	EdgeListIter edge_iter_;

};

class CstResNeighbConstIterator : public ResidueNeighborConstIterator {
public:
	typedef graph::Node::EdgeListConstIter EdgeListConstIter;

public:
	CstResNeighbConstIterator( Size focused_node, EdgeListConstIter edge_iter);

	virtual ~CstResNeighbConstIterator();

	virtual ResidueNeighborConstIterator const & operator = ( ResidueNeighborConstIterator const & );
	virtual ResidueNeighborConstIterator const & operator ++ ();
	virtual bool operator == ( ResidueNeighborConstIterator const & ) const;
	virtual bool operator != ( ResidueNeighborConstIterator const & ) const;

	virtual Size upper_neighbor_id() const;
	virtual Size lower_neighbor_id() const;

	virtual Size residue_iterated_on() const;
	virtual Size neighbor_id() const;

	virtual void retrieve_energy( EnergyMap & ) const;
	virtual void accumulate_energy( EnergyMap & ) const;

	virtual bool energy_computed() const;

protected:
	static ConstraintEdge const * downcast_cstedge( graph::Edge const * );

private:
	Size focused_node_;
	graph::Node::EdgeListConstIter edge_iter_; // no need to store a const iterator if this class guarantees that no non-const methods are called

};

class CstEnergyContainer : public LREnergyContainer {

public:
	CstEnergyContainer();

	CstEnergyContainer( pose::Pose const & );

	virtual
	~CstEnergyContainer();

	virtual
	bool empty() const;

	virtual
	LREnergyContainerOP clone() const;

	virtual
	void
	set_num_nodes( Size );

	virtual
	bool
	any_neighbors_for_residue( int resid ) const;

	virtual
	bool
	any_upper_neighbors_for_residue( int resid ) const;

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_begin( int resid ) const;

	virtual
	ResidueNeighborConstIteratorOP
	const_neighbor_iterator_end( int resid ) const;

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_begin( int resid ) const;

	virtual
	ResidueNeighborConstIteratorOP
	const_upper_neighbor_iterator_end( int resid ) const;

	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_begin( int resid );

	virtual
	ResidueNeighborIteratorOP
	neighbor_iterator_end( int resid );

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_begin( int resid );

	virtual
	ResidueNeighborIteratorOP
	upper_neighbor_iterator_end( int resid );

	/// Does the constraint graph that this CEC defines match up with the constraint set stored
	/// in the pose?
	bool
	matches( ConstraintSetCOP cst_set );

private:
	ConstraintGraphOP cst_graph_;

	Size cst_set_revision_id_;
	ConstraintSetCOP constraint_set_;

};

}
}
}

#endif
