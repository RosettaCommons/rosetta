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

#ifndef INCLUDED_core_scoring_constraints_ConstraintGraph_hh
#define INCLUDED_core_scoring_constraints_ConstraintGraph_hh

// Unit headers
#include <core/scoring/constraints/ConstraintGraph.fwd.hh>

// Project headers
#include <core/graph/Graph.hh>
#include <core/types.hh>

namespace core {
namespace scoring {
namespace constraints {

class ConstraintNode : public graph::Node
{

public:
	typedef graph::Node parent;
	typedef graph::Node Node;

public:
	ConstraintNode( graph::Graph*, Size node_id );
	virtual ~ConstraintNode();

	virtual void copy_from( Node const * source );

	virtual Size count_static_memory() const;
	virtual Size count_dynamic_memory() const;

};

class ConstraintEdge : public graph::Edge
{

public:
	typedef graph::Edge parent;
	typedef graph::Edge Edge;

public:
	ConstraintEdge( graph::Graph* owner, Size first_node_ind, Size second_node_ind);
	ConstraintEdge( graph::Graph* owner, ConstraintEdge const & example_edge );
	virtual ~ConstraintEdge();

	virtual void copy_from( Edge const * source );

	virtual Size count_static_memory() const;
	virtual Size count_dynamic_memory() const;

	void rna_bond_geometry_energy( Energy );
	void atom_pair_constraint_energy( Energy );
	void coordinate_constraint_energy( Energy );
	void angle_constraint_energy( Energy );
	void dihedral_constraint_energy( Energy );
	void backbone_stub_constraint_energy( Energy );
	void backbone_stub_linear_constraint_energy( Energy );
	void res_type_linking_constraint_energy( Energy );


	Energy rna_bond_geometry_energy() const;
	Energy atom_pair_constraint_energy() const;
	Energy coordinate_constraint_energy() const;
	Energy angle_constraint_energy() const;
	Energy dihedral_constraint_energy() const;
	Energy backbone_stub_constraint_energy() const;
	Energy backbone_stub_linear_constraint_energy() const;
	Energy res_type_linking_constraint_energy() const;

	void energy_computed( bool setting );
	bool energy_computed() const;

private:

	Energy rna_bond_geometry_energy_;
	Energy atom_pair_constraint_energy_;
	Energy coordinate_constraint_energy_;
	Energy angle_constraint_energy_;
	Energy dihedral_constraint_energy_;
	Energy backbone_stub_constraint_energy_;
	Energy backbone_stub_linear_constraint_energy_;
	Energy res_type_linking_constraint_energy_;


	bool energy_computed_;

};

class ConstraintGraph : public graph::Graph
{

public:
	typedef graph::Graph parent;
	typedef graph::Graph Graph;

public:
	ConstraintGraph();
	ConstraintGraph(Size num_nodes);
	ConstraintGraph( ConstraintGraph const & source );
	ConstraintGraph & operator = ( ConstraintGraph const & source );

	virtual ~ConstraintGraph();

	virtual void delete_edge( graph::Edge * edge );

	//virtual void copy_from( Graph const * source ); //? why haven't I implemented something like this in the base class?

protected:

	virtual Size count_static_memory() const;
	virtual Size count_dynamic_memory() const;

	virtual graph::Node * create_new_node( Size node_index );
	virtual graph::Edge * create_new_edge( Size index1, Size index2);
	virtual graph::Edge * create_new_edge( graph::Edge const * example_edge );


};

}
}
}

#endif
