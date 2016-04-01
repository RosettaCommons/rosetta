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

// Unit headers
#include <core/scoring/constraints/ConstraintGraph.hh>


namespace core {
namespace scoring {
namespace constraints {

ConstraintNode::ConstraintNode( graph::Graph* owner, Size node_id )
:
	parent( owner, node_id )
{}


ConstraintNode::~ConstraintNode()
{}

void
ConstraintNode::copy_from( Node const * source )
{
	parent::copy_from( source ); // I can't remember if this is supposed to be recursive
}

Size
ConstraintNode::count_static_memory() const
{
	return sizeof ( ConstraintNode );
}

Size
ConstraintNode::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}

ConstraintEdge::~ConstraintEdge()
{}

ConstraintEdge::ConstraintEdge( graph::Graph * owner, Size first_node_ind, Size second_node_ind)
:
	parent( owner, first_node_ind, second_node_ind ),
	bond_geometry_energy_( 0.0 ),
	rna_bond_geometry_energy_( 0.0 ),
	atom_pair_constraint_energy_( 0.0 ),
	coordinate_constraint_energy_( 0.0 ),
	angle_constraint_energy_( 0.0 ),
	dihedral_constraint_energy_( 0.0 ),
	backbone_stub_constraint_energy_( 0.0 ),
	backbone_stub_linear_constraint_energy_( 0.0 ),
	res_type_linking_constraint_energy_( 0.0 ),
	metalbinding_constraint_energy_( 0.0 ),
	energy_computed_( false )
{
}

ConstraintEdge::ConstraintEdge( graph::Graph * owner, ConstraintEdge const & example_edge )
:
	parent( owner, example_edge.get_first_node_ind(), example_edge.get_second_node_ind() ),
	bond_geometry_energy_( example_edge.bond_geometry_energy_ ),
	rna_bond_geometry_energy_( example_edge.rna_bond_geometry_energy_ ),
	atom_pair_constraint_energy_( example_edge.atom_pair_constraint_energy_ ),
	coordinate_constraint_energy_( example_edge.coordinate_constraint_energy_ ),
	angle_constraint_energy_( example_edge.angle_constraint_energy_ ),
	dihedral_constraint_energy_( example_edge.dihedral_constraint_energy_ ),
	backbone_stub_constraint_energy_( example_edge.backbone_stub_constraint_energy_ ),
	backbone_stub_linear_constraint_energy_( example_edge.backbone_stub_linear_constraint_energy_ ),
	res_type_linking_constraint_energy_( example_edge.res_type_linking_constraint_energy_ ),
	metalbinding_constraint_energy_( example_edge.metalbinding_constraint_energy_ ),
	energy_computed_( example_edge.energy_computed_ )
{}

void
ConstraintEdge::copy_from( graph::Edge const * source )
{
	ConstraintEdge const * cst_source = static_cast< ConstraintEdge const * > ( source );

	bond_geometry_energy_ = cst_source->bond_geometry_energy_;
	rna_bond_geometry_energy_ = cst_source->rna_bond_geometry_energy_;
	atom_pair_constraint_energy_ = cst_source->atom_pair_constraint_energy_;
	coordinate_constraint_energy_ = cst_source->coordinate_constraint_energy_;
	angle_constraint_energy_ = cst_source->angle_constraint_energy_;
	dihedral_constraint_energy_ = cst_source->dihedral_constraint_energy_;
	backbone_stub_constraint_energy_ = cst_source->backbone_stub_constraint_energy_;
	backbone_stub_linear_constraint_energy_ = cst_source->backbone_stub_linear_constraint_energy_;
	res_type_linking_constraint_energy_ = cst_source->res_type_linking_constraint_energy_;
	metalbinding_constraint_energy_ = cst_source->metalbinding_constraint_energy_;
	energy_computed_ = cst_source->energy_computed_;
}

Size
ConstraintEdge::count_static_memory() const
{
	return sizeof ( ConstraintEdge );
}

Size
ConstraintEdge::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}

void
ConstraintEdge::bond_geometry_energy( Energy setting )
{
	bond_geometry_energy_ = setting;
}

void
ConstraintEdge::rna_bond_geometry_energy( Energy setting )
{
	rna_bond_geometry_energy_ = setting;
}


void
ConstraintEdge::atom_pair_constraint_energy( Energy setting )
{
	atom_pair_constraint_energy_ = setting;
}

void
ConstraintEdge::coordinate_constraint_energy( Energy setting )
{
	coordinate_constraint_energy_ = setting;
}

void
ConstraintEdge::angle_constraint_energy( Energy setting )
{
	angle_constraint_energy_ = setting;
}

void
ConstraintEdge::dihedral_constraint_energy( Energy setting )
{
	dihedral_constraint_energy_ = setting;
}

void
ConstraintEdge::backbone_stub_constraint_energy( Energy setting )
{
	backbone_stub_constraint_energy_ = setting;
}

void
ConstraintEdge::backbone_stub_linear_constraint_energy( Energy setting )
{
	backbone_stub_linear_constraint_energy_ = setting;
}

void
ConstraintEdge::res_type_linking_constraint_energy( Energy setting )
{
	res_type_linking_constraint_energy_ = setting;
}

void
ConstraintEdge::metalbinding_constraint_energy( Energy setting )
{
	metalbinding_constraint_energy_ = setting;
}

Energy
ConstraintEdge::bond_geometry_energy() const
{
	return bond_geometry_energy_;
}

Energy
ConstraintEdge::rna_bond_geometry_energy() const
{
	return rna_bond_geometry_energy_;
}


Energy
ConstraintEdge::atom_pair_constraint_energy() const
{
	return atom_pair_constraint_energy_;
}

Energy
ConstraintEdge::coordinate_constraint_energy() const
{
	return coordinate_constraint_energy_;
}

Energy
ConstraintEdge::angle_constraint_energy() const
{
	return angle_constraint_energy_;
}

Energy
ConstraintEdge::dihedral_constraint_energy() const
{
	return dihedral_constraint_energy_;
}

Energy
ConstraintEdge::backbone_stub_constraint_energy() const
{
	return backbone_stub_constraint_energy_;
}

Energy
ConstraintEdge::backbone_stub_linear_constraint_energy() const
{
	return backbone_stub_linear_constraint_energy_;
}

Energy
ConstraintEdge::res_type_linking_constraint_energy() const
{
	return res_type_linking_constraint_energy_;
}

Energy
ConstraintEdge::metalbinding_constraint_energy() const
{
	return metalbinding_constraint_energy_;
}


void ConstraintEdge::energy_computed( bool setting )
{
	energy_computed_ = setting;
}

bool ConstraintEdge::energy_computed() const
{
	return energy_computed_;
}


ConstraintGraph::ConstraintGraph() : parent() {}

ConstraintGraph::ConstraintGraph( Size num_nodes )
:
	parent()
{
	set_num_nodes( num_nodes );
}

ConstraintGraph::ConstraintGraph( ConstraintGraph const & source ) : parent() {
	parent::operator = ( source );
}

ConstraintGraph &
ConstraintGraph::operator = ( ConstraintGraph const & source )
{
	parent::operator = ( source );
	return *this;
}

ConstraintGraph::~ConstraintGraph() { delete_everything(); }

void ConstraintGraph::delete_edge( graph::Edge * edge )
{
	delete edge;
}


Size
ConstraintGraph::count_static_memory() const
{
	return sizeof( ConstraintGraph );
}

Size
ConstraintGraph::count_dynamic_memory() const
{
	return parent::count_dynamic_memory();
}

graph::Node*
ConstraintGraph::create_new_node( Size node_index )
{
	return new ConstraintNode( this, node_index );
}

graph::Edge*
ConstraintGraph::create_new_edge( Size index1, Size index2)
{
	return new ConstraintEdge( this, index1, index2 );
}


graph::Edge*
ConstraintGraph::create_new_edge( graph::Edge const * example_edge )
{
	return new ConstraintEdge(
		this,
		* ( static_cast< ConstraintEdge const * > (example_edge) )
	);
}


}
}
}
