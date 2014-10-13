// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/ConstraintsEnergyContainer.cc
/// @brief  Constraints Energy Container class implementation
/// @author Andrew Leaver-Fay

// Unit headers
#include <core/scoring/constraints/ConstraintEnergyContainer.hh>

// Package headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintGraph.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/pose/Pose.hh>

// STL Headers
#include <cassert>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace constraints {

CstResNeighbIterator::CstResNeighbIterator(
	Size focused_node,
	CstResNeighbIterator::EdgeListIter edge_iter
)
:
	focused_node_( focused_node ),
	edge_iter_( edge_iter )
{}


CstResNeighbIterator::~CstResNeighbIterator()
{}

ResidueNeighborIterator const &
CstResNeighbIterator::operator = ( ResidueNeighborIterator const & rhs)
{
	assert( &(dynamic_cast< CstResNeighbIterator const & > ( rhs )) );
	CstResNeighbIterator const & crni_rhs = static_cast< CstResNeighbIterator const & > ( rhs );

	focused_node_ = crni_rhs.focused_node_;
	edge_iter_ = crni_rhs.edge_iter_;
	return *this;
}

ResidueNeighborIterator const &
CstResNeighbIterator::operator ++ ()
{
	++edge_iter_;
	return *this;
}

bool
CstResNeighbIterator::operator == ( ResidueNeighborIterator const & rhs ) const
{
	assert( &( dynamic_cast< CstResNeighbIterator const & > ( rhs )) );
	CstResNeighbIterator const & crni_rhs = static_cast< CstResNeighbIterator const & > ( rhs );

	return ( edge_iter_ == crni_rhs.edge_iter_ );
}

bool
CstResNeighbIterator::operator != ( ResidueNeighborIterator const & rhs ) const
{
	assert( &( dynamic_cast< CstResNeighbIterator const & > ( rhs )) );
	CstResNeighbIterator const & crni_rhs = static_cast< CstResNeighbIterator const & > ( rhs );
	return ( edge_iter_ != crni_rhs.edge_iter_ );
}


Size
CstResNeighbIterator::upper_neighbor_id() const
{
	return (Size) (*edge_iter_)->get_second_node_ind();
}

Size
CstResNeighbIterator::lower_neighbor_id() const
{
	return (Size) (*edge_iter_)->get_first_node_ind();
}

Size
CstResNeighbIterator::residue_iterated_on() const
{
	return focused_node_;
}

Size
CstResNeighbIterator::neighbor_id() const
{
	return (Size) (*edge_iter_)->get_other_ind( (int) focused_node_ );
}


void
CstResNeighbIterator::save_energy( EnergyMap const & emap )
{
	downcast_cstedge(*edge_iter_)->bond_geometry_energy( emap[ bond_geometry ] );
	downcast_cstedge(*edge_iter_)->rna_bond_geometry_energy( emap[ rna_bond_geometry ] );
	downcast_cstedge(*edge_iter_)->atom_pair_constraint_energy( emap[ atom_pair_constraint ] );
	downcast_cstedge(*edge_iter_)->coordinate_constraint_energy( emap[ coordinate_constraint ] );
	downcast_cstedge(*edge_iter_)->angle_constraint_energy( emap[ angle_constraint ] );
	downcast_cstedge(*edge_iter_)->dihedral_constraint_energy( emap[ dihedral_constraint ] );
	downcast_cstedge(*edge_iter_)->backbone_stub_constraint_energy( emap[ backbone_stub_constraint ] );
	downcast_cstedge(*edge_iter_)->backbone_stub_linear_constraint_energy( emap[ backbone_stub_linear_constraint ] );
	downcast_cstedge(*edge_iter_)->res_type_linking_constraint_energy( emap[ res_type_linking_constraint ]);

}


void
CstResNeighbIterator::retrieve_energy( EnergyMap & emap ) const
{
	emap[ bond_geometry ] = downcast_cstedge(*edge_iter_)->bond_geometry_energy();
	emap[ rna_bond_geometry ] = downcast_cstedge(*edge_iter_)->rna_bond_geometry_energy();
	emap[ atom_pair_constraint ] = downcast_cstedge(*edge_iter_)->atom_pair_constraint_energy();
	emap[ coordinate_constraint ] = downcast_cstedge(*edge_iter_)->coordinate_constraint_energy();
	emap[ angle_constraint ]     = downcast_cstedge(*edge_iter_)->angle_constraint_energy();
	emap[ dihedral_constraint ]  = downcast_cstedge(*edge_iter_)->dihedral_constraint_energy();
	emap[ backbone_stub_constraint ]  = downcast_cstedge(*edge_iter_)->backbone_stub_constraint_energy();
	emap[ res_type_linking_constraint ]  = downcast_cstedge(*edge_iter_)->res_type_linking_constraint_energy();

}


void
CstResNeighbIterator::accumulate_energy( EnergyMap & emap ) const
{
	emap[ bond_geometry ] += downcast_cstedge(*edge_iter_)->bond_geometry_energy();
	emap[ rna_bond_geometry ] += downcast_cstedge(*edge_iter_)->rna_bond_geometry_energy();
	emap[ atom_pair_constraint ] += downcast_cstedge(*edge_iter_)->atom_pair_constraint_energy();
	emap[ coordinate_constraint ] += downcast_cstedge(*edge_iter_)->coordinate_constraint_energy();
	emap[ angle_constraint ]     += downcast_cstedge(*edge_iter_)->angle_constraint_energy();
	emap[ dihedral_constraint ]  += downcast_cstedge(*edge_iter_)->dihedral_constraint_energy();
	emap[ backbone_stub_constraint ]  += downcast_cstedge(*edge_iter_)->backbone_stub_constraint_energy();
	emap[ backbone_stub_linear_constraint ]  += downcast_cstedge(*edge_iter_)->backbone_stub_linear_constraint_energy();
	emap[ res_type_linking_constraint ]  += downcast_cstedge(*edge_iter_)->res_type_linking_constraint_energy();

}

void CstResNeighbIterator::mark_energy_computed()
{
	downcast_cstedge(*edge_iter_)->energy_computed( true );
}

void CstResNeighbIterator::mark_energy_uncomputed()
{
	downcast_cstedge(*edge_iter_)->energy_computed( false );
}


bool
CstResNeighbIterator::energy_computed() const
{
	return downcast_cstedge(*edge_iter_)->energy_computed();
}

ConstraintEdge *
CstResNeighbIterator::downcast_cstedge( graph::Edge * edge )
{
	assert( dynamic_cast< ConstraintEdge * > ( edge ) );
	return static_cast< ConstraintEdge * > (edge);
}

////////////////////////////////////////////////////////////////////////
///// Constraint Residue Neighbor Constant Iterator class implementation
////////////////////////////////////////////////////////////////////////

CstResNeighbConstIterator::CstResNeighbConstIterator(
	Size focused_node,
	CstResNeighbConstIterator::EdgeListConstIter edge_iter
)
:
	focused_node_( focused_node ),
	edge_iter_( edge_iter )
{}

CstResNeighbConstIterator::~CstResNeighbConstIterator()
{}

ResidueNeighborConstIterator const &
CstResNeighbConstIterator::operator = ( ResidueNeighborConstIterator const & rhs )
{
	//assert( dynamic_cast< CstResNeighbConstIterator const & > ( rhs ) );
	CstResNeighbConstIterator const & crnci_rhs = static_cast< CstResNeighbConstIterator const & > ( rhs );

	focused_node_ = crnci_rhs.focused_node_;
	edge_iter_ = crnci_rhs.edge_iter_;
	return *this;
}

ResidueNeighborConstIterator const &
CstResNeighbConstIterator::operator ++ ()
{
	++edge_iter_;
	return *this;
}

/// @brief returns true if the two edge-list iterators are equal
bool
CstResNeighbConstIterator::operator == ( ResidueNeighborConstIterator const & rhs ) const
{
	assert( &( dynamic_cast< CstResNeighbConstIterator const & > ( rhs )) );
	CstResNeighbConstIterator const & crnci_rhs = static_cast< CstResNeighbConstIterator const & > ( rhs );
	return ( edge_iter_ == crnci_rhs.edge_iter_ );
}


/// @brief returns true if the two edge-list iterators are not equal
bool
CstResNeighbConstIterator::operator != ( ResidueNeighborConstIterator const & rhs ) const
{
	assert( &( dynamic_cast< CstResNeighbConstIterator const & > ( rhs )) );
	CstResNeighbConstIterator const & crnci_rhs = static_cast< CstResNeighbConstIterator const & > ( rhs );
	return ( edge_iter_ != crnci_rhs.edge_iter_ );
}

Size
CstResNeighbConstIterator::upper_neighbor_id() const
{
	return (Size) (*edge_iter_)->get_second_node_ind();
}

Size
CstResNeighbConstIterator::lower_neighbor_id() const
{
	return (Size) (*edge_iter_)->get_first_node_ind();
}


Size
CstResNeighbConstIterator::residue_iterated_on() const
{
	return focused_node_;
}

Size
CstResNeighbConstIterator::neighbor_id() const
{
	return (Size) (*edge_iter_)->get_other_ind( (int) focused_node_ );
}

/// @brief overwrites the three constraint-energy positions in the emap with
/// the three contraint energies stored on the edge pointed to by the edge iter.
/// Does not zero out the other positions in the emap.
void
CstResNeighbConstIterator::retrieve_energy( EnergyMap & emap ) const
{
	emap[ bond_geometry ] = downcast_cstedge(*edge_iter_)->bond_geometry_energy();
	emap[ rna_bond_geometry ] = downcast_cstedge(*edge_iter_)->rna_bond_geometry_energy();
	emap[ atom_pair_constraint ] = downcast_cstedge(*edge_iter_)->atom_pair_constraint_energy();
	emap[ coordinate_constraint ] = downcast_cstedge(*edge_iter_)->coordinate_constraint_energy();
	emap[ angle_constraint ]     = downcast_cstedge(*edge_iter_)->angle_constraint_energy();
	emap[ dihedral_constraint ]  = downcast_cstedge(*edge_iter_)->dihedral_constraint_energy();
	emap[ backbone_stub_constraint ]  = downcast_cstedge(*edge_iter_)->backbone_stub_constraint_energy();
	emap[ backbone_stub_linear_constraint ]  = downcast_cstedge(*edge_iter_)->backbone_stub_linear_constraint_energy();
	emap[ res_type_linking_constraint ]  = downcast_cstedge(*edge_iter_)->res_type_linking_constraint_energy();

}

/// @brief accumulates the three constraint-energy positions in the emap with
/// the three contraint energies stored on the edge pointed to by the edge iter.
/// Does not touch the other positions in the emap.
void
CstResNeighbConstIterator::accumulate_energy( EnergyMap & emap ) const
{
	emap[ bond_geometry ] += downcast_cstedge(*edge_iter_)->bond_geometry_energy();
	emap[ rna_bond_geometry ] += downcast_cstedge(*edge_iter_)->rna_bond_geometry_energy();
	emap[ atom_pair_constraint ] += downcast_cstedge(*edge_iter_)->atom_pair_constraint_energy();
	emap[ coordinate_constraint ] += downcast_cstedge(*edge_iter_)->coordinate_constraint_energy();
	emap[ angle_constraint ]     += downcast_cstedge(*edge_iter_)->angle_constraint_energy();
	emap[ dihedral_constraint ]  += downcast_cstedge(*edge_iter_)->dihedral_constraint_energy();
	emap[ backbone_stub_constraint ]  += downcast_cstedge(*edge_iter_)->backbone_stub_constraint_energy();
	emap[ backbone_stub_linear_constraint ]  += downcast_cstedge(*edge_iter_)->backbone_stub_linear_constraint_energy();
	emap[ res_type_linking_constraint ]  += downcast_cstedge(*edge_iter_)->res_type_linking_constraint_energy();

}

bool
CstResNeighbConstIterator::energy_computed() const
{
	return downcast_cstedge(*edge_iter_)->energy_computed();
}

ConstraintEdge const *
CstResNeighbConstIterator::downcast_cstedge( graph::Edge const * edge )
{
	assert( dynamic_cast< ConstraintEdge const * > ( edge ) );
	return static_cast< ConstraintEdge const * > (edge);
}

/////////////////////////////////////////////////////
/// Constraints Energy Container Class Implementation
/////////////////////////////////////////////////////

CstEnergyContainer::CstEnergyContainer() : cst_graph_( /* 0 */ ), cst_set_revision_id_( 0 ), constraint_set_( 0 ) {}

bool
CstEnergyContainer::empty() const {
	return cst_graph_ == 0;
}


CstEnergyContainer::CstEnergyContainer( pose::Pose const & pose )
{
	cst_graph_ = ConstraintGraphOP( new ConstraintGraph() );
	cst_graph_->set_num_nodes( pose.total_residue() );

	constraint_set_ = pose.constraint_set();
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( ConstraintSet::ResiduePairConstraintsIterator
				rpc_iter = constraint_set_->residue_pair_constraints_begin( ii ),
				rpc_end = constraint_set_->residue_pair_constraints_end( ii );
				rpc_iter != rpc_end; ++rpc_iter ) {
			if ( ii < rpc_iter->first ) {
				cst_graph_->add_edge( ii, rpc_iter->first );
			}
		}
	}
	cst_set_revision_id_ = constraint_set_->revision_id();
}



CstEnergyContainer::~CstEnergyContainer()
{}

LREnergyContainerOP
CstEnergyContainer::clone() const
{
	CstEnergyContainerOP cstec( new CstEnergyContainer );
	if ( !empty() ) {
		cstec->cst_graph_ = ConstraintGraphOP( new ConstraintGraph( *cst_graph_ ) );
		cstec->cst_set_revision_id_ =  cst_set_revision_id_;
		cstec->constraint_set_ = constraint_set_;
	}
	return cstec;
}

void
CstEnergyContainer::set_num_nodes( Size newsize )
{
	if ( !empty() ) {
		cst_graph_->set_num_nodes( newsize );
		cst_set_revision_id_ = 0;
		constraint_set_.reset(); // flag that the CstEnergyContainer needs to be recreated.
	}
}

bool
CstEnergyContainer::any_neighbors_for_residue( int resid ) const
{
	assert( !empty() );
	return cst_graph_->get_node( resid )->num_edges() != 0;
}

bool
CstEnergyContainer::any_upper_neighbors_for_residue( int resid ) const
{
	assert( !empty() );
	return cst_graph_->get_node( resid )->get_num_edges_to_larger_indexed_nodes() != 0;
}


ResidueNeighborConstIteratorOP
CstEnergyContainer::const_neighbor_iterator_begin( int resid ) const
{
	assert( !empty() );
	return ResidueNeighborConstIteratorOP( new CstResNeighbConstIterator( resid, cst_graph_->get_node( resid )->const_edge_list_begin() ) );
}

ResidueNeighborConstIteratorOP
CstEnergyContainer::const_neighbor_iterator_end( int resid ) const
{
	assert( !empty() );
	return ResidueNeighborConstIteratorOP( new CstResNeighbConstIterator( resid, cst_graph_->get_node( resid )->const_edge_list_end() ) );
}

ResidueNeighborConstIteratorOP
CstEnergyContainer::const_upper_neighbor_iterator_begin( int resid ) const
{
	assert( !empty() );
	return ResidueNeighborConstIteratorOP( new CstResNeighbConstIterator( resid, cst_graph_->get_node( resid )->const_upper_edge_list_begin() ) );
}

ResidueNeighborConstIteratorOP
CstEnergyContainer::const_upper_neighbor_iterator_end( int resid ) const
{
	assert( !empty() );
	return ResidueNeighborConstIteratorOP( new CstResNeighbConstIterator( resid, cst_graph_->get_node( resid )->const_upper_edge_list_end() ) );
}

ResidueNeighborIteratorOP
CstEnergyContainer::neighbor_iterator_begin( int resid )
{
	assert( !empty() );
	return ResidueNeighborIteratorOP( new CstResNeighbIterator( resid, cst_graph_->get_node( resid )->edge_list_begin() ) );
}

ResidueNeighborIteratorOP
CstEnergyContainer::neighbor_iterator_end( int resid )
{
	assert( !empty() );
	return ResidueNeighborIteratorOP( new CstResNeighbIterator( resid, cst_graph_->get_node( resid )->edge_list_end() ) );
}

ResidueNeighborIteratorOP
CstEnergyContainer::upper_neighbor_iterator_begin( int resid )
{
	assert( !empty() );
	return ResidueNeighborIteratorOP( new CstResNeighbIterator( resid, cst_graph_->get_node( resid )->upper_edge_list_begin() ) );
}

ResidueNeighborIteratorOP
CstEnergyContainer::upper_neighbor_iterator_end( int resid )
{
	assert( !empty() );
	return ResidueNeighborIteratorOP( new CstResNeighbIterator( resid, cst_graph_->get_node( resid )->upper_edge_list_end() ) ); // pbmod
}

bool
CstEnergyContainer::matches( ConstraintSetCOP cst_set )
{
	assert( !empty() );
	return ( cst_set.get() == constraint_set_.get() && cst_set_revision_id_ == cst_set->revision_id() );
}


}
}
}

