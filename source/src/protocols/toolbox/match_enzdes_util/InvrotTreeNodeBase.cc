// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.cc
/// @brief  Forward declaration for inverse rotamer tree node base
/// @author Florian Richter, flosopher@gmail.com, mar 2012

/// unit headers
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.hh>

//project headers
#include <core/conformation/Residue.hh>

//utility headers
#include <utility/exit.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


InvrotCollector::InvrotCollector(
  Size num_residue_lists )
{
  invrots_.clear();
  invrots_.resize( num_residue_lists );
  owner_nodes_and_locations_.clear();
}

InvrotCollector::InvrotCollector( InvrotCollector const & other )
  : ReferenceCount(), invrots_(other.invrots_), owner_nodes_and_locations_(other.owner_nodes_and_locations_)
{}

InvrotCollector::~InvrotCollector(){}

InvrotCollectorOP
InvrotCollector::clone() const {
  return InvrotCollectorOP( new InvrotCollector( *this ) );
}

void
InvrotCollector::set_invrots_for_listnum(
  Size listnum,
  std::list<core::conformation::ResidueCOP> const & invrots,
  InvrotTreeNodeBaseCOP tree_node,
  Size location_in_node
)
{
  runtime_assert( listnum < invrots_.size() );
  //if( owner_nodes_and_locations_.find( tree_node ) != owner_nodes_and_locations_.end() ) utility_exit_with_message("Trying to add stuff for a node that has already been added.");

  //if( invrots_[listnum].size() != 0) utility_exit_with_message("Trying to overwrite rotamers in invrot collector");
  invrots_[listnum] = invrots;
  if( owner_nodes_and_locations_.find( tree_node ) != owner_nodes_and_locations_.end() ) owner_nodes_and_locations_.erase( tree_node );
  owner_nodes_and_locations_.insert( std::pair< InvrotTreeNodeBaseCOP, Size >( tree_node, location_in_node ) );
}



InvrotTreeNodeBase::InvrotTreeNodeBase(
  InvrotTreeNodeBaseCAP parent_node )
  : ReferenceCount(), parent_node_(parent_node ), location_in_parent_node_(1)
{}

InvrotTreeNodeBase::~InvrotTreeNodeBase(){}

}
}
}
