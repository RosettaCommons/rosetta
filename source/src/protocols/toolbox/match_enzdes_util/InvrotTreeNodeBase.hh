// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.hh
/// @brief  Forward declaration for inverse rotamer tree node base
/// @author Florian Richter, flosopher@gmail.com, mar 2012

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTreeNodeBase_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTreeNodeBase_hh

/// unit headeers
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.fwd.hh>

/// package headers
#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.fwd.hh>

/// project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/types.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.fwd.hh>

//c++ heades
#include <list>
#include <map>

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


///@ brief helper class for collecting all different definitions of invrots from a tree
class InvrotCollector : public utility::pointer::ReferenceCount {

public:
	InvrotCollector( Size num_residue_lists );

	InvrotCollector( InvrotCollector const & other );

	virtual ~InvrotCollector();

	InvrotCollectorOP
	clone() const;

	void
	set_invrots_for_listnum(
		Size listnum,
		std::list<core::conformation::ResidueCOP> const & invrots,
		InvrotTreeNodeBaseCOP tree_node,
		Size location_in_node
	);

	std::vector< std::list<core::conformation::ResidueCOP> > const &
	invrots() const {
		return invrots_; }

	std::map< InvrotTreeNodeBaseCOP, Size > const &
	owner_nodes_and_locations() const {
		return owner_nodes_and_locations_; }

private:

	//we're using standard vectors bc the targets put into 0th element
	std::vector< std::list<core::conformation::ResidueCOP> > invrots_;
	std::map< InvrotTreeNodeBaseCOP, Size > owner_nodes_and_locations_;
};



/// @brief abstract base class for an inverse rotamer tree node.
/// Exists so that the target node (InvrotTreeTarget) and the different
/// geom cst nodes (InvrotTreeNode) have a common base class.
/// This is necessary so that a node can point at its parent
/// node in a tree without having to worry about whether that
/// is a target or a regular node
class InvrotTreeNodeBase : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< InvrotTreeNodeBase >
{
public:

  InvrotTreeNodeBase(
    InvrotTreeNodeBaseCAP parent_node
  );

  virtual ~InvrotTreeNodeBase();

	/// self pointers
	inline InvrotTreeNodeBaseCOP get_self_ptr() const { return shared_from_this(); }
	inline InvrotTreeNodeBaseOP get_self_ptr() { return shared_from_this(); }
	inline InvrotTreeNodeBaseCAP get_self_weak_ptr() const { return InvrotTreeNodeBaseCAP( shared_from_this() ); }
	inline InvrotTreeNodeBaseAP get_self_weak_ptr() { return InvrotTreeNodeBaseAP( shared_from_this() ); }

  InvrotTreeNodeBaseCAP
  parent_node() const {
    return parent_node_; }

	void
	set_location_in_parent_node( Size location ) {
		location_in_parent_node_ = location; }

	Size location_in_parent_node() const {
		return location_in_parent_node_; }

  /// @brief nodes need to be able to generate constraints
  virtual
  core::scoring::constraints::ConstraintCOP
  generate_constraints(
    core::pose::Pose const & pose,
    AllowedSeqposForGeomCstCOP geomcst_seqpos
  ) const = 0;


  /// @brief this function traverses up the tree
  /// and adds the target residue for every geomcst
  /// in the branch toward this node
  /// used for clash checking, i.e. in case where
  /// the interaction is lig<-geomcst1-<geomcst2-<geomcst3,
  /// we don't want geomcst2 rots that clash with
  /// the ligand, and we don't wand geomcst3 rots
  /// that clash with the ligand or geomcst1 res
  /// the child node argument usually represents the node
  /// that is asking for the target residues, i.e. the node
  /// that's calling this function
  virtual
  utility::vector1< std::list< core::conformation::ResidueCOP > >
  all_target_residues( InvrotTreeNodeBaseCAP child_node ) const = 0;


  /// @brief convenience funtion to get all
  /// inverse rotamers in the tree
  /// puts all the inverse rotamers associated
  /// with this node into vector, and should
  /// call this function on daughter nodes
	/// uses std:vector because the targets will be put into
	/// 0th element
	/// needs vector of vector bc there can be different definitions
	/// of the tree
  virtual
  void
  collect_all_inverse_rotamers(
    utility::vector1< InvrotCollectorOP > & invrot_collectors
  ) const = 0;


private:

  InvrotTreeNodeBaseCAP parent_node_; //pointer to parent node
	Size location_in_parent_node_;
};


}
}
}

#endif
