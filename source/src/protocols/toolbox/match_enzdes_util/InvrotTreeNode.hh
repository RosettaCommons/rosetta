// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTreeNode.hh
/// @brief  .hh file for inverse rotamer tree node
/// @author Florian Richter, flosopher@gmail.com, mar 2012

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTreeNode_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTreeNode_hh


//unit headers
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNode.fwd.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>

//project headers
#include <core/id/AtomID.fwd.hh>

// Utility headers
//#include <util

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


class InvrotTreeNode : public InvrotTreeNodeBase {

	// what's better for this? list or vector1?
	typedef std::pair< std::list<core::conformation::ResidueCOP >, utility::vector1< InvrotTreeNodeBaseOP > > invrots_node_ptrs_pair;

	typedef core::Size Size;

public:

	InvrotTreeNode( InvrotTreeNodeBaseCAP parent );

	~InvrotTreeNode();

	/// @brief generate invrots according to
	/// an enzdes cstfile
	/// returns true if it was possible to
	/// build invrots, i.e. geometry correctly
	/// defined in cstfile and no clashes with
	/// any parent inverse rotamers
	bool
	initialize_from_enzcst_io(
		core::conformation::Residue const & target_residue,
		EnzConstraintIOCOP enzcst_io,
		Size invrot_geomcst,
		core::pose::PoseCOP pose = NULL
	);

	bool
	initialize_from_enzcst_io_and_invrots(
		std::list< core::conformation::ResidueCOP > const & all_invrots,
		EnzConstraintIOCOP enzcst_io,
		Size invrot_geomcst,
		core::pose::PoseCOP pose = NULL
	);

	/// @brief
	/// 1. for each invrots/node pointer pair:
	///     a.generate an ambigous constraint for the inverse rots
	///     for each node pointer to child nodes (only if there are child nodes, of course)
	///         b. call the generate_constraints function on the pointers to child nodes
	///     c. package the ambig constraint generated for the inverse rots
	///     and the constraints that come from the child nodes into one multicst
	///     caveat: if the child nodes send up a multi constraint and not an ambiguous
	///     constraint, 'unpackage' those multicst and generate a new multicst
	///     containing this nodes ambigcst and the member_csts from the
	///     child nodes' multicsts
	///
	/// 2. if there is more than one invrots/node pointer pair,
	///    throw all the generated multicsts into one ambig cst
	/// 3. return either the multicst (one pair) or the ambigcst
	/// 4. hope this works
	core::scoring::constraints::ConstraintCOP
	generate_constraints(
		core::pose::Pose const & pose,
		AllowedSeqposForGeomCstCOP geomcst_seqpos
	) const;


	/// @brief this function returns the AtomID for an atom
	/// in the pose that's supposed to stay fixed during
	/// folding, i.e. the neighbor atom of the first target
	/// needed to generate the right backbone_stub constraints
	core::id::AtomID
	get_fixed_pt( core::pose::Pose const & pose ) const;

	utility::vector1< std::list< core::conformation::ResidueCOP > >
	all_target_residues( InvrotTreeNodeBaseCAP child_node ) const;


	void
	remove_invrots_clashing_with_parent_res(
		std::list< core::conformation::ResidueCOP > & invrots,
		bool covalent
	) const;

	void
	collect_all_inverse_rotamers(
		utility::vector1< InvrotCollectorOP > & invrot_collectors
	) const;

private:

	Size geom_cst_; //which geometric constraint the invrots belong to

	utility::vector1< invrots_node_ptrs_pair > invrots_and_next_nodes_;

	bool generate_invrot_csts_; //if this node should actually generate csts according to the invrots stored in it. initialized to true, but gets set to fals if inverse rotamer trees for partial matches are constructed
};

}
}
}

#endif
