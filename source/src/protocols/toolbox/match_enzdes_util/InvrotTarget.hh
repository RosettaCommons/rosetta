// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/InvrotTarget.hh
/// @brief  .hh file for inverse rotamer target
/// @author Florian Richter, flosopher@gmail.com, mar 2012

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTarget_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_InvrotTarget_hh


//unit headers
#include <protocols/toolbox/match_enzdes_util/InvrotTarget.fwd.hh>

//package headers
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNodeBase.hh>

#ifdef PYROSETTA
	#include <core/conformation/Residue.hh>
#endif

#ifdef WIN32
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNode.hh>
#else
#include <protocols/toolbox/match_enzdes_util/InvrotTreeNode.fwd.hh>
#endif

// Utility headers
//#include <util

namespace protocols {
namespace toolbox {
namespace match_enzdes_util {


/// @brief the 'root' target against which the inverse rotamers are built
/// abstract base class to allow for invrots being built against
/// any sort of entity
class InvrotTarget : public InvrotTreeNodeBase {

public:

	InvrotTarget();

	~InvrotTarget();

	core::conformation::ResidueCOP
	target_res_for_geom_cst( core::Size geom_cst ) const;

	std::list< core::conformation::ResidueCOP >
	all_target_res() const;

	/// @brief generate constraints against a certain pose
	core::scoring::constraints::ConstraintCOP
	generate_constraints(
		core::pose::Pose const & pose,
		AllowedSeqposForGeomCstCOP geomcst_seqpos
	) const;

	/// @brief
	/// can initialize tree nodes according to
	/// an enzcst io
	/// note that this function presumes that
	/// representative_target_res_for_geom_cst_ and
	/// all_target_res_ have been set through calling the
	/// generate_representative_target_res_for_geom_cst()
	/// function as implemented by the child class
	bool
	initialize_tree_nodes_from_enzcst_io(
		EnzConstraintIOCOP enzcst_io );

	utility::vector1< std::list< core::conformation::ResidueCOP > >
	all_target_residues( InvrotTreeNodeBaseCAP child_node ) const;

	virtual
	void
	collect_all_inverse_rotamers(
		utility::vector1< InvrotCollectorOP > & invrot_collector
	) const;

protected:

	/// @brief this function figures out the coordinates
	/// (in Residue format) of the target residue for
	/// each geomcst.
	/// for a ligand, it's simply that ligand
	virtual
	void
	generate_representative_target_res_for_geom_cst( Size const num_geom_cst ) = 0;

	void
	set_all_target_res( std::list< core::conformation::ResidueCOP > const & all_target_res );

	void
	set_representative_target_res_for_geom_cst(
		utility::vector1< core::conformation::ResidueCOP > const & representative_res );

private:

	utility::vector1< core::conformation::ResidueCOP > representative_target_res_for_geom_cst_;
	std::list< core::conformation::ResidueCOP > all_target_res_;

	utility::vector1< InvrotTreeNodeOP > next_nodes_;

};


class SingleResidueInvrotTarget : public InvrotTarget {

public:

	SingleResidueInvrotTarget(
		utility::vector1< core::conformation::ResidueCOP > const & single_res
	);

	~SingleResidueInvrotTarget();

private:

	void
	generate_representative_target_res_for_geom_cst( Size const num_geom_cst );


};


}
}
}

#endif
