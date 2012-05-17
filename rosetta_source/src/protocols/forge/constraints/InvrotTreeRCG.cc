// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/constraints/InvrotTreeRCG.cc
///
/// @brief
/// @author Florian Richter, floric@u.washington.edu, may 2012

#include <protocols/forge/constraints/InvrotTreeRCG.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/constraints/Constraint.hh>

#include <protocols/toolbox/match_enzdes_util/AllowedSeqposForGeomCst.hh>
#include <protocols/toolbox/match_enzdes_util/InvrotTree.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.forge.constraints.InvrotTreeRCG" );

namespace protocols{
namespace forge{
namespace constraints{


InvrotTreeRCG::InvrotTreeRCG(
	toolbox::match_enzdes_util::InvrotTreeOP invrot_tree,
	toolbox::match_enzdes_util::AllowedSeqposForGeomCstCOP geomcst_seqpos
	) : invrot_tree_(invrot_tree), geomcst_seqpos_(geomcst_seqpos)
{}

InvrotTreeRCG::~InvrotTreeRCG(){}

void
InvrotTreeRCG::generate_remodel_constraints(
		core::pose::Pose const & pose )
{
	//for now
	runtime_assert( invrot_tree_->num_target_states() == 1 );
	Size target_state = 1;

	invrot_tree_->generate_inverse_rotamer_constraints( pose, geomcst_seqpos_ );

	this->add_constraint( invrot_tree_->get_constraint_for_target_state( target_state ) );

}

} //namespace remodel
} //namespace forge
} //namespace protocols
