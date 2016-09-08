// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/constraint_util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/constraint_util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/id/SequenceMapping.hh>
#include <basic/Tracer.hh>

using namespace core;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.constraint_util" );

///////////////////////////////////////////////////////////////////////////////
// these functions are pretty old -- will be souped up when I revive constraints
// in stepwise monte carlo.
///////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace modeler {

///////////////////////////////////////////////////////////////////////////////
// Currently only handles atom pair constraints.
///////////////////////////////////////////////////////////////////////////////
core::scoring::constraints::ConstraintSetOP
constraint_set_slice( core::scoring::constraints::ConstraintSetOP & cst_set,
	utility::vector1< core::Size > const & slice_res,
	pose::Pose const & pose,
	pose::Pose const & full_pose )
{
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::id;

	ConstraintSetOP cst_set_new( new scoring::constraints::ConstraintSet );

	ConstraintCOPs csts( cst_set->get_all_constraints() );

	//  std::map< Size, Size > slice_map;
	//  for (Size i = 1; i <= slice_res.size(); i++) slice_map[ slice_res[ i ] ] = i;
	utility::vector1< Size > mapping( full_pose.size(), 0);
	for ( Size i = 1; i <= slice_res.size(); i++ ) mapping[ slice_res[ i ] ] = i;
	SequenceMappingOP smap( new SequenceMapping( mapping ) );

	for ( Size n = 1; n <= csts.size(); n++ ) {

		ConstraintCOP const & cst( csts[n] );
		ConstraintOP cst_new = cst->remapped_clone( full_pose, pose, smap );
		if ( cst_new ) {
			cst_set_new->add_constraint( cst_new );
		}
	}

	std::cout << "NUM CONSTRAINTS " << cst_set_new->get_all_constraints().size() << " out of " <<
		csts.size() << std::endl;

	return cst_set_new;
}

///////////////////////////////////////////////////////////////////////
void
check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose::Pose const & pose,
	core::scoring::ScoreFunctionOP & scorefxn ){

	using namespace core::scoring;

	if ( pose.constraint_set()->has_constraints() ) {
		if ( scorefxn->has_zero_weight( atom_pair_constraint ) ||
				scorefxn->has_zero_weight( coordinate_constraint ) ) {
			utility_exit_with_message( "Since we want constraints, need to use a scorefunction with non-zero atom_pair_constraint and coordinate_constraint weight");
		}
	}
}


} //modeler
} //stepwise
} //protocols
