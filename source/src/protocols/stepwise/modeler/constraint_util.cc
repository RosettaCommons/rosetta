// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/constraint_util.cc
/// @brief
/// @detailed
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
#include <basic/Tracer.hh>

using namespace core;

static basic::Tracer TR( "protocols.stepwise.modeler.constraint_util" );

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

		//		std::map< Size, Size > slice_map;
		//		for (Size i = 1; i <= slice_res.size(); i++) slice_map[ slice_res[ i ] ] = i;
		utility::vector1< Size > mapping( full_pose.total_residue(), 0);
		for (Size i = 1; i <= slice_res.size(); i++) mapping[ slice_res[ i ] ] = i;
		SequenceMappingOP smap = new SequenceMapping( mapping );

		for ( Size n = 1; n <= csts.size(); n++ ) {

			ConstraintCOP const & cst( csts[n] );
			ConstraintOP cst_new = cst->remapped_clone( full_pose, pose, smap );
			if ( cst_new ) {
				cst_set_new->add_constraint( cst_new );
				//				std::cout << "HEY CONSTRAINTS!!! "
				//									<< cst_new->atom(1).rsd() << " " << pose.residue_type( cst_new->atom(1).rsd() ).atom_name( cst_new->atom(1).atomno() )
				//									<< " to "
				//									<< cst_new->atom(2).rsd() << " " << pose.residue_type( cst_new->atom(2).rsd() ).atom_name( cst_new->atom(2).atomno() ) << std::endl;
			}

			// currently only defined for pairwise distance constraints,
			//  and coordinate constraints
			//			if ( cst->score_type() == atom_pair_constraint)  {

// 				Size const i = cst->atom( 1 ).rsd();
// 				Size const j = cst->atom( 2 ).rsd();
// 				//			Size const dist( shortest_path_in_fold_tree.dist( i , j ) );
// 				//			if ( dist  > separation_cutoff ) continue;

// 				if ( slice_map.find( i ) == slice_map.end()  ) continue;
// 				if ( slice_map.find( j ) == slice_map.end()  ) continue;

// 				std::cout << "CST MAP: " << i << " " << slice_map[ i] << "          " << j << " " << slice_map[ j ] << std::endl;

// 				std::string const & atom_name1 = full_pose.residue_type( i ).atom_name( cst->atom(1).atomno() );
// 				std::string const & atom_name2 = full_pose.residue_type( j ).atom_name( cst->atom(2).atomno() );

// 				AtomID atom1_new( named_atom_id_to_atom_id( NamedAtomID( atom_name1, slice_map[ i ] ), pose );
// 				AtomID atom2_new( named_atom_id_to_atom_id( NamedAtomID( atom_name2, slice_map[ j ] ), pose );

// 				ConstraintOP cst_new = new AtomPairConstraint( atom1_new, atom2_new,
// 																											 cst->get_func().clone() /*is this defined?*/, cst->score_type() );

//			if ( cst_new ) cst_set_new->add_constraint( cst_new );


				//			} else if ( cst->score_type() == coordinate_constraint)  {




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

		if ( pose.constraint_set()->has_constraints() )	{
			if ( scorefxn->has_zero_weight( atom_pair_constraint ) ||
					 scorefxn->has_zero_weight( coordinate_constraint ) ) {
				utility_exit_with_message( "Since we want constraints, need to use a scorefunction with non-zero atom_pair_constraint and coordinate_constraint weight");
			}
		}

	}



} //modeler
} //stepwise
} //protocols
