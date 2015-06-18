// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/MgMinimizer.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/magnesium/MgMinimizer.hh>
#include <protocols/magnesium/util.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/ScoreFunction.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.magnesium.MgMinimizer" );

using namespace core;
using utility::vector1;

namespace protocols {
namespace magnesium {

	//Constructor
	MgMinimizer::MgMinimizer():
		minimize_scorefxn_( get_mg_scorefxn() /* can be replaced by user*/ ),
		mg_coord_cst_dist_( 0.2 )
	{}

	//Destructor
	MgMinimizer::~MgMinimizer()
	{}

	/////////////////////////////////////////////////////
	void
	MgMinimizer::apply( core::pose::Pose & pose ) {

		using namespace core::optimization;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::scoring::func;
		using namespace core::scoring::constraints;
		using namespace protocols::magnesium;

		vector1< Size > pose_mg_res = mg_res_;
		if ( pose_mg_res.size() == 0 ) pose_mg_res = get_mg_res( pose );

		//	update_mg_hoh_fold_tree( pose ); // trust input fold tree...
		kinematics::MoveMap const mm = get_mg_hoh_minimize_move_map( pose, pose_mg_res );

		ScoreFunctionOP minimize_scorefxn_working = minimize_scorefxn_->clone();
		bool const coord_cst_in_scorefxn = minimize_scorefxn_working->has_nonzero_weight( coordinate_constraint );
		if ( !coord_cst_in_scorefxn ) minimize_scorefxn_working->set_weight( coordinate_constraint, 1.0 );

		vector1< ConstraintCOP > coord_csts;
		if ( mg_coord_cst_dist_ > 0.0 ) {
			for ( Size n = 1; n <= pose_mg_res.size(); n++ ) {
				coord_csts.push_back( pose.add_constraint( ConstraintOP(
																																new CoordinateConstraint( AtomID( 1, pose_mg_res[ n ] ),
																																													AtomID( 1, 1 ),
																																													pose.residue( pose_mg_res[ n ] ).xyz( 1 ),
																																													FuncOP( new HarmonicFunc( 0.0, mg_coord_cst_dist_ ) ) ) ) ) );
			}
		}

		std::string const min_type = "dfpmin_armijo_nonmonotone";
		Real const min_tolerance = 0.000025;
		bool const use_nblist( true );
		MinimizerOptionsOP minimizer_options_( new MinimizerOptions( min_type, min_tolerance, use_nblist, false, false ) );
		minimizer_options_->nblist_auto_update( true );
		AtomTreeMinimizerOP atom_tree_minimizer_( new AtomTreeMinimizer );
		atom_tree_minimizer_->run( pose, mm, *minimize_scorefxn_working, *minimizer_options_ );

		if ( coord_csts.size() > 0 ) {
			for ( Size n = 1; n <= coord_csts.size(); n++ ) pose.remove_constraint( coord_csts[ n ] );
		}
		if ( !coord_cst_in_scorefxn ) minimize_scorefxn_working->set_weight( coordinate_constraint, 0.0 );

		(*minimize_scorefxn_working)( pose );

	}

	//////////////////////////////////////////
	// @brief move all jumps associated with Mg(2+)
	core::kinematics::MoveMap
	MgMinimizer::get_mg_hoh_minimize_move_map( core::pose::Pose const & pose,
																						 utility::vector1< Size > const & mg_res ) const {
		using namespace core::kinematics;
		FoldTree const & f( pose.fold_tree() );
		MoveMap mm;
		mm.set_bb( false );
		mm.set_chi( false );
		mm.set_jump( false );
		for ( Size n = 1; n <= f.num_jump(); n++ ) {
			if ( mg_res.has_value( Size( f.upstream_jump_residue( n ) ) ) ||
					 mg_res.has_value( Size( f.downstream_jump_residue( n ) ) ) ) {
				TR.Debug << "Going to minimize jump: " << n << std::endl;
				mm.set_jump( n, true );
			}
		}
		return mm;
	}


} //magnesium
} //protocols
