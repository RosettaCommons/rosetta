// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author

#include <core/util/cryst_util.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDB_Info.hh>

#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/LineMinimizer.hh>

#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/AtomTreeMultifunc.hh>
#include <core/optimization/AtomTreeMinimizer.hh>

#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/optimization/CartesianMinimizer.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh>
#include <core/optimization/symmetry/SymAtomTreeMultifunc.hh>

#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/symmetry/util.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <iostream>

namespace core {
namespace util {

static basic::Tracer TS( "core.scoring.cryst.util" );

core::Real getMLweight(core::scoring::ScoreFunction & scorefxn,  core::pose::Pose & pose ) {
	// if no movemap is specified assume everything can move
	core::kinematics::MoveMap mm;
	mm.set_bb  ( true );
	mm.set_chi ( true );
	mm.set_jump( true );

	return getMLweight( scorefxn, pose, mm );
}

core::Real getMLweight( core::scoring::ScoreFunction & scorefxn, core::pose::Pose &pose_orig , core::kinematics::MoveMap &move_map) {
	using namespace core::optimization;
	using namespace core::optimization::symmetry;

	// create two scorefunctions, one experimental only, one rosetta only
	core::scoring::ScoreFunctionOP rosetta_scorefxn = scorefxn.clone(); //core::scoring::get_score_function();
	core::scoring::ScoreFunctionOP xtal_scorefxn = new core::scoring::ScoreFunction();

	rosetta_scorefxn->set_weight( core::scoring::xtal_ml, 0.0 );
	xtal_scorefxn->set_weight( core::scoring::xtal_ml, 1.0 );

	core::Real grad2_ros=0, grad2_xtal=0;

	core::pose::Pose pose = pose_orig; // copy the pose since we will perturb it

	// symmetrize?
	if (core::pose::symmetry::is_symmetric(pose)) {
		core::conformation::symmetry::SymmetricConformation const & symm_conf (
				dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation() ) );
		core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

		// symmetrize scorefunct & movemap
		rosetta_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *rosetta_scorefxn );
		xtal_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *xtal_scorefxn );
		core::pose::symmetry::make_symmetric_movemap( pose, move_map );

		// compute gradients using both scorefunctions
		SymMinimizerMap min_map( pose, move_map, symm_info );
		Multivec vars( min_map.nangles() ), dEros_dvars, dExtal_dvars;
		min_map.copy_dofs_from_pose( pose, vars );

		(*rosetta_scorefxn)(pose);  // score pose first
		rosetta_scorefxn->setup_for_minimizing( pose, min_map );
		SymAtomTreeMultifunc f_ros( pose, min_map, *rosetta_scorefxn, false, false );
		f_ros.dfunc( vars, dEros_dvars );

		(*xtal_scorefxn)(pose);  // score pose first
		xtal_scorefxn->setup_for_minimizing( pose, min_map );
		SymAtomTreeMultifunc f_xtal( pose, min_map, *xtal_scorefxn, false, false );
		f_xtal.dfunc( vars, dExtal_dvars );

		// sum
		for (int i=1; i<=(int)vars.size(); ++i) {
			grad2_xtal += dExtal_dvars[ i ]*dExtal_dvars[ i ];
			grad2_ros += dEros_dvars[ i ]*dEros_dvars[ i ];
		}

		core::Real w_xtal1 = (grad2_xtal != 0) ?  1 * sqrt( grad2_ros / grad2_xtal ) : 1;
		TS << " Gradient ratio = " << w_xtal1 << " = sqrt( " << grad2_ros << " / " << grad2_xtal << " )" << std::endl;

		return w_xtal1;

		// now take a step in the grad direction and recompute
		// 		{
		// 			(*rosetta_scorefxn)(pose);  // score pose first
		// 			rosetta_scorefxn->setup_for_minimizing( pose, min_map );
		// 			BrentLineMinimization test_brent( f_ros, vars.size() );
		// 			core::Real newval = test_brent( vars, dEros_dvars );
		//
		// 			f_ros.dfunc( vars, dEros_dvars );
		//
		// 			min_map.copy_dofs_to_pose( pose, vars );
		// 			(*xtal_scorefxn)(pose);  // score pose first
		// 			xtal_scorefxn->setup_for_minimizing( pose, min_map );
		// 			f_xtal.dfunc( vars, dExtal_dvars );
		//
		// 			// sum
		// 			grad2_xtal = grad2_ros = 0;
		// 			for (int i=1; i<=(int)vars.size(); ++i) {
		// 				grad2_xtal += dExtal_dvars[ i ]*dExtal_dvars[ i ];
		// 				grad2_ros += dEros_dvars[ i ]*dEros_dvars[ i ];
		// 			}
		// 		}
		//
		// 		core::Real w_xtal2 = (grad2_xtal != 0) ?  1 * sqrt( grad2_ros / grad2_xtal ) : 1;
		// 		std::cerr << " w2 = " << w_xtal2 << " = sqrt( " << grad2_ros << " / " << grad2_xtal << " )" << std::endl;
		//
		// 		return 0.5*(w_xtal1+w_xtal2);
	} else {
		// compute gradients using both scorefunctions
		MinimizerMap min_map;
		min_map.setup( pose, move_map );
		Multivec vars( min_map.nangles() ), dEros_dvars, dExtal_dvars;
		min_map.copy_dofs_from_pose( pose, vars );

		(*rosetta_scorefxn)(pose);  // score pose first
		rosetta_scorefxn->setup_for_minimizing( pose, min_map );
		AtomTreeMultifunc f_ros( pose, min_map, *rosetta_scorefxn, false, false );
		f_ros.dfunc( vars, dEros_dvars );

		(*xtal_scorefxn)(pose);  // score pose first
		xtal_scorefxn->setup_for_minimizing( pose, min_map );
		AtomTreeMultifunc f_xtal( pose, min_map, *xtal_scorefxn, false, false );
		f_xtal.dfunc( vars, dExtal_dvars );

		// sum
		for (int i=1; i<=(int)vars.size(); ++i) {
			grad2_xtal += dExtal_dvars[ i ]*dExtal_dvars[ i ];
			grad2_ros += dEros_dvars[ i ]*dEros_dvars[ i ];
		}

		core::Real w_xtal;
		if (grad2_xtal != 0) {
			w_xtal = 1 * sqrt( grad2_ros / grad2_xtal );
		} else {
			w_xtal = 1;
		}

		TS << " Gradient ratio = " << w_xtal << " = sqrt( " << grad2_ros << " / " << grad2_xtal << " )" << std::endl;

		return w_xtal;
	}
}


/////////////////////////////
/////////////////////////////


core::Real getMLweight_cart( core::scoring::ScoreFunction & scorefxn, core::pose::Pose & pose ) {
	// if no movemap is specified assume everything can move
	core::kinematics::MoveMap mm;
	mm.set_bb  ( true );
	mm.set_chi ( true );
	mm.set_jump( true );

	return getMLweight_cart( scorefxn, pose, mm );
}

core::Real getMLweight_cart( core::scoring::ScoreFunction & scorefxn, core::pose::Pose &pose , core::kinematics::MoveMap &move_map) {
	using namespace core::optimization;
	using namespace core::optimization::symmetry;

	// create two scorefunctions, one experimental only, one rosetta only
	core::scoring::ScoreFunctionOP rosetta_scorefxn  = scorefxn.clone();//core::scoring::get_score_function();
	core::scoring::ScoreFunctionOP xtal_scorefxn = new core::scoring::ScoreFunction();

	rosetta_scorefxn->set_weight( core::scoring::xtal_ml, 0.0 );
	xtal_scorefxn->set_weight( core::scoring::xtal_ml, 1.0 );

	// dont use cart bonded term to fit weights
	//rosetta_scorefxn->set_weight( core::scoring::cart_bonded, 0.0 );

	core::Real grad2_ros=0, grad2_xtal=0;

	// symmetrize?
	if (core::pose::symmetry::is_symmetric(pose)) {
		core::conformation::symmetry::SymmetricConformation const & symm_conf (
				dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation() ) );
		core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

		// symmetrize scorefunct & movemap
		rosetta_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *rosetta_scorefxn );
		xtal_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *xtal_scorefxn );
		core::pose::symmetry::make_symmetric_movemap( pose, move_map );
	}

	// compute gradients using both scorefunctions
	CartesianMinimizerMap min_map;
	min_map.setup( pose, move_map );
	Multivec vars( min_map.ndofs() ), dEros_dvars, dExtal_dvars;
	min_map.copy_dofs_from_pose( pose, vars );

	(*rosetta_scorefxn)(pose);  // score pose first
	rosetta_scorefxn->setup_for_minimizing( pose, min_map );
	CartesianMultifunc f_ros( pose, min_map, *rosetta_scorefxn, false, false );
	f_ros.dfunc( vars, dEros_dvars );

	(*xtal_scorefxn)(pose);  // score pose first
	xtal_scorefxn->setup_for_minimizing( pose, min_map );
	CartesianMultifunc f_xtal( pose, min_map, *xtal_scorefxn, false, false );
	f_xtal.dfunc( vars, dExtal_dvars );

	// sum
	for (int i=1; i<=(int)vars.size(); ++i) {
		grad2_xtal += dExtal_dvars[ i ]*dExtal_dvars[ i ];
		grad2_ros += dEros_dvars[ i ]*dEros_dvars[ i ];
	}

	core::Real w_xtal1;
	if (grad2_xtal != 0)
		w_xtal1 = 1 * sqrt( grad2_ros / grad2_xtal );
	else
		w_xtal1 = 1;
	std::cerr << " w0 = " << w_xtal1 << " = sqrt( " << grad2_ros << " / " << grad2_xtal << " )" << std::endl;

	return 0.5*(w_xtal1+w_xtal1);
}

}
} // namespace core
