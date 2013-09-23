// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/AtomTreeMinimizer.cc
/// @brief  High-level atom tree minimizer class
/// @author Phil Bradley


// Unit headers
#include <core/optimization/AtomTreeMinimizer.hh>

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>
#include <core/optimization/AtomTreeMultifunc.hh>

// Project headers
#include <core/kinematics/MoveMap.fwd.hh>
// AUTO-REMOVED #include <core/id/AtomID_Mask.hh>
// AUTO-REMOVED #include <core/id/AtomID_Map.Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>


#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


using namespace ObjexxFCL::format;

namespace core {
namespace optimization {

AtomTreeMinimizer::AtomTreeMinimizer()
{}

AtomTreeMinimizer::~AtomTreeMinimizer() {}

///////////////////////////////////////////////////////////////////////////////
Real
AtomTreeMinimizer::run(
	pose::Pose & pose,
	kinematics::MoveMap const & move_map,
	scoring::ScoreFunction const & scorefxn,
	MinimizerOptions const & options
) /*const*/
{
	if ( options.deriv_check() ) {
		deriv_check_result_ = new NumericalDerivCheckResult;
		deriv_check_result_->send_to_stdout( options.deriv_check_to_stdout() );
	}

	bool const use_nblist( options.use_nblist() );

	//fpd fail if this is called on a symmetric pose
	runtime_assert( !core::pose::symmetry::is_symmetric( pose ) );

	// it's important that the structure be scored prior to nblist setup
	Real const start_score( scorefxn( pose ) );

	//std::cout << "start_score:" << start_score << std::endl;
	//pose.energies().show( std::cout );

	// setup the map of the degrees of freedom
	MinimizerMap min_map;
	min_map.setup( pose, move_map );

	// if we are using the nblist, set it up
	if ( use_nblist ) {
		// setup a mask of the moving dofs
		pose.energies().set_use_nblist( pose, min_map.domain_map(), options.nblist_auto_update() );
	}

	scorefxn.setup_for_minimizing( pose, min_map );

	// setup the function that we will pass to the low-level minimizer
	AtomTreeMultifunc f( pose, min_map, scorefxn,
		options.deriv_check(), options.deriv_check_verbose() );

	if ( deriv_check_result_ ) f.set_deriv_check_result( deriv_check_result_ );

	// starting position -- "dofs" = Degrees Of Freedom
	Multivec dofs( min_map.nangles() );
	min_map.copy_dofs_from_pose( pose, dofs );

	Real const start_func( f( dofs ) );

	//std::cout << "start_func: " << start_func <<  std::endl;
	//pose.energies().show( std::cout );

	// now do the optimization with the low-level minimizer function
	Minimizer minimizer( f, options );
	minimizer.run( dofs );

	Real const end_func( f( dofs ) );

	//std::cout << "end_func: " << end_func << std::endl;
	//pose.energies().show( std::cout );

	// turn off nblist
	if ( use_nblist ) pose.energies().reset_nblist();

	// if we were doing rigid-body minimization, fold the rotation and
	// translation offsets into the jump transforms
	//
	// also sets rb dofs to 0.0, so in principle func value should be the same
	//
	min_map.reset_jump_rb_deltas( pose, dofs );

	// rescore
	Real const end_score( scorefxn( pose ) );

	//std::cout << "end_score:" << std::endl;
	//pose.energies().show( std::cout );

	// we may not really need all these extra function evaluations
	// good for diagnostics though

	static basic::Tracer core_optimize( "core.optimize",  basic::t_debug);
	core_optimize << "AtomTreeMinimizer::run: nangles= " << min_map.nangles() <<
		" start_score: " << F(12,3,start_score) <<
		" start_func: "  << F(12,3,start_func ) <<
		" end_score: "   << F(12,3,end_score  ) <<
		" end_func: "    << F(12,3,end_func   ) << std::endl;


	return end_score;
}

NumericalDerivCheckResultOP
AtomTreeMinimizer::deriv_check_result() const
{
	return deriv_check_result_;
}

} // namespace optimization
} // namespace core
