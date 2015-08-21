// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/CartesianMinimizer.cc
/// @brief  High-level Cartesian minimizer class
/// @author Frank DiMaio


// Unit headers
#include <core/optimization/CartesianMinimizer.hh>

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>
#include <core/optimization/CartesianMultifunc.hh>

// Project headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>


#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

#include <core/kinematics/Jump.hh>
#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/DOF_Node.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


using namespace ObjexxFCL::format;

namespace core {
namespace optimization {

static thread_local basic::Tracer TR( "core.optimization.CartesianMinimizer" );

CartesianMinimizer::CartesianMinimizer()
{}

CartesianMinimizer::~CartesianMinimizer() {}

///////////////////////////////////////////////////////////////////////////////
Real
CartesianMinimizer::run(
	pose::Pose & pose,
	kinematics::MoveMap const & move_map,
	scoring::ScoreFunction const & scorefxn,
	MinimizerOptions const & options
) /*const*/
{
	if ( ! scorefxn.ready_for_nonideal_scoring() ) {
		utility_exit_with_message( "Scorefunction not set up for nonideal/Cartesian scoring" );
	}
	if ( options.min_type() != "lbfgs_armijo_nonmonotone" && options.min_type() != "lbfgs_armijo" && options.min_type() != "linmin" ) {
		TR.Warning << "WARNING: Use of the 'lbfgs_armijo_nonmonotone' minimizer with Cartesian minimization is recommended " <<
			"for better runtime performance. (Using '" << options.min_type() << "' minimizer instead.)" << std::endl;
	}

	if ( options.deriv_check() ) {
		deriv_check_result_ = NumericalDerivCheckResultOP( new NumericalDerivCheckResult );
		deriv_check_result_->send_to_stdout( options.deriv_check_to_stdout() );
	}

	bool const use_nblist( options.use_nblist() );

	//fpd this works with symm!
	//runtime_assert( !core::pose::symmetry::is_symmetric( pose ) );

	// it's important that the structure be scored prior to nblist setup
	Real const start_score( scorefxn( pose ) );

	// setup the map of the degrees of freedom
	CartesianMinimizerMap min_map;
	min_map.setup( pose, move_map );

	// if we are using the nblist, set it up
	if ( use_nblist ) {
		// setup a mask of the moving dofs
		pose.energies().set_use_nblist( pose, min_map.domain_map(), options.nblist_auto_update() );
	}

	scorefxn.setup_for_minimizing( pose, min_map );

	// setup the function that we will pass to the low-level minimizer
	CartesianMultifunc f( pose, min_map, scorefxn,
		options.deriv_check(), options.deriv_check_verbose() );

	if ( deriv_check_result_ ) f.set_deriv_check_result( deriv_check_result_ );

	// starting position -- "dofs" = Degrees Of Freedom
	Multivec dofs( min_map.ndofs() );
	min_map.copy_dofs_from_pose( pose, dofs );

	Real const start_func( f( dofs ) );

	// now do the optimization with the low-level minimizer function
	Minimizer minimizer( f, options );
	minimizer.run( dofs );

	Real const end_func( f( dofs ) );

	// turn off nblist
	if ( use_nblist ) pose.energies().reset_nblist();

	// rescore
	Real const end_score( scorefxn( pose ) );

	// we may not really need all these extra function evaluations
	// good for diagnostics though
	basic::Tracer core_optimize( "core.optimize", basic::t_debug );
	core_optimize << "CartesianMinimizer::run: ndofs= " << min_map.ndofs() <<
		" start_score: " << F(12,3,start_score) <<
		" start_func: "  << F(12,3,start_func ) <<
		" end_score: "   << F(12,3,end_score  ) <<
		" end_func: "    << F(12,3,end_func   ) << std::endl;

	return end_score;
}

NumericalDerivCheckResultOP
CartesianMinimizer::deriv_check_result() const
{
	return deriv_check_result_;
}

} // namespace optimization
} // namespace core
