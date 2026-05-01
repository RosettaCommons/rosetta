// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/optimization/ParametricAtomTreeMultifunc.hh>
#include <core/optimization/parametric_minimize_util.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/Pose.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers

// External headers
#include <ObjexxFCL/format.hh>


using namespace ObjexxFCL::format;

namespace core {
namespace optimization {

static basic::Tracer TR( "core.optimization.AtomTreeMinimizer" );

AtomTreeMinimizer::AtomTreeMinimizer() = default;

AtomTreeMinimizer::~AtomTreeMinimizer() = default;

///////////////////////////////////////////////////////////////////////////////
Real
AtomTreeMinimizer::run(
	pose::Pose & pose,
	kinematics::MoveMap const & move_map,
	scoring::ScoreFunction const & scorefxn,
	MinimizerOptions const & options
) /*const*/
{
	check_setup( pose, move_map, scorefxn, options);

	if ( options.deriv_check() ) {
		deriv_check_result_ = utility::pointer::make_shared< NumericalDerivCheckResult >();
		deriv_check_result_->send_to_stdout( options.deriv_check_to_stdout() );
	}

	bool const use_nblist( options.use_nblist() );

	//fpd fail if this is called on a symmetric pose
	runtime_assert( !core::pose::symmetry::is_symmetric( pose ) );

	// it's important that the structure be scored prior to nblist setup
	Real const start_score( scorefxn( pose ) );

	if ( TR.Debug.visible() ) {
		TR.Debug << "start_score:" << start_score << std::endl;
		pose.energies().show( TR.Trace );
	}

	// Check whether parametric DOFs should be included
	bool const use_parametric = move_map.get_parametric() && pose.conformation().n_parameters_sets() > 0;
	utility::vector1< ParametricDOFInfo > parametric_dofs;
	std::set< Size > parametric_residues;
	if ( use_parametric ) {
		enumerate_parametric_dofs( pose, parametric_dofs );
		parametric_residues = get_parametric_residues( pose );
		if ( TR.Debug.visible() ) {
			TR.Debug << "Parametric minimization: " << parametric_dofs.size() << " parametric DOFs, "
				<< parametric_residues.size() << " residues under parametric control." << std::endl;
		}
	}

	// Setup the map of the degrees of freedom.
	// If parametric minimization is active, we need a MoveMap that excludes backbone
	// torsions of parametric residues (to prevent redundant DOF control).
	kinematics::MoveMap effective_move_map( move_map );
	if ( use_parametric ) {
		for ( Size resid : parametric_residues ) {
			effective_move_map.set_bb( resid, false );
		}
	}

	MinimizerMap min_map;
	min_map.setup( pose, effective_move_map );

	// if we are using the nblist, set it up
	if ( use_nblist ) {
		// setup a mask of the moving dofs
		pose.energies().set_use_nblist( pose, min_map.domain_map(), options.nblist_auto_update() );
	}

	scorefxn.setup_for_minimizing( pose, min_map );

	// Setup the multifunc — use parametric version if parametric DOFs are present
	Multivec dofs;
	Real start_func, end_func;
	if ( use_parametric && !parametric_dofs.empty() ) {
		ParametricAtomTreeMultifunc f( pose, min_map, scorefxn, parametric_dofs,
			options.deriv_check(), options.deriv_check_verbose() );
		dofs.resize( f.total_dofs() );
		min_map.copy_dofs_from_pose( pose, dofs );
		for ( Size p = 1; p <= parametric_dofs.size(); ++p ) {
			dofs[ min_map.nangles() + p ] = get_parametric_dof_value( pose, parametric_dofs[p] );
		}

		start_func = f( dofs );
		if ( TR.Trace.visible() && !use_nblist ) {
			pose.energies().show( TR.Trace );
		}

		Minimizer minimizer( f, options );
		minimizer.run( dofs );
		end_func = f( dofs );
	
		if ( TR.Debug.visible() && ! use_nblist ) {
			TR.Debug << "end_func: " << end_func << std::endl;
			pose.energies().show( TR.Trace );
		}

		// turn off nblist
		if ( use_nblist ) pose.energies().reset_nblist();
	
		// if we were doing rigid-body minimization, fold the rotation and
		// translation offsets into the jump transforms
		//
		// also sets rb dofs to 0.0, so in principle func value should be the same
		//
		min_map.reset_jump_rb_deltas( pose, dofs );

	} else {
		AtomTreeMultifunc f( pose, min_map, scorefxn,
			options.deriv_check(), options.deriv_check_verbose() );

		if ( deriv_check_result_ ) f.set_deriv_check_result( deriv_check_result_ );

		dofs.resize( min_map.nangles() );
		min_map.copy_dofs_from_pose( pose, dofs );

		start_func = f( dofs );
		if ( TR.Trace.visible() && !use_nblist ) {
			pose.energies().show( TR.Trace );
		}

		Minimizer minimizer( f, options );
		minimizer.run( dofs );
		end_func = f( dofs );
	
		if ( TR.Debug.visible() && ! use_nblist ) {
			TR.Debug << "end_func: " << end_func << std::endl;
			pose.energies().show( TR.Trace );
		}

		// turn off nblist
		if ( use_nblist ) pose.energies().reset_nblist();

		// if we were doing rigid-body minimization, fold the rotation and
		// translation offsets into the jump transforms
		//
		// also sets rb dofs to 0.0, so in principle func value should be the same
		//
		min_map.reset_jump_rb_deltas( pose, dofs );

	}

	// rescore
	Real const end_score( scorefxn( pose ) );

	if ( TR.Debug.visible() && ! use_nblist ) {
		TR.Debug << "end_score:" << std::endl;
		pose.energies().show( TR.Trace );
	}

	// we may not really need all these extra function evaluations
	// good for diagnostics though

	scorefxn.finalize_after_minimizing( pose );

	basic::TracerImpl core_optimize( "core.optimize", basic::t_debug );
	core_optimize << "AtomTreeMinimizer::run: nangles= " << min_map.nangles() <<
		" start_score: " << F(12,3,start_score) <<
		" start_func: "  << F(12,3,start_func ) <<
		" end_score: "   << F(12,3,end_score  ) <<
		" end_func: "    << F(12,3,end_func   ) << std::endl;


	return end_score;
}

/// @brief Do consistency checks for minimizer setup.
void
AtomTreeMinimizer::check_setup(pose::Pose const &,
	kinematics::MoveMap const & move_map,
	scoring::ScoreFunction const & scorefxn,
	MinimizerOptions const & options) const {

	// Do we have any nonideal bond length/bond angle settings in the movemap?
	bool nonideal( move_map.get( core::id::D ) || move_map.get( core::id::THETA ) );
	if ( ! nonideal ) {
		for ( auto iter( move_map.dof_id_begin() ); iter != move_map.dof_id_end(); ++iter ) {
			if ( iter->second && ( iter->first.type() == core::id::D || iter->first.type() == core::id::THETA ) ) {
				nonideal = true;
				break;
			}
		}
	}

	if ( nonideal && ! scorefxn.ready_for_nonideal_scoring() ) {
		// It would be good to have this be a hard error, but several protocols are currently
		// almost-but-not-quite ready_for_nonideal_scoring().
		//TR.Error << "** Non-ideal Movemap used in minimization:" << std::endl;
		//move_map.show( TR.Error );
		//TR.Error << "** Scorefunction not ready for nonideal scoring:" << std::endl;
		//scorefxn.show( TR.Error );
		//TR.Error << std::endl;
		//utility_exit_with_message( "Scorefunction not set up for nonideal scoring" );
		TR.Warning << "***************************************************************" << std::endl;
		TR.Warning << "** WARNING: Non-ideal minimization used with a ScoreFunction **" << std::endl;
		TR.Warning << "**    which isn't set up for non-ideal minimization          **" << std::endl;
		TR.Warning << "***************************************************************" << std::endl;
	}

	if ( nonideal && options.min_type() != "lbfgs_armijo_nonmonotone" ) {
		TR.Warning << "Use of the 'lbfgs_armijo_nonmonotone' minimizer with nonideal minimization is recommended " <<
			"for better runtime performance. (Using '" << options.min_type() << "' minimizer instead.)" << std::endl;
	}
}

NumericalDerivCheckResultOP
AtomTreeMinimizer::deriv_check_result() const
{
	return deriv_check_result_;
}

} // namespace optimization
} // namespace core
