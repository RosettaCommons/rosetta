// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/normalmode/NormalModeMinimizer.cc
/// @brief  High-level atom tree minimizer class
/// @author Phil Bradley


// Unit headers
#include <protocols/normalmode/NormalMode.hh>
#include <protocols/normalmode/NormalModeMinimizer.hh>
#include <protocols/normalmode/NormalModeMinimizerCreator.hh>
#include <protocols/normalmode/NormalModeMultiFunc.hh>

// Package headers
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>

// Project headers
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <iostream>
#include <protocols/rosetta_scripts/util.hh>

#include <boost/foreach.hpp>
#include <basic/Tracer.hh>

#include <ObjexxFCL/format.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


using namespace ObjexxFCL::format;
using namespace core;
using namespace kinematics;
using namespace optimization;
using namespace scoring;

namespace protocols {
namespace normalmode {

std::string
NormalModeMinimizerCreator::keyname() const
{
	return NormalModeMinimizerCreator::mover_name();
}

protocols::moves::MoverOP
NormalModeMinimizerCreator::create_mover() const {
	return protocols::moves::MoverOP( new NormalModeMinimizer );
}

std::string
NormalModeMinimizerCreator::mover_name()
{
	return "NormalModeMinimizer";
}

////////////

static THREAD_LOCAL basic::Tracer nma_tr( "protocols.normalmode", basic::t_debug );


NormalModeMinimizer::NormalModeMinimizer() {
	options_ = MinimizerOptionsOP( new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.01, true, false, false ) );
	options_->max_iter(200);
}

NormalModeMinimizer::~NormalModeMinimizer() {}

void
NormalModeMinimizer::apply( pose::Pose & pose ) {
	if ( ! movemap_ ) movemap_ = core::kinematics::MoveMapOP( new MoveMap );
	if ( ! scorefxn_ ) scorefxn_ = get_score_function(); // get a default (INITIALIZED!) ScoreFunction

	if ( options_->deriv_check() ) {
		deriv_check_result_ = NumericalDerivCheckResultOP( new NumericalDerivCheckResult );
		deriv_check_result_->send_to_stdout( options_->deriv_check_to_stdout() );
	}

	bool const use_nblist( options_->use_nblist() );

	//fpd << fix this!
	runtime_assert( !core::pose::symmetry::is_symmetric( pose ) );

	Real const start_score( (*scorefxn_)( pose ) );
	MinimizerMap min_map;
	min_map.setup( pose, *movemap_ );

	if ( use_nblist ) {
		pose.energies().set_use_nblist( pose, min_map.domain_map(), options_->nblist_auto_update() );
	}

	scorefxn_->setup_for_minimizing( pose, min_map );

	// Normalmode setup: solve for starting pose
	// We can put in more options for NormalMode here, e.g. distance cutoffs
	// but let's just use default options right now
	protocols::normalmode::NormalMode normalmode;
	normalmode.torsion( true );
	normalmode.solve( pose );

	// setup the function that we will pass to the low-level minimizer
	protocols::normalmode::NormalModeMultifunc
		f( pose, min_map, *scorefxn_, normalmode,
		false, //omega
		options_->deriv_check(), options_->deriv_check_verbose() );

	// Reset modes to be minimized in MultiFunc when defined by user
	if ( modes_using_.size() > 0 ) {
		f.set_modes( modes_using_ );
	}

	if ( deriv_check_result_ ) deriv_check_local( pose, f );

	// starting position -- "dofs" = Degrees Of Freedom
	Multivec dofs( min_map.nangles() );
	min_map.copy_dofs_from_pose( pose, dofs );

	scorefxn_->show( pose );

	// Convert into vars
	Multivec vars( f.dofs_to_vars( dofs ) );

	Real const start_func( f( vars ) );

	nma_tr << "start_func: " << start_func <<  std::endl;

	// now do the optimization with the low-level minimizer function
	Minimizer minimizer( f, *options_ );
	minimizer.run( vars );

	Real const end_func( f( vars ) );

	nma_tr << "end_func: " << end_func << std::endl;

	// turn off nblist
	if ( use_nblist ) pose.energies().reset_nblist();

	dofs = f.vars_to_dofs( vars );
	min_map.reset_jump_rb_deltas( pose, dofs );

	Real const end_score( (*scorefxn_)( pose ) );

	nma_tr << "NormalModeMinimizer::run: nangles= " << min_map.nangles() <<
		" start_score: " << F(12,3,start_score) <<
		" start_func: "  << F(12,3,start_func ) <<
		" end_score: "   << F(12,3,end_score  ) <<
		" end_func: "    << F(12,3,end_func   ) << std::endl;
}

std::string
NormalModeMinimizer::get_name() const {
	return NormalModeMinimizerCreator::mover_name();
}


void NormalModeMinimizer::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	Filters_map const &,
	protocols::moves::Movers_map const &,
	Pose const & pose )
{
	if ( ! movemap_ ) movemap_ = core::kinematics::MoveMapOP( new MoveMap );

	std::string const scorefxn_name( tag->getOption< std::string >( "scorefxn", "score12" ) );
	scorefxn_ = data.get_ptr<ScoreFunction>( "scorefxns", scorefxn_name );

	// high-level movemap
	bool const chi( tag->getOption< bool >( "chi" ) ), bb( tag->getOption< bool >( "bb" ) );
	movemap_->set_chi( chi );
	movemap_->set_bb( bb );
	if ( tag->hasOption("jump") ) {
		if ( ! movemap_ ) movemap_ = core::kinematics::MoveMapOP( new MoveMap );
		utility::vector1<std::string> jumps = utility::string_split( tag->getOption<std::string>( "jump" ), ',' );

		// string 'ALL' makes all jumps movable
		if ( jumps.size() == 1 && (jumps[1] == "ALL" || jumps[1] == "All" || jumps[1] == "all") ) {
			movemap_->set_jump( true );
		} else if ( tag->getOption< core::Size > ( "jump" ) == 0 ) {
			movemap_->set_jump( false );
		} else {
			BOOST_FOREACH ( std::string jump, jumps ) {
				Size const value = std::atoi( jump.c_str() );
				movemap_->set_jump( value, true );
			}
		}
	}

	// minimizer
	options_->max_iter( tag->getOption< int >( "max_iter", 200 ) );
	options_->min_type( tag->getOption< std::string >( "type", "lbfgs_armijo_nonmonotone" ) );
	options_->minimize_tolerance( tag->getOption< core::Real >( "tolerance", 0.01 ) );

	// nma-min specific
	dampen_ = tag->getOption< core::Real >( "dampen", 0.05 );

	// fine-grained movemap control
	protocols::rosetta_scripts::parse_movemap( tag, pose, movemap_, data, false );
}


NumericalDerivCheckResultOP
NormalModeMinimizer::deriv_check_result() const {
	return deriv_check_result_;
}

void
NormalModeMinimizer::deriv_check_local( pose::Pose const &pose,
	protocols::normalmode::NormalModeMultifunc f ) const {

	pose::Pose pose_work( pose );
	Real const dvar( 1.0e-2 );

	Multivec const vars0( f.nvar(), 0.0 );
	Multivec analytical_deriv( f.nvar(), 0.0 );

	// First calculate function and derivative for starting pose
	f.dfunc( vars0, analytical_deriv );

	// Then enumerate over degree of freedom
	Multivec vars;
	nma_tr << LJ(3,"DOF") << " "
		<< LJ(6,"dval") << " "
		<< LJ(10,"g_ana") << " "
		<< LJ(10,"g_num") << " "
		<< LJ(10,"dg") << " "
		<< LJ(10,"dg/g_num") << std::endl;

	for ( Size i_var = 1; i_var <= f.nvar(); ++i_var ) {

		Real f1, f2;
		// +dvar
		vars = vars0; vars[i_var] += dvar;
		f1 = f( vars );

		// -dvar
		vars = vars0; vars[i_var] -= dvar;
		f2 = f( vars );

		Real numeric_deriv = 0.5*(f1 - f2) /dvar;

		// Report
		nma_tr << I(3,i_var) << " "
			<< E(6,1,dvar) << " "
			<< F(10,5,analytical_deriv[i_var]) << " "
			<< F(10,5,numeric_deriv ) << " "
			<< F(10,5,analytical_deriv[i_var] - numeric_deriv  ) << " "
			<< F(10,5,(analytical_deriv[i_var] - numeric_deriv) / numeric_deriv   ) << std::endl;
	}

}

} // namespace normalmode
} // namespace protocols
