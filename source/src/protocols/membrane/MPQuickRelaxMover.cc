// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief      Does a quick relax of a membrane protein
/// @details Uses SmallMover and ShearMover with adjustable maximum dihedral
///    angle changes, then repacking and a single round of minimization
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_MPQuickRelaxMover_cc
#define INCLUDED_protocols_membrane_MPQuickRelaxMover_cc

// Unit Headers
#include <protocols/membrane/MPQuickRelaxMover.hh>
#include <protocols/membrane/MPQuickRelaxMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/kinematics/MoveMap.hh>
#include <core/import_pose/import_pose.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <protocols/membrane/util.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>


// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.MPQuickRelaxMover" );

namespace protocols {
namespace membrane {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Defaults: dih angle_max = 1, nmoves = nres, movemap = bb and all chi
MPQuickRelaxMover::MPQuickRelaxMover() :
	protocols::moves::Mover(),
	movemap_( new core::kinematics::MoveMap() )
{
	set_defaults();
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor
MPQuickRelaxMover::MPQuickRelaxMover( core::Real angle_max, std::string nmoves ) :
	protocols::moves::Mover(),
	movemap_( new core::kinematics::MoveMap() )
{
	set_defaults();
	register_options();
	init_from_cmd();

	angle_max_ = angle_max;
	moves_ = nmoves;
}

/// @brief Custom constructor
MPQuickRelaxMover::MPQuickRelaxMover(
	core::Real angle_max,
	std::string nmoves,
	core::kinematics::MoveMapOP movemap ) :
	protocols::moves::Mover(),
	movemap_( new core::kinematics::MoveMap() )
{
	set_defaults();
	register_options();
	init_from_cmd();

	angle_max_ = angle_max;
	moves_ = nmoves;
	movemap_ = movemap;
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
MPQuickRelaxMover::MPQuickRelaxMover( MPQuickRelaxMover const & src ) : protocols::moves::Mover( src ),
	native_( src.native_ ),
	angle_max_( src.angle_max_ ),
	moves_( src.moves_ ),
	nmoves_( src.nmoves_ ),
	movemap_( src.movemap_ ),
	sfxn_( src.sfxn_ ),
	cst_file_( src.cst_file_ ),
	cst_weight_( src.cst_weight_ ),
	addmem_( src.addmem_ ),
	mem_from_topo_ ( src.mem_from_topo_ ),
	opt_mem_( src.opt_mem_ )
{}

/// @brief Assignment Operator
MPQuickRelaxMover & MPQuickRelaxMover::operator = ( MPQuickRelaxMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new MPQuickRelaxMover( *this ) );
}

/// @brief Destructor
MPQuickRelaxMover::~MPQuickRelaxMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPQuickRelaxMover::clone() const {
	return ( protocols::moves::MoverOP( new MPQuickRelaxMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPQuickRelaxMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MPQuickRelaxMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MPQuickRelaxMover::parse_my_tag(
	utility::tag::TagCOP /*tag*/,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// TODO: implement this

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
MPQuickRelaxMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MPQuickRelaxMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MPQuickRelaxMoverCreator::keyname() const {
	return MPQuickRelaxMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MPQuickRelaxMoverCreator::mover_name() {
	return "MPQuickRelaxMover";
}

/// @brief Get the name of this Mover (MPQuickRelaxMover)
std::string
MPQuickRelaxMover::get_name() const {
	return "MPQuickRelaxMover";
}

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Run a quick relax on a membrane protein
/// @details Requires AddMembraneMover to be run beforehand
void MPQuickRelaxMover::apply( core::pose::Pose & pose ) {

	using namespace protocols::membrane;
	using namespace protocols::simple_moves;
	using namespace protocols::moves;
	using namespace core::pack::task;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::kinematics;

	TR << "Running MPQuickRelax protocol..." << std::endl;

	// add membrane to pose
	if ( addmem_ == true ) {
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );
	}

	// starting foldtree
	TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// get number of residues and number of moves
	core::Size nres( nres_protein( pose ) );
	if ( moves_ == "nres" ) {
		nmoves_ = nres;
	} else {
		nmoves_ = utility::string2Size( moves_ );
	}

	// find membrane position around the protein
	if ( mem_from_topo_ ) {
		MembranePositionFromTopologyMoverOP mempos( new MembranePositionFromTopologyMover() );
		mempos->anchor_at_res1( false );
		mempos->apply( pose );
	}

	// optimize membrane position using smooth scorefunction
	if ( opt_mem_ ) {
		OptimizeMembranePositionMoverOP optmem( new OptimizeMembranePositionMover() );
		optmem->apply( pose );
	}

	// creating constraint set
	if ( cst_file_.size() > 0 && cst_weight_ > 0 ) {

		ConstraintSetOP cst( ConstraintIO::get_instance()->read_constraints( cst_file_, ConstraintSetOP( new ConstraintSet() ), pose ) );
		pose.constraint_set( cst );

		// set constraints in scorefunction
		sfxn_->set_weight( atom_pair_constraint, cst_weight_ );
		sfxn_->set_weight( angle_constraint, cst_weight_ );
		sfxn_->set_weight( dihedral_constraint, cst_weight_ );
	}

	// starting position: Shake up the protein
	TR << "Choosing starting position: shake up the protein - Waka waka ..." << std::endl;

	// set small and shearmover
	core::Real kT = 1.0;
	SmallMoverOP small( new SmallMover( movemap_, kT, nmoves_ ) );
	small->angle_max( angle_max_ );
	small->apply( pose );

	ShearMoverOP shear( new ShearMover( movemap_, kT, nmoves_ ) );
	shear->angle_max( angle_max_ );
	shear->apply( pose );

	// create MC object
	MonteCarloOP mc( new MonteCarlo( pose, *sfxn_, 1.0 ) );

	// initialize AtomTreeMinimizer
	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
	core::optimization::AtomTreeMinimizer atm;

	//    fa_rep min_tol cst_wts
	// repeat 5
	// ramp_repack_min 0.02  0.01     1.0
	// ramp_repack_min 0.250 0.01     0.5
	// ramp_repack_min 0.550 0.01     0.0
	// ramp_repack_min 1     0.00001  0.0

	// sfxn_->set_weight( fa_rep, 0.5 );

	// run small and shearmover again
	TR << "SmallMover and ShearMover - Waka waka ..." << std::endl;
	small->apply( pose );
	shear->apply( pose );

	// do this for a certain number of iterations
	// since both the packer and especially minimization takes a while, just do one for now
	core::Size breakpoint( 0 );

	// score the pose
	core::Real tot_score = ( *sfxn_ )( pose );
	core::Real fa_rep = pose.energies().total_energies()[ core::scoring::fa_rep ];
	core::Real fa_atr = pose.energies().total_energies()[ core::scoring::fa_atr ];

	// if fa_rep exceeds threshold (number of residues of the protein, as empirically determined)
	// or score is greater than 0, keep continuing
	while ( tot_score > 0 || fa_rep > nres_protein( pose ) + 100 ) {

		pose.energies().show_total_headers( TR );
		TR << std::endl;
		pose.energies().show_totals( TR );
		TR << std::endl;
		TR << "tot: " << tot_score << " fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;

		++breakpoint;

		// ramp repulsive up and down
		if ( breakpoint % 2 == 1 ) {
			sfxn_->set_weight( core::scoring::fa_rep, 1.0);
		} else {
			sfxn_->set_weight( core::scoring::fa_rep, 0.1);
		}

		// if breakpoint, fail_retry
		//  if ( breakpoint >= 10 && breakpoint <= 20 ) {
		//   TR << "Counter > 10, score: " << tot_score << ", fa_rep: " << fa_rep << ", FAIL_RETRY..." << std::endl;
		//   set_last_move_status(protocols::moves::FAIL_RETRY);
		//  }
		//  else if ( breakpoint > 20 ){
		//   TR << "Counter > 20, score: " << tot_score << ", fa_rep: " << fa_rep << ", FAIL_DO_NOT_RETRY..." << std::endl;
		//   set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
		//   return;
		//  }

		if ( breakpoint == 5 ) {
			TR << "Counter == 6, score: " << tot_score << ", fa_rep: " << fa_rep << ", FAIL_DO_NOT_RETRY..." << std::endl;
			set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
			return;
		}

		// packing
		TR << "Packing rotamers..." << std::endl;
		PackerTaskOP repack = TaskFactory::create_packer_task( pose );
		repack->restrict_to_repacking();
		core::pack::pack_rotamers( pose, *sfxn_, repack );

		// minimize
		TR << "Minimizing..." << std::endl;
		atm.run( pose, *movemap_, *sfxn_, min_opts );

		// evaluate Boltzmann
		mc->boltzmann( pose );
		TR << "accepted? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;

		// set pose to lowest scoring pose from MC simulation
		pose = mc->lowest_score_pose();

		// score the pose
		tot_score = ( *sfxn_ )( pose );
		fa_rep = pose.energies().total_energies()[ core::scoring::fa_rep ];
		pose.energies().show_total_headers( TR );
		TR << std::endl;
		pose.energies().show_totals( TR );
		TR << std::endl;
		TR << "tot: " << tot_score << "fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;
		TR << "Iteration: " << breakpoint << " score: " << tot_score << " fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;

		pose.dump_pdb("relax_test_" + utility::to_string( breakpoint ) + ".pdb");

	} // number of iterations for search

	// reset the weight in the scorefunction to what it was before (for mpsmooth)
	sfxn_->set_weight( core::scoring::fa_rep, 0.44);
	tot_score = ( *sfxn_ )( pose );
	fa_rep = pose.energies().total_energies()[ core::scoring::fa_rep ];
	TR << "Final score: " << tot_score << " final fa_rep: " << fa_rep << std::endl;

	// superimpose poses with native
	SuperimposeMoverOP super( new SuperimposeMover( *native_, 1, nres, 1, nres, true ) );
	super->apply( pose );

	// get job for adding rmsd to scorefile
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	// calculate and store the rmsd in the score file
	job->add_string_real_pair("rms", core::scoring::bb_rmsd( pose, *native_ ));

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

}// apply

////////////////////////////////////////////////////////////////////////////////

/// @brief Run AddMembraneMover again?
/// @details If you want to keep your anchor point for MEM, then pick no
void MPQuickRelaxMover::add_membrane_again( bool yesno ) {
	addmem_ = yesno;
}

/// @brief Run MembranePositionFromTopology again?
/// @details Will change the starting membrane position
void MPQuickRelaxMover::membrane_from_topology( bool yesno ) {
	mem_from_topo_ = yesno;
}

/// @brief Optimize membrane position before relax?
void MPQuickRelaxMover::optimize_membrane( bool yesno ) {
	opt_mem_ = yesno;
}

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
void MPQuickRelaxMover::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::mp::quickrelax::angle_max );
	option.add_relevant( OptionKeys::mp::quickrelax::nmoves );
	option.add_relevant( OptionKeys::constraints::cst_file );
	option.add_relevant( OptionKeys::constraints::cst_weight );

}

////////////////////////////////////////////////////////////////////////////////

/// @brief Set default values
void MPQuickRelaxMover::set_defaults() {

	// maximum change in dihedral angles
	angle_max_ = 1.0;

	// number of moves
	moves_ = "nres";

	// Movemap
	movemap_->set_bb( true );
	movemap_->set_chi( true );

	// create scorefunction
	sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );

	// default constraint weight of 1
	cst_weight_ = 1.0;

	// run AddMembraneMover again?
	addmem_ = true;

	// starting MembranePositionFromTopology?
	mem_from_topo_ = false;

	// optimize membrane before relax?
	opt_mem_ = false;

}// set_defaults

////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize from commandline
void MPQuickRelaxMover::init_from_cmd() {

	using namespace basic::options;

	// read native and attach membrane to it
	if ( option[ OptionKeys::in::file::native ].user() ) {
		native_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( *native_ );
	}

	if ( option[ OptionKeys::mp::quickrelax::angle_max ].user() ) {
		angle_max_ = option[ OptionKeys::mp::quickrelax::angle_max ]();
	}

	if ( option[ OptionKeys::mp::quickrelax::nmoves ].user() ) {
		if ( option[ OptionKeys::mp::quickrelax::nmoves ]() == "nres" ) {
			moves_ = utility::to_string( nres_protein( *native_ ) );
		} else {
			moves_ = option[ OptionKeys::mp::quickrelax::nmoves ]();
		}
	}

	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		cst_file_ = option[ OptionKeys::constraints::cst_file ]()[1];
	}

	if ( option[ OptionKeys::constraints::cst_weight ].user() ) {
		cst_weight_ = option[ OptionKeys::constraints::cst_weight ]();
	}

}// init from cmdline


} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_MPQuickRelaxMover_cc
