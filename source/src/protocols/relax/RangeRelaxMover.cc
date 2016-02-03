// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/relax/RangeRelaxMover.hh
/// @brief      Relaxes the protein by relaxing in ranges
/// @details Relaxes a protein by iteratively relaxing ranges of the protein;
///    No ramping required. Much faster than FastRelax and good for
///    large to very large proteins (tested up to 5250 residues);
///    For the membrane version, use MPRangeRelax which runs this
///    Mover in the underneath
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_relax_RangeRelaxMover_cc
#define INCLUDED_protocols_relax_RangeRelaxMover_cc

// Unit Headers
#include <protocols/relax/RangeRelaxMover.hh>
//#include <protocols/relax/RangeRelaxMoverCreator.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>

// Project Headers
#include <core/id/AtomID.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MonteCarlo.hh>
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
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <core/conformation/util.hh>
#include <numeric/random/random.hh>
#include <protocols/membrane/util.hh>
#include <core/scoring/rms_util.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.relax.RangeRelaxMover" );

namespace protocols {
namespace relax {

using namespace core;
using namespace core::pose;
using namespace protocols::moves;
using namespace protocols::membrane;
using namespace core::optimization;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Defaults: dih angle_max = 1, nmoves = nres, movemap = bb and all chi
RangeRelaxMover::RangeRelaxMover() :
	protocols::moves::Mover(),
	movemap_( new MoveMap() )
{
	set_defaults();
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor
RangeRelaxMover::RangeRelaxMover( core::Size center_resnumber ) :
	protocols::moves::Mover(),
	movemap_( new MoveMap() )
{
	set_defaults();
	register_options();
	init_from_cmd();

	center_resnumber_ = center_resnumber;
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
RangeRelaxMover::RangeRelaxMover( RangeRelaxMover const & src ) :
	protocols::moves::Mover( src ),
	native_( src.native_ ),
	center_resnumber_( src.center_resnumber_ ),
	angle_max_( src.angle_max_ ),
	moves_( src.moves_ ),
	nmoves_( src.nmoves_ ),
	movemap_( src.movemap_ ),
	sfxn_( src.sfxn_ ),
	sfxn_nocst_( src.sfxn_nocst_ ),
	cst_file_( src.cst_file_ ),
	cst_weight_( src.cst_weight_ ),
	cst_to_start_( src.cst_to_start_ ),
	cst_to_native_( src.cst_to_native_ ),
	addmem_( src.addmem_ ),
	optmem_( src.optmem_ ),
	repack_again_( src.repack_again_ ),
	cycles_ ( src.cycles_ ),
	min_cycles_ ( src.min_cycles_ ),
	spherical_wave_( src.spherical_wave_ ),
	idealize_( src.idealize_ )
{}

/// @brief Assignment Operator
RangeRelaxMover & RangeRelaxMover::operator = ( RangeRelaxMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new RangeRelaxMover( *this ) );
}

/// @brief Destructor
RangeRelaxMover::~RangeRelaxMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
//protocols::moves::MoverOP
//RangeRelaxMover::clone() const {
// return ( protocols::moves::MoverOP( new RangeRelaxMover( *this ) ) );
//}
//
// /// @brief Create a Fresh Instance of this Mover
//protocols::moves::MoverOP
//RangeRelaxMover::fresh_instance() const {
// return protocols::moves::MoverOP( new RangeRelaxMover() );
//}
//
// /// @brief Pase Rosetta Scripts Options for this Mover
//void
//RangeRelaxMover::parse_my_tag(
//        utility::tag::TagCOP /*tag*/,
//        basic::datacache::DataMap &,
//        protocols::filters::Filters_map const &,
//        protocols::moves::Movers_map const &,
//        core::pose::Pose const &
//        ) {
//
//  // TODO: implement this
//
//}
//
// /// @brief Create a new copy of this mover
//protocols::moves::MoverOP
//RangeRelaxMoverCreator::create_mover() const {
// return protocols::moves::MoverOP( new RangeRelaxMover() );
//}
//
// /// @brief Return the Name of this mover (as seen by Rscripts)
//std::string
//RangeRelaxMoverCreator::keyname() const {
// return RangeRelaxMoverCreator::mover_name();
//}
//
// /// @brief Mover name for Rosetta Scripts
//std::string
//RangeRelaxMoverCreator::mover_name() {
// return "RangeRelaxMover";
//}
//

/// @brief Get the name of this Mover (RangeRelaxMover)
std::string
RangeRelaxMover::get_name() const {
	return "RangeRelaxMover";
}

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Run RangeRelax
void RangeRelaxMover::apply( Pose & pose ) {

	using namespace basic::options;
	using namespace numeric::random;
	using namespace protocols::simple_moves;
	using namespace protocols::membrane;
	using namespace protocols::moves;
	using namespace core::pack::task;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace core::scoring::constraints;

	TR << "Running RangeRelax protocol..." << std::endl;

	FoldTree orig_ft = finalize_setup( pose );

	// switching to highres pose
	if ( pose.is_centroid() ) {
		TR << "Pose is centroid. Switching to full-atom." << std::endl;
		SwitchResidueTypeSetMoverOP highres( new SwitchResidueTypeSetMover() );
		highres->apply( pose );
	}

	// starting position: Shake up the protein
	TR << "Choosing starting position: shake up the protein - Waka waka ..." << std::endl;

	// if not constraining to starting coords, do small initial perturbation of the pose
	Real kT = option[ OptionKeys::relax::range::kT ]();
	SmallMoverOP small( new SmallMover( movemap_, kT, nmoves_ ) );
	ShearMoverOP shear( new ShearMover( movemap_, kT, nmoves_ ) );
	// if ( cst_to_start_ == false ) {

	small->angle_max( 0.05 );
	small->apply( pose );
	shear->angle_max( 0.05 );
	shear->apply( pose );
	// }

	// create MC object
	MonteCarloOP mc( new MonteCarlo( pose, *sfxn_, kT ) );

	// initialize AtomTreeMinimizer
	core::optimization::MinimizerOptions min_opts( "lbfgs_armijo_atol", 0.01, true );
	min_opts.max_iter( min_cycles_ );
	core::optimization::AtomTreeMinimizer atm;

	// run small and shearmover with defined angle max
	// from then on it's only repacking and minimization
	if ( cst_to_start_ == false ) {

		TR << "SmallMover and ShearMover - Waka waka ..." << std::endl;

		// get a random number up to angle_max, this should randomize the sampling
		// range a little more than it used to
		TR << "angle max: " << angle_max_ << std::endl;

		small->angle_max( angle_max_ );
		small->apply( pose );
		shear->angle_max( angle_max_ );
		shear->apply( pose );
	}

	// only go through a certain number of iterations to keep runtimes short
	for ( Size i = 1; i <= cycles_; ++i ) {

		// if total score < 0 and fa_rep < 1.3 * nres, stop iterations to keep speed;
		// for this, score the pose first without cst
		Real score_nocst = ( *sfxn_nocst_ )( pose );
		Real fa_rep = pose.energies().total_energies()[ scoring::fa_rep ];
		if ( score_nocst < 0 && fa_rep <= 1.3 * nres_protein( pose ) ) {
			TR << "Total score and fa_rep already good, skipping iteration " << i << std::endl;
			continue;
		}
		( *sfxn_ )( pose );

		////////////////////////////////////////////////////////////////////////

		// do a repack run, depending on which algorithm to use: spherical
		// or window repack along the sequence, starting from a center residue
		if ( spherical_wave_ == true ) {
			repack_spherical_range( pose, mc );
		} else {
			repack_sequence_window( pose, mc );
		}

		////////////////////////////////////////////////////////////////////////

		// packing again
		if ( repack_again_ == true ) {
			repack_all( pose, mc );
		}

		////////////////////////////////////////////////////////////////////////

		// minimize all
		TR << "Minimizing round " << i << std::endl;
		atm.run( pose, *movemap_, *sfxn_, min_opts );

		// evaluate Boltzmann
		mc->boltzmann( pose );
		TR << "accepted minimization? " << mc->mc_accepted() << " pose energy: " << pose.energies().total_energy() << std::endl;

		TR << "iteration " << i << std::endl;
		print_score( pose );

	} // number of iterations for search

	// idealize pose
	if ( idealize_ == true ) {
		idealize_pose( pose, min_opts, atm );
	}

	// if membrane protein, optimize membrane position (MEM is already at root)
	// TODO: make MPRangeRelaxMover inherit from RangeRelax!
	if ( pose.conformation().is_membrane() && optmem_ == true ) {
		OptimizeMembranePositionMoverOP optmem( new OptimizeMembranePositionMover() );
		optmem->apply( pose );
	}

	// score the pose and print the scores
	TR << "Final score: " << std::endl;
	print_score( pose );

	// superimpose poses with native
	core::Size nres( nres_protein( pose ) );
	SuperimposeMoverOP super( new SuperimposeMover( *native_, 1, nres, 1, nres, true ) );
	super->apply( pose );

	// get job for adding rmsd to scorefile
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

	// calculate and store the rmsd in the score file
	job->add_string_real_pair("rms", core::scoring::bb_rmsd( pose, *native_ ));

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	pose.fold_tree().show( TR );

}// apply

////////////////////////////////////////////////////////////////////////////////

/// @brief Run AddMembraneMover again?
/// @details If you want to keep your anchor point for MEM, then pick no
void RangeRelaxMover::add_membrane_again( bool yesno ) {
	addmem_ = yesno;
}

/// @brief Optimize membrane
void RangeRelaxMover::optimize_membrane( bool yesno ) {
	optmem_ = yesno;
}

/// @brief Idealize pose after run?
/// @details Might lead to better decoy but takes longer
void RangeRelaxMover::idealize( bool yesno ) {
	idealize_ = yesno;
}

/// @brief Set maximum dihedral angle perturbation
void RangeRelaxMover::set_angle_max( core::Real angle_max ) {
	angle_max_ = angle_max;
}

/// @brief Set number of moves, can be "nres" or a number
void RangeRelaxMover::set_nmoves( std::string nmoves ) {
	moves_ = nmoves;
}

/// @brief Use membrane scorefunction
void RangeRelaxMover::set_scorefunction( ScoreFunctionOP sfxn ) {
	sfxn_ = sfxn;
	sfxn_nocst_ = sfxn->clone();
}

/// @brief Set native
void RangeRelaxMover::set_native( PoseOP pose ) {
	native_ = pose;
}

////////////////////////////////////////////////////////////////////////////////

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
void RangeRelaxMover::register_options() {

	using namespace basic::options;

	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::constraints::cst_file );
	option.add_relevant( OptionKeys::constraints::cst_weight );
	option.add_relevant( OptionKeys::constraints::cst_fa_file );
	option.add_relevant( OptionKeys::constraints::cst_fa_weight );
	option.add_relevant( OptionKeys::relax::constrain_relax_to_native_coords );
	option.add_relevant( OptionKeys::relax::constrain_relax_to_start_coords );
	option.add_relevant( OptionKeys::relax::range::angle_max );
	option.add_relevant( OptionKeys::relax::range::nmoves );
	option.add_relevant( OptionKeys::relax::range::spherical_wave );
	option.add_relevant( OptionKeys::relax::range::repack_again );
	option.add_relevant( OptionKeys::relax::range::cycles );
	option.add_relevant( OptionKeys::relax::range::min_cycles );
	option.add_relevant( OptionKeys::relax::range::idealize );
	option.add_relevant( OptionKeys::mp::transform::optimize_embedding );

	option.add_relevant( OptionKeys::relax::range::kT );

}

////////////////////////////////////////////////////////////////////////////////

/// @brief Set default values
void RangeRelaxMover::set_defaults() {

	// native
	native_ = 0;

	// center residue number
	center_resnumber_ = 0;

	// maximum change in dihedral angles
	angle_max_ = 0.1;

	// number of moves
	moves_ = "nres";

	// Movemap
	movemap_->set_bb( true );
	movemap_->set_chi( true );
	movemap_->set_jump( true );

	// create scorefunction
	sfxn_ = get_score_function();

	// scorefunction without constraints
	sfxn_nocst_ = get_score_function();

	// constraint file
	cst_file_ = "";

	// default constraint weight of 1
	cst_weight_ = 1.0;

	// constrain to starting coordinates?
	cst_to_start_ = false;

	// constraint to native?
	cst_to_native_ = false;

	// run AddMembraneMover again?
	addmem_ = false;

	// additional round of repack?
	repack_again_ = false;

	// cycles for repacking and minimization
	cycles_ = 3;

	// maximum number of minimization cycles, 2k is the default for lbfgs, but for
	// some membrane proteins this wasn't enough (however it didn't seem to
	// compromise the outcome)
	min_cycles_ = 2000;

	// relax in spherical wave pattern
	spherical_wave_ = false;

	// idealize after run? increases sampling range but scores typically get worse
	idealize_ = false;

	// optimize membrane?
	optmem_ = true;

}// set_defaults

////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize from commandline
void RangeRelaxMover::init_from_cmd() {

	using namespace basic::options;

	// read native and attach membrane to it
	if ( option[ OptionKeys::in::file::native ].user() ) {
		native_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);

		if ( addmem_ == true ) {
			AddMembraneMoverOP addmem( new AddMembraneMover() );
			addmem->apply( *native_ );
		}
	}

	if ( option[ OptionKeys::constraints::cst_file ].user() || option[ OptionKeys::constraints::cst_weight ].user() ) {
		utility_exit_with_message( "You are trying to read in centroid constraints but want to do high-res refinement. Convert your restraints to full-atom mode and use the options -constraints:cst_fa_file and -constraints:cst_fa_weight instead!" );
	}

	if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
		cst_file_ = option[ OptionKeys::constraints::cst_fa_file ]()[1];

		if ( !option[ OptionKeys::constraints::cst_fa_weight ].user() ) {
			TR << "No cst_fa_weight given. Using default option of " << cst_weight_ << std::endl;
		}
	}

	if ( option[ OptionKeys::constraints::cst_fa_weight ].user() ) {
		cst_weight_ = option[ OptionKeys::constraints::cst_fa_weight ]();

		if ( !option[ OptionKeys::constraints::cst_fa_file ].user() ) {
			utility_exit_with_message( "No constraint file given. Use option -constraints:cst_fa_file unless you are trying to constraint to the native or the starting structure. In this case, remove the -cst_fa_weight flag altogether. Quitting." );
		}
	}

	if ( option[ OptionKeys::relax::constrain_relax_to_start_coords ].user() ) {
		cst_to_start_ = option[ OptionKeys::relax::constrain_relax_to_start_coords ]();

		if ( cst_file_.size() != 0 ) {
			utility_exit_with_message( "Can't constrain to starting coordinates AND use additional constraints. It's complicated. Sorry." );
		}
	}

	if ( option[ OptionKeys::relax::constrain_relax_to_native_coords ].user() ) {
		cst_to_native_ = option[ OptionKeys::relax::constrain_relax_to_native_coords ]();

		if ( cst_to_native_ == true && native_ == 0 ) {
			utility_exit_with_message( "You want to constrain the coordinates to the native structure without a native given. Please provide a native structure with -in:file:native. Quitting." );
		}
		if ( cst_file_.size() != 0 ) {
			utility_exit_with_message( "Can't constrain to the native structure AND use additional constraints. It's complicated. I would suggest running with constraints first, then select the best-scoring decoy and relax again with constraining to starting coordinates." );
		}
	}

	if ( option[ OptionKeys::relax::range::angle_max ].user() ) {
		angle_max_ = option[ OptionKeys::relax::range::angle_max ]();
	}

	if ( option[ OptionKeys::relax::range::nmoves ].user() ) {
		moves_ = option[ OptionKeys::relax::range::nmoves ]();
	}

	if ( option[ OptionKeys::relax::range::repack_again ].user() ) {
		repack_again_ = option[ OptionKeys::relax::range::repack_again ]();
	}

	if ( option[ OptionKeys::relax::range::cycles ].user() ) {
		cycles_ = option[ OptionKeys::relax::range::cycles ]();
	}

	if ( option[ OptionKeys::relax::range::min_cycles ].user() ) {
		min_cycles_ = option[ OptionKeys::relax::range::min_cycles ]();
	}

	if ( option[ OptionKeys::relax::range::spherical_wave ].user() ) {
		spherical_wave_ = option[ OptionKeys::relax::range::spherical_wave ]();
	}

	if ( option[ OptionKeys::relax::range::idealize ].user() ) {
		idealize_ = option[ OptionKeys::relax::range::idealize ]();
	}

	if ( option[ OptionKeys::mp::transform::optimize_embedding ].user() ) {
		optmem_ = option[ OptionKeys::mp::transform::optimize_embedding ]();
	}

}// init from cmdline

//////////////////////////////////////////////////////////////////////

/// @brief Finalize setup
core::kinematics::FoldTree RangeRelaxMover::finalize_setup( Pose & pose ) {

	using namespace core::scoring::constraints;
	using namespace core::kinematics;
	using namespace protocols::membrane;
	using namespace protocols::simple_moves;
	using namespace protocols::moves;

	// add membrane to pose
	if ( addmem_ == true ) {
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );

		// center residue number is pose TM COM; membrane foldtree is set in
		// MPRangeRelaxMover, here we only provide the center residue number
		// to relax from
		center_resnumber_ = rsd_closest_to_pose_tm_com( pose );

		TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	}

	// remember starting foldtree
	FoldTree orig_ft = pose.fold_tree();

	// setup soluble foldtree with rsd close to pose COM being the root
	if ( ! pose.conformation().is_membrane() ) {
		FoldTree ft = pose.fold_tree();
		core::Size anchor( setup_foldtree_pose_com( pose ) );
		ft.reorder( anchor );
		pose.fold_tree( ft );
	}

	// show foldtree
	pose.fold_tree().show( TR );

	// if not set, get center residue number from the pose, which is the COM of the pose
	// this is where we start relaxing from
	if ( center_resnumber_ == 0 ) {
		center_resnumber_ = static_cast< core::Size >( residue_center_of_mass( pose, 1, nres_protein( pose ) ) );
	}

	// get number of residues and number of moves, this doesn't take into
	// account virtual residues
	if ( moves_ == "nres" ) {
		nmoves_ = nres_protein( pose );
	} else {
		nmoves_ = utility::string2Size( moves_ );
	}

	// creating constraint set
	if ( cst_file_.size() > 0 && cst_weight_ > 0 ) {

		TR << "Adding constraints to the pose and scorefunction." << std::endl;
		ConstraintSetOP cst( ConstraintIO::get_instance()->read_constraints( cst_file_, ConstraintSetOP( new ConstraintSet() ), pose ) );
		pose.constraint_set( cst );
		pose.constraint_set()->show_definition( TR, pose );

		// set constraints in scorefunction
		sfxn_->set_weight( atom_pair_constraint, cst_weight_ );
		sfxn_->set_weight( angle_constraint, cst_weight_ );
		sfxn_->set_weight( dihedral_constraint, cst_weight_ );

		// remove constraints in sfxn_nocst_
		TR << "Note: only dealing with constraints for coordinates, atom pairs, angles, and dihedral angles." << std::endl;
		sfxn_nocst_->set_weight( atom_pair_constraint, 0 );
		sfxn_nocst_->set_weight( angle_constraint, 0 );
		sfxn_nocst_->set_weight( dihedral_constraint, 0 );
		sfxn_nocst_->set_weight( coordinate_constraint, 0 );
	}

	// create constraints to native
	if ( cst_to_native_ == true ) {

		if ( nres_protein( pose ) != nres_protein( *native_ ) ) {
			utility_exit_with_message( "Trying to use option -relax::constrain_relax_to_native_coords. The native structure does not have the same number of residues as the Pose. Can't create constraints for poses of different lengths!" );
		}

		// superimpose pose to native so that the coordinate constraint score won't go
		// through the roof
		core::Size nres( nres_protein( pose ) );
		SuperimposeMoverOP super( new SuperimposeMover( *native_, 1, nres, 1, nres, true ) );
		super->apply( pose );

		constrain_to_reference( pose, *native_ );
	}

	// create constraints to starting model
	if ( cst_to_start_ == true ) {

		if ( native_ == 0 ) {
			TR << "WARNING: Setting native to input structure because you chose the option -relax:constrain_relax_to_start_coords!" << std::endl;
			native_ = pose.clone();
		}

		constrain_to_reference( pose, pose );
	}

	return orig_ft;

} // finalize setup

//////////////////////////////////////////////////////////////////////

/// @brief Create constraints to reference pose
void RangeRelaxMover::constrain_to_reference( Pose & pose, Pose & ref_pose ) {

	using namespace basic::options;
	using namespace core::id;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	// go through residues and only constrain CA coordinates
	// we only use CA to make it faster, otherwise it would be too many constraints
	// and lots of time needed for scoring, which is bad for large structures
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {

		// get AtomID: atomno=2 for CA, resno
		AtomID aid( 2, i );
		AtomID aid_fixed( 1, 1 );

		// get "fixed" AtomID, this is only to know when to rescore the pose
		if ( i == 1 ) {
			AtomID aid_fixed( 1, nres_protein( pose ) );
		}

		// get coords to which we want to constraint to
		core::Vector coord = ref_pose.residue( i ).atom( "CA" ).xyz();

		// get constraint function: harmonic around 0 with stdev
		FuncOP harmonic( new HarmonicFunc( 0.0, 1.0 ) );

		// add constraint to pose
		pose.add_constraint( ConstraintCOP( ConstraintOP( new CoordinateConstraint( aid, aid_fixed, coord, harmonic ) ) ) );

	}

	// set weight in sxfn
	sfxn_->set_weight( coordinate_constraint, 0.3 );
	sfxn_nocst_->set_weight( coordinate_constraint, 0 );

} // constrain to ref pose

//////////////////////////////////////////////////////////////////////

/// @brief Repack in a sequence window
void RangeRelaxMover::repack_sequence_window( Pose & pose, MonteCarloOP mc ){

	using namespace core::pack::task;

	// PACKING SIDECHAINS FROM RESIDUE CLOSEST TO COM OUTWARDS
	// start 5 residues lower and repack for 10 residues
	core::SSize halfrange = 4;
	core::SSize range = 2 * halfrange;
	core::SSize nres = static_cast< core::SSize >( pose.total_residue() );

	// create packer task - will be re-used
	PackerTaskOP repack = TaskFactory::create_packer_task( pose );
	repack->restrict_to_repacking();

	// pack the center
	TR << "Packing center..." << std::endl;
	utility::vector1< bool > repack_center( get_window_repack_residues( pose, center_resnumber_, center_resnumber_, halfrange ) );
	repack->restrict_to_residues( repack_center );
	core::pack::pack_rotamers( pose, *sfxn_, repack );

	// go backwards and forwards by 10 residues and repack and minimize those
	// get the center points first for the lower and upper ranges, then get
	// the repack residues from that
	for ( core::SSize i = center_resnumber_, j = center_resnumber_;
			i >= 1 || j <= nres;
			i -= range, j += range ) {

		// pack both lower and upper ranges
		TR << "Packing..." << std::endl;
		utility::vector1< bool > repack_res( get_window_repack_residues( pose, i, j, halfrange ) );
		repack->restrict_to_residues( repack_res );
		core::pack::pack_rotamers( pose, *sfxn_, repack );

	}

	// score pose and print scores
	print_score( pose );

	// evaluate Boltzmann
	mc->boltzmann( pose );

} // repack sequence window

//////////////////////////////////////////////////////////////////////

/// @brief Repack in spherical range
void RangeRelaxMover::repack_spherical_range( Pose & pose, MonteCarloOP mc ) {

	using namespace core::pack::task;

	// PACKING SIDECHAINS FROM CENTER RESIDUE OUTWARDS
	// create packer task - will be re-used
	PackerTaskOP repack = TaskFactory::create_packer_task( pose );
	repack->restrict_to_repacking();

	// find maximum radius in pose
	core::Vector center_coords = pose.residue( center_resnumber_ ).xyz( "CA" );
	core::Real max_radius( 0 );
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		core::Vector coords = pose.residue( i ).xyz( "CA" );
		core::Real radius = ( coords - center_coords ).length();

		if ( radius > max_radius ) {
			max_radius = radius;
		}
	}

	// go outwards in pose in 5 A shells and pack residues
	for ( core::Real i = 0; i < max_radius; i += 5.0 ) {

		TR << "Packing range from " << i << std::endl;
		TR << "center residue " << center_resnumber_ << std::endl;
		utility::vector1< bool > repack_res( get_spherical_repack_residues( pose, i, i+5.0 ) );
		repack->restrict_to_residues( repack_res );
		core::pack::pack_rotamers( pose, *sfxn_, repack );
	}

	// score pose and print scores
	print_score( pose );

	// evaluate Boltzmann
	mc->boltzmann( pose );

} // repack spherical range

//////////////////////////////////////////////////////////////////////

/// @brief Repack all residues again
void RangeRelaxMover::repack_all( Pose & pose, MonteCarloOP mc ) {

	using namespace core::pack::task;

	TR << "Packing all rotamers one more time..." << std::endl;
	PackerTaskOP repack_all = TaskFactory::create_packer_task( pose );
	repack_all->restrict_to_repacking();
	core::pack::pack_rotamers( pose, *sfxn_, repack_all );

	// score pose and print scores
	print_score( pose );

	// evaluate Boltzmann and reset pose
	mc->boltzmann( pose );

	// score pose and print scores
	print_score( pose );

} // repack all residues

//////////////////////////////////////////////////////////////////////

utility::vector1< bool > RangeRelaxMover::get_window_repack_residues( Pose & pose, core::SSize center1, core::SSize center2, core::SSize halfrange ){

	core::SSize nres = static_cast< core::SSize >( pose.total_residue() );
	core::SSize m, n;

	// initialize vector with false
	utility::vector1< bool > repack_res( nres, false );

	TR << "center1 " << center1 << " center2 " << center2 << std::endl;

	// go through residues in the lower range
	for ( core::SSize i = center1, j = center1; i >= center1-halfrange || j <= center1+halfrange; --i, ++j ) {

		// reset counters when out of the pose
		i <= 1 || i >= nres-1 ? m = 1 : m = i;
		j >= nres-1 || j <= 1 ? n = nres-1 : n = j;

		repack_res[ m ] = true;
		repack_res[ n ] = true;
	}

	// go through residues in the upper range
	for ( core::SSize i = center2, j = center2; i >= center2-halfrange || j <= center2+halfrange; --i, ++j ) {

		// reset counters when out of the pose
		i <= 1 || i >= nres-1 ? m = 1 : m = i;
		j >= nres-1 || j <= 1 ? n = nres-1 : n = j;

		repack_res[ m ] = true;
		repack_res[ n ] = true;
	}

	return repack_res;

}// get pack residues

//////////////////////////////////////////////////////////////////////

/// @brief Initialize from commandline
utility::vector1< bool > RangeRelaxMover::get_spherical_repack_residues( Pose & pose, core::Real inner_radius, core::Real outer_radius ){

	core::Size nres = pose.total_residue();

	// initialize vector with false
	utility::vector1< bool > repack_res( nres, false );

	// get CA coordinates of center residue
	core::Vector center_coords = pose.residue( center_resnumber_ ).xyz( "CA" );

	// go through residues in the pose, skip virtuals
	// if CA is between inner and outer radius to the center residue CA, then
	//   set bool to true
	for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {

		// get coordinates of this residue
		core::Vector coords = pose.residue( i ).xyz( "CA" );
		core::Real distance = ( coords - center_coords ).length();

		// if residue within shell, set repack residue to true
		if ( distance >= inner_radius && distance < outer_radius ) {
			TR << "res " << i << " set to true" << std::endl;
			repack_res[ i ] = true;
		}
	}

	return repack_res;

}// get pack residues

//////////////////////////////////////////////////////////////////////

/// @brief Idealize pose
void RangeRelaxMover::idealize_pose( Pose & pose, MinimizerOptions minopts, AtomTreeMinimizer atm  ){

	using namespace core::id;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;

	TR << "Idealizing relaxed pose..." << std::endl;

	// remember original constraint set:
	// we don't remove constraints here because even when we idealize, we might
	// want to keep the constraints that were set by the user because idealize
	// can mess up models quite a bit
	ConstraintSetOP original_cst_set( pose.constraint_set()->clone() );

	// make a copy of the scorefunction and use the copy for minimization with
	// idealize constraints; original sfxn won't be touched!
	ScoreFunctionOP cst_sfxn( sfxn_->clone() );

	core::Size nres( nres_protein( pose ) );
	core::Real dist_dev( 0.25 );
	core::Real const heavyatom_dis2_threshold( 5.5 * 5.5 );
	core::Real const polarH_dis2_threshold( 2.5 * 2.5 );

	// go through residues
	for ( core::Size i = 1; i <= nres; ++i ) {

		Residue i_rsd( pose.residue( i ) );

		// go through residues again
		for ( core::Size j = i + 1; j <= nres; ++j ) {

			Residue j_rsd( pose.residue( j ) );

			// go through atoms of residue 1
			for ( core::Size ii = 1; ii <= i_rsd.natoms(); ++ii ) {

				chemical::AtomType it( i_rsd.atom_type( ii ) );

				// go through atoms of residue 2
				for ( core::Size jj = 1; jj <= j_rsd.natoms(); ++jj ) {

					chemical::AtomType jt( j_rsd.atom_type( jj ) );
					Real dis2( i_rsd.xyz( ii ).distance_squared( j_rsd.xyz( jj ) ) );

					// only look at atoms within a certain threshold
					if ( ( it.is_polar_hydrogen() && jt.is_acceptor()  && dis2 <    polarH_dis2_threshold ) ||
							( jt.is_polar_hydrogen() && it.is_acceptor()  && dis2 <    polarH_dis2_threshold ) ||
							( it.is_heavyatom()      && jt.is_heavyatom() && dis2 < heavyatom_dis2_threshold ) ) {

						FuncOP harmonic( new HarmonicFunc( std::sqrt( dis2 ), dist_dev ) );
						AtomPairConstraintOP pair( new AtomPairConstraint( AtomID(ii,i), AtomID(jj,j), harmonic ) );
						pose.add_constraint( ConstraintCOP( ConstraintOP( pair ) ) );

					}
				} // jj
			} // ii
		} // j>=i
	} // i

	// add constraint weight to scorefunction
	cst_sfxn->set_weight( atom_pair_constraint, 0.03 );

	// keep prolines closed during idealizations.
	cst_sfxn->set_weight( pro_close, 0.5 );

	// keep disulphides together.
	cst_sfxn->set_weight( dslf_ss_dst, 0.5 );//SG-SG bond length
	cst_sfxn->set_weight( dslf_cs_ang, 2.0 );//CB-SG-SG covalent bond angle

	// idealize the pose
	for ( core::Size i = 1; i <= nres; ++i ) {
		idealize_position( i, pose.conformation() );
	}

	// minimize all
	TR << "Minimizing after idealizing decoy..." << std::endl;
	atm.run( pose, *movemap_, *cst_sfxn, minopts );

	// find non-ideal positions
	// for ( core::Size i = 1; i <= nres; ++i ){
	//  TR << "ideal? " << i << " " << is_ideal_position( i, pose.conformation() ) << std::endl;
	// }

	// remove idealize constraints and restore the original constraint set
	pose.remove_constraints();
	pose.constraint_set( original_cst_set );

} // idealize pose

//////////////////////////////////////////////////////////////////////

/// @brief Print score to cout
void RangeRelaxMover::print_score( Pose & pose ) {

	// print energies and iteration
	Real tot_score = ( *sfxn_ )( pose );
	Real fa_rep = pose.energies().total_energies()[ scoring::fa_rep ];
	Real fa_atr = pose.energies().total_energies()[ scoring::fa_atr ];

	pose.energies().show_total_headers( TR );
	TR << std::endl;
	pose.energies().show_totals( TR );
	TR << std::endl;

	TR << "tot: " << tot_score << " fa_rep: " << fa_rep << " fa_atr: " << fa_atr << std::endl;

} // print score

} // relax
} // protocols

#endif // INCLUDED_protocols_relax_RangeRelaxMover_cc
