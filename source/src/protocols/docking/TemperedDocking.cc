// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TemperedDocking
/// @brief low resolution docking using parallel or simulated tempering
/// @detailed
/// @author Oliver Lange

//Unit Headers
#include <protocols/docking/TemperedDocking.hh>
#include <protocols/docking/TemperedDockingCreator.hh>

//Package Headers
#include <protocols/docking/util.hh>

//Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/canonical_sampling/SimulatedTempering.hh>
#include <protocols/canonical_sampling/MpiParallelTempering.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/canonical_sampling/SilentTrajectoryRecorder.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tools/make_vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

//#include <numeric/trig.functions.hh>
//#include <numeric/xyzMatrix.fwd.hh>

// C++ Headers
#include <string>

#include <basic/Tracer.hh>

//Auto Headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <core/scoring/ScoreFunctionFactory.hh>


// cmd-line options


bool protocols::docking::TemperedDocking::options_registered_( false );

void protocols::docking::TemperedDocking::register_options() {
  using namespace basic::options;
  using namespace OptionKeys;
  if ( options_registered_ ) return;
  options_registered_ = true;
	protocols::canonical_sampling::SimulatedTempering::register_options();
	protocols::canonical_sampling::MpiParallelTempering::register_options();
	protocols::canonical_sampling::SilentTrajectoryRecorder::register_options();
	OPT( docking::partners );
	OPT( score::patch );
	OPT( constraints::cst_file );
	OPT( run::n_cycles );
	OPT( run::n_replica );
	OPT( rigid::translation );
	OPT( rigid::rotation );
}

static basic::Tracer tr("protocols.docking.TemperedDocking");

using namespace core;

namespace protocols {
namespace docking {

//// ----------------------------------- BEGIN CONSTRUCTORS --------------------------------------------------
// Constructors
TemperedDocking::TemperedDocking() {
	init( utility::tools::make_vector1<core::SSize>(1), NULL, NULL);
}

void TemperedDocking::init(
	DockJumps const movable_jumps,
	core::scoring::ScoreFunctionCOP docking_score_low,
	core::scoring::ScoreFunctionCOP docking_score_high
) {
	Mover::type( "TemperedDocking" );
	movable_jumps_ = movable_jumps;
	docking_scorefxn_low_ = docking_score_low;
	docking_scorefxn_high_= docking_score_high;

	// setup all the booleans with default values
	// they will get overwritten by the options and/or passed values
	set_defaults();
	init_from_options();

	// set up objects based on the boolean values defined above
	setup_objects();
}

///@brief clone operator, calls the copy constructor
protocols::moves::MoverOP
TemperedDocking::clone() const {
	return new TemperedDocking(*this);
}

///@brief copy ctor
TemperedDocking::TemperedDocking( TemperedDocking const & rhs ) :
	Mover(rhs)
{
	copy(*this, rhs);
}

///@brief assignment operator
TemperedDocking & TemperedDocking::operator=( TemperedDocking const & rhs ){
	//abort self-assignment
	if (this == &rhs) return *this;
	Mover::operator=(rhs);
	copy(*this, rhs);
	return *this;
}

void TemperedDocking::copy(TemperedDocking & lhs, TemperedDocking const & rhs){

	lhs.autofoldtree_ = rhs.autofoldtree_;
	lhs.flags_and_objects_are_in_sync_ = rhs.flags_and_objects_are_in_sync_;
	lhs.first_apply_with_current_setup_ = rhs.first_apply_with_current_setup_;
	lhs.sc_min_ = rhs.sc_min_;

	lhs.use_csts_ = rhs.use_csts_;

	lhs.fold_tree_ = rhs.fold_tree_;
	lhs.partners_ = rhs.partners_;
	lhs.previous_sequence_ = rhs.previous_sequence_;

	lhs.movable_jumps_ = rhs.movable_jumps_;
	lhs.docking_scorefxn_low_ = rhs.docking_scorefxn_low_->clone();
	lhs.docking_scorefxn_high_ = rhs.docking_scorefxn_high_->clone();

	lhs.to_centroid_ = rhs.to_centroid_->clone();

	if( rhs.docking_constraint_ )	lhs.docking_constraint_ = rhs.docking_constraint_->clone();
}

//// ----------------------------------- END CONSTRUCTORS --------------------------------------------------

///@brief setup steps when new pose comes in (this is not repeated for every reply, ... only when new sequence.. )
void
TemperedDocking::finalize_setup( pose::Pose & pose ) //setup objects requiring pose during apply
{
	if( autofoldtree_ ) {
		tr.Info << "Setting docking foldtree" << std::endl;
		tr.Info << "old fold tree: " << pose.fold_tree() << std::endl;
		setup_foldtree( pose, partners_, movable_jumps_ );
		tr.Info << "new fold tree: " << pose.fold_tree() << std::endl;
	}
	fold_tree_ = pose.fold_tree();

	rb_mover_->clear_jumps();
	for ( Size i=1; i<=pose.num_jump(); ++i ) {
		// should honor movemap or movable_jumps_ here...
		rb_mover_->add_jump( i );
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @begin Docking apply
///
/// @brief main docking protocol
///
/// @detailed
///
void
TemperedDocking::apply( pose::Pose & pose )
{
	/// -------------- any changes ? ------------------------------------------------------
	if ( !flags_and_objects_are_in_sync_ ){
		sync_objects_with_flags();
	}

	/// ---------------- do we have to call finalize_setup() again ? -------------------------////
	if ( previous_sequence_.compare( pose.sequence() ) != 0 ){
		first_apply_with_current_setup_ = true;
		previous_sequence_ = pose.sequence();
	}

	if ( first_apply_with_current_setup_ ){
		finalize_setup( pose );
		first_apply_with_current_setup_ = false;
	}

	pose.fold_tree( fold_tree_ );
	show(tr.Info);

	/// go into centroid mode
	to_centroid_->apply( pose );

	/// apply constraints
	if ( docking_constraint_ ) {
		tr.Info << "setting up the constraint set mover" << std::endl;
		docking_constraint_->apply( pose );
	}

	sampler_->apply( pose );


}

//// --------------------------------- Setters -------------------------------------------------


void TemperedDocking::set_lowres_scorefxn( core::scoring::ScoreFunctionOP docking_scorefxn_low )
{
	docking_scorefxn_low_ = docking_scorefxn_low;
}

void TemperedDocking::set_highres_scorefxn( core::scoring::ScoreFunctionOP scfxn )
{
	set_highres_scorefxn( scfxn, scfxn, scfxn );
}

void TemperedDocking::set_highres_scorefxn( // delete
	core::scoring::ScoreFunctionCOP high,
	core::scoring::ScoreFunctionCOP pack )
{
	set_highres_scorefxn( high, pack, high );
}

void TemperedDocking::set_highres_scorefxn(
	core::scoring::ScoreFunctionCOP docking_scorefxn_high,
	core::scoring::ScoreFunctionCOP docking_scorefxn_pack,
	core::scoring::ScoreFunctionCOP docking_scorefxn_output )
{
	docking_scorefxn_high_ = docking_scorefxn_high;
	docking_scorefxn_pack_ = docking_scorefxn_pack;
	docking_scorefxn_output_ = docking_scorefxn_output;
}

void TemperedDocking::set_sc_min( bool ) {}

void TemperedDocking::set_use_constraints( bool setting )
{
	if ( use_csts_ != setting ) flags_and_objects_are_in_sync_ = false;
	use_csts_ = setting;
}

/// ----------------------------------------- Getters ------------------------------------------

//getters for const access to movers and data of docking protocol
protocols::moves::MoverCOP TemperedDocking::to_centroid() const {
	return to_centroid_;
}



/// ---------------------------- diagnostic output ------------------------------------------

/// @details  Show the complete setup of the docking protocol
void
TemperedDocking::show( std::ostream & out ) const {
	/*if ( !flags_and_objects_are_in_sync_ ){
		sync_objects_with_flags();
	}*/  // show() should be const.
	out << *this;
}

std::ostream & operator<<(std::ostream& out, const TemperedDocking & dp )
{
	using namespace ObjexxFCL::format;

	// All output will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << line_marker << A( 47, "Rosetta 3 Docking Protocol" ) << space( 27 ) << line_marker << std::endl;
	out << line_marker << space( 74 ) << line_marker << std::endl;
	// Display the movable jumps that will be used in docking
	out << line_marker << " Dockable Jumps: ";

	core::Size spaces_so_far = 23;
	bool first = true;
	for ( DockJumps::const_iterator it = dp.movable_jumps_.begin() ; it != dp.movable_jumps_.end() ; ++it ){
		if (!first) {
			out << ", ";
			spaces_so_far += 2;
		}
		else first = false;

		out << I( 1, *it );
		spaces_so_far += 1;
	}
	core::Size remaining_spaces = 80 - spaces_so_far;
	if ( remaining_spaces > 0 )	out << space( 80 - spaces_so_far );
	out << line_marker << std::endl;

	// Close the box I have drawn
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	return out;
}

/// ------------------ initializing helpers ---------------------------------

void
TemperedDocking::set_defaults()
{
	//using namespace basic::options;
	using namespace core::scoring;
	autofoldtree_ = true;

	sc_min_ = false;
	partners_ = "_";
	n_cycles_ = 10000;

}

void TemperedDocking::setup_objects()
{
	// initialize constraint set mover
	docking_constraint_ = NULL;

	// stores the sequence of the previous pose, so that the TemperedDocking can re setup the fold tree
	previous_sequence_ = "";

	fold_tree_ = core::kinematics::FoldTree();

	// Residue movers
	to_centroid_ = new protocols::simple_moves::SwitchResidueTypeSetMover( core::chemical::CENTROID );

	// correctly set up the score functions from either passed in values or defaults

	sampler_ = new protocols::canonical_sampling::MetropolisHastingsMover();
	rb_mover_ = new rigid::RigidBodyPerturbNoCenterMover();

	rb_mover_->rot_magnitude( rigid_rot_mag_ );
	rb_mover_->trans_magnitude( rigid_trans_mag_ );

	moves::MonteCarloOP mc_object = new moves::MonteCarlo( *docking_scorefxn_low_, 0.6 );

	sampler_->set_monte_carlo( mc_object );
	sampler_->set_ntrials( n_cycles_ );
	tempering_->set_monte_carlo( mc_object );

	sampler_->add_mover( rb_mover_, 1.0 );
	sampler_->set_tempering( tempering_ );
	sampler_->add_observer( new protocols::canonical_sampling::SilentTrajectoryRecorder );
	sync_objects_with_flags();
}

void TemperedDocking::sync_objects_with_flags()
{
	// set up constraint set mover
	if ( !use_csts_ ) {
		docking_constraint_ = NULL;
	} else {
		if ( !docking_constraint_ ) {
			docking_constraint_ = new protocols::simple_moves::ConstraintSetMover();
		}
	}

	flags_and_objects_are_in_sync_ = true;
	first_apply_with_current_setup_ = true;
}

void
TemperedDocking::init_from_options()
{
	using namespace basic::options;

	// check for n_cycle option
	if ( option[ OptionKeys::run::n_cycles ].user() ) {
		n_cycles_ = option[ OptionKeys::run::n_cycles ]();
	}

	// check for movable jumps option
	if( option[ OptionKeys::docking::multibody ].user() ) {
		set_movable_jumps( option[ OptionKeys::docking::multibody ]() );
	}

	// Sets the member directly so sync_objects_with_flags won't be triggered prematurely.
	// A public setter exists for this member.
	if( option[ OptionKeys::docking::sc_min ].user() ) {
		sc_min_ = option[ OptionKeys::docking::sc_min ]();
	}

	// This defaults to "_"
	if( option[ OptionKeys::docking::partners ].user() ) {
		set_partners(option[ OptionKeys::docking::partners ]());
	}

	// Defaults to false
	set_use_constraints( option[ OptionKeys::constraints::cst_file ].user() );

	if ( !docking_scorefxn_low_ ) {
		//docking_scorefxn_low_ = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen" );
		docking_scorefxn_low_ = core::scoring::ScoreFunctionFactory::create_score_function( "interchain_cen",
			option[ OptionKeys::score::patch ]()
		);
	}
	if ( !docking_scorefxn_high_ ) {
		docking_scorefxn_high_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
		docking_scorefxn_pack_ = core::scoring::getScoreFunctionLegacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS );
		docking_scorefxn_output_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking" );
	}

	rigid_rot_mag_ = option[ OptionKeys::rigid::rotation ]();
	rigid_trans_mag_ = option[ OptionKeys::rigid::translation ]();

	if ( option[ OptionKeys::run::n_replica ]() > 1 ) {
		tempering_ = new protocols::canonical_sampling::MpiParallelTempering();
	} else {
		tempering_ = new protocols::canonical_sampling::SimulatedTempering();
	}

}


void
TemperedDocking::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	using namespace core::scoring;

	if( tag->hasOption("docking_score_low" ) ){
		std::string const score_low( tag->getOption<std::string>( "docking_score_low" ) );
		ScoreFunctionOP scorelo = data.get< ScoreFunction * >( "scorefxns", score_low )->clone();
		set_lowres_scorefxn(scorelo);
	}
	if( tag->hasOption("docking_score_high" ) ){
		std::string const score_high( tag->getOption<std::string>( "docking_score_high" ) );
		ScoreFunctionOP scorehi = data.get< ScoreFunction * >( "scorefxns", score_high )->clone();
		set_highres_scorefxn(scorehi);
	}
	//get through partners
	if( tag->hasOption( "partners" ) ){
		std::string const partners( tag->getOption<std::string>( "partners") );
		set_partners(partners);
	}

}//end parse_my_tag


/// ------------------ End initializing helpers ---------------------------------

std::string
TemperedDockingCreator::keyname() const
{
	return TemperedDockingCreator::mover_name();
}

protocols::moves::MoverOP
TemperedDockingCreator::create_mover() const {
	return new TemperedDocking();
}

std::string
TemperedDockingCreator::mover_name()
{
	return "TemperedDocking";
}
} //docking
} //protocols
