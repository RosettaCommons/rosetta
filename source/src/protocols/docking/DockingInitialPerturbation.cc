// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DockingInitialPerturbation.cc
/// @brief initial position functions
/// @detailed
///		This contains the functions that create initial positions for docking
///		You can either randomize partner 1 or partner 2, spin partner 2, or
///		perform a simple perturbation.
/// 	Also contains docking mcm protocol
/// @author Monica Berrondo

#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/DockingInitialPerturbationCreator.hh> // zhe

// Rosetta Headers
#include <protocols/moves/Mover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/RigidBodyInfo.hh> // zhe
#include <basic/datacache/DataMap.hh> // zhe

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/options/option.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers
#include <string>

//Utility Headers
// AUTO-REMOVED #include <numeric/conversions.hh>

#include <numeric/trig.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <basic/Tracer.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/tag/Tag.hh>
using basic::T;

// option key includes

#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.docking.DockingInitialPerturbation");
static core::Size trans ( 1 ), rot ( 2 );

using namespace core;

namespace protocols {
namespace docking {


// Creator part for DockingInitialPerturbation, used in scripts
std::string
DockingInitialPerturbationCreator::keyname() const
{
	return DockingInitialPerturbationCreator::mover_name();
}

protocols::moves::MoverOP
DockingInitialPerturbationCreator::create_mover() const
{
	return new DockingInitialPerturbation;
}

std::string
DockingInitialPerturbationCreator::mover_name()
{
	return "DockingInitialPerturbation";
}


// initial perturbation on one of the partners
// the type of perturbation is defined in the options
// some of the options are randomize1 (randomizes the first partner)
// randomize2 (randomizes the second partner), dock_pert, and spin
//------------------------------------------------------------------------------
//
//     there are several ways to perturb the structure before beginning
//     the search; they are controlled through command-line flags
//
//     at the end, partners are slid into contact and scored
//


//constructors
DockingInitialPerturbation::DockingInitialPerturbation()
 :
	Mover(),
	slide_( true ),
	rigid_body_info_( NULL )
{
	Mover::type( "DockingInitialPerturbation" );
	movable_jumps_ = utility::tools::make_vector1<core::Size>(1);
	multiple_jumps_ = false;
	init();
}

DockingInitialPerturbation::DockingInitialPerturbation(
		core::Size const rb_jump,
		bool const slide
) :
	Mover(),
	slide_(slide),
	rigid_body_info_( NULL )
{
	Mover::type( "DockingInitialPerturbation" );
	movable_jumps_ = utility::tools::make_vector1<core::Size>(rb_jump);
	multiple_jumps_ = false;
	init();
}

DockingInitialPerturbation::DockingInitialPerturbation(
		DockJumps const movable_jumps,
		bool const slide
) :
	Mover(),
	slide_( slide ),
	movable_jumps_( movable_jumps ),
	rigid_body_info_( NULL )
{
	Mover::type( "DockingInitialPerturbation" );
	if ( movable_jumps_.size() > 1 ) multiple_jumps_ = true;
	else multiple_jumps_ = false;
	init();
}

void
DockingInitialPerturbation::init()
{
	set_default();
 	init_from_options();  // put this into apply in case scripts is used, then this will not be needed.
}

void
DockingInitialPerturbation::set_default()
{
	randomize1_ = false;
	randomize2_ = false;
	if_dock_pert_ = false;
	if_uniform_trans_ = false;
	spin_ = false;
	center_at_interface_ = false;
	//	dock_pert_ = new utility::vector1< Real >(NULL);
	//	uniform_trans_ = NULL;
}

protocols::moves::MoverOP
DockingInitialPerturbation::clone() const {
	return new DockingInitialPerturbation( *this );
}

//destructor
DockingInitialPerturbation::~DockingInitialPerturbation() {}

// ALERT!
// register_options() and init_from_options() are not called anywhere yet!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!11

void
DockingInitialPerturbation::init_from_options()
{
	using namespace basic::options;
	TR << "Reading options..." << std::endl;

	if ( option[ OptionKeys::docking::randomize1 ].user() )
		set_randomize1(option[ OptionKeys::docking::randomize1 ]());

	if ( option[ OptionKeys::docking::randomize2 ].user() )
		set_randomize2(option[ OptionKeys::docking::randomize2 ]());

	if ( option[ OptionKeys::docking::dock_pert ].user() )
		set_dock_pert(option[ OptionKeys::docking::dock_pert ]());


	if ( option[ OptionKeys::docking::uniform_trans ].user() )
		set_uniform_trans(option[ OptionKeys::docking::uniform_trans ]());

	if ( option[ OptionKeys::docking::spin ].user() )
		set_spin(option[ OptionKeys::docking::spin ]());

	if ( option[ OptionKeys::docking::center_at_interface ].user() )
		set_center(option[ OptionKeys::docking::center_at_interface ]());
}

void
DockingInitialPerturbation::register_options()
{
	using namespace basic::options;

	option.add_relevant( OptionKeys::docking::randomize1 );
	option.add_relevant( OptionKeys::docking::randomize2 );
	option.add_relevant( OptionKeys::docking::dock_pert );
	option.add_relevant( OptionKeys::docking::uniform_trans );
	option.add_relevant( OptionKeys::docking::spin );
	option.add_relevant( OptionKeys::docking::center_at_interface );
}


////////////////////////////////////////////////////////////////////////////////
/// @begin initial_perturbation
///
/// @brief   Make starting perturbations for rigid body moves
///
/// @detailed    There are several ways to perturb the structure before beginning
///     the search; they are controlled through command-line flags
///     At the end, partners are slid into contact and scored (including
///     mc_reset).
///     Also, they are tested for passing the FAB filter.
///
///
/// @references see dock_structure or pose_docking from rosetta++
///
/// @authors Monica Berrondo June 14 2007
///
/// @last_modified October 17 2007
/////////////////////////////////////////////////////////////////////////////////
void DockingInitialPerturbation::apply( core::pose::Pose & pose )
{
	if ( rigid_body_info_ ) {
		movable_jumps_ = rigid_body_info_->movable_jumps();
		TR.Debug << "finished reading movable_jumps_ from RigidBodyInfo" << std::endl;
		if ( movable_jumps_.empty() ) {
			utility_exit_with_message( "DockSetupMover has to be applied before DockingInitialPerturbation !" );
		}
	}


	runtime_assert( !movable_jumps_.empty() );

	for ( DockJumps::const_iterator it=movable_jumps_.begin(); it != movable_jumps_.end(); ++it){
		apply_body( pose, *it );
		TR.Debug <<"movable_jumps_ value in apply:" << *it << std::endl;
	}
}

void
DockingInitialPerturbation::apply_body(core::pose::Pose & pose, core::Size jump_number )
{
	using namespace moves;

	if ( randomize1_ ) {
		TR << "randomize1: true" << std::endl;
		rigid::RigidBodyRandomizeMover mover( pose, jump_number, rigid::partner_upstream );
		mover.apply( pose );
	}
	if ( randomize2_ ) {
		TR << "randomize2: true" << std::endl;
		rigid::RigidBodyRandomizeMover mover( pose, jump_number, rigid::partner_downstream );
		mover.apply( pose );
	}
	if ( if_dock_pert_ ) {
		// DO NOT supply default values for this option -- reasonable values differ for protein and ligand protocols.
		// Also, default values will cause a perturbation to *always* occur, even with no command line flags -- very surprising.
		// Adding defaults WILL BREAK existing protocols in unexpected ways.
		// Decided by Jeff, Monica, Ian, and Sid in March 2008.
		//
		// Jeff notes that eventually there should be 3 parameters, like Rosetta++:
		// rotation, normal translation, and perpendicular translation.
		TR << "dock_pert: true" << std::endl;
		/// read in dock_pert options from commandline.  the first value is the
		/// rotation magnitude and the second value is the translational value
		utility::vector1< Real > pert_mags = dock_pert_;
		//TR << "option[ docking::rotational ]()" << rot << "\n";
		//TR << "option[ docking::parallel ]()" << trans << "\n";
		TR << "option[ docking::dock_pert ]()" << pert_mags[rot] << ' ' << pert_mags[trans] << std::endl;
		rigid::RigidBodyPerturbMoverOP mover;
		if (center_at_interface_) mover = new rigid::RigidBodyPerturbMover( jump_number, pert_mags[rot], pert_mags[trans], rigid::partner_downstream, true );
		else mover = new rigid::RigidBodyPerturbMover( jump_number, pert_mags[rot], pert_mags[trans] );
		mover->apply( pose );
	}

	if ( if_uniform_trans_ ) {
		core::Real mag( uniform_trans_ );
		TR << "uniform_trans: " << mag << std::endl;
		rigid::UniformSphereTransMover mover( jump_number, mag );
		mover.apply( pose );
	}
	if ( spin_ ) {
		TR << "axis_spin: true" << std::endl;
		rigid::RigidBodySpinMover mover( jump_number );
		mover.apply( pose );
	}
	// DO NOT do this for e.g. ligand docking
	if ( slide_ && !pose.is_fullatom() ) {
		DockingSlideIntoContact slide( jump_number );
		slide.apply( pose );
		TR.Debug << "centroid mode, DockingSlideIntoContact applied" << std::endl;
	} else if ( slide_ ) {
		FaDockingSlideIntoContact slide( jump_number );
		slide.apply( pose );
		TR.Debug << "fa-standard mode, FaDockingSlideIntoContact applied" << std::endl;
	}


}
std::string
DockingInitialPerturbation::get_name() const {
	return "DockingInitialPerturbation";
}

void
DockingInitialPerturbation::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	if ( !data_map.has( "RigidBodyInfo", "docking_setup" ) ) {
		TR << "RigidBodyInfo not found in basic::datacache::DataMap" << std::endl;
		rigid_body_info_ = new protocols::docking::RigidBodyInfo;
		data_map.add( "RigidBodyInfo", "docking_setup", rigid_body_info_ );
		//		throw utility::excn::EXCN_RosettaScriptsOption( "RigidBodyInfo not found in basic::datacache::DataMap, DockingInitialPerturbation can not be done, so exit here!" );
	} else {
		rigid_body_info_ = data_map.get< protocols::docking::RigidBodyInfo* >( "RigidBodyInfo", "docking_setup" );
		TR.Debug << "get RigidBodyInfo pointer from basic::datacache::DataMap" << std::endl;
	}

	if ( tag->hasOption( "randomize1" ) ) {
		set_randomize1( tag->getOption<bool>( "randomize1" ) );
	}

	if ( tag->hasOption( "randomize2" ) ) {
		set_randomize2( tag->getOption<bool>( "randomize2" ) );
	}

	if ( tag->hasOption( "dock_pert" ) && tag->getOption<bool>( "dock_pert" ) ) {
		dock_pert_.push_back( tag->getOption<core::Real>("trans" ) );
		dock_pert_.push_back( tag->getOption<core::Real>( "rot" ) );
		set_dock_pert( dock_pert_ );
	}

	if ( tag->hasOption( "uniform_trans" ) ) {
		set_uniform_trans( tag->getOption<core::Real>( "uniform_trans" ) );
	}

	if ( tag->hasOption( "spin" ) ) {
		set_spin( tag->getOption<bool>( "spin" ) );
	}

	if ( tag->hasOption( "center_at_interface" ) ) {
		set_center( tag->getOption<bool>( "center_at_interface" ) );
	}

	slide_ = tag->getOption<bool>( "slide", true );
}


////////////////////////////////////////// DockingSlideIntoContact ////////////////////////////////

// default constructor
DockingSlideIntoContact::DockingSlideIntoContact() : Mover()
{
	using namespace core::scoring;
	Mover::type( "DockingSlideIntoContact" );
	rb_jump_ = 1;
	scorefxn_ = ScoreFunctionFactory::create_score_function( CENTROID_WTS, DOCK_LOW_PATCH );
	scorefxn_ = ScoreFunctionFactory::create_score_function( "interchain_cen" );
}

//constructor
DockingSlideIntoContact::DockingSlideIntoContact(
		core::Size const rb_jump
) : Mover(), rb_jump_(rb_jump)
{
	using namespace core::scoring;
	Mover::type( "DockingSlideIntoContact" );
	scorefxn_ = ScoreFunctionFactory::create_score_function( CENTROID_WTS, DOCK_LOW_PATCH );
	scorefxn_ = ScoreFunctionFactory::create_score_function( "interchain_cen" );
}

//destructor
DockingSlideIntoContact::~DockingSlideIntoContact() {}


void DockingSlideIntoContact::apply( core::pose::Pose & pose )
{
	using namespace moves;

	rigid::RigidBodyTransMover mover( pose, rb_jump_ );
	( *scorefxn_ )( pose );

	TR << "sliding into contact" << std::endl;
	TR << "Moving away" << std::endl;
	core::Size const counter_breakpoint( 500 );
	core::Size counter( 0 );
	// first try moving away from each other
	while ( pose.energies().total_energies()[ scoring::interchain_vdw ] > 0.1 && counter <= counter_breakpoint ) {
		mover.apply( pose );
		( *scorefxn_ )( pose );
		++counter;
	}
	if( counter > counter_breakpoint ){
		TR<<"failed moving away with original vector. Aborting DockingSlideIntoContact."<<std::endl;
		set_current_tag( "fail" );
		return;
	}
	counter = 0;
	// then try moving towards each other
	TR << "Moving together" << std::endl;
	mover.trans_axis().negate();
	while ( counter <= counter_breakpoint && pose.energies().total_energies()[ scoring::interchain_vdw ] < 0.1 ) {
		mover.apply( pose );
		( *scorefxn_ )( pose );
		++counter;
	}
	if( counter > counter_breakpoint ){
		TR<<"moving together failed. Aborting DockingSlideIntoContact."<<std::endl;
		set_current_tag( "fail" );
		return;
	}
	// move away again until just touching
	mover.trans_axis().negate();
	mover.apply( pose );
}

std::string
DockingSlideIntoContact::get_name() const {
	return "DockingSlideIntoContact";
}

void
DockingSlideIntoContact::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Jump number: " << get_jump_num() << std::endl;
}

std::ostream &operator<< ( std::ostream & os, DockingSlideIntoContact const & mover )
{
	mover.show(os);
	return os;
}


////////////////////////////////////////// FaDockingSlideIntoContact ////////////////////////////////

// default constructor
FaDockingSlideIntoContact::FaDockingSlideIntoContact()
{
	Mover::type( "FaDockingSlideIntoContact" );
	scorefxn_ = new core::scoring::ScoreFunction();
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
}


//constructor
FaDockingSlideIntoContact::FaDockingSlideIntoContact(
		core::Size const rb_jump
) : Mover(), rb_jump_(rb_jump), tolerance_(0.2)
{
	Mover::type( "FaDockingSlideIntoContact" );
	scorefxn_ = new core::scoring::ScoreFunction();
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
}

FaDockingSlideIntoContact::FaDockingSlideIntoContact( utility::vector1<core::Size> rb_jumps):
		Mover(), rb_jumps_(rb_jumps),tolerance_(0.2){
	Mover::type( "FaDockingSlideIntoContact" );
	scorefxn_ = new core::scoring::ScoreFunction();
	scorefxn_->set_weight( core::scoring::fa_rep, 1.0 );
}

//destructor
FaDockingSlideIntoContact::~FaDockingSlideIntoContact() {}

void FaDockingSlideIntoContact::apply( core::pose::Pose & pose )
{
	using namespace core::scoring;

	// A very hacky way of guessing whether the components are touching:
	// if pushed together by 1A, does fa_rep change at all?
	// (The docking rb_* score terms aren't implemented as of this writing.)
	(*scorefxn_)( pose );
	core::Real const initial_fa_rep = pose.energies().total_energies()[ fa_rep ];
	bool are_touching = false;

	utility::vector1<rigid::RigidBodyTransMover> trans_movers;

	if(rb_jumps_.size()<1){
		trans_movers.push_back( rigid::RigidBodyTransMover(pose,rb_jump_));
	}
	else{
		for(
				utility::vector1<core::Size>::iterator jump_idx = rb_jumps_.begin(),
				end = rb_jumps_.end();
				jump_idx != end;
				jump_idx++
		){
			trans_movers.push_back( rigid::RigidBodyTransMover(pose, *jump_idx));
		}
	}

	utility::vector1< rigid::RigidBodyTransMover >::iterator const end(trans_movers.end());

	//int i=1;
	// Take 2A steps till clash, then back apart one step.  Now you're within 2A of touching.
	// Repeat with 1A steps, 0.5A steps, 0.25A steps, etc until you're as close are you want.
	for( core::Real stepsize = 2.0; stepsize > tolerance_; stepsize /= 2.0 ) {
		for(
				utility::vector1< rigid::RigidBodyTransMover >::iterator trans_mover(trans_movers.begin());
				trans_mover != end;
				trans_mover++
		){
			trans_mover->trans_axis( trans_mover->trans_axis().negate() ); // now move together
			trans_mover->step_size(stepsize);
		}
		core::Size const counter_breakpoint( 500 );
		core::Size counter( 0 );
		do
		{
			for(
					utility::vector1< rigid::RigidBodyTransMover >::iterator trans_mover(trans_movers.begin());
					trans_mover != end;
					trans_mover++
			){
				trans_mover->apply( pose );
			}
			(*scorefxn_)( pose );
			core::Real const push_together_fa_rep = pose.energies().total_energies()[ fa_rep ];
			//std::cout << "fa_rep = " << push_together_fa_rep << std::endl;
			are_touching = (std::abs(initial_fa_rep - push_together_fa_rep) > 1e-4);
			//std::ostringstream s;
			//s << "snapshot" << i << ".pdb";
			//pose.dump_pdb(s.str());
			//i += 1;
			++counter;
		} while( counter <= counter_breakpoint && !are_touching );
		if( counter > counter_breakpoint ){
			TR<<"Failed Fadocking Slide Together. Aborting."<<std::endl;
			set_current_tag( "fail" );
		}
		for(
				utility::vector1< rigid::RigidBodyTransMover >::iterator trans_mover(trans_movers.begin());
				trans_mover != end;
				trans_mover++
		){
			trans_mover->trans_axis( trans_mover->trans_axis().negate() ); // now move apart
			trans_mover->apply( pose );
		}
	}
}

std::string
FaDockingSlideIntoContact::get_name() const {
	return "FaDockingSlideTogether";
}

void
FaDockingSlideIntoContact::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Jump number: " << get_jump_num() << "\nTolerance:   " << get_tolerance() << std::endl;
}

std::ostream &operator<< ( std::ostream &os, FaDockingSlideIntoContact const &fadock )
{
	fadock.show(os);
	return os;
}

} // namespace docking
} // namespace protocols
