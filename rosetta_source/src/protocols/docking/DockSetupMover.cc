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
#include <protocols/docking/DockSetupMover.hh>
#include <protocols/docking/DockSetupMoverCreator.hh>

//Package Headers
#include <protocols/docking/util.hh>
#include <protocols/moves/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
//Project Headers
#include <core/pose/Pose.hh>

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


// cmd-line options
#include <basic/options/option_macros.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

static basic::Tracer tr("protocols.docking.DockSetupMover");

using namespace core;

namespace protocols {
namespace docking {

//// ----------------------------------- BEGIN CONSTRUCTORS --------------------------------------------------
// Constructors
DockSetupMover::DockSetupMover() {
	set_defaults();
}

///@brief clone operator, calls the copy constructor
protocols::moves::MoverOP
DockSetupMover::clone() const {
	return new DockSetupMover(*this);
}

///@brief copy ctor
DockSetupMover::DockSetupMover( DockSetupMover const & rhs ) : Mover(rhs) {
	copy(*this, rhs);
}

///@brief assignment operator
DockSetupMover & DockSetupMover::operator=( DockSetupMover const & rhs ) {
	//abort self-assignment
	if (this == &rhs) return *this;
	Mover::operator=(rhs);
	copy(*this, rhs);
	return *this;
}

void DockSetupMover::copy(DockSetupMover & lhs, DockSetupMover const & rhs) {
	lhs.partners_ = rhs.partners_;
	//	lhs.previous_sequence_ = rhs.previous_sequence_;
	//	lhs.movable_jumps_ = rhs.movable_jumps_;
}

//// ----------------------------------- END CONSTRUCTORS --------------------------------------------------
void
DockSetupMover::apply( pose::Pose & pose ) {
	docking::setup_foldtree( pose, partners_, movable_jumps_ );
	rb_mover_->clear_jumps();
	for ( Size i=1; i<=pose.num_jump(); ++i ) {
		// should honor movemap or movable_jumps_ here...
		rb_mover_->add_jump( i );
	}
}

//// --------------------------------- Setters -------------------------------------------------



/// ----------------------------------------- Getters ------------------------------------------


/// ---------------------------- diagnostic output ------------------------------------------

/// @details  Show the complete setup of the docking protocol
void
DockSetupMover::show( std::ostream & out ) {
	out << *this;
}

std::ostream & operator<<(std::ostream& out, const DockSetupMover & dp )
{
	using namespace ObjexxFCL::fmt;

	// All output will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << line_marker << A( 47, "Rosetta 3 DockSetupMover" ) << space( 27 ) << line_marker << std::endl;
	out << line_marker << space( 74 ) << line_marker << std::endl;
	// Close the box I have drawn
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	return out;
}

/// ------------------ initializing helpers ---------------------------------

void
DockSetupMover::set_defaults()
{
	using namespace basic::options;
	partners_ = option[ OptionKeys::docking::partners ]();
}

void
DockSetupMover::parse_my_tag(
  TagPtr const tag,
	moves::DataMap&,
	protocols::filters::Filters_map const&,
	moves::Movers_map const& movers,
	core::pose::Pose const&
) {
	using namespace core::scoring;
	//get through partners
	if( tag->hasOption( "partners" ) ){
		std::string const partners( tag->getOption<std::string>( "partners") );
		set_partners(partners);
	}
	moves::MoverOP mover = rosetta_scripts::parse_mover( tag->getOption< std::string >( "rb_mover", "null" ), movers );
	rb_mover_ = dynamic_cast< moves::RigidBodyPerturbNoCenterMover* >( mover() );
	if ( !rb_mover_ ) {
		utility_exit_with_message( "DockSetupMover requires an rb_mover argument" );
	}
	movable_jumps_.clear();
	if ( tag->hasOption( "moveable_jump" ) ) {
		movable_jumps_.push_back( tag->getOption< core::Size >( "moveable_jump" ));
	}
	if (partners_ == "_" and movable_jumps_.size()<1 ) {
		movable_jumps_.push_back( 1 );
	}
}//end parse_my_tag


/// ------------------ End initializing helpers ---------------------------------

std::string
DockSetupMoverCreator::keyname() const
{
	return DockSetupMoverCreator::mover_name();
}

protocols::moves::MoverOP
DockSetupMoverCreator::create_mover() const {
	return new DockSetupMover();
}

std::string
DockSetupMoverCreator::mover_name()
{
	return "DockSetupMover";
}
} //docking
} //protocols
