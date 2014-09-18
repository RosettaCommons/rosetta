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
/// @detailed
/// @author Oliver Lange

//Unit Headers
#include <protocols/docking/DockSetupMover.hh>
#include <protocols/docking/DockSetupMoverCreator.hh>

//Package Headers
#include <protocols/docking/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/docking/RigidBodyInfo.hh>
#include <basic/datacache/DataMap.hh>
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
#include <utility/excn/Exceptions.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

static thread_local basic::Tracer tr( "protocols.docking.DockSetupMover" );

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
	lhs.rb_mover_ = rhs.rb_mover_;
	lhs.movable_jumps_ = rhs.movable_jumps_;
	lhs.rigid_body_info_ = rhs.rigid_body_info_;
}

//// ----------------------------------- END CONSTRUCTORS --------------------------------------------------
void
DockSetupMover::apply( pose::Pose & pose ) {
	docking::setup_foldtree( pose, partners_, movable_jumps_ );
	//	runtime_assert( rb_mover_ );
	//	rb_mover_->clear_jumps(); //this doesn't work because of cloning --- have to communicate via data-map !
	//	for ( Size i=1; i<=pose.num_jump(); ++i ) {
	  // should honor movemap or movable_jumps_ here...
	for ( Size i=1; i<=movable_jumps_.size(); ++i ) {
		rigid_body_info_->add_jump( movable_jumps_[i] );
		if ( rb_mover_ ) {
			rb_mover_->add_jump( movable_jumps_[i] );
		}
	}
 	if ( rigid_body_info_->movable_jumps().empty() ) {
 		utility_exit_with_message( "RigidBodyInfo was not correctly set!" );
 	}
}

//// --------------------------------- Setters -------------------------------------------------



/// ----------------------------------------- Getters ------------------------------------------


/// ---------------------------- diagnostic output ------------------------------------------

/// @details  Show the complete setup of the docking protocol
void
DockSetupMover::show( std::ostream & out ) const {
	using namespace ObjexxFCL::format;
// All output will be 80 characters - 80 is a nice number, don't you think?
	std::string line_marker = "///";
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << line_marker << A( 47, "Rosetta 3 DockSetupMover" ) << space( 27 ) << line_marker << std::endl;
	out << line_marker << space( 74 ) << line_marker << std::endl;
	out << line_marker << A( 47, "partners: ") << A(10, partners_ ) << space(14) << line_marker << std::endl;
	if ( movable_jumps_.size() ) {
		out << line_marker << A( 47, "movable jumps: " ) << movable_jumps_.front() << space( 26 ) << line_marker<< std::endl;
	}
	// Close the box I have drawn
	out << "////////////////////////////////////////////////////////////////////////////////" << std::endl;
}

std::ostream & operator<<(std::ostream& out, const DockSetupMover & dp ) {

	dp.show( out );
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
  TagCOP const tag,
	basic::datacache::DataMap& data_map,
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
	rb_mover_ = dynamic_cast< rigid::RigidBodyPerturbNoCenterMover* >( mover() );
// 	if ( !rb_mover_ ) {
// 		throw utility::excn::EXCN_RosettaScriptsOption( "DockSetupMover requires an rb_mover argument" );
// 	}
	movable_jumps_.clear();
	if ( tag->hasOption( "moveable_jump" ) ) {
		movable_jumps_.push_back( tag->getOption< core::Size >( "moveable_jump" ));
	}
	if (partners_ == "_" && movable_jumps_.size()<1 ) {
		movable_jumps_.push_back( 1 );
	}
	// using RigidBodyInfo to store movable_jumps, then rb_mover is free from DockSetupMover
	if ( !data_map.has( "RigidBodyInfo", "docking_setup" ) ) {
		// as member variable: RigidBodyInfoOP rigid_body_info_;
		rigid_body_info_ = new protocols::docking::RigidBodyInfo;
		data_map.add( "RigidBodyInfo", "docking_setup", rigid_body_info_ );
		tr.Debug << "added RigidBodyInfo into basic::datacache::DataMap" << std::endl;
	} else {
		rigid_body_info_ = data_map.get< protocols::docking::RigidBodyInfo* >( "RigidBodyInfo", "docking_setup" );
		assert( rigid_body_info_ );
		tr.Debug << "RigidBodyInfo supposed to be in basic::datacache::DataMap, but somehow failed to get it from basic::datacache::DataMap" << std::endl;
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
