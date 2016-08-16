// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/Mover.cc
/// @brief Method code and full headers for Mover--
/// keeps heavily-included Mover.hh small and concise to maximize compiling
/// efficiency and to make the class definitions easier to read.
/// @author Monica Berrondo
/// @author Jeff Gray
/// @author Steven Lewis
/// @author Sarel Fleishman

// Unit Headers
#include <protocols/moves/Mover.hh>

// Package headers

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/jobdist/Jobs.hh>
#include <utility/string_util.hh>

// tracer
#include <basic/Tracer.hh>

#include <protocols/jobdist/Jobs.hh>
#include <utility/exit.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace moves {

using namespace core;
using namespace pose;

using basic::T;
using basic::Error;
using basic::Warning;
static THREAD_LOCAL basic::Tracer TR( "protocols.moves.Mover" );

Mover::Mover()
:
	utility::pointer::ReferenceCount(),
	type_( "MoverBase" ),
	current_tag_( "NoTag" ),
	input_pose_(/* 0 */),
	native_pose_(/* 0 */),
	last_status_( MS_SUCCESS )
{}

Mover::~Mover(){}

Mover::Mover( std::string const & type_name ) :
	utility::pointer::ReferenceCount(),
	type_( type_name ),
	current_tag_( "NoTag" ),
	input_pose_(/* 0 */),
	native_pose_(/* 0 */),
	last_status_( MS_SUCCESS )
{}

Mover::Mover( Mover const & other ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< Mover >(),
	type_( other.type_ ),
	current_tag_( other.current_tag_ ),
	input_pose_(other.input_pose_),
	native_pose_(other.native_pose_),
	last_status_( other.last_status_ ),
	current_job_( other.current_job_ )
{}


/// @brief assignment operator
Mover& Mover::operator=( Mover const & rhs ) {
	//abort self-assignment
	if ( this == &rhs ) return *this;
	type_ = rhs.type_;
	current_tag_ = rhs.current_tag_;
	input_pose_ = rhs.input_pose_;
	native_pose_ = rhs.native_pose_;
	last_status_ = rhs.last_status_;
	current_job_ = rhs.current_job_;
	return *this;
}

PoseCOP
Mover::get_input_pose() const { return input_pose_; }

void
Mover::set_input_pose( PoseCOP pose ) { input_pose_ = pose; }

PoseCOP
Mover::get_native_pose() const { return native_pose_; }

void
Mover::set_native_pose( PoseCOP pose ) { native_pose_ = pose; }

// Factory<Mover> functions
MoverOP Mover::create() {
	return fresh_instance();
}

// end mpr support
/// @details Some movers need not be parsed, so we shouldn't stop executions. This, however, calls attention to the lack of this method, which could be due to something as silly as a wrong parameters definition.
void Mover::parse_my_tag(
	TagCOP const,
	basic::datacache::DataMap &,
	Filters_map const &,
	Movers_map const &,
	Pose const &
)
{
	TR << "***WARNING!!!! WARNING!!!*** parse_my_tag has been invoked for this mover but it hasn't been defined. Are you sure this is appropriate?" << std::endl;
}

void
Mover::set_current_job( protocols::jobdist::BasicJobCOP job )
{
	current_job_ = job;
}


jobdist::BasicJobCOP
Mover::get_current_job() const {
	return current_job_;
}


/// @details Mechanism by which a mover may return multiple output poses from a single input pose.  After apply is called, calling this function will return subsequent output poses.  NULL is returned if the mover has no more output poses remaining.  The base class implementation returns NULL always; multi-output movers must override this.  RosettaScripts uses this interface to produce a 'branched' protocol.
core::pose::PoseOP
Mover::get_additional_output() {
	return NULL;
}

//////////////////////////start Job Distributor interface//////////////////////////////////
/// @details used by job distributor to get status
MoverStatus Mover::get_last_move_status() const { return last_status_; }

/// @details The job distributor (august 08 vintage) uses this to ensure non-accumulation of status across apply()s.
void Mover::reset_status() { last_status_ = MS_SUCCESS; }

/// @details use this function for implementing filtering in your protocol - failed jobs should set the status to something other than "SUCCESS".
void Mover::set_last_move_status( MoverStatus status ) { last_status_ = status; }

/// @details Movers default to not regenerating
bool Mover::reinitialize_for_each_job() const { return false; }

/// @details Movers default to not regenerating
bool Mover::reinitialize_for_new_input() const { return false; }

/// @details fresh_instance is meant to return a new object of this class, created with the default constructor.  This really should be a pure virtual in the base class, but adding pure virtuals to Mover would massively disrupt the code.  This default implementation crashes at runtime instead of compiletime if you try to call it.  If this code is causing you problems, your Mover needs to override this function.  This is used by the August 08 job distributor.
MoverOP Mover::fresh_instance() const
{
	utility_exit_with_message("fresh_instance has been called on a Mover which has not overridden the base class implementation.  Probably you tried to pass a Mover to the job distributor which does not have fresh_instance implemented.  Implement the function and try again.\n");
	return MoverOP(NULL);
	//return MoverOP( new Mover ); //this is what your Mover should return - it's illegal here because Mover does not define the pure virtual function apply().
}

/// @details clone is meant to return an OP'ed deep copy of this object.  This really should be a pure virtual in the base class, but adding pure virtuals to Mover would massively disrupt the code.  This default implementation crashes at runtime instead of compiletime if you try to call it.  If this code is causing you problems, your Mover needs to override this function.
MoverOP Mover::clone() const {
	utility_exit_with_message( "clone has been called on a Mover which has not overridden the base class implementation.  Probably you tried to pass a Mover to the job distributor or parser which does not have clone implemented.  Implement the function and try again.\n");
	return MoverOP(NULL);
}

// Outputs details about the Mover, including current settings.
/// @details Ideally, a child Mover should call Mover.show() and add additional information particular to that Mover.
void
Mover::show(std::ostream & output) const
{
	output << "Mover name: " << get_name();
	output << ", Mover type: " << get_type();
	output << ", Mover current tag:" << get_current_tag() << std::endl;
}


} // moves

// Insertion operator (overloaded so that the Mover can be "printed").
/// @details This helper function is located in protocols so that all Movers can access it. If a child Mover does not
/// specifically implement the overloading of operator<<, one will still be able to "print" that Mover, because the
/// implementation below will call the child Mover's show() method (which may or may not be inherited from the parent
/// Mover.
// Currently, placing this here triggers other errors that I need to investigate. ~Labonte
//std::ostream &
//operator<<(std::ostream & output, moves::Mover const & mover)
//{
// mover.show(output);
// return output;
//}

} // protocols
