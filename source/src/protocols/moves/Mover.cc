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
/// @author Vikram K. Mulligan (vmullig@uw.edu) -- Added citation manager functions.

// Unit Headers
#include <protocols/moves/Mover.hh>

// Package headers

// Project headers
#include <core/pose/Pose.hh>
#include <protocols/jobdist/Jobs.hh>
#include <utility/string_util.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/filters/Filter.hh>

// Basic Headers
#include <basic/citation_manager/Citation.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

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

static basic::Tracer TR( "protocols.moves.Mover" );

Mover::Mover()
:
	utility::VirtualBase(),
	type_( "MoverBase" ),
	current_tag_( "NoTag" ),
	input_pose_(/* 0 */),
	native_pose_(/* 0 */),
	last_status_( MS_SUCCESS )
{}

Mover::Mover( std::string const & type_name ) :
	utility::VirtualBase(),
	type_( type_name ),
	current_tag_( "NoTag" ),
	input_pose_(/* 0 */),
	native_pose_(/* 0 */),
	last_status_( MS_SUCCESS )
{}


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

/// @details Some movers need not be parsed, so we shouldn't force people to reimplement this method.
/// However, we should be chatty about the fact that someone is using a RosettaScripts interface to a mover which didn't define parse_my_tag()
void Mover::parse_my_tag(
	TagCOP,
	basic::datacache::DataMap &
) {
	TR.Warning << "The parse_my_tag() for mover " << name() << " has been invoked, but it hasn't been overridden. Are you sure this is appropriate?" << std::endl;
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
	return nullptr;
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
	return MoverOP(nullptr);
	//return utility::pointer::make_shared< Mover >(); //this is what your Mover should return - it's illegal here because Mover does not define the pure virtual function apply().
}

/// @details clone is meant to return an OP'ed deep copy of this object.  This really should be a pure virtual in the base class, but adding pure virtuals to Mover would massively disrupt the code.  This default implementation crashes at runtime instead of compiletime if you try to call it.  If this code is causing you problems, your Mover needs to override this function.
MoverOP Mover::clone() const {
	utility_exit_with_message( "clone has been called on a Mover which has not overridden the base class implementation.  Probably you tried to pass a Mover to the job distributor or parser which does not have clone implemented.  Implement the function and try again.\n");
	return MoverOP(nullptr);
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

//////////////////FUNCTIONS NEEDED FOR THE CITATION MANAGER/////////////////////////////

/// @brief Does this mover provide information about how to cite it?
/// @details Defaults to false.  Derived classes may override this to provide citation info.  If set to
/// true, the provide_citation_info() override should also be provided.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
Mover::mover_provides_citation_info() const { return false; }

/// @brief Provide the citation.
/// @returns A vector of citation collections.  This allows the mover to provide citations for
/// itself and for any modules that it invokes.
/// @details The default implementation of this function provides an empty vector.  It may be
/// overriden by Movers wishing to provide citation information.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::CitationCollectionCOP >
Mover::provide_citation_info() const {
	return utility::vector1< basic::citation_manager::CitationCollectionCOP >();
}

/// @brief Does this mover indicate that it is unpublished (and, by extension, that the author should be
/// included in publications resulting from it)?
/// @details Defaults to false.  Derived classes may override this to provide authorship info.  If set to
/// true, the provide_authorship_info_for_unpublished() override should also be provided.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
Mover::mover_is_unpublished() const { return false; }

/// @brief Provide a list of authors and their e-mail addresses, as strings.
/// @returns A list of pairs of (author, e-mail address).  Empty list if not unpublished.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
Mover::provide_authorship_info_for_unpublished() const {
	return utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >();
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
