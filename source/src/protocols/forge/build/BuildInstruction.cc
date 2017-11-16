// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/components/instructions/BuildInstruction.cc
/// @brief tracks modifications to be made and is capable of residue length changes on a Pose
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/build/BuildInstruction.hh>

// project headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <basic/Tracer.hh>

// utility headers
#include <utility/signals/Link.hh>

// C++ headers
#include <algorithm>

#ifdef WIN32
#include <iterator>
#endif

#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace forge {
namespace build {


// Tracer instance for this file
// Named after the original location of this code
static basic::Tracer TR( "protocols.forge.build.BuildInstruction" );


/// @brief default constructor
BuildInstruction::BuildInstruction() :
	Super(),
	state_( BuildInstructionState::READY ),
	original_interval_( Interval( 0, 0 ) ),
	detach_after_modify_( true )
{}


/// @brief interval constructor
/// @param[in] i the residue range to operate on
/// @param[in] rts the residue type set to use, default FA_STANDARD
BuildInstruction::BuildInstruction(
	Interval const & i,
	ResidueTypeSetCAP rts
) :
	Super(),
	state_( BuildInstructionState::READY ),
	original_interval_( i ),
	rts_( rts ),
	detach_after_modify_( true )
{}


/// @brief copy constructor
BuildInstruction::BuildInstruction( BuildInstruction const & rval ) :
	Super( rval ),
	state_( rval.state_ ),
	original_interval_( rval.original_interval_ ),
	rts_( rval.rts_ ),
	detach_after_modify_( rval.detach_after_modify_ )
{
	// length_obs_link_ not copied!
	// dependencies_ not copied!  See comments in header.
}


/// @brief default destructor
BuildInstruction::~BuildInstruction() {
	detach_from();
}


/// @brief copy assignment
BuildInstruction & BuildInstruction::operator =( BuildInstruction const & rval ) {
	if ( this != &rval ) {
		Super::operator =( rval );
		state_ = rval.state_;
		original_interval_ = rval.original_interval_;
		rts_ = rval.rts_;
		// length_obs_link_ not copied on purpose!
		detach_after_modify_ = rval.detach_after_modify_;
		// dependencies_ not copied on purpose!  See comments in header.
	}
	return *this;
}


/// @brief modify this pose
BuildInstructionState::Enum BuildInstruction::modify( Pose & pose ) {
	// Can't call modify() more than once without invoking reset_accounting().
	// Can consider replacing this with a call to reset_accounting() itself...
	runtime_assert( !modify_was_successful() );

	// Before doing anything, check to make sure the necessary dependencies
	// are satisfied so that the BuildInstruction has the information it
	// needs to complete successfully.
	if ( !dependencies_satisfied() ) {
		state( BuildInstructionState::WAITING_ON_DEPENDENCIES );
		return state_;
	}

	// we're mangling the structure so must fully clear Energies, otherwise
	// things may break in users' protocols
	pose.energies().clear();

	attach_to( pose );

	// do the actual work
	modify_impl( pose );

	if ( detach_after_modify_ ) {
		detach_from();
	}

	// mark state
	state( BuildInstructionState::MODIFY_WAS_SUCCESSFUL );

	return state_;
}


/// @brief Is the BuildInstruction's state at <em>READY</em>?
/// @remarks <em>READY</em> indicates the BuildInstruction has been reset
///  and is ready to modify a Pose.
bool BuildInstruction::ready() const {
	return state_ == BuildInstructionState::READY;
}


/// @brief Is the BuildInstruction's state at <em>WAITING_ON_DEPENDENCIES</em>?
/// @remarks <em>WAITING_ON_DEPENDENCIES</em> indicates
///  the BuildInstruction is waiting for its dependencies to be satisfied
///  before allowing modifications to proceed.
bool BuildInstruction::waiting_on_dependencies() const {
	return state_ == BuildInstructionState::WAITING_ON_DEPENDENCIES;
}


/// @brief Is the BuildInstruction's state at <em>MODIFY_WAS_SUCCESSFUL</em>?
/// @remarks <em>MODIFY_WAS_SUCCESSFUL</em> indicates the BuildInstruction has
///  finished modifications to the Pose, and its residue indexing is now
///  consistent with the newly modified Pose.
bool BuildInstruction::modify_was_successful() const {
	return state_ == BuildInstructionState::MODIFY_WAS_SUCCESSFUL;
}


/// @brief attach to a Pose's conformation
void BuildInstruction::attach_to( Pose & pose ) {
	detach_from();
	length_obs_link_ = pose.conformation().attach_length_obs( &BuildInstruction::on_length_change, this );
}


/// @brief detach from a Pose's conformation
void BuildInstruction::detach_from() {
	length_obs_link_.invalidate();
}


/// @brief update any indexing wrt length change to Pose/Conformation being watched
void BuildInstruction::on_length_change( LengthEvent const & event ) {
	switch ( event.tag ) {
	case ( LengthEvent::INVALIDATE ) : {
		detach_from();
		break;
	}
	case ( LengthEvent::RESIDUE_APPEND ) : {
		on_residue_append( event );
		break;
	}
	case ( LengthEvent::RESIDUE_PREPEND ) : {
		on_residue_prepend( event );
		break;
	}
	case ( LengthEvent::RESIDUE_DELETE ) : {
		on_residue_delete( event );
		break;
	}
	default : { // do nothing, fall through
		break;
	}
	}
}


/// @brief Reset intervals, positions, etc to initial state and drop
///  observers.  State set to READY.
void BuildInstruction::reset_accounting() {
	detach_from();
	reset_accounting_impl();
	state( BuildInstructionState::READY );
}


/// @brief are dependencies satisfied so that modify_impl() can complete
///  successfully?
/// @return True if modify_impl() can be called, False otherwise.
/// @remarks This function will automatically be checked within modify()
///  prior to running modify_impl().
bool BuildInstruction::dependencies_satisfied() const {
	for ( BuildInstructionCAPs::const_iterator i = dependencies_.begin(), ie = dependencies_.end(); i != ie; ++i ) {
		protocols::forge::build::BuildInstructionCOP instruction(*i);
		if ( !instruction->modify_was_successful() ) { // some instruction has not completed yet
			return false;
		}
	}

	// all instructions have completed, everything ok
	return true;
}


/// @brief add an instruction to this BuildInstruction's dependency list
/// @remarks before this instruction's modify() may complete properly, all
///  instructions in the dependency list must have completed modify() first
void BuildInstruction::add_dependency_to( BuildInstructionCAP i ) {
	if ( !dependent_on( i ) ) {
		dependencies_.push_back( i );
	}
}


/// @brief is this instruction dependent upon the given instruction?
bool BuildInstruction::dependent_on( BuildInstructionCAP i ) const {
	// Doesn't work with weak_ptr:
	// return std::find( dependencies_.begin(), dependencies_.end(), i ) != dependencies_.end();
	for ( BuildInstructionCAP d : dependencies_ ) {
		if ( utility::pointer::equal(d, i) ) {
			return true;
		}
	}
	return false;
}


/// @brief compares fixed and mutable positions to determine compatibility with
///  another instruction
/// @remarks override this to obtain custom compatibility check
bool BuildInstruction::compatible_with( BuildInstruction const & rval ) const {
	using std::insert_iterator;
	using std::inserter;

	// First check residue type sets are equivalent via their names.
	// This check should be performed unless mixing-matching of different
	// residue type sets within a Conformation is implemented.

	// Compare by name:
	// if ( rts_.lock()->name() != rval.rts_.lock()->name() )

	// Compare pointers (faster):
	if ( !utility::pointer::equal(rts_, rval.rts_ ) ) {
		TR.Error << "ResidueTypeSet incompatibility!" << std::endl;
		return false;
	}

	Positions fixed = original_fixed_positions();
	Positions muta = original_mutable_positions();

	Positions rval_fixed = rval.original_fixed_positions();
	Positions rval_muta = rval.original_mutable_positions();

	Positions iset; // holds the intersections of two sets

	// check mutable regions do not overlap
	set_intersection( muta.begin(), muta.end(), rval_muta.begin(), rval_muta.end(),
		inserter( iset, iset.begin() ) );
	if ( !iset.empty() ) {
		TR.Error << "mutable regions intersect; incompatibility!" << std::endl;
		return false;
	}
	iset.clear();

	// check fixed region and mutable region from right does not overlap
	set_intersection( fixed.begin(), fixed.end(), rval_muta.begin(), rval_muta.end(),
		inserter( iset, iset.begin() ) );
	if ( !iset.empty() ) {
		TR.Error << "fixed and mutable regions intersect; incompatibility!" << std::endl;
		return false;
	}
	iset.clear();

	// check fixed regions from right and mutable region does not overlap
	set_intersection( rval_fixed.begin(), rval_fixed.end(), muta.begin(), muta.end(),
		inserter( iset, iset.begin() ) );
	if ( !iset.empty() ) {
		TR.Error << "mutable and fixed regions intersect; incompatibility!" << std::endl;
		return false;
	}
	iset.clear();

	return true;
}


} // namespace build
} // namespace forge
} // namespace protocols
