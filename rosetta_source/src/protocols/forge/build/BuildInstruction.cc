// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <iterator>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/types.hh>
//XRW_B_T1
//#include <core/coarse/Translator.fwd.hh>
//XRW_E_T1
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/graph/ArrayPool.hh>
#include <core/graph/Graph.fwd.hh>
#include <core/graph/Graph.hh>
#include <core/conformation/PointGraph.fwd.hh>
#include <core/conformation/PointGraphData.fwd.hh>
#include <core/graph/UpperEdgeGraph.fwd.hh>
#include <core/graph/unordered_object_pool.fwd.hpp>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/AtomID_Mask.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomPointer.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/AtomWithDOFChange.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/DomainMap.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/types.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
//#include <core/optimization/MinimizerMap.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/scoring/ContextGraph.fwd.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/scoring/EnergyGraph.fwd.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/LREnergyContainer.fwd.hh>
#include <core/scoring/NeighborList.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/TenANeighborGraph.fwd.hh>
#include <core/scoring/TwelveANeighborGraph.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.fwd.hh>
#include <core/scoring/methods/Methods.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <basic/MetricValue.fwd.hh>
// AUTO-REMOVED #include <basic/OStream.fwd.hh>
#include <utility/stream_util.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableData.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <basic/datacache/DataCache.fwd.hh>
#include <basic/datacache/DataCache.hh>
#include <protocols/forge/build/BuildInstruction.fwd.hh>
#include <protocols/forge/build/Interval.fwd.hh>
#include <protocols/forge/build/Interval.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector0.fwd.hh>
#include <utility/vector0.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/io/all.fwd.hh>
#include <utility/io/icstream.fwd.hh>
#include <utility/io/irstream.fwd.hh>
#include <utility/io/izstream.fwd.hh>
#include <utility/io/ocstream.fwd.hh>
#include <utility/io/orstream.fwd.hh>
#include <utility/io/ozstream.fwd.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/BodyPosition.fwd.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/Quaternion.fwd.hh>
#include <numeric/all.fwd.hh>
#include <numeric/conversions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzTriple.fwd.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray.all.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/config.hpp>
#include <boost/function.hpp>
#include <boost/pool/detail/mutex.hpp>
#include <boost/pool/poolfwd.hpp>



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
		if ( !(**i).modify_was_successful() ) { // some instruction has not completed yet
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
		dependencies_.push_back( i.get() );
	}
}


/// @brief is this instruction dependent upon the given instruction?
bool BuildInstruction::dependent_on( BuildInstructionCAP i ) const {
	return std::find( dependencies_.begin(), dependencies_.end(), i ) != dependencies_.end();
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
	if ( rts_->name() != rval.rts_->name() ) {
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
