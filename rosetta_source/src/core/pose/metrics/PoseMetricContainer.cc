// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/metrics/PoseMetricContainer.cc
/// @brief  PoseMetricContainer class
/// @author John Karanicolas

// Unit headers
#include <core/pose/metrics/PoseMetricContainer.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>

#include <core/pose/Pose.hh>
#include <core/pose/signals/ConformationEvent.hh>
#include <core/pose/signals/DestructionEvent.hh>
#include <core/pose/signals/EnergyEvent.hh>
#include <core/conformation/Conformation.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.fwd.hh>
#include <utility/exit.hh>

// C++ Headers
#include <map>

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
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
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
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
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
#include <core/pose/metrics/PoseMetricCalculatorBase.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreFunctionInfo.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Constraints.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
// AUTO-REMOVED #include <basic/OStream.fwd.hh>
#include <utility/stream_util.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <utility/down_cast.hh>
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
#include <utility/signals/Link.hh>
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
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/function.hpp>


namespace core {
namespace pose {
namespace metrics {


/// @brief default destructor
PoseMetricContainer::PoseMetricContainer() :
	pose_ptr_( NULL ),
	structure_is_outdated_( true ),
	energies_are_outdated_( true )
{}


/// @brief copy constructor
/// @warning observation of subject Pose in src not duplicated
PoseMetricContainer::PoseMetricContainer( PoseMetricContainer const & src ) :
	ReferenceCount(),
	pose_ptr_( NULL ), // this is correct behavior
	structure_is_outdated_( src.structure_is_outdated_ ),
	energies_are_outdated_( src.energies_are_outdated_ )
{
	clone_calculators( src.metric_calculators_ );
}


/// @brief default destructor
PoseMetricContainer::~PoseMetricContainer() {
	detach_from();
}


/// @brief copy assignment
/// @warning container keeps observing whatever Pose it was originally observing
PoseMetricContainer & PoseMetricContainer::operator =( PoseMetricContainer const & src ) {
	if ( this != &src ) {
		// subject Pose in src is *not* copied, we retain whatever Pose we were originally
		// observing
		clear(); // flags reset here, we'll need an update after this
		clone_calculators( src.metric_calculators_ );
	}
	return *this;
}


/// @brief clear the list of metric calculators, reset flags
void PoseMetricContainer::clear() {
	metric_calculators_.clear();
	structure_is_outdated_ = true;
	energies_are_outdated_ = true;
}


/// @brief get a value out of a PoseMetricCalculator
// TODO: should we really be passing a Pose pointer in...?  Consider using the internal pose_ptr_;
void PoseMetricContainer::get( std::string const & calculator_name, std::string const & key, basic::MetricValueBase & val, Pose const & this_pose ) {
	get_calculator(calculator_name)->get( key, val, this_pose );
}


/// @brief return a string with the results of a PoseMetricCalculator
std::string PoseMetricContainer::print( std::string const & calculator_name, std::string const & key, Pose const & this_pose ) {
	return get_calculator(calculator_name)->get( key, this_pose );
}


/// @brief attach to a Pose
void PoseMetricContainer::attach_to( core::pose::Pose & pose ) {
	if ( pose_ptr_ ) {
		detach_from();
	}

	pose.attach_conformation_obs( &PoseMetricContainer::on_pose_conf_change, this );
	pose.attach_energy_obs( &PoseMetricContainer::on_pose_energy_change, this );
	pose.attach_destruction_obs( &PoseMetricContainer::on_pose_destruction_change, this );
	pose_ptr_ = &pose;

	// new pose, so assume we need a full update
	structure_is_outdated_ = true;
	energies_are_outdated_ = true;
}


/// @brief detach from Pose
void PoseMetricContainer::detach_from() {
	if ( pose_ptr_ ) {
		pose_ptr_->detach_conformation_obs( &PoseMetricContainer::on_pose_conf_change, this );
		pose_ptr_->detach_energy_obs( &PoseMetricContainer::on_pose_energy_change, this );
		pose_ptr_->detach_destruction_obs( &PoseMetricContainer::on_pose_destruction_change, this );
	}

	pose_ptr_ = NULL;
}


/// @brief is observing a Pose?
/// @return the Pose if observing, otherwise NULL
Pose const * PoseMetricContainer::is_observing() {
	return pose_ptr_;
}


/// @brief upon receiving a pose::signals::DestructionEvent, detaches
void PoseMetricContainer::on_pose_destruction_change( core::pose::signals::DestructionEvent const & ) {
	detach_from();
}


/// @brief upon receiving pose::signals::ConformationEvent, sets flag telling
///  calculators to refresh structure based calculations
void PoseMetricContainer::on_pose_conf_change( core::pose::signals::ConformationEvent const & ) {
	structure_is_outdated_ = true;
}


/// @brief upon receiving a pose:signals::EnergyEvent, sets flag telling
///  calculators to refresh energy based calculations
void PoseMetricContainer::on_pose_energy_change( core::pose::signals::EnergyEvent const & ) {
	energies_are_outdated_ = true;
}


/// @brief set PoseMetricCalculator structure changed flags
void PoseMetricContainer::process_structure_change() {
	// notify all calculators that the structure is outdated
	std::map< std::string, PoseMetricCalculatorOP >::iterator iter_calculators;
	for ( iter_calculators = metric_calculators_.begin(); iter_calculators != metric_calculators_.end(); ++iter_calculators ) {
		iter_calculators->second->notify_structure_change();
	}

	// reset global flag
	structure_is_outdated_ = false;
}


/// @brief set PoseMetricCalculator energy changed flags
void PoseMetricContainer::process_energy_change() {

	// notify all calculators that the energy is outdated
	Name2Calculator::iterator iter_calculators;
	for ( iter_calculators = metric_calculators_.begin(); iter_calculators != metric_calculators_.end(); ++iter_calculators ) {
		iter_calculators->second->notify_energy_change();
	}

	// reset global flag
	energies_are_outdated_ = false;
}


/// @brief get a PoseMetricCalculator by name
PoseMetricCalculatorOP PoseMetricContainer::get_calculator( std::string const & calculator_name ) {

	if ( !pose_ptr_ ) {
		basic::Error() << "This PoseMetricContainer is not observing a Pose" << std::endl;
		utility_exit();
	}

	// propagate news of the outdated structure/sequence/energy to the metrics and calculators
	if ( structure_is_outdated_ ) {
		process_structure_change();
	}

	//	if ( sequence_is_outdated_ ) process_sequence_change();

	if ( energies_are_outdated_ ) {
		process_energy_change();
	}

	// call the "get" function for the appropriate metric
	Name2Calculator::const_iterator calculator_iter;
	calculator_iter = metric_calculators_.find( calculator_name );
	if (calculator_iter == metric_calculators_.end() ) {
		// this calculator has not yet been setup, add it
		add_calculator( calculator_name );
		calculator_iter = metric_calculators_.find( calculator_name );
		if (calculator_iter == metric_calculators_.end() ) {
			basic::Error() << "Could not lookup calculator " << calculator_name << " despite trying to add it" << std::endl;
			utility_exit();
		}
	}

	return calculator_iter->second;
}


/// @brief add a PoseMetricCalculator by name
void PoseMetricContainer::add_calculator( std::string const & calculator_name ) {
	PoseMetricCalculatorOP new_calculator = CalculatorFactory::Instance().retrieve_calculator( calculator_name );
	metric_calculators_.insert( std::make_pair( calculator_name, new_calculator ) );
}


/// @brief clone calculators, primarily for copy
void PoseMetricContainer::clone_calculators( Name2Calculator const & calcs ) {
	for ( Name2Calculator::const_iterator i = calcs.begin(), ie = calcs.end(); i != ie; ++i ) {
		metric_calculators_.insert( std::make_pair( i->first, i->second->clone() ) );
	}
}

/*
void PoseMetricContainer::process_sequence_change() {
	// notify all calculators that the sequence is outdated
	Name2Calculator:iterator iter_calculators;
	for ( iter_calculators = metric_calculators_.begin(); iter_calculators != metric_calculators_.end(); ++iter_calculators ) {
		iter_calculators->second->notify_sequence_change();
	}
	// reset global flag
	sequence_is_outdated_ = false;
	return;
}
*/


} // metrics
} // pose
} // core
