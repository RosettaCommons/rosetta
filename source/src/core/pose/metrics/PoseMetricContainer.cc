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

// Utility headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.fwd.hh>
#include <utility/exit.hh>

// C++ Headers
#include <map>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/kinematics/Jump.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/map.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace metrics {


/// @brief default constructor
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

	// if ( sequence_is_outdated_ ) process_sequence_change();

	if ( energies_are_outdated_ ) {
		process_energy_change();
	}

	// call the "get" function for the appropriate metric
	Name2Calculator::const_iterator calculator_iter;
	calculator_iter = metric_calculators_.find( calculator_name );
	if ( calculator_iter == metric_calculators_.end() ) {
		// this calculator has not yet been setup, add it
		add_calculator( calculator_name );
		calculator_iter = metric_calculators_.find( calculator_name );
		if ( calculator_iter == metric_calculators_.end() ) {
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

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::metrics::PoseMetricContainer::save( Archive & arc ) const {

	// Do not serialize the raw pointer to the Pose being observed; it will
	// be the responsibility of the deserialized Pose to reestablish the
	// observation link.
	// arc( CEREAL_NVP( pose_ptr_ ) ); // const core::pose::Pose *; raw pointer: const core::pose::Pose *
	// EXEMPT pose_ptr_
	arc( CEREAL_NVP( structure_is_outdated_ ) ); // _Bool
	arc( CEREAL_NVP( energies_are_outdated_ ) ); // _Bool
	arc( CEREAL_NVP( metric_calculators_ ) ); // Name2Calculator
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::metrics::PoseMetricContainer::load( Archive & arc ) {
	// Do not deserialize the raw pointer to the Pose being observed
	// arc( pose_ptr_ ); // const core::pose::Pose *; raw pointer: const core::pose::Pose *
	// EXEMPT pose_ptr_
	arc( structure_is_outdated_ ); // _Bool
	arc( energies_are_outdated_ ); // _Bool
	arc( metric_calculators_ ); // Name2Calculator
}
SAVE_AND_LOAD_SERIALIZABLE( core::pose::metrics::PoseMetricContainer );
CEREAL_REGISTER_TYPE( core::pose::metrics::PoseMetricContainer )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_metrics_PoseMetricContainer )
#endif // SERIALIZATION
