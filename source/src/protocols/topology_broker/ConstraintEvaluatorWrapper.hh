// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @details responsibilities:
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_ConstraintEvaluatorWrapper_hh
#define INCLUDED_protocols_topology_broker_ConstraintEvaluatorWrapper_hh


// Unit Headers
#include <protocols/topology_broker/ConstraintEvaluatorWrapper.fwd.hh>

// Package Headers
#include <protocols/topology_broker/ConstraintClaimer.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers
#include <protocols/evaluation/PoseEvaluator.hh>
// Utility headers

//// C++ headers

// option key includes


namespace protocols {
namespace topology_broker {

class ConstraintEvaluatorWrapper : public evaluation::PoseEvaluator {

public:
	ConstraintEvaluatorWrapper( std::string const& name, ConstraintClaimerCOP claimer ); //for factory

	using evaluation::PoseEvaluator::apply;

	//sets xxx_cst and xxx_viol columns
	virtual void apply( core::pose::Pose& pose, std::string tag, core::io::silent::SilentStruct &pss) const;

	//returns constraint score
	virtual core::Real apply( core::pose::Pose& pose ) const;

	virtual core::Size size() const { return 1; }
	virtual std::string name( core::Size i ) const;

private:
	std::string name_;
	ConstraintClaimerCOP claimer_;
}; //class ConstraintEvaluatorWrapper

}
}

#endif
