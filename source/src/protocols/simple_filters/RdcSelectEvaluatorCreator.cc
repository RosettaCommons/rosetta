// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/RdcSelectEvaluatorCreator.hh
/// @brief  Header for RdcSelectEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/RdcSelectEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <core/io/silent/silent.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>

// due to template function
#include <core/io/silent/SilentStruct.hh>
#include <utility/vector0.hh>

//Auto Headers
#include <basic/options/keys/OptionKeys.hh>


#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh>
#endif


static THREAD_LOCAL basic::Tracer tr( "protocols.evalution.RdcSelectEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

RdcSelectEvaluatorCreator::~RdcSelectEvaluatorCreator() {}

void RdcSelectEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

}

void RdcSelectEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & ) const {

}

std::string RdcSelectEvaluatorCreator::type_name() const {
	return "RdcSelectEvaluatorCreator";
}

} //namespace
} //namespace
