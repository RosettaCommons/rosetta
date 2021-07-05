// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/JumpNrEvaluatorCreator.hh
/// @brief  Header for JumpNrEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/JumpNrEvaluatorCreator.hh>

// Package Headers

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/JumpEvaluator.hh>

// ObjexxFCL Headers

// Utility headers

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// due to template function


// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>


//Auto Headers


#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh>
#endif


static basic::Tracer tr( "protocols.evalution.JumpNrEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

JumpNrEvaluatorCreator::~JumpNrEvaluatorCreator() = default;

void JumpNrEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::jump_nr );

}

void JumpNrEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {

	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;


	if ( option[ OptionKeys::evaluation::jump_nr ]() ) {
		eval.add_evaluation( utility::pointer::make_shared< JumpNrEvaluator >() );
	}

}

std::string JumpNrEvaluatorCreator::type_name() const {
	return "JumpNrEvaluatorCreator";
}

} //namespace
} //namespace
