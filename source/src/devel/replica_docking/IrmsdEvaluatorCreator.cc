// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Zhe Zhang

// Unit Headers
#include <devel/replica_docking/WrapFilterAsEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/SingleValuePoseEvaluator.fwd.hh>
#include <protocols/evaluation/SingleValuePoseEvaluator.hh>
#include <devel/replica_docking/IrmsdEvaluator.hh>
#include <core/io/silent/silent.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// due to template function
#include <core/io/silent/SilentStruct.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <utility/vector0.hh>

//Auto Headers


static basic::Tracer tr( "devel.replica_docking.IrmsdEvaluatorCreator" );

namespace devel {
namespace replica_docking {

IrmsdEvaluatorCreator::~IrmsdEvaluatorCreator() {}

void IrmsdEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::Irms );

}

void IrmsdEvaluatorCreator::add_evaluators( protocols::evaluation::SingleValuePoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	if ( option[ OptionKeys::evaluation::Irms ]() ) {

		tr << "Add evaluator of IRMSD " << std::endl;
		eval.add_evaluation( new IrmsdEvaluator( "Irms") );

	}
}

std::string IrmsdEvaluatorCreator::type_name() const {
	return "IrmsdEvaluator";
}

} //namespace
} //namespace
