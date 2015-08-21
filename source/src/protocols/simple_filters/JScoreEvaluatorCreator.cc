// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/JScoreEvaluatorCreator.hh
/// @brief  Header for JScoreEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/JScoreEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/JScoreEvaluator.hh>

#include <core/io/silent/silent.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

//// C++ headers

// due to template function
#include <core/io/silent/SilentStruct.hh>


// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <utility/vector0.hh>

#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh>
#endif


static thread_local basic::Tracer tr( "protocols.evalution.JScoreEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

JScoreEvaluatorCreator::~JScoreEvaluatorCreator() {}

void JScoreEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::jscore_evaluator );

}

void JScoreEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;


	if ( option[ OptionKeys::evaluation::jscore_evaluator ].user() ) {
		using std::string;
		using utility::vector1;
		using core::scoring::ScoreFunctionOP;
		vector1< string > const & tags( option[ OptionKeys::evaluation::jscore_evaluator ]() );

		for ( Size ii = 1; ii <= tags.size() - 1; ii += 2 ) {
			string scorefxn_name( tags[ii]   );
			string rsd_set_name ( tags[ii+1] );

			eval.add_evaluation( PoseEvaluatorOP( new JScoreEvaluator( scorefxn_name, rsd_set_name ) ) );
		}
	}


}

std::string JScoreEvaluatorCreator::type_name() const {
	return "JScoreEvaluatorCreator";
}

} //namespace
} //namespace
