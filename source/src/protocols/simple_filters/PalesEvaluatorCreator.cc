// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/PalesEvaluatorCreator.hh
/// @brief  Header for PalesEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/PalesEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

#include <protocols/simple_filters/PalesEvaluator.hh>

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

//Auto Headers


#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh>
#endif


static THREAD_LOCAL basic::Tracer tr( "protocols.evalution.PalesEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

PalesEvaluatorCreator::~PalesEvaluatorCreator() {}

void PalesEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::pales );

}

void PalesEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;


	if ( option[ OptionKeys::evaluation::pales ].user() ) {
		typedef utility::vector1< std::string > CSVector;
		CSVector const& pales( option[ OptionKeys::evaluation::pales ]() );

		for ( CSVector::const_iterator it=pales.begin(); it!=pales.end(); ++it ) {
			std::string fname( *it );
			std::string column;
			++it;
			if ( it != pales.end() ) {
				column = *it;
			} else {
				utility_exit_with_message(
					"need to specify dupletss <pales_rdcs> <column>  with option -evaluation:pales   last read: "+fname );
			}
			eval.add_evaluation( PoseEvaluatorOP( new PalesEvaluator( column, fname ) ) );
		}
	}

}

std::string PalesEvaluatorCreator::type_name() const {
	return "PalesEvaluatorCreator";
}

} //namespace
} //namespace
