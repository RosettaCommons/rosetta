// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/RdcEvaluatorCreator.hh
/// @brief  Header for RdcEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/RdcEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/RDC_Evaluator.hh>

#include <core/chemical/ResidueType.hh>
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


static THREAD_LOCAL basic::Tracer tr( "protocols.evalution.RdcEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

RdcEvaluatorCreator::~RdcEvaluatorCreator() {}

void RdcEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::rdc );

}

void RdcEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;


	if ( option[ OptionKeys::evaluation::rdc ].user() ) {
		typedef utility::vector1< std::string > RdcVector;
		RdcVector const& rdc( option[ OptionKeys::evaluation::rdc ]() );
		utility::vector1< core::Size> empty_selection;
		for ( RdcVector::const_iterator it=rdc.begin(); it!=rdc.end(); ++it ) {
			std::string fname( *it );
			std::string column;
			++it;
			if ( it != rdc.end() ) {
				column = *it;
			} else {
				utility_exit_with_message(
					"need to specify dupletts <rdcs> <column> with option -evaluation:rdc   last read: "+fname );
			}
			eval.add_evaluation( PoseEvaluatorOP( new simple_filters::SelectRDC_Evaluator( empty_selection, column, fname ) ) );
		} // iterate over tripletts in option -rmsd
	}
	if ( option[ OptionKeys::evaluation::built_in_rdc ].user() ) {
		eval.add_evaluation( PoseEvaluatorOP( new simple_filters::RDC_Evaluator(option[ OptionKeys::evaluation::built_in_rdc ]()) ) );
	} // iterate over tripletts in option -rmsd
}

std::string RdcEvaluatorCreator::type_name() const {
	return "RdcEvaluatorCreator";
}

} //namespace
} //namespace
