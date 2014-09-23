// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraints_additional/ConstraintsEvaluatorCreator.hh
/// @brief  Header for ConstraintsEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/topology_broker/ConstraintEvaluatorWrapperCreator.hh>
#include <protocols/topology_broker/ConstraintEvaluatorWrapper.hh>
// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>
#include <protocols/topology_broker/ConstraintClaimer.hh>
// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/constraints_additional/ConstraintEvaluator.hh>
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

#ifdef WIN32
	#include <core/scoring/constraints/Constraint.hh>
#endif


static thread_local basic::Tracer tr( "protocols.topology_broker.ConstraintEvaluatorWrapperCreator" );

namespace protocols {
namespace topology_broker {

ConstraintEvaluatorWrapperCreator::~ConstraintEvaluatorWrapperCreator() {}

void ConstraintEvaluatorWrapperCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::combined_constraints );
	OPT( evaluation::combined_constraints_column );

}

void ConstraintEvaluatorWrapperCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace evaluation;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

  if ( option[ OptionKeys::evaluation::combined_constraints ].user() ) {
    /*
      this creates Evaluators to evaluate different constraint sets against your decoys
      pls: provide also as many column-names to match your constraint sets
      ---
    */
    utility::vector1< std::string > const& cst_target( option[ OptionKeys::evaluation::combined_constraints ]() );
    utility::vector1< std::string > const& cst_col_name( option[ OptionKeys::evaluation::combined_constraints_column ]() );
    for ( Size ct = 1; ct <= cst_target.size(); ct ++ ) {
      std::string tag( ObjexxFCL::string_of( ct ) );
			if ( cst_col_name.size() >= ct ) tag = cst_col_name[ ct ];
			topology_broker::ConstraintClaimerOP cst =
				new topology_broker::ConstraintClaimer( cst_target[ ct ], tag );
			cst->set_combine_ratio( 2 );
			cst->set_fullatom( true );
			cst->set_centroid( false );
			cst->set_filter_weight( 1 );
			cst->set_skip_redundant( 1 );
      eval.add_evaluation( new ConstraintEvaluatorWrapper( cst->tag(), cst ) );
    }
  }
}

std::string ConstraintEvaluatorWrapperCreator::type_name() const {
	return "ConstraintEvaluator";
}

} //namespace
} //namespace
