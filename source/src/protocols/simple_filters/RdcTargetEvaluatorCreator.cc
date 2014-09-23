// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/RdcTargetEvaluatorCreator.hh
/// @brief  Header for RdcTargetEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/RdcTargetEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/RDC_Evaluator.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.hh>

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
#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>

#ifdef WIN32
	#include <core/scoring/constraints/Constraint.hh>
#endif


static thread_local basic::Tracer tr( "protocols.evalution.RdcTargetEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

RdcTargetEvaluatorCreator::~RdcTargetEvaluatorCreator() {}

void RdcTargetEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::rdc_target );
	OPT( evaluation::rdc_column );
}

void RdcTargetEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;


  if ( option[ OptionKeys::evaluation::rdc_target ].user() ) {
    /*
      this creates Evaluators to get RMSD and GDTMM for as many structures as you specify in this FileVector option.
      pls: provide also as many column-names to match your rmsd-targets.
      ---
      missing density in the input file will be ignored, e.g., this is your way to signal on which residues the rmsd should be computed
    */
    utility::vector1< std::string > const& rdc_target( option[ OptionKeys::evaluation::rdc_target ]() );
    utility::vector1< std::string > const& rdc_col_name( option[ OptionKeys::evaluation::rdc_column ]() );
    for ( Size ct = 1; ct <= rdc_target.size(); ct ++ ) {
      pose::PoseOP rdc_pose = new pose::Pose;
      core::import_pose::pose_from_pdb( *rdc_pose, rdc_target[ ct ] );
      std::string tag( ObjexxFCL::string_of( ct ) );
      if ( rdc_col_name.size() >= ct ) tag = rdc_col_name[ ct ];
      eval.add_evaluation( new simple_filters::SelectRDC_Evaluator( rdc_pose, tag ) );
    }
  }

}

std::string RdcTargetEvaluatorCreator::type_name() const {
	return "RdcTargetEvaluatorCreator";
}

} //namespace
} //namespace
