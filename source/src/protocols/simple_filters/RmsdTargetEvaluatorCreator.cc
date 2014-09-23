// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/RmsdTargetEvaluatorCreator.hh
/// @brief  Header for RmsdTargetEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/RmsdTargetEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>
#include <protocols/evaluation/util.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/simple_filters/ScoreEvaluator.hh>

#include <core/chemical/ResidueType.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.hh>


#include <core/scoring/ScoreFunction.hh>

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

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>

//Auto Headers


#ifdef WIN32
	#include <core/scoring/constraints/Constraint.hh>
#endif


static thread_local basic::Tracer tr( "protocols.evalution.RmsdTargetEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

RmsdTargetEvaluatorCreator::~RmsdTargetEvaluatorCreator() {}

void RmsdTargetEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::rmsd_target );
	OPT( evaluation::rmsd_target );
	OPT( evaluation::rmsd_column );
	OPT( evaluation::gdtmm );
	OPT( evaluation::score_with_rmsd );
	OPT( evaluation::symmetric_rmsd );
}

void RmsdTargetEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;


  if ( option[ OptionKeys::evaluation::rmsd_target ].user() ) {
    /*
      this creates Evaluators to get RMSD and GDTMM for as many structures as you specify in this FileVector option.
      pls: provide also as many column-names to match your rmsd-targets.
      ---
      missing density in the input file will be ignored, e.g., this is your way to signal on which residues the rmsd should be computed
    */
    utility::vector1< std::string > const & rmsd_target  ( option[ OptionKeys::evaluation::rmsd_target ]() );
		utility::vector1< std::string >         rmsd_col_name;

		if ( option[ OptionKeys::evaluation::rmsd_column ].user() ){
			rmsd_col_name = option[ OptionKeys::evaluation::rmsd_column ]();
		}else{
			// make up default names: firone is just empty, leading to the columns being correctly named "rms", "gdtmm" etc.. the
			// subsequent columns then become "_2", "_3" etc..
			rmsd_col_name.push_back("");
			for( core::Size j=2; j<=rmsd_target.size(); ++j){
				rmsd_col_name.push_back("_" + ObjexxFCL::string_of(j));
			}
		}
    for ( Size ct = 1; ct <= rmsd_target.size(); ct ++ ) {
      pose::PoseOP rmsd_pose( new pose::Pose );
      core::import_pose::pose_from_pdb( *rmsd_pose, rmsd_target[ ct ] );
      std::string tag( ObjexxFCL::string_of( ct ) );
      if ( rmsd_col_name.size() >= ct ) tag = rmsd_col_name[ ct ];
      eval.add_evaluation( PoseEvaluatorOP( new simple_filters::SelectRmsdEvaluator( rmsd_pose, tag ) ) );
      if ( option[ OptionKeys::evaluation::gdtmm ]() ) eval.add_evaluation( PoseEvaluatorOP( new simple_filters::SelectGdtEvaluator( rmsd_pose, tag ) ) );
			if ( option[ OptionKeys::evaluation::score_with_rmsd ]() ){
				core::scoring::ResidueSelection selection;
				evaluation::find_existing_residues( rmsd_pose, tag, selection );
				core::scoring::ResidueSelectionVector vector;
				std::copy( selection.begin(), selection.end(), std::back_inserter( vector ) );
				eval.add_evaluation( PoseEvaluatorOP( new simple_filters::TruncatedScoreEvaluator( tag, vector ) ) );
			}
			if ( option[ OptionKeys::evaluation::symmetric_rmsd ]() ) {
				eval.add_evaluation( PoseEvaluatorOP( new simple_filters::SymmetricRmsdEvaluator( rmsd_pose, tag ) ) );
			}
    }
  }


}

std::string RmsdTargetEvaluatorCreator::type_name() const {
	return "RmsdTargetEvaluatorCreator";
}

} //namespace
} //namespace
