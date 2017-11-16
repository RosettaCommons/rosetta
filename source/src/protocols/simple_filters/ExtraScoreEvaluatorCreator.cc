// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/ExtraScoreEvaluatorCreator.hh
/// @brief  Header for ExtraScoreEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/ExtraScoreEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/ScoreEvaluator.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>
#include <core/scoring/constraints/util.hh>

#include <core/chemical/ResidueType.hh>

#include <core/io/silent/silent.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


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


static basic::Tracer tr( "protocols.evaluation.ExtraScoreEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

ExtraScoreEvaluatorCreator::~ExtraScoreEvaluatorCreator() = default;

void ExtraScoreEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::extra_score );
	OPT( evaluation::extra_score_column );
	OPT( evaluation::extra_score_patch );
	OPT( evaluation::extra_score_select );

}

void ExtraScoreEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;


	if ( option[ OptionKeys::evaluation::extra_score ].user() ) {
		using namespace core::scoring;
		utility::vector1< std::string > const& extra_scores( option[ OptionKeys::evaluation::extra_score ]() );
		utility::vector1< std::string > const& extra_score_names( option[ OptionKeys::evaluation::extra_score_column]() );
		if ( extra_scores.size() != extra_score_names.size() ) {
			utility_exit_with_message("-extra_score: you need to provide as much extra_score_names as extra_scores! ");
		}
		for ( Size ct = 1; ct <= extra_scores.size(); ct ++ ) {
			std::string const& tag = extra_score_names[ ct ];
			std::string patch( "NOPATCH" );
			if ( option[ OptionKeys::evaluation::extra_score_patch ].user() ) {
				if ( option[ OptionKeys::evaluation::extra_score_patch ]().size() != extra_scores.size() ) {
					utility_exit_with_message("-extra_score: you need to provide as much extra_score_patch(es) as "
						"extra_scores! use NOPATCH as placeholder");
				}
				patch = option[ OptionKeys::evaluation::extra_score_patch ]()[ ct ];
			}
			ScoreFunctionOP scfxn( nullptr );

			if ( patch != "NOPATCH" ) {
				scfxn = ScoreFunctionFactory::create_score_function( extra_scores[ ct ], patch );
			} else {
				scfxn = ScoreFunctionFactory::create_score_function( extra_scores[ ct ] );
			}

			std::string name( extra_scores[ ct ] );
			if ( (name == "score0") ||
					(name == "score2") ||
					(name == "score3") ||
					(name == "score5") ) {
				core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *scfxn );
			} else {
				core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *scfxn );
			}

			std::string select_string( "SELECT_ALL" );
			if ( option[ OptionKeys::evaluation::extra_score_select ].user() ) {
				if ( option[ OptionKeys::evaluation::extra_score_select ]().size() != extra_scores.size() ) {
					utility_exit_with_message("-extra_score: you need to provide as much extra_score_patch(es) as "
						"extra_scores! use SELECT_ALL as placeholder");
				}
				select_string = option[ OptionKeys::evaluation::extra_score_select ]()[ ct ];
			}
			if ( select_string != "SELECT_ALL" ) {
				std::ifstream is( select_string.c_str() );

				if ( !is.good() ) {
					utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + select_string + "'" );
				}

				loops::PoseNumberedLoopFileReader reader;
				reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
				loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file( is, select_string, false /*no strict checking */ );
				loops::Loops core( loops );

				utility::vector1< Size> selection;
				core.get_residues( selection );
				eval.add_evaluation( PoseEvaluatorOP( new simple_filters::TruncatedScoreEvaluator( tag, selection, scfxn ) ) );
			} else {
				eval.add_evaluation( PoseEvaluatorOP( new ScoreEvaluator( tag, scfxn ) ) );
			}
		}
	}


}

std::string ExtraScoreEvaluatorCreator::type_name() const {
	return "ExtraScoreEvaluatorCreator";
}

} //namespace
} //namespace
