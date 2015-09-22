// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/StructureSimilarityEvaluatorCreator.hh
/// @brief  Header for StructureSimilarityEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/StructureSimilarityEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/StructuralSimilarityEvaluator.hh>

// Project Headers
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/chemical/util.hh>

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
#include <utility/vector0.hh>

#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh>
#endif


static THREAD_LOCAL basic::Tracer tr( "protocols.evalution.StructureSimilarityEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

StructureSimilarityEvaluatorCreator::~StructureSimilarityEvaluatorCreator() {}

void StructureSimilarityEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::structural_similarity );

}

void StructureSimilarityEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;

	if ( option[ OptionKeys::evaluation::structural_similarity ].user() ) {
		using std::string;
		using core::import_pose::pose_stream::SilentFilePoseInputStream;

		SilentFilePoseInputStream silent_input(
			option[ OptionKeys::evaluation::structural_similarity ]()
		);
		core::chemical::ResidueTypeSetCOP rsd_set( core::chemical::rsd_set_from_cmd_line() );
		utility::vector1< core::pose::Pose > poses;
		while ( silent_input.has_another_pose() ) {
			core::pose::Pose pose;
			silent_input.fill_pose( pose, *rsd_set );
			poses.push_back( pose );
		}
		eval.add_evaluation( PoseEvaluatorOP( new StructuralSimilarityEvaluator(poses) ) );
	}

}

std::string StructureSimilarityEvaluatorCreator::type_name() const {
	return "StructureSimilarityEvaluatorCreator";
}

} //namespace
} //namespace
