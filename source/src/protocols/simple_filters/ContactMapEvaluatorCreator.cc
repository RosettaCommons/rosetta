// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/ContactMapEvaluatorCreator.hh
/// @brief  Header for ContactMapEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/ContactMapEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>

#include <protocols/simple_filters/ContactMapEvaluator.hh>

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

//// C++ headers

// due to template function
#include <core/io/silent/SilentStruct.hh>


// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>

#ifdef WIN32
	#include <core/scoring/constraints/Constraint.hh>
#endif


static thread_local basic::Tracer tr( "protocols.evalution.ContactMapEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

ContactMapEvaluatorCreator::~ContactMapEvaluatorCreator() {}

void ContactMapEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::contact_map );
	OPT( in::file::native );
}

void ContactMapEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;


	if ( option[ OptionKeys::evaluation::contact_map ] ) {
		core::pose::PoseOP native_pose = NULL;
		if ( option[ in::file::native ].user() ) {
			native_pose = core::pose::PoseOP( new core::pose::Pose );
			core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]() );
		}

		if ( !native_pose ) {
			tr.Error << "Error: -evaluation::contact_map must be specified with a native!\n";
		} else {
			core::Real max_dist(12);
			core::Size min_seqsep(5);
			eval.add_evaluation(
				PoseEvaluatorOP( new ContactMapEvaluator( *native_pose, max_dist, min_seqsep ) )
			);
		}
	}


}

std::string ContactMapEvaluatorCreator::type_name() const {
	return "ContactMapEvaluatorCreator";
}

} //namespace
} //namespace
