// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/comparative_modeling/AlignRmsdTargetEvaluatorCreator.hh
/// @brief  Header for AlignRmsdTargetEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/comparative_modeling/AlignRmsdTargetEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/comparative_modeling/Align_RmsdEvaluator.hh>
#include <protocols/comparative_modeling/Align_RotamerEvaluator.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.hh>

#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>

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


static basic::Tracer tr( "protocols.comparative_modeling.AlignRmsdTargetEvaluatorCreator" );

namespace protocols {
namespace comparative_modeling {

AlignRmsdTargetEvaluatorCreator::~AlignRmsdTargetEvaluatorCreator() = default;

void AlignRmsdTargetEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::align_rmsd_target );
	OPT( evaluation::align_rmsd_column );
	OPT( evaluation::align_rmsd_fns );
	OPT( evaluation::align_rmsd_format );

}

void AlignRmsdTargetEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;

	if ( option[ OptionKeys::evaluation::align_rmsd_target ].user() ) {
		using std::string;
		using utility::vector1;
		vector1< string > const & align_rmsd_target(
			option[ OptionKeys::evaluation::align_rmsd_target ]()
		);
		vector1< string > const & align_rmsd_col_names(
			option[ OptionKeys::evaluation::align_rmsd_column ]()
		);

		vector1< string > align_rmsd_fns;
		if ( option[ OptionKeys::evaluation::align_rmsd_fns ].user() ) {
			align_rmsd_fns = option[ OptionKeys::evaluation::align_rmsd_fns ]();
		}
		runtime_assert( align_rmsd_target.size() == align_rmsd_col_names.size() );

		// evaluation::gdtmm will conflict, but lets override
		bool gdt_by_TM( option[ OptionKeys::evaluation::gdttm ]() );

		for ( Size ii = 1; ii <= align_rmsd_target.size(); ++ii ) {
			pose::PoseOP rmsd_pose( new pose::Pose );
			core::import_pose::pose_from_file( *rmsd_pose, align_rmsd_target[ii] , core::import_pose::PDB_file);
			//string const tag( align_rmsd_target[ii] );
			string const tag( align_rmsd_col_names[ii] );
			core::sequence::SequenceAlignmentOP aln(nullptr);
			if ( align_rmsd_fns.size() >= ii ) {
				*aln = core::sequence::read_aln(
					align_rmsd_fns[ii], option[ OptionKeys::evaluation::align_rmsd_format ]()
					).front();
			}
			eval.add_evaluation( PoseEvaluatorOP( new Align_RmsdEvaluator(rmsd_pose,tag,true,aln,gdt_by_TM) ) );
			eval.add_evaluation( PoseEvaluatorOP( new Align_RotamerEvaluator(rmsd_pose,tag,true,aln) ) );
		}
	}

}

std::string AlignRmsdTargetEvaluatorCreator::type_name() const {
	return "AlignRmsdTargetEvaluator";
}

} //namespace
} //namespace
