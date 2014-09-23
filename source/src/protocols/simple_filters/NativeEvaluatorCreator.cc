// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/simple_filters/NativeEvaluatorCreator.hh
/// @brief  Header for NativeEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/NativeEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/simple_filters/ScoreEvaluator.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>

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

// due to template function
#include <core/io/silent/SilentStruct.hh>

// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>

#ifdef WIN32
	#include <core/scoring/constraints/Constraint.hh>
#endif

#include <utility/exit.hh>
#include <fstream>

static thread_local basic::Tracer tr( "protocols.evalution.NativeEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

NativeEvaluatorCreator::~NativeEvaluatorCreator() {}

void NativeEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( in::file::native );
	OPT( in::file::native_exclude_res );
	OPT( evaluation::gdtmm );
	OPT( abinitio::rmsd_residues );

}

void NativeEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;

	core::pose::PoseOP native_pose = NULL;
	if ( option[ in::file::native ].user() ) {
		native_pose = core::pose::PoseOP( new core::pose::Pose );
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]() );
	}


	// set rmsd native
	if ( native_pose && option[ in::file::native ].user() ) {
		if ( option[ OptionKeys::symmetry::symmetric_rmsd ].user() ){
				eval.add_evaluation( PoseEvaluatorOP( new SymmetricRmsdEvaluator( native_pose, "" ) ) );
		}

		if ( option[ in::file::native_exclude_res ].user() ) {
			eval.add_evaluation( PoseEvaluatorOP( new SelectRmsdEvaluator(
				native_pose,
				core::scoring::invert_exclude_residues(
				native_pose->total_residue(), option[ in::file::native_exclude_res ]()), "" ) )
			);
			if ( option[ OptionKeys::evaluation::gdtmm ]() ) {
				eval.add_evaluation( PoseEvaluatorOP( new SelectGdtEvaluator(
					native_pose,
					core::scoring::invert_exclude_residues( native_pose->total_residue(),
					option[ in::file::native_exclude_res ]()), "" ) )
				);
			}
		} else if ( option[ OptionKeys::abinitio::rmsd_residues ].user() ){
			core::Size start = option[ OptionKeys::abinitio::rmsd_residues ]()[ 1 ];
			Size end = option[ OptionKeys::abinitio::rmsd_residues ]()[ 2 ];
			eval.add_evaluation( PoseEvaluatorOP( new RmsdEvaluator( native_pose, start, end,  "", option[ OptionKeys::abinitio::bGDT ]() ) ) );
		} else {
			eval.add_evaluation( PoseEvaluatorOP( new SelectRmsdEvaluator( native_pose, "" ) ) );
			if ( option[ OptionKeys::evaluation::gdtmm ]() ) eval.add_evaluation( PoseEvaluatorOP( new SelectGdtEvaluator( native_pose, "" ) ) );
			eval.add_evaluation( PoseEvaluatorOP( new SelectMaxsubEvaluator( native_pose, "" ) ) );
		}
	} // if ( native_pose_ )
	
	if ( option[ OptionKeys::evaluation::rmsd_select ].user() ) {
		utility::vector1< utility::file::FileName > const& rmsd_core( option[ OptionKeys::evaluation::rmsd_select ]() );
		if ( !option[ in::file::native ].user() ) utility_exit_with_message( "need to specify in:file:native together with rmsd_select " );
		
		for ( Size ct = 1; ct <= rmsd_core.size(); ct ++ ) {
			std::ifstream is( rmsd_core[ ct ].name().c_str() );
			
			if (!is.good()) {
				utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + rmsd_core[ ct ].name() + "'" );
			}
			
			loops::PoseNumberedLoopFileReader reader;
			reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
			loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file( is, rmsd_core[ ct ].name(), false /*no strict checking */ );

			loops::Loops core( loops );
			utility::vector1< Size> selection;
			core.get_residues( selection );
			if ( native_pose ) eval.add_evaluation( PoseEvaluatorOP( new simple_filters::SelectRmsdEvaluator( native_pose, selection, rmsd_core[ ct ].base() ) ) );
			if ( option[ OptionKeys::evaluation::score_with_rmsd ]() ) {
				eval.add_evaluation( PoseEvaluatorOP( new simple_filters::TruncatedScoreEvaluator( rmsd_core[ ct ].base(), selection ) ) );
			}
		}
	}


}

std::string NativeEvaluatorCreator::type_name() const {
	return "NativeEvaluatorCreator";
}

} //namespace
} //namespace
