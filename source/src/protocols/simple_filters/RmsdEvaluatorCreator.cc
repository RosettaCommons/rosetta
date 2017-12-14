// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/RmsdEvaluatorCreator.hh
/// @brief  Header for RmsdEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/RmsdEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>

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


static basic::Tracer tr( "protocols.evaluation.RmsdEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

RmsdEvaluatorCreator::~RmsdEvaluatorCreator() = default;

void RmsdEvaluatorCreator::register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;

	OPT( evaluation::rmsd );
	OPT( evaluation::gdtmm );
	OPT( in::file::native );

}

void RmsdEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;

	if ( option[ OptionKeys::evaluation::rmsd ].user() ) {

		core::pose::PoseOP native_pose = nullptr;
		if ( option[ in::file::native ].user() ) {
			native_pose = core::pose::PoseOP( new core::pose::Pose );
			core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
		}

		using RmsdVector = utility::vector1<std::string>;
		RmsdVector const& rmsd( option[ OptionKeys::evaluation::rmsd ]() );

		//the ultimate in RMSD technology:
		//you can have any number of tripletts:
		// each triplett contains  the
		//    target: file or special tags: NATIVE (use in:file:native ) IRMS (use structure in Job object)
		///   columnname: XXX creates: rms_XXX  gdtmm_XXX
		///   selection: file or special tag: FULL
		///   modifier: EXCLUDE ( take the inverse of the selection )
		///                INLINE r1 r2 r3 ... rm END_INLINE   -- select residues r1, r2, r3, ..., rm
		for ( auto it=rmsd.begin(); it!=rmsd.end(); ++it ) {
			core::pose::PoseOP target_pose = nullptr;
			std::string fname( *it );
			std::string column;
			std::string selection_file;
			bool CA( true );
			bool invert( false );
			bool loop_rms( false );
			bool superimpose_for_looprms( false );
			++it;
			loops::Loops core;
			while ( it != rmsd.end() && (*it).find('=')!=std::string::npos ) {
				Size pos( it->find('=') );
				std::string key=it->substr(0,pos);
				std::string value=it->substr(pos+1);
				if ( key=="heavy" && value=="yes" ) {
					tr.Info << "switching to full-atom RMSD " << std::endl;
					CA=false;
				} else if ( key == "loop" ) {
					loop_rms = true;
				} else if ( key == "superimpose" ) {
					superimpose_for_looprms = (value=="yes");
				} else if ( key == "core" ) {
					std::ifstream is( value.c_str() );

					if ( !is.good() ) {
						utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + value + "'" );
					}
					loops::PoseNumberedLoopFileReader reader;
					reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
					loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file(is, value, false );
					core = loops::Loops( loops );
				} else {
					utility_exit_with_message( "key not recognized: "+key+" possible keys: { heavy }" );
				}
				tr.Info << "detected option: " << key << " : " << value << std::endl;
				++it;
			}
			if ( it != rmsd.end() ) {
				column = *it;
			} else {
				utility_exit_with_message(
					"need to specify tripletts <target> [key/value pairs] <column> <selection/FULL> with option -evaluation:rmsd   last read: "+fname );
			}
			++it;
			if ( it != rmsd.end() ) {
				selection_file = *it;
			} else {
				utility_exit_with_message(
					"need to specify tripletts <target> <column> <selection/FULL> with option -evaluation:rmsd   last read: "+column );
			}
			if ( fname == "NATIVE" ) target_pose = native_pose;
			else if ( fname != "IRMS" ) {
				target_pose = core::pose::PoseOP( new pose::Pose );
				core::import_pose::pose_from_file( *target_pose, fname , core::import_pose::PDB_file);
			}

			if ( selection_file == "EXCLUDE" ) {
				invert = true;
				++it;
				if ( it != rmsd.end() ) {
					selection_file = *it;  //read next tag
				} else {
					utility_exit_with_message(
						"need to specify a <selection/FULL> after 'EXCLUDE' with option -evaluation:rmsd    last read: "+column );
				}
			}
			loops::Loops loops;
			utility::vector1< Size> selection; //figure out selection
			if ( selection_file == "INLINE" ) {
				std::string next_tag( "" );
				++it;
				if ( it != rmsd.end() ) {
					next_tag = *it;
				} else {
					utility_exit_with_message(
						"need to find END_INLINE after INLINE in option -evaluation:rmsd    last read: "+column );
				} //error condition
				while ( next_tag != "END_INLINE" ) {
					selection.push_back( utility::string2int( next_tag ) );
					++it;
					if ( it != rmsd.end() ) {
						next_tag = *it;
					} else {
						utility_exit_with_message(
							"need to find END_INLINE after INLINE in option -evaluation:rmsd    last read: "+column );
					} //error condition
				}
			} else if ( selection_file != "FULL" ) {
				std::ifstream is( selection_file.c_str() );

				if ( !is.good() ) {
					utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + selection_file + "'" );
				}
				loops::SerializedLoopList list_of_loops;
				if ( loop_rms ) {
					loops::PoseNumberedLoopFileReader reader;
					list_of_loops = reader.read_pose_numbered_loops_file(is, selection_file, false );
					loops = loops::Loops( list_of_loops );
				} else {
					loops::PoseNumberedLoopFileReader reader;
					reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
					list_of_loops = reader.read_pose_numbered_loops_file(is, selection_file, false );
					loops::Loops core = loops::Loops( list_of_loops );
					core.get_residues( selection );
				}
			}
			if ( invert ) {
				utility::vector1< Size > inverted_selection;
				for ( Size i = 1; i<=target_pose->size(); ++i ) {
					bool found( false );
					for ( Size j = 1; j<=selection.size() && !found; ++j ) {
						if ( selection[ j ]==i ) {
							found = true;
						}
					}
					if ( !found ) inverted_selection.push_back( i );
				}
				selection = inverted_selection;
			}
			if ( selection_file != "FULL" ) {
				if ( !loop_rms ) {
					eval.add_evaluation( PoseEvaluatorOP( new simple_filters::SelectRmsdEvaluator( target_pose, selection, column, CA ) ) );
				} else {
					eval.add_evaluation( PoseEvaluatorOP( new simple_filters::LoopRmsdEvaluator( target_pose, loops, core, column, CA, superimpose_for_looprms ) ) );
				}
				if ( option[ OptionKeys::evaluation::gdtmm ]() ) {
					eval.add_evaluation( PoseEvaluatorOP( new simple_filters::SelectGdtEvaluator( target_pose, selection, column) ) );
				}
			} else {
				eval.add_evaluation( PoseEvaluatorOP( new simple_filters::SelectRmsdEvaluator( target_pose, column, CA ) ) );
				if ( option[ OptionKeys::evaluation::gdtmm ]() ) {
					eval.add_evaluation( PoseEvaluatorOP( new simple_filters::SelectGdtEvaluator( target_pose, column ) ) );
				}
			} // no selection
		} // iterate over tripletts in option -rmsd
	} // option -rmsd

}

std::string RmsdEvaluatorCreator::type_name() const {
	return "RmsdEvaluatorCreator";
}

} //namespace
} //namespace
