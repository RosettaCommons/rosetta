// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/simple_filters/ChiWellRmsdEvaluatorCreator.hh
/// @brief  Header for ChiWellRmsdEvaluatorCreator
/// @author Matthew O'Meara

// Unit Headers
#include <protocols/simple_filters/ChiWellRmsdEvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/EvaluatorCreator.hh>

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/simple_filters/ChiWellRmsdEvaluator.hh>

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


static THREAD_LOCAL basic::Tracer tr( "protocols.evaluation.ChiWellRmsdEvaluatorCreator" );

namespace protocols {
namespace simple_filters {

ChiWellRmsdEvaluatorCreator::~ChiWellRmsdEvaluatorCreator() {}

void ChiWellRmsdEvaluatorCreator::add_evaluators( evaluation::MetaPoseEvaluator & eval ) const {
	using namespace core;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using protocols::evaluation::PoseEvaluatorOP;
	tr.Info << "Evaluation Creator active ... " << std::endl;
	if ( option[ OptionKeys::evaluation::chirmsd ].user() ) {

		core::pose::PoseOP native_pose = NULL;
		if ( option[ in::file::native ].user() ) {
			native_pose = core::pose::PoseOP( new core::pose::Pose );
			core::import_pose::pose_from_file( *native_pose, option[ in::file::native ]() , core::import_pose::PDB_file);
		}

		typedef utility::vector1< std::string > RmsdVector;
		RmsdVector const& rmsd( option[ OptionKeys::evaluation::chirmsd ]() );

		//the ultimate in RMSD technology:
		//you can have any number of tripletts:
		// each triplett contains  the
		//    target: file or special tags: NATIVE (use in:file:native ) IRMS (use structure in Job object)
		///   columnname: XXX creates: rms_XXX  gdtmm_XXX
		///   selection: file or special tag: FULL
		///   modifier: EXCLUDE ( take the inverse of the selection )
		///                INLINE r1 r2 r3 ... rm END_INLINE   -- select residues r1, r2, r3, ..., rm
		for ( RmsdVector::const_iterator it=rmsd.begin(); it!=rmsd.end(); ++it ) {
			core::pose::PoseOP target_pose = NULL;
			std::string fname( *it );
			std::string column;
			std::string selection_file;
			bool invert( false );
			core::Size nchi_max( 4 );
			core::Real sasa_max( 10000.0 );
			++it;
			while ( it != rmsd.end() && (*it).find('=')!=std::string::npos ) {
				Size pos( it->find('=') );
				std::string key=it->substr(0,pos);
				std::string value=it->substr(pos+1);
				if ( key=="nchi" ) {
					nchi_max=utility::string2int( value );
				} else if ( key =="sasa" ) {
					sasa_max=utility::string2float( value );
				} else {
					utility_exit_with_message( "key not recognized: "+key+" possible keys: { heavy, sasa }" );
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

				loops::PoseNumberedLoopFileReader reader;
				reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
				loops::SerializedLoopList loops = reader.read_pose_numbered_loops_file( is, selection_file, false /*no strict checking */ );
				loops::Loops core( loops );
				core.get_residues( selection );
			}
			if ( invert ) {
				utility::vector1< Size > inverted_selection;
				for ( Size i = 1; i<=target_pose->total_residue(); ++i ) {
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
				eval.add_evaluation( PoseEvaluatorOP( new simple_filters::ChiWellRmsdEvaluator( target_pose, nchi_max, sasa_max, selection, column) ) );
			} else {
				eval.add_evaluation( PoseEvaluatorOP( new simple_filters::ChiWellRmsdEvaluator( target_pose, nchi_max, sasa_max, column ) ) );
			}
		} // iterate over tripletts in option -rmsd
	} // option -rmsd
}

std::string ChiWellRmsdEvaluatorCreator::type_name() const {
	return "ChiWellRmsdEvaluatorCreator";
}

} //namespace
} //namespace
