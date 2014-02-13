// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/util.cc
/// @brief utility functions for silent-file classes.
/// @author James Thompson

// C++ Headers
#include <string>
#include <map>

// mini headers
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/types.hh>

// Core Headers
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/Residue.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/LazySilentFilePoseInputStream.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>

#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
#include <utility/file/gzip_util.hh>
#include <utility/file/file_sys_util.hh>

// option key includes

#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/boinc.OptionKeys.gen.hh>
#include <utility/vector1.hh>
#include <stdio.h>

namespace core {
namespace io {
namespace silent {

static basic::Tracer tr("core.io.silent");

void
gzip() {
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	// gzip output file if desired
	std::string const basename = option[ out::file::silent ]();
	if ( option[ out::silent_gz ]() ) {
		utility::vector1< std::string > file_list;
		file_list.push_back( basename );

		for ( utility::vector1< std::string >::const_iterator
				fn = file_list.begin(), end = file_list.end(); fn != end; ++fn
		) {
			utility::io::izstream in_stream( *fn );
			if ( in_stream ) {
				in_stream.close();
				in_stream.clear();
				tr.Info << "GZIP SILENT FILE: " << *fn << std::endl;
				utility::file::gzip( *fn, true );
			} // file exists
		} // loop over each file
	} // if ( option[ out::silent_gz ]() )
} // gzip


/////////////////////////////////////////////////////////////////
std::map< std::string, bool >
initialize_tag_is_done( std::string const & silent_file ){

	std::map< std::string, bool > tag_is_done;

	utility::vector1< std::string > tags_done;

	SilentFileData silent_file_data;
	if ( utility::file::file_exists( silent_file ) ) {
		tags_done = silent_file_data.read_tags_fast( silent_file );
		for ( utility::vector1< std::string >::const_iterator iter = tags_done.begin(); iter != tags_done.end(); iter++ ) {
			std::cout << "Already done: " << *iter << std::endl;
			tag_is_done[ *iter ] = true;
		}
	}

	return tag_is_done;
}

/////////////////////////////////////////////////////////////////
//@details Takes the output and filters to limit the number of structures.  This code's primary use is boinc where you don't want to generate more than 1 structure every 60 seconds. 
// ***This also does a core cut. You want all boinc servers to get the same amount of credit for work-done. So the PRIMARY CUT IS SCORE. The time cut is mandatory for VERY fast computers that
// generate more than 1 model every 60 seconds even with the score cut.  NOTE IF CENTOID be sure you use -in:file:centroid_input 
void boincOutputFilter(Real runTime, Real minTimePerModel){
	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using core::import_pose::pose_from_pdb;
	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using utility::vector1;
	using utility::file::file_exists;
	using std::string;
	string tempSilentLocation = option[boinc::score_cut_fl]();
	Real scoreCutPct = option[boinc::score_cut_pct]();
	vector1<Real> scores;
	//step1: read in poses and output them to temp file.
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	SilentStructOP localOutputSS(new core::io::silent::BinaryProteinSilentStruct);
	SilentFileData sfd;
	string outputFileName = option[out::file::silent]();
	if(option[out::silent_gz]()){
		string tmpOutputFileName = outputFileName + ".gz";
		if(!file_exists(tmpOutputFileName))
			utility_exit_with_message( "could not find " + tmpOutputFileName);
		utility::file::gunzip(tmpOutputFileName,true);
	}
	if(!file_exists(outputFileName))
		utility_exit_with_message( "could not find " + outputFileName);
	vector1< string> tempSilentVectorOrig;
	tempSilentVectorOrig.push_back(outputFileName);
	LazySilentFilePoseInputStreamOP tempSilentStreamOrig(
			new LazySilentFilePoseInputStream(tempSilentVectorOrig)
		);
	MetaPoseInputStream input;
	input.add_pose_input_stream(tempSilentStreamOrig);
	while(input.has_another_pose()){
		core::pose::PoseOP input_poseOP;
		input_poseOP = new core::pose::Pose();
		input.fill_pose(*input_poseOP,*rsd_set);
		std::string tag = core::pose::tag_from_pose(*input_poseOP);
		if("fa_standard" == rsd_set->name()){
			core::scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(TALARIS_2013 ));
			Real score = scorefxn->score(*input_poseOP);
			scores.push_back(score);
		}else{
			core::scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function("score3"));
			Real score = scorefxn->score(*input_poseOP);
			scores.push_back(score);
		}
		localOutputSS->fill_struct(*input_poseOP);
		localOutputSS->decoy_tag(tag);
		sfd.write_silent_struct(*localOutputSS,tempSilentLocation);
	}
	//sort scores
	std::sort(scores.begin(), scores.end());
	//figure out the score cut
	Size keepNumb = (Size)floor( ( (Real) scores.size() )*(1-scoreCutPct) );
	Size maxByTime = keepNumb;
	//in boinc maxtime should be set.  
	if(runTime!=-1)
		 maxByTime= (Size)floor(runTime/minTimePerModel);
	if(maxByTime < keepNumb)
		keepNumb = maxByTime;
	if(keepNumb < 1)
		keepNumb = 1;
	Real scoreThreshold = scores[keepNumb];


	//setup rd2 input stream
	vector1< string> tempSilentVector;
	tempSilentVector.push_back(tempSilentLocation);
	LazySilentFilePoseInputStreamOP tempSilentStream(
			new LazySilentFilePoseInputStream(tempSilentVector)
		);
	//wipe out input
	remove(outputFileName.c_str());
	remove(tempSilentLocation.c_str());
	//filter round 2
	MetaPoseInputStream input_rd2;
	input_rd2.add_pose_input_stream(tempSilentStream);
	while(input_rd2.has_another_pose()){
		core::pose::PoseOP input_poseOP;
		input_poseOP = new core::pose::Pose();
		input_rd2.fill_pose(*input_poseOP,*rsd_set);
		std::string tag = core::pose::tag_from_pose(*input_poseOP);
		Real tmpScore;
		if("fa_standard" == rsd_set->name()){
			core::scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function(TALARIS_2013 ));
			tmpScore = scorefxn->score(*input_poseOP);

		}else{
			core::scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function("score3"));
			tmpScore = scorefxn->score(*input_poseOP);
		}
		if(tmpScore<=scoreThreshold){
			localOutputSS->fill_struct(*input_poseOP);
			localOutputSS->decoy_tag(tag);
			sfd.write_silent_struct(*localOutputSS,outputFileName);
		}
	}
	if(option[out::silent_gz]())
		utility::file::gzip(outputFileName,true);
}

} // namespace silent
} // namespace io
} // namespace core
