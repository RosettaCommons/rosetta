// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/brunette/unalignedEvaluate.cc
///
/// @brief  This takes an alignment file, fasta file, and a silent file.  It then evaluates the qualities of the loops.

/// @usage: -in:file:alignent <alignment> -in:file:s <decoy.out> -in:file:native <native pose> -in:file:fasta <fasta> -unalignedEvaluate:pair_silent_all_aln -cm:aln_format grishin -in:file:template_pdb <templates>

/// @author TJ Brunette
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <basic/Tracer.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/id/SequenceMapping.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>


#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/Sequence.hh>
#include <core/sequence/SequenceAlignment.hh>

#include <core/types.hh>

#include <devel/init.hh>
// option key includes
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <utility/options/keys/BooleanOptionKey.fwd.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/io/ozstream.hh>

#include <protocols/comparative_modeling/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <numeric/xyzVector.hh>

#include <apps/pilot/brunette/tj_util.hh>
#include <map>
#include <iostream>
#include <fstream>

static THREAD_LOCAL basic::Tracer tr( "filterAlignments" );
using std::string;
using std::map;
using std::multimap;
using std::vector;
using core::pose::Pose;
using protocols::loops::Loops;
using namespace core::sequence;
using namespace core::id;


SequenceAlignment filterNonCore(SequenceAlignment inptAln, SequenceAlignment partialThreadAln, Size startRes, Size endRes, Pose partialThreadPose, vector1<core::pose::PoseOP> poses){
	Real PERCENT_RES_OVERLAP = .70;
 	Real MIN_RMSD_MATCH = 2;
	Real PERCENT_TOP_LOOPS = .05;
	const Size nres1 = inptAln.sequence(1)->ungapped_sequence().size();
	const Size nres2 = inptAln.sequence(2)->ungapped_sequence().size();
	SequenceOP origAlnSequence1_ = new Sequence(inptAln.sequence(1)->ungapped_sequence(),inptAln.sequence(1)->id(),inptAln.sequence(1)->start());
	SequenceOP origAlnSequence2_ = new Sequence(inptAln.sequence(2)->ungapped_sequence(),inptAln.sequence(2)->id(),inptAln.sequence(2)->start());
	SequenceMapping inptMapping( inptAln.sequence_mapping( 1, 2 ) );
	SequenceMapping partialThreadMapping(partialThreadAln.sequence_mapping(1,2));
	//step 1 see if alignment has atleast PERCENT_RES_OVERLAP overlap with loop
	Size loopLength = (endRes-startRes)+1;
	Real reqResOverlap = floor((Real)loopLength*PERCENT_RES_OVERLAP);
	Size definedRes = 0;
	for ( Size ii = startRes; ii <= endRes; ++ii ) {
		std::cout <<"inptMapping" <<  ii <<"-" << inptMapping[ii] << std::endl;
		std::cout <<"partialThreadMap" << ii << "-" << partialThreadMapping[ii] << std::endl;
		if(inptMapping[ii] != 0)
			definedRes ++;
	}
	std::cout << "definedRes" << definedRes << "req" << reqResOverlap <<"length" << loopLength << std::endl;
	if(definedRes < reqResOverlap){
		SequenceMapping filtMapping(nres1,nres2);
		SequenceAlignment filtAln = mapping_to_alignment(filtMapping,origAlnSequence1_,origAlnSequence2_);
		return(filtAln);
	}
	//step 2 select top loops
	multimap<Real,core::pose::PoseOP> rankedPoses;
	for(vector1<core::pose::PoseOP>::iterator pose_itr = poses.begin(); pose_itr != poses.end(); ++pose_itr){
		Real tmp_score = region_score12(**pose_itr,startRes,endRes);
		rankedPoses.insert(std::pair<Real,core::pose::PoseOP>(tmp_score,*pose_itr));
	}
	//step 3 see if largest contiguous chunk within partial thread and loop matches the residues < MIN_RMSD_MATCH.  If yes, copy over aligned residues  to output alignment
	Size numbTopPoses = (Size)round(poses.size()*PERCENT_TOP_LOOPS);
	if(numbTopPoses == 0)
		numbTopPoses = 1;
	//get maximal overlap from start or end of loop
	Size currentOverlapStart = startRes;
	Size currentOverlapEnd = startRes;
	Size currentOverlapLength = 0;
	Size maxOverlapStart = startRes;
	Size maxOverlapEnd = startRes;
	Size maxOverlapLength = 0;
	bool previousResInOverlap = false;
	for( Size ii = startRes; ii <= endRes; ++ii ){
		if(inptMapping[ii]!=0){
			if(!previousResInOverlap){
				currentOverlapStart = ii;
				currentOverlapEnd = ii;
				currentOverlapLength  = 1;
			}
			else{
				currentOverlapEnd = ii;
				currentOverlapLength = (currentOverlapEnd-currentOverlapStart)+1;
			}
			if(currentOverlapLength>maxOverlapLength){
				maxOverlapLength = currentOverlapLength;
				maxOverlapStart = currentOverlapStart;
				maxOverlapEnd = currentOverlapEnd;
			}
			previousResInOverlap = true;
		}
		else
			previousResInOverlap = false;
	}

	std::cout << "start" << maxOverlapStart << "end" << maxOverlapEnd <<"," << maxOverlapLength << std::endl;
	Size poseCount = 0;
	bool goodLoopExists = false;
	//get correct aln
	for(multimap<Real,core::pose::PoseOP>::iterator score_itr=rankedPoses.begin(); score_itr != rankedPoses.end() &&(poseCount<numbTopPoses) ; ++score_itr){
		poseCount++;
		Real rmsd = region_rmsd(*score_itr->second,partialThreadPose,partialThreadAln,maxOverlapStart,maxOverlapEnd);
		if(rmsd<MIN_RMSD_MATCH)
			goodLoopExists = true;
	}
	if(goodLoopExists == true){//allow alignment to pass through
		SequenceMapping filtMapping(nres1,nres2);
		for ( Size ii = startRes; ii <= endRes; ++ii ) {
			if(inptMapping[ii] != 0)
				filtMapping[ii] = inptMapping[ii];
		}
		SequenceAlignment filtAln = mapping_to_alignment(filtMapping,origAlnSequence1_,origAlnSequence2_);
		return(filtAln);
	}
	else{
		SequenceMapping filtMapping(nres1,nres2);
		SequenceAlignment filtAln = mapping_to_alignment(filtMapping,origAlnSequence1_,origAlnSequence2_);
		return(filtAln);
	}
}

SequenceAlignment filterCore(SequenceAlignment inptAln, SequenceAlignment coreAln){
	SequenceMapping inptMapping( inptAln.sequence_mapping( 1, 2 ) );
	SequenceMapping coreMapping( coreAln.sequence_mapping( 1, 2 ) );
	const Size nres1 = inptAln.sequence(1)->ungapped_sequence().size();
	const Size nres2 = inptAln.sequence(2)->ungapped_sequence().size();
	SequenceMapping filtMapping(nres1,nres2);
	for ( Size ii = 1; ii <= nres1; ++ii ) {
		if(coreMapping[ii] != 0)
			filtMapping[ii] = inptMapping[ii];
	}
	SequenceOP origAlnSequence1_ = new Sequence(inptAln.sequence(1)->ungapped_sequence(),inptAln.sequence(1)->id(),inptAln.sequence(1)->start());
	SequenceOP origAlnSequence2_ = new Sequence(inptAln.sequence(2)->ungapped_sequence(),inptAln.sequence(2)->id(),inptAln.sequence(2)->start());
	SequenceAlignment filtAln = mapping_to_alignment(filtMapping,origAlnSequence1_,origAlnSequence2_);
 	return(filtAln);
}

SequenceAlignment combineIdenticalLengthAln(vector1<SequenceAlignment> subsetAlns){
	SequenceAlignment firstAln = subsetAlns[1];
	const Size nres1 = firstAln.sequence(1)->ungapped_sequence().size();
	const Size nres2 = firstAln.sequence(2)->ungapped_sequence().size();
	SequenceOP origAlnSequence1_ = new Sequence(firstAln.sequence(1)->ungapped_sequence(),firstAln.sequence(1)->id(),firstAln.sequence(1)->start());
	SequenceOP origAlnSequence2_ = new Sequence(firstAln.sequence(2)->ungapped_sequence(),firstAln.sequence(2)->id(),firstAln.sequence(2)->start());
	SequenceMapping finalMapping(nres1,nres2);
	vector1<SequenceMapping> subsetMappings;
	for(Size jj=1; jj<=subsetAlns.size(); ++jj){
		SequenceMapping tempMapping( subsetAlns[jj].sequence_mapping( 1, 2 ) );
		subsetMappings.push_back(tempMapping);
	}
	for ( Size ii = 1; ii <= nres1; ++ii ) {
		for( Size jj= 1; jj <= subsetMappings.size(); ++jj){
			if(subsetMappings[jj][ii] != 0)
				finalMapping[ii] = subsetMappings[jj][ii];
		}
	}
	SequenceAlignment finalAln = mapping_to_alignment(finalMapping,origAlnSequence1_,origAlnSequence2_);
	return finalAln;
}

SequenceAlignment filterAlignment(SequenceAlignment inptAln, SequenceAlignment coreAln,protocols::loops::LoopsOP unconvergedLoops,Pose partialThreadPose,vector1<core::pose::PoseOP> poses){
	std::cout << "filterAlignmentA" << std::endl;
	vector1<SequenceAlignment> subsetAlns;
	//Evaluate each loop then conserved region then combine them.

	SequenceOP fastaSequence = new Sequence(*poses[1]);
	SequenceOP partialThreadSequence = new Sequence(partialThreadPose);
	//SequenceAlignment partialThreadAln =  align_naive(partialThreadSequence,fastaSequence);
	SequenceMapping inptMapping( inptAln.sequence_mapping( 1, 2 ) );
	const Size nres1 = partialThreadSequence->ungapped_sequence().size();
	const Size nres2 = fastaSequence->ungapped_sequence().size();
	SequenceMapping partialThreadMap(nres1,nres2);
	Size partialThreadPos= 1;
	for(int ii=1; ii<=nres2; ++ii){
		if(inptMapping[ii] != 0){
			partialThreadMap[partialThreadPos] = ii;
			partialThreadPos++;
		}
	}
	SequenceAlignment partialThreadAln = mapping_to_alignment(partialThreadMap,partialThreadSequence,fastaSequence);
	for(Loops::const_iterator itr_loop = unconvergedLoops->begin(); itr_loop != unconvergedLoops->end(); ++itr_loop){
		std::cout << "filterAlignmentA" << std::endl;
		Size start_res = itr_loop->start();
		Size end_res = itr_loop->stop();
		std::cout << "start_res" <<start_res <<"end_res" << end_res << std::endl;
		SequenceAlignment tmpAln = filterNonCore(inptAln,partialThreadAln,start_res,end_res, partialThreadPose, poses);
		subsetAlns.push_back(tmpAln);
	}
	std::cout << "filterCore" << std::endl;
	SequenceAlignment tmpAln = filterCore(inptAln,coreAln);
	subsetAlns.push_back(tmpAln);
	std::cout << "combiningaln" << std::endl;
	SequenceAlignment combinedAln = combineIdenticalLengthAln(subsetAlns);
	std::cout << "allcombined" << std::endl;
	return(combinedAln);
}

vector1 <SequenceAlignment> filterAllAlignments(map<string,SequenceAlignment> alignDataMapped,SequenceAlignment coreAln,protocols::loops::LoopsOP unconvergedLoops ,map< string, Pose> partialThreadsMapped, vector1<core::pose::PoseOP> poses){
	vector1 <SequenceAlignment> outAlns;
	map<string,SequenceAlignment>::iterator itr;
	for(map<string,SequenceAlignment>::iterator itr = alignDataMapped.begin(); itr != alignDataMapped.end(); ++itr){
		std::cout << "*******filtering " << itr->second << std::endl;
		SequenceAlignment tmp_aln = filterAlignment(itr->second,coreAln,unconvergedLoops,partialThreadsMapped[itr->first],poses);
		outAlns.push_back(tmp_aln);
	}
	return(outAlns);
}

void output_alignments(vector1 <SequenceAlignment> alns, std::ostream & out){
  for ( Size ii = 1; ii <= alns.size(); ++ii ) {
    alns[ii].printGrishinFormat(out);
  }
}

namespace filterAlignments {
basic::options::FileOptionKey core_pdb("filterAlignments:core_pdb");
}

int main( int argc, char * argv [] ) {
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::sequence;
	using namespace core::id;
	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using namespace protocols::comparative_modeling;
	using core::import_pose::pose_from_pdb;
	using core::sequence::read_fasta_file;
	using utility::file_basename;
	option.add (filterAlignments::core_pdb,"core pdb");
	devel::init(argc, argv);
	string query_sequence (
			read_fasta_file( option[ in::file::fasta ]()[1])[1]->sequence()	);
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	//get unconverged loops
	Pose core_pose;
	pose_from_pdb(
			core_pose,
			*rsd_set,
			option[filterAlignments::core_pdb]()
			);
	SequenceOP core_sequence = new Sequence(core_pose);
	SequenceOP fasta_sequence = new Sequence(query_sequence,"fasta",1);
	SequenceAlignment coreAln =  align_naive( fasta_sequence, core_sequence);
	protocols::loops::LoopsOP unconvergedLoops = loops_from_alignment(query_sequence.size(),coreAln,1);
	//get partial threads
	map< string, Pose > templateData = poses_from_cmd_line(
			option[ in::file::template_pdb ]());
	map<string,SequenceAlignment> alnDataMapped = input_alignmentsMapped();
	map< string, Pose> partialThreadsMapped = generate_partial_threads(alnDataMapped,templateData,query_sequence,false);
	//create vector of input poses.
	MetaPoseInputStream input = streams_from_cmd_line();
	vector1<core::pose::PoseOP> poses;
	while(input.has_another_pose()){
			core::pose::PoseOP input_poseOP;
			input_poseOP = new core::pose::Pose();
			input.fill_pose(*input_poseOP,*rsd_set);
			poses.push_back(input_poseOP);
	}
	//filter alignments
	vector1 <SequenceAlignment> out_alns = filterAllAlignments(alnDataMapped,coreAln,unconvergedLoops,partialThreadsMapped,poses);
	//output alignments
	string out_filename = option[out::file::alignment ]();
	std::ofstream out_aln_stream( out_filename.c_str() );
	output_alignments(out_alns,out_aln_stream);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

