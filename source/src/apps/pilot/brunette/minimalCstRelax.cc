// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.AlignmentCluster

/// @file   apps/pilot/brunette/minimalCstRelax.cc
///
/// @brief  For the first structure in the alignment compare the deviation between native and the relaxed native structure. Then constraints are modified and coordinate constraints are added. If an alignment is given only the residues from the beginning to end of the alignment are used.

/// @usage: -in:file:s <pdb files> [options: -minimalCstRelax:coordinate_cst_gap in first round of relax gap between CA coordinate constraints] -in::file::alignmen [alignment.filt file]
/// @author TJ Brunette


#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <core/types.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/MissingDensityToJumpMover.hh>
#include <protocols/simple_moves/MissingDensityToJumpMover.fwd.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <basic/Tracer.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>


#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/cst_util.hh>


//utilities
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/options/keys/FileOptionKey.fwd.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/FileVectorOptionKey.fwd.hh>
#include <utility/options/keys/FileVectorOptionKey.hh>


#include <devel/init.hh>
//#include <devel/init.hh>
// option key includes
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>

#include <ObjexxFCL/format.hh>

#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include <devel/cstEnergyBalance/minimalCstRelaxUtil.hh>

basic::Tracer tr( "minimalCstRelax" );

using core::Size;
using std::map;
using std::string;
using namespace core::sequence;

namespace minimalCstRelax {
	basic::options::IntegerOptionKey coordinate_cst_gap("minimalCstRelax:coordinate_cst_gap");
	basic::options::BooleanOptionKey relax_pdb("minimalCstRelax:relax_pdb");
}

///@brief: gets the pose to relax.Returns either the input pose or if alignments exists a  subset of the input pose
void get_poseToRelax(core::pose::PoseOP & toRelax_poseOP, core::pose::PoseOP & input_poseOP,  core::import_pose::pose_stream::MetaPoseInputStream input, std::map<string, SequenceAlignment> alns,Size & offset){
	using namespace core::chemical;
	using utility::file_basename;
	using namespace devel::cstEnergyBalance;
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	input_poseOP = new core::pose::Pose();
	input.fill_pose(*input_poseOP,*rsd_set);
	toRelax_poseOP = input_poseOP;
	string const pose_tag( file_basename( core::pose::tag_from_pose( *input_poseOP ) ) );
	string pdbid = pose_tag.substr(0,pose_tag.length()-4);
	if (alns.size() == 0){//the case where there is no alignment file
		offset = 0;
	}
	else{
		Size firstRes,lastRes;
		if(alns.find(pdbid) == alns.end())
			utility_exit_with_message( "[ERROR] Unable to find alignment for " + pdbid);
		SequenceAlignment aln = alns.find(pdbid)->second;
		get_terminal_aln_res(aln,2,firstRes,lastRes);
		toRelax_poseOP = new core::pose::Pose(*input_poseOP,firstRes,lastRes);
		offset = firstRes-1;
	}
}

///@brief: outputs the ca atoms to constrain either with taking account the offset
void output_caAtomsToConstraint(const std::set< Size> caAtomsToConstrain,const Size offset,const string coordCstFile, const core::pose::PoseOP input_poseOP){
	using namespace devel::cstEnergyBalance;
	std::set < Size > adjusted_caAtomsToConstrain;
	std::set < Size >::iterator caAtomsToConstrain_iter;
	std::ofstream out( coordCstFile.c_str() );
	if (offset == 0){
		output_coordCsts(caAtomsToConstrain,out,*input_poseOP);
	}
	else{
		caAtomsToConstrain_iter = caAtomsToConstrain.begin();
		while(caAtomsToConstrain_iter != caAtomsToConstrain.end()){
			adjusted_caAtomsToConstrain.insert(*caAtomsToConstrain_iter+offset);
			caAtomsToConstrain_iter++;
		}
		output_coordCsts(adjusted_caAtomsToConstrain,out,*input_poseOP);
	}
	out.close();
}

int main( int argc, char * argv [] ) {
	try {

	using namespace protocols;
	using namespace relax;
	using namespace devel::cstEnergyBalance;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::io::silent;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::import_pose::pose_stream;
	using namespace protocols::moves;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
 	using utility::file_basename;
	option.add (minimalCstRelax::coordinate_cst_gap, "gap between residues with a coordinate constraint").def(4);
	option.add (minimalCstRelax::relax_pdb, "should a final relax be completed").def(false);
	devel::init(argc, argv);
	const Real COORDINATE_CST_WT = 1;
	SilentFileData sfd;
	SilentStructOP ss(new core::io::silent::ScoreFileSilentStruct);
	//Size gapBtwConstrainedResidues = option[ minimalCstRelax::coordinate_cst_gap]();  // unused ~Labonte
	//Get input pdbs---------------------------------------------------------------
	MetaPoseInputStream input = streams_from_cmd_line();
	std::map<string,SequenceAlignment> alns = input_alignmentsMapped(true);
	while(input.has_another_pose()){
		core::pose::PoseOP toRelax_poseOP, input_poseOP;
		Size offset;
		get_poseToRelax(toRelax_poseOP,input_poseOP,input,alns,offset);
		string const pose_tag( file_basename( core::pose::tag_from_pose( *input_poseOP ) ) );
		string pdbid = pose_tag.substr(0,pose_tag.length()-4);
		std::set< Size> caAtomsToConstrain= get_residuesToConstrain(*toRelax_poseOP);
		string coordCstFile  = pdbid+".coordCsts";
		//input_pose is used so that the coord-csts are in the correct position on the template
		output_caAtomsToConstraint(caAtomsToConstrain,offset,coordCstFile,input_poseOP);
		if( option[ minimalCstRelax::relax_pdb]()){
			core::pose::PoseOP finalRelax_poseOP =  new core::pose::Pose(*toRelax_poseOP);
			protocols::simple_moves::MissingDensityToJumpMoverOP fixMissingDensityMover (new protocols::simple_moves::MissingDensityToJumpMover());
			core::scoring::ScoreFunctionOP scorefxn_w_csts_ = get_score_function();
			scorefxn_w_csts_->set_weight( coordinate_constraint, COORDINATE_CST_WT );
			fixMissingDensityMover->apply(*finalRelax_poseOP);
			add_virtual_residue_to_cterm(*finalRelax_poseOP);
			ConstraintSetOP cst_set = convert_caAtomsToConstrain_to_coordCsts(caAtomsToConstrain,*toRelax_poseOP);
			finalRelax_poseOP->constraint_set(cst_set);
			FastRelax relaxer3( scorefxn_w_csts_,5,"NO CST RAMPING" );
			relaxer3.apply(*finalRelax_poseOP);
			delete_virtual_residues(*finalRelax_poseOP);
			ss->fill_struct(*finalRelax_poseOP);
			ss->decoy_tag(pose_tag);
			Real gdtmm = core::scoring::CA_gdtmm(*finalRelax_poseOP,*toRelax_poseOP);
			ss->add_energy("gdtmm",gdtmm);
			string outputName = pose_tag.substr(0,pose_tag.length()-4);
			outputName = outputName+"_0001.pdb";
			finalRelax_poseOP->dump_pdb(outputName);
			sfd.write_silent_struct(*ss,option[ out::file::silent ]());
		}
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

