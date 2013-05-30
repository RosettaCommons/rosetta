// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.AlignmentCluster

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

#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/electron_density/util.hh>
#include <core/scoring/electron_density/FastDensEnergy.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

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
#include <iostream>
#include <fstream>

basic::Tracer tr( "unalignedEvaluate" );
using std::string;
using std::map;
using std::vector;
using core::pose::Pose;
using protocols::loops::Loops;
using namespace core::scoring::constraints;


Real region_density_score(Pose & pose,Pose const templatePose, SequenceAlignment aln, Size start_res,Size end_res){
	using namespace core::scoring;
	using namespace core::scoring::methods;
	using namespace core::sequence;
	using namespace core::id;
	EnergyMethodOptions eopts;
	superimpose_pose_using_aln(pose,templatePose,aln,2);
	core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
	EnergyMap weights( pose.energies().weights() );
	for ( int ii = 1; ii <= n_score_types; ++ii )
		weights[ ScoreType(ii) ] = 0.0;
	core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn );
	scorefxn->score(pose);
	protocols::moves::SequenceMoverOP seqmov = new protocols::moves::SequenceMover;
	seqmov->add_mover( new protocols::electron_density::SetupForDensityScoringMover );
	seqmov->apply(pose);
	Real density_score = 0;
	for(Size ii = start_res; ii <= end_res; ++ii){
		Real score = pose.energies().residue_total_energies(ii)[core::scoring::elec_dens_fast];
		density_score += score;
	}
	//Real density_score = scorefxn->get_sub_score( pose, mask );
		/*core::scoring::EnergyMap emap;
	emap.zero();
	Real total_score = 0;
	core::scoring::methods::EnergyMethodOptionsOP emopts(
			new core::scoring::methods::EnergyMethodOptions( scorefxn->energy_method_options() )
	);
	for(Size ii = start_res; ii <= end_res; ++ii){
	  Size jj = pose.total_residue();
		emopts->residue_pair_energy(pose.residue(ii), pose.residue(jj), pose, *scorefxn, emap);
		total_score += emap[ core::scoring::elec_dens_fast ];
	}
	*/
	return (density_score);
}


Real region_constraint_score(Pose & pose, Size start_res, Size end_res,ConstraintSetOP cstSet_){
	using namespace core::scoring::constraints;
	core::scoring::ScoreFunctionOP scorefxn( new core::scoring::ScoreFunction );
	core::scoring::EnergyMap emap;
	emap.zero();
	Real total_score = 0;
	for(Size ii = start_res; ii <= end_res; ++ii){
		for(Size jj=1; jj <= pose.total_residue(); ++jj){
			cstSet_->residue_pair_energy(pose.residue(ii), pose.residue(jj), pose, *scorefxn, emap);
			total_score += emap[ core::scoring::atom_pair_constraint ];
		}
	}
	return (total_score);
}

Real region_score(Pose & pose, Size start_res, Size end_res){
	using namespace core::scoring::constraints;
	using namespace core::scoring;
	core::scoring::ScoreFunctionOP scorefxn( getScoreFunction() );
	Real all_res_score = scorefxn->score(pose);
	utility::vector1< bool > mask( pose.total_residue(), true );
	for(Size ii= start_res; ii<=end_res; ++ii)
		mask[ii] = false;
	Real excluded_region_score = scorefxn->get_sub_score( pose, mask );
	return(all_res_score-excluded_region_score);
}

Real region_rmsd_structures_aligned(Pose & mod_pose, Pose const & ref_pose, SequenceAlignment aln, Size  start_res, Size end_res){
	using namespace core::sequence;
	using namespace core::id;
	superimpose_pose_using_aln(mod_pose,ref_pose,aln,2);
	SequenceMapping map = aln.sequence_mapping(2,1);
	Real sum = 0;
	for (Size ii =start_res; ii<= end_res; ++ii){
		Size ref_res = map[ii];
		sum += mod_pose.residue(ii).xyz("CA").distance(ref_pose.residue(ref_res).xyz("CA"));
	}
	return(sum/(end_res-start_res+1));
}

Real two_region_rmsd(Pose & mod_pose, Pose const & ref_pose, SequenceAlignment aln, Size start_res, Size end_res, Size region2_start_res, Size region2_end_res){
	using namespace core::sequence;
	using namespace core::scoring;
	using namespace core::id;
	vector1<Size> mod_pose_positions;
	vector1<Size> ref_pose_positions;
	SequenceMapping map = aln.sequence_mapping(2,1);
	for (Size ii=start_res; ii<= end_res; ++ii){
		mod_pose_positions.push_back(ii);
		ref_pose_positions.push_back(map[ii]);
	}
	if(region2_start_res != 0){
		for(Size ii=region2_start_res; ii<=region2_end_res; ++ii){
			mod_pose_positions.push_back(ii);
			ref_pose_positions.push_back(map[ii]);
		}
	}
	core::kinematics::FoldTree f_mod(mod_pose_positions.size());
	core::kinematics::FoldTree f_ref(mod_pose_positions.size());
	Pose sub_mod_pose;
	Pose sub_ref_pose;
	core::pose::create_subpose(mod_pose,mod_pose_positions,f_mod,sub_mod_pose);
	core::pose::create_subpose(ref_pose,ref_pose_positions,f_ref,sub_ref_pose);
	return(CA_rmsd(sub_mod_pose,sub_ref_pose));
}

Real region_rmsd(Pose & mod_pose, Pose const & ref_pose, SequenceAlignment aln, Size start_res, Size end_res){
	return(two_region_rmsd(mod_pose,ref_pose,aln,start_res,end_res,0,0));
}

bool present_in_native(Size start_res, Size end_res, SequenceAlignment native_aln){
	using namespace core::sequence;
	using namespace core::id;
	SequenceMapping map = native_aln.sequence_mapping(2,1);
	bool allResiduesFound = true;
	if(start_res <= 0)
		return(false);
	for(Size ii=start_res;ii<=end_res; ii++)
		if(map[ii] == 0)
			allResiduesFound = false;
	return(allResiduesFound);
}

void score_loop(Pose & pose, Pose const & native_pose, SequenceAlignment native_aln,SequenceAlignment template_aln, Pose native_relaxed_pose,  Size  start_res, Size end_res,ConstraintSetOP cstSet_, ConstraintSetOP lookBackCstSet_, Pose const templatePose, std::ostream& out){
	using namespace core::sequence;
	using namespace core::id;
	using namespace core::scoring;
	Size STEM_RESIDUES = 5;
	//note checking if start_stem exists is done in loop present in native
	Size start_stem = start_res-STEM_RESIDUES-1;
	Size end_stem = end_res+STEM_RESIDUES+1;
	SequenceMapping map = native_aln.sequence_mapping(2,1);
	Real rmsd_stem = 9999;
	bool tail = false;
	bool disorded_stem = false;
	if((present_in_native(start_stem,start_res-1,native_aln))&&(present_in_native(end_res+1,end_stem,native_aln))){
		rmsd_stem = two_region_rmsd(pose,native_pose,native_aln,start_stem,start_res-1,end_res+1,end_stem);
	}
	else if(present_in_native(start_stem,start_res-1,native_aln)){
		rmsd_stem = region_rmsd(pose,native_pose,native_aln,start_stem,start_res-1);
		if (start_stem <= 0)
			tail = true;
		else
			disorded_stem = true;
	}
 	else if(present_in_native(end_res+1,end_stem,native_aln)){
		rmsd_stem = region_rmsd(pose,native_pose,native_aln,end_res+1,end_stem);
		if(end_res >= map.size2())
			tail = true;
		else
			disorded_stem = true;
	}
	Real rmsd_loop_wStruct = region_rmsd_structures_aligned(pose,native_pose,native_aln,start_res,end_res);
	Real rmsd_loop = region_rmsd(pose,native_pose,native_aln,start_res,end_res);
	Size length = end_res-start_res+1;
	Real constraint_score = region_constraint_score(pose,start_res,end_res,cstSet_);
	Real score12 = region_score(pose,start_res,end_res);
	Real native_constraint_score = region_constraint_score(native_relaxed_pose,start_res,end_res,cstSet_);
	Real native_score12 = region_score(native_relaxed_pose,start_res,end_res);
	core::scoring::ScoreFunctionOP scorefxn( getScoreFunction() );
	Real all_res_score = scorefxn->score(pose);
	Real density_score = region_density_score(pose,templatePose, template_aln, start_res,end_res);
	Real lookBack_constraint_score = region_constraint_score(pose,start_res,end_res,lookBackCstSet_);
	out <<std::setprecision(2) << std::fixed;
	out <<"score \t" << start_res << "\t" << end_res << "\t" << length << "\t" << rmsd_loop << "\t" << rmsd_loop_wStruct << "\t"<< rmsd_stem << "\t" << constraint_score <<"\t" << score12 <<"\t" << native_constraint_score << "\t" << native_score12 << "\t" << tail << "\t" << disorded_stem << "\t" << all_res_score << "\t" << density_score <<"\t" << lookBack_constraint_score <<  std::endl;
}

void score_loops(Pose pose,	Loops const unalignedLoops, Pose native_pose,SequenceAlignment native_alignment, SequenceAlignment template_alignment, Pose native_relaxed_pose, Pose const templatePoseData,ConstraintSetOP cstSet_, ConstraintSetOP lookBackCstSet_, std::ostream& out){

	for(Loops::const_iterator ii = unalignedLoops.begin(); ii != unalignedLoops.end(); ++ii){
		Size start_res = ii->start();
		Size end_res = ii->stop();
		if(present_in_native(start_res,end_res,native_alignment)){
			score_loop(pose,native_pose,native_alignment, template_alignment, native_relaxed_pose, start_res,end_res,cstSet_,lookBackCstSet_, templatePoseData, out);
		}
	}
}


namespace unalignedEvaluate {
basic::options::FileOptionKey relaxed_native("unalignedEvaluate:relaxed_native");
basic::options::FileOptionKey core_pdb("unalignedEvaluate:core_pdb");

basic::options::FileOptionKey lookback_csts("unalignedEvaluate:lookback_csts");
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
	option.add (unalignedEvaluate::relaxed_native, "relaxed native");
	option.add (unalignedEvaluate::core_pdb,"core pdb probably want to use theseus");
	option.add (unalignedEvaluate::lookback_csts,"lookback_csts");
	devel::init(argc, argv);
	map< string, SequenceAlignment> alnDataMapped = input_alignmentsMapped(true);
	map< string, SequenceAlignment>::iterator alnDataMapped_itr;
	string query_sequence (
			read_fasta_file( option[ in::file::fasta ]()[1])[1]->sequence()	);
	map< string, const Loops > unalignedLoopsMapped = get_unalignedLoopsMapped(alnDataMapped,query_sequence.size());
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	Pose native_pose;
	pose_from_pdb(
			native_pose,
			*rsd_set,
			option[ in::file::native ]()
			);
	Pose native_relaxed_pose;
	pose_from_pdb(
			native_relaxed_pose,
			*rsd_set,
			option[unalignedEvaluate::relaxed_native]()
			);
	Pose core_pose;
	pose_from_pdb(
			core_pose,
			*rsd_set,
			option[unalignedEvaluate::core_pdb]()
			);
	map< string, Pose > templatePoseData = poses_from_cmd_line_noPDBtag(
			option[ in::file::template_pdb ]());
	using namespace core::scoring::constraints;
	std::string cstfile = core::scoring::constraints::get_cst_file_option();
	ConstraintSetOP cstSet_ = ConstraintIO::read_constraints( cstfile, new ConstraintSet, native_relaxed_pose );
	ConstraintSetOP lookBackCstSet_ = ConstraintIO::read_constraints(option[unalignedEvaluate::lookback_csts]() , new ConstraintSet, native_relaxed_pose );
	MetaPoseInputStream input = streams_from_cmd_line();
	SequenceOP native_sequence = new Sequence(native_pose);
	SequenceOP fasta_sequence = new Sequence(query_sequence,"fasta",1);
	SequenceAlignment native_alignment =  align_naive(native_sequence, fasta_sequence);
	std::string out_nametag = option[ out::file::o ];
	Pose const templatePose = templatePoseData.begin()->second;
	SequenceOP template_sequence = new Sequence(templatePose);
	SequenceAlignment template_aln =  align_naive(template_sequence, fasta_sequence);
	SequenceOP core_sequence = new Sequence(core_pose);
	SequenceAlignment core_aln =  align_naive( fasta_sequence, core_sequence);
	SequenceMapping map = core_aln.sequence_mapping(1,2);
	protocols::loops::LoopsOP unconverged_loops = loops_from_alignment(query_sequence.size(),core_aln,1);
	std::string out_file_name_str( out_nametag + ".scores" );
	utility::io::ozstream output(out_file_name_str);
	output << "labels \t start_res \t end_res \t length \t rmsd_loop \t rmsd_loop_struct_aln \t rmsd_stem \t constraint_score \t score12 \t native_constraint_score \t native_score12 \t disorded_stem \t all_res_score \t density_score \t lookback_cst_score" << std::endl;
	while(input.has_another_pose()){
		core::pose::PoseOP input_poseOP;
		input_poseOP = new core::pose::Pose();
		input.fill_pose(*input_poseOP,*rsd_set);
		score_loops(*input_poseOP,*unconverged_loops,native_pose,native_alignment,template_aln,native_relaxed_pose,templatePose,cstSet_, lookBackCstSet_, output);
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
}
