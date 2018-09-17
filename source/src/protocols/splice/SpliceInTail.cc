// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/splice/SpliceInTail.cc
/// @brief
/// @author Gideon Lapidoth (glapidoth@gmail.com)

// Unit headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/splice/SpliceSegment.hh>
#include <protocols/splice/TailSegmentMover.hh>
#include <protocols/splice/SpliceInTail.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/splice/SpliceInTailCreator.hh>
#include <utility/string_util.hh>
#include <utility/exit.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/AA.hh>
#include <protocols/protein_interface_design/filters/TorsionFilter.hh>
#include <protocols/protein_interface_design/util.hh>
#include <boost/algorithm/string/predicate.hpp>//for comparing string case insensitive
#include <protocols/task_operations/RestrictChainToRepackingOperation.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/AtomType.fwd.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <protocols/enzdes/AddorRemoveCsts.hh>
#include <protocols/splice/util.hh>
#include <protocols/moves/Mover.hh>
#include <core/id/SequenceMapping.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/datacache/DataMapObj.hh>
#include <protocols/minimization_packing/RotamerTrialsMinMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <protocols/protein_interface_design/movers/AddChainBreak.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <string>
#include <boost/algorithm/string.hpp>
//Auto Headers
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueProperty.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/task_operations/DesignAroundOperation.hh>
#include <protocols/task_operations/ProteinInterfaceDesignOperation.hh>
#include <protocols/task_operations/ThreadSequenceOperation.hh>
#include <protocols/task_operations/SeqprofConsensusOperation.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <numeric/xyzVector.hh>
#include <protocols/loops/FoldTreeFromLoopsWrapper.hh>
#include <protocols/protein_interface_design/movers/LoopLengthChange.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <numeric/constants.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/util.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/scoring/Energies.hh>
#include <numeric/xyz.functions.hh>
#include <protocols/simple_moves/CutChainMover.hh>
//////////////////////////////////////////////////
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
///////////////////////////////////////////////////
#include <fstream>
#include <ctime>
#include <protocols/splice/RBInMover.hh>
#include <protocols/splice/RBOutMover.hh>
#include <protocols/toolbox/superimpose.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/extra_pose_info_util.hh>

namespace protocols {
namespace splice {

using namespace core::conformation;
static  basic::Tracer TR_constraints("protocols.splice.Splice_constraints");
static  basic::Tracer TR("protocols.splice.SpliceInTail");
static  basic::Tracer TR_pssm("protocols.splice.Splice_pssm");
std::string SpliceInTailCreator::keyname() const {
	return SpliceInTailCreator::mover_name();
}

protocols::moves::MoverOP SpliceInTailCreator::create_mover() const {
	return protocols::moves::MoverOP(new SpliceInTail);
}

std::string SpliceInTailCreator::mover_name() {
	return "SpliceInTail";
}

SpliceInTail::SpliceInTail() : Mover(SpliceInTailCreator::mover_name())
{
	splicemanager.profile_weight_away_from_interface(1.0);
	splicemanager.design_shell(6.0);
	splicemanager.repack_shell(8.0);
	basic::options::option[basic::options::OptionKeys::out::file::pdb_comments].value(
		true);
	splicemanager.rb_sensitive(false);
	tolerance_ = 0.23;
	allowed_cuts_ = 1;
	tail_dbase_subset_.clear();
	end_tail_dbase_subset_ = DataccacheBoolDataOP( new basic::datacache::DataMapObj<bool> );

}
SpliceInTail::~SpliceInTail() = default;


void SpliceInTail::apply(core::pose::Pose & pose) {
	using namespace protocols::rosetta_scripts;
	using core::chemical::DISULFIDE;
	set_last_move_status(protocols::moves::MS_SUCCESS);
	splicemanager.scorefxn(scorefxn());
	//Get entry from dbase file

	TR<<database_pdb_entry()<<std::endl;
	core::Size const dbase_entry(find_dbase_entry(pose));
	if ( dbase_entry == 0 ) { // failed to read entry
		TR << "Did not find the request source conformation;Should we fail loudly if this happens??"<< std::endl;
		return;
	}
	splicemanager.mm(core::kinematics::MoveMapOP(new core::kinematics::MoveMap));

	runtime_assert(dbase_entry <= tail_torsion_database_.size());
	splicemanager.dofs() = tail_torsion_database_[dbase_entry];
	TR<<"Size of dofs:"<<splicemanager.dofs().size()<<std::endl;

	if ( splicemanager.dofs().cut_site() ) {
		splicemanager.tail_seg("c");
	} else {
		splicemanager.tail_seg("n");
	}
	splicemanager.source_pdb_name(splicemanager.dofs().source_pdb());
	protocols::splice::load_pdb_segments_from_pose_comments(pose,splicemanager.pdb_segments());

	splicemanager.modify_pdb_segments_with_current_segment(splicemanager.source_pdb_name());
	/// change the pose segment length to that of the source segment length
	core::pose::add_comment(pose, "segment_" + splicemanager.segment_type(), splicemanager.source_pdb_name());//change correct association between current loop and pdb file

	assign_from_res_to_res(pose);

	TR<<"pose from res:"<<splicemanager.pose_from_res()<<std::endl;
	TR<<"pose to res:"<<splicemanager.pose_to_res()<<std::endl;
	splicemanager.residue_diff(splicemanager.dofs().size() - (splicemanager.pose_to_res() - splicemanager.pose_from_res())-1);


	if ( splicemanager.dofs().cut_site() ) { //if c-ter
		splicemanager.cut_site(1);

	} else { //if n-ter
		splicemanager.cut_site(pose.total_residue()+splicemanager.residue_diff()-1);
	}

	if ( !(splicemanager.pose_from_res() && splicemanager.pose_to_res() && splicemanager.cut_site()) ) {
		TR << "WARNING something is messed up and I will fail because of:" << std::endl;
		TR << "from_res() is " << splicemanager.pose_from_res() << " to_res is " <<splicemanager.pose_to_res() << " cut_site() " << splicemanager.cut_site() << std::endl;
		TR << "pdb dumped as mess.pdb" << std::endl;
		pose.dump_pdb("mess.pdb");
		runtime_assert(splicemanager.pose_from_res() && splicemanager.pose_to_res() && splicemanager.cut_site());
	}
	TR<<"Source name is:"<<splicemanager.source_pdb_name()<<std::endl;

	TR<<"Residue diff: "<<splicemanager.residue_diff()<<std::endl;
	protocols::protein_interface_design::movers::LoopLengthChange llc;
	set_loop_length_change(llc);
	TR << "Foldtree before loop length change: " << pose.fold_tree()<< std::endl;
	TR << "cut_site:" << splicemanager.cut_site() << std::endl;
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name() + "_before_2ndllc_test.pdb");
	}
	llc.apply(pose);
	splicemanager.pose_to_res(splicemanager.pose_to_res() + splicemanager.residue_diff());
	build_ideal_segment(pose);
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name() + "_after_2ndllc_test.pdb");
	}
	TR << "Before setting the new fold_tree after loop length change:"<< std::endl;
	TR << "pose_from_res():" << splicemanager.pose_from_res() << std::endl;
	TR << "to_res():" << splicemanager.pose_to_res() << std::endl;
	TR << "residue_diff:" << splicemanager.residue_diff() << std::endl;
	splicemanager.anchor_res(set_anchor_res());
	set_fold_tree_nodes(pose);
	splicemanager.fold_tree(pose, fold_tree_nodes_);
	TR << "Foldtree after loop length change: " << pose.fold_tree() << std::endl;
	/// set torsions
	core::Size const total_residue_new(splicemanager.dofs().size()-1); //how long is the introduced segment, I subtract 1 because the size of dofs includes +1 for the name of the segment (see db file)
	//Set BB psi/phi/omega angles on pose segment
	TR << "Changing dofs" << std::endl;
	splicemanager.set_BB_dofs(pose);
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name() + "_after_change_BB_dofs.pdb");
	}
	TR << "Ended change dofs. Now testing chain breaks." << std::endl;
	core::Real tolerance = tolerance_; // 0.23 angstrom deviation from mean peptide bond length
	bool fail_retry_if_found = true;
	bool crash_if_found = false;
	splicemanager.chainbreak_check( pose , tolerance, fail_retry_if_found , crash_if_found );
	splicemanager.add_sequence_constraints(pose);
	////////////////////////////////////////////////////////////////////////////////////////////
	///Apply user defined design shell and repack shell around spliced segment
	core::pack::task::TaskFactoryOP tf;
	if ( splicemanager.task_factory() == nullptr ) {
		tf = core::pack::task::TaskFactoryOP(new core::pack::task::TaskFactory());
	} else {
		tf = core::pack::task::TaskFactoryOP(new core::pack::task::TaskFactory(*splicemanager.task_factory()));
	}
	using namespace protocols::task_operations;
	DesignAroundOperationOP dao(new DesignAroundOperation);
	dao->design_shell((splicemanager.task_factory() == NULL ? 0.0 : splicemanager.design_shell())); // threaded sequence operation needs to design, and will restrict design to the loop, unless design_task_factory is defined, in which case a larger shell can be defined
	dao->repack_shell(splicemanager.repack_shell());
	for ( core::Size i = splicemanager.pose_from_res(); i <= splicemanager.pose_to_res(); ++i ) {
		if ( !pose.residue(i).has_variant_type(DISULFIDE) ) {
			dao->include_residue(i);
		}
	}
	tf->push_back(dao);
	//splicemanager.task_factory(tf);

	//////////////////Add chainbreak attribute to the cut site residues.////////////////////////
	//Add chainbreak attribute to the cut site resiudes.
	protocols::protein_interface_design::movers::AddChainBreak acb;
	acb.resnum(utility::to_string(splicemanager.cut_site()));
	acb.find_automatically(false);
	acb.change_foldtree(false);
	TR << "Adding ccd chainbreak at: " << splicemanager.cut_site() << std::endl;
	acb.apply(pose);
	////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////Perform Repacking and Rotamer Trial min mover here on the segment to improve energy ////////////////
	core::pack::task::PackerTaskOP ptask = tf->create_task_and_apply_taskoperations(pose);
	protocols::minimization_packing::PackRotamersMover prm(scorefxn(), ptask);
	prm.apply(pose);
	rtmin(pose,tf);
	////////////////////////////////////////////////////////////////////////////////////////////




	splicemanager.mm()->set_chi(false);
	splicemanager.mm()->set_bb(false);
	// splicemanager.mm()->set_jump(false);
	core::kinematics::FoldTree ft;
	ft = pose.fold_tree();
	/* if (conf.num_chains() > 1) {//if ligand is present we need to add edge between receptor and ligand
	core::Size const jump_num = ft.num_jump();
	TR << "Number of jumps in the fold tree is:" << jump_num << std::endl;
	splicemanager.mm()->set_jump(jump_num, true); /// 1Feb13 for cases in which the we're splicing in the presence of a ligand

	}*/
	for ( core::Size i = splicemanager.pose_from_res(); i <= splicemanager.pose_to_res(); ++i ) {
		splicemanager.mm()->set_chi(i, true);
		splicemanager.mm()->set_bb(i, true);
	}
	//Apply coordinate and dihedral constraints to spliced segment
	core::Size anchor_res(set_anchor_res());

	splicemanager.residue_diff(0);//after changing the segment length, the residue diff is null
	std::vector<std::string> BB_atoms = { "CA", "O", "N", "CB" };
	for ( auto BB_atom : BB_atoms ) { // access by const reference
		splicemanager.add_coordinate_constraints(pose, pose, splicemanager.pose_from_res(),splicemanager.pose_from_res() + total_residue_new - 1,anchor_res, BB_atom);
	}
	//add coordinate constraints to the entire side chain. I can only apply this in cases where the
	splicemanager.add_dihedral_constraints(pose, pose,splicemanager.pose_from_res()+1,splicemanager.pose_from_res() + total_residue_new - 1); //The plus 1 for from res is bcause dihedral constraints are added between i and i-1
	////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1<core::Size> designable_residues = residue_packer_states(pose, tf, true, false);
	TR << "Residues Allowed to Design: " << std::endl;
	for ( utility::vector1<core::Size>::const_iterator i(designable_residues.begin()); i != designable_residues.end(); ++i ) {
		TR << *i << "+";
	}
	TR << std::endl;

	utility::vector1<core::Size> packable_residues = residue_packer_states(pose,tf, true, true);
	TR << "Residues Allowed to Repack: " << std::endl;
	for ( utility::vector1<core::Size>::const_iterator i(
			packable_residues.begin()); i != packable_residues.end(); ++i ) {
		TR << *i << "+";
	}
	TR << std::endl;
	//TR << splicemanager <<std::endl;
	//splicemanager.mm()->show();
	//mm=splicemanager.mm();
	minimize_segment(pose);
	splicemanager.add_sequence_constraints(pose);
	acb.apply(pose);
	TR << "fold tree before chain break val score: " << pose.fold_tree()<< std::endl;
	acb.remove(true);//remove cutpoint variant, this will cause trouble in itterative splicein moves. GDL sep17
	acb.apply(pose);
} //apply

std::string SpliceInTail::get_name() const {
	return SpliceInTailCreator::mover_name();
}

void SpliceInTail::parse_my_tag(TagCOP const tag, basic::datacache::DataMap &data,protocols::filters::Filters_map const & /*filters*/,protocols::moves::Movers_map const &/* movers*/,core::pose::Pose const & /*pose*/) {
	typedef utility::vector1<std::string> StringVec;
	utility::vector1<TagCOP> const sub_tags(tag->getTags());
	splicemanager.parse_tags(tag,data);
	tolerance_ = tag->getOption < core::Real > ("tolerance", 0.23); //for debugging purposes
	splicemanager.allowed_cuts(tag->getOption < core::Size > ("allowed_cuts", 1));
	dbase_iterate(tag->getOption<bool>("dbase_iterate", false));
	min_seg(tag->getOption<bool>("min_seg", true));
	//TR<<"READING TORSION DB "<<splicemanager.dbase_file_name()<<std::endl;

	std::string delta;
	if ( tag->hasOption("delta_lengths") ) {
		delta = tag->getOption<std::string>("delta_lengths");
		//check in put is valid, either : "1,2,3,..." or "1:10", Gideon 19may15
		StringVec const lengths_keys(utility::string_split(delta, ','));
		for ( std::string const & delta: lengths_keys ) {
			if ( delta == "" ) continue;
			int const delta_i( 1 * atoi( delta.c_str() ) );
			delta_lengths_.push_back( delta_i );
		}//for
	} else { //fi  (tag->hasOption("delta_lengths"))
		delta_lengths_.push_back(0);
	}
	std::sort(delta_lengths_.begin(), delta_lengths_.end());
	std::unique(delta_lengths_.begin(), delta_lengths_.end());


	//Read map between segments and PSSMs
	splicemanager.parse_segments(sub_tags, tag, data);
	scorefxn(protocols::rosetta_scripts::parse_score_function(tag, data));



	database_entry(tag->getOption<core::Size>("database_entry", 0));
	database_pdb_entry(tag->getOption<std::string>("database_pdb_entry", ""));
	runtime_assert(!(tag->hasOption("database_entry") && tag->hasOption("database_pdb_entry")));
	read_torsion_database();
	TR << "torsion_database: " << splicemanager.dbase_file_name() << " ";
	if ( database_entry() == 0 ) {
		if ( database_pdb_entry_ == "" ) {
			TR << " database entry will be randomly picked at run time. ";
		} else {
			TR << " picking database entry " << database_pdb_entry() << std::endl;
		}
	} else {
		TR << " database_entry: " << database_entry() << " ";
		runtime_assert(database_entry() <= tail_torsion_database_.size());
	}
}

protocols::moves::MoverOP SpliceInTail::clone() const {
	return (protocols::moves::MoverOP(new SpliceInTail(*this)));
}

core::Size SpliceInTail::find_dbase_entry(core::pose::Pose const & pose) {
	core::Size dbase_entry(database_entry());
	TR<<"The size of tail_torsion_database_ is:"<<tail_torsion_database_.size()<<std::endl;
	if ( first_pass_ ) { /// setup the dbase subset where loop lengths fit the selection criteria
		for ( core::Size i = 1; i <= tail_torsion_database_.size(); ++i ) {// find entries that fit the length criteria
			using namespace protocols::rosetta_scripts;
			ResidueBBDofs const & dofs(tail_torsion_database_[i]);

			splicemanager.dofs(tail_torsion_database_[i]);
			core::pose::Pose tmp_pose = pose;
			assign_from_res_to_res(tmp_pose);
			//TR<<"splicemanager.pose_from_res()" <<splicemanager.pose_from_res()<<std::endl;
			//TR<<"splicemanager.pose_to_res()" <<splicemanager.pose_to_res()<<std::endl;
			core::Size const nearest_to_entry_start_on_pose(splicemanager.pose_from_res());
			core::Size const nearest_to_entry_stop_on_pose(splicemanager.pose_to_res());
			core::Size const pose_residues = nearest_to_entry_stop_on_pose - nearest_to_entry_start_on_pose + 1;

			//   //TR<<"Dofs start loop: "<<dofs.start_loop()<< " Dofs stop loop:" << dofs.stop_loop()<<"," <<dofs.source_pdb() <<std::endl;
			//   core::Size nearest_to_entry_start_on_pose =0;
			//   core::Size nearest_to_entry_stop_on_pose = 0;
			//
			//   if (dofs.cut_site()==0){ //n-ter tail
			//    nearest_to_entry_start_on_pose = 1;
			//    nearest_to_entry_stop_on_pose= find_nearest_res(pose, *splicemanager.template_pose(),dofs.start_loop(), 0/*chain*/);
			//   }
			//   else if (dofs.cut_site()==1){//c-ter tail
			//    nearest_to_entry_start_on_pose = find_nearest_res(pose, *splicemanager.template_pose(),dofs.start_loop(), 0/*chain*/);
			//    nearest_to_entry_stop_on_pose = pose.conformation().chain_end(1);
			//    //TR<<" nearest_to_entry_start_on_pose:"<< nearest_to_entry_start_on_pose<<std::endl;
			//    //TR<<" dofs from res:"<< dofs.start_loop()<<std::endl;
			//    //TR<<" nearest_to_entry_stop_on_pose:"<< nearest_to_entry_stop_on_pose<<std::endl;
			//   }
			//   //TR<<" nearest_to_entry_start_on_pose:"<< nearest_to_entry_start_on_pose<<std::endl;
			//
			//   //Gideon, 5May15: So nearest stop entry on pose is a bit problemtic.  In tail segments find_nearest won't always work.
			//   // So instead I will do this case by case:
			//
			//   core::Size const pose_residues = nearest_to_entry_stop_on_pose - nearest_to_entry_start_on_pose + 1;
			//   //TR<<"Size of curent loop:"<<pose_residues<<std::endl;

			int const delta(dofs.size() - pose_residues);
			// TR << "The size of pose is " << pose_residues << std::endl;
			// TR << "The size of possible insert is " << dofs.size() << std::endl;
			// TR << "The delta is " << delta << std::endl;

			bool const fit = std::find(delta_lengths_.begin(), delta_lengths_.end(), delta) != delta_lengths_.end();
			if ( fit || database_pdb_entry_ != "" || dbase_entry != 0 ) {
				tail_dbase_subset_.push_back(i);
			}
		}   //for (core::Size i = 1; i <= tail_torsion_database_.size(); ++i)
		if ( tail_dbase_subset_.empty() ) {
			TR << "Loop of appropriate length not found in database. Returning"
				<< std::endl;
			return 0;
		}
		//TR<<tail_dbase_subset_[0]<<std::endl;
		//get the previous and next segments given the current segment, in case there are other segments
		numeric::random::random_permutation(tail_dbase_subset_.begin(),
			tail_dbase_subset_.end(), numeric::random::rg());
		current_dbase_entry_ = tail_dbase_subset_.begin();
		load_from_checkpoint();
		first_pass_ = false;
	} // fi first_pass
	if ( dbase_iterate() ) {
		load_from_checkpoint();
		if ( current_dbase_entry_ == dbase_end() ) {
			TR<< "Request to read past end of dbase. Splice returns without doing anything."<< std::endl;
			return 0;
		}
		dbase_entry = *current_dbase_entry_;
		if ( !first_pass_ ) {
			current_dbase_entry_++;
		}
		if ( current_dbase_entry_ == dbase_end() ) {
			TR << "Reached last dbase entry" << std::endl;
			end_tail_dbase_subset_->obj = true;
		}
	} else if ( dbase_entry == 0 ) { // fi dbase_iterate
		if ( database_pdb_entry_ == "" ) { //randomize dbase entry
			TR << "The dbase_subset size is " << tail_dbase_subset_.size()<< std::endl;
			core::Size entry_no_to_choose = (core::Size) (numeric::random::rg().uniform()* tail_dbase_subset_.size() + 1);
			if ( tail_dbase_subset_.size() == 0 ) {
				//if ( debug_)
				pose.dump_pdb(splicemanager.mover_name() + "db_size_0.pdb");
				utility_exit_with_message("The database size is "+ utility::to_string(tail_dbase_subset_.size())+ ". Dumped " + splicemanager.mover_name()+ "db_size_0.pdb to inspect");
			}

			TR << "trying to pick entry " << entry_no_to_choose << std::endl;
			dbase_entry = tail_dbase_subset_[entry_no_to_choose];
		} else { // look for the pdb_entry name
			//dbase_entry = ( core::Size )( RG.uniform() * tail_dbase_subset_.size() + 1 );
			for ( core::Size count = 1; count <= tail_dbase_subset_.size(); ++count ) {
				//TR<<tail_torsion_database_[tail_dbase_subset_[count]].source_pdb()<<std::endl;
				if ( tail_torsion_database_[tail_dbase_subset_[count]].source_pdb()== database_pdb_entry_ ) {
					TR << "Found entry for " << database_pdb_entry_
						<< " at number " << tail_dbase_subset_[count]
						<< std::endl;
					dbase_entry = tail_dbase_subset_[count];
					break;
				}
			}
			//This assertion does not make sense. Need to fix this. Gideon 4May15
			runtime_assert(dbase_entry <= tail_dbase_subset_.size());
		} //else
	} //fi dbase_entry==0
	return dbase_entry;
}



std::string
SpliceInTail::mover_name()  {
	return "SpliceInTail";
}

std::string SpliceInTail_complex_type_name_for_subsubtag( std::string const & foo ) {
	return "subsubtag_SpliceInTail_" + foo + "_type";
}

std::string SpliceInTail_complex_type_name_for_subtag( std::string const & foo ) {
	return "subtag_SpliceInTail_" + foo + "_type";
}


void SpliceInTail::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;

	XMLSchemaRestriction n_or_c;
	n_or_c.name( "n_or_c" );
	n_or_c.base_type( xs_string );
	n_or_c.add_restriction( xsr_enumeration, "n" );
	n_or_c.add_restriction( xsr_enumeration, "c" );
	xsd.add_top_level_element( n_or_c );

	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "debug", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "min_seg", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "CG_const", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rb_sensitive", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "chain_num", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute( "segment", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "superimposed", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_n", xsct_non_negative_integer, "XRW TO DO", "4" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_c", xsct_non_negative_integer, "XRW TO DO", "13" )
		+ XMLSchemaAttribute( "tail_segment", "n_or_c", "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "skip_alignment", xsct_rosetta_bool, "XRW TO DO", "false" );

	// The "Segments" subtag
	AttributeList segments_subtag_attlist;
	segments_subtag_attlist + XMLSchemaAttribute::attribute_w_default( "profile_weight_away_from_interface", xsct_real, "XRW TO DO", "1.0" )
		+ XMLSchemaAttribute( "current_segment", xs_string, "XRW TO DO" );

	// The "segment" sub-subtag"
	AttributeList subtag_segments_subtag_attlist;
	subtag_segments_subtag_attlist + XMLSchemaAttribute( "pdb_profile_match", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "profiles", xs_string, "XRW TO DO" );

	XMLSchemaComplexTypeGenerator segment_subsubtag_gen;
	segment_subsubtag_gen.complex_type_naming_func( & SpliceInTail_complex_type_name_for_subsubtag )
		.element_name( "segment" )
		.description( "individual segment tag" )
		.add_attributes( subtag_segments_subtag_attlist )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList subsubelements;
	subsubelements.add_already_defined_subelement( "segment", SpliceInTail_complex_type_name_for_subsubtag/*, 0*/ );

	XMLSchemaComplexTypeGenerator segments_subtag_gen;
	segments_subtag_gen.complex_type_naming_func( & SpliceInTail_complex_type_name_for_subtag )
		.element_name( "Segments" )
		.description( "Wrapper for multiple segments tags" )
		.add_attributes( segments_subtag_attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subsubelements )
		.write_complex_type_to_schema( xsd );


	XMLSchemaSimpleSubelementList subelements;
	subelements.add_already_defined_subelement( "Segments", SpliceInTail_complex_type_name_for_subtag/*, 0*/ );

	attlist + XMLSchemaAttribute( "use_sequence_profile", xsct_rosetta_bool, "XRW TO DO" );
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "add_sequence_constraints_only", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "template_file", xs_string, "XRW TO DO" );
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist

		+ XMLSchemaAttribute( "residue_numbers_setter", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "torsion_database", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "database_entry", xsct_non_negative_integer, "XRW TO DO" )
		+ XMLSchemaAttribute( "database_pdb_entry", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "design_shell", xsct_real, "XRW TO DO", "6.0" )
		+ XMLSchemaAttribute::attribute_w_default( "repack_shell", xsct_real, "XRW TO DO", "8.0" )
		+ XMLSchemaAttribute::attribute_w_default( "thread_ala", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute( "delta_lengths", xsct_int_cslist, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "thread_original_sequence", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rtmin", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "dbase_iterate", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "checkpointing_file", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "use_sequence_profiles", xsct_rosetta_bool, "XRW TO DO", "true" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, subelements );
}

void SpliceInTailCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SpliceInTail::provide_xml_schema( xsd );
}
void SpliceInTail::set_loop_length_change( protocols::protein_interface_design::movers::LoopLengthChange & llc){
	TR<<"loop start:"<<splicemanager.pose_from_res()<<std::endl;
	TR<<"loop end:"<<splicemanager.pose_to_res()<<std::endl;
	llc.loop_start(splicemanager.pose_from_res());
	llc.loop_end(splicemanager.pose_to_res());
	llc.tail(1);
	llc.delta(splicemanager.residue_diff());
	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		llc.direction(1);
	} else if ( boost::iequals(splicemanager.tail_seg(), "c") ) {
		llc.direction(0);
	}

}
void SpliceInTail::set_fold_tree_nodes(core::pose::Pose const & pose){
	splicemanager.mm()->set_jump(false);

	core::conformation::Conformation const & conf(pose.conformation());
	TR<<"Num of chains: "<<conf.num_chains()<<std::endl;
	fold_tree_nodes_.clear();
	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		fold_tree_nodes_.push_back({{(int) splicemanager.pose_to_res()+1,(int) splicemanager.pose_from_res(),-1}});
		fold_tree_nodes_.push_back({{(int) splicemanager.pose_to_res()+1,(int)conf.chain_end(1),-1}});
	} else if ( boost::iequals(splicemanager.tail_seg(), "c") ) {
		fold_tree_nodes_.push_back({{(int) splicemanager.pose_from_res()-1,(int) splicemanager.pose_to_res(),-1}});
		fold_tree_nodes_.push_back({{(int) splicemanager.pose_from_res()-1,1,-1}});
	}
	if ( conf.num_chains()>1 ) { //if ligand is present we need to add edge between receptor and ligand
		core::Size CoM = (core::Size ) core::pose::residue_center_of_mass( pose, conf.chain_begin(2), conf.chain_end(2) );
		fold_tree_nodes_.push_back({{(int) conf.chain_end(1),(int) CoM,1}});
		splicemanager.mm()->set_jump(1, true);

	}
}

core::Size SpliceInTail::set_anchor_res(){
	if ( boost::iequals(splicemanager.tail_seg(), "n") ) {
		return splicemanager.pose_to_res()+1;
	} else if ( boost::iequals(splicemanager.tail_seg(), "c") ) {
		return splicemanager.pose_from_res()-1;
	}
	return 0;
}
void SpliceInTail::assign_from_res_to_res(core::pose::Pose const & pose){
	using namespace protocols::rosetta_scripts;
	TR<<"I am here"<<std::endl;
	if ( splicemanager.dofs().cut_site() ) { //if true this a C-ter tail
		splicemanager.pose_from_res(find_nearest_res(pose,*splicemanager.template_pose(),splicemanager.dofs().start_loop(),0));
	} else { //if n-ter
		splicemanager.pose_from_res(1);
	}

	if ( splicemanager.dofs().cut_site() ) { //if true this a C-ter tail
		TR<<"I am here"<<std::endl;
		core::conformation::Conformation const & conf(pose.conformation());
		splicemanager.pose_to_res(conf.chain_end(1));
		TR<<"splicemanager.pose_to_res()"<<splicemanager.pose_to_res()<<std::endl;
	} else { //if n-ter
		splicemanager.pose_to_res(find_nearest_res(pose,*splicemanager.template_pose(),splicemanager.dofs().start_loop(),0));
	}
}
void SpliceInTail::build_ideal_segment(core::pose::Pose & pose){
	using namespace core::chemical;
	using namespace core::conformation;
	ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose() );
	TR<<"Building ideal segment"<<std::endl;
	TR<<"build_ideal_segment:"<<splicemanager.pose_from_res()<<std::endl;
	TR<<"build_ideal_segment:"<<splicemanager.pose_to_res()<<std::endl;


	ResidueCOP new_res = ResidueFactory::create_residue( residue_set->name_map( name_from_aa( aa_from_oneletter_code( 'A' ) ) ) );
	core::Size new_segment_length(splicemanager.pose_to_res()-splicemanager.pose_from_res());

	if ( (boost::iequals(splicemanager.tail_seg(), "n")) ) {
		core::Size new_to_res = splicemanager.pose_to_res()-(splicemanager.pose_to_res()-splicemanager.pose_from_res());
		pose.delete_residue_range_slow( splicemanager.pose_from_res(), splicemanager.pose_to_res()-1);
		//pose.dump_pdb("after_delete_segment.pdb");
		for ( core::Size res=1; res<=new_segment_length; res++ ) {
			pose.conformation().safely_prepend_polymer_residue_before_seqpos(*new_res,new_to_res,true);
			//TR<<"Appending residue before residue: "<<res+new_to_res-1<<std::endl;
			if ( splicemanager.debug() ) {
				pose.dump_pdb(splicemanager.mover_name() +"_build_ideal_segment.pdb");
			}
		}
	} else if ( (boost::iequals(splicemanager.tail_seg(), "c")) ) {
		pose.delete_residue_range_slow( splicemanager.pose_from_res()+1, splicemanager.pose_to_res());
		if ( splicemanager.debug() ) {
			pose.dump_pdb(splicemanager.mover_name() +"after_delete_segment.pdb");
		}
		for ( core::Size res=splicemanager.pose_from_res()+1; res<=splicemanager.pose_to_res(); res++ ) {
			pose.append_polymer_residue_after_seqpos(*new_res,res-1,true);
			if ( splicemanager.debug() ) {
				pose.dump_pdb(splicemanager.mover_name() +"_build_ideal_segment.pdb");
			}
		}
	}
}

} //splice
} //protocols
