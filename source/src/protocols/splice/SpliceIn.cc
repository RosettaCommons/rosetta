// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/splice/SpliceIn.cc
/// @brief
/// @author Gideon Lapidoth (glapidoth@gmail.com)

// Unit headers
#include <core/pose/extra_pose_info_util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/splice/SpliceSegment.hh>
#include <protocols/splice/TailSegmentMover.hh>
#include <protocols/splice/SpliceIn.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/splice/SpliceInCreator.hh>
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


namespace protocols {
namespace splice {

using namespace core::conformation;
static  basic::Tracer TR_constraints("protocols.splice.Splice_constraints");
static  basic::Tracer TR("protocols.splice.SpliceIn");
static  basic::Tracer TR_pssm("protocols.splice.Splice_pssm");
std::string SpliceInCreator::keyname() const {
	return SpliceInCreator::mover_name();
}

protocols::moves::MoverOP SpliceInCreator::create_mover() const {
	return protocols::moves::MoverOP(new SpliceIn);
}

std::string SpliceInCreator::mover_name() {
	return "SpliceIn";
}

SpliceIn::SpliceIn() : Mover(SpliceInCreator::mover_name()), dbase_iterate_(false), checkpointing_file_(""),min_seg_(true),first_pass_(true)
{
	splicemanager.profile_weight_away_from_interface(1.0);
	splicemanager.design_shell(6.0);
	splicemanager.repack_shell(8.0);
	basic::options::option[basic::options::OptionKeys::out::file::pdb_comments].value(
		true);
	splicemanager.rb_sensitive(false);
	tolerance_ = 0.23;
	allowed_cuts_ = 1;
	dbase_subset_.clear();
	end_dbase_subset_ = DataccacheBoolDataOP( new basic::datacache::DataMapObj<bool> );

}
SpliceIn::~SpliceIn() = default;

void SpliceIn::assign_from_res_to_res(core::pose::Pose const & pose) {
	//pose.dump_pdb("set_from__res.pdb");
	using namespace protocols::rosetta_scripts;
	splicemanager.pose_from_res(find_nearest_res(pose,*splicemanager.template_pose(),splicemanager.dofs().start_loop(),0));
	splicemanager.pose_to_res (find_nearest_res(pose,*splicemanager.template_pose(),splicemanager.dofs().stop_loop(),0));
	splicemanager.cut_site (splicemanager.dofs().cut_site() - splicemanager.dofs().start_loop() + splicemanager.pose_from_res());
	splicemanager.residue_diff(splicemanager.dofs().size() - (splicemanager.pose_to_res() - splicemanager.pose_from_res()) -1);
	if ( !(splicemanager.pose_from_res() && splicemanager.pose_to_res() && splicemanager.cut_site()) ) {
		TR << "WARNING something is messed up and I will fail because of:" << std::endl;
		TR << "from_res() is " << splicemanager.pose_from_res() << " to_res is " <<splicemanager.pose_to_res() << " cut_site() " << splicemanager.cut_site() << std::endl;
		TR << "pdb dumped as mess.pdb" << std::endl;
		pose.dump_pdb("mess.pdb");
		runtime_assert(splicemanager.pose_from_res() && splicemanager.pose_to_res() && splicemanager.cut_site());
	}
}

void SpliceIn::apply(core::pose::Pose & pose) {
	using namespace protocols::rosetta_scripts;
	using core::chemical::DISULFIDE;
	set_last_move_status(protocols::moves::MS_SUCCESS);
	splicemanager.scorefxn(scorefxn());
	//Get entry from dbase file
	//assign_from_res_to_res(pose);
	core::Size const dbase_entry(find_dbase_entry(pose));
	if ( dbase_entry == 0 ) { // failed to read entry
		TR << "Did not find the request source conformation;Should we fail loudly if this happens??"<< std::endl;
		return;
	}
	splicemanager.mm(core::kinematics::MoveMapOP(new core::kinematics::MoveMap));
	runtime_assert(dbase_entry <= torsion_database_.size());
	splicemanager.dofs(torsion_database_[dbase_entry]);
	splicemanager.source_pdb_name(splicemanager.dofs().source_pdb());
	protocols::splice::load_pdb_segments_from_pose_comments(pose,splicemanager.pdb_segments());
	splicemanager.modify_pdb_segments_with_current_segment(splicemanager.source_pdb_name());
	/// change the pose segment length to that of the source segment length
	core::pose::add_comment(pose, "segment_" + splicemanager.segment_type(), splicemanager.source_pdb_name());//change correct association between current loop and pdb file
	TR<<"DOFS cut site:"<<splicemanager.dofs().cut_site()<<std::endl;
	assign_from_res_to_res(pose);

	TR<<"The dofs size is:"<<splicemanager.dofs().size()<<std::endl;
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
	TR << "pose_from_res()" << splicemanager.pose_from_res() << std::endl;
	TR << "to_res()" << splicemanager.pose_to_res() << std::endl;
	TR << "residue_diff" << splicemanager.residue_diff() << std::endl;
	splicemanager.anchor_res( set_anchor_res() );
	set_fold_tree_nodes(pose);
	splicemanager.fold_tree(pose, fold_tree_nodes_);
	TR << "Foldtree after loop length change: " << pose.fold_tree() << std::endl;
	/// set torsions
	core::Size const total_residue_new(splicemanager.dofs().size()-1); //how long is the introduced segment, I subtract 1 because the size of dofs includes +1 for the name of the segment (see db file)
	//Set BB psi/phi/omega angles on pose segment
	TR << "Changing dofs" << std::endl;
	splicemanager.set_BB_dofs(pose);
	fix_chain_break_residue(pose,splicemanager);

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

	/////////////////////Perform Repacking and Rotamer Trial min mover here on the segment to improve energy ////////////////
	pack(pose,tf);
	rtmin(pose,tf);
	////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////Perform Minimization on the segment ////////////////

	splicemanager.mm()->set_chi(false);
	splicemanager.mm()->set_bb(false);
	//splicemanager.mm()->set_jump(false);
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


	std::vector<std::string> BB_atoms = { "CA", "O", "N", "CB" };
	for ( auto BB_atom : BB_atoms ) { // access by const reference
		splicemanager.add_coordinate_constraints(pose, pose, splicemanager.pose_from_res(),splicemanager.pose_from_res() + total_residue_new - 1,anchor_res, BB_atom);
	}
	//add coordinate constraints to the entire side chain. I can only apply this in cases where the
	splicemanager.add_dihedral_constraints(pose, pose,splicemanager.pose_from_res()+1,splicemanager.pose_from_res() + total_residue_new - 1); //The plus 1 for from res is bcause dihedral constraints are added between i and i-1

	minimize_segment(pose);
	fix_chain_break_residue(pose,splicemanager);
	splicemanager.add_sequence_constraints(pose);
	acb.remove(true);//remove cutpoint variant, this will cause trouble in itterative splicein moves. GDL sep17
	acb.apply(pose);

} //apply

std::string SpliceIn::get_name() const {
	return SpliceInCreator::mover_name();
}

void SpliceIn::parse_my_tag(TagCOP const tag, basic::datacache::DataMap &data,protocols::filters::Filters_map const & /*filters*/,protocols::moves::Movers_map const &/* movers*/,core::pose::Pose const & /*pose*/) {
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
		runtime_assert(database_entry() <= torsion_database_.size());
	}


}

protocols::moves::MoverOP SpliceIn::clone() const {
	return (protocols::moves::MoverOP(new SpliceIn(*this)));
}

core::Size SpliceIn::find_dbase_entry(core::pose::Pose const & pose) {
	core::Size dbase_entry(database_entry());
	if ( first_pass_ ) { /// setup the dbase subset where loop lengths fit the selection criteria
		for ( core::Size i = 1; i <= torsion_database_.size(); ++i ) {// find entries that fit the length criteria
			using namespace protocols::rosetta_scripts;
			ResidueBBDofs const & dofs(torsion_database_[i]);
			splicemanager.dofs(torsion_database_[i]);
			core::pose::Pose tmp_pose = pose;
			assign_from_res_to_res(tmp_pose);
			core::Size const nearest_to_entry_start_on_pose(splicemanager.pose_from_res());
			core::Size const nearest_to_entry_stop_on_pose(splicemanager.pose_to_res());
			core::Size const pose_residues = nearest_to_entry_stop_on_pose - nearest_to_entry_start_on_pose + 1;

			int const delta(dofs.size() - pose_residues);
			//   TR << "The size of pose is " << pose_residues << std::endl;
			//   TR << "Allowed deltas:" << delta_lengths_ << std::endl;
			//   TR << "The current delta is " << delta << std::endl;
			//   TR<< "nearest_to_entry_start_on_pose: "<<nearest_to_entry_start_on_pose<<std::endl;
			//   TR<< "nearest_to_entry_stop_on_pose: "<<nearest_to_entry_stop_on_pose<<std::endl;
			//   TR<< "dofs.start_loop(): "<<dofs.start_loop()<<std::endl;
			//   TR<< "dofs.stop_loop(): "<<dofs.stop_loop()<<std::endl;
			//   TR<< "dofs.size():"<<dofs.size()<<std::endl;
			bool const fit = std::find(delta_lengths_.begin(), delta_lengths_.end(), delta) != delta_lengths_.end();
			if ( fit || database_pdb_entry_ != "" || dbase_entry != 0 ) {
				dbase_subset_.push_back(i);
			}
		}   //for (core::Size i = 1; i <= db.size(); ++i)
		if ( dbase_subset_.empty() ) {
			TR << "Loop of appropriate length not found in database. Returning"
				<< std::endl;
			return 0;
		}
		//TR<<dbase_subset_[0]<<std::endl;
		//get the previous and next segments given the current segment, in case there are other segments
		numeric::random::random_permutation(dbase_subset_.begin(),
			dbase_subset_.end(), numeric::random::rg());
		current_dbase_entry_ = dbase_subset_.begin();
		load_from_checkpoint();
		first_pass_ = false;
	} // fi first_pass
	if ( dbase_iterate() ) {
		TR<<"Iterating through DB"<<std::endl;
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
			end_dbase_subset_->obj = true;
		}
	} else if ( dbase_entry == 0 ) { // fi dbase_iterate
		if ( database_pdb_entry_ == "" ) { //randomize dbase entry
			TR << "The dbase_subset size is " << dbase_subset_.size()<< std::endl;
			core::Size entry_no_to_choose = (core::Size) (numeric::random::rg().uniform()* dbase_subset_.size() + 1);
			if ( dbase_subset_.size() == 0 ) {
				//if ( debug_)
				pose.dump_pdb(splicemanager.mover_name() + "db_size_0.pdb");
				utility_exit_with_message("The database size is "+ utility::to_string(dbase_subset_.size())+ ". Dumped " + splicemanager.mover_name()+ "db_size_0.pdb to inspect");
			}

			TR << "trying to pick entry " << entry_no_to_choose << std::endl;
			dbase_entry = dbase_subset_[entry_no_to_choose];
		} else { // look for the pdb_entry name
			//dbase_entry = ( core::Size )( RG.uniform() * dbase_subset_.size() + 1 );
			for ( core::Size count = 1; count <= dbase_subset_.size(); ++count ) {
				if ( torsion_database_[dbase_subset_[count]].source_pdb()== database_pdb_entry_ ) {
					TR << "Found entry for " << database_pdb_entry_ << " at number " << dbase_subset_[count]  << std::endl;
					dbase_entry = dbase_subset_[count];
					break;
				}
			}
			//This assertion does not make sense. Need to fix this. Gideon 4May15
			//runtime_assert(dbase_entry <= dbase_subset_.size());
		} //else
	} //fi dbase_entry==0
	TR << "Chose db entry: " << dbase_entry << std::endl;
	return dbase_entry;
}

void SpliceIn::load_from_checkpoint() {
	using namespace std;
	if ( checkpointing_file_ == "" ) {
		return;
	}
	utility::io::izstream data(checkpointing_file_);
	if ( !data ) {
		return;
	}
	TR << "Loading from checkpoint" << std::endl;
	/// first read the dbase_subset from the checkpointing file
	{
		string line;
		getline(data, line);
		if ( line.length() == 0 ) {
			TR << "Checkpointing file empty or corrupted. Not loading."
				<< std::endl;
			return;
		}
		istringstream line_stream(line);
		dbase_subset_.clear();
		while ( !line_stream.eof() ) {
			core::Size entry;
			line_stream >> entry;
			dbase_subset_.push_back(entry);
		}
	}
	TR << "dbase subset order loaded from checkpoint is: ";
	for ( core::Size const i : dbase_subset_ ) {
		TR << i << ' ';
	}
	{
		std::string line;
		getline(data, line);
		istringstream line_stream(line);
		core::Size entry;
		line_stream >> entry;
		current_dbase_entry_ = std::find(dbase_subset_.begin(),
			dbase_subset_.end(), entry);
	}
	TR << "current dbase entry loaded from checkpoint is: "
		<< *current_dbase_entry_ << std::endl;
}
void SpliceIn::minimize_segment(core::pose::Pose & pose){
	if ( !min_seg() ) {
		return;
	}

	splicemanager.mm()->show();
	TR<<splicemanager.cut_site()<<std::endl;
	TR<<"Fold tree before minimization: "<<pose.fold_tree()<<std::endl;
	//core::scoring::ScoreFunctionOP scorefxn_with_chainbrk = scorefxn()->clone();
	//scorefxn_with_chainbrk->set_weight( core::scoring::chainbreak, 1.0 );
	protocols::minimization_packing::MinMover min_mover( splicemanager.mm(), scorefxn(), "dfpmin_armijo_nonmonotone", 0.01, true /*use_nblist*/ );

	TR << "scorefxn before min mover " << std::endl;
	scorefxn()->show(pose);
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"before_min_seg.pdb");
	}
	min_mover.apply(pose);
	TR << "scorefxn after min mover " << std::endl;
	scorefxn()->show(pose);
	//remove coordinate constraints post minimization
	pose.remove_constraints();
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"after_min_seg.pdb");
	}
}

void SpliceIn::read_torsion_database() {
	using namespace std;
	TR << "Reading torsion database "<<splicemanager.dbase_file_name() << std::endl;
	utility::io::izstream data(splicemanager.dbase_file_name());
	if ( !data ) {
		utility_exit_with_message("cannot open torsion database " + splicemanager.dbase_file_name() + "\n");
	}
	std::string line;//go through torsion DB line be line
	std::string line_tail;
	while ( getline(data, line) ) {
		utility::vector1<std::string> elements_in_line(utility::string_split(line, ' '));
		if ( elements_in_line.size() % 4 != 0 ) {
			utility_exit_with_message("While reading torsion database " + splicemanager.dbase_file_name()+ " found a line where the number of elements is not divisible by 4. This likely stems from an error in the database:\n"+ line);
		}
		//if we find tail segment dofs we break the line into 2.
		using namespace boost::algorithm;
		trim(line);
		std::istringstream line_stream(line);
		ResidueBBDofs bbdof_entry;
		bbdof_entry.clear();
		//TR<<"Segment dofs:"<<std::endl;
		//TR<<line<<std::endl;
		while ( !line_stream.eof() ) {
			core::Real phi, psi, omega;
			std::string resn;
			line_stream >> phi >> psi >> omega >> resn;
			if ( line_stream.eof() ) { // the end of the line signifies that we're reading the start, stop, cut, source_pdb fields
				bbdof_entry.start_loop((core::Size) phi);
				bbdof_entry.stop_loop((core::Size) psi);
				bbdof_entry.cut_site((core::Size) omega);
				bbdof_entry.source_pdb(resn);
			} else {
				bbdof_entry.push_back(BBDofs(0/*resstd::map<std::string, core::Size> cys_pos;//store all cysteine positions in the AB chainid*/,phi, psi, omega, resn)); /// resid may one day be used. Currently it isn't
			}
		}
		if ( bbdof_entry.stop_loop()!=0 ) { //I am using 0 to mark that this entry is a tail segment entry
			torsion_database_.push_back(bbdof_entry);
			//TR<<"The size of torsion_database_ is:"<<torsion_database_.size()<<std::endl;

		} else {
			tail_torsion_database_.push_back(bbdof_entry);//This is a tail segment
			//TR<<"segment bbdof_entry size is:"<<bbdof_entry.size()<<std::endl;
			//TR<<"The size of tail_torsion_database_ is:"<<tail_torsion_database_.size()<<std::endl;
		}
	}
	TR << "Finished reading torsion database with " << torsion_database_.size() << " entries" << std::endl;
	TR << "Finished reading tail torsion database with " << tail_torsion_database_.size() << " entries" << std::endl;
}

std::string
SpliceIn::mover_name()  {
	return "SpliceIn";
}

std::string SpliceIn_complex_type_name_for_subsubtag( std::string const & foo ) {
	return "subsubtag_splicein_" + foo + "_type";
}

std::string SpliceIn_complex_type_name_for_subtag( std::string const & foo ) {
	return "subtag_splicein_" + foo + "_type";
}


void SpliceIn::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
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
		+ XMLSchemaAttribute::attribute_w_default( "tolerance", xsct_real, "tolernace of bond length std when checking for chain breaks", "0.23" )
		+ XMLSchemaAttribute::attribute_w_default( "ignore_chain_break", xsct_rosetta_bool, "if there is a chain break don't fail", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "debug", xsct_rosetta_bool, "more verbose during run", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "min_seg", xsct_rosetta_bool, "minimize segment after splice in?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "CG_const", xsct_rosetta_bool, "apply coordinate constraints on C-gamma atoms", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rb_sensitive", xsct_rosetta_bool, "apply rigid body addaptations", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "chain_num", xsct_non_negative_integer, "which chain number to apply splice on", "1" )
		+ XMLSchemaAttribute( "segment", xs_string, "segment name" )
		+ XMLSchemaAttribute::attribute_w_default( "superimposed", xsct_rosetta_bool, "are thr structures superimposed", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin", xsct_rosetta_bool, "delete hairpin segment?", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_n", xsct_non_negative_integer, "how many residues to delete on N-ter hairpin", "4" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_c", xsct_non_negative_integer, "how many residues to delete on C-ter hairpin", "13" )
		+ XMLSchemaAttribute( "tail_segment", "n_or_c", "what direction is the tail segment" )
		+ XMLSchemaAttribute( "source_pdb_to_res", xsct_refpose_enabled_residue_number, "residue number of last residue on source segment" )
		+ XMLSchemaAttribute::attribute_w_default( "skip_alignment", xsct_rosetta_bool, "whether to align pose to source pdb", "false" );

	// The "Segments" subtag
	AttributeList segments_subtag_attlist;
	segments_subtag_attlist + XMLSchemaAttribute::attribute_w_default( "profile_weight_away_from_interface", xsct_real, "XRW TO DO", "1.0" )
		+ XMLSchemaAttribute( "current_segment", xs_string, "XRW TO DO" );

	// The "segment" sub-subtag"
	AttributeList subtag_segments_subtag_attlist;
	subtag_segments_subtag_attlist + XMLSchemaAttribute( "pdb_profile_match", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "profiles", xs_string, "XRW TO DO" );

	XMLSchemaComplexTypeGenerator segment_subsubtag_gen;
	segment_subsubtag_gen.complex_type_naming_func( & SpliceIn_complex_type_name_for_subsubtag )
		.element_name( "segment" )
		.description( "individual segment tag" )
		.add_attributes( subtag_segments_subtag_attlist )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList subsubelements;
	subsubelements.add_already_defined_subelement( "segment", SpliceIn_complex_type_name_for_subsubtag/*, 0*/ );

	XMLSchemaComplexTypeGenerator segments_subtag_gen;
	segments_subtag_gen.complex_type_naming_func( & SpliceIn_complex_type_name_for_subtag )
		.element_name( "Segments" )
		.description( "Wrapper for multiple segments tags" )
		.add_attributes( segments_subtag_attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subsubelements )
		.write_complex_type_to_schema( xsd );


	XMLSchemaSimpleSubelementList subelements;
	subelements.add_already_defined_subelement( "Segments", SpliceIn_complex_type_name_for_subtag/*, 0*/ );

	attlist + XMLSchemaAttribute( "use_sequence_profile", xsct_rosetta_bool, "Use sequence profiles in design?" );
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute( "template_file", xs_string, "What is the template file to use during design" );
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist
		+ XMLSchemaAttribute( "design_task_operations", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "residue_numbers_setter", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "torsion_database", xs_string, "torsion db file name" )
		+ XMLSchemaAttribute( "database_entry", xsct_non_negative_integer, "choose db entry by number" )
		+ XMLSchemaAttribute( "database_pdb_entry", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "design_shell", xsct_real, "design shell radius", "6.0" )
		+ XMLSchemaAttribute::attribute_w_default( "repack_shell", xsct_real, "packing shell radius", "8.0" )
		+ XMLSchemaAttribute::attribute_w_default( "randomize_cut", xsct_rosetta_bool, "should the cut be put randomly in the designed segment ", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "cut_secondarystruc", xsct_rosetta_bool, "put cut into ss", "false" )
		+ XMLSchemaAttribute( "delta_lengths", xsct_int_cslist, "comma separated list of allowed segment delta lengths" )
		+ XMLSchemaAttribute::attribute_w_default( "thread_original_sequence", xsct_rosetta_bool, "thread sequence from source pdb onto design", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rtmin", xsct_rosetta_bool, "should we apply rtmin after design", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "dbase_iterate", xsct_rosetta_bool, "should we iterate through the db file", "false" )
		+ XMLSchemaAttribute( "checkpointing_file", xs_string, "name of checkpointing file" );


	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, subelements );
}

void SpliceInCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SpliceIn::provide_xml_schema( xsd );
}
void SpliceIn::set_loop_length_change( protocols::protein_interface_design::movers::LoopLengthChange & llc){
	TR<<"loop start:"<<splicemanager.pose_from_res()<<std::endl;
	TR<<"loop end:"<<splicemanager.cut_site()<<std::endl;
	llc.loop_start(splicemanager.pose_from_res());
	llc.loop_end(splicemanager.cut_site());
	llc.delta(splicemanager.residue_diff());
	llc.tail(false);;
}
void SpliceIn::set_fold_tree_nodes(core::pose::Pose const & pose){
	fold_tree_nodes_.clear();
	splicemanager.mm()->set_jump(false);
	core::conformation::Conformation const & conf(pose.conformation());
	fold_tree_nodes_.push_back({{1,(int) splicemanager.pose_from_res()-1,-1}});
	fold_tree_nodes_.push_back({{(int)splicemanager.pose_from_res()-1,(int) splicemanager.cut_site(),-1}});
	fold_tree_nodes_.push_back({{(int)splicemanager.pose_to_res()+1,(int)conf.chain_end(1),-1}});
	fold_tree_nodes_.push_back({{(int)splicemanager.pose_to_res()+1,(int)splicemanager.cut_site()+1,-1}});
	fold_tree_nodes_.push_back({{(int)splicemanager.pose_from_res()-1,(int)splicemanager.pose_to_res()+1,1}});
	if ( conf.num_chains()>1 ) { //if ligand is present we need to add edge between receptor and ligand
		core::Size CoM = (core::Size ) core::pose::residue_center_of_mass( pose, conf.chain_begin(2), conf.chain_end(2) );
		fold_tree_nodes_.push_back({{(int) conf.chain_end(1),(int) CoM,2}});
		splicemanager.mm()->set_jump(2, true);
	}
}

core::Size SpliceIn::set_anchor_res(){
	return splicemanager.pose_from_res()-1;
}

void SpliceIn::build_ideal_segment(core::pose::Pose & pose){
	using namespace core::chemical;
	using namespace core::conformation;
	ResidueTypeSetCOP residue_set( pose.residue_type_set_for_pose() );
	TR<<splicemanager.pose_from_res()<<std::endl;
	TR<<splicemanager.pose_to_res()<<std::endl;

	pose.delete_residue_range_slow( splicemanager.pose_from_res()+1, splicemanager.pose_to_res()-1);
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_after_delete_segment.pdb");
	}
	core::Size new_to_res = splicemanager.pose_to_res()-(splicemanager.pose_to_res()-splicemanager.pose_from_res());
	ResidueCOP new_res = ResidueFactory::create_residue( residue_set->name_map( name_from_aa( aa_from_oneletter_code( 'A' ) ) ) );
	for ( core::Size res=splicemanager.pose_from_res()+1; res<=splicemanager.cut_site(); ++res ) {
		pose.append_polymer_residue_after_seqpos(*new_res,res-1,true);
		new_to_res++;
	}
	//TR<<"new to res:"<<new_to_res<<std::endl;
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_build_ideal_segment_n-ter.pdb");
	}

	core::Size new_segment_length(splicemanager.pose_to_res()-splicemanager.cut_site()-1);
	for ( core::Size res=1; res<=new_segment_length; res++ ) {
		pose.conformation().safely_prepend_polymer_residue_before_seqpos(*new_res,new_to_res+1,true);
		//TR<<"Appending residue before residue: "<<res+new_to_res-1<<std::endl;
		//pose.dump_pdb("build_ideal_segment.pdb");
	}
}

void SpliceIn::rtmin( core::pose::Pose & pose,core::pack::task::TaskFactoryOP tf){
	using namespace protocols::rosetta_scripts;
	//tf->push_back(core::pack::task::operation::TaskOperationCOP(new core::pack::task::operation::RestrictToRepacking ));
	protocols::minimization_packing::RotamerTrialsMinMover rtmin(scorefxn(), tf);

	utility::vector1<core::Size> packable_residues = residue_packer_states(pose, tf, false, true);

	report_des_pack_res(pose,tf);

	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_before_rtmin.pdb");
	}
	TR << "scorefxn before rtmin " << std::endl;
	scorefxn()->show(pose);
	rtmin.apply(pose);
	TR << "scorefxn after rtmin " << std::endl;
	scorefxn()->show(pose);
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_after_rtmin.pdb");
	}
}


//Since applying the dofs does not actually do any design to the segment I am adding a packer mover//gideonla DEC17
void SpliceIn::pack( core::pose::Pose & pose,core::pack::task::TaskFactoryOP tf){
	using namespace protocols::rosetta_scripts;
	core::pack::task::PackerTaskOP ptask = tf->create_task_and_apply_taskoperations(pose);
	protocols::minimization_packing::PackRotamersMover prm(scorefxn(),ptask);
	prm.apply(pose);
	//pose.dump_pdb(splicemanager.mover_name()+"_after_pack.pdb");
}

} //splice
} //protocols
