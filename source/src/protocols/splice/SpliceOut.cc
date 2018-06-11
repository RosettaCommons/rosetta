// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/splice/SpliceOut.cc
/// @brief
/// @author Gideon Lapidoth (glapidoth@gmail.com)


// Unit headers
#include <core/pose/extra_pose_info_util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/splice/SpliceSegment.hh>
#include <protocols/splice/TailSegmentMover.hh>
#include <protocols/splice/SpliceOut.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <protocols/task_operations/PreventChainFromRepackingOperation.hh>
#include <protocols/splice/SpliceOutCreator.hh>
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
#include <core/id/SequenceMapping.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <basic/Tracer.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/minimization_packing/RotamerTrialsMinMover.hh>
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
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>


namespace protocols {
namespace splice {

using namespace core::conformation;

static  basic::Tracer TR( "protocols.splice.SpliceOut" );
static  basic::Tracer TR_ccd( "protocols.splice.Splice_ccd" );
static  basic::Tracer TR_constraints( "protocols.splice.Splice_constraints" );
static  basic::Tracer TR_pssm( "protocols.splice.Splice_pssm" );
static  basic::Tracer TR_min( "protocols.splice.Splice_min" );
static  basic::Tracer TR_tail( "protocols.splice.Splice_Tail" );
std::string SpliceOutCreator::keyname() const {
	return SpliceOutCreator::mover_name();
}

protocols::moves::MoverOP SpliceOutCreator::create_mover() const {
	return protocols::moves::MoverOP( new SpliceOut );
}

std::string SpliceOutCreator::mover_name() {
	return "SpliceOut";
}


SpliceOut::SpliceOut() :
	Mover(SpliceOutCreator::mover_name()), saved_from_res_(0), saved_to_res_(0), source_pdb_(""), scorefxn_(
	/* NULL */), rms_cutoff_(999999), randomize_cut_(false), cut_secondarystruc_(false), task_factory_( /* NULL */),
	source_pose_(nullptr),saved_fold_tree_( /* NULL */),thread_original_sequence_(false),rtmin_(true), dbase_file_name_(""), splice_filter_( NULL),source_from_res_(0), source_to_res_(0),write_to_database_(true)
{
	splicemanager.profile_weight_away_from_interface(1.0);
	design_shell_ = 6.0;
	repack_shell_ = 8.0;
	basic::options::option[basic::options::OptionKeys::out::file::pdb_comments].value(true);
	splicemanager.rb_sensitive(false);
	tolerance_=0.23;
	delete_hairpin( false );
	delete_hairpin_c( 0 );
	delete_hairpin_n( 0 );
	mover_type_={"MinMover","LoopMover_Refine_CCD","TailSegmentMover"};
	call_mover=nullptr;
	splicemanager.template_file("");
	splicemanager.adjust_stem_positions_by_template_=true;
}


SpliceOut::~SpliceOut() = default;

void SpliceOut::apply(core::pose::Pose & pose) {
	using namespace protocols::rosetta_scripts;
	using core::chemical::DISULFIDE;
	set_last_move_status(protocols::moves::MS_SUCCESS);
	splicemanager.scorefxn(scorefxn());
	///The template pdb is always static. so the segment boundary residue numbers are given according to the template pdb.
	///However, the pose (what we are splicing onto) might have a different length than the template (usually because it
	///already went through a splice mover) so we need to change the user given "from_res" "to_res" (which match those of the template) to match those of the pose.

	splicemanager.adjust_stem_positions_by_template(pose);
	//get the name of the source pdb
	splicemanager.source_pdb_name(parse_pdb_code(source_pdb_));

	//Set source from_res and to_res
	set_source_from_to_res();



	if ( ((source_from_res() == 0) || (source_to_res() == 0)) ) {
		std::ostringstream os; os << "source_from_res: " << source_from_res() << " source_to_res: " << source_to_res() << std::endl;
		utility_exit_with_message(os.str());
	}
	/*before using the conformation from the source pdb we first check that it does not have
	any chain breaks, if it does we kill the run and not use that source*/
	check_source_segment();

	TR << "source_from_res: " << source_from_res() << " source_to_res: " << source_to_res() << std::endl;
	TR << "From_res on pose: "<<splicemanager.pose_from_res()<<",to_res on pose: "<<splicemanager.pose_to_res()<<std::endl;

	splicemanager.residue_diff(source_to_res() - source_from_res() - (splicemanager.pose_to_res() - splicemanager.pose_from_res()));
	TR<<"Residue diff is: "<<splicemanager.residue_diff()<<std::endl;

	//copy the BB phi,psi,omega from the source pdb to the DB
	copy_dofs_from_source();

	/// find a position on the source segment to cut the segment for CCD
	place_cut_site_in_segment(pose);

	/// change the pose segment length to that of the source segment length
	splicemanager.anchor_res( set_anchor_res() );

	superimpose_source_on_pose(pose,*source_pose_);
	protocols::protein_interface_design::movers::LoopLengthChange llc;
	set_loop_length_change(llc);
	TR << "Foldtree before loop length change: " << pose.fold_tree() << std::endl;
	TR<<"cut_site:"<<splicemanager.cut_site()<<std::endl;
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_before_2ndllc_test.pdb");
	}
	llc.apply(pose);
	splicemanager.update_pose_stem_positions();

	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_after_2ndllc_test.pdb");
	}

	build_ideal_segment(pose);

	set_fold_tree_nodes(pose);
	splicemanager.fold_tree(pose, fold_tree_nodes_);

	TR << "Foldtree after loop length change: " << pose.fold_tree() << std::endl;

	/// set torsions
	core::Size const total_residue_new(splicemanager.dofs().size()); //how long is the introduced segment
	//Set BB psi/phi/omega angles on pose segment
	TR << "Changing dofs"<<std::endl;
	splicemanager.set_BB_dofs(pose);
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_after_change_BB_dofs.pdb");
	}

	splicemanager.add_sequence_constraints(pose);

	////////////////////////////////////////////////////////////////////////////////////////////
	///Apply user defined design shell and repack shell around spliced segment
	core::pack::task::TaskFactoryOP tf;
	if ( splicemanager.task_factory() == nullptr ) {
		tf = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory );
	} else {
		tf = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory(*splicemanager.task_factory()) );
	}
	using namespace protocols::task_operations;
	DesignAroundOperationOP dao( new DesignAroundOperation );
	dao->design_shell((splicemanager.task_factory() == NULL ? 0.0 : splicemanager.design_shell())); // threaded sequence operation needs to design, and will restrict design to the loop, unless design_task_factory is defined, in which case a larger shell can be defined
	dao->repack_shell(splicemanager.repack_shell());
	for ( core::Size i = splicemanager.pose_from_res(); i <= splicemanager.pose_to_res(); ++i ) {
		if ( !pose.residue(i).has_variant_type(DISULFIDE) ) {
			dao->include_residue(i);
		}
	}
	tf->push_back(dao);
	////////////////////////////////////////////////////////////////////////////////////////////
	//Add chainbreak attribute to the cut site resiudes.
	protocols::protein_interface_design::movers::AddChainBreak acb;
	acb.resnum(utility::to_string(splicemanager.cut_site()));
	acb.find_automatically(false);
	acb.change_foldtree(false);
	TR << "Adding ccd chainbreak at: " << splicemanager.cut_site() << std::endl;
	acb.apply(pose);
	////////////////////////////////////////////////////////////////////////////////////////////
	core::kinematics::MoveMapOP mm;
	mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap );
	mm->set_chi(false);
	mm->set_bb(false);
	mm->set_jump(false);
	core::conformation::Conformation const & conf(pose.conformation());
	core::kinematics::FoldTree ft;
	ft=pose.fold_tree();
	if ( conf.num_chains()>1 ) { //if ligand is present we need to add edge between receptor and ligand
		core::Size const jump_num= ft.num_jump();
		TR<<"Number of jumps in the fold tree is:"<<jump_num<<std::endl;
		mm->set_jump(jump_num, true); /// 1Feb13 for cases in which the we're splicing in the presence of a ligand

	}
	for ( core::Size i = splicemanager.pose_from_res(); i <= splicemanager.pose_to_res(); ++i ) {
		mm->set_chi(i, true);
		mm->set_bb(i, true);
	}
	splicemanager.mm(mm);
	//Apply coordinate and dihedral constraints to spliced segment
	TR<<"Anchor res before applying coordinate constraints: "<<splicemanager.anchor_res()<<std::endl;
	std::vector<std::string> BB_atoms = {"CA","O","N","CB"};
	for ( auto BB_atom : BB_atoms ) { // access by const reference
		splicemanager.add_coordinate_constraints(pose, *source_pose_, splicemanager.pose_from_res(), splicemanager.pose_from_res()+total_residue_new-1, splicemanager.anchor_res(),BB_atom);
	}
	//add coordinate constraints to the entire side chain. I can only apply this in cases where the
	if ( CG_const_ ) {
		core::pack::task::PackerTaskOP ptask_for_coor_const(tf->create_task_and_apply_taskoperations(pose));
		splicemanager.add_coordinate_constraints(pose, *source_pose_, splicemanager.pose_from_res(), splicemanager.pose_from_res()+total_residue_new-1, splicemanager.anchor_res(),"CG",ptask_for_coor_const); //add coordiante constraints to loop
		TR_constraints<<"score function after CG constraints"<<std::endl;
	}
	splicemanager.add_dihedral_constraints(pose, *source_pose_, splicemanager.pose_from_res()+1, splicemanager.pose_from_res()+total_residue_new-1); //The plus 1 for from res is bcause dihedral constraints are added between i and i-1
	////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1<core::Size> designable_residues = residue_packer_states(pose, tf, true, false);
	TR << "Residues Allowed to Design: "<< std::endl;
	for ( utility::vector1<core::Size>::const_iterator i(designable_residues.begin());
			i != designable_residues.end(); ++i ) {
		TR << *i << "+";
	}
	TR << std::endl;

	utility::vector1<core::Size> packable_residues = residue_packer_states(pose, tf, true, true);
	TR << "Residues allowed to repack: "<< std::endl;
	for ( utility::vector1<core::Size>::const_iterator i(packable_residues.begin());
			i != packable_residues.end(); ++i ) {
		TR << *i << "+";
	}
	TR << std::endl;

	TR<<"Using "<< TR.Red << submover_->get_name()<<TR.Reset<<" to close loop"<<std::endl;
	splicemanager.task_factory(tf);//update the private member taskoperation
	(this->*call_mover)(pose,mm);
	acb.resnum(utility::to_string(splicemanager.cut_site()));
	acb.apply(pose);
	minimize_segment(pose);
	pose.remove_constraints();//remove coor/dih constraints before rtmin
	splicemanager.add_sequence_constraints(pose);

	if ( (SpliceOutRMSDFilter(&pose)) || SpliceOutFilter(&pose) ) {
		return;
	}

	/// tell us what the torsions of the new (closed) loop are.
	core::pack::task::TaskFactoryOP tf_dofs( new core::pack::task::TaskFactory );
	protocols::task_operations::DesignAroundOperationOP dao_dofs( new protocols::task_operations::DesignAroundOperation );
	for ( core::Size i = splicemanager.pose_from_res(); i <= splicemanager.pose_to_res(); ++i ) {
		dao_dofs->include_residue(i);
	}
	dao_dofs->design_shell(0);  /// only include the loop residues
	tf_dofs->push_back(dao_dofs);
	protocols::protein_interface_design::filters::Torsion torsion;
	torsion.task_factory(tf_dofs);
	torsion.task_factory_set(true);
	torsion.apply(pose);

	/// Now write to dbase disk file
	write_database_to_file(pose);
	//TR<<core::pose::get_all_comments(pose)<<std::endl;
	// remove_hairpin( pose );
	saved_fold_tree_ = core::kinematics::FoldTreeOP( new core::kinematics::FoldTree(pose.fold_tree()) );
	//retrieve_values();

}//apply

void
SpliceOut::write_database_to_file(core::pose::Pose const & pose){
	if ( !(write_to_database_) ) {
		return;
	}
	TR<<"Writing to database"<<std::endl;
	std::ofstream dbase_file;
	std::stringstream torsion_stringstream; // We cannot write directly into the dbase_file, as the assert statement might stop the run, and thus leave the db poluted, Chris Sep23
	dbase_file.open(splicemanager.dbase_file_name().c_str(), std::ios::app);
	for ( core::Size i = splicemanager.pose_from_res(); i <= splicemanager.pose_to_res(); ++i ) {
		if ( i != splicemanager.pose_to_res()&&i != splicemanager.pose_from_res() ) { // don't do the check for the last residue, which would have psi and omega equals 0 for tail segments
			//TR<<"Checking residue:"<<i<<std::endl;
			runtime_assert( pose.phi(i) && pose.psi(i) && pose.omega(i) );// Make sure that all dihedral angles have non zero values,
		}
		//0 value dihedral angle is probably an error, Gideon Sep14
		torsion_stringstream << pose.phi(i) << ' ' << pose.psi(i) << ' ' << pose.omega(i) << ' ' << pose.residue(i).name3() << ' ';
	}
	torsion_stringstream << splicemanager.template_from_res() << ' ' << splicemanager.template_to_res() << ' ' << splicemanager.cut_site() << ' ';
	torsion_stringstream << splicemanager.source_pdb_name()<<"\n";
	dbase_file << torsion_stringstream.str();
	dbase_file.close();
}


void
SpliceOut::remove_hairpin( core::pose::Pose & pose ) const{
	if ( !delete_hairpin() ) {
		return;
	}

	using namespace protocols::rosetta_scripts;
	using core::Size;
	using core::Real;

	Size const start = find_nearest_disulfide( pose, pose.conformation().chain_end( 1 ) ) + delete_hairpin_n();
	Size const end   = pose.conformation().chain_end( 1 ) - delete_hairpin_c();
	runtime_assert( start < end );
	for ( Size i = start + 1; i <= end - 1; ++i ) {
		pose.delete_polymer_residue( start + 1 );
	}
}

/// splice apply might change the from_res/to_res internals since they sometimes refer to the template file. If that happens, we want the values to
/// revert to their original values before the end of the apply function (so retrieve_values) below must be called before return.
void SpliceOut::save_values() {
	saved_from_res_ = splicemanager.pose_from_res();
	saved_to_res_ = splicemanager.pose_to_res();
}

void SpliceOut::retrieve_values() {
	splicemanager.pose_from_res(saved_from_res_);
	splicemanager.pose_to_res(saved_to_res_);
}

std::string SpliceOut::get_name() const {
	return SpliceOutCreator::mover_name();
}

//All tags that would be used by all SpliceOut derived classes are here
void SpliceOut::parse_SpliceOut_tags(TagCOP const tag,protocols::moves::Movers_map const & movers,protocols::filters::Filters_map const & filters){
	CG_const_=tag->getOption<bool>("CG_const", false);
	superimposed(tag->getOption< bool >( "superimposed", true ) );
	if ( delete_hairpin() ) {
		delete_hairpin_n( tag->getOption< core::Size >( "delete_hairpin_n", 4 ) );
		delete_hairpin_c( tag->getOption< core::Size >( "delete_hairpin_c", 13 ) );
		TR<<"deleting the hairpin with parameters delete_hairpin_n: "<<delete_hairpin_n()<<" delete_hairpin_c: "<<delete_hairpin_c()<<std::endl;
	}
	handle_mover_tag(tag,movers);
	if ( tag->hasOption("source_pdb") ) {
		source_pdb(tag->getOption<std::string>("source_pdb"));
		source_pose_ = core::pose::PoseOP( new core::pose::Pose );
		source_pose_ = core::import_pose::pose_from_file( source_pdb_);
	}


	delete_hairpin( tag->getOption< bool >( "delete_hairpin", false ) );
	if ( tag->hasOption("splice_filter") ) {
		splice_filter(protocols::rosetta_scripts::parse_filter(tag->getOption<std::string>("splice_filter"), filters));
	}
	rms_cutoff(tag->getOption<core::Real>("rms_cutoff", 999999));
	rms_cutoff_loop(tag->getOption<core::Real>("rms_cutoff_loop", -1));//Added by gideonla Sep15, used in concatenation with the "rms_cutoff" sets a different rms cutoff for loop segments
	ignore_chain_break_=tag->getOption<bool>("ignore_chain_break",false);//for debugging purposes
}

//Modify this function in derived classes to change the specific flavour of parse_my_tag
void SpliceOut::abstract_parse_tag(TagCOP const tag){
	randomize_cut(tag->getOption<bool>("randomize_cut", false));
	splicemanager.cut_site(tag->getOption<core::Size>("cut_site",0));
	runtime_assert((tag->hasOption("randomize_cut") && tag->hasOption("source_pose")) || !tag->hasOption("source_pose"));
	cut_secondarystruc(tag->getOption<bool>("cut_secondarystruc", false));
	tolerance_=tag->getOption<core::Real>("tolerance",0.23);//for debugging purposes
	splicemanager.allowed_cuts(tag->getOption<core::Size>("allowed_cuts",1));
	splicemanager.template_from_res(core::pose::parse_resnum(tag->getOption<std::string>("from_res", "0"), *(splicemanager.template_pose())));
	splicemanager.template_to_res(core::pose::parse_resnum(tag->getOption<std::string>("to_res", "0"), *(splicemanager.template_pose())));

}


void SpliceOut::parse_my_tag(TagCOP const tag, basic::datacache::DataMap &data, protocols::filters::Filters_map const & filters, protocols::moves::Movers_map const & movers, core::pose::Pose const & /*pose*/) {
	utility::vector1<TagCOP> const sub_tags(tag->getTags());

	parse_SpliceOut_tags(tag,movers,filters);
	splicemanager.parse_segments(sub_tags,tag, data);
	splicemanager.parse_tags(tag,data);
	abstract_parse_tag(tag);
	scorefxn(protocols::rosetta_scripts::parse_score_function(tag, data));
	if ( !superimposed() ) {
		source_from_res(core::pose::parse_resnum(tag->getOption<std::string>("source_from_res", "0"), *source_pose_));
		source_to_res(core::pose::parse_resnum(tag->getOption<std::string>("source_to_res", "0"), *source_pose_));
		utility_exit_with_message("You have noted that structures are not aligned and have not given source_from_res and source_to_res, fix XML \n");
	}
}

protocols::moves::MoverOP SpliceOut::clone() const {
	return (protocols::moves::MoverOP(new SpliceOut(*this)));
}

void SpliceOut::scorefxn(core::scoring::ScoreFunctionOP sf) {
	scorefxn_ = sf;
}

core::scoring::ScoreFunctionOP SpliceOut::scorefxn() const {
	return scorefxn_;
}


/// the torsion dbase should have the following structure:
/// each line represents a single loop. Each four values represent <phi> <psi> <omega> <3-let resid>; the last entry in a line represents <loop start> <loop stop> <cut site> cut; where cut signifies that this is the loop designator

protocols::filters::FilterOP SpliceOut::splice_filter() const {
	return splice_filter_;
}

void SpliceOut::splice_filter(protocols::filters::FilterOP f) {
	splice_filter_ = f;
}



void SpliceOut::load_pdb_segments_from_pose_comments(core::pose::Pose const & pose) {
	// if(use_sequence_profiles_){
	//If we are using sequence profiles then the condition is true and function can run
	using namespace std;
	map<string, string> const comments = core::pose::get_all_comments(pose);
	TR << "The size of comments is: " << comments.size() << std::endl;
	for ( std::map<string, string>::const_iterator i = comments.begin(); i != comments.end(); ++i ) {
		//TR<<"the size of j is: "<<j<<std::endl;
		std::string const key(i->first);
		//TR<<"the size of j after i->first is: "<<j<<std::endl;
		std::string const val(i->second);
		//TR<<"the size of j after i->second is: "<<j<<std::endl;
		if ( key.substr(0, 7) != "segment" ) { /// the expected format is segment_??, where we're interested in ??
			continue;
		}
		std::string const short_key(key.substr(8, 1000));
		splicemanager.pdb_segments()[short_key] = val;
		TR << "recording segment/pdb pair: " << short_key << '/' << val << std::endl;
	}
	// }
}



void SpliceOut::check_pose_pssm_match(core::pose::PoseOP pose, utility::vector1<core::sequence::SequenceProfileOP> profile_vector) {
	core::Size aapos = 0;
	//TR<<"TESTING PSSMs"<<std::endl;
	for ( core::Size seg = 1; seg <= profile_vector.size(); seg++ ) { //go over all the PSSM segments provided by the user
		for ( core::Size pos /*go over profile ids*/= 1; pos <= profile_vector[seg]->size(); ++pos ) {
			++aapos;  //go over pose residue
			TR_pssm << pose->residue(aapos).name1() << aapos << ","<< profile_vector[seg]->prof_row(pos) << std::endl;
			if ( (profile_vector[seg]->prof_row(pos)[2]) > 8 ) { //If the profile vector holds a disulfide Cys it will have a pssm score over 8
				std::stringstream ss;
				std::string s;
				ss << pose->residue(aapos).name1();
				ss >> s;
				//TR<<"found a dis cys="<<s<<std::endl;
				if ( s.compare("C") != 0 ) {
					std::string seqpos;
					std::ostringstream convert;
					convert << aapos; // insert the textual representation of 'Number' in the characters in the stream
					seqpos = convert.str();
					pose->dump_pdb( splicemanager.source_pdb_name()+"_align_problem.pdb");
					utility_exit_with_message(" PSSM and pose might be misaligned, position " + s + seqpos + " should be a CYS\n");
				} //fi
			} //fi
		} //end inner segment for
	} //end pssm segment for
}


///@brief apply coordinate constraints on the segment being inserted. "to" and "from" are residue number of the pose(!), anchor residue number is also on the pose



void SpliceOut::dbase_file_name(std::string const & s) {
	dbase_file_name_ = s;
}

std::string SpliceOut::dbase_file_name() const {
	return dbase_file_name_;
}



/// @brief Since we want to minimally perturb the active site confirmation (whether it be binding or catalytic) the cut site should be placed farthest away.
/// For each protein family we will probably need special definitions unless we find a more general way. For antibodies I go from the conserved Trp res at the base of CDR1
/// and then continue going to the C-ter until I find the distal loops.
core::Size SpliceOut::find_non_active_site_cut_site(core::pose::Pose const & pose) {
	using namespace core::sequence;
	SequenceProfileOP profile;
	TR<<"Placing cut away from functional site"<<std::endl;
	std::string const source_pdb_name(parse_pdb_code(pose.pdb_info()->name()));
	//use pssm to find conserved trp
	if ( splicemanager.splice_segments()[ splicemanager.segment_type() ]->pdb_profile( source_pdb_name )==0 ) {
		utility_exit_with_message(" could not find the source pdb name: "+ source_pdb_name + ", in pdb_profile_match file."+splicemanager.segment_type()+" or PSSM file is missing\n");
	} else {
		profile=splicemanager.splice_segments()[ splicemanager.segment_type() ]->pdb_profile( source_pdb_name );
	}

	core::Size aapos=1;
	for ( core::Size pos /*go over profile ids*/= 1; pos <= profile->size(); ++pos ) {
		TR_pssm << pose.residue(aapos).name1() << aapos << ","<< profile->prof_row(pos) << std::endl;
		if ( (profile->prof_row(pos)[19]) > 9 ) { //Conserved CDR1 stem will have a pssm score over 8
			break;
		}
		++aapos;//found conserved trp
	} //for
	core::scoring::dssp::Dssp dssp(pose);
	dssp.dssp_reduced();// switch to simplified H E L notation
	for ( core::Size pos /*go over profile ids*/= aapos; pos <= pose.total_residue(); ++pos ) {
		if ( dssp.get_dssp_secstruct(pos) == 'L' ) { // allow site for cutting if it's either in a loop or if cutting secondary structure is allowed
			while ( pose.residue(pos).name3()=="PRO"||pose.residue(pos+1).name3()=="PRO" )
					pos=pos+1;//Can't place cut site after proline
			TR<<"Found cut site at:"<<pos<<pose.residue(pos).name1()<<std::endl;
			return pos+2;
		}
	}

	return 0;
}

void SpliceOut::check_source_segment() {
	TR<<"Checking source segment for breaks"<<std::endl;
	protocols::simple_moves::CutChainMover ccm;
	ccm.bond_length(2.5);
	core::Size source_pdb_cut(ccm.chain_cut(*source_pose_,source_from_res(), source_to_res()));
	//TR<<"Source PDB cut: "<<source_pdb_cut<<std::endl;
	//TR<<"Ignore chain break val: "<<ignore_chain_break_<<std::endl;
	if ( source_pdb_cut!=0 && !ignore_chain_break_ ) { //if there is a chain break in the source segment we don't use that source segment
		std::ostringstream os; os <<source_pdb_cut;
		utility_exit_with_message("found chain break in source PDB "+source_pdb_+" at "+os.str()+", exiting\n");
	}
}

void SpliceOut::copy_dofs_from_source(){
	splicemanager.dofs().clear();
	for ( core::Size i = source_from_res(); i <=  source_to_res(); ++i ) {
		/// Feed the source_pose dofs into the BBDofs array
		BBDofs residue_dofs;
		residue_dofs.resid(i); /// resid is probably never used
		residue_dofs.phi(source_pose_->phi(i));
		residue_dofs.psi(source_pose_->psi(i));
		residue_dofs.omega(source_pose_->omega(i));
		/// convert 3let residue code to 1let code
		std::stringstream ss;
		std::string s;
		ss << source_pose_->residue(i).name1();
		ss >> s;
		residue_dofs.resn(s);
		splicemanager.dofs().push_back(residue_dofs);
	} // for i source_from_res..source_to_res
}

/// choose cutsite randomly within loop residues on the loop (no 2ary structure)
void SpliceOut::place_cut_site_in_segment(core::pose::Pose const & pose){
	if ( splicemanager.cut_site()!=0 ) {
		return;// If cut_site is already defined than don't randomly set
	}
	core::scoring::dssp::Dssp dssp(*source_pose_);
	dssp.dssp_reduced(); // switch to simplified H E L notation
	std::vector<core::Size> loop_positions_in_source;
	loop_positions_in_source.clear();
	TR << "DSSP of source segment: ";
	core::Size last_res( std::min(source_to_res(), splicemanager.pose_to_res() - splicemanager.pose_from_res() + source_from_res()) );
	core::Size first_res( source_from_res() );

	for ( core::Size i = first_res; i <= last_res; ++i ) {
		if ( ( dssp.get_dssp_secstruct(i) == 'L' || cut_secondarystruc() ) && pose.conformation().residue(i).name1() != 'P'  && pose.conformation().residue(i+1).name1() != 'P' ) { // allow site for cutting if it's either in a loop or if cutting secondary structure is allowed. Don't allow cuts next to prolines... cause trouble in calculating proline closure score.
			loop_positions_in_source.push_back(i);
		}
		TR << dssp.get_dssp_secstruct(i);
	}
	TR << std::endl;
	//New test to see what is the sequence of the new loop
	TR.Debug << "The sequence of the source loop is: ";
	for ( core::Size i = first_res; i <= last_res; ++i ) {
		TR.Debug << source_pose_->residue(i).name1() << " ";
	}
	TR << std::endl;

	splicemanager.cut_site(loop_positions_in_source[(core::Size) (numeric::random::rg().uniform() * loop_positions_in_source.size())]- source_from_res() + splicemanager.pose_from_res());
	TR << "Cut site on source PDB: "<< splicemanager.cut_site()+ source_from_res() -splicemanager.pose_from_res()<<std::endl;
	TR << "Cut placed at: " << splicemanager.cut_site() << std::endl;
}

void SpliceOut::superimpose_source_on_pose( core::pose::Pose const & pose, core::pose::Pose & source_pose ){
	using namespace protocols::toolbox;

	if ( superimposed() ) { //Structures are superimposed. Nothing to do
		return;
	}
	utility::vector1< core::Size > pose_positions, template_positions;
	pose_positions.clear(); template_positions.clear();
	pose_positions.push_back( source_from_res() -1 );
	pose_positions.push_back( source_from_res() );
	pose_positions.push_back( source_from_res() +1 );
	pose_positions.push_back( source_to_res() -1 );
	pose_positions.push_back( source_to_res() );
	pose_positions.push_back( source_to_res() +1 );
	template_positions.push_back( splicemanager.template_from_res() -1 );
	template_positions.push_back( splicemanager.template_from_res() );
	template_positions.push_back( splicemanager.template_from_res() +1 );
	template_positions.push_back( splicemanager.template_to_res() -1 );
	template_positions.push_back( splicemanager.template_to_res() );
	template_positions.push_back( splicemanager.template_to_res() +1 );
	TR<<" template scafold_res: "<<splicemanager.template_from_res() -1<<",source scafold_res: "<< source_from_res() -1<<std::endl;
	utility::vector1< numeric::xyzVector< core::Real > > init_coords( coords( source_pose, pose_positions ) ), ref_coords( coords( pose/*this is the starting pose, the structure on which we want to graft the loop from the source protein*/, template_positions ));
	TR<<"template ref coords: "<<ref_coords[1][1]<<std::endl;
	TR<<"source ref coords: "<<init_coords[1][1]<<std::endl;

	numeric::xyzMatrix< core::Real > rotation;
	numeric::xyzVector< core::Real > to_init_center, to_fit_center;

	superposition_transform( init_coords, ref_coords, rotation, to_init_center, to_fit_center );

	apply_superposition_transform( source_pose, rotation, to_init_center, to_fit_center );
	/// DEBUGGING
	if ( splicemanager.debug() ) {
		pose.dump_pdb("pose_pdb.pdb");
		source_pose.dump_pdb("source_pose_pdb.pdb");
	}
}
void SpliceOut::ccd_mover( core::pose::Pose & pose,core::kinematics::MoveMapOP mm){
	using namespace protocols::loops;
	protocols::loops::loop_mover::refine::LoopMover_Refine_CCD*  ccd_mover = dynamic_cast< protocols::loops::loop_mover::refine::LoopMover_Refine_CCD* >( submover_.get() );
	TR_ccd<<"Loop definition before ccd:startn, startc, cut_site"<<splicemanager.pose_from_res()<<","<< splicemanager.pose_to_res()<<","<< splicemanager.cut_site()<<std::endl;
	Loop loop(splicemanager.pose_from_res(), splicemanager.pose_to_res(), splicemanager.cut_site()); /// Gideon & Sarel (8Jul13): we're now respecting the user's choice of from_res to_res and not melting the framework
	LoopsOP loops( new Loops() );
	loops->push_back(loop);
	ccd_mover->loops(loops);
	TR_ccd<<"fold tree before ccd"<<pose.fold_tree()<<std::endl;
	core::pack::task::TaskFactoryOP tf;
	tf = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory(*splicemanager.task_factory()) );
	ccd_mover->set_task_factory(tf);
	ccd_mover->move_map(mm);
	//ccd_mover->set_scorefxn(scorefxn());
	mm->show();
	TR_ccd<<"Score function before CCD"<<std::endl;
	scorefxn()->show(pose);
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_before_ccd.pdb");
	}
	ccd_mover->apply(pose);
	TR << "Score function after CCD" << std::endl;
	scorefxn()->show(pose);
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_after_ccd.pdb");
	}

}
void SpliceOut::min_mover( core::pose::Pose & pose,core::kinematics::MoveMapOP mm){
	protocols::minimization_packing::MinMover*  min_mover = dynamic_cast< protocols::minimization_packing::MinMover* >( submover_.get() );
	min_mover->set_movemap(mm);
	TR_min<<"fold tree before MinMover"<<pose.fold_tree()<<std::endl;
	TR_min<<"Score function before MinMover"<<std::endl;
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_before_minmover.pdb");
	}
	//min_mover->score_function(scorefxn());
	min_mover->apply(pose);//if tail segment is on we don't need ccd
	TR_min << "Score function after MinMover" << std::endl;
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_after_minmover.pdb");
	}
}

void SpliceOut::tail_mover( core::pose::Pose & pose,core::kinematics::MoveMapOP mm){
	TailSegmentMover*  tail_mover = dynamic_cast< TailSegmentMover* >( submover_.get() );
	mm->show();
	tail_mover->set_movemap(mm);
	scorefxn()->show(pose);
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_before_tailmover.pdb");
	}
	core::pack::task::TaskFactoryOP tf;
	tf = core::pack::task::TaskFactoryOP( new core::pack::task::TaskFactory(*splicemanager.task_factory()) );
	tail_mover->set_task_factory(tf);
	//tail_mover->set_fa_scorefxn(scorefxn());
	tail_mover->get_scorefunction();
	tail_mover->apply(pose);//if tail segment is on we don't need ccd
	TR_tail << "Score function after tail segment mover" << std::endl;
	scorefxn()->show(pose);
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_after_tailmover.pdb");
	}
}

void SpliceOut::handle_mover_tag(TagCOP const tag,protocols::moves::Movers_map const & movers){
	std::string mover_name("");
	if ( tag->hasOption("mover") ) {
		mover_name = tag->getOption< std::string >( "mover" );
	} else {
		utility_exit_with_message("Must specify mover to copy source conformation onto template segment\n");
	}
	protocols::moves::Movers_map::const_iterator  find_mover ( movers.find( mover_name ));
	if ( find_mover == movers.end() ) {
		utility_exit_with_message("Mover \""+mover_name+"\" not found in mover section in the XML, check XML\n");
	}
	submover_ = find_mover->second;
	TR<<submover_->get_name()<<std::endl;
	if ( std::find(mover_type_.begin(), mover_type_.end(), submover_->get_name()) == mover_type_.end() ) {
		utility_exit_with_message("Please choose only \"Minmover\" or \"LoopMover_Refine_CCD\" for \"mover=\" tag \n");
	} else {
		if ( submover_->get_name()=="MinMover" ) {
			call_mover =&SpliceOut::min_mover;
		} else if ( submover_->get_name()=="LoopMover_Refine_CCD" ) {
			call_mover = &SpliceOut::ccd_mover;
		} else if ( submover_->get_name()=="TailSegmentMover" ) {
			call_mover = &SpliceOut::tail_mover;
		}
	}

}

void SpliceOut::set_fold_tree_nodes(core::pose::Pose const & pose){
	fold_tree_nodes_.clear();
	core::conformation::Conformation const & conf(pose.conformation());
	fold_tree_nodes_.push_back({{1,(int) splicemanager.pose_from_res()-1,-1}});
	fold_tree_nodes_.push_back({{(int)splicemanager.pose_from_res()-1,(int) splicemanager.cut_site(),-1}});
	fold_tree_nodes_.push_back({{(int)splicemanager.pose_to_res()+1,(int)conf.chain_end(1),-1}});
	fold_tree_nodes_.push_back({{(int)splicemanager.pose_to_res()+1,(int)splicemanager.cut_site()+1,-1}});
	fold_tree_nodes_.push_back({{(int)splicemanager.pose_from_res()-1,(int)splicemanager.pose_to_res()+1,1}});
	if ( conf.num_chains()>1 ) { //if ligand is present we need to add edge between receptor and ligand
		core::Size CoM = (core::Size ) core::pose::residue_center_of_mass( pose, conf.chain_begin(2), conf.chain_end(2) );
		fold_tree_nodes_.push_back({{(int) conf.chain_end(1),(int) CoM,2}});
	}
}


void SpliceOut::set_source_from_to_res(){
	TR<<"set source from"<<std::endl;
	using namespace protocols::rosetta_scripts;
	source_from_res(find_nearest_res(*source_pose_, *splicemanager.template_pose(), splicemanager.template_from_res(), 0/*chain*/));
	source_to_res(find_nearest_res(*source_pose_, *splicemanager.template_pose(), splicemanager.template_to_res(), 0/*chain*/));
}


core::Size SpliceOut::set_anchor_res(){
	return splicemanager.pose_from_res()-1;
}
void SpliceOut::set_loop_length_change( protocols::protein_interface_design::movers::LoopLengthChange & llc){
	llc.loop_start(splicemanager.pose_from_res());
	llc.loop_end(splicemanager.cut_site());
	llc.delta(splicemanager.residue_diff());
}

std::string SpliceOut_complex_type_name_for_subsubtag( std::string const & foo ) {
	return "subsubtag_spliceout_" + foo + "_type";
}

std::string SpliceOut_complex_type_name_for_subtag( std::string const & foo ) {
	return "subtag_spliceout_" + foo + "_type";
}


void SpliceOut::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;


	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::attribute_w_default( "tolerance", xsct_real, "XRW TO DO", "0.23" )
		+ XMLSchemaAttribute::attribute_w_default( "ignore_chain_break", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "debug", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "CG_const", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rb_sensitive", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "chain_num", xsct_non_negative_integer, "XRW TO DO", "1" )
		+ XMLSchemaAttribute::attribute_w_default( "cut_site", xsct_non_negative_integer, "residue number of where to place cut", "1" )
		+ XMLSchemaAttribute( "Segment", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "superimposed", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_n", xsct_non_negative_integer, "XRW TO DO", "4" )
		+ XMLSchemaAttribute::attribute_w_default( "delete_hairpin_c", xsct_non_negative_integer, "XRW TO DO", "13" )
		+ XMLSchemaAttribute( "source_pdb_to_res", xsct_refpose_enabled_residue_number, "XRW TO DO" );

	// The "Segments" subtag
	AttributeList segments_subtag_attlist;
	segments_subtag_attlist + XMLSchemaAttribute( "current_segment", xs_string, "XRW TO DO" );

	// The "segment" sub-subtag"
	AttributeList subtag_segments_subtag_attlist;
	subtag_segments_subtag_attlist + XMLSchemaAttribute( "pdb_profile_match", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "profiles", xs_string, "XRW TO DO" );

	XMLSchemaComplexTypeGenerator segment_subsubtag_gen;
	segment_subsubtag_gen.complex_type_naming_func( & SpliceOut_complex_type_name_for_subsubtag )
		.element_name( "Segment" )
		.description( "individual segment tag" )
		.add_attributes( subtag_segments_subtag_attlist )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );

	XMLSchemaSimpleSubelementList subsubelements;
	subsubelements.add_already_defined_subelement( "Segment", SpliceOut_complex_type_name_for_subsubtag/*, 0*/ );

	XMLSchemaComplexTypeGenerator segments_subtag_gen;
	segments_subtag_gen.complex_type_naming_func( & SpliceOut_complex_type_name_for_subtag )
		.element_name( "Segments" )
		.description( "Wrapper for multiple segments tags" )
		.add_attributes( segments_subtag_attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( subsubelements )
		.write_complex_type_to_schema( xsd );


	XMLSchemaSimpleSubelementList subelements;
	subelements.add_already_defined_subelement( "Segments", SpliceOut_complex_type_name_for_subtag/*, 0*/ );

	attlist + XMLSchemaAttribute( "use_sequence_profile", xsct_rosetta_bool, "XRW TO DO" );
	protocols::rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "add_sequence_constraints_only", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "template_file", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "set_fold_tree_only", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute( "source_pdb", xs_string, "XRW TO DO");
	protocols::rosetta_scripts::attributes_for_parse_task_operations( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "from_res", xsct_refpose_enabled_residue_number, "XRW TO DO", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "to_res", xsct_refpose_enabled_residue_number, "XRW TO DO", "0" )
		+ XMLSchemaAttribute( "design_task_operations", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "residue_numbers_setter", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "torsion_database", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute::attribute_w_default( "design_shell", xsct_real, "XRW TO DO", "6.0" )
		+ XMLSchemaAttribute::attribute_w_default( "repack_shell", xsct_real, "XRW TO DO", "8.0" )
		+ XMLSchemaAttribute::attribute_w_default( "rms_cutoff", xsct_real, "XRW TO DO", "999999" )
		+ XMLSchemaAttribute::attribute_w_default( "rms_cutoff_loop", xsct_real, "XRW TO DO", "999999" )
		+ XMLSchemaAttribute::attribute_w_default( "res_move", xsct_non_negative_integer, "XRW TO DO", "1000" )
		+ XMLSchemaAttribute::attribute_w_default( "randomize_cut", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "cut_secondarystruc", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "thread_ala", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "design", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "thread_original_sequence", xsct_rosetta_bool, "XRW TO DO", "false" )
		+ XMLSchemaAttribute::attribute_w_default( "rtmin", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute( "locked_residue", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "checkpointing_file", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "splice_filter", xs_string, "XRW TO DO" )
		+ XMLSchemaAttribute( "mover", xs_string, "Which mover to use to close the segment" )
		+ XMLSchemaAttribute::attribute_w_default( "restrict_to_repacking_chain2", xsct_rosetta_bool, "XRW TO DO", "true" )
		+ XMLSchemaAttribute::attribute_w_default( "use_sequence_profiles", xsct_rosetta_bool, "XRW TO DO", "true" );

	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "XRW TO DO", attlist, subelements );

}
void SpliceOutCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SpliceOut::provide_xml_schema( xsd );
}
std::string
SpliceOut::mover_name()  {
	return "SpliceOut";
}

int
SpliceOut::SpliceOutFilter(core::pose::Pose * pose)  {
	if ( splice_filter()==NULL ) {
		return 0;
	}
	protocols::protein_interface_design::movers::AddChainBreak acb;
	acb.resnum(utility::to_string(splicemanager.cut_site()));
	acb.find_automatically(false);
	acb.change_foldtree(false);
	TR << "Adding ccd chainbreak at: " << splicemanager.cut_site() << std::endl;
	std::string Result_filter;
	std::ostringstream convert_filter;
	acb.apply(*pose);
	TR<<"grep  "<<pose->fold_tree()<<std::endl;
	convert_filter << splice_filter()->score(*pose);
	convert_filter << splice_filter()->report_sm(*pose);
	Result_filter = convert_filter.str();
	TR<< name_for_filter()+"Filter Val:"<< Result_filter<<std::endl;
	core::pose::add_comment(*pose, name_for_filter()+"Filter Val:", Result_filter);
	if ( !splice_filter()->apply(*pose) ) {
		//  pose.dump_pdb("failed_filter_pose.pdb");
		TR << "Failing because filter fails" << std::endl;
		set_last_move_status(protocols::moves::FAIL_RETRY);
		retrieve_values();
		return 1;
	}
	return 0;
}

int
SpliceOut::SpliceOutRMSDFilter(core::pose::Pose * pose)  {
	bool pass_rmsd_filter = calculate_rmsd(*pose,*source_pose_,splicemanager.dofs().size(),source_from_res(), splicemanager.pose_from_res(),rms_cutoff(), rms_cutoff_loop()  );

	if ( !pass_rmsd_filter ) {
		set_last_move_status(protocols::moves::FAIL_RETRY);
		retrieve_values();
		// pose.dump_pdb("failed_rmsd_pose.pdb");
		return 1;
	}
	return 0;
}
void SpliceOut::build_ideal_segment(core::pose::Pose & pose){
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
	}
	if ( splicemanager.debug() ) {
		pose.dump_pdb(splicemanager.mover_name()+"_build_ideal_segment.pdb");
	}
}
void SpliceOut::minimize_segment(core::pose::Pose & pose){
	splicemanager.mm()->show();
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
std::string SpliceOut::name_for_filter(){
	return "SpliceOut";
}

} //splice
} //protocols
