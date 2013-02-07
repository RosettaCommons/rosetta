// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LoopCreationMover.cc
///
/// @brief
/// @author Tim Jacobs

// Unit Headers
#include <devel/loop_creation/LoopCreationMover.hh>
#include <devel/loop_creation/LoopCreationMoverCreator.hh>

//Basic
#include <basic/MetricValue.hh>
#include <basic/datacache/DiagnosticData.hh>

//Core
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

//Protocols
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/analysis/LoopAnalyzerMover.hh>
#include <protocols/loops/loops_main.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>


#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <protocols/filters/BasicFilters.hh>

#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>

//Devel
#include <devel/loop_creation/LoopInserter.hh>
#include <devel/loop_creation/LoopCloser.hh>
#include <devel/loop_creation/KICLoopCloser.hh>

//Utility
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/sort_predicates.hh>

//Basic
#include <basic/Tracer.hh>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif


namespace devel {
namespace loop_creation {

static basic::Tracer TR( "devel.loop_creation.LoopCreationMover" );

using namespace core;

//****CREATOR METHODS****//
std::string
LoopCreationMoverCreator::keyname() const
{
	return LoopCreationMoverCreator::mover_name();
}

protocols::moves::MoverOP
LoopCreationMoverCreator::create_mover() const {
	return new LoopCreationMover;
}

std::string
LoopCreationMoverCreator::mover_name()
{
	return "LoopCreationMover";
}
//****END CREATOR METHODS****//

///@brief default constructor
LoopCreationMover::LoopCreationMover():
	loop_inserter_(NULL),
	loop_closer_(NULL),
	loop_filter_(NULL),
	num_insertions_(5),
	closures_per_insertion_(10),
	refine_(true),
	design_loops_(false),
	include_neighbors_(false),
	minimize_loops_(true),
	loop_anchor_(0),
	dump_pdbs_(false)
{
	init();
}

///@brief explicit constructor
LoopCreationMover::LoopCreationMover(
	LoopInserterOP loop_inserter,
	LoopCloserOP loop_closer,
	core::Size num_insertions,
	core::Size closures_per_insertion,
	bool refine,
	bool design_loops,
	bool include_neighbors,
	bool minimize_loops
):
	Mover("LoopCreationMover"),
	loop_inserter_(loop_inserter),
	loop_closer_(loop_closer),
	num_insertions_(num_insertions),
	closures_per_insertion_(closures_per_insertion),
	refine_(refine),
	design_loops_(design_loops),
	include_neighbors_(include_neighbors),
	minimize_loops_(minimize_loops)
{
	init();
}

void
LoopCreationMover::init()
{
	loop_filter_=new protocols::filters::TrueFilter;
}

protocols::moves::MoverOP
LoopCreationMover::clone() const {
	return( protocols::moves::MoverOP( new LoopCreationMover( *this ) ) );
}
protocols::moves::MoverOP
LoopCreationMover::fresh_instance() const {
	return protocols::moves::MoverOP( new LoopCreationMover );
}
	
std::string
LoopCreationMover::get_name() const {
	return "LoopCreationMover";
}

protocols::loops::Loop
LoopCreationMover::get_last_created_loop() const
{
	return last_created_loop_;
}

//Handy structure for passing around data accrued during loop creation
struct loop_creation_data{
	core::pose::Pose pose;
	protocols::loops::Loop loop;
};
	
//Attempt to build residues to connect residues anchor point to anchor point+1.
//Do this num_iterations time and score each one that passes the loop filter.
//Set the best scoring to the current pose.
void
LoopCreationMover::apply(
	pose::Pose & pose
){
	using namespace std;
	using namespace core;
	using utility::vector1;
	
	//Sanity checks
	if(loop_closer_ == 0)
	{
		utility_exit_with_message("No LoopCloser given to LoopCreationMover!");
	}
	if(loop_inserter_ == 0)
	{
		utility_exit_with_message("No LoopInserter given to LoopCreationMover!");
	}
	if(loop_filter_ == 0)
	{
		loop_filter_ = new protocols::filters::TrueFilter;//default to true filter
	}
	TR.Debug << "Fold tree prior to LoopCreationMover: " << pose.fold_tree() << endl;

	protocols::jd2::JobOP const job_me ( protocols::jd2::JobDistributor::get_instance()->current_job() );
	std::string const job_name ( protocols::jd2::JobDistributor::get_instance()->job_outputter()->output_name(job_me) );

	//If a loop anchor was set, pass it off to the inserter
	if(loop_anchor_)
	{
		loop_inserter_->loop_anchor(loop_anchor_);
	}
	
	bool loop_closed = false;
	bool loop_passed = false;
	protocols::loops::Loop created_loop;
	
	vector1<pair<Real, loop_creation_data> > loop_scores;
	pose::Pose edit_pose;
	for(Size j=1; j<=num_insertions_; ++j)
	{
		edit_pose = pose;
		
		//Insert new loop residues
		TR.Debug << "Atempting insertion " << j << " of " << num_insertions_ << std::endl;
		clock_t start_time = clock();
		loop_inserter_->apply(edit_pose);
		clock_t insert_time = clock();
		TR.Debug << "Clocks - loop inserter run " << j << " of " << num_insertions_
			<< " - total time: " << insert_time - start_time << std::endl;
		
		created_loop = loop_inserter_->get_created_loop();
		TR.Debug << "New loop: " << created_loop << endl;

		//Save the inserted-only pose for input into loop closure
		pose::Pose inserted_pose = edit_pose;

		//DEBUG
		if(dump_pdbs_){
			std::stringstream filename;
			filename << job_name << "_insertion_" << j << ".pdb";
			inserted_pose.dump_pdb(filename.str());
		}
		//END DEBUG

		for(Size k=1; k<=closures_per_insertion_; ++k)
		{
			edit_pose = inserted_pose;
			
			start_time = clock();
			loop_closer_->loop(created_loop);
			loop_closer_->apply(edit_pose);
			clock_t closure_time = clock();
			TR.Debug << "Clocks - loop closer run " << k << " of " << closures_per_insertion_
				<< " - total time: " << closure_time - start_time << std::endl;
			if(loop_closer_->success())
			{
				TR.Debug << "Loop insertion " << j << ", closer iteration " << k << " succeeded" << std::endl;

				//ensure that a chemical bond exists now that we've closed the loop
				edit_pose.conformation().declare_chemical_bond(created_loop.cut(), "C", created_loop.cut()+1,"N");

				//For some reason the creation/closing of loops really f***s the hydrogen placement.
				//This is an attempt to fix that
				for(core::Size i=created_loop.start()-1; i<=created_loop.stop()+1; ++i)
				{
					core::conformation::ResidueOP new_res = edit_pose.residue(i).clone();
					core::conformation::idealize_hydrogens(*new_res, edit_pose.conformation());
					edit_pose.conformation().replace_residue(i, *new_res, false);
				}

				//DEBUG
				if(dump_pdbs_){
					std::stringstream filename;
					filename << job_name << "_insertion_" << j << "_closure_" << k << ".pdb";
					edit_pose.dump_pdb(filename.str());
				}
				//END DEBUG

				loop_closed=true;
				Real score=0.0;
				if(refine_){
					start_time = clock();
					score = refine_loop(edit_pose, created_loop);
					clock_t refine_time = clock();
					TR.Debug << "Clocks - loop refinement completed in total time: " << refine_time - start_time << std::endl;

					//DEBUG
					if(dump_pdbs_){
						std::stringstream filename;
						filename << job_name << "_insertion_" << j << "_closure_" << k << "_refined.pdb";
						edit_pose.dump_pdb(filename.str());
					}
					//END DEBUG
				}

				loop_creation_data lcd;
				lcd.pose = edit_pose;
				lcd.loop = created_loop;
				if(loop_filter_->apply(edit_pose))
				{
					loop_passed=true;
					loop_scores.push_back(make_pair(score, lcd));
			
					TR.Debug << "Loop closure iteration " << k << ", of insertion " << j
						<< ", completed succesfully with score: " << score << std::endl;
				}
			}
		}
	}
	if(!loop_closed)
	{
		stringstream err;
		TR << "Unable to create and close a loop between residues "
			<< loop_inserter_->loop_anchor() << " and " << loop_inserter_->loop_anchor()+1
			<< ", consider longer loop lengths or a different LoopInserter and/or LoopCloser." << endl;
		utility_exit_with_message(err.str());
//		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
//		return;
	}
	if(!loop_passed)
	{
		stringstream err;
		TR << "No closed loops between residues "
			<< loop_inserter_->loop_anchor() << " and " << loop_inserter_->loop_anchor()+1
			<< " passed filters. Try relaxing filters." << endl;
		utility_exit_with_message(err.str());
//		set_last_move_status(protocols::moves::FAIL_DO_NOT_RETRY);
//		return;
	}
	
	sort(loop_scores.begin(), loop_scores.end(), utility::SortFirst<Real, loop_creation_data>());
	loop_creation_data const & best_lcd = loop_scores[1].second;
	TR << "Best loop from has a score of: " << loop_scores[1].first << std::endl;
	
	last_created_loop_=best_lcd.loop;
	pose=best_lcd.pose;
	
	protocols::loops::Loops loops;
	loops.add_loop(last_created_loop_);
	protocols::analysis::LoopAnalyzerMoverOP la = new protocols::analysis::LoopAnalyzerMover(loops, true);
	la->apply(pose);
}

void
LoopCreationMover::copy_last_loop_to_new_anchor(
	core::pose::Pose & pose,
	core::Size new_anchor
){
	using namespace core;
	using namespace core::chemical;
	
	TR.Debug << "removing termini types" << std::endl;
	pose::remove_upper_terminus_type_from_pose_residue(pose, new_anchor);
	pose::remove_lower_terminus_type_from_pose_residue(pose, new_anchor+1);

	core::Size old_anchor = loop_inserter_->loop_anchor();
	core::Size old_cut = last_created_loop_.cut();
	
	TR.Debug << "Loop to copy: " << last_created_loop_ << std::endl;
	
	core::Size modifications_begin = loop_inserter_->modified_range().first;
	core::Size mod_begin_offset = old_anchor-modifications_begin;
	core::Size new_modifications_begin=new_anchor-mod_begin_offset;
	
	core::Size modifications_end = loop_inserter_->modified_range().second;
	core::Size mod_end_offset = old_anchor-modifications_end;
	core::Size new_modifications_end=new_anchor-mod_end_offset;
	new_modifications_end-=last_created_loop_.size();
	TR.Debug << "Old modifications range: " << modifications_begin << " " << modifications_end << std::endl;
	
	TR.Debug << "new mod begin " << new_modifications_begin << std::endl;
	TR.Debug << "new mod end " << new_modifications_end << std::endl;
	
//	core::Size new_loop_start_offset = old_anchor-last_created_loop_.start();
//	core::Size new_loop_start=new_anchor-new_loop_start_offset;
//	TR.Debug << "new loop start " << new_loop_start << std::endl;
	
	protocols::loops::Loop extended_loop(new_modifications_begin, new_modifications_end, new_anchor);
	protocols::loops::set_single_loop_fold_tree(pose, extended_loop);
	
	//copy the loop residues
	core::Size anchor_pos = new_anchor;
	for(core::Size i=last_created_loop_.start(); i<=last_created_loop_.cut(); ++i)
	{
		TR.Debug << "Append old resnum " << i << "(" << oneletter_code_from_aa(pose.aa(i))
			<< ") to new anchor " << anchor_pos << "(" <<oneletter_code_from_aa(pose.aa(anchor_pos)) << ")" << std::endl;
		pose.append_polymer_residue_after_seqpos(pose.residue(i), anchor_pos, true);
		++anchor_pos;
		TR.Debug << "Fold tree after attach: " << pose.fold_tree() << std::endl;
	}
	
	++anchor_pos;
	for(core::Size i=last_created_loop_.stop(); i>=last_created_loop_.cut()+1; --i)
	{
		TR.Debug << "Prepending old resnum " << i << "(" << oneletter_code_from_aa(pose.aa(i))
			<< ") to new anchor " << anchor_pos << "(" <<oneletter_code_from_aa(pose.aa(anchor_pos)) << ")" << std::endl;
		pose.prepend_polymer_residue_before_seqpos(pose.residue(i), anchor_pos, true);
		TR.Debug << "Fold tree after prepend: " << pose.fold_tree() << std::endl;
	}
	
	for(core::Size resnum=modifications_begin; resnum<=old_cut; ++resnum)
	{
		core::Size offset = old_anchor-resnum;
		TR.Debug << "Copying torsions from " << resnum << " to " << new_anchor-offset << std::endl;
		pose.set_phi( new_anchor-offset, pose.phi(resnum) );
		pose.set_psi( new_anchor-offset, pose.psi(resnum) );
		pose.set_omega( new_anchor-offset, pose.omega(resnum) );
	}
	for(core::Size resnum=modifications_end; resnum>=old_cut+1; --resnum)
	{
		core::Size offset = old_anchor-resnum;
		TR.Debug << "Copying torsions from " << resnum << " to " << new_anchor-offset << std::endl;
		pose.set_phi( new_anchor-offset, pose.phi(resnum) );
		pose.set_psi( new_anchor-offset, pose.psi(resnum) );
		pose.set_omega( new_anchor-offset, pose.omega(resnum) );
	}
//	pose.conformation().declare_chemical_bond(created_loop.cut(), "C", created_loop.cut()+1,"N");
}

core::Real
LoopCreationMover::refine_loop(
	core::pose::Pose & pose,
	protocols::loops::Loop loop
){
	using namespace std;
	using namespace core;
	using utility::vector1;
	
	scoring::ScoreFunctionOP scorefxn = scoring::getScoreFunction();
	scorefxn->set_weight( scoring::chainbreak, 100.0 );//loop should already be closed, so this will just prevent re-breaking by minimization
	
	
	//Initialize a mover if we haven't already
	protocols::simple_moves::PackRotamersMoverOP pack_mover =
		new protocols::simple_moves::PackRotamersMover;
	
	
	//Setup task factory for minimization and packing
	pack::task::TaskFactoryOP task_factory = new pack::task::TaskFactory;
	
	set<Size> loop_residues;
	for(Size i=loop.start(); i<=loop.stop(); ++i)
	{
		loop_residues.insert(i);
	}
	
	//If we're including loop neighbors in packing/redesign then calculate them
	set<Size> residues_to_pack;
	if(include_neighbors_)
	{
		string const nb_calc("neighbor_calculator");
		pose::metrics::CalculatorFactory::Instance().register_calculator( nb_calc,
			new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( loop_residues ) );
			
		basic::MetricValue< set< Size > > neighbor_mv;
		pose.metric( nb_calc, "neighbors", neighbor_mv);
		residues_to_pack = neighbor_mv.value();
		pose::metrics::CalculatorFactory::Instance().remove_calculator(nb_calc);
	}
	else
	{
		residues_to_pack=loop_residues;
	}

	//Restrict all non-packable residues and allow movements at all packable residues
	pack::task::operation::PreventRepackingOP prevent_packing =
		new pack::task::operation::PreventRepacking;
	pack::task::operation::RestrictResidueToRepackingOP repack_res =
		new pack::task::operation::RestrictResidueToRepacking();
	for(Size i=1; i<=pose.total_residue(); ++i)
	{
		if(residues_to_pack.find(i)==residues_to_pack.end())
		{
			prevent_packing->include_residue(i);
		}
		else
		{
			//If we're not designing then restrict packable residues to repack-only
			if(!design_loops_)
			{
				repack_res->include_residue(i);
			}
		}
	}
	task_factory->push_back(prevent_packing);
	task_factory->push_back(repack_res);
	
	pack_mover->task_factory(task_factory);
	pack_mover->apply(pose);
	
	if(minimize_loops_)
	{
		core::Size modifications_begin = loop_inserter_->modified_range().first;
		core::Size modifications_end = loop_inserter_->modified_range().second;
		protocols::loops::Loop extended_loop(modifications_begin, modifications_end, loop.cut());
		protocols::loops::set_single_loop_fold_tree(pose, extended_loop);
		
		kinematics::MoveMapOP movemap = new kinematics::MoveMap;
		for(core::Size i=loop.start(); i<=loop.stop(); ++i)
		{
			movemap->set_bb(i, true);
			movemap->set_chi(i, true);
		}

//		movemap->set_bb(true);
//		movemap->set_chi(true);
			
		TR.Debug << "Movemap for minimization: " << *movemap << std::endl;
			
		protocols::simple_moves::MinMoverOP min_mover =
			new protocols::simple_moves::MinMover(movemap, scorefxn, "dfpmin_armijo_nonmonotone", 0.01, false );
			
		TR << "Score prior to minimization: " << scorefxn->score(pose) << std::endl;
//		pose.dump_pdb("pre_minimization.pdb");
		
		min_mover->apply(pose);
	}
	
//	utility::vector1< bool > residues_to_score(pose.total_residue(), false);
//	for(core::Size i=loop.start(); i<=loop.stop(); ++i)
//	{
//		residues_to_score[i]=true;
//	}
//	
//	core::Real loop_score = scorefxn->get_sub_score(pose, residues_to_score);
	core::Real loop_score = scorefxn->score(pose);
	return (loop_score/loop.size());
}
	
LoopInserterOP
LoopCreationMover::loop_inserter() const
{
	return loop_inserter_;
}

LoopCloserOP
LoopCreationMover::loop_closer() const
{
	return loop_closer_;
}

void
LoopCreationMover::parse_my_tag(
	TagPtr const tag,
	protocols::moves::DataMap & data,
	Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	Pose const & pose )
{
	using namespace std;
	
	//****REQUIRED TAGS****//
	if(tag->hasOption("loop_inserter")){
		string const loop_inserter_name( tag->getOption< string >( "loop_inserter" ) );
		protocols::moves::Movers_map::const_iterator find_mover( movers.find( loop_inserter_name ) );
		bool const mover_found( find_mover != movers.end() );
		if( mover_found )
			loop_inserter_ = dynamic_cast< LoopInserter* > (find_mover->second());
		else
			utility_exit_with_message( "Mover " + loop_inserter_name + " not found" );
	}
	else{
		utility_exit_with_message("loop_inserter tag is required");
	}
	
	if(tag->hasOption("loop_closer")){
		string const loop_closer_name( tag->getOption< string >( "loop_closer" ) );
		protocols::moves::Movers_map::const_iterator find_mover( movers.find( loop_closer_name ) );
		bool const mover_found( find_mover != movers.end() );
		if( mover_found )
			loop_closer_ = dynamic_cast< LoopCloser* > (find_mover->second());
		else
			utility_exit_with_message( "Mover " + loop_closer_name + " not found" );
	}
	else{
		utility_exit_with_message("loop_closer tag is required");
	}
	
	//****OPTIONAL TAGS****//

	if(tag->hasOption("dump_pdbs")){
		dump_pdbs_ = tag->getOption<bool>("dump_pdbs");
	}

	if(tag->hasOption("refine")){
		refine_ = tag->getOption<bool>("refine");
	}
	
	if(tag->hasOption("num_insertions")){
		num_insertions_ = tag->getOption< core::Size >("num_insertions");
	}
	
	if(tag->hasOption("closures_per_insertion")){
		closures_per_insertion_ = tag->getOption< core::Size >("closures_per_insertion");
	}
	
	if(tag->hasOption("minimize_loops")){
		minimize_loops_ = tag->getOption<bool>("minimize_loops");
	}
	
	if(tag->hasOption("design_loops")){
		design_loops_ = tag->getOption<bool>("design_loops");
	}
	
	if(tag->hasOption("include_neighbors")){
		include_neighbors_ = tag->getOption<bool>("include_neighbors");
	}
	
	//Useful for setting loop anchor from loop creation mover
	if(tag->hasOption("loop_anchor")){
		loop_anchor_ = tag->getOption< core::Size >("loop_anchor");
	}
}

} //loop creation
} //devel
