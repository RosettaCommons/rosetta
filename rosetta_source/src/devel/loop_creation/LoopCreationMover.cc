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
#include <devel/loop_creation/LoopInserter.hh>
#include <devel/loop_creation/LoopCloser.hh>
#include <devel/loop_creation/KICLoopCloser.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <protocols/filters/BasicFilters.hh>

#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>

//Utility
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/vector0.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/sort_predicates.hh>

//Basic
#include <basic/Tracer.hh>

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
loop_inserter_(0),
loop_closer_(0),
loop_filter_(0),
num_iterations_(10),
design_loops_(false),
include_neighbors_(false),
minimize_loops_(true)
{
	init();
}

///@brief explicit constructor
LoopCreationMover::LoopCreationMover(
	LoopInserterOP loop_inserter,
	LoopCloserOP loop_closer,
	core::Size num_iterations,
	bool design_loops,
	bool include_neighbors,
	bool minimize_loops
):
Mover("LoopCreationMover"),
loop_inserter_(loop_inserter),
loop_closer_(loop_closer),
num_iterations_(num_iterations),
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
	
//Attempt to build residues to connect residues anchor point to anchor point+1
void
LoopCreationMover::apply(
	pose::Pose & pose
){
	using namespace std;
	using namespace core;
	using utility::vector1;
	
	calc_counter_=0;

	//Sanity checks
	if(loop_closer_ == 0)
	{
		utility_exit_with_message("No LoopCloser given to LoopCreationMover!");
	}
	if(loop_inserter_ == 0)
	{
		utility_exit_with_message("No LoopInserter given to LoopCreationMover!");
	}
	TR.Debug << "Fold tree prior to LoopCreationMover: " << pose.fold_tree() << endl;
	
	bool loop_closed = false;
	bool loop_passed = false;
	protocols::loops::Loop created_loop;
	
	vector1<pair<Real, loop_creation_data> > loop_scores;
	pose::PoseOP edit_pose;
	for(Size j=1; j<=num_iterations_; ++j)
	{
		pose::Pose edit_pose = pose;
		
		//Insert new loop residues
		loop_inserter_->apply(edit_pose);
		created_loop = loop_inserter_->get_created_loop();
		
		TR.Debug << "New loop: " << created_loop << endl;
		TR.Debug << "New ft: " << edit_pose.fold_tree() << endl;
		
		//DEBUGGING
//		edit_pose.dump_pdb("after_insert.pdb");
		
		loop_closer_->loop(created_loop);
		loop_closer_->apply(edit_pose);
		if(loop_closer_->success())
		{
			loop_closed=true;
			if(loop_filter_->apply(edit_pose))
			{
				loop_passed=true;
				Real score = refine_loop(edit_pose, created_loop);
				
				loop_creation_data lcd;
				lcd.pose = edit_pose;
				lcd.loop = created_loop;
				
				loop_scores.push_back(make_pair(score, lcd));
			}
		}
	}
	if(!loop_closed)
	{
		stringstream err;
		err << "Unable to create and close a loop between residues "
			<< loop_inserter_->loop_anchor() << " and " << loop_inserter_->loop_anchor()+1
			<< ", consider longer loop lengths or a different LoopInserter and/or LoopCloser." << endl;
		utility_exit_with_message(err.str());
	}
	if(!loop_passed)
	{
		stringstream err;
		err << "No closed loops between residues "
			<< loop_inserter_->loop_anchor() << " and " << loop_inserter_->loop_anchor()+1
			<< " passed filters. Try relaxing filters." << endl;
		utility_exit_with_message(err.str());
	}
	
	sort(loop_scores.begin(), loop_scores.end(), utility::SortFirst<Real, loop_creation_data>());
	loop_creation_data const & best_lcd = loop_scores[1].second;
	
	last_created_loop_=best_lcd.loop;
	pose=best_lcd.pose;
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
	scorefxn->set_weight( scoring::chainbreak, 1.0 );
	
	//Initialize a mover if we haven't already
	if(pack_mover_==0){
		pack_mover_ = new protocols::simple_moves::PackRotamersMover;
	}
	
	kinematics::MoveMapOP movemap = new kinematics::MoveMap;
	
	//Setup task factory for minimization and packing
	pack::task::TaskFactoryOP task_factory = new pack::task::TaskFactory;
	
	set<Size> loop_residues;
	for(Size i=loop.start(); i<=loop.stop(); ++i)
	{
		loop_residues.insert(i);
	}
	
	set<Size> residues_to_pack;
	if(include_neighbors_)
	{
		++calc_counter_;
		string const nb_calc("neighbor_calculator_"+utility::to_string(calc_counter_));
		pose::metrics::CalculatorFactory::Instance().register_calculator( nb_calc,
			new protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator( loop_residues ) );
			
		basic::MetricValue< set< Size > > neighbor_mv;
		pose.metric( nb_calc, "neighbors", neighbor_mv);
		residues_to_pack = neighbor_mv.value();
		
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
			movemap->set_bb(i, true);
			movemap->set_chi(i, true);
			//If we're not designing then restrict packable residues to repack-only
			if(!design_loops_)
			{
				repack_res->include_residue(i);
			}
		}
	}
	task_factory->push_back(prevent_packing);
	task_factory->push_back(repack_res);
	
	pack_mover_->task_factory(task_factory);
	pack_mover_->apply(pose);
	
	if(minimize_loops_)
	{
		kinematics::FoldTree ft;
		
		ft.add_edge(1, loop.start()-1, kinematics::Edge::PEPTIDE);
		
		ft.add_edge(1, loop.start(), 1);
		ft.add_edge(loop.start(), loop.cut(), kinematics::Edge::PEPTIDE);
		
		ft.add_edge(1, loop.stop(), 2);
		ft.add_edge(loop.stop(), loop.cut()+1, kinematics::Edge::PEPTIDE);
		
		ft.add_edge(1, loop.stop()+1, 3);
		ft.add_edge(loop.stop()+1, pose.total_residue(), kinematics::Edge::PEPTIDE);
		pose.fold_tree(ft);
	
			
		pose::add_variant_type_to_pose_residue(
			pose, chemical::CUTPOINT_LOWER, loop.cut());
		pose::add_variant_type_to_pose_residue(
			pose, chemical::CUTPOINT_UPPER, loop.cut()+1);
			
		protocols::simple_moves::MinMoverOP base_min_mover =
			new protocols::simple_moves::MinMover(movemap, scorefxn, "dfpmin_armijo_nonmonotone", 0.01, false );
			
		pose.dump_pdb("pre_minimization.pdb");
		
		base_min_mover->apply(pose);
	}
	return scorefxn->score(pose);
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
	if(tag->hasOption("num_iterations")){
		num_iterations_ = tag->getOption< core::Size >("num_iterations");
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
}

} //loop creation
} //devel
