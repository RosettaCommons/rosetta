// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/remodel/RemodelLoopMover.cc
/// @brief  Loop modeling protocol based on routines from Remodel and EpiGraft
///         packages in Rosetta++.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)
/// @author Possu Huang (possu@u.washington.edu)


// numeric headers
#include <numeric/random/random.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>

#include <core/types.hh>
#include <core/conformation/util.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/ResfileReader.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/util/disulfide_util.hh>

// unit headers
#include <protocols/forge/remodel/RemodelDesignMover.hh>
#include <protocols/forge/remodel/RemodelRotamerLinks.hh>
#include <protocols/forge/methods/util.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/simple_moves/MutateResidue.hh>

// C++ headers

using namespace basic::options;

namespace protocols {
namespace forge{
namespace remodel{

// Tracer instance for this file
// Named after the original location of this code
static thread_local basic::Tracer TR( "protocols.forge.remodel.RemodelDesignMover" );

// RNG


// @brief default constructor
RemodelDesignMover::RemodelDesignMover()
{
	// has to reinitialize state before apply
	state_.clear();
}

/// @brief value constructor
RemodelDesignMover::RemodelDesignMover( RemodelData const & remodel_data,
																				RemodelWorkingSet const & working_model,
																				ScoreFunctionOP const & sfxn )
{

	using core::pose::metrics::CalculatorFactory;
	using core::pose::metrics::PoseMetricCalculatorOP;
	using protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator;

	remodel_data_ = remodel_data;
	working_model_ = working_model;
	archived_starting_task_ = working_model.task;
	score_fxn_ = sfxn->clone();

	// setup calculators
	CalculatorFactory::Instance().remove_calculator( "neighborhood_calc" );

	std::set< Size > und_pos;

	std::set< core::Size > uup = working_model.manager.union_of_intervals_containing_undefined_positions();
	//for ( std::set<core::Size>::iterator i = uup.begin(); i!=uup.end(); i++){
	//	TR << *i <<  " UUP in DesignMover" <<  std::endl;
	//}

	if (option[OptionKeys::remodel::repeat_structure].user() ) {
		Size repeatCount =option[OptionKeys::remodel::repeat_structure];
		for ( Size rep = 0; rep < repeatCount ; rep++ ) {
			for ( std::set< Size >::iterator it = uup.begin(); it != uup.end(); ++it ) {
			//DEBUG
				//std::cout << *it + remodel_data.blueprint.size()*rep << std::endl;
				//std::cout << "manger size"  << working_model.manager.union_of_intervals_containing_undefined_positions().size() <<  std::endl;
				//std::cout << *it  << std::endl;
				if ( !(*it+remodel_data.blueprint.size()*rep > (remodel_data.blueprint.size()*repeatCount)) ){ //Extrapolation of positions shouldn't go beyond the length of pose
					und_pos.insert(*it + remodel_data.blueprint.size()*rep);
				}
			}
		}

	} else {
		und_pos = working_model.manager.union_of_intervals_containing_undefined_positions();
	}

	TR << "Creating NeighborhoodByDistanceCalculator using und_pos: [ ";
	for ( std::set< Size >::iterator itr = und_pos.begin(); itr != und_pos.end(); itr++ ) {
		TR << *itr << " ";
	}
	TR << "]" << std::endl;

	if ( und_pos.empty() ) {
		TR << "Warning: union_of_intervals_containing_undefined_positions() returned empty set. NeighborhoodByDistanceCalculator could return undefined results." << std::endl;
	}

	CalculatorFactory::Instance().register_calculator( "neighborhood_calc", PoseMetricCalculatorOP( new NeighborhoodByDistanceCalculator( und_pos ) ) );

}

/// @brief default destructor
RemodelDesignMover::~RemodelDesignMover(){}

/// @brief clone this object
RemodelDesignMover::MoverOP RemodelDesignMover::clone() const {
  return RemodelDesignMover::MoverOP( new RemodelDesignMover( *this ) );
}

/// @brief create this type of object
RemodelDesignMover::MoverOP RemodelDesignMover::fresh_instance() const {
  return RemodelDesignMover::MoverOP( new RemodelDesignMover() );
}

/// @brief packer task accessor
core::pack::task::PackerTaskOP & RemodelDesignMover::task(){
	return working_model_.task;
}

/// @brief score function setter
void RemodelDesignMover::scorefunction( ScoreFunctionOP const & sfxn ) {
	score_fxn_ = sfxn->clone();
}

std::string RemodelDesignMover::get_name() const {
	return "RemodelDesignMover";
}

bool RemodelDesignMover::check_state() {
	if (state_.empty()) {
		return false;
		TR << "state tag not set " << std::endl;
	} else {
		TR << "Design Mover state: " << state_ << std::endl;
		return true;
	}
}

void RemodelDesignMover::set_state( std::string state_tag ){
	state_ = state_tag;
	// reset the task
	working_model_.task = archived_starting_task_;
}



///
/// @begin RemodelDesignMover::apply
///
/// @brief
/// Apply method for Mover.
/// Checks value of option -remodel::design::no_design
/// -remodel::design::find_neighbors
/// -remodel::design::design_neigbors
/// -remodel::design::skip_partial
///
void RemodelDesignMover::apply( Pose & pose ) {

	if (option[OptionKeys::remodel::design::no_design].user() ) {
		TR << "bypassing design due to invokation of -no_design" << std::endl;
		return;
	}

	// make decision as to which mode to apply
	bool manual = remodel_data_.has_design_info_;
	bool neighbor =option[OptionKeys::remodel::design::find_neighbors].user();
	bool design =option[OptionKeys::remodel::design::design_neighbors].user();

	if ( !check_state() ) {
		basic::Error() << "check_state failed, has to set_state first " << std::endl;
	}

	// based on the values of the variables manual, neighbor and design, figure out which "mode" of design we're doing
	// and call the appropriate packertask function
	if ( manual ){
		if (neighbor){
			if (design){
				mode5_packertask(pose);
			} else {
				mode6_packertask(pose);
			}
		} else {
			mode4_packertask(pose);
		}
	} else {
		if (neighbor){
			if (design){
				mode2_packertask(pose);
			} else {
				mode3_packertask(pose);
			}
		} else {
			if(option[OptionKeys::remodel::design::design_all].user()){
				mode1_1_packertask(pose);
			}
			else {
				mode1_packertask(pose);
			}
		}
	}

	if (option[OptionKeys::remodel::repeat_structure].user() ) {
		//try turning off bumpcheck
		//TR << "bumpcheck off" << std::endl;
		working_model_.task->set_bump_check( true );

		// make rotamer links
		RemodelRotamerLinksOP linkOP( new RemodelRotamerLinks );
		linkOP->apply( pose, *working_model_.task );
	}

	if (!strcmp(state_.c_str(), "stage")){
		if (manual){
			//do nothing
		}
		else {  //auto build always reduce task

			if (!option[OptionKeys::remodel::design::skip_partial].user()){
				reduce_task(pose, working_model_.task, true, true, false);
			}
			else {
				reduce_task(pose, working_model_.task, true, true, true);
			}
		}
	}
	else if (!strcmp(state_.c_str(), "finish")){
		if (manual){
			//do nothing
		}
		else {
		// if finishing design, no need to reduce, but require resetting positios
    //if (!option[OptionKeys::remodel::design::skip_partial].user()){
			reduce_task(pose, working_model_.task, true, true, true);
		}
	}

	//debug
	//TR.Debug << working_model_.task->task_string(pose) << std::endl;
	//TR.Debug << *working_model_.task << std::endl;

	core::pack::pack_rotamers(pose, *score_fxn_ , working_model_.task);
	score_fxn_->show(TR, pose);
	TR << std::endl;

// pose.dump_pdb("junkCheck.pdb");
//	core::pack::pack_rotamers(pose, *scorefxn , working_model_.task);
//	core::pack::pack_rotamers(pose, *scorefxn , working_model_.task);
// core::pack::task::TaskFactoryOP TF = new core::pack::task::TaskFactory;
//TF->create_packer_task(pose);
//	core::pack::pack_rotamers(pose, *scorefxn , TF->create_packer_task(pose));
//	core::pack::pack_rotamers(pose, *scorefxn , TF->create_packer_task(pose));
//	core::pack::pack_rotamers(pose, *scorefxn , TF->create_packer_task(pose));
}


void RemodelDesignMover::run_calculator( core::pose::Pose const & pose, std::string const & calculator, std::string const & calculation, utility::vector1_bool & residues ) {

	runtime_assert( residues.size() == pose.total_residue() );

	// find the set of residues
	typedef std::set< core::Size > SizeSet;
	basic::MetricValue< SizeSet > mv_sizeset;
	pose.metric( calculator, calculation, mv_sizeset );
	SizeSet const & sizeset( mv_sizeset.value() );
	//TR << "runCalculator " << std::endl;

	// insert this into the vector
	for( SizeSet::const_iterator it(sizeset.begin()), end(sizeset.end()) ; it != end; ++it ) {
		//TR << *it <<  " debug run_calc " << std::endl;
		residues[*it] = true;
	}

	return;
}

void RemodelDesignMover::reduce_task( Pose & pose, core::pack::task::PackerTaskOP &task, bool core, bool boundary, bool surface){

  // setup calculators
	using core::pose::metrics::CalculatorFactory;
	using core::pose::metrics::PoseMetricCalculatorOP;
  using protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator;
	using core::pack::task::ResfileCommandOP;

	//std::cout << "START REDUCE" << std::endl;
	//std::cout << *task << std::endl;

  // need to collect all positions before reduction
  CalculatorFactory::Instance().remove_calculator( "reducetask_calc" );
	utility::vector1_bool boollist(pose.total_residue());
  run_calculator(pose, "neighborhood_calc", "neighbors", boollist);
	std::set<Size> positionList;
	for (Size i = 1; i <= pose.total_residue(); i++){
		if (boollist[i] && !option[OptionKeys::remodel::design::design_all].user()){
			positionList.insert(i);
		}
		else { // in case of design all flag, take all positions
			positionList.insert(i);
		}
	}

  CalculatorFactory::Instance().register_calculator(
    "reducetask_calc",
    PoseMetricCalculatorOP( new NeighborhoodByDistanceCalculator( positionList ) )
  );

	//compute the bsasa values for each position
	utility::vector1< core::Real > sasa_list;
	sasa_list = protocols::forge::methods::calc_rsd_sasa( pose );

	utility::vector1_bool resclass(pose.total_residue(),false);
	utility::vector1_bool neighbor_count(pose.total_residue(),false);
	utility::vector1<Size> corePos;
	utility::vector1<Size> boundaryPos;
	utility::vector1<Size> surfPos;

	Size CORE_CUTOFF =option[OptionKeys::remodel::core_cutoff];
	Size BOUNDARY_CUTOFF =option[OptionKeys::remodel::boundary_cutoff];

	//get num_neighbors for each position
	basic::MetricValue< std::map< core::Size, core::Size > >nbr_map;
	pose.metric("reducetask_calc", "num_neighbors_map", nbr_map);
  std::map< core::Size, core::Size > sizemap;

	if(option[OptionKeys::remodel::resclass_by_sasa].user()){
		//simply repackage the values so either metric can feed into the code
		int count = 1 ;
		for ( utility::vector1< core::Real >::iterator it = sasa_list.begin(), ite = sasa_list.end(); it != ite ; it++){
				sizemap[count] = Size(*it);
				count++;
		}
	} else {
		sizemap = nbr_map.value();
	}

		TR.Debug << "sizemap content " << sizemap.size() << std::endl;

	for (std::map< core::Size, core::Size>::iterator it= sizemap.begin(); it!=sizemap.end(); it++){
		TR.Debug << "neighborlist " << (*it).first << " " <<  (*it).second << std::endl;
	}

	if(option[OptionKeys::remodel::repeat_structure].user()){

		utility::vector1<bool> visited(pose.total_residue(),false);

		for (Size i = 1; i<= resclass.size(); i++){ //check everyposition in the packertask

			if (visited[i]){
				continue;
			}
			// process linkage info
			utility::vector1<int> copies = task->rotamer_links()->get_equiv( i );

			int coreCount = 0;
			int boundaryCount = 0;
			int surfCount = 0;

			for (unsigned jj = 1; jj <= copies.size(); ++jj){

					visited[ copies[jj] ] = true;

					TR.Debug << "sizemap in repeat decision " << copies[jj] << " " << sizemap[ copies[jj] ] << std::endl;

				if(option[OptionKeys::remodel::resclass_by_sasa].user()){
					//take the counts for each set
					if ( sizemap[ copies[jj] ] <= CORE_CUTOFF){
						coreCount++;
					} else if ( sizemap[ copies[jj] ] > CORE_CUTOFF && sizemap[ copies[jj] ] <= BOUNDARY_CUTOFF ){
						boundaryCount++;
					} else if ( sizemap[ copies[jj] ] > BOUNDARY_CUTOFF){
						surfCount++;
					} else {
						TR << "RESCLASS ERROR" << sizemap[i] << std::endl;
					}
			}else{
			//take the counts for each set
					if ( sizemap[ copies[jj] ] >= CORE_CUTOFF){
						coreCount++;
					} else if ( sizemap[ copies[jj] ] < CORE_CUTOFF && sizemap[ copies[jj] ] >= BOUNDARY_CUTOFF ){
						boundaryCount++;
					} else if ( sizemap[ copies[jj] ] < BOUNDARY_CUTOFF){
						surfCount++;
					} else {
						TR << "RESCLASS ERROR" << sizemap[i] << std::endl;
					}
			}

			}

			//assign
			if (coreCount > 0) { //if any of them is core, turn everythign to core
				//std::cout << "core: ";
				for (unsigned jj = 1; jj <= copies.size(); ++jj){
					corePos.push_back(copies[jj]);
				//std::cout << copies[jj] << " " ;
				}
				//std::cout << std::endl;
			}
			else if (coreCount == 0 && boundaryCount > 0){
				//std::cout << "boundary: ";
				for (unsigned jj = 1; jj <= copies.size(); ++jj){
					boundaryPos.push_back(copies[jj]);
				//std::cout << copies[jj] << " " ;
				}
				//std::cout << std::endl;
			}
			else if (coreCount == 0 && boundaryCount == 0){
				//std::cout << "surf: ";
				for (unsigned jj = 1; jj <= copies.size(); ++jj){
					surfPos.push_back(copies[jj]);
				//std::cout << copies[jj] << " " ;
				}
				//std::cout << std::endl;
			}
			else {
				TR << "no idea what kind of scenario this would be" << std::endl;
			}

		}
	} //if repeat
	else {
		for (Size i = 1; i<= resclass.size(); i++){ //check everyposition in the packertask
			if (task->nonconst_residue_task(i).being_packed() && sizemap[i]){
				TR << "touch position " << i << std::endl;
				if(option[OptionKeys::remodel::resclass_by_sasa].user()){
								if (sizemap[i] <= CORE_CUTOFF){
									corePos.push_back(i);
								} else if ( sizemap[i] > CORE_CUTOFF && sizemap[i] <= BOUNDARY_CUTOFF){
									boundaryPos.push_back(i);
								} else if ( sizemap[i] > BOUNDARY_CUTOFF){
									surfPos.push_back(i);
								} else {
									TR << "RESCLASS ERROR" << sizemap[i] << std::endl;
								}
				}else{
								if (sizemap[i] >= CORE_CUTOFF){
									corePos.push_back(i);
								} else if ( sizemap[i] < CORE_CUTOFF && sizemap[i] >= BOUNDARY_CUTOFF){
									boundaryPos.push_back(i);
								} else if ( sizemap[i] < BOUNDARY_CUTOFF){
									surfPos.push_back(i);
								} else {
									TR << "RESCLASS ERROR" << sizemap[i] << std::endl;
								}
				}
			}
		}
	}

	//std::cout << "MID REDUCE" << std::endl;
	//std::cout << *task << std::endl;

/*
		//debug:
		for (utility::vector1<Size>::iterator it= corePos.begin(), end=corePos.end(); it!=end; it++){
			TR.Debug << "DEBUG: core positions:" << *it << std::endl;
		}
		for (utility::vector1<Size>::iterator it= boundaryPos.begin(), end=boundaryPos.end(); it!=end; it++){
			TR.Debug << "DEBUG: boundary positions:" << *it << std::endl;
		}
		for (utility::vector1<Size>::iterator it= surfPos.begin(), end=surfPos.end(); it!=end; it++){
			TR.Debug << "DEBUG: surface positions:" << *it << std::endl;
		}
*/
		//build new reduced task
//		core::pack::task::TaskFactoryOP TF = protocols::forge::methods::remodel_generic_taskfactory();
//		task=TF->create_task_and_apply_taskoperations( pose );

		//borrow RESFILE command for the job -- save a lot of coding effort.
		if ( core ){
			for (utility::vector1<Size>::iterator it= corePos.begin(), end=corePos.end(); it!=end; it++){
				ResfileCommandOP command( new core::pack::task::APOLAR );
				utility::vector1<std::string> decoy;
				decoy.push_back("APOLAR");
				Size resid = *it;
				TR << "APOLAR " << resid << std::endl;
				Size whichtoken = 1;
				command->initialize_from_tokens(decoy, whichtoken, resid );
				command->residue_action(*task, resid); //not sure about this token thing...
			}
		} else {
			for (utility::vector1<Size>::iterator it= corePos.begin(), end=corePos.end(); it!=end; it++){
					ResfileCommandOP command( new core::pack::task::NATRO );
					utility::vector1<std::string> decoy;
					decoy.push_back("NATRO");
					Size resid = *it;
					TR << "NATRO " <<  "(" << pose.residue(resid).name() << ")"<<  std::endl;
					Size whichtoken = 1;
					command->initialize_from_tokens(decoy, whichtoken, resid );
					command->residue_action(*task, resid); //not sure about this token thing...
				}
		}

		if ( boundary ){
			for (utility::vector1<Size>::iterator it= boundaryPos.begin(), end=boundaryPos.end(); it!=end; it++){
				ResfileCommandOP command( new core::pack::task::ALLAAxc ); //note! no cys
				utility::vector1<std::string> decoy;
				decoy.push_back("ALLAAxc");
				Size resid = *it;
				TR << "ALLAAxc " << resid << std::endl;
				Size whichtoken = 1;
				command->initialize_from_tokens( decoy, whichtoken, resid);
				command->residue_action(*task, resid); //not sure about this token thing...
			}
		} else {
		for (utility::vector1<Size>::iterator it= boundaryPos.begin(), end=boundaryPos.end(); it!=end; it++){
					ResfileCommandOP command( new core::pack::task::NATRO );
					utility::vector1<std::string> decoy;
					decoy.push_back("NATRO");
					Size resid = *it;
					TR << "NATRO " << resid << "(" << pose.residue(resid).name() << ")" << std::endl;
					Size whichtoken = 1;
					command->initialize_from_tokens( decoy, whichtoken, resid);
					command->residue_action( *task, resid); //not sure about this token thing...
				}
		}


		if ( surface ){
			for (utility::vector1<Size>::iterator it= surfPos.begin(), end=surfPos.end(); it!=end; it++){
				ResfileCommandOP command( new core::pack::task::POLAR );
				utility::vector1<std::string> decoy;
				decoy.push_back("POLAR");
				Size resid = *it;
				TR << "POLAR " << resid << std::endl;
				Size whichtoken = 1;
				command->initialize_from_tokens(decoy, whichtoken, resid);
				command->residue_action( *task, resid); //not sure about this token thing...
			}
		} else {
			for (utility::vector1<Size>::iterator it= surfPos.begin(), end=surfPos.end(); it!=end; it++){
						ResfileCommandOP command( new core::pack::task::NATRO );
						utility::vector1<std::string> decoy;
						decoy.push_back("NATRO");
						Size resid = *it;
						TR << "NATRO " << resid << " (" << pose.residue(resid).name() << ")" << std::endl;
						Size whichtoken = 1;
						command->initialize_from_tokens(decoy, whichtoken, resid);
						command->residue_action(*task, resid); //not sure about this token thing...
			}
		}

	//std::cout << "END REDUCE" << std::endl;
	//std::cout << *task << std::endl;

}

core::Real build_and_score_disulfide(core::pose::Pose & blank_pose, core::scoring::ScoreFunctionOP sfxn, const bool relax_bb, core::Size const res1, core::Size const res2) {

	
	//blank_pose.dump_pdb("pre_disulf.pdb");

	core::conformation::Residue old_res1 = blank_pose.residue(res1);
	core::conformation::Residue old_res2 = blank_pose.residue(res2);

	core::kinematics::MoveMapOP mm = core::kinematics::MoveMapOP( new core::kinematics::MoveMap());
	mm->set_bb(relax_bb);
	mm->set_chi(true);
	//utility::vector1<std::pair<core::Size, core::Size> > disulf;
	//disulf.push_back(std::pair<core::Size, core::Size>(res1, res2));

	//make_disulfide(blank_pose, disulf, mm);
	
	core::conformation::form_disulfide(blank_pose.conformation(), res1, res2);
	core::util:: rebuild_disulfide(blank_pose, res1,res2, NULL /*task*/, NULL /*scfxn*/, mm, NULL /*min scfxn*/);

	core::Real const score = (*sfxn)(blank_pose) * 0.50;
	//sfxn->show( TR, blank_pose );
	//blank_pose.dump_pdb("with_disulf.pdb");


	blank_pose.replace_residue(res1, old_res1, false);
	blank_pose.replace_residue(res2, old_res2, false);
	//blank_pose.dump_pdb("post_disulf.pdb");

	return score;
}

bool RemodelDesignMover::find_disulfides_in_the_neighborhood(Pose & pose, utility::vector1<std::pair<Size,Size> > & disulf_partners, const core::Energy & match_rt_limit, const Size & rosetta_scripts_min_loop,
															const bool & rosetta_scripts_include_current_ds, const bool & rosetta_scripts_keep_current_ds,
															const bool & relax_bb_for_disulf, const bool & use_match_rt, const bool & use_disulf_fa_score, const core::Real disulf_fa_max) {

	using core::scoring::disulfides::DisulfideMatchingPotential;
	using namespace core::chemical;

	bool pass=0;

	DisulfideMatchingPotential disulfPot;
	core::Energy  match_t;
	core::Energy  match_r;
	core::Energy  match_rt;

	// initialize default
	Size landingRangeStart = 1;
	Size landingRangeStop = pose.total_residue();

	// alternatively via blueprint
	if (remodel_data_.disulfLandingRange.size() != 0){
		landingRangeStart = remodel_data_.disulfLandingRange[0];
		landingRangeStop = remodel_data_.disulfLandingRange[1];
		TR << "Assigning Landing Range by Blueprint: " << landingRangeStart << " to " << landingRangeStop << std::endl;
	}

	if(option[OptionKeys::remodel::disulf_landing_range].user()){ //overwrite ranges if existed
		landingRangeStart =option[OptionKeys::remodel::disulf_landing_range][1];
		landingRangeStop =option[OptionKeys::remodel::disulf_landing_range][2];
		TR << "Assigning Landing Range by Arguments: " << landingRangeStart << " to " << landingRangeStop << std::endl;
	}


	TR << "FINDING DISULF" << std::endl;
	utility::vector1< bool > modeled_clusters( pose.total_residue(), false );
	utility::vector1< bool > residue_clusters( pose.total_residue(), false );
	run_calculator( pose, "neighborhood_calc", "neighbors", residue_clusters );
	run_calculator( pose, "neighborhood_calc", "central_residues", modeled_clusters );

	TR << "residue_clusters: [ ";
	for ( Size ii = 1; ii <= residue_clusters.size(); ++ii ) {
		if ( residue_clusters[ ii ] == true ) {
			TR << ii << " ";
		}
	}
	TR << "]" << std::endl;

	TR << "modeled_clusters: [ ";
	for ( Size ii = 1; ii <= modeled_clusters.size(); ++ii ) {
		if ( modeled_clusters[ ii ] == true ) {
			TR << ii << " ";
		}
	}
	TR << "]" << std::endl;

	// manual overwrite of the disulfide mobile range
	if ( remodel_data_.disulfMobileRange.size() != 0 ) {
		Size i = 1;
		for ( utility::vector1_bool::iterator itr = modeled_clusters.begin(); itr != modeled_clusters.end(); itr++) {

			*itr = false;
			if ( i == remodel_data_.disulfMobileRange[0] ) {
				*itr = true;
				TR << "Use disulf mobile range start: " << i << std::endl;
		 	} else if ( i > remodel_data_.disulfMobileRange[0] && i < remodel_data_.disulfMobileRange[1] ) {
				*itr = true;
			}
			else if ( i == remodel_data_.disulfMobileRange[1] ){
				*itr = true;
				TR << "Use disulf mobile range stop: " << i << std::endl;
			}
			i++;
		}
	} else {
		TR << "RemodelData disulfMobileRange not overwritten because it was not initialized previously." << std::endl;
	}

	// figure out which positions are "central" positions - I presume these are positions from which DS bonds can emanate.
	// then figure out which positions are not "central" but still "modeled". I assume these are the disulfide landing range
	// positions.
	utility::vector1<Size> cen_res;
	utility::vector1<Size> nbr_res;
	for ( Size ii = 1; ii <= modeled_clusters.size(); ++ii ) {
		if ( modeled_clusters[ ii ] == 1 ) {
			TR << "central " << ii <<  std::endl;
			cen_res.push_back( ii );
		}
	//	if ( modeled_clusters[ ii ] == 0 && residue_clusters[ ii ] == 1 ) {
	if ( residue_clusters[ii] == 1 ) {
			TR << "neighbor " << ii <<  std::endl;
			nbr_res.push_back( ii );
		}
	}

	TR << "central residues: ";
	for ( utility::vector1<Size>::iterator itr = cen_res.begin(), end=cen_res.end(); itr!=end; itr++ ) {
		TR << *itr << ",";
	}
	TR << std::endl;

	TR << "neighbor residues: ";
	for ( utility::vector1<Size>::iterator itr = nbr_res.begin(), end=nbr_res.end(); itr!=end; itr++ ) {
		TR <<  *itr << ",";
	}
	TR << std::endl;

	core::pose::Pose pose_copy = pose;
	for (core::Size i = 1; i <= pose_copy.total_residue(); ++i) {
		if (pose_copy.residue(i).name3() != "GLY") {
			protocols::simple_moves::MutateResidue make_ala(i,"ALA");
			make_ala.apply(pose_copy);
			TR << "Mutating residue " << i << " to ALA" << std::endl;
		}
	}
	core::scoring::ScoreFunctionOP sfxn_disulfide_only = core::scoring::ScoreFunctionOP(new core::scoring::ScoreFunction());
	sfxn_disulfide_only->set_weight(core::scoring::dslf_fa13, 1.0);


	for ( utility::vector1<Size>::iterator itr = cen_res.begin(), end=cen_res.end(); itr!=end; itr++ ) {
		for ( utility::vector1<Size>::iterator itr2 = nbr_res.begin(), end2=nbr_res.end(); itr2!=end2 ; itr2++ ) {
			if ((nbr_res != cen_res || (*itr2 > (*itr + rosetta_scripts_min_loop))) &&
					std::abs( core::SSize(*itr2 - *itr) )  > std::abs( core::SSize(rosetta_scripts_min_loop) ) &&
				(*itr2) <= landingRangeStop && (*itr2) >= landingRangeStart) {
				TR << "DISULF trying disulfide between " << *itr << " and " << *itr2 << std::endl;
				// distance check
				if ( pose.residue(*itr).aa() != aa_gly && pose.residue(*itr2).aa() != aa_gly ) {
					//TR << "DISULF " <<  *itr << "x" << *itr2 << std::endl;
					if (!rosetta_scripts_keep_current_ds || ( pose.residue(*itr).name() != "CYD" && pose.residue(*itr2).name() != "CYD")) {
						Real dist_squared = pose.residue(*itr).xyz("CB").distance_squared(pose.residue(*itr2).xyz("CB"));
						if ( dist_squared > 25 ) {
							TR << "DISULF \tTOO FAR. CB-CB distance squared: " << dist_squared << std::endl;

						} else {
							//if ( match_rt < option[OptionKeys::remodel::match_rt_limit] && std::abs( seqGap ) > 1 && (*itr2) <= landingRangeStop && (*itr2) >= landingRangeStart ) {
							if (use_disulf_fa_score == true) {
								core::Real const disulfide_fa_score = build_and_score_disulfide(pose_copy, sfxn_disulfide_only, relax_bb_for_disulf, *itr, *itr2);
								TR << "DISULF FA SCORE RES " << *itr << " " << *itr2 << " " << disulfide_fa_score << std::endl;
								if ( disulfide_fa_score < disulf_fa_max   ) {
									TR << "DISULF possible " << dist_squared << std::endl;
									TR << "DISULF " <<  *itr << "x" << *itr2 << std::endl;
									//TR << " " << match_rt << std::endl;
									std::pair< Size, Size > temp_pair;
									std::pair< Size, Size > alt_pair;
									
									temp_pair = std::make_pair( *itr, *itr2 );
									alt_pair = std::make_pair( *itr2, *itr );
									
									if (std::find(disulf_partners.begin(), disulf_partners.end(), temp_pair) == disulf_partners.end() &&
										std::find(disulf_partners.begin(), disulf_partners.end(), alt_pair) == disulf_partners.end()) {
										disulf_partners.push_back( temp_pair );									
									}

										


									pass = 1;
								} else {
									if (rosetta_scripts_include_current_ds && pose.residue(*itr).is_bonded(pose.residue(*itr2))) {
										TR << "DISULF \tIncluding pre-existing disulfide despite failed disulf_fa_max check." << std::endl;
										std::pair< Size, Size > temp_pair;
										std::pair< Size, Size > alt_pair;
										
										temp_pair = std::make_pair( *itr, *itr2 );
										alt_pair = std::make_pair( *itr2, *itr );
										
										if (std::find(disulf_partners.begin(), disulf_partners.end(), temp_pair) == disulf_partners.end() &&
											std::find(disulf_partners.begin(), disulf_partners.end(), alt_pair) == disulf_partners.end()) {
											disulf_partners.push_back( temp_pair );									
										}
									} else {
										TR << "DISULF \tFailed disulf_fa_max check." << std::endl;
									}
								}							


							}



							if (use_match_rt == true) {
								disulfPot.score_disulfide( pose.residue(*itr), pose.residue(*itr2), match_t, match_r, match_rt );
								TR << "DISULF \tmatch_t: " << match_t << ", match_r: " << match_r << ", match_rt: " << match_rt << std::endl;
								if ( match_rt < match_rt_limit  ) {
									TR << "DISULF possible " << dist_squared << std::endl;
									TR << "DISULF " <<  *itr << "x" << *itr2 << std::endl;
									TR << "match_rt " << match_rt << std::endl;
									std::pair< Size, Size > temp_pair;
									std::pair< Size, Size > alt_pair;

									temp_pair = std::make_pair( *itr, *itr2 );
									alt_pair = std::make_pair( *itr2, *itr );

									if (std::find(disulf_partners.begin(), disulf_partners.end(), temp_pair) == disulf_partners.end() &&
										std::find(disulf_partners.begin(), disulf_partners.end(), alt_pair) == disulf_partners.end()) {
										disulf_partners.push_back( temp_pair );
									}

										


									pass = 1;
								} else {
									if (rosetta_scripts_include_current_ds && pose.residue(*itr).is_bonded(pose.residue(*itr2))) {
										TR << "DISULF \tIncluding pre-existing disulfide despite failed match_rt_limit check." << std::endl;
										std::pair< Size, Size > temp_pair;
										std::pair< Size, Size > alt_pair;
										
										temp_pair = std::make_pair( *itr, *itr2 );
										alt_pair = std::make_pair( *itr2, *itr );
										
										if (std::find(disulf_partners.begin(), disulf_partners.end(), temp_pair) == disulf_partners.end() &&
											std::find(disulf_partners.begin(), disulf_partners.end(), alt_pair) == disulf_partners.end()) {
											disulf_partners.push_back( temp_pair );									
										}
									} else {
										TR << "DISULF \tFailed match_rt_limit check." << std::endl;
									}
								}
							}

						}

					} else {
						TR <<"DISULF \tkeep_current_ds set to True, skipping residues that are already in disulfides." << std::endl;
					}
				} else {
					TR << "DISULF \tFailed glycine check. Cysteines cannot replace glycine positions in the native pose." << std::endl;
				}
			}
		}
	}

	return pass;

}

void RemodelDesignMover::make_disulfide(Pose & pose, utility::vector1<std::pair<Size, Size> > & disulf_partners, core::kinematics::MoveMapOP mm){
	//utility::vector1<std::pair<Size,Size>> dummy_vector;
	for (utility::vector1<std::pair<Size,Size> >::iterator itr = disulf_partners.begin(); itr != disulf_partners.end(); itr++){
		core::conformation::form_disulfide(pose.conformation(), (*itr).first, (*itr).second);
		core::util:: rebuild_disulfide(pose, (*itr).first,(*itr).second, NULL /*task*/, NULL /*scfxn*/, mm, NULL /*min scfxn*/);
		TR << "build_disulf between " << (*itr).first << " and " << (*itr).second << std::endl;
	//	pose.dump_pdb("disulf.pdb");
		}
}

void RemodelDesignMover::make_disulfide_fast(Pose & pose, utility::vector1<std::pair<Size, Size> > & disulf_partners){
	//utility::vector1<std::pair<Size,Size>> dummy_vector;
	for (utility::vector1<std::pair<Size,Size> >::iterator itr = disulf_partners.begin(); itr != disulf_partners.end(); itr++){
		core::conformation::form_disulfide(pose.conformation(), (*itr).first, (*itr).second);
		TR << "build_disulf between " << (*itr).first << " and " << (*itr).second << std::endl;
		}
}

/// these are split up for convenience reasons, so one can bypass blueprint setting if needed be
void RemodelDesignMover::mode1_packertask(Pose & pose){ // auto loop only
	TR << "MODE 1: AUTO DESIGN of remodeled regions only" << std::endl;

  core::pack::task::TaskFactoryOP TF = protocols::forge::methods::remodel_generic_taskfactory();

	//if need more operations added, put them here.

	//create the real task
  working_model_.task = TF->create_task_and_apply_taskoperations( pose );
	utility::vector1_bool additional_sites(pose.total_residue(), false);

	run_calculator(pose, "neighborhood_calc", "central_residues", additional_sites);

	working_model_.task->restrict_to_residues( additional_sites );

	TR << "number to be packed after adding sites: " << working_model_.task->num_to_be_packed() << std::endl;

}

/// these are split up for convenience reasons, so one can bypass blueprint setting if needed be
void RemodelDesignMover::mode1_1_packertask(Pose & pose){ // auto loop only
	TR << "MODE 1.1: AUTO DESIGN everything -- denovo" << std::endl;

  core::pack::task::TaskFactoryOP TF = protocols::forge::methods::remodel_generic_taskfactory();

	//if need more operations added, put them here.

	//create the real task
  working_model_.task = TF->create_task_and_apply_taskoperations( pose );
	utility::vector1_bool additional_sites(pose.total_residue(), true);

	working_model_.task->restrict_to_residues( additional_sites );

	TR << "number to be packed after adding sites: " << working_model_.task->num_to_be_packed() << std::endl;

}

void RemodelDesignMover::mode2_packertask(Pose & pose){ // auto loop with design neighbor
	TR << "MODE 2: AUTO DESIGN of remodeled regions and neighbors" << std::endl;

  core::pack::task::TaskFactoryOP TF = protocols::forge::methods::remodel_generic_taskfactory();

	//if need more operations added, put them here.

	//create the real task
  working_model_.task = TF->create_task_and_apply_taskoperations( pose );
	utility::vector1_bool additional_sites(pose.total_residue(), false);

	run_calculator(pose, "neighborhood_calc", "neighbors", additional_sites);

	working_model_.task->restrict_to_residues( additional_sites );

	TR << "number to be packed after adding sites: " << working_model_.task->num_to_be_packed() << std::endl;

}

void RemodelDesignMover::mode3_packertask(Pose & pose){ // auto loop with repack neighbor
  TR << "MODE 3: AUTO DESIGN of remodeled regions and repack neighbors only" << std::endl;

  core::pack::task::TaskFactoryOP TF = protocols::forge::methods::remodel_generic_taskfactory();

	//if need more operations added, put them here.

	//create the real task
  working_model_.task = TF->create_task_and_apply_taskoperations( pose );
	utility::vector1_bool additional_sites(pose.total_residue(), false);
	utility::vector1_bool core_sites(pose.total_residue(), false);

	run_calculator(pose, "neighborhood_calc", "neighbors", additional_sites);
	run_calculator(pose, "neighborhood_calc", "central_residues", core_sites);

	working_model_.task->restrict_to_residues( additional_sites );

	//lock down the rest of the positions to repack only
	for (Size ii = 1; ii <= core_sites.size();++ii){
		if ( !core_sites[ii] ){
			if (working_model_.task->nonconst_residue_task(ii).being_packed()){
				working_model_.task->nonconst_residue_task(ii).restrict_to_repacking();
				//debug	TR << "restrict position " << ii << " to repack only" << std::endl;
			}
		}
	}

	TR << "number to be packed after adding sites: " << working_model_.task->num_to_be_packed() << std::endl;

}


void RemodelDesignMover::mode4_packertask(Pose & pose){ // full manuaing namespace core::scoring;

	TR << "MODE 4: Manual DESIGN REMODEL" << std::endl;
	working_model_.manualPackerTaskGen(pose, remodel_data_);

}


void RemodelDesignMover::mode5_packertask(Pose & pose){ // manual with auto design neighbor

  TR << "MODE 5: Manual DESIGN REMODEL with DESIGN NEIGHBOR" << std::endl;

  core::pack::task::TaskFactoryOP TF = protocols::forge::methods::remodel_generic_taskfactory();

	//if need more operations added, put them here.

	//create the real task
  working_model_.task = TF->create_task_and_apply_taskoperations( pose );

	utility::vector1_bool additional_sites(pose.total_residue(), false);

	//debug	run_calculator(pose, "neighborhood_calc", "central_residues", additional_sites);
	run_calculator(pose, "neighborhood_calc", "neighbors", additional_sites);

	// identify all the positions to change -- this includes central_residues and
	// neighbors from calculator runs
	working_model_.task->restrict_to_residues( additional_sites );

	// process the information in blueprint and save the positions touched.
  non_default_positions_ = protocols::forge::methods::parse_resfile_string_with_no_lockdown(pose, *working_model_.task, remodel_data_.parsed_string_for_resfile );

	TR << "number to be packed after adding sites: " << working_model_.task->num_to_be_packed() << std::endl;
}


void RemodelDesignMover::mode6_packertask(Pose & pose){ // manual with auto repack neighbor

  //build the base task with positions designed
  mode5_packertask(pose);

  TR << "MODE 6: Manual DESIGN REMODEL with REPACK NEIGHBOR" << std::endl;

	//lock down the rest of the positions to repack only
	for (Size ii = 1; ii <= non_default_positions_.size();++ii){
		if ( !non_default_positions_[ii] ){
			if (working_model_.task->nonconst_residue_task(ii).being_packed()){
				working_model_.task->nonconst_residue_task(ii).restrict_to_repacking();
				//debug	TR << "restrict position " << ii << " to repack only" << std::endl;
			}
		}
	}
}



} // remodel
} // forge
} // protocol
