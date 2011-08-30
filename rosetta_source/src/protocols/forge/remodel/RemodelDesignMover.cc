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

// unit headers
#include <protocols/forge/remodel/RemodelDesignMover.hh>
#include <protocols/forge/remodel/RemodelMover.hh>
#include <protocols/forge/methods/util.hh>

// package headers

// project headers
#include <basic/Tracer.hh>
#include <basic/MetricValue.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/disulfides/DisulfideMatchingPotential.hh>
#include <core/util/disulfide_util.hh>
#include <core/conformation/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/pack/task/PackerTask_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborhoodByDistanceCalculator.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>
#include <core/pack/task/ResfileReader.hh>
#include <protocols/forge/build/BuildInstruction.hh> // REQUIRED FOR WINDOWS

// numeric headers
#include <numeric/random/random.hh>

// boost headers

// C++ headers


namespace protocols {
namespace forge{
namespace remodel{

// Tracer instance for this file
// Named after the original location of this code
static basic::Tracer TR( "protocols.forge.remodel.RemodelDesignMover" );

// RNG
static numeric::random::RandomGenerator RG( 2342342 ); // magic number, don't change


// @brief default constructor
RemodelDesignMover::RemodelDesignMover(){
// has to reinitialize state before apply
	state_.clear();
}

/// @brief value constructor
RemodelDesignMover::RemodelDesignMover(RemodelData remodel_data, RemodelWorkingSet working_model, ScoreFunctionOP sfxn)
{

	using core::pose::metrics::CalculatorFactory;
	using protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator;

  remodel_data_ = remodel_data;
	working_model_ = working_model;
	archived_starting_task_ = working_model.task;
	score_fxn_ = sfxn->clone();

  // setup calculators
  CalculatorFactory::Instance().remove_calculator( "neighborhood_calc" );
  CalculatorFactory::Instance().register_calculator(
    "neighborhood_calc",
    new NeighborhoodByDistanceCalculator( working_model.manager.union_of_intervals_containing_undefined_positions())
  );

}
/// @brief copy constructor

/// @brief default destructor
  RemodelDesignMover::~RemodelDesignMover(){}

/// @brief clone this object
  RemodelDesignMover::MoverOP RemodelDesignMover::clone() {
  return new RemodelDesignMover( *this );
}


/// @brief create this type of object
  RemodelDesignMover::MoverOP RemodelDesignMover::fresh_instance() {
  return new RemodelDesignMover();
}

void RemodelDesignMover::apply( Pose & pose )
{
	// make decision as to which mode to apply
	bool manual = remodel_data_.has_design_info_;
	bool neighbor = basic::options::option[basic::options::OptionKeys::remodel::design::find_neighbors].user();
	bool design = basic::options::option[basic::options::OptionKeys::remodel::design::design_neighbors].user();

	if (!check_state()){
		basic::Error() << "check_state failed, has to set_state first " << std::endl;
	}

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
			mode1_packertask(pose);
		}
	}
	if (!strcmp(state_.c_str(), "stage")){
    if (!basic::options::option[basic::options::OptionKeys::remodel::design::skip_partial].user()){
			reduce_task(pose, working_model_.task, true, true, false);
		}
	}
	else if (!strcmp(state_.c_str(), "finish")){
		// if finishing design, no need to reduce, but require resetting positios
		//reduce_task(pose, working_model_.task, true, true, true);
	}

/*
	//mode4_packertask(pose);
	//mode5_packertask(pose);
	//mode6_packertask(pose);
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH );
	core::pack::pack_rotamers(pose, *scorefxn , working_model_.task);

	mode5_packertask(pose);
	core::pack::pack_rotamers(pose, *scorefxn , working_model_.task);
	mode5_packertask(pose);
	core::pack::pack_rotamers(pose, *scorefxn , working_model_.task);
*/
/*
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::STANDARD_WTS, core::scoring::SCORE12_PATCH );
	if (basic::options::option[basic::options::OptionKeys::packing::soft_rep_design]){
			TR << "SWITCHING REMODEL DESIGN MOVER SCOREFUNCTION TO SOFT_REP_DESIGN" << std::endl;
			scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::SOFT_REP_DESIGN_WTS);
	}
	*/
	//debug
//TR <<  *working_model_.task << std::endl;
	core::pack::pack_rotamers(pose, *score_fxn_ , working_model_.task);
// pose.dump_pdb("junkCheck.pdb");
//	core::pack::pack_rotamers(pose, *scorefxn , working_model_.task);
//	core::pack::pack_rotamers(pose, *scorefxn , working_model_.task);
 // core::pack::task::TaskFactoryOP TF = new core::pack::task::TaskFactory;
  //TF->create_packer_task(pose);
//	core::pack::pack_rotamers(pose, *scorefxn , TF->create_packer_task(pose));
//	core::pack::pack_rotamers(pose, *scorefxn , TF->create_packer_task(pose));
//	core::pack::pack_rotamers(pose, *scorefxn , TF->create_packer_task(pose));
}

std::string
RemodelDesignMover::get_name() const {
	return "RemodelDesignMover";
}

bool RemodelDesignMover::check_state(){
	if (state_.empty()) {
		return false;
		TR << "state tag not set " << std::endl;
	} else {
		TR << "Design Mover State: " << state_ << std::endl;
		return true;
	}
}

void RemodelDesignMover::set_state( std::string state_tag ){
		state_ = state_tag;
		//reset the task
		working_model_.task = archived_starting_task_;
}


void
run_calculator(
									core::pose::Pose const & pose,
									std::string const & calculator,
									std::string const & calculation,
									utility::vector1_bool & residues )
{
  runtime_assert(residues.size() == pose.total_residue());

  //find the set of residues
  typedef std::set< core::Size > SizeSet;
  basic::MetricValue< SizeSet > mv_sizeset;
  pose.metric(calculator, calculation, mv_sizeset);
  SizeSet const & sizeset(mv_sizeset.value());
	//TR << "runCAlculator " << std::endl;
  //insert this into the vector
  for(SizeSet::const_iterator it(sizeset.begin()), end(sizeset.end()) ; it != end; ++it){
   // TR << *it <<  " debug run_calc " << std::endl;
    residues[*it] = true;  }

  return;
}

void RemodelDesignMover::reduce_task( Pose & pose, core::pack::task::PackerTaskOP &task, bool core, bool boundary, bool surface){

  // setup calculators
	using core::pose::metrics::CalculatorFactory;
  using protocols::toolbox::pose_metric_calculators::NeighborhoodByDistanceCalculator;
	using core::pack::task::ResfileCommandOP;

  // need to collect all positions before reduction
  CalculatorFactory::Instance().remove_calculator( "reducetask_calc" );
	utility::vector1_bool boollist(pose.total_residue());
  run_calculator(pose, "neighborhood_calc", "neighbors", boollist);
	std::set<Size> positionList;
	for (Size i = 1; i <= pose.total_residue(); i++){
		if (boollist[i]){
			positionList.insert(i);
		}
	}

  CalculatorFactory::Instance().register_calculator(
    "reducetask_calc",
    new NeighborhoodByDistanceCalculator( positionList )
  );

	utility::vector1_bool resclass(pose.total_residue(),false);
	utility::vector1_bool neighbor_count(pose.total_residue(),false);
	utility::vector1<Size> corePos;
	utility::vector1<Size> boundaryPos;
	utility::vector1<Size> surfPos;

	Size CORE_CUTOFF=15;
	Size BOUNDARY_CUTOFF=10;

	//get num_neighbors for each position
	basic::MetricValue< std::map< core::Size, core::Size > >nbr_map;
	pose.metric("reducetask_calc", "num_neighbors_map", nbr_map);
  std::map< core::Size, core::Size > sizemap=nbr_map.value();
	for (Size i = 1; i<= resclass.size(); i++){ //check everyposition in the packertask
		if (task->nonconst_residue_task(i).being_packed() && sizemap[i]){
			TR << "touch position " << i << std::endl;
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
/*
		//debug:
		for (utility::vector1<Size>::iterator it= corePos.begin(), end=corePos.end(); it!=end; it++){
			TR << "DEBUG: core positions:" << *it << std::endl;
		}
		for (utility::vector1<Size>::iterator it= boundaryPos.begin(), end=boundaryPos.end(); it!=end; it++){
			TR << "DEBUG: boundary positions:" << *it << std::endl;
		}
		for (utility::vector1<Size>::iterator it= surfPos.begin(), end=surfPos.end(); it!=end; it++){
			TR << "DEBUG: surface positions:" << *it << std::endl;
		}
*/
		//build new reduced task
//		core::pack::task::TaskFactoryOP TF = protocols::forge::methods::remodel_generic_taskfactory();
//		task=TF->create_task_and_apply_taskoperations( pose );

		//borrow RESFILE command for the job -- save a lot of coding effort.
		if ( core ){
			for (utility::vector1<Size>::iterator it= corePos.begin(), end=corePos.end(); it!=end; it++){
				ResfileCommandOP command  = new core::pack::task::APOLAR;
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
					ResfileCommandOP command  = new core::pack::task::NATRO;
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
				ResfileCommandOP command  = new core::pack::task::ALLAAxc; //note! no cys
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
					ResfileCommandOP command  = new core::pack::task::NATRO;
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
				ResfileCommandOP command  = new core::pack::task::POLAR;
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
						ResfileCommandOP command  = new core::pack::task::NATRO;
						utility::vector1<std::string> decoy;
						decoy.push_back("NATRO");
						Size resid = *it;
						TR << "NATRO " << resid << resid << "(" << pose.residue(resid).name() << ")" << std::endl;
						Size whichtoken = 1;
						command->initialize_from_tokens(decoy, whichtoken, resid);
						command->residue_action(*task, resid); //not sure about this token thing...
			}
		}


}


bool RemodelDesignMover::find_disulfides_in_the_neighborhood(Pose & pose, utility::vector1<std::pair<Size,Size> > & disulf_partners){

	using core::scoring::disulfides::DisulfideMatchingPotential;
	using namespace basic::options;
	using namespace core::chemical;

	bool pass=0;

	DisulfideMatchingPotential disulfPot;
	core::Energy  match_t;
	core::Energy  match_r;
	core::Energy  match_rt;

//initialize default
	Size landingRangeStart = 1;
	Size landingRangeStop = pose.total_residue();

//alternatively via blueprint
	if (remodel_data_.disulfLandingRange.size() != 0){
		landingRangeStart = remodel_data_.disulfLandingRange[0];
		landingRangeStop = remodel_data_.disulfLandingRange[1];
		TR << "Assignging Landing Range by Blueprint: " << landingRangeStart << " to " << landingRangeStop << std::endl;
	}

	if (option[  OptionKeys::remodel::disulf_landing_range].user()){ //overwrite ranges if existed
		landingRangeStart = option[ OptionKeys::remodel::disulf_landing_range]()[1];
		landingRangeStop = option[ OptionKeys::remodel::disulf_landing_range]()[2];
		TR << "Assignging Landing Range by Arguments: " << landingRangeStart << " to " << landingRangeStop << std::endl;
	}


	TR << "FINDING DISULF" << std::endl;
	utility::vector1_bool modeled_clusters(pose.total_residue(), false);
	utility::vector1_bool residue_clusters(pose.total_residue(), false);
	utility::vector1<Size> cen_res;
	utility::vector1<Size> nbr_res;
	run_calculator(pose, "neighborhood_calc", "neighbors", residue_clusters);
	run_calculator(pose, "neighborhood_calc", "central_residues", modeled_clusters);
//debug
	Size i=1;
	for (utility::vector1_bool::iterator itr=modeled_clusters.begin(), end=modeled_clusters.end();  itr !=end; itr++){
	//TR << *itr<< std::endl;
		if (modeled_clusters[i] == 0 && residue_clusters[i] == 1) {
	//		TR << "neighbor " << i <<  std::endl;
			nbr_res.push_back(i);
		}
		if (modeled_clusters[i] == 1){
	//		TR << "central " << i <<  std::endl;
			cen_res.push_back(i);
		}
	i++;
	}
	TR << "central residues: ";
	for (utility::vector1<Size>::iterator itr = cen_res.begin(), end=cen_res.end(); itr!=end; itr++){
		TR << *itr << ",";
	}
	TR << std::endl;

	TR << "neighbor residues: ";
	for (utility::vector1<Size>::iterator itr = nbr_res.begin(), end=nbr_res.end(); itr!=end; itr++){
		TR <<  *itr << ",";
	}
	TR << std::endl;

	for (utility::vector1<Size>::iterator itr = cen_res.begin(), end=cen_res.end(); itr!=end; itr++){
		for (utility::vector1<Size>::iterator itr2 = nbr_res.begin(), end2=nbr_res.end(); itr2!=end2 ; itr2++){
			//TR << "DISULF " <<  *itr << "x" << *itr2 << std::endl;
			//disance check
			if (pose.residue(*itr).aa() != aa_gly && pose.residue(*itr2).aa() != aa_gly){
			//TR << "DISULF " <<  *itr << "x" << *itr2 << std::endl;
			Real dist = pose.residue(*itr).xyz("CB").distance_squared(pose.residue(*itr2).xyz("CB"));
			if (dist > 25){
			//TR << "TOO FAR " << dist << std::endl;
				}
				else {
					disulfPot.score_disulfide(pose.residue(*itr), pose.residue(*itr2), match_t, match_r, match_rt);
					//TR << "match_t " << match_t << std::endl;
					//TR << "match_r " << match_r << std::endl;
					int seqGap = (int)*itr2-(int)*itr;
					if (match_rt < option[OptionKeys::remodel::match_rt_limit] && std::abs(seqGap)>1 && (*itr2) <= landingRangeStop && (*itr2) >= landingRangeStart){
					TR << "DISULF possible " << dist << std::endl;
					TR << "DISULF " <<  *itr << "x" << *itr2 << std::endl;
					TR << "match_rt " << match_rt << std::endl;
						std::pair<Size, Size> temp_pair;
						temp_pair = std::make_pair(*itr, *itr2);
						disulf_partners.push_back(temp_pair);
						pass = 1;
					}
				}
			}
		}
	}

	return pass;

}

void RemodelDesignMover::make_disulfide(Pose & pose, utility::vector1<std::pair<Size, Size> > & disulf_partners){
	//utility::vector1<std::pair<Size,Size>> dummy_vector;
	for (utility::vector1<std::pair<Size,Size> >::iterator itr = disulf_partners.begin(); itr != disulf_partners.end(); itr++){
		core::conformation::form_disulfide(pose.conformation(), (*itr).first, (*itr).second);
		core::util:: rebuild_disulfide(pose, (*itr).first,(*itr).second);
	TR << "build_disulf between " << (*itr).first << " and " << (*itr).second << std::endl;
	//	pose.dump_pdb("disulf.pdb");
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
core::pack::task::PackerTaskOP &
RemodelDesignMover::task(){
	return working_model_.task;
}

void RemodelDesignMover::scorefunction( ScoreFunctionOP sfxn) {
	score_fxn_ = sfxn->clone();
}




} // remodel
} // forge
} // protocol
