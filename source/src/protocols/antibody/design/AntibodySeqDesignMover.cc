// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/AntibodySeqDesignMover.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/AntibodySeqDesignMover.hh>
#include <protocols/antibody/design/ConservativeDesignOperation.hh>
#include <protocols/antibody/design/ResidueProbDesignOperation.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/design/AntibodySeqDesignTFCreator.hh>

#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/clusters/util.hh>
#include <protocols/antibody/constraints/util.hh>
#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.hh>
#include <protocols/antibody/util.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/util.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/toolbox/task_operations/LimitAromaChi2Operation.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/flxbb/FlxbbDesign.hh>
#include <protocols/flxbb/DesignTask.hh>
#include <protocols/moves/MonteCarlo.hh>

// Basic Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.antibody.design.AntibodySeqDesignMover" );

namespace protocols {
namespace antibody {
namespace design {
	using namespace protocols::antibody;
	using namespace protocols::antibody::clusters;
	using namespace protocols::antibody::constraints;

	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;
	using namespace core::kinematics;


AntibodySeqDesignMover::AntibodySeqDesignMover(AntibodyInfoOP ab_info):
	protocols::moves::Mover("AntibodySeqDesignMover"),
		ab_info_(ab_info),
		scorefxn_(/* NULL */)

{
	if (ab_info_->get_current_AntibodyNumberingScheme()!="AHO_Scheme" && ab_info_->get_current_CDRDefinition() != "North"){
		utility_exit_with_message("Antibody Design Protocol requires AHO_scheme and North definitions");
	}

	set_defaults();
	read_command_line_options();

	//Modeler only currently used to get distance for adding paratope-epitope constraints.
	modeler_ = AntibodyDesignModelerOP( new AntibodyDesignModeler(ab_info_) );

}

AntibodySeqDesignMover::AntibodySeqDesignMover(AntibodyInfoOP ab_info, std::string instruction_path):
	protocols::moves::Mover("AntibodySeqDesignMover"),
	ab_info_(ab_info),
	scorefxn_(/* NULL */),
	instruction_path_(instruction_path)

{
	ab_info_ = ab_info;
	if (ab_info_->get_current_AntibodyNumberingScheme()!="AHO_Scheme" && ab_info_->get_current_CDRDefinition() != "North"){
		utility_exit_with_message("Antibody Design Protocol requires AHO_scheme and North definitions");
	}

	set_defaults();
	read_command_line_options();

	//Modeler only currently used to get distance for adding paratope-epitope constraints.
	modeler_ = AntibodyDesignModelerOP( new AntibodyDesignModeler(ab_info_) );

}

void
AntibodySeqDesignMover::set_defaults() {

	set_use_cluster_constraints(true);
	paratope_cdrs_.clear();
	paratope_cdrs_.resize(6, true);
	epitope_residues_.clear();
	cdr_seq_design_options_.clear();


}

AntibodySeqDesignMover::~AntibodySeqDesignMover(){}

void
AntibodySeqDesignMover::read_command_line_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	default_instruction_path_ = basic::options::option [basic::options::OptionKeys::antibody::design::base_instructions]();
	neighbor_detection_dis(option [OptionKeys::antibody::design::neighbor_dis]());

	set_basic_design(option [OptionKeys::antibody::design::benchmark_basic_design]());
	set_design_method(design_type_from_string(option [OptionKeys::antibody::design::design_method]()));
	set_rounds(option [OptionKeys::antibody::design::design_rounds]());
	use_epitope_constraints_ = basic::options::option [basic::options::OptionKeys::antibody::design::use_epitope_constraints]();
	full_cdr_min_ = !basic::options::option [basic::options::OptionKeys::antibody::design::disable_full_cdr_min]();

}

void
AntibodySeqDesignMover::set_seq_designer_options(AntibodyCDRSeqDesignOptions options) {
	cdr_seq_design_options_ = options;
}

void
AntibodySeqDesignMover::setup_options_classes(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if (cdr_seq_design_options_.size() == 0){
		cdr_seq_design_options_ = protocols::antibody::design::get_seq_design_options();
	}
}

std::string
AntibodySeqDesignMover::get_name() const {
	return "AntibodySeqDesignMover";
}


void
AntibodySeqDesignMover::set_basic_design(const bool setting){
	basic_design_ = setting;
}

void
AntibodySeqDesignMover::neighbor_detection_dis(core::Real const neighbor_distance){
	neighbor_dis_ = neighbor_distance;
}

void
AntibodySeqDesignMover::set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn) {
	scorefxn_ = scorefxn->clone();
}

void
AntibodySeqDesignMover::set_rounds(const core::Size rounds){
	rounds_ = rounds;
}

void
AntibodySeqDesignMover::set_design_method(DesignTypeEnum const design_method){
	design_method_ = design_method;
}

void
AntibodySeqDesignMover::set_use_cluster_constraints(const bool setting) {
	use_cluster_constraints_ = setting;
}

void
AntibodySeqDesignMover::set_paratope_cdrs(const vector1<bool>& cdrs){
	if (cdrs.size() != 6){
		utility_exit_with_message("Passed paratope cdrs does not equal the total number of cdrs!");
	}
	paratope_cdrs_ = cdrs;
}

void
AntibodySeqDesignMover::set_epitope_residues(vector1<PDBNumbering> epitope_residues){
	epitope_residues_ = epitope_residues;
}


void
AntibodySeqDesignMover::setup_cdr_constraints(core::pose::Pose & pose){


	using namespace protocols::antibody::constraints;

	if (use_cluster_constraints_){

		for (core::Size i = 1; i<= core::Size(CDRNameEnum_total); ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);

			if (cdr_has_res_constraints(ab_info_, pose, cdr, "DihedralConstraint")){
				continue;
			} else if (cdr_has_res_constraints(ab_info_, pose, cdr, "CoordinateConstraint")){
				continue;
			} else {
				add_harmonic_cluster_cst_or_dihedral_cst(ab_info_, pose, cdr);
			}

		}
	}

}

void
AntibodySeqDesignMover::setup_paratope_epitope_constraints(core::pose::Pose & pose){

		if (use_epitope_constraints_){
			paratope_epitope_cst_mover_ = protocols::antibody::constraints::ParatopeEpitopeSiteConstraintMoverOP( new ParatopeEpitopeSiteConstraintMover(ab_info_) );
			paratope_epitope_cst_mover_->constrain_to_paratope_cdrs(paratope_cdrs_);
			paratope_epitope_cst_mover_->constrain_to_epitope_residues(epitope_residues_, pose);
			paratope_epitope_cst_mover_->set_interface_distance(modeler_->interface_detection_dis());
			paratope_epitope_cst_mover_->apply(pose);
		}
		else {
			paratope_cst_mover_ = protocols::antibody::constraints::ParatopeSiteConstraintMoverOP( new ParatopeSiteConstraintMover(ab_info_) );
			paratope_cst_mover_->constrain_to_paratope_cdrs(paratope_cdrs_);
			paratope_cst_mover_->set_interface_distance(modeler_->interface_detection_dis());
			paratope_cst_mover_->apply(pose);
		}

}

void
AntibodySeqDesignMover::setup_scorefxn() {

	if (! scorefxn_ ){
		scorefxn_ = core::scoring::get_score_function();
		scorefxn_->apply_patch_from_file("antibody_design");
	}
}

void
AntibodySeqDesignMover::apply(core::pose::Pose& pose){

		using namespace protocols::antibody;
		using namespace protocols::antibody::design;
		using namespace core::pack::task;
		using namespace protocols::toolbox::task_operations;

		setup_options_classes();
		if (! protocols::antibody::clusters::check_if_pose_renumbered_for_clusters(pose)){
			utility_exit_with_message("PDB must be numbered correctly to identify North CDR clusters.  Please see Antibody Design documentation.");
		}


		//Move to show
		utility::vector1<std::string> strategies(SeqDesignStrategyEnum_total, "");
		strategies[ seq_design_profiles ] = "Profile-Based";
		strategies[ seq_design_conservative ] = "Conservative";
		strategies[ seq_design_basic ] = "Basic";

		TR<< "//////////////////////////////////////////////////////////////////////////////////////////////////" <<std::endl;
		for (core::Size i = 1; i <= core::Size(CDRNameEnum_total); ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			CDRSeqDesignOptionsCOP options = cdr_seq_design_options_[cdr];

			TR<<"//////// "<< ab_info_->get_CDR_name(cdr) <<" /////////////////////////////////////////////// "<<std::endl;
			TR<<"////" << std::endl;
			TR<< "///   Design? "<<std::boolalpha <<  options->design() << std::endl;
			if (options->design()){
				TR<< "////    Strategy: " << strategies[ options->design_strategy() ] << std::endl;
			}
			TR<<"///" << std::endl;
		}
		TR<< "//////////////////////////////////////////////////////////////////////////////////////////////////" <<std::endl;


		//////////Setup//////////////////////////


		ab_info_->setup_CDR_clusters(pose);

		AntibodySeqDesignTFCreatorOP tf_creator( new AntibodySeqDesignTFCreator(ab_info_, cdr_seq_design_options_) );
		tf_creator->neighbor_detection_dis(neighbor_dis_);
		tf_creator->set_basic_design(basic_design_);

		core::pack::task::TaskFactoryOP tf = tf_creator->generate_tf_seq_design(pose);
		protocols::simple_moves::PackRotamersMover packer = protocols::simple_moves::PackRotamersMover(scorefxn_);
		packer.task_factory(tf);


		setup_scorefxn();
		setup_cdr_constraints(pose);
		setup_paratope_epitope_constraints(pose);

		modeler_->set_scorefunction(scorefxn_);
		scorefxn_->show(pose);

		core::kinematics::MoveMapOP movemap( new MoveMap() );
		for (core::Size i = 1; i<=CDRNameEnum_total; ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			CDRSeqDesignOptionsCOP options = cdr_seq_design_options_[cdr];
			if (full_cdr_min_ || options->design()){
				ab_info_->add_CDR_to_MoveMap(pose, movemap, cdr, false, true, neighbor_dis_);
			}
		}

		//Setup MonteCarlo for >1 round
		protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo(pose, *scorefxn_, 1.0) );

		//Design methods
		if (design_method_ == fixbb){
			for (core::Size i = 1; i <= rounds_; ++i){
				packer.apply(pose);
				TR << (*scorefxn_)(pose) <<std::endl;
				mc->boltzmann(pose);
			}
		}
		else if( design_method_ == flxbb) {
			protocols::flxbb::FlxbbDesignOP flx( new protocols::flxbb::FlxbbDesign() );
			protocols::flxbb::DesignTaskOP des( new protocols::flxbb::DesignTask_Normal() );
			flx->set_movemap(movemap);
			core::pack::task::TaskFactory::const_iterator it;
			for (it = tf->begin(); it != tf->end(); ++it){
				des->add_task_operation(*it);
			}

			flx->add_design_task(des);
			flx->set_scorefxn_design(scorefxn_);
			flx->set_scorefxn_relax(scorefxn_);

			for (core::Size i =1; i <= rounds_; ++i){
				flx->apply(pose);
				TR << (*scorefxn_)(pose) <<std::endl;
				mc->boltzmann(pose);
			}


		}
		else if ( design_method_ == relaxed_design) {

			protocols::relax::FastRelaxOP rel( new protocols::relax::FastRelax(scorefxn_) );
			rel->set_movemap(movemap);

			//Optionally minimize bond length + angles - Test using cmd-line first!
			rel->set_task_factory(tf);
			for (core::Size i = 1; i<=rounds_; ++i){
				rel->apply(pose);
				TR << (*scorefxn_)(pose) <<std::endl;
				mc->boltzmann(pose);
			}

		}
		else if ( design_method_ == docked_design){
			TR << "Docked design not currently implemented. " << std::endl;
		}
		else{
			utility_exit_with_message("Antibody Design method not recognized. ");
		}


		mc->recover_low(pose);
		TR << "AntibodySeqDesignMover complete." << std::endl;
}


} //design
} //antibody
} //protocols
