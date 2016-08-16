// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/AntibodySeqDesignTFCreator.cc
/// @brief Class for generating TaskFactories and TaskOperations for antibody design from a set of options
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/AntibodySeqDesignTFCreator.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/task_operations/AddCDRProfilesOperation.hh>
#include <protocols/toolbox/task_operations/ResidueProbDesignOperation.hh>
#include <protocols/toolbox/task_operations/ConservativeDesignOperation.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>

#include <protocols/loops/Loops.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh>

#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.design.AntibodySeqDesignTFCreator" );

namespace protocols {
namespace antibody {
namespace design {

using namespace core::pack::task;
using namespace core::pack::task::operation;
using namespace protocols::toolbox::task_operations;
using namespace protocols::loops;
using namespace protocols::antibody::task_operations;

AntibodySeqDesignTFCreator::AntibodySeqDesignTFCreator( AntibodyInfoCOP ab_info, bool force_north_paper_db ):
	utility::pointer::ReferenceCount(),
	ab_info_( ab_info ),
	force_north_paper_db_(force_north_paper_db)
{
	stem_size_ = 2;
	setup_default_options();
	read_command_line_options();
}

AntibodySeqDesignTFCreator::AntibodySeqDesignTFCreator( AntibodyInfoCOP ab_info,
	const utility::vector1<CDRSeqDesignOptionsOP> design_options,
	bool force_north_paper_db, core::Size stem_size):
	utility::pointer::ReferenceCount(),
	ab_info_( ab_info ),
	stem_size_(stem_size),
	force_north_paper_db_(force_north_paper_db)

{
	assert( design_options.size() == 8 );

	setup_default_options();
	read_command_line_options();
	cdr_design_options_ = design_options;

}

AntibodySeqDesignTFCreator::~AntibodySeqDesignTFCreator() {}

void
AntibodySeqDesignTFCreator::setup_default_options(){
	cdr_design_options_.clear();
	for ( core::Size i = 1; i <= 8; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		CDRSeqDesignOptionsOP options( new CDRSeqDesignOptions( cdr ) );
		cdr_design_options_.push_back( options );
	}
	design_proline_ = true; //Proline conservation in clusters is very good but not 100 percent.
	design_antigen_ = false;
	design_framework_ = false;

	for ( core::Size i = 1; i <= 8; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		if ( ! ab_info_->has_CDR( cdr ) ) {
			cdr_design_options_[ cdr ]->design( false);
		}
	}

	no_data_cdrs_.clear();
	no_data_cdrs_.resize(8, false);
}

void
AntibodySeqDesignTFCreator::read_command_line_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	neighbor_detection_dis( option [OptionKeys::antibody::design::neighbor_dis]() );
	set_probability_data_cutoff( option [OptionKeys::antibody::design::seq_design_stats_cutoff]() );
	zero_prob_weight_ = option [ OptionKeys::antibody::design::sample_zero_probs_at]();
	design_antigen_ = option [ OptionKeys::antibody::design::design_antigen]();
	design_framework_ = option [ OptionKeys::antibody::design::design_framework]();
	design_framework_conservative_ = option [ OptionKeys::antibody::design::conservative_framework_design]();
	design_framework_conserved_res_ = false; // Place cmd-line option here if needed.
	design_h3_stem_ = option[ OptionKeys::antibody::design::design_H3_stem ]();
	profile_picking_rounds_ = option[ OptionKeys::antibody::design::seq_design_profile_samples]();
	use_outliers_ = option[OptionKeys::antibody::design::use_outliers]();

}

void
AntibodySeqDesignTFCreator::set_cdr_design_options( CDRNameEnum cdr, CDRSeqDesignOptionsCOP design_options ){
	cdr_design_options_[ cdr ] = design_options->clone();
	if ( ! ab_info_->has_CDR( cdr ) ) {
		cdr_design_options_[ cdr ]->design( false);
	}
}

void
AntibodySeqDesignTFCreator::set_cdr_design_options( const utility::vector1<CDRSeqDesignOptionsOP> design_options ){
	cdr_design_options_ = design_options;
	for ( core::Size i = 1; i <= 8; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		if ( ! ab_info_->has_CDR( cdr ) ) {
			cdr_design_options_[ cdr ]->design( false);
		}
	}
}

utility::vector1<CDRSeqDesignOptionsOP>
AntibodySeqDesignTFCreator::get_cdr_design_options(){
	return cdr_design_options_;
}

CDRSeqDesignOptionsOP
AntibodySeqDesignTFCreator::get_cdr_design_options(CDRNameEnum cdr) {
	return cdr_design_options_[cdr];
}

void
AntibodySeqDesignTFCreator::neighbor_detection_dis( core::Real const neighbor_distance ){
	neighbor_dis_ = neighbor_distance;
}

void
AntibodySeqDesignTFCreator::set_probability_data_cutoff( core::Size const cutoff ){
	prob_cutoff_ = cutoff;
}

void
AntibodySeqDesignTFCreator::set_design_H3_stem(bool design_H3_stem) {
	design_h3_stem_ = design_H3_stem;
}

void
AntibodySeqDesignTFCreator::set_zero_prob_weight_at( const core::Real weight ){
	zero_prob_weight_ = weight;
}

void
AntibodySeqDesignTFCreator::design_proline( const bool setting ) {
	design_proline_ = setting;
}


void
AntibodySeqDesignTFCreator::design_antigen(bool antigen_design) {
	design_antigen_ = antigen_design;
}

void
AntibodySeqDesignTFCreator::design_framework(bool framework_design) {
	design_framework_ = framework_design;
}

void
AntibodySeqDesignTFCreator::set_design_framework_conservative(bool design_framework_conservative) {
	design_framework_conservative_ = design_framework_conservative;
}

void
AntibodySeqDesignTFCreator::set_design_framework_conserved_res(bool design_framework_conserved_res) {
	design_framework_conserved_res_ = design_framework_conserved_res;
}

core::pack::task::TaskFactoryOP
AntibodySeqDesignTFCreator::generate_tf_seq_design( const core::pose::Pose & pose, bool disable_non_designing_cdrs /* False */){

	TaskFactoryOP tf( new TaskFactory() );

	/*
	//Setup Basic TaskOP
	tf->push_back(TaskOperationCOP( new InitializeFromCommandline() ));

	add_extra_restrict_operations(tf, pose);

	//Setup Prob TaskOp.
	TR << "Adding ResidueProbDesignOp " << std::endl;
	ResidueProbDesignOperationOP prob_task = generate_task_op_cdr_profile( pose );
	tf->push_back( prob_task );

	//Use conservative mutations for non-cluster positions + Optionally H3.
	ConservativeDesignOperationOP cons_task = generate_task_op_cdr_conservative( pose );

	*/

	add_extra_restrict_operations(tf, pose);

	AddCDRProfilesOperationOP profile_strategy_task = this->generate_task_op_cdr_profiles(pose);
	tf->push_back( profile_strategy_task );

	if ( design_framework_conservative_ ) {
		ConservativeDesignOperationOP cons_task = get_framework_conservative_op(pose);
		tf->push_back(cons_task);
	}

	if ( disable_non_designing_cdrs ) {
		disable_design_for_non_designing_cdrs(tf, pose);
	}


	//disable_design_for_no_fallback_cdrs(tf, pose);

	return tf;
}

core::pack::task::TaskFactoryOP
AntibodySeqDesignTFCreator::generate_tf_seq_design_graft_design(
	const core::pose::Pose& pose,
	CDRNameEnum cdr,
	utility::vector1<bool> const & neighbor_cdr_min)
{

	assert( neighbor_cdr_min.size() == 8 );

	//Design CDR of interest, including neighbor CDRs, framework, +/or antigen.

	TaskFactoryOP tf( new TaskFactory() );

	//Make sure to add current cdr to design options.
	utility::vector1<CDRSeqDesignOptionsOP> local_options = cdr_design_options_;
	local_options[ cdr ]->design( true );

	//Create a vector of CDRs that should only min.
	utility::vector1<bool> only_min_cdrs(core::Size(CDRNameEnum_proto_total), false);
	for ( core::Size i = 1; i <= core::Size(CDRNameEnum_proto_total); ++i ) {
		CDRNameEnum cdr_i = static_cast<CDRNameEnum>( i );
		if ( neighbor_cdr_min[ i ] && (! local_options[ i ]->design()) && ab_info_->has_CDR( cdr_i ) ) {
			only_min_cdrs[ i ] = true;
		}
	}


	//Setup Basic TaskOP
	tf->push_back( TaskOperationCOP( new InitializeFromCommandline() ) );

	//Limit packing and design to only CDRs we are minimizing.
	RestrictToLoopsAndNeighborsOP restrict_to_min_cdrs = this->generate_task_op_cdr_design(pose, neighbor_cdr_min, true);
	tf->push_back(restrict_to_min_cdrs);

	//Limit design to only CDRs in options and neighbors
	RestrictToLoopsAndNeighborsOP restrict_design_to_set_cdrs = this->generate_task_op_cdr_design(pose, true);
	restrict_design_to_set_cdrs->set_restrict_only_design_to_loops(true); //Don't worry about the wrongly named option.  It will restrict design to loops and neighbors.
	tf->push_back(restrict_design_to_set_cdrs);

	//Restrict Design to only CDRs enabled in instructions - aka if we are not design L2, do not try and design L2.
	//Still pack L2 though
	for ( core::Size i = 1; i <= core::Size( ab_info_->get_total_num_CDRs( true /* Include DE loops */) ); ++i ) {
		if ( !cdr_design_options_[ i ]->design() ) {
			CDRNameEnum cdr_i = static_cast<CDRNameEnum>( i );
			//TR << "Disabling "<< ab_info_->get_CDR_name( cdr_i ) << std::endl;
			RestrictResidueToRepackingOP restrict= disable_design_cdr( ab_info_, cdr_i, pose );
			tf->push_back( restrict );
		}
	}

	add_extra_restrict_operations(tf, pose);

	//Now that we have the residues we are designing, we now control what we design:
	AddCDRProfilesOperationOP profile_strategy_task = this->generate_task_op_cdr_profiles(pose);
	tf->push_back( profile_strategy_task );

	if ( design_framework_conservative_ ) {
		ConservativeDesignOperationOP cons_task = get_framework_conservative_op(pose);
		tf->push_back(cons_task);
	}

	/*
	//Setup Prob TaskOp.
	ResidueProbDesignOperationOP prob_task = generate_task_op_cdr_profile( pose );
	tf->push_back( prob_task );

	//Use conservative mutations for non-cluster positions + Optionally H3.
	ConservativeDesignOperationOP cons_task = generate_task_op_cdr_conservative( pose );
	tf->push_back( cons_task );

	disable_design_for_no_fallback_cdrs(tf, pose);
	*/

	return tf;

}

AddCDRProfilesOperationOP
AntibodySeqDesignTFCreator::generate_task_op_cdr_profiles(core::pose::Pose const & pose){

	AddCDRProfilesOperationOP profile_op = AddCDRProfilesOperationOP(new AddCDRProfilesOperation(ab_info_));
	profile_op->set_use_outliers(use_outliers_);
	profile_op->set_stats_cutoff(prob_cutoff_);
	profile_op->set_picking_rounds(profile_picking_rounds_);
	profile_op->set_force_north_paper_db(force_north_paper_db_);
	profile_op->set_sample_zero_probs_at(zero_prob_weight_);
	profile_op->set_design_options(cdr_design_options_);
	profile_op->pre_load_data(pose);

	return profile_op;
}

ConservativeDesignOperationOP
AntibodySeqDesignTFCreator::get_framework_conservative_op(const core::pose::Pose& pose){

	utility::vector1<core::Size> conservative_positions;

	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( ab_info_->get_region_of_residue(pose, i, false /* DE as framework */) ==  framework_region ) {
			conservative_positions.push_back( i );
		}
	}

	ConservativeDesignOperationOP cons_task( new ConservativeDesignOperation() );
	if ( !conservative_positions.empty() ) {
		TR << "Adding ConservativeDesignOp for Framework Residues"<<std::endl;

		cons_task->limit_to_positions( conservative_positions );
		cons_task->include_native_aa( true );
		cons_task->add_to_allowed_aas( false );
	}

	if ( has_native_sequence( pose ) ) {
		TR <<"Using original bb sequence for conservative design" << std::endl;
		std::string native_seq = get_native_sequence( pose );
		cons_task->set_native_sequence( native_seq );
	}

	return cons_task;
}

void
AntibodySeqDesignTFCreator::disable_design_for_non_designing_cdrs(
	core::pack::task::TaskFactoryOP tf,
	const core::pose::Pose& pose) {
	for ( core::Size i =1; i <= core::Size( ab_info_->get_total_num_CDRs( true /* include DE loops */) ); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );


		if ( ! cdr_design_options_[ i ]->design() ) {
			//TR << "Disabling2 " << ab_info_->get_CDR_name( cdr ) << std::endl;
			tf->push_back(protocols::antibody::design::disable_design_cdr(ab_info_, cdr, pose));
		}
	}
}

void
AntibodySeqDesignTFCreator::disable_proline_design(core::pack::task::TaskFactoryOP tf, const core::pose::Pose& pose) const{


	//One restrict op per CDR.  That way we can pop them off the TF  individually if we need to.
	core::pack::task::operation::RestrictResidueToRepackingOP restrict( new core::pack::task::operation::RestrictResidueToRepacking() );
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( pose.residue( i ).aa() == core::chemical::aa_pro ) {
			restrict->include_residue( i );
		}

	}
	tf->push_back( restrict );
}

protocols::loops::LoopsOP
AntibodySeqDesignTFCreator::get_design_cdr_loops(const core::pose::Pose& pose, core::Size stem_size /* 0 */ ) const{
	utility::vector1<bool> design_cdrs( 8, false );

	for ( core::Size i = 1; i <= CDRNameEnum_proto_total; ++i ) {
		if ( cdr_design_options_[ i ]->design() ) {
			CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
			if ( TR.Debug.visible() ) {
				TR.Debug << "Design on: " << ab_info_->get_CDR_name( cdr ) << std::endl;
			}
			design_cdrs[ i ] = true;
		}
	}
	return protocols::antibody::get_cdr_loops( ab_info_, pose, design_cdrs, stem_size );
}

protocols::loops::LoopsOP
AntibodySeqDesignTFCreator::get_design_cdr_loops_with_stem( const core::pose::Pose& pose ) const {
	return get_design_cdr_loops( pose, stem_size_ );
}

protocols::toolbox::task_operations::RestrictToLoopsAndNeighborsOP
AntibodySeqDesignTFCreator::generate_task_op_cdr_design( core::pose::Pose const & pose, bool design_neighbors /*true*/) const{

	protocols::loops::LoopsOP cdr_loops = get_design_cdr_loops( pose );
	return this->get_general_loop_task_op(cdr_loops, design_neighbors);
}

protocols::toolbox::task_operations::RestrictToLoopsAndNeighborsOP
AntibodySeqDesignTFCreator::generate_task_op_cdr_design(const core::pose::Pose& pose, utility::vector1<bool> cdrs, bool design_neighbors /* true */) const {

	LoopsOP cdr_loops = ab_info_->get_CDR_loops( pose, cdrs, stem_size_ );
	return this->get_general_loop_task_op(cdr_loops, design_neighbors);
}

protocols::toolbox::task_operations::RestrictToLoopsAndNeighborsOP
AntibodySeqDesignTFCreator::generate_task_op_all_cdr_design( const core::pose::Pose& pose, bool design_neighbors /* true */ ) const {
	///Limit Packing and Design to CDR loops and neighbors

	utility::vector1< bool > cdrs(8, true);
	for ( core::Size i = 1; i <= 8; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		if ( ! ab_info_->has_CDR( cdr ) ) {
			cdrs[ i ] = false;
		}
	}
	LoopsOP all_cdrs = ab_info_->get_CDR_loops( pose, cdrs, stem_size_ );
	return this->get_general_loop_task_op(all_cdrs, design_neighbors);
}

protocols::toolbox::task_operations::RestrictToLoopsAndNeighborsOP
AntibodySeqDesignTFCreator::get_general_loop_task_op(LoopsOP loops, bool design_neighbors /* true */) const {

	RestrictToLoopsAndNeighborsOP loop_task( new RestrictToLoopsAndNeighbors() );
	loop_task->set_loops( loops );
	loop_task->set_design_loop( true );
	loop_task->set_include_neighbors( true );
	loop_task->set_cutoff_distance( neighbor_dis_ );
	loop_task->set_design_neighbors( design_neighbors );
	return loop_task;
}


void
AntibodySeqDesignTFCreator::add_extra_restrict_operations(core::pack::task::TaskFactoryOP tf, const core::pose::Pose & pose) const{

	if ( ! design_antigen_ ) {
		RestrictResidueToRepackingOP restrict_antigen = disable_design_region( ab_info_, pose, antigen_region );
		tf->push_back( restrict_antigen );
	}

	if ( ! design_framework_ ) {
		RestrictResidueToRepackingOP restrict_framework = disable_design_region( ab_info_, pose, framework_region, false /* cdr4 as framework */ );
		tf->push_back( restrict_framework );
	}

	if ( ! design_framework_conserved_res_ ) {
		RestrictResidueToRepackingOP restrict_conserved = disable_conserved_framework_positions( ab_info_, pose );
		tf->push_back( restrict_conserved );
	}

	if ( ! design_h3_stem_ ) {
		RestrictResidueToRepackingOP restrict_h3_stem = disable_h3_stem_positions(ab_info_, pose);
		tf->push_back( restrict_h3_stem );
	}

	tf->push_back(TaskOperationCOP( new operation::NoRepackDisulfides() ));

	//Optionally disable Proline design
	if ( !design_proline_ ) {
		disable_proline_design( tf, pose );
	}

}


}
}
}
