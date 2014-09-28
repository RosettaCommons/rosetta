// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/AntibodySeqDesignTFCreator.cc
/// @brief Class for generating TaskFactories and TaskOperations for antibody design from a set of options
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/AntibodySeqDesignTFCreator.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/design/ResidueProbDesignOperation.hh>
#include <protocols/antibody/design/ConservativeDesignOperation.hh>
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

static thread_local basic::Tracer TR( "protocols.antibody.design.AntibodySeqDesignTFCreator" );

namespace protocols {
namespace antibody {
namespace design {

	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using namespace protocols::toolbox::task_operations;
	using namespace protocols::loops;

AntibodySeqDesignTFCreator::AntibodySeqDesignTFCreator( AntibodyInfoCOP ab_info, bool force_north_paper_db ):
	utility::pointer::ReferenceCount(),
	ab_info_( ab_info ),
	force_north_paper_db_(force_north_paper_db)
{


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
	assert( design_options.size() == 6 );

	setup_default_options();
	read_command_line_options();
	cdr_design_options_ = design_options;

}

AntibodySeqDesignTFCreator::~AntibodySeqDesignTFCreator() {}

void
AntibodySeqDesignTFCreator::setup_default_options(){
	cdr_design_options_.clear();
	for (core::Size i = 1; i <= 6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		CDRSeqDesignOptionsOP options( new CDRSeqDesignOptions( cdr ) );
		cdr_design_options_.push_back( options );
	}
	design_proline_ = true; //Proline conservation in clusters is very good but not 100 percent.
	design_antigen_ = false;
	design_framework_ = false;
}

void
AntibodySeqDesignTFCreator::read_command_line_options(){
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	neighbor_detection_dis( option [OptionKeys::antibody::design::neighbor_dis]() );
	set_probability_data_cutoff( option [OptionKeys::antibody::design::stats_cutoff]() );
	set_basic_design( option [OptionKeys::antibody::design::benchmark_basic_design]() );
	zero_prob_weight_ = option [ OptionKeys::antibody::design::sample_zero_probs_at]();
	design_antigen_ = option [ OptionKeys::antibody::design::design_antigen]();
	design_framework_ = option [ OptionKeys::antibody::design::design_framework]();
	design_framework_conservative_ = option [ OptionKeys::antibody::design::conservative_framework_design]();
	design_framework_conserved_res_ = false; // Place cmd-line option here if needed.

}

void
AntibodySeqDesignTFCreator::set_cdr_design_options( CDRNameEnum cdr, CDRSeqDesignOptionsCOP design_options ){
	cdr_design_options_[ cdr ] = design_options->clone();
}

void
AntibodySeqDesignTFCreator::set_cdr_design_options( const utility::vector1<CDRSeqDesignOptionsOP> design_options ){
	cdr_design_options_ = design_options;
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
AntibodySeqDesignTFCreator::set_zero_prob_weight_at( const core::Real weight ){
	zero_prob_weight_ = weight;
}

void
AntibodySeqDesignTFCreator::design_proline( const bool setting ) {
	design_proline_ = setting;
}

void
AntibodySeqDesignTFCreator::set_basic_design( const bool setting ){
	basic_design_ = setting;
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

std::map< core::Size, std::map< core::chemical::AA, core::Real > >
AntibodySeqDesignTFCreator::setup_probability_data( const core::pose::Pose& pose ){
	no_data_cdrs_.clear();
	AntibodyDatabaseManager manager = AntibodyDatabaseManager( ab_info_, force_north_paper_db_ );
	std::map< core::Size, std::map< core::chemical::AA, core::Real > > prob_set;
	vector1<CDRNameEnum> no_data_cdrs_ = manager.load_cdr_design_data( cdr_design_options_, pose, prob_set, prob_cutoff_ );

	return prob_set;
}

void
AntibodySeqDesignTFCreator::remove_conservative_design_residues_from_prob_set(vector1<core::Size> const & positions, std::map< core::Size, std::map< core::chemical::AA, core::Real > > & prob_set){
	for (core::Size i = 1; i<= positions.size(); ++i){
		std::map< core::Size, std::map<core::chemical::AA, core::Real > >::iterator it = prob_set.find(positions[ i ]);
		if ( it != prob_set.end() ){
			prob_set.erase( it );
			TR << "Removing "<< i << "  from probabilistic design" <<std::endl;
		}
	}
}

vector1<core::Size>
AntibodySeqDesignTFCreator::get_conservative_design_residues(const core::pose::Pose& pose){

	vector1<core::Size> conservative_positions;

	//If NO  or sparse data then default to conservative mutations.
	for (core::SSize i = CDRNameEnum_start; i <= CDRNameEnum_total; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		bool no_data_cdr = no_data_cdrs_.has_value(cdr);
		if ( ( cdr_design_options_[ cdr ]->design_strategy() ==  seq_design_conservative ) || no_data_cdr   ) {

			core::Size start = ab_info_->get_CDR_start( cdr, pose );
			core::Size end = ab_info_->get_CDR_end( cdr, pose );

			for ( core::Size res = start; res <= end; ++res ){
				conservative_positions.push_back( res );
				//TR << "Treating "<< res << " as conservative " << std::endl;
			}


		}
	}
	if ( design_framework_conservative_ ){
		for ( core::Size i = 1; i <= pose.total_residue(); ++i ){
			if ( ab_info_->get_region_of_residue(pose, i) ==  framework_region ){
				conservative_positions.push_back( i );
			}
		}
	}
	return conservative_positions;
}

void
AntibodySeqDesignTFCreator::disable_design_for_non_designing_cdrs(
	core::pack::task::TaskFactoryOP tf,
	const core::pose::Pose& pose) {
	for (core::Size i =1; i <= 6; ++i){
		if ( ! cdr_design_options_[ i ]->design() ){
			CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
			tf->push_back(protocols::antibody::design::disable_design_cdr(ab_info_, cdr, pose));
		}
	}
}

void
AntibodySeqDesignTFCreator::disable_proline_design(core::pack::task::TaskFactoryOP tf, const core::pose::Pose& pose){


	//One restrict op per CDR.  That way we can pop them off the TF  individually if we need to.
	core::pack::task::operation::RestrictResidueToRepackingOP restrict( new core::pack::task::operation::RestrictResidueToRepacking() );
	for (core::Size i = 1; i <= pose.total_residue(); ++i){
		if (pose.residue( i ).aa() == core::chemical::aa_pro){
			restrict->include_residue( i );
		}

	}
	tf->push_back( restrict );
}

protocols::loops::LoopsOP
AntibodySeqDesignTFCreator::get_design_cdr_loops(const core::pose::Pose& pose, core::Size stem_size /* 0 */ ) const{
	utility::vector1<bool> design_cdrs( 6, false );

	for (core::Size i = 1; i <= core::Size(CDRNameEnum_total); ++i){
		if (cdr_design_options_[ i ]->design()){
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
	LoopsOP all_cdrs = ab_info_->get_CDR_loops( pose, stem_size_ );
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

ResidueProbDesignOperationOP
AntibodySeqDesignTFCreator::generate_task_op_cdr_profile( core::pose::Pose const & pose ){

	ResidueProbDesignOperationOP prob_task( new ResidueProbDesignOperation() );
	std::map< core::Size, std::map< core::chemical::AA, core::Real > > prob_set = setup_probability_data( pose );
	vector1<core::Size> conservative_positions = get_conservative_design_residues( pose );

	remove_conservative_design_residues_from_prob_set( conservative_positions, prob_set );

	prob_task->set_aa_probability_set( prob_set );
	prob_task->set_keep_task_allowed_aas( false );
	prob_task->set_include_native_restype( true );
	prob_task->set_sample_zero_probs_at( zero_prob_weight_ ); //Control by cmd line flag. We do want to have some variability that isn't known.;
	return prob_task;
}

ConservativeDesignOperationOP
AntibodySeqDesignTFCreator::generate_task_op_cdr_conservative( core::pose::Pose const & pose ){
	vector1<core::Size> conservative_positions = get_conservative_design_residues( pose );
	ConservativeDesignOperationOP cons_task( new ConservativeDesignOperation() );
	if (!conservative_positions.empty()){
		TR << "Adding ConservativeDesignOp "<<std::endl;

		cons_task->limit_to_positions( conservative_positions );
		cons_task->include_native_aa( true );
		cons_task->add_to_allowed_aas( false );
	}
	else{
		TR << "No conservative positions to add!" << std::endl;
	}
	return cons_task;
}

core::pack::task::TaskFactoryOP
AntibodySeqDesignTFCreator::generate_tf_seq_design( const core::pose::Pose & pose){

	TaskFactoryOP tf( new TaskFactory() );


	//Setup Basic TaskOP
	tf->push_back(TaskOperationCOP( new InitializeFromCommandline() ));

	if (! design_antigen_  ){
		RestrictResidueToRepackingOP restrict_antigen = disable_design_region( ab_info_, pose, antigen_region );
		tf->push_back( restrict_antigen );
	}

	if ( ! design_framework_ ) {
		RestrictResidueToRepackingOP restrict_framework = disable_design_region( ab_info_, pose, framework_region );
		tf->push_back( restrict_framework );
	}

	if ( ! design_framework_conserved_res_ ) {
		RestrictResidueToRepackingOP restrict_conserved = disable_conserved_framework_positions( ab_info_, pose );
		tf->push_back( restrict_conserved );
	}

	tf->push_back(TaskOperationCOP( new operation::NoRepackDisulfides() ));

	//Optionally disable Proline design
	if (!design_proline_){
		disable_proline_design( tf, pose );
	}

	//For benchmarking
	if (basic_design_){
		return tf;
	}
	//Setup Prob TaskOp.
	TR << "Adding ResidueProbDesignOp " << std::endl;
	ResidueProbDesignOperationOP prob_task = generate_task_op_cdr_profile( pose );
	tf->push_back( prob_task );

	//Use conservative mutations for non-cluster positions + Optionally H3.
	ConservativeDesignOperationOP cons_task = generate_task_op_cdr_conservative( pose );
	tf->push_back( cons_task );

	return tf;
}

core::pack::task::TaskFactoryOP
AntibodySeqDesignTFCreator::generate_tf_seq_design_graft_design(
		const core::pose::Pose& pose,
		CDRNameEnum cdr,
		utility::vector1<bool> const & neighbor_cdr_min)
{

	assert( neighbor_cdr_min.size() == 6 );

	//Design CDR of interest, including neighbor CDRs, framework, +/or antigen.

	TaskFactoryOP tf( new TaskFactory() );

	//Make sure to add current cdr to design options.
	utility::vector1<CDRSeqDesignOptionsOP> local_options = cdr_design_options_;
	local_options[ cdr ]->design( true );

	//Create a vector of CDRs that should only min.
	utility::vector1<bool> only_min_cdrs(6, false);
	for ( core::Size i = 1; i <= 6; ++i ){
		if (neighbor_cdr_min[ i ] && (! local_options[ i ]->design()) ){
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
	restrict_design_to_set_cdrs->set_restrict_only_design_to_loops(true);
	tf->push_back(restrict_design_to_set_cdrs);

	//Restrict Design to only CDRs enabled in instructions - aka if we are not design L2, do not try and design L2.
	//Still pack L2 though
	for ( core::Size i = 1; i <= 6; ++i ){
		if (!cdr_design_options_[ i ]->design()){
			CDRNameEnum cdr_i = static_cast<CDRNameEnum>( i );
			RestrictResidueToRepackingOP restrict= disable_design_cdr( ab_info_, cdr_i, pose );
			tf->push_back( restrict );
		}
	}

	if (! design_antigen_ ){
		RestrictResidueToRepackingOP restrict_antigen = disable_design_region( ab_info_, pose, antigen_region );
		tf->push_back( restrict_antigen );
	}

	if ( ! design_framework_ ){
		RestrictResidueToRepackingOP restrict_framework = disable_design_region( ab_info_, pose, framework_region );
		tf->push_back( restrict_framework );
	}

	if ( ! design_framework_conserved_res_ ) {
		RestrictResidueToRepackingOP restrict_conserved = disable_conserved_framework_positions( ab_info_, pose );
		tf->push_back( restrict_conserved );
	}

	tf->push_back(TaskOperationCOP( new operation::NoRepackDisulfides() ));

	//Optionally disable Proline design
	if ( !design_proline_ ){
		disable_proline_design( tf, pose );
	}

	//For benchmarking
	if ( basic_design_ ){
		return tf;
	}

	//Now that we have the residues we are designing, we now control what we design:
	//Setup Prob TaskOp.
	ResidueProbDesignOperationOP prob_task = generate_task_op_cdr_profile( pose );
	tf->push_back( prob_task );

	//Use conservative mutations for non-cluster positions + Optionally H3.
	ConservativeDesignOperationOP cons_task = generate_task_op_cdr_conservative( pose );
	tf->push_back( cons_task );

	return tf;

}


}
}
}










