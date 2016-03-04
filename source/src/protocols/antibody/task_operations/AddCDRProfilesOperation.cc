// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocol/antibody/task_operations/AddCDRProfilesOperation.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/task_operations/AddCDRProfilesOperation.hh>
#include <protocols/antibody/task_operations/AddCDRProfilesOperationCreator.hh>

#include <protocols/antibody/task_operations/AddCDRProfileSetsOperation.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/design/AntibodyDesignEnumManager.hh>

#include <protocols/toolbox/task_operations/ResidueProbDesignOperation.hh>
#include <protocols/toolbox/task_operations/ConservativeDesignOperation.hh>

#include <core/pack/task/operation/TaskOperations.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.task_operations.AddCDRProfilesOperation");

namespace protocols {
namespace antibody {
namespace task_operations {

using namespace core::pack::task::operation;
using namespace protocols::antibody::design;
using namespace protocols::antibody::clusters;
using namespace protocols::toolbox::task_operations;
using namespace basic::options;
using namespace basic::options::OptionKeys;

AddCDRProfilesOperation::AddCDRProfilesOperation():
	TaskOperation(),
	ab_info_(/* NULL */),
	profile_sets_task_(/* NULL */)
{
	set_defaults();
}

AddCDRProfilesOperation::AddCDRProfilesOperation(AntibodyInfoCOP ab_info):
	TaskOperation(),
	ab_info_(ab_info),
	profile_sets_task_(/* NULL */)
{
	set_defaults();
}

AddCDRProfilesOperation::AddCDRProfilesOperation(AntibodyInfoCOP ab_info, const utility::vector1<bool>& cdrs):
	TaskOperation(),
	ab_info_(ab_info),
	profile_sets_task_(/* NULL */)
{
	set_defaults();
	set_cdrs(cdrs);
}

void
AddCDRProfilesOperation::set_defaults(){

	seq_design_options_.clear();

	pre_loaded_data_ = false;

	//Profile Options
	picking_rounds_ = 1;
	keep_task_allowed_aas_ = false;
	include_native_restype_ = true;
	force_north_paper_db_ = false;
	use_outliers_ = false;
	stats_cutoff_ = 10;
	zero_prob_sample_ = 0;
	cons_design_data_source_ = "blosum62";
	pre_loaded_data_ = false;
	prob_set_.clear();
	seq_design_options_.clear();

	AntibodyEnumManager manager = AntibodyEnumManager();
	std::string numbering_scheme = option [OptionKeys::antibody::numbering_scheme]();
	numbering_scheme_ = manager.numbering_scheme_string_to_enum(numbering_scheme);

	for ( core::Size i = 1; i <=6; ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		CDRSeqDesignOptionsOP opt = CDRSeqDesignOptionsOP(new CDRSeqDesignOptions(cdr));
		seq_design_options_.push_back(opt);
	}

	cons_task_ = ConservativeDesignOperationOP( new ConservativeDesignOperation());

}



AddCDRProfilesOperation::AddCDRProfilesOperation(AddCDRProfilesOperation const & src):
	TaskOperation(src),
	ab_info_(src.ab_info_),
	seq_design_options_(src.seq_design_options_),
	cons_task_(src.cons_task_),
	profile_sets_task_(src.profile_sets_task_),
	picking_rounds_(src.picking_rounds_),
	keep_task_allowed_aas_(src.keep_task_allowed_aas_),
	include_native_restype_(src.include_native_restype_),
	force_north_paper_db_(src.force_north_paper_db_),
	use_outliers_(src.use_outliers_),
	stats_cutoff_(src.stats_cutoff_),
	zero_prob_sample_(src.zero_prob_sample_),
	cons_design_data_source_(src.cons_design_data_source_),
	pre_loaded_data_(src.pre_loaded_data_),
	prob_set_(src.prob_set_),
	no_profile_data_cdrs_(src.no_profile_data_cdrs_),
	no_profile_sets_data_cdrs_(src.no_profile_sets_data_cdrs_),
	numbering_scheme_(src.numbering_scheme_)


{

}

TaskOperationOP
AddCDRProfilesOperation::clone() const{
	return TaskOperationOP( new AddCDRProfilesOperation( *this ));
}

AddCDRProfilesOperation::~AddCDRProfilesOperation(){}

void
AddCDRProfilesOperation::parse_tag(utility::tag::TagCOP tag, basic::datacache::DataMap&){
	if ( tag->hasOption("cdrs") ) {
		TR << "Setting CDRs from settings" << std::endl;
		set_cdrs(get_cdr_bool_from_tag(tag, "cdrs"));
	}

	AntibodyDesignEnumManager manager = AntibodyDesignEnumManager();
	set_fallback_strategy(manager.seq_design_strategy_string_to_enum(tag->getOption< std::string >("fallback_strategy", "seq_design_conservative")));
	keep_task_allowed_aas_ = tag->getOption< bool>("add_to_current", keep_task_allowed_aas_);
	include_native_restype_ = tag->getOption< bool>("include_native_restype", include_native_restype_);
	picking_rounds_ = tag->getOption< core::Size >("picking_rounds", picking_rounds_);
	force_north_paper_db_ = tag->getOption< bool >("force_north_paper_db", force_north_paper_db_);
	use_outliers_ = tag->getOption< bool >("use_outliers", use_outliers_);
	stats_cutoff_ = tag->getOption< core::Size >("stats_cutoff", stats_cutoff_);
	zero_prob_sample_ = tag->getOption< core::Real >("sample_zero_probs_at", zero_prob_sample_);
	set_cons_design_data_source( tag->getOption< std::string >("cons_design_data_source", cons_design_data_source_) );

	if ( tag->hasOption("numbering_scheme") ) {


		AntibodyEnumManager manager = AntibodyEnumManager();
		numbering_scheme_ = manager.numbering_scheme_string_to_enum(tag->getOption<std::string>("numbering_scheme"));

	}
	if ( tag->hasOption("cdr_definition") && tag->getOption<std::string>("cdr_definition") != "North" ) {
		TR <<"This operation only works with the North CDR definition." <<std::endl;
	}



}

void
AddCDRProfilesOperation::set_cdr_only(CDRNameEnum cdr){
	for ( core::Size i = 1; i <= 6; ++i ) {
		CDRNameEnum current_cdr = static_cast<CDRNameEnum>( i );

		if ( current_cdr == cdr ) {
			seq_design_options_[ i ]->design(true);
		} else {
			seq_design_options_[ i ]->design(false);
		}
	}
}

void
AddCDRProfilesOperation::set_cdrs(const utility::vector1<bool>& cdrs) {
	for ( core::Size i = 1; i <= 6; ++i ) {
		if ( cdrs[ i ] ) {
			seq_design_options_[ i ]->design(true);
		} else {
			seq_design_options_[ i ]->design(false);
		}
	}
}

void
AddCDRProfilesOperation::set_design_options(design::AntibodyCDRSeqDesignOptions seq_design_options){
	seq_design_options_ = seq_design_options;
}

void
AddCDRProfilesOperation::set_add_to_current(bool add_to_current){
	keep_task_allowed_aas_ = add_to_current;
}

void
AddCDRProfilesOperation::set_fallback_strategy(SeqDesignStrategyEnum fallback_strategy){
	for ( core::Size i = 1; i <= 6; ++i ) {
		seq_design_options_[ i ]->fallback_strategy(fallback_strategy);
	}
}

void
AddCDRProfilesOperation::set_primary_strategy(SeqDesignStrategyEnum primary_strategy){
	for ( core::Size i = 1; i <= 6; ++i ) {
		seq_design_options_[ i ]->design_strategy(primary_strategy);
	}
}

void
AddCDRProfilesOperation::set_picking_rounds(core::Size rounds){
	picking_rounds_ = rounds;
}

void
AddCDRProfilesOperation::set_include_native_type(bool use_native) {
	include_native_restype_ = use_native;
}

void
AddCDRProfilesOperation::set_force_north_paper_db(bool force_north_db){
	force_north_paper_db_ = force_north_db;
}

void
AddCDRProfilesOperation::set_stats_cutoff(core::Size stats_cutoff) {
	stats_cutoff_ = stats_cutoff;
}

void
AddCDRProfilesOperation::set_use_outliers(bool use_outliers){
	use_outliers_ = use_outliers;
}

void
AddCDRProfilesOperation::set_sample_zero_probs_at(core::Real zero_prob_sample){
	zero_prob_sample_ = zero_prob_sample;
}

void
AddCDRProfilesOperation::set_cons_design_data_source(std::string data_source){

	//Re-load conservative data
	if ( cons_design_data_source_ != data_source ) {
		cons_task_->set_data_source(data_source);
	}
	cons_design_data_source_ = data_source;
}

utility::vector1<bool>
AddCDRProfilesOperation::get_profile_and_design_cdrs() const {
	utility::vector1<bool> profiles(6, false);
	for ( core::Size i = 1; i <= 6; ++i ) {
		if ( ! seq_design_options_[ i ]->design() ) continue;
		SeqDesignStrategyEnum strat = seq_design_options_[ i ]->design_strategy();
		if ( strat == seq_design_profiles || strat == seq_design_profile_sets_combined ) {
			profiles[ i ] = true;
		}
	}
	return profiles;
}

utility::vector1<bool>
AddCDRProfilesOperation::get_profile_set_and_design_cdrs() const {
	utility::vector1<bool> profiles(6, false);
	for ( core::Size i = 1; i <= 6; ++i ) {
		if ( ! seq_design_options_[ i ]->design() ) continue;

		SeqDesignStrategyEnum strat = seq_design_options_[ i ]->design_strategy();
		if ( strat == seq_design_profile_sets || strat == seq_design_profile_sets_combined ) {
			profiles[ i ] = true;
		}
	}
	return profiles;
}

utility::vector1<bool>
AddCDRProfilesOperation::get_design_cdrs() const{
	utility::vector1<bool> designing(6, false);
	for ( core::Size i = 1; i <= 6; ++i ) {
		if ( seq_design_options_[ i ]->design() ) {
			designing[ i ] = true;
		}
	}
	return designing;
}

void
AddCDRProfilesOperation::pre_load_data(const core::pose::Pose& pose){
	utility::vector1<bool> profile_set_cdrs = this->get_profile_set_and_design_cdrs();
	utility::vector1<bool> profile_cdrs = this->get_profile_and_design_cdrs();

	core::Size n_profile_set_cdrs = count(profile_set_cdrs.begin(), profile_set_cdrs.end(), true);
	core::Size n_profile_cdrs = count(profile_cdrs.begin(), profile_cdrs.end(), true);

	if ( ! ab_info_ ) {
		ab_info_ = AntibodyInfoOP(new AntibodyInfo(pose, numbering_scheme_, North));
	}

	if ( n_profile_set_cdrs > 0 ) {

		profile_sets_task_ = AddCDRProfileSetsOperationOP( new AddCDRProfileSetsOperation(ab_info_));
		profile_sets_task_->set_cdrs(profile_set_cdrs);
		profile_sets_task_->set_force_north_paper_db(force_north_paper_db_);
		profile_sets_task_->set_use_outliers(use_outliers_);
		profile_sets_task_->set_cutoff(stats_cutoff_);

		no_profile_sets_data_cdrs_ = profile_sets_task_->pre_load_data(pose);

	}

	if ( n_profile_cdrs > 0 ) {

		prob_set_ = get_cluster_profile_probability_data(
			ab_info_,
			pose,
			profile_cdrs,
			no_profile_data_cdrs_,
			stats_cutoff_,
			use_outliers_,
			force_north_paper_db_);

	}

	pre_loaded_data_ = true;
	TR << "Data pre-loaded." << std::endl;
}

void
AddCDRProfilesOperation::apply(const core::pose::Pose& pose, core::pack::task::PackerTask& task) const{

	//This is due to const apply and no pose in parse_my_tag.
	AntibodyInfoOP local_ab_info;
	if ( ! ab_info_ ) {
		local_ab_info = AntibodyInfoOP(new AntibodyInfo(pose, numbering_scheme_, North));
	} else {
		local_ab_info = ab_info_->clone();
	}


	utility::vector1<bool> no_profile_set_cdrs;
	utility::vector1<bool> no_profile_cdrs;

	utility::vector1<bool> profile_set_cdrs = this->get_profile_set_and_design_cdrs();
	utility::vector1<bool> profile_cdrs = this->get_profile_and_design_cdrs();

	core::Size n_profile_set_cdrs = count(profile_set_cdrs.begin(), profile_set_cdrs.end(), true);
	core::Size n_profile_cdrs = count(profile_cdrs.begin(), profile_cdrs.end(), true);


	/// Load, setup, and run ProfileSets operation.
	if ( n_profile_set_cdrs > 0 ) {
		AddCDRProfileSetsOperationOP profile_sets_task; //yay const apply
		if ( pre_loaded_data_ ) {
			no_profile_set_cdrs = no_profile_sets_data_cdrs_;
			profile_sets_task = AddCDRProfileSetsOperationOP( new AddCDRProfileSetsOperation( *profile_sets_task_ ) );
		} else {
			profile_sets_task = AddCDRProfileSetsOperationOP( new AddCDRProfileSetsOperation( local_ab_info));

			profile_sets_task->set_cdrs(profile_set_cdrs);
			profile_sets_task->set_force_north_paper_db(force_north_paper_db_);
			profile_sets_task->set_use_outliers(use_outliers_);
			profile_sets_task->set_cutoff(stats_cutoff_);

			no_profile_set_cdrs = profile_sets_task->pre_load_data(pose);
		}
		TR << "applying profile sets op" << std::endl;

		profile_sets_task->set_picking_rounds(picking_rounds_);
		profile_sets_task->set_include_native_type(include_native_restype_);
		profile_sets_task->set_add_to_current(keep_task_allowed_aas_);
		profile_sets_task->apply(pose, task);
	}

	if ( n_profile_cdrs > 0 ) {
		std::map< core::Size, std::map< core::chemical::AA, core::Real > > prob_set;

		if ( pre_loaded_data_ ) {
			prob_set = prob_set_;
			no_profile_cdrs = no_profile_data_cdrs_;
		} else {
			prob_set = get_cluster_profile_probability_data(
				local_ab_info,
				pose,
				profile_cdrs,
				no_profile_cdrs,
				stats_cutoff_,
				use_outliers_,
				force_north_paper_db_);
		}

		ResidueProbDesignOperation prob_task = ResidueProbDesignOperation();

		prob_task.set_picking_rounds(picking_rounds_);
		prob_task.set_aa_probability_set( prob_set );

		if ( n_profile_set_cdrs > 0 ) {
			prob_task.set_keep_task_allowed_aas(true);
		} else {
			prob_task.set_keep_task_allowed_aas( keep_task_allowed_aas_);
		}
		prob_task.set_include_native_restype( include_native_restype_ );
		prob_task.set_sample_zero_probs_at( zero_prob_sample_ );
		TR << "applying prob task op"<<std::endl;
		prob_task.apply(pose, task);
	}

	utility::vector1< bool > no_data_cdrs(6, false);

	if ( n_profile_cdrs > 0 ) {
		no_data_cdrs = no_profile_cdrs;
	} else if ( n_profile_set_cdrs > 0 ) {
		no_data_cdrs = no_profile_set_cdrs;
	}

	/// Add the conservative design op.
	cons_task_->add_to_allowed_aas(keep_task_allowed_aas_);
	cons_task_->include_native_aa(include_native_restype_);
	core::Size cons_task_residues = 0;

	for ( core::Size i = 1; i <= core::Size(local_ab_info->get_total_num_CDRs()); ++i ) {
		CDRNameEnum cdr = static_cast< CDRNameEnum >( i );
		if ( (no_data_cdrs[ i ]  && seq_design_options_[ i ]->fallback_strategy() == seq_design_conservative) || (seq_design_options_[ i ]->design() && seq_design_options_[ i ]->design_strategy() == seq_design_conservative) ) {

			TR << "Using conservative op for " << local_ab_info->get_CDR_name(cdr) << std::endl;
			for ( core::Size resnum = local_ab_info->get_CDR_start(cdr, pose, North); resnum <= local_ab_info->get_CDR_end(cdr, pose, North); ++resnum ) {
				cons_task_->include_residue( resnum );
				cons_task_residues += 1;
			}
		}
	}

	if ( cons_task_residues > 0 ) {
		if ( has_native_sequence( pose ) ) {
			TR << "Using original bb sequence for conservative design." << std::endl;
			std::string native_seq = get_native_sequence( pose );
			cons_task_->set_native_sequence( native_seq );
		} else {
			cons_task_->use_pose_sequence_as_native( pose );
		}

		cons_task_->apply(pose, task);

	}



	/// Disable design if seq design is set as none.
	RestrictResidueToRepacking restrict = RestrictResidueToRepacking();
	core::Size restrict_residues = 0;
	//TR << "No data cdrs: " << utility::to_string(no_data_cdrs) << std::endl;
	for ( core::Size i = 1; i <= core::Size(local_ab_info->get_total_num_CDRs()); ++i ) {
		CDRNameEnum cdr = static_cast< CDRNameEnum >( i );
		if ( (no_data_cdrs[ cdr ] && seq_design_options_[ i ]->fallback_strategy() == seq_design_none) || (seq_design_options_[ i ]->design() && seq_design_options_[ i ]->design_strategy() == seq_design_none) ) {
			TR << "Disabling design for " << local_ab_info->get_CDR_name(cdr) << std::endl;

			for ( core::Size resnum = local_ab_info->get_CDR_start(cdr, pose, North); resnum <= local_ab_info->get_CDR_end(cdr, pose, North); ++resnum ) {
				restrict.include_residue( resnum );
				restrict_residues += 1;
			}
		}

	}
	if ( restrict_residues > 0 ) {
		restrict.apply(pose, task);
	}
}

core::pack::task::operation::TaskOperationOP
AddCDRProfilesOperationCreator::create_task_operation() const
{
	return core::pack::task::operation::TaskOperationOP( new AddCDRProfilesOperation );
}

} //task_operations
} //antibody
} //protocols











