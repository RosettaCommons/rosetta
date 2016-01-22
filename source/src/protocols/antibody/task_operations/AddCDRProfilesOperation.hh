// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/task_operations/AddCDRProfilesOperation.hh
/// @brief Add CDR Cluster-based profile ResidueProbDesignOperation to sample
///  within a CDR cluster profile
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_task_operations_AddCDRProfilesOperation_hh
#define INCLUDED_protocols_antibody_task_operations_AddCDRProfilesOperation_hh

#include <protocols/antibody/task_operations/AddCDRProfilesOperation.fwd.hh>
#include <protocols/antibody/task_operations/AddCDRProfileSetsOperation.fwd.hh>

#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.fwd.hh>
#include <protocols/toolbox/task_operations/ConservativeDesignOperation.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperationCreator.hh>
#include <core/chemical/AA.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <map>

namespace protocols {
namespace antibody {
namespace task_operations {


/// @brief Add Cluster-based CDR Profiles as the task operation for the set of CDRs by default.
/// See protocols/toolbox/task_operations/ResidueProbTaskOperation for more.
///
/// CDR definitions used are North/Dunbrack as the clusters are defined using it.
///
/// @details If Cluster-based profiles cannot be used, will use the fallback strategy.
/// This can happen if the the CDR is of an unknown cluster or there is too little data
/// about the cluster to use profiles.
///
///  See protocols/antibody/design/AntibodyDesignEnum; SeqDesignStrategyEnum
///  for possible fallback strategies. Right now, only conservative or basic/none are implemented.
///
/// FALLBACK STRATEGIES:
///    seq_design_conservative adds a conservative mutation set to the possible residue types (blosum62 default),
///    seq_design_basic will do nothing (as the default for design is to allow all residue positions);
///    seq_design_none will disable design for that CDR (essentially your saying that if it doesn't have profiles, don't design it)
///
/// Due to constness of the apply method, cannot store which CDRs used the fall back strategy.
/// Functions in antibody/database/util can be used to query the database for the number of datapoints
/// Functions in AntibodyInfo can query what the CDR cluster is (NA for unknown clusters)
///
/// This TaskOperation is not currently recommended for H3 as it does not cluster well.
///
/// Optionally sample whole CDR sequences via the primary strategy of:
/// seq_design_profile_sets (use sets instead of profile probability)
///     seq_design_profile_sets_combined (use profile sets and profile probability)
///
class AddCDRProfilesOperation : public core::pack::task::operation::TaskOperation {
public:

	AddCDRProfilesOperation();

	AddCDRProfilesOperation(AntibodyInfoCOP ab_info);

	AddCDRProfilesOperation(AntibodyInfoCOP ab_info, utility::vector1< bool > const & cdrs);

	AddCDRProfilesOperation(AddCDRProfilesOperation const & src);

	virtual ~AddCDRProfilesOperation();

	core::pack::task::operation::TaskOperationOP
	clone() const;

	/// @brief Configure from a RosettaScripts XML tag.
	virtual void
	parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & );

	////////////////////

	/// @brief Pre load the profile data with this function instead of loading it when we apply.
	/// Use this function after all settings are ready to save time with each apply.
	/// This is needed due to const apply for TaskOps
	void
	pre_load_data(core::pose::Pose const & pose);

	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;



	/// @brief Set the CDR-level options as opposed to the defaults.
	void
	set_design_options(design::AntibodyCDRSeqDesignOptions seq_design_options);

	void
	set_cdr_only(CDRNameEnum cdr);

	void
	set_cdrs(utility::vector1< bool > const & cdrs);




	//////// Set stategy options in CDRSeqDesignOptions //////

	/// @brief Set the primary strategy for all CDRs.  See AntibodyDesignEnum.hh for more
	/// Default is Cluster-based Sequence Probabilities - seq_design_profiles
	void
	set_primary_strategy(design::SeqDesignStrategyEnum primary_strategy);

	/// @brief Set the fallback strategy for all CDRs. If the primary strategy could not be
	/// completed due to lack of data, will use this fallback strategy.
	/// See AntibodyDesignEnum.hh for more
	///
	/// Default is Design using the set of Conservative Mutations for each position (Blosum-62) - seq_design_conservative
	/// Set to seq_design_none to disable CDRs as the fallback instead of designing them.
	///
	void
	set_fallback_strategy(design::SeqDesignStrategyEnum fallback_strategy);


	////////

	/// @brief Set the number of times a sequence each chosen.  Increase this number to increase variability of design.
	/// Default 1 round
	void
	set_picking_rounds(core::Size rounds);

	/// @brief Add to the current set of amino acids in the task or replace them?
	/// Default False
	void
	set_add_to_current(bool add_to_current);

	/// @brief Include the poses current residue type in the allowed amino acids.
	/// Default True.
	void
	set_include_native_type(bool use_native);

	/// @brief Use cluster outliers as defined using DihedralDistance and RMSD.
	/// Default false.
	void
	set_use_outliers( bool use_outliers);

	void
	set_defaults();

	/// @brief Force the use of the north paper DB.  Used for benchmarking and Unit Tests.
	void
	set_force_north_paper_db(bool force_north_db);

	/// @brief Set the cutoff.  Will not add the profile set if the total is less than or equal to this number.
	/// Default is 10.
	void
	set_stats_cutoff( core::Size stats_cutoff);

	/// @brief For residue types that have a probability of 0, use this setting
	/// to give a probability to them that is not zero.
	/// Used to increase variability of designs.
	void
	set_sample_zero_probs_at( core::Real zero_prob_sample);

	/// @brief Set the data source for conservative design.  Default is blosum62.  Increased blosum are more stringent, more conservative design.
	///
	/// @details
	///legal = 'chothia_1976', 'BLOSUM30', 'blosum30', 'BLOSUM35', 'blosum35', 'BLOSUM40', 'blosum40',
	///  'BLOSUM45', 'blosum45', 'BLOSUM50', 'blosum50', 'BLOSUM55', 'blosum55', 'BLOSUM60', 'blosum60',
	/// 'BLOSUM62', 'blosum62', 'BLOSUM65', 'blosum65', 'BLOSUM70', 'blosum70', 'BLOSUM75', 'blosum75',
	///'BLOSUM80', 'blosum80', 'BLOSUM85', 'blosum85', 'BLOSUM90', 'blosum90', 'BLOSUM100', 'blosum100'
	void
	set_cons_design_data_source(std::string data_source);


private:

	utility::vector1<bool>
	get_profile_and_design_cdrs() const;

	utility::vector1<bool>
	get_profile_set_and_design_cdrs() const;

	utility::vector1<bool>
	get_design_cdrs() const;

private:

	AntibodyInfoCOP ab_info_;
	design::AntibodyCDRSeqDesignOptions seq_design_options_;
	toolbox::task_operations::ConservativeDesignOperationOP cons_task_;
	AddCDRProfileSetsOperationOP profile_sets_task_;

	//Profile Options
	core::Size picking_rounds_;
	bool keep_task_allowed_aas_;
	bool include_native_restype_;
	bool force_north_paper_db_;
	bool use_outliers_;

	core::Size stats_cutoff_;
	core::Real zero_prob_sample_;
	std::string cons_design_data_source_;

	bool pre_loaded_data_;
	std::map< core::Size, std::map< core::chemical::AA, core::Real > > prob_set_;
	utility::vector1<bool> no_profile_data_cdrs_;
	utility::vector1<bool> no_profile_sets_data_cdrs_;

	///Needed for default and RS constructor.
	AntibodyNumberingSchemeEnum numbering_scheme_;

};

class AddCDRProfilesOperationCreator : public core::pack::task::operation::TaskOperationCreator {
public:
	virtual core::pack::task::operation::TaskOperationOP create_task_operation() const;
	virtual std::string keyname() const { return "AddCDRProfilesOperation"; }
	//core::pack::task::operation::TaskOperationOP clone() const;
};

} //task_operations
} //antibody
} //protocols

#endif //INCLUDED_protocols_antibody_task_operations_AddCDRProfilesOperation_hh



