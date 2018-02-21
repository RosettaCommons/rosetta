// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/multistage_rosetta_scripts/MRSJobQueen.hh
/// @brief Derived Job Queen for Multistage Rosetta Scripts
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_protocols_multistage_rosetta_scripts_MRSJobQueen_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_MRSJobQueen_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/multistage_rosetta_scripts/MRSJobQueen.fwd.hh>
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.fwd.hh>
#include <protocols/multistage_rosetta_scripts/TagManager.hh>

//JD3
#include <protocols/jd3/JobQueen.hh>
#include <protocols/jd3/standard/StandardJobQueen.hh>
#include <protocols/jd3/JobGenealogist.hh>
#include <protocols/jd3/Job.fwd.hh>
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <protocols/jd3/JobOutputIndex.fwd.hh>
#include <protocols/jd3/JobResult.fwd.hh>
#include <protocols/jd3/JobSummary.fwd.hh>
#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/dag_node_managers/NodeManager.hh>
#include <protocols/jd3/standard/MoverAndPoseJob.hh>
#include <protocols/jd3/output/OutputSpecification.fwd.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <basic/datacache/DataMap.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/XMLSchemaValidation.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/options/OptionCollection.fwd.hh>

#include <unordered_map>
#include <boost/functional/hash/hash.hpp>

namespace protocols {
namespace multistage_rosetta_scripts {

struct JobResultID_hash {
	std::size_t operator () ( jd3::JobResultID const & id ) const {
		std::size_t seed = 0;
		boost::hash_combine( seed, id.first );
		boost::hash_combine( seed, id.second );
		return seed;
	}
};

///@brief We do not want to load all of the input poses into memory at once. Instead we use this struct to keep track of the most recent pose loaded (which is assumed to be the most likely one we are going to need next)
struct PoseForPoseID
{
	PoseForPoseID( core::Size pose_id_arg, core::pose::PoseOP const & pose_arg ) :
		pose_id( pose_id_arg ),
		pose( pose_arg )
	{}
	core::Size pose_id;
	core::pose::PoseOP pose;
};

struct SortByLowEnergy
{
	bool operator() (
		const std::pair< core::Size, jd3::standard::EnergyJobSummaryOP > & left,
		const std::pair< core::Size, jd3::standard::EnergyJobSummaryOP > & right
	) const {
		return left.second->energy() < right.second->energy();
	}
};

class MRSJobQueen: public jd3::standard::StandardJobQueen {

public:

	//constructor
	MRSJobQueen();

	//destructor
	~MRSJobQueen() override;

public:

	jd3::JobDigraphOP
	initial_job_dag()
	override;

	std::list< jd3::LarvalJobOP > determine_job_list(
		Size job_dag_node_index,
		Size max_njobs
	) override;

	jd3::JobOP
	complete_larval_job_maturation(
		jd3::LarvalJobCOP larval_job,
		utility::options::OptionCollectionCOP job_options,
		utility::vector1< jd3::JobResultCOP > const & input_job_results
	) override ;

	///@brief this was only created for the unit test. Please do not call this.
	void note_job_completed(
		core::Size job_id,
		jd3::JobStatus status,
		core::Size nresults,
		bool are_you_a_unit_test
	);

	void note_job_completed(
		core::Size,
		jd3::JobStatus,
		core::Size
	) override {
		runtime_assert( false );
	}

	void note_job_completed(
		jd3::LarvalJobCOP job,
		jd3::JobStatus status,
		core::Size nresults
	) override;

	void completed_job_summary(
		core::Size job_id,
		core::Size result_index,
		jd3::JobSummaryOP summary
	) override;

	void completed_job_summary(
		jd3::LarvalJobCOP job,
		core::Size result_index,
		jd3::JobSummaryOP summary
	) override {
		completed_job_summary( job->job_index(), result_index, summary );
	}

	core::Size stage_for_global_job_id( core::Size global_job_id ) const;

	std::list< jd3::JobResultID > job_results_that_should_be_discarded() override;

	std::list< jd3::output::OutputSpecificationOP > jobs_that_should_be_output() override;

	TagManager const & tag_manager() const {
		return tag_manager_;
	}

	void determine_validity_of_stage_tags();

	core::Size num_input_structs() const { return num_input_structs_; }
	core::Size num_stages() const { return num_stages_; }
	core::Size num_results_to_keep_for_stage( core::Size stage ) const {
		return num_results_to_keep_for_stage_[ stage ];
	}

protected:

	core::pose::PoseOP pose_for_job_derived(
		jd3::LarvalJobCOP job,
		utility::options::OptionCollection const & options
	);


	core::pose::PoseOP pose_for_inner_job_derived(
		jd3::standard::StandardInnerLarvalJobCOP inner_job,
		utility::options::OptionCollection const & options
	);

	void all_settings_to_defualt();

	jd3::LarvalJobOP get_nth_job_for_initial_stage( core::Size local_job_id );
	jd3::LarvalJobOP get_nth_job_for_noninitial_stage( core::Size stage, core::Size local_job_id );

	void
	append_common_tag_subelements(
		utility::tag::XMLSchemaDefinition & xsd,
		utility::tag::XMLSchemaComplexTypeGenerator & ct_gen
	) const override ;

	void
	parse_job_definition_tags(
		utility::tag::TagCOP common_block_tags,
		utility::vector1< jd3::standard::PreliminaryLarvalJob > const &
	) override;

	void parse_common_tag( utility::tag::TagCOP common_tag );
	void parse_single_stage_tag( utility::tag::TagCOP subprotocol_tag );

	void
	append_job_tag_subelements(
		utility::tag::XMLSchemaDefinition & job_definition_xsd,
		utility::tag::XMLSchemaComplexTypeGenerator & job_ct_gen
	) const override ;

	void parse_single_job_tag( jd3::standard::PreliminaryLarvalJob const & prelim_larval_job, core::Size input_pose_id );

	static std::string protocols_subelement_mangler( std::string const & name ){
		return "protocols_multistage_rosetta_scripts_" + name + "_complex_type";
	}

	static std::string job_subelement_mangler( std::string const & name ){
		return "protocols_multistage_rosetta_scripts_job_" + name + "_complex_type";
	}

	core::Size input_pose_id_for_jobid( core::Size global_jobid ) const;

	void print_job_lineage() const;

	void assign_output_index(
		jd3::LarvalJobCOP larval_job,
		Size result_index_for_job,
		Size n_results_for_job,
		jd3::JobOutputIndex & output_index
	) override;

	void assign_output_index(
		core::Size global_job_id,
		core::Size local_result_id,
		jd3::JobOutputIndex & output_index
	);

	void cluster( core::Size stage_about_to_start );

private:
	bool has_been_initialized_;

	SortByLowEnergy sorter_;

	core::Size num_input_structs_;
	core::Size num_stages_;

	utility::vector1< core::Size > num_total_jobs_for_stage_;
	utility::vector1< core::Size > num_results_to_keep_for_stage_;
	utility::vector1< core::Size > max_num_results_to_keep_per_instance_for_stage_;
	utility::vector1< core::Size > max_num_results_to_keep_per_input_struct_for_stage_;
	utility::vector1< core::Size > num_jobs_per_input_for_stage_;
	utility::vector1< core::Size > result_cutoffs_for_stage_;
	utility::vector1< core::Size > job_results_have_been_discarded_for_stage_;

	utility::vector1< jd3::dag_node_managers::NodeManagerOP > node_managers_;

	//XML data
	utility::vector1< utility::tag::TagCOP > tag_for_stage_;

	TagManager tag_manager_;

	utility::tag::XMLSchemaValidatorOP validator_;

	utility::vector1< std::string > outputters_;
	utility::vector1< std::string > input_job_tags_;
	utility::vector1< core::Size > num_structs_output_for_input_job_tag_;

	utility::vector1< std::pair< jd3::InnerLarvalJobOP, core::Size > > current_inner_larval_job_for_stage_;//second element is counter
	jd3::JobGenealogistOP job_genealogist_;

	std::map< jd3::JobResultID, jd3::output::OutputSpecificationOP > pose_output_specification_for_job_result_id_;

	PoseForPoseID most_recent_pose_id_;

	utility::vector1< std::pair< core::Size, utility::tag::TagCOP > > checkpoints_;


	//Data for clustering:
	utility::vector1< utility::tag::TagCOP > cluster_metric_tag_for_stage_;
	utility::vector1< core::Size > num_results_to_keep_after_clustering_for_stage_;
	utility::vector1< std::unordered_map< jd3::JobResultID, cluster::ClusterMetricCOP, JobResultID_hash > > cluster_data_for_results_of_stage_;

	///@brief No need to store this data for more than one stage at a time
	utility::vector1< jd3::JobResultID > most_recent_cluster_results_;
};

inline core::Size MRSJobQueen::stage_for_global_job_id( core::Size global_job_id ) const {
	for ( core::Size stage_ii = 2; stage_ii <= num_stages_; ++stage_ii ) {
		if ( node_managers_[ stage_ii ]->job_offset() >= global_job_id ) {
			return stage_ii - 1;
		}
	}
	return num_stages_;
}

inline core::Size MRSJobQueen::input_pose_id_for_jobid( core::Size global_jobid ) const{
	core::Size const stage = stage_for_global_job_id( global_jobid );
	return job_genealogist_->input_source_for_job( stage, global_jobid );
}


} //multistage_rosetta_scripts
} //protocols

#endif
