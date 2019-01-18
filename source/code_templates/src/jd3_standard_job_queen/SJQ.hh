// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.hh
/// @brief --brief--
/// @author --name-- (--email--)

#ifndef INCLUDED_--path_underscore--_--class--_HH
#define INCLUDED_--path_underscore--_--class--_HH

// Unit headers
#include <--path--/--class--.fwd.hh>
#include <protocols/jd3/standard/StandardJobQueen.hh>

// protocol headers
#include <protocols/jd3/deallocation/DeallocationMessage.fwd.hh>
#include <protocols/jd3/output/OutputSpecification.fwd.hh>
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <protocols/jd3/dag_node_managers/NodeManager.fwd.hh>
#include <protocols/jd3/dag_node_managers/SimpleNodeManager.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>
#include <protocols/jd3/JobGenealogist.fwd.hh>

// utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/resource_manager/ResourceManager.fwd.hh>

// utility headers
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>

// C++ headers
#include <list>
#include <map>
#include <set>

	///////////////////// How to Write an SJQ //////////////////////////
	//
	// https://www.rosettacommons.org/docs/wiki/development_documentation/tutorials/jd3_derived_jq/jd3_derived_jq_home
	// https://wiki.rosettacommons.org/index.php/JD3FAQ
	//
	/////////////////////////////////////////////////////////////////////





--namespace--

///@brief --brief--
class --class-- : public protocols::jd3::standard::StandardJobQueen {


public:
	///@brief --brief--
	--class--();

	~--class--();

	/// These are common functions to override.  This is not ALL funcitons that can be overriden!
	// https://www.rosettacommons.org/docs/latest/development_documentation/tutorials/jd3_derived_jq/checklist
	//

public:

	////////////////////////
	/// Starting the Job ///
	////////////////////////

	/// This is the Job graph.  Here is where you can create what your initial job(s) will look like.
	/// This is stored in the JQ by the JD.
	///
	protocols::jd3::JobDigraphOP
	initial_job_dag() override;

	///@brief Get the Jobs ready to go.
	///
	/// It is recommended to override next_batch_of_larval_jobs_for_job_node instead of this method
	///  However, if you do override this - especially for more complicated protocols - be sure to call similar functions here
	///
	std::list< protocols::jd3::LarvalJobOP >
	determine_job_list (
        	core::Size job_dag_node_index,
        	core::Size max_njobs
	) override;

	///@brief Called during the SJQs determine_job_list function.
	/// Override this function if you are calling the SJQ determine job list.
	/// This is used to to create LarvalJobs for nodes other than preliminary job nodes.
	protocols::jd3::LarvalJobs
	next_batch_of_larval_jobs_for_job_node( Size job_dag_node_index, Size max_njobs ) override;

protected:

	///@brief /// This is the core of your code.
	protocols::jd3::JobOP
	complete_larval_job_maturation(
	        protocols::jd3::LarvalJobCOP larval_job,
	        utility::options::OptionCollectionCOP job_options,
	        utility::vector1< protocols::jd3::JobResultCOP > const & result
	) override;

public:

	///////////////////////////
	//// Completing the Job ///
	///////////////////////////

	///@brief The JD calls this function on completion of a LarvalJob, after updating the JobTracker
	///
	///@details
	///  If you override this method, call the SJQs version first for PJN tracking and output.
	void
	note_job_completed( protocols::jd3::LarvalJobCOP job, protocols::jd3::JobStatus status, core::Size nresults ) override;

	///@brief As each job completes, this function is called for each result coming from a LarvalJob.
	///
	///@details
    /// Override this method to store any results you need.
	void completed_job_summary(
		    jd3::LarvalJobCOP job,
 		    core::Size result_index,
 		    jd3::JobSummaryOP summary
	) override;

	///@brief Note which jobs we don't need for the next set of Job Nodes.
	std::list< jd3::JobResultID >
	job_results_that_should_be_discarded() override;

	///@brief Figure out output.  Not entirely needed.
	std::list< jd3::output::OutputSpecificationOP >
	jobs_that_should_be_output() override;

protected:

	//////////////////
	/// XML Schema ///
	//////////////////

	///@brief Here, you specify which tags that are defined in the <common> element.
	void
	append_common_tag_subelements(
	        utility::tag::XMLSchemaDefinition & job_definition_xsd,
	        utility::tag::XMLSchemaComplexTypeGenerator & ct_gen
	) const override;

	///@brief Here you can specify Schema's to use when creating your XSD for a specific JOB.
	/// AKA - you can add elements such as ResidueSelectors, task operations, etc. to parse.
	void
	append_job_tag_subelements(
	        utility::tag::XMLSchemaDefinition & job_definition_xsd,
	        utility::tag::XMLSchemaComplexTypeGenerator & job_ct_gen
	) const override;








private:


};



--end_namespace--

#endif //--path_underscore--_--class--_HH


