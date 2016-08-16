// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/resource_manager/planner/JD2ResourceManagerJobInputter.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_jd2_JD2ResourceManagerJobInputter_hh
#define INCLUDED_protocols_jd2_JD2ResourceManagerJobInputter_hh

// Package headers
#include <protocols/jd2/JobInputter.hh>

// Basic headers
#include <basic/options/keys/OptionKeys.hh>
#include <basic/resource_manager/JobOptions.fwd.hh>
#include <protocols/jd2/JobsContainer.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

//C++ headers
#include <istream>
#include <string>

namespace protocols {
namespace jd2 {

/// @brief The %JD2ResourceManagerJobInputter works with the JobDistributor to define
/// the set of jobs that are to be executed.  It reads the list of jobs and the set
/// of resouces they require from one or more XML files.
///
/// @details The %JD2ResourceManagerJobInputter, in its fill_jobs() method, reads
/// a set of XML files that are provided on the command line with the
/// "jd2::resource_definition_files" flag.  Multiple files can be given with this
/// flag, and each will be read in the order that they are given. e.g.
/// -jd2::resource_definition_files file1.xml file2.xml
///
/// The XML file format for defining jobs is structured into blocks
/// of declarations of objects of the same type.  The blocks may be specified
/// in any order and may be repeated as many times as is necessary.
/// For example, it is useful to first declare a DatabaseConnection resource
/// in a Resources block before then declaring a DatabaseResourceLocator object
/// (which requires a database connection) in a ResourceLocators block, before
/// then declaring resources in a second Resources block which may rely on the
/// DatabaseResourceLocator for their creation.  The key is that the data is
/// processed in the order in which it is declared.  If a Resource lists a
/// ResourceLocator that isn't declared until later in the file, the
/// %JD2ResourceManagerJobInputter will exit with an error message; it will not
/// wait until the entire file has been read to try and resolve whether the
/// requested ResourceLocator will eventually be declared.
///
/// The structure of this XML file is given below.
/// \verbatim
/* (<--- please ignore this; it is merely an artifact of writing doxygen)

<JD2ResourceManager>
<ResourceLocators>
... ResourceLocator definitions go in this block
</ResourceLocators>
<ResourceOptions>
... ResourceOption definitions go in this block
</ResourceOptions>
<Resources>
... Resource definitions go in this block
</Resources>
<Jobs>
Job definitions go here.  Currently the structure is:
List a) options and b) resources that are used by all jobs
List c) individual jobs
List options and resources for an individual job
or,
Provide d) a table of jobs and the resources that go with them
from a database.

a) Give a command line option that will apply to all jobs.
<Option optionname=value/>
Each option's name must be given with full namespacing.  (This hopefully will change)
E.g.: if you want to set the -ex1 flag to true, you must pass it as
"packing:ex1=1"
with exactly one colon between the namespace (packing) and the flag (ex1).

NOTE: Not all command line options are currently controlled by
the ResourceManager.  Check first if the option you're interested
in providing is controlled by the ResouceManager by looking...
(APL note to self: figure out where someone should look to tell
if an option is controlled by the ResourceManager)

b) Give a resource that should be used by all jobs
<Data desc=resource_description resource_tag=resource_name/>
* where resource_description is the string that a protocol will use to find the desired resource
e.g. "native", and
* where resource_name is the string that was given as a name to some previously-declared
resource.
or
<Data desc=startstruct pdb=pdbfilename/>
* where pdbfilename is the string that identifies a particular pdb located on the file system
or
<Data desc=startstruct pdb=pdbidentifier locator=locator_tag/>
* where pdbidentifier is a string that identifies a particular pdb, and
* where locator_tag is a string that was given as a name to some previously-declared
ResourceLocator.
E.g., if the given locator_tag names a DatabaseResourceLocator, then the pdbidentifier
will be used to construct an SQL query that rerieves a pdb-file string from a database.

c) List individual jobs
<Job name=jobname nstruct=number>
* where name=jobname is optional and allows you to control the name given to a job
and thus to the name of the output structures the job creates.  If no name is given
then the name will be taken from the "startstruct" data for this job.
And,
* where nstruct=number is optional and allows you to control the number of times this
job is to be performed.  If this is not provided, then the number of iterations is
set to 1.  number should be an integer string, e.g. "nstruct=100"

Note: All jobs must have a "startstruct" specified, either individually for the
given job, or through the Data block above that applies to all jobs.  If the
startstruct is given in the Data block above that applies to all jobs, then each
individual job needs to be given a different jobname.

provide options for this job
<Option optionname=value/>
* where the format for optionname and value are exactly as described above

provide the resources for this job
<Data desc=resource_description resource_tag=resource_name/> or
<Data desc=startstruct pdb=pdfilename/> or
<Data desc=startstruct pdb=pdbidentifier locator=locator_tag/>
* where the formats for these three types of resource statements are
exactly as described above
</Job>

d) Provide a table of jobs from a database
<JobsTable sql_command=command database_resource=dbresource/>
* where command is a string representing an SQL query, and
* where dbresource is the name for a previously-declared database-connection resource
connection resource.
Each row of the jobs table should have one of the following formats
job_name, 'Resource', desc, resource_tag
job_name, 'Option', option_key, option_value
* desc: A job-agnostic description for a resource like 'native' or 'symm_data'
that can be referenced in the protocol.
* resource_tag: The tag of a resource described in the <Resources/> block.
* option_key: An optionally namespaced option key (string) for the options
system like 'in:file:native'
* option_value: A value or list of values that processed into the option system
</Jobs>
</JD2ResourceManager>
\endverbatim
*/
/// See the documentation in the JD2ResourceManager for a description of the file format
/// for the ResourceLocators block (JD2ResourceManager::read_resource_locators_tags),
/// the ResourceOptions block (JD2ResourceManager::read_resource_options_tags) and the
/// Resources block (JD2ResourceManager::read_resources_tags).
class JD2ResourceManagerJobInputter : public JobInputter
{
public:
	JD2ResourceManagerJobInputter();
	virtual ~JD2ResourceManagerJobInputter();

	/// @brief this function is responsible for filling the pose reference with the pose indicated by the job.  The Job object (within its InnerJob) contains a PoseCOP.  This function needs to either fill the pose reference from the InnerJob or, on first demand of a pose from that InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill the reference.
	virtual void pose_from_job( core::pose::Pose & pose, JobOP job );

	/// @brief this function determines what jobs exist.  This function neither knows nor cares what jobs are already complete on disk/memory - it just figures out what ones should exist given the input.  NOTE: your JobInputter should order Job objects in the Jobs vector to have as few "transitions" between inputs as possible (group all Jobs of the same input next to each other).  This improves efficiency of the "FAIL_BAD_INPUT" functionality.  Note I said "should", not "must".
	virtual void fill_jobs( JobsContainer & jobs );

	/// @brief return the type of input source that the JobInputter is currently
	///  using
	virtual JobInputterInputSource::Enum input_source() const;

	void
	fill_jobs_from_stream(
		std::istream & instream,
		Jobs & jobs
	);

private:

	/// @brief Delete the Resources associated with a job
	virtual
	void
	cleanup_after_job_completion(
		std::string const & job_tag);

	void
	parse_jobs_tags(
		utility::tag::TagCOP tag,
		Jobs & jobs
	);

	void
	parse_job_tag(
		utility::tag::TagCOP jobs_tags,
		std::map< std::string, std::string > const & generic_resources_for_job,
		basic::resource_manager::JobOptions const & generic_job_options,
		Jobs & jobs
	);

	void
	parse_jobs_table_tag(
		utility::tag::TagCOP tag,
		std::map< std::string, std::string > const & generic_resources_for_job,
		basic::resource_manager::JobOptions const & generic_job_options,
		Jobs & jobs
	);

	void
	record_job(
		std::string const & job_name,
		std::map< std::string, std::string > const & resources_for_job,
		basic::resource_manager::JobOptionsOP job_options,
		Jobs & jobs
	);

	void
	read_Option_subtag_for_job(
		utility::tag::TagCOP options_tag,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	parse_options_name_and_value(
		std::string const & optname,
		std::string const & value,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_BooleanOption_subtag_for_job(
		basic::options::BooleanOptionKey const & boolopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_FileOption_subtag_for_job(
		basic::options::FileOptionKey const & fileopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_IntegerOption_subtag_for_job(
		basic::options::IntegerOptionKey const & intopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_PathOption_subtag_for_job(
		basic::options::PathOptionKey const & pathopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_RealOption_subtag_for_job(
		basic::options::RealOptionKey const & realopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_StringOption_subtag_for_job(
		basic::options::StringOptionKey const & stringopt,
		std::string const & optname,
		std::string const & val,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_BooleanVectorOption_subtag_for_job(
		basic::options::BooleanVectorOptionKey const & boolvectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_FileVectorOption_subtag_for_job(
		basic::options::FileVectorOptionKey const & filevectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_IntegerVectorOption_subtag_for_job(
		basic::options::IntegerVectorOptionKey const & intvectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_PathVectorOption_subtag_for_job(
		basic::options::PathVectorOptionKey const & pathvectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_RealVectorOption_subtag_for_job(
		basic::options::RealVectorOptionKey const & realvectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);

	void
	read_StringVectorOption_subtag_for_job(
		basic::options::StringVectorOptionKey const & strinvectopt,
		std::string const & optname,
		std::string const & val,
		utility::vector1< std::string > const & vals,
		basic::resource_manager::JobOptionsOP job_options
	);


	void
	read_Data_for_subtag(
		utility::tag::TagCOP options_tag,
		std::string const & jobname,
		std::string & input_tag,
		std::map< std::string, std::string > & resources_for_job
	);

	void
	read_ResidueType_for_subtag(
		utility::tag::TagCOP options_tag,
		std::map< std::string, std::string > & resources_for_job
	);

	void
	check_each_job_has_startstruct(
		Jobs const & jobs
	) const;

private:
	/// @brief save the last input tag so the resources loaded for it can
	///be unloaded when we see a new input tag
	std::string last_input_tag_;

};

} // namespace jd2
} // namespace protocols

#endif
