// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/resource_manager/planner/JD2ResourceManagerJobInputter.cc
/// @brief
/// @author

// Unit headers
#include <protocols/jd2/JD2ResourceManagerJobInputter.hh>
#include <protocols/jd2/JD2ResourceManagerJobInputterCreator.hh>

// Package headers
#include <protocols/jd2/JD2ResourceManager.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobsContainer.hh>
#include <basic/resource_manager/JobOptions.hh>
#include <protocols/jd2/InnerJob.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/PoseResidueTypeSet.hh>

// Utility Headers
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.hh>
#include <utility/file/FileName.hh>
#include <utility/options/keys/BooleanOptionKey.hh>
#include <utility/options/keys/IntegerOptionKey.hh>
#include <utility/options/keys/FileOptionKey.hh>
#include <utility/options/keys/PathOptionKey.hh>
#include <utility/options/keys/RealOptionKey.hh>
#include <utility/options/keys/StringOptionKey.hh>
#include <utility/string_util.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>

// Boost headers
#include <boost/lexical_cast.hpp>

// External headers
#include <cppdb/frontend.h>

//C++ headers
#include <istream>
#include <string>
#include <sstream>

namespace protocols {
namespace jd2 {

using utility::excn::EXCN_Msg_Exception;
using basic::resource_manager::ResourceOP;
using basic::resource_manager::ResourceManager;
using basic::resource_manager::ResourceTag;
using basic::resource_manager::JobOptions;
using basic::resource_manager::JobOptionsOP;
using core::pose::PoseOP;
using core::pose::Pose;

using std::endl;
using std::string;
using std::stringstream;

static THREAD_LOCAL basic::Tracer tr( "protocols.jd2.JD2ResourceManagerJobInputter" );


//CREATOR SECTION
std::string
JD2ResourceManagerJobInputterCreator::keyname() const
{
	return "JD2ResourceManagerJobInputter";
}

protocols::jd2::JobInputterOP
JD2ResourceManagerJobInputterCreator::create_JobInputter() const {
	return protocols::jd2::JobInputterOP( new JD2ResourceManagerJobInputter );
}


JD2ResourceManagerJobInputter::JD2ResourceManagerJobInputter() :
	last_input_tag_("")
{}

JD2ResourceManagerJobInputter::~JD2ResourceManagerJobInputter() = default;

void
JD2ResourceManagerJobInputter::pose_from_job(
	Pose & pose,
	JobOP job )
{
	tr.Debug << "JD2ResourceManagerJobInputter::pose_from_job" << endl;

	std::string const & input_tag(job->inner_job()->input_tag());

	if ( last_input_tag_ == "" ) {
		last_input_tag_ = input_tag;
	} else {
		if ( last_input_tag_ != input_tag ) {
			cleanup_after_job_completion(last_input_tag_);
			last_input_tag_ = input_tag;
		}
	}

	if ( !job->inner_job()->get_pose() ) {
		tr.Debug
			<< "Retrieving pose from ResourceManager (tag = " << input_tag << ")" << endl;
		pose.clear();
		JD2ResourceManager * jd2_resource_manager(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		debug_assert( jd2_resource_manager != nullptr );

		ResourceOP resource;

		// Check to see if we have a Residue resource, if so load it into the chemical manager if it hasn't already been loaded
		// APL: Sam, is there anywhere else this code could live?
		// RM: It's been moved to PoseFromPDBLoader::create_resource() - not an ideal place, but one that works.

		try {
			tr << "Loading startstruct " << jd2_resource_manager->find_resource_tag_by_job_tag( "startstruct", input_tag ) << " for job " <<
				input_tag <<std::endl;
			resource = jd2_resource_manager->get_resource_by_job_tag("startstruct", input_tag);
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::ostringstream err;
			err << e.msg() << std::endl;
			err << "Failed to access 'startstruct' resource from the JD2ResourceManager for job '";
			err << input_tag << "' with nstruct index " << job->nstruct_index();
			err << "\n" << "Exception caught and re-thrown from JD2ResourceManagerJobInputter::pose_from_job\n";
			throw utility::excn::EXCN_Msg_Exception( err.str() );
		}

		PoseOP resource_pose( utility::pointer::dynamic_pointer_cast< core::pose::Pose > ( resource ) );

		/// make sure the resource that we requested from the resource manager is in fact a pose.
		if ( ! resource_pose ) {
			std::ostringstream err;
			err
				<< "Error trying to access starting structure (startstruct) "
				<< "for job with tag '";
			err << input_tag << "'.";
			if ( jd2_resource_manager->has_resource_tag_by_job_tag("startstruct", input_tag) ) {
				ResourceTag rt(jd2_resource_manager->find_resource_tag_by_job_tag("startstruct", input_tag));
				err
					<< " The resource tag '" << rt << "' "
					<< "for this job is not for a Pose object" << endl;
			} else {
				err << " This starting structure for this job is critically absent\n";
			}
			throw EXCN_Msg_Exception( err.str() );
		}

		/// Save a copy of the resource_pose into the inner job.
		load_pose_into_job(resource_pose, job);
		pose = *resource_pose;
	} else {
		pose.clear();
		pose = *(job->inner_job()->get_pose());
		tr.Debug
			<< "filling pose from saved copy (tag = "
			<< job->input_tag() << ")" << endl;
	}

}

void
JD2ResourceManagerJobInputter::fill_jobs( JobsContainer & jobs )
{
	tr << "Initializing jobs from resource definition files" << std::endl;
	utility::vector1< utility::file::FileName > resource_definition_files
		= basic::options::option[ basic::options::OptionKeys::jd2::resource_definition_files ];
	for ( core::Size ii = 1; ii <= resource_definition_files.size(); ++ii ) {

		tr.Debug << "Reading from resource definition file " << resource_definition_files[ii] << std::endl;
		utility::io::izstream instream( resource_definition_files[ii]().c_str() );
		utility::vector1 < JobOP > jobs_list;
		fill_jobs_from_stream( instream, jobs_list );
		for ( core::Size jj=1; jj<=jobs_list.size(); ++jj ) {
			jobs.push_back(jobs_list[jj]);
		}
	}
}

void
JD2ResourceManagerJobInputter::cleanup_after_job_completion(
	std::string const & job_tag)
{
	JD2ResourceManager * jd2_resource_manager(
		JD2ResourceManager::get_jd2_resource_manager_instance());

	jd2_resource_manager->mark_job_tag_as_complete(job_tag);
	std::list<ResourceTag> resources_for_job(jd2_resource_manager->get_resource_tags_for_job_tag(job_tag));

	auto resource_list_it = resources_for_job.begin();
	for ( ; resource_list_it != resources_for_job.end(); ++resource_list_it ) {
		core::Size job_count(jd2_resource_manager->get_count_of_jobs_associated_with_resource_tag(*resource_list_it));
		if ( job_count == 0 ) {
			if ( !jd2_resource_manager->ResourceManager::has_resource(*resource_list_it) ) {
				tr << "Skipping deleting '" << *resource_list_it << "' because it doesn't exist in the ResourceManager." << std::endl;
				continue;
			}

			// We no longer need to explicitly delete ResidueTypes - Since we load them into the pose, they should be freed when the pose is freed.

			tr << "Deleting resource " << *resource_list_it <<std::endl;

			jd2_resource_manager->free_resource_by_tag(*resource_list_it);
		}
	}
}

void
JD2ResourceManagerJobInputter::fill_jobs_from_stream( std::istream & instream, Jobs & jobs )
{
	JD2ResourceManager * jd2_resource_manager(
		JD2ResourceManager::get_jd2_resource_manager_instance());

	utility::tag::TagCOP resource_tags = utility::tag::Tag::create( instream );

	for ( auto
			tag_iter = resource_tags->getTags().begin(),
			tag_iter_end = resource_tags->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter ) {
		std::string const & tagname = (*tag_iter)->getName();
		if ( tagname == "ResourceLocators" ) {
			jd2_resource_manager->read_resource_locators_tags( *tag_iter );
		} else if ( tagname == "ResourceOptions" ) {
			jd2_resource_manager->read_resource_options_tags( *tag_iter );
		} else if ( tagname == "Resources" ) {
			jd2_resource_manager->read_resources_tags( *tag_iter );
		} else if ( tagname == "Jobs" ) {
			parse_jobs_tags( *tag_iter, jobs );
		}
	}
}

void
JD2ResourceManagerJobInputter::parse_jobs_tags(
	utility::tag::TagCOP jobs_tags,
	Jobs & jobs
)
{
	//identify options and resources that apply to each job
	std::map< std::string, std::string > generic_resources_for_job;
	JobOptionsOP generic_job_options( new JobOptions() );

	for ( auto
			tag_iter = jobs_tags->getTags().begin(),
			tag_iter_end = jobs_tags->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter ) {
		std::string const & tagname = (*tag_iter)->getName();
		if ( tagname == "Option" ) {
			read_Option_subtag_for_job( *tag_iter, generic_job_options );
		} else if ( tagname == "Data" ) {
			std::string dummy_input_tag;
			read_Data_for_subtag(
				*tag_iter, "generic", dummy_input_tag,
				generic_resources_for_job );
		}
	}


	for ( auto
			tag_iter = jobs_tags->getTags().begin(),
			tag_iter_end = jobs_tags->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter ) {
		std::string const & tagname = (*tag_iter)->getName();
		if ( tagname == "Job" ) {
			parse_job_tag( *tag_iter, generic_resources_for_job, *generic_job_options, jobs );
		} else if ( tagname == "JobsTable" ) {
			parse_jobs_table_tag( *tag_iter, generic_resources_for_job, *generic_job_options, jobs );
		} else {
			std::ostringstream err;
			err << "Error parsing jobs tags in JD2ResourceManagerJobInputter: unrecognized tag '" << tagname << "'";
			EXCN_Msg_Exception( err.str() );
		}
	}

	check_each_job_has_startstruct(jobs);
}

void
JD2ResourceManagerJobInputter::parse_job_tag(
	utility::tag::TagCOP jobs_tag,
	std::map< std::string, std::string > const & generic_resources_for_job,
	JobOptions const & generic_job_options,
	Jobs & jobs
)
{

	//InnerJobOP inner_job = new InnerJob;
	std::string jobname;
	std::string input_tag; // <--- pdb name
	std::map< std::string, std::string > resources_for_job(generic_resources_for_job);


	/// if the name is not given
	if ( jobs_tag->hasOption( "name" ) ) {
		jobname = jobs_tag->getOption< std::string >( "name" );
	}

	JobOptionsOP job_options( new JobOptions(generic_job_options) );
	for ( auto
			tag_iter = jobs_tag->getTags().begin(),
			tag_iter_end = jobs_tag->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter ) {
		std::string const & tagname = (*tag_iter)->getName();
		if ( tagname == "Option" ) {
			read_Option_subtag_for_job( *tag_iter, job_options );
		} else if ( tagname == "Data" ) {
			read_Data_for_subtag( *tag_iter, jobname, input_tag, resources_for_job );
		} else if ( tagname == "ResidueType" ) {
			read_ResidueType_for_subtag(*tag_iter, resources_for_job);
		}// are there other kinds of tags?
	}

	if ( jobs_tag->hasOption("nstruct") ) {
		// if it's not specified here, it can be specified in an <Option/>
		// tag or default to 1.
		parse_options_name_and_value(
			"nstruct", jobs_tag->getOption<string>("nstruct"), job_options);
	}

	if ( jobname.size() == 0 ) {
		JD2ResourceManager * jd2rm( JD2ResourceManager::get_jd2_resource_manager_instance());
		if ( jd2rm->has_resource_configuration( input_tag ) ) {
			std::ostringstream err;
			err << "Repeat input startstruct for a job which has not been given a name. Ordinarily, if a job is\n";
			err << "given without a name, it is assigned a name based on its startstruct.  If, however, the same\n";
			err << "startstruct is given multiple times, then we cannot ensure that job names are unique.  Please\n";
			err << "give each job with a shared startstrct a different name.\n";
			err << "Offeding startstruct: '" << input_tag << "\n";
			throw EXCN_Msg_Exception( err.str() );
		}

		jobname = input_tag;
	}

	record_job(jobname, resources_for_job, job_options, jobs);
}

void
JD2ResourceManagerJobInputter::parse_jobs_table_tag(
	utility::tag::TagCOP tag,
	std::map< std::string, std::string > const & generic_resources_for_job,
	JobOptions const & generic_job_options,
	Jobs & jobs
) {
	using namespace basic::database;

	utility::sql_database::sessionOP db_session = parse_database_connection(tag);

	string sql_command;
	if ( tag->hasOption("sql_command") ) {
		sql_command = tag->getOption<string>("sql_command");
		check_statement_sanity(sql_command);
	} else {
		stringstream err_msg;
		err_msg
			<< "The JobsTable tag requires a 'sql_command' tag that "
			<< "is an SQL SELECT statement that returns the following column formats" << endl
			<< "ordered by <job_name>:" << endl
			<< "\t(<job_name>, 'Resource', <desc>, <resource_tag>)" << endl
			<< "\t(<job_name>, 'Option', <option_key>, <option_value>)" << endl;
		throw utility::excn::EXCN_Msg_Exception(err_msg.str());
	}

	cppdb::statement select_stmt(safely_prepare_statement(sql_command, db_session));
	cppdb::result res(safely_read_from_database(select_stmt));

	std::string job_table_schema =
		"Each row of the jobs table should have one of the following formats\n"
		"\n"
		"    job_name, 'Resource', desc, resource_tag\n"
		"    job_name, 'Option', option_key, option_value\n"
		"\n"
		" * desc: A job-agnostic description for a resource like 'native' or 'symm_data'\n"
		"   that can be referenced in the protocol.\n"
		"\n"
		" * resource_tag: The tag of a resource described in the <Resources/> block.\n"
		"\n"
		" * option_key: An optionally namespaced option key (string) for the options\n"
		"   system like 'in:file:native'\n"
		"\n"
		" * option_value: A value or list of values that processed into the option system\n";

	if ( res.cols() != 4 ) {
		stringstream err_msg;
		err_msg
			<< "The JobsTable tag requires a 'sql_command' tag" << endl
			<< job_table_schema << endl
			<< "Instead, the query returned " << res.cols() << ":" << endl
			<< "SQL query:" << endl
			<< sql_command << endl;
		throw utility::excn::EXCN_Msg_Exception(err_msg.str());
	}

	Size row_number(0);
	string previous_job_name("");
	string input_tag;
	std::map< string, string > resources_for_job(generic_resources_for_job);
	JobOptionsOP job_options( new JobOptions(generic_job_options) );
	while ( res.next() ) {
		row_number++;

		string job_name, resource_type, key, value;
		res >> job_name >> resource_type >> key >> value;

		if ( row_number != 1 && previous_job_name != job_name ) {
			// we've just finished collecting information for the job,
			// record job and reset the invariants
			record_job(previous_job_name, resources_for_job, job_options, jobs);
			resources_for_job = generic_resources_for_job;
			job_options = JobOptionsOP( new JobOptions(generic_job_options) );
			previous_job_name = job_name;

		}


		if ( resource_type == "Resource" ) {
			resources_for_job[ key ] = value;
		} else if ( resource_type == "Option" ) {
			parse_options_name_and_value(key, value, job_options);
		} else {
			stringstream err_msg;
			err_msg
				<< "Unrecognized data type '" << resource_type << "' for job '" << job_name << "'" << endl
				<< "The JobsTable tag requires a 'sql_command' tag" << endl
				<< job_table_schema << endl
				<< endl
				<< "SQL query:" << endl
				<< sql_command << endl
				<< endl
				<< "Resources:" << endl;
			for ( std::map<string,string>::const_iterator i=resources_for_job.begin(), ie=resources_for_job.end(); i != ie; ++i ) {
				err_msg << "\t '" << i->first << "' <- '" << i->second << "'" << endl;
			}
			err_msg
				<< "Options:" << endl
				<< job_options << endl;

			throw EXCN_Msg_Exception(err_msg.str());
		}

		if ( previous_job_name.empty() ) {
			previous_job_name = job_name;
		}
	}

	if ( row_number == 0 ) {
		tr.Warning << "JobsTable returned no rows." << endl;
	}

	record_job(previous_job_name, resources_for_job, job_options, jobs);
}

void
JD2ResourceManagerJobInputter::record_job(
	string const & job_name,
	std::map< string, string > const & resources_for_job,
	JobOptionsOP job_options,
	Jobs & jobs
) {

	JD2ResourceManager * jd2rm(
		JD2ResourceManager::get_jd2_resource_manager_instance());

	using namespace basic::options;

	for ( auto const & iter : resources_for_job ) {
		jd2rm->add_resource_tag_by_job_tag( iter.first, job_name, iter.second );
	}

	JD2ResourceManager::get_jd2_resource_manager_instance()->add_job_options(
		job_name, job_options);

	//Have the jobs with this job_name already been created?
	core::Size n_jobs(0);
	for ( Jobs::const_iterator ii = jobs.begin(), iie = jobs.end(); ii != iie; ++ii ) {
		if ( (*ii)->input_tag() == job_name ) {
			n_jobs++;
		}
	}

	//If not, then add them
	if ( n_jobs==0 ) {
		Size nstruct(1);
		if ( job_options->has_option(OptionKeys::out::nstruct) ) {
			nstruct = job_options->get_option(OptionKeys::out::nstruct);
		}

		InnerJobOP inner_job( new InnerJob( job_name, nstruct ) );
		for ( Size ii = 1; ii <= nstruct; ++ii ) {
			jobs.push_back( protocols::jd2::JobOP( new Job( inner_job, ii ) ));
		}
	} else {
		if ( job_options->has_option(OptionKeys::out::nstruct) ) {
			Size requested_nstruct = job_options->get_option(OptionKeys::out::nstruct);
			if ( requested_nstruct != n_jobs ) {
				std::stringstream err_msg;
				err_msg
					<< "Conflicting specification for nstruct for job with name '" << job_name << "'" << std::endl
					<< "Previous nstruct=" << n_jobs << ", new nstruct=" << requested_nstruct << std::endl;
				throw EXCN_Msg_Exception( err_msg.str() );
			}
		}
	}
}

void
JD2ResourceManagerJobInputter::read_Option_subtag_for_job(
	utility::tag::TagCOP options_tag,
	JobOptionsOP job_options
)
{
	for ( auto const & opt_iter : options_tag->getOptions() ) {
		std::string const & optname = opt_iter.first;
		std::string const & val = opt_iter.second;
		parse_options_name_and_value(optname, val, job_options);
	}
}


void
JD2ResourceManagerJobInputter::parse_options_name_and_value(
	std::string const & optname,
	std::string const & val,
	JobOptionsOP job_options
)
{
	using namespace basic::options;

	std::string full_key;
	try{
		full_key = option.find_key_cl(optname, "", true);
	} catch (...) {
		std::stringstream err_msg;
		err_msg
			<< "Error: Option key '" << optname << "' not found. Please remember to use only one colon when giving options." << endl;
		throw EXCN_Msg_Exception( err_msg.str() );
	}

	if ( basic::options::OptionKeys::has( full_key ) ) {
		OptionKey const & opt( basic::options::OptionKeys::key( full_key ));
		if ( opt.scalar() ) {
			// scalar options
			if ( dynamic_cast< BooleanOptionKey const * > (&opt) ) {
				BooleanOptionKey const & boolopt( static_cast< BooleanOptionKey const & > (opt) );
				read_BooleanOption_subtag_for_job( boolopt, full_key, val, job_options );
			} else if ( dynamic_cast< FileOptionKey const * > (&opt) ) {
				FileOptionKey const & fileopt( static_cast< FileOptionKey const & > (opt) );
				read_FileOption_subtag_for_job( fileopt, full_key, val, job_options );
			} else if ( dynamic_cast< IntegerOptionKey const * > (&opt) ) {
				IntegerOptionKey const & iopt( static_cast< IntegerOptionKey const & > (opt) );
				read_IntegerOption_subtag_for_job( iopt, full_key, val, job_options );
			} else if ( dynamic_cast< PathOptionKey const * > (&opt) ) {
				PathOptionKey const & pathopt( static_cast< PathOptionKey const & > (opt) );
				read_PathOption_subtag_for_job( pathopt, full_key, val, job_options );
			} else if ( dynamic_cast< RealOptionKey const * > (&opt) ) {
				RealOptionKey const & ropt( static_cast< RealOptionKey const & > (opt) );
				read_RealOption_subtag_for_job( ropt, full_key, val, job_options );
			} else if ( dynamic_cast< StringOptionKey const * > (&opt) ) {
				StringOptionKey const & stopt( static_cast< StringOptionKey const & > (opt) );
				read_StringOption_subtag_for_job( stopt, full_key, val, job_options );
			}

		} else {
			/// vector option
			utility::vector1< std::string > vals = utility::string_split( val, ',' );
			if ( dynamic_cast< BooleanVectorOptionKey const * > (&opt) ) {
				BooleanVectorOptionKey const & boolvectopt( static_cast< BooleanVectorOptionKey const & > (opt) );
				read_BooleanVectorOption_subtag_for_job( boolvectopt, full_key, val, vals, job_options );
			} else if ( dynamic_cast< FileVectorOptionKey const * > (&opt) ) {
				FileVectorOptionKey const & filevectopt( static_cast< FileVectorOptionKey const & > (opt) );
				read_FileVectorOption_subtag_for_job( filevectopt, full_key, val, vals, job_options );
			} else if ( dynamic_cast< IntegerVectorOptionKey const * > (&opt) ) {
				IntegerVectorOptionKey const & ivectopt( static_cast< IntegerVectorOptionKey const & > (opt) );
				read_IntegerVectorOption_subtag_for_job( ivectopt, full_key, val, vals, job_options );
			} else if ( dynamic_cast< PathVectorOptionKey const * > (&opt) ) {
				PathVectorOptionKey const & pathvectopt( static_cast< PathVectorOptionKey const & > (opt) );
				read_PathVectorOption_subtag_for_job( pathvectopt, full_key, val, vals, job_options );
			} else if ( dynamic_cast< RealVectorOptionKey const * > (&opt) ) {
				RealVectorOptionKey const & rvectopt( static_cast< RealVectorOptionKey const & > (opt) );
				read_RealVectorOption_subtag_for_job( rvectopt, full_key, val, vals, job_options );
			} else if ( dynamic_cast< StringVectorOptionKey const * > (&opt) ) {
				StringVectorOptionKey const & stvectopt( static_cast< StringVectorOptionKey const & > (opt) );
				read_StringVectorOption_subtag_for_job( stvectopt, full_key, val, vals, job_options );
			}
		}
	} else {
		std::ostringstream err;
		err << "Error: option '" << optname << "' corresponding to the full key '" << full_key << "' does not match any existing option in Rosetta.\n";
		err << "Thrown from JD2ResourceManagerJobInputter::parse_job_tag\n";
		throw EXCN_Msg_Exception( err.str() );
	}
}

void
JD2ResourceManagerJobInputter::read_BooleanOption_subtag_for_job(
	basic::options::BooleanOptionKey const & boolopt,
	std::string const & optname,
	std::string const & val,
	basic::resource_manager::JobOptionsOP job_options
)
{
	bool boolval;
	try {
		boolval = boost::lexical_cast< bool >( val );
	} catch ( boost::bad_lexical_cast const& ) {
		std::ostringstream err;
		err << "Error converting value '" << val << "' given for option '" << optname << "' to a boolean from within JD2ResourceManagerJobInputter::parse_job_tag\n Boolean options must be given either a '1' or a '0'";
		throw EXCN_Msg_Exception( err.str() );
	}
	job_options->add_option( boolopt, boolval );

}

void
JD2ResourceManagerJobInputter::read_FileOption_subtag_for_job(
	basic::options::FileOptionKey const & fileopt,
	std::string const &,
	std::string const & val,
	basic::resource_manager::JobOptionsOP job_options
)
{
	job_options->add_option( fileopt, val );
}


void
JD2ResourceManagerJobInputter::read_IntegerOption_subtag_for_job(
	basic::options::IntegerOptionKey const & intopt,
	std::string const & optname,
	std::string const & val,
	basic::resource_manager::JobOptionsOP job_options
)
{
	int intval;
	try {
		intval = boost::lexical_cast< int >( val );
	} catch ( boost::bad_lexical_cast const& ) {
		std::ostringstream err;
		err << "Error converting value '" << val << "' given for option '"
			<< optname << "' to an integer from within JD2ResourceManagerJobInputter::parse_job_tag\n";
		throw EXCN_Msg_Exception( err.str() );
	}
	job_options->add_option( intopt, intval );

}

void
JD2ResourceManagerJobInputter::read_PathOption_subtag_for_job(
	basic::options::PathOptionKey const & pathopt,
	std::string const &,
	std::string const & val,
	basic::resource_manager::JobOptionsOP job_options
)
{
	job_options->add_option( pathopt, val );
}


void
JD2ResourceManagerJobInputter::read_RealOption_subtag_for_job(
	basic::options::RealOptionKey const & realopt,
	std::string const & optname,
	std::string const & val,
	basic::resource_manager::JobOptionsOP job_options
)
{
	core::Real realval;
	try {
		realval = boost::lexical_cast< core::Real >( val );
	} catch ( boost::bad_lexical_cast const& ) {
		std::ostringstream err;
		err << "Error converting value '" << val << "' given for option '" << optname << "' to a floating point number from within JD2ResourceManagerJobInputter::parse_job_tag\n";
		throw EXCN_Msg_Exception( err.str() );
	}
	job_options->add_option( realopt, realval );
}

void
JD2ResourceManagerJobInputter::read_StringOption_subtag_for_job(
	basic::options::StringOptionKey const & stringopt,
	std::string const & ,
	std::string const & val,
	basic::resource_manager::JobOptionsOP job_options
)
{
	job_options->add_option( stringopt, val );
}


void
JD2ResourceManagerJobInputter::read_BooleanVectorOption_subtag_for_job(
	basic::options::BooleanVectorOptionKey const & boolvectopt,
	std::string const & optname,
	std::string const & val,
	utility::vector1< std::string > const & vals,
	basic::resource_manager::JobOptionsOP job_options
)
{
	utility::vector1< bool > boolvect( vals.size() );
	for ( Size ii = 1; ii <= vals.size(); ++ii ) {
		bool boolval;
		try {
			boolval = boost::lexical_cast< bool >( vals[ii] );
		} catch ( boost::bad_lexical_cast const& ) {
			std::ostringstream err;
			err << "Error converting value '" << vals[ii] << "', option # " << ii
				<< ", given for the comma-separated vector option '"
				<< optname << "' to a boolean\nfrom within JD2ResourceManagerJobInputter::parse_job_tag\n"
				<< "Original value string: '" << val << "'\n"
				<< "Boolean options must be given either a '1' or a '0'\n";
			throw EXCN_Msg_Exception( err.str() );
		}
		boolvect[ ii ] = boolval;
	}
	job_options->add_option( boolvectopt, boolvect );
}

void
JD2ResourceManagerJobInputter::read_FileVectorOption_subtag_for_job(
	basic::options::FileVectorOptionKey const & filevectopt,
	std::string const & ,
	std::string const & ,
	utility::vector1< std::string > const & vals,
	basic::resource_manager::JobOptionsOP job_options
)
{
	job_options->add_option( filevectopt, vals );
}

void
JD2ResourceManagerJobInputter::read_IntegerVectorOption_subtag_for_job(
	basic::options::IntegerVectorOptionKey const & intvectopt,
	std::string const & optname,
	std::string const & val,
	utility::vector1< std::string > const & vals,
	basic::resource_manager::JobOptionsOP job_options
)
{
	utility::vector1< int > intvect( vals.size() );
	for ( Size ii = 1; ii <= vals.size(); ++ii ) {
		bool intval;
		try {
			intval = boost::lexical_cast< int >( vals[ii] );
		} catch ( boost::bad_lexical_cast const& ) {
			std::ostringstream err;
			err << "Error converting value '" << vals[ii] << "', option # " << ii
				<< ", given for the comma-separated vector option '"
				<< optname << "' to an integer\nfrom within JD2ResourceManagerJobInputter::parse_job_tag\n"
				<< "Original value string: '" << val << "'\n";
			throw EXCN_Msg_Exception( err.str() );
		}
		intvect[ ii ] = intval;
	}
	job_options->add_option( intvectopt, intvect );
}

void
JD2ResourceManagerJobInputter::read_PathVectorOption_subtag_for_job(
	basic::options::PathVectorOptionKey const & pathvectopt,
	std::string const & ,
	std::string const & ,
	utility::vector1< std::string > const & vals,
	basic::resource_manager::JobOptionsOP job_options
)
{
	job_options->add_option( pathvectopt, vals );
}

void
JD2ResourceManagerJobInputter::read_RealVectorOption_subtag_for_job(
	basic::options::RealVectorOptionKey const & realvectopt,
	std::string const & optname,
	std::string const & val,
	utility::vector1< std::string > const & vals,
	basic::resource_manager::JobOptionsOP job_options
)
{
	utility::vector1< core::Real > realvect( vals.size() );
	for ( Size ii = 1; ii <= vals.size(); ++ii ) {
		core::Real realval;
		try {
			realval = boost::lexical_cast< core::Real >( vals[ii] );
		} catch ( boost::bad_lexical_cast const& ) {
			std::ostringstream err;
			err << "Error converting value '" << vals[ii] << "', option # " << ii << ", given for the comma-separated vector option '"
				<< optname << "' to a boolean\nfrom within JD2ResourceManagerJobInputter::parse_job_tag\n"
				<< "Original value string: '" << val << "'\n";
			throw EXCN_Msg_Exception( err.str() );
		}
		realvect[ ii ] = realval;
	}
	job_options->add_option( realvectopt, realvect );
}

void
JD2ResourceManagerJobInputter::read_StringVectorOption_subtag_for_job(
	basic::options::StringVectorOptionKey const & strinvectopt,
	std::string const & ,
	std::string const & ,
	utility::vector1< std::string > const & vals,
	basic::resource_manager::JobOptionsOP job_options
)
{
	job_options->add_option( strinvectopt, vals );
}


void
JD2ResourceManagerJobInputter::read_ResidueType_for_subtag(
	utility::tag::TagCOP options_tag,
	std::map< std::string, std::string > & resources_for_job
){
	std::string rname;

	std::ostringstream err;
	if ( !options_tag->hasOption("resource_tag") ) {
		err << "you must specify the resource_tag option when using a ResidueType tag";
		err <<std::endl;
		throw EXCN_Msg_Exception( err.str() );
	}
	for ( auto const & opt_iter : options_tag->getOptions() ) {
		if ( opt_iter.first == "resource_tag" ) {
			rname = opt_iter.second;
			resources_for_job["residue"] = rname;
		}
	}


}

/// @details Read the <Data/> subtag of the <Job/>
///
///   Form 1: <Data desc="resource_description" resource="resource_tag"/>
///      Precondition requirements:
///         has_resource_configuration(resource_tag)
///      Postcondition guarentees:
///         resources_for_job["resource_description"] == "resource_tag"
///
void
JD2ResourceManagerJobInputter::read_Data_for_subtag(
	utility::tag::TagCOP data_tag,
	std::string const & jobname,
	std::string & input_tag,
	std::map< std::string, std::string > & resources_for_job
)
{
	JD2ResourceManager * jd2rm( JD2ResourceManager::get_jd2_resource_manager_instance());

	bool desc_found( false );
	bool resource_found( false );
	bool pdb_found( false );
	bool local_startstruct_found( false );
	std::string desc;
	std::string rname;
	std::string locator = "";
	for ( auto const & opt_iter : data_tag->getOptions() ) {
		if ( opt_iter.first == "desc" ) {
			if ( opt_iter.second == "startstruct" ) {
				local_startstruct_found = true;
			}
			desc_found = true;
			desc = opt_iter.second;
		} else if ( opt_iter.first == "resource_tag" ) {
			resource_found = true;
			rname = opt_iter.second;
		} else if ( opt_iter.first == "pdb" ) {
			pdb_found = true;
			rname = opt_iter.second;
		} else if ( opt_iter.first == "locator" ) {
			locator = opt_iter.second;
		}
	}

	/// now let's make sure the provided data is consistent and complete
	bool in_error = false;
	std::ostringstream err;
	if ( ! desc_found ) {
		in_error = true;
		err << "Failed to find a data description (desc) amongst the options pairs listed reading a 'Data' tag in a Job tag.\n A desc option must always be given.\n";
	}
	if ( rname.size() == 0 ) {
		in_error = true;
		err << "Failed to find a resource name or a pdb name reading a 'Data' tag in a Job tag.  Either a 'resource_tag' or a 'pdb' option must be provided.\n";
	}
	if ( resource_found && pdb_found ) {
		in_error = true;
		err << "Error: Both a 'resource_tag' and a 'pdb' tag were found for a 'Data' tag in the Job tag." << std::endl;
	} else if ( resource_found && ! jd2rm->has_resource_configuration( rname ) ) {
		in_error = true;
		err << "Error: In Data subtag with descr='" << desc << "', the ResourceManager has no resource configuration for ResourceTag '" << rname << "'.\n";
	}
	if ( ! local_startstruct_found && pdb_found ) {
		in_error = true;
		err << "Error:  'pdb' tag given for a non-'startstruct' option in the 'Data' tag of a Job tag." << std::endl;
	}
	if ( resources_for_job.find( desc ) != resources_for_job.end() ) {
		in_error = true;
		err << "Error: description '" << desc << "' appears twice";
		err << ".\nFirst specified as '" << resources_for_job[ desc ] << "', and now as '" << rname << "'.";
		err << "A 'desc' may only be given once\n";
	}

	if ( in_error ) {
		err << "Problem encountered ";
		if ( jobname.size() != 0 ) {
			err << "for job named '" << jobname << "'";
		} else if ( desc != "startstruct" && resources_for_job.find("startstruct") != resources_for_job.end() ) {
			err << "for job whose starstruct is given as '" << resources_for_job[ "startstruct" ] << "'";
		}
		err << ".\nOptions given:\n";
		for ( auto const & opt_iter : data_tag->getOptions() ) {
			err << "\t(" << opt_iter.first << ", " << opt_iter.second << ")\n";
		}
		err << "Thrown from protocols::jd2::JD2ResourceManagerJobInputter::parse_job_tag\n";
		throw EXCN_Msg_Exception( err.str() );
	}

	if ( local_startstruct_found && pdb_found ) {
		/// add a pose to the ResourceManager with a canonical name, and then add that pose
		/// as a resource for this job.
		/// Then, when the
		input_tag = rname;
		std::string pdb_resource_name = "pdb_resource_" + rname;
		if ( jd2rm->has_resource_configuration( pdb_resource_name ) ) {
			// TODO: make sure the other resource has the right type
		} else {
			tr << "Adding implicit resource '" << pdb_resource_name << "' for job whose startstruct is given as pdb='" << rname << "'." << std::endl;
			basic::resource_manager::ResourceConfiguration rconfig;
			rconfig.resource_tag         = pdb_resource_name;
			rconfig.locator_tag          = locator; // <--- default
			rconfig.locator_id           = rname;
			rconfig.loader_type          = "PoseFromPDB";
			rconfig.resource_options_tag = ""; // <-- default PoseFromPDB options
			jd2rm->add_resource_configuration( pdb_resource_name, rconfig );
		}
		resources_for_job[ desc ] = pdb_resource_name;
	} else {
		/// put this resource name in a list for later telling the resource manager
		/// that this job (whose name may not yet be resolved) requires a particular resource
		resources_for_job[ desc ] = rname;
	}

}

void
JD2ResourceManagerJobInputter::check_each_job_has_startstruct(
	Jobs const & jobs
) const {
	JD2ResourceManager * jd2rm(
		JD2ResourceManager::get_jd2_resource_manager_instance());

	for ( auto const & job : jobs ) {
		if ( !(jd2rm->has_resource_tag_by_job_tag("startstruct", job->input_tag())) ) {
			std::stringstream errmsg;
			errmsg
				<< "Error: Job '" << job->input_tag() << "' given without a 'startstruct'";
			throw EXCN_Msg_Exception( errmsg.str() );
		}
	}
}

JobInputterInputSource::Enum
JD2ResourceManagerJobInputter::input_source() const
{
	return JobInputterInputSource::RESOURCE_MANAGED_JOB;
}


} // namespace jd2
} // namespace protocols

