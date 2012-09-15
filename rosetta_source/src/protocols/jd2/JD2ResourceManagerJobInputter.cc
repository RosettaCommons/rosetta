// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/resource_manager/planner/JD2ResourceManagerJobInputter.cc
/// @brief
/// @author

// Unit headers
#include <protocols/jd2/JD2ResourceManagerJobInputter.hh>
#include <protocols/jd2/JD2ResourceManagerJobInputterCreator.hh>

// Package headers
#include <protocols/jd2/JD2ResourceManager.hh>
#include <protocols/jd2/Job.hh>
#include <basic/resource_manager/JobOptions.hh>
#include <protocols/jd2/InnerJob.hh>

// Project headers
#include <core/pose/Pose.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>

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

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>

// Boost headers
#include <boost/lexical_cast.hpp>

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


static basic::Tracer tr("protocols.jd2.JD2ResourceManagerJobInputter");



//CREATOR SECTION
std::string
JD2ResourceManagerJobInputterCreator::keyname() const
{
	return "JD2ResourceManagerJobInputter";
}

protocols::jd2::JobInputterOP
JD2ResourceManagerJobInputterCreator::create_JobInputter() const {
	return new JD2ResourceManagerJobInputter;
}



JD2ResourceManagerJobInputter::~JD2ResourceManagerJobInputter() {}

void
JD2ResourceManagerJobInputter::pose_from_job(
	Pose & pose,
	JobOP job )
{
	tr.Debug << "JD2ResourceManagerJobInputter::pose_from_job" << endl;

	if ( !job->inner_job()->get_pose() ) {
		tr.Debug
			<< "Retrieving pose from ResourceManager (tag = "
			<< job->inner_job()->input_tag() << ")" << endl;
		pose.clear();
		JD2ResourceManager * jd2_resource_manager(
			JD2ResourceManager::get_jd2_resource_manager_instance());
		ResourceOP resource;

		//Check to see if we have a Residue resource, if so load it into the chemical manager if it hasn't already been loaded
		if(jd2_resource_manager->has_resource_tag_by_job_tag("residue",job->inner_job()->input_tag()))
		{
			ResourceOP residue_resource = jd2_resource_manager->get_resource_by_job_tag("residue",job->inner_job()->input_tag());

			core::chemical::ResidueTypeOP new_residue(dynamic_cast<core::chemical::ResidueType *>(residue_resource()));
			std::string type_set_name(new_residue->residue_type_set().name());
			if(!core::chemical::ChemicalManager::get_instance()->residue_type_set(type_set_name)->has_name(new_residue->name()))
			{
				tr << "loading residue " << new_residue->name() << " into " << type_set_name <<" residue_type_set" <<std::endl;
				core::chemical::ChemicalManager::get_instance()->nonconst_residue_type_set(type_set_name).add_residue_type(new_residue);
			}
		}

		try {
			tr << "Loading startstruct " << jd2_resource_manager->find_resource_tag_by_job_tag( "startstruct", job->inner_job()->input_tag() ) << " for job " <<
				job->inner_job()->input_tag() <<std::endl;
			 resource = jd2_resource_manager->get_resource_by_job_tag(
			"startstruct", job->inner_job()->input_tag() );
		} catch ( utility::excn::EXCN_Msg_Exception const & e ) {
			std::ostringstream err;
			err << e.msg() << std::endl;
			err << "Failed to access 'startstruct' resource from the JD2ResourceManager for job '";
			err << job->inner_job()->input_tag() << "' with nstruct index " << job->nstruct_index();
			err << "\n" << "Exception caught and re-thrown from JD2ResourceManagerJobInputter::pose_from_job\n";
			throw utility::excn::EXCN_Msg_Exception( err.str() );
		}

		PoseOP resource_pose( dynamic_cast< Pose * > ( resource() ) );

		/// make sure the resource that we requested from the resource manager is in fact a pose.
		if ( ! resource_pose ) {
			std::ostringstream err;
			err
				<< "Error trying to access starting structure (startstruct) "
				<< "for job with tag '";
			err << job->inner_job()->input_tag() << "'.";
			if ( jd2_resource_manager->has_resource_tag_by_job_tag(
					"startstruct", job->inner_job()->input_tag() ) ) {
				ResourceTag rt =
					jd2_resource_manager->find_resource_tag_by_job_tag(
					"startstruct", job->inner_job()->input_tag() );
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
JD2ResourceManagerJobInputter::fill_jobs( Jobs & jobs )
{
	tr << "Initializing jobs from resource definition files" << std::endl;
	utility::vector1< utility::file::FileName > resource_definition_files
		= basic::options::option[ basic::options::OptionKeys::jd2::resource_definition_files ];
	for ( core::Size ii = 1; ii <= resource_definition_files.size(); ++ii ) {

		tr.Debug << "Reading from resource definition file " << resource_definition_files[ii] << std::endl;
		utility::io::izstream instream( resource_definition_files[ii]().c_str() );
		fill_jobs_from_stream( instream, jobs );
	}
}

void JD2ResourceManagerJobInputter::cleanup_input_after_job_completion(JobOP current_job)
{
	JD2ResourceManager * jd2_resource_manager(
				JD2ResourceManager::get_jd2_resource_manager_instance());
	std::string job_tag(current_job->inner_job()->input_tag());
	jd2_resource_manager->mark_job_tag_as_complete(job_tag);
	std::list<ResourceTag> resources_for_job(jd2_resource_manager->get_resource_tags_for_job_tag(job_tag));

	std::list<ResourceTag>::iterator resource_list_it = resources_for_job.begin();
	for(; resource_list_it != resources_for_job.end(); ++resource_list_it)
	{
		core::Size job_count(jd2_resource_manager->get_count_of_jobs_associated_with_resource_tag(*resource_list_it));
		if(job_count == 0)
		{
			//This might be a ResidueType.  if it is, we should delete the residue from the resource from the ChemicalManager
			ResourceOP current_residue = jd2_resource_manager->find_resource(*resource_list_it);
			core::chemical::ResidueTypeOP new_residue(dynamic_cast<core::chemical::ResidueType *>(current_residue()));

			if(new_residue)
			{
				std::string residue_type_set = new_residue->residue_type_set().name();
				core::chemical::ChemicalManager::get_instance()->
					ChemicalManager::nonconst_residue_type_set(residue_type_set).remove_residue_type(new_residue->name());
			}

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

	utility::tag::TagPtr resource_tags = utility::tag::Tag::create( instream );

	for ( utility::tag::Tag::tags_t::const_iterator
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
	utility::tag::TagPtr jobs_tags,
	Jobs & jobs
)
{
	for ( utility::tag::Tag::tags_t::const_iterator
			tag_iter = jobs_tags->getTags().begin(),
			tag_iter_end = jobs_tags->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter ) {
		std::string const & tagname = (*tag_iter)->getName();
		if ( tagname == "Job" ) {
			parse_job_tag( *tag_iter, jobs );
		} else {
			std::ostringstream err;
			err << "Error parsing jobs tags in JD2ResourceManagerJobInputter: unrecognized tag '" << tagname << "'";
			EXCN_Msg_Exception( err.str() );
		}
	}
}

void
JD2ResourceManagerJobInputter::parse_job_tag(
	utility::tag::TagPtr jobs_tag,
	Jobs & jobs
)
{
  JD2ResourceManager * jd2rm( JD2ResourceManager::get_jd2_resource_manager_instance());

	//InnerJobOP inner_job = new InnerJob;
	bool startstruct_found( false );
	std::string jobname;
	std::string input_tag; // <--- pdb name
	std::map< std::string, std::string > resources_for_job;


	/// if the name is not given
	if ( jobs_tag->hasOption( "name" )) {
		jobname = jobs_tag->getOption< std::string >( "name" );
	}

	JobOptionsOP job_options = new JobOptions;
	for ( utility::tag::Tag::tags_t::const_iterator
			tag_iter = jobs_tag->getTags().begin(),
			tag_iter_end = jobs_tag->getTags().end();
			tag_iter != tag_iter_end; ++tag_iter ) {
		std::string const & tagname = (*tag_iter)->getName();
		if ( tagname == "Option" ) {
			read_Option_subtag_for_job( *tag_iter, job_options );
		} else if ( tagname == "Data" ) {
			read_Data_for_subtag( *tag_iter, jobname, input_tag, startstruct_found, resources_for_job );
		} // are there other kinds of tags?
	}

	if ( ! startstruct_found ) {
		std::string errmsg( "Error: Job given without a 'startstruct'" );
		throw EXCN_Msg_Exception( errmsg );
	}

	if ( jobname.size() == 0 ) {
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

	int nstruct( 1 );
	if ( jobs_tag->hasOption( "nstruct" )) {
		std::string nstruct_string( jobs_tag->getOption< std::string >( "nstruct" ));
		try {
			nstruct = boost::lexical_cast< int >( nstruct_string );
		} catch ( boost::bad_lexical_cast const& ) {
			std::ostringstream err;
			err << "Error converting value '" << nstruct_string << "' given for job '" << jobname << "' for the nstruct option of the Job tag from JD2ResourceManagerJobInputter::parse_job_tag\n";
			throw EXCN_Msg_Exception( err.str() );
		}
	}

	/// OK: let's inform the JD2ResourceManager about all the resources required by this job
	for ( std::map< std::string, std::string >::const_iterator
			iter = resources_for_job.begin(), iter_end = resources_for_job.end();
			iter != iter_end; ++iter ) {
		jd2rm->add_resource_tag_by_job_tag( iter->first, jobname, iter->second );
	}

	InnerJobOP inner_job = new InnerJob( jobname, nstruct );
	JD2ResourceManager::get_jd2_resource_manager_instance()->add_job_options(
		jobname, job_options);
	for ( int ii = 1; ii <= nstruct; ++ii ) {
		jobs.push_back( new Job( inner_job, ii ));
	}

}

void
JD2ResourceManagerJobInputter::read_Option_subtag_for_job(
	utility::tag::TagPtr options_tag,
	JobOptionsOP job_options
)
{
	using namespace basic::options;
	for ( utility::tag::Tag::options_t::const_iterator
			opt_iter = options_tag->getOptions().begin(),
			opt_iter_end = options_tag->getOptions().end();
			opt_iter != opt_iter_end; ++opt_iter ) {
		std::string const & optname = opt_iter->first;
		std::string const & val = opt_iter->second;
		if ( basic::options::OptionKeys::has( optname ) ) {
			OptionKey const & opt( basic::options::OptionKeys::key( optname ));
			if ( opt.scalar() ) {
				// scalar options
				if ( dynamic_cast< BooleanOptionKey const * > (&opt) ) {
					BooleanOptionKey const & boolopt( static_cast< BooleanOptionKey const & > (opt) );
					read_BooleanOption_subtag_for_job( boolopt, optname, val, job_options );
				} else if ( dynamic_cast< FileOptionKey const * > (&opt) ) {
					FileOptionKey const & fileopt( static_cast< FileOptionKey const & > (opt) );
					read_FileOption_subtag_for_job( fileopt, optname, val, job_options );
				} else if ( dynamic_cast< IntegerOptionKey const * > (&opt) ) {
					IntegerOptionKey const & iopt( static_cast< IntegerOptionKey const & > (opt) );
					read_IntegerOption_subtag_for_job( iopt, optname, val, job_options );
				} else if ( dynamic_cast< PathOptionKey const * > (&opt) ) {
					PathOptionKey const & pathopt( static_cast< PathOptionKey const & > (opt) );
					read_PathOption_subtag_for_job( pathopt, optname, val, job_options );
				} else if ( dynamic_cast< RealOptionKey const * > (&opt) ) {
					RealOptionKey const & ropt( static_cast< RealOptionKey const & > (opt) );
					read_RealOption_subtag_for_job( ropt, optname, val, job_options );
				} else if ( dynamic_cast< StringOptionKey const * > (&opt) ) {
					StringOptionKey const & stopt( static_cast< StringOptionKey const & > (opt) );
					read_StringOption_subtag_for_job( stopt, optname, val, job_options );
				}

			} else {
				/// vector option
				utility::vector1< std::string > vals = utility::string_split( val, ',' );
				if ( dynamic_cast< BooleanVectorOptionKey const * > (&opt) ) {
					BooleanVectorOptionKey const & boolvectopt( static_cast< BooleanVectorOptionKey const & > (opt) );
					read_BooleanVectorOption_subtag_for_job( boolvectopt, optname, val, vals, job_options );
				} else if ( dynamic_cast< FileVectorOptionKey const * > (&opt) ) {
					FileVectorOptionKey const & filevectopt( static_cast< FileVectorOptionKey const & > (opt) );
					read_FileVectorOption_subtag_for_job( filevectopt, optname, val, vals, job_options );
				} else if ( dynamic_cast< IntegerVectorOptionKey const * > (&opt) ) {
					IntegerVectorOptionKey const & ivectopt( static_cast< IntegerVectorOptionKey const & > (opt) );
					read_IntegerVectorOption_subtag_for_job( ivectopt, optname, val, vals, job_options );
				} else if ( dynamic_cast< PathVectorOptionKey const * > (&opt) ) {
					PathVectorOptionKey const & pathvectopt( static_cast< PathVectorOptionKey const & > (opt) );
					read_PathVectorOption_subtag_for_job( pathvectopt, optname, val, vals, job_options );
				} else if ( dynamic_cast< RealVectorOptionKey const * > (&opt) ) {
					RealVectorOptionKey const & rvectopt( static_cast< RealVectorOptionKey const & > (opt) );
					read_RealVectorOption_subtag_for_job( rvectopt, optname, val, vals, job_options );
				} else if ( dynamic_cast< StringVectorOptionKey const * > (&opt) ) {
					StringVectorOptionKey const & stvectopt( static_cast< StringVectorOptionKey const & > (opt) );
					read_StringVectorOption_subtag_for_job( stvectopt, optname, val, vals, job_options );
				}
			}
		} else {
			std::ostringstream err;
			err << "Error: option '" << optname << "' does not match any existing option in Rosetta.\n";
			err << "Options must be fully namespaced (e.g. 'ex1' will not work, but 'packing:ex1' would work.\n";
			err << "Please remember to use only one colon when giving options.\n";
			err << "Thrown from JD2ResourceManagerJobInputter::parse_job_tag\n";
			throw EXCN_Msg_Exception( err.str() );
		}
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


///@details Read the <Data/> subtag of the <Job/>
///
///   Form 1: <Data desc="resource_description" resource="resource_tag"/>
///      Precondition requirements:
///         has_resource_configuration(resource_tag)
///      Postcondition guarentees:
///         resources_for_job["resource_description"] == "resource_tag"
///
void
JD2ResourceManagerJobInputter::read_Data_for_subtag(
	utility::tag::TagPtr data_tag,
	std::string const & jobname,
	std::string & input_tag,
	bool & startstruct_found,
	std::map< std::string, std::string > & resources_for_job
)
{
  JD2ResourceManager * jd2rm( JD2ResourceManager::get_jd2_resource_manager_instance());

	bool desc_found( false );
	bool resource_found( false );
	bool local_startstruct_found( false );
	bool pdb_found( false );
	std::string desc;
	std::string rname;
	for ( utility::tag::Tag::options_t::const_iterator
			opt_iter = data_tag->getOptions().begin(),
			opt_iter_end = data_tag->getOptions().end();
			opt_iter != opt_iter_end; ++opt_iter ) {
		if ( opt_iter->first == "desc" ) {
			if ( opt_iter->second == "startstruct" ) {
				startstruct_found = true;
				local_startstruct_found = true;
			}
			desc_found = true;
			desc = opt_iter->second;
		} else if ( opt_iter->first == "resource_tag" ) {
			resource_found = true;
			rname = opt_iter->second;
		} else if ( opt_iter->first == "pdb" ) {
			pdb_found = true;
			rname = opt_iter->second;
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
	} else if ( resource_found && ! jd2rm->has_resource_configuration( rname )) {
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
		} else if ( desc != "startstruct" && startstruct_found ) {
			err << "for job whose starstruct is given as '" << resources_for_job[ "startstruct" ] << "'";
		}
		err << ".\nOptions given:\n";
		for ( utility::tag::Tag::options_t::const_iterator
				opt_iter = data_tag->getOptions().begin(),
				opt_iter_end = data_tag->getOptions().end();
				opt_iter != opt_iter_end; ++opt_iter ) {
			err << "\t(" << opt_iter->first << ", " << opt_iter->second << ")\n";
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
			tr << "Adding implicit resource '" << pdb_resource_name << "' for job whose startstruct is given as pdb='" << rname << "'.";
			basic::resource_manager::ResourceConfiguration rconfig;
			rconfig.resource_tag         = pdb_resource_name;
			rconfig.locator_tag          = ""; // <--- default
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

JobInputterInputSource::Enum
JD2ResourceManagerJobInputter::input_source() const
{
	return JobInputterInputSource::RESOURCE_MANAGED_JOB;
}


} // namespace jd2
} // namespace protocols

