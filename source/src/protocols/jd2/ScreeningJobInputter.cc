// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/ScreeningJobInputter.cc
/// @brief  The ScreeningJobInputter
/// @detail The ScreeningJobInputter simplifies the process of inputting and handling input to a vHTS protocol.
/// The ScreeningJobInputter takes in a JSON file using the -in:file:screening_job_file option
/// This file is in the following format:
/// {
///  "jobs" :
///   [
///    {
///     "group_name" : "group_name_a",
///     "proteins" :
///      ["protein_a.pdb", "protein_b.pdb"],
///     "ligands" :
///      ["ligand_a.pdb", "ligand_b.pdb"]
///    }
///   ]
///  "params :
///   [
///      "aaa.param",
///      "aab.param"
//    ]
/// }
/// In this scheme, jobs will be created for all combinations of rpoteins and ligands, and the job will be tagged with
/// "group_name" using a job string/real pair, in the form "input_group_name group_name_a"
/// @author Sam DeLuca <samuel.l.deluca@vanderbilt.edu>

#include <protocols/jd2/ScreeningJobInputter.hh>
#include <protocols/jd2/ScreeningJobInputterCreator.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/InnerJob.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/io/izstream.hh>

#include <basic/Tracer.hh>
#include <algorithm>
namespace protocols {
namespace jd2 {

static THREAD_LOCAL basic::Tracer TR( "protocols.jd2.ScreeningJobInputter" );

//CREATOR SECTION
std::string
ScreeningJobInputterCreator::keyname() const
{
	return "ScreeningJobInputter";
}

protocols::jd2::JobInputterOP
ScreeningJobInputterCreator::create_JobInputter() const {
	return protocols::jd2::JobInputterOP( new ScreeningJobInputter );
}


ScreeningJobInputter::ScreeningJobInputter()
{

}

ScreeningJobInputter::~ScreeningJobInputter()
{

}

void ScreeningJobInputter::pose_from_job(core::pose::Pose & pose, JobOP job)
{
	TR << "ScreeningJobInputter::pose_from_job" << std::endl;

	if ( !job->inner_job()->get_pose() ) {
		TR << "filling pose from PDB " << job->input_tag() << std::endl;
		core::import_pose::pose_from_file( pose, job->input_tag() , core::import_pose::PDB_file);
		load_pose_into_job(pose, job);
	} else {
		TR << "filling pose from saved copy " << job->input_tag() << std::endl;
		pose = *(job->inner_job()->get_pose());
	}
}

void ScreeningJobInputter::fill_jobs(JobsContainer & jobs)
{
	std::string file_name(basic::options::option[basic::options::OptionKeys::in::file::screening_job_file]());
	utility::io::izstream data(file_name.c_str(),std::ifstream::in);
	if ( !data.good() ) {
		utility_exit_with_message("Unable to open Screening Job File: " + file_name);
	}
	core::Size const nstruct( get_nstruct() );

	utility::json_spirit::mValue job_json_data;
	utility::json_spirit::mObject job_object_data;
	try
{
		utility::json_spirit::read(data,job_json_data);
		job_object_data = job_json_data.get_obj();
	}catch(std::runtime_error &)
{
		throw utility::excn::EXCN_BadInput(
			"screening file " + file_name + "is incorrectly formatted. "
			"it must be a dict with two keys: 'params' containing a list of needed params, and "
			"'jobs' which must be a list of dicts. each dict should contain keys 'group_name', 'proteins' and 'ligands'");
	}


	utility::json_spirit::mArray param_group_data;
	try
{
		param_group_data = job_object_data["params"].get_array();
	}catch(std::runtime_error &)
{
		throw utility::excn::EXCN_BadInput("the screening file " + file_name + " does not contain a 'params' section");
	}

//parse params files and insert them into the chemical manager
	if ( param_group_data.size() > 0 ) {
		for ( core::Size i = 0; i < param_group_data.size(); ++i ) {
			std::string param_name = param_group_data[i].get_str();
			core::chemical::ChemicalManager::get_instance()->
				nonconst_residue_type_set(core::chemical::FA_STANDARD).add_custom_residue_type(param_name);

		}
	}

	utility::json_spirit::mArray job_group_data;
	try
{
		job_group_data = job_object_data["jobs"].get_array();
	}catch(std::runtime_error & )
{
		throw utility::excn::EXCN_BadInput("the screening file " + file_name + " does not contain a 'jobs' section");
	}
	for ( core::Size i = 0; i < job_group_data.size(); ++i ) {
		utility::json_spirit::mObject group_map(job_group_data[i].get_obj());

		std::string group_name;
		utility::json_spirit::mArray protein_path_data;
		utility::json_spirit::mArray ligand_path_data;

		//The exceptions thrown by json_spirit are incredibly generic and meaningless. we catch them and throw more useful ones
		try
{
			group_name = group_map["group_name"].get_str();
		}catch(std::runtime_error &)
{
			throw utility::excn::EXCN_BadInput("a group in screening file " + file_name + " does not contain the element 'group_name' or is misformatted");
		}

		utility::json_spirit::mArray startfrom_data;
		bool startfrom_present = false;
		if ( group_map.find("startfrom") != group_map.end() ) {
			startfrom_data = group_map["startfrom"].get_array();
			startfrom_present = true;
		}

		try
{
			protein_path_data = group_map["proteins"].get_array();
		}catch(std::runtime_error &)
{
			throw utility::excn::EXCN_BadInput("the group " + group_name +" in screening file " + file_name + " does not contain the element 'proteins' or is misformatted");
		}

		try
{
			ligand_path_data = group_map["ligands"].get_array();
		}catch(std::runtime_error &)
{
			throw utility::excn::EXCN_BadInput("the group " + group_name +" in screening file " + file_name + " does not contain the element 'ligands' or is misformatted");
		}

//If we specify a native structure, store it
		bool native_present = false;
		std::string native_string;
		if ( group_map.find("native") != group_map.end() ) {
			native_string = group_map["native"].get_str();
			native_present = true;
		}
		//Make a job for each combination of a protein and a ligand defined in the group
		for ( core::Size protein_path_index = 0; protein_path_index < protein_path_data.size(); ++protein_path_index ) {
			std::string protein_path(protein_path_data[protein_path_index].get_str());

			for ( core::Size ligand_path_index = 0; ligand_path_index < ligand_path_data.size(); ++ligand_path_index ) {
				std::string ligand_path(ligand_path_data[ligand_path_index].get_str());
				//multiple pdb paths combined with a space get concatenated into a single pose
				std::string input_tag(protein_path + " " + ligand_path);

				InnerJobOP ijob( new InnerJob( input_tag,nstruct) );
				for ( core::Size index(1); index <= nstruct; ++index ) {
					JobOP current_job( new Job(ijob,index) );
					//Tag the current job with the group name so that we can keep track of what is is
					current_job->add_string_string_pair("input_group_name",group_name);
					//Add the path to the native structure so that we can compute RMS values later
					//This has to be explicitly supported at the protocol level.  It is not equivilent
					//to -in:file:native
					if ( native_present ) {
						current_job->add_string_string_pair("native_path",native_string);
					}

					if ( startfrom_present ) {
						core::Real xcoord = startfrom_data[0].get_real();
						core::Real ycoord = startfrom_data[1].get_real();
						core::Real zcoord = startfrom_data[2].get_real();

						current_job->add_string_real_pair("start_x", xcoord);
						current_job->add_string_real_pair("start_y", ycoord);
						current_job->add_string_real_pair("start_z", zcoord);

					}
					jobs.push_back( current_job );
					TR << "pushing " << input_tag << " nstruct index " << index << std::endl;
				}//loop over nstruct
			}
		}

	}

	if ( basic::options::option[basic::options::OptionKeys::in::file::shuffle_screening_jobs]() ) {
		jobs.shuffle();
	}
}

JobInputterInputSource::Enum ScreeningJobInputter::input_source() const
{
	return JobInputterInputSource::SCREENING_FILE;

}

}
}
