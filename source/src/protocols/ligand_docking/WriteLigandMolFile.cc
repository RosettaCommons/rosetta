// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/WriteLigandMolFile.cc
/// @brief  Ligand MolFile writer
/// @author Sam DeLuca

// MPI headers
#ifdef USEMPI
#include <mpi.h> //keep this first
#endif

#include <protocols/ligand_docking/WriteLigandMolFile.hh>
#include <protocols/ligand_docking/WriteLigandMolFileCreator.hh>
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/file/file_sys_util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <core/chemical/sdf/mol_writer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <core/pose/util.hh>
#include <map>
#include <basic/Tracer.hh>

namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer write_ligand_tracer( "protocols.ligand_docking.WriteLigandMolFile" );

std::string
WriteLigandMolFileCreator::keyname() const
{
	return WriteLigandMolFileCreator::mover_name();
}

protocols::moves::MoverOP
WriteLigandMolFileCreator::create_mover() const {
	return new WriteLigandMolFile;
}

std::string
WriteLigandMolFileCreator::mover_name()
{
	return "WriteLigandMolFile";
}

WriteLigandMolFile::WriteLigandMolFile() :
	Mover("WriteLigandMolFile"), chain_(""),directory_(""),prefix_(""),hash_file_names_(false)
{

}
WriteLigandMolFile::~WriteLigandMolFile()
{

}

WriteLigandMolFile::WriteLigandMolFile(WriteLigandMolFile const & that) :
	Mover("WriteLigandMolFile"),chain_(that.chain_),
    directory_(that.directory_),
    prefix_(that.prefix_),
    hash_file_names_(that.hash_file_names_)
{

}

protocols::moves::MoverOP WriteLigandMolFile::clone() const
{
	return new WriteLigandMolFile(*this);
}

protocols::moves::MoverOP WriteLigandMolFile::fresh_instance() const
{
	return new WriteLigandMolFile;
}

std::string WriteLigandMolFile::get_name() const
{
	return "WriteLigandMolFile";
}

void WriteLigandMolFile::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{
	if(!tag->hasOption("chain"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("'WriteLigandMolFile' requires the option 'chain'");
	}

	if(!tag->hasOption("directory"))
	{
		throw utility::excn::EXCN_RosettaScriptsOption("'WriteLigandMolFile' requires the option 'directory'");
	}
    if(!tag->hasOption("prefix"))
    {
        throw utility::excn::EXCN_RosettaScriptsOption("'WriteLigandMolFile' requires the option 'prefix'");
    }

	std::string hash_status = tag->getOption<std::string>("hash_file_names","false");
	if(hash_status == "true")
	{
		hash_file_names_ = true;
	}else
	{
		hash_file_names_ = false;
	}

	chain_ = tag->getOption<std::string>("chain");
	directory_ = tag->getOption<std::string>("directory");
    prefix_ = tag->getOption<std::string>("prefix");

}

void WriteLigandMolFile::apply(core::pose::Pose & pose)
{

    
	jd2::JobDistributor* job_dist(jd2::JobDistributor::get_instance());
	jd2::JobOP current_job(job_dist->current_job());

    core::chemical::sdf::MolWriter mol_writer("V2000");
    
	jd2::Job::StringStringPairs string_string_data(job_dist->current_job()->get_string_string_pairs());
	jd2::Job::StringRealPairs string_real_data(job_dist->current_job()->get_string_real_pairs());
    
	std::map<std::string,std::string> job_data;
    
	for(jd2::Job::StringStringPairs::const_iterator it = string_string_data.begin(); it != string_string_data.end();++it )
	{
		job_data.insert(std::make_pair(it->first,it->second));
	}
    
	for(jd2::Job::StringRealPairs::const_iterator it = string_real_data.begin(); it != string_real_data.end();++it )
	{
		job_data.insert(std::make_pair(it->first,utility::to_string(it->second)));
	}
    
	mol_writer.set_job_data(job_data);
    
	core::Size chain_id = core::pose::get_chain_id_from_chain(chain_,pose);
	core::Size residue_index = pose.conformation().chain_begin(chain_id);
	core::conformation::ResidueCOP ligand_residue(pose.conformation().residue(residue_index).get_self_ptr());
    
    utility::file::create_directory(directory_);
    
#ifdef USEMPI
	int mpi_rank;
	MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);/* get current process id */
    std::string mpi_rank_string(utility::to_string(mpi_rank));
    std::string output_file = directory_+"/"+prefix_+"_"+mpi_rank_string+".sdf.gz";
#else
    std::string output_file = directory_+"/"+prefix_+".sdf.gz";
#endif
    utility::io::ozstream output;
    std::stringstream header;
    output.open_append_if_existed( output_file, header);
    
    write_ligand_tracer << "Writing ligand to " << output_file <<std::endl;
	mol_writer.output_residue(output,ligand_residue);
    
}

}
}

