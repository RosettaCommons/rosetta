// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com)
///

// Unit Headers
#include <protocols/ligand_docking/StartFrom.hh>
#include <protocols/ligand_docking/StartFromCreator.hh>

#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

//project headers
#include <numeric/random/random.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

// Utility Headers

#include <utility/tag/Tag.hh>
#include <utility/exit.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>

#include <utility/json_spirit/json_spirit_reader.h>

// Boost headers
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

#include <utility/excn/Exceptions.hh>
#include <fstream>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static basic::Tracer start_from_tracer("protocols.ligand_docking.ligand_options.Start_from", basic::t_debug);

std::string
StartFromCreator::keyname() const
{
	return StartFromCreator::mover_name();
}

protocols::moves::MoverOP
StartFromCreator::create_mover() const {
	return new StartFrom;
}

std::string
StartFromCreator::mover_name()
{
	return "StartFrom";
}

///@brief
StartFrom::StartFrom():
		//utility::pointer::ReferenceCount(),
		Mover("StartFrom"),
		chain_(""),
		starting_points_(),
		use_file_name_(false){}

StartFrom::StartFrom(StartFrom const & that):
		//utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		chain_(that.chain_),
		starting_points_(that.starting_points_),
		potential_starting_positions_(that.potential_starting_positions_),
		use_file_name_(that.use_file_name_)
{}

StartFrom::~StartFrom() {}

protocols::moves::MoverOP StartFrom::clone() const {
	return new StartFrom( *this );
}

protocols::moves::MoverOP StartFrom::fresh_instance() const {
	return new StartFrom;
}

std::string StartFrom::get_name() const{
	return "StartFrom";
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
StartFrom::parse_my_tag(
		utility::tag::TagCOP const tag,
		basic::datacache::DataMap & /*datamap*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "StartFrom" ){
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover requires chain tag");

	chain_ = tag->getOption<std::string>("chain");

	foreach(utility::tag::TagCOP child_tag, tag->getTags()){
		std::string name= child_tag->getName();
		if( name == "features"){
			std::cout << "found features tag with type '" << child_tag->getOption<std::string>("type") << "'" << std::endl;

		} else if( name == "Coordinates")
		{
			if ( ! child_tag->hasOption("x") ) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover Coordinates tag requires 'x' coordinates option");
			if ( ! child_tag->hasOption("y") ) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover Coordinates tag requires 'y' coordinates option");
			if ( ! child_tag->hasOption("z") ) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover Coordinates tag requires 'z' coordinates option");

			std::string pdb_tag = "default";
			if(child_tag->hasOption("pdb_tag"))
			{
				pdb_tag = child_tag->getOption<std::string>("pdb_tag");
			}

			core::Vector v(
					child_tag->getOption<core::Real>("x"),
					child_tag->getOption<core::Real>("y"),
					child_tag->getOption<core::Real>("z")
			);

			coords(v,pdb_tag);
		}else if(name == "File")
		{
			std::string location_id_mode = child_tag->getOption<std::string>("struct_identifier","hash");
			if(location_id_mode =="hash")
			{
				use_file_name_ = false;
			}else if(location_id_mode =="file")
			{
				use_file_name_ = true;
			}else
			{
				throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' Mover option struct_identifier can only be 'hash' or 'file'");
			}

			if(!child_tag->hasOption("filename")) throw utility::excn::EXCN_RosettaScriptsOption("'StartFrom' mover File tag requires 'filename' coordinates option");
			parse_startfrom_file(child_tag->getOption<std::string>("filename"));

		}
	}
}

void StartFrom::coords(core::Vector const & coords,std::string const & pdb_tag)
{
	std::map< std::string, utility::vector1<core::Vector> >::iterator start_point_it = starting_points_.find(pdb_tag);

	if(start_point_it == starting_points_.end())
	{
		utility::vector1<core::Vector> new_point_set;
		new_point_set.push_back(coords);
		starting_points_.insert(std::make_pair(pdb_tag,new_point_set));
	}else
	{
		start_point_it->second.push_back(coords);
	}
}

void StartFrom::chain(std::string const & chain)
{
	chain_ = chain;
}

void StartFrom::apply(core::pose::Pose & pose){
	assert(!starting_points_.empty() || !potential_starting_positions_.empty());
	int const starting_point_index= numeric::random::RG.random_range(1, starting_points_.size());

	if(!core::pose::has_chain(chain_,pose))
	{
		utility_exit_with_message("StartFrom mover cannot find the chain " +chain_+ " in the current pose.");
	}

	if(!starting_points_.empty())
	{
		std::string input_tag(protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag());
		std::list<std::string> component_tags(utility::split_to_list(input_tag));

		utility::vector1<core::Vector> centroid_points;

		for(std::list<std::string>::iterator tag_it = component_tags.begin(); tag_it != component_tags.end();++tag_it)
		{
			std::map< std::string, utility::vector1<core::Vector> >::iterator tag_points( starting_points_.find(*tag_it));
			if(tag_points != starting_points_.end())
			{
				centroid_points = tag_points->second;
				break;
			}else
			{
				tag_points = starting_points_.find("default");
				if(tag_points == starting_points_.end())
				{
					utility_exit_with_message("There are no default starting coordinates specified in the StartFrom mover, and none of the specified coordinates match the current tag");
				}
				centroid_points = tag_points->second;
				break;
			}
		}
		assert(centroid_points.size());

		core::Vector desired_centroid = centroid_points[starting_point_index];

		core::Size jump_id = core::pose::get_jump_id_from_chain(chain_, pose);
		move_ligand_to_desired_centroid(jump_id, desired_centroid, pose);
	}else if(!potential_starting_positions_.empty())
	{
		std::map<std::string,core::Vector >::iterator position_id;
		if(use_file_name_)
		{
			std::string job_tag(jd2::JobDistributor::get_instance()->current_job()->input_tag());
			utility::vector1<std::string> input_filenames(utility::split(job_tag));
			foreach(std::string filename, input_filenames)
			{
				utility::file::FileName file_data(filename);
				std::string base_name(file_data.base());
				position_id = potential_starting_positions_.find(base_name);
				if(position_id != potential_starting_positions_.end())
				{
					break;
				}
			}
		}else
		{
			std::string hash = core::pose::get_sha1_hash_excluding_chain(chain_[0],pose);
			position_id = potential_starting_positions_.find(hash);
		}
		if(position_id != potential_starting_positions_.end())
		{
			core::Size jump_id = core::pose::get_jump_id_from_chain(chain_, pose);
			move_ligand_to_desired_centroid(jump_id, position_id->second , pose);
		}else
		{
			start_from_tracer << "cannot find structure with position id and tag " <<
					protocols::jd2::JobDistributor::get_instance()->current_job()->input_tag() << std::endl;
			utility_exit_with_message("the current structure is not in the startfrom_file");
		}
	}else
	{
		utility_exit_with_message("You must specify either a Coordinates or a File tag in the StartFrom mover");
	}
}

void StartFrom::parse_startfrom_file(std::string filename)
{
    if(!utility::file::file_exists(filename))
    {
        utility_exit_with_message("cannot parse "+filename+" because it does not exist");
    }
    
	utility::io::izstream infile;
	infile.open(filename.c_str(),std::ifstream::in);
	utility::json_spirit::mValue startfrom_data;
    if(!utility::json_spirit::read(infile,startfrom_data))
    {
        infile.close();
        utility_exit_with_message("cannot parse JSON in file "+ filename);
    }

	infile.close();

	//The format is something like this:
	/*
	[
	    {
	        "input_tag" : "path/infile.pdb",
	        "file_name" : "infile.pdb",
	        "x" : 0.0020,
	        "y" : -0.004,
	        "z" : 0.0020,
	        "hash" : "aa2aff055d19bc32e483df7ff4ae08361a768931"
	    }
	]
	*/

	utility::json_spirit::mArray start_positions = startfrom_data.get_array();
	for(utility::json_spirit::mArray::iterator start_it = start_positions.begin(); start_it != start_positions.end();++start_it)
	{
		utility::json_spirit::mObject position_data(start_it->get_obj());

		std::string identifier;
		if(use_file_name_)
		{
			identifier = position_data["file_name"].get_str();
		}else
		{
			identifier = position_data["hash"].get_str();
		}

		core::Real x = position_data["x"].get_real();
		core::Real y = position_data["y"].get_real();
		core::Real z = position_data["z"].get_real();

		core::Vector coords(x,y,z);
		if(potential_starting_positions_.find(identifier) !=potential_starting_positions_.end() )
		{

			start_from_tracer << "WARNING: There is more than one entry in the startfrom_file with the hash " << identifier <<std::endl;
			//utility_exit_with_message("hashes in startfrom files must all be unique");
		}
		potential_starting_positions_[identifier] = coords;

	}

}

void
move_ligand_to_desired_centroid(
		core::Size const jump_id,
		core::Vector const desired_centroid,
		core::pose::Pose & pose
){
	core::Vector const ligand_centroid = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
	core::Vector const trans_vec = desired_centroid - ligand_centroid;
	core::Real const trans_len = trans_vec.length();
	if (trans_len > 1e-3) { // otherwise we get NaNs
		protocols::rigid::RigidBodyTransMover mover(pose, jump_id);
		mover.step_size(trans_len);
		mover.trans_axis(trans_vec);
		mover.apply(pose);
	}
}


} //namespace ligand_docking
} //namespace protocols
