// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/core/chemical/sdf/MolData.cc
/// @author Sam DeLuca

#include <core/chemical/sdf/MolData.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <map>
#include <string>

// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector1.hh>

namespace core {
namespace chemical {
namespace sdf {
static basic::Tracer MolDataTracer("core.io.sdf.MolData");

MolData::MolData()
{

}

MolData::~MolData()
{

}

core::Size MolData::size()
{
	return mol_data_map_.size();
}

void MolData::clear()
{
	mol_data_map_.clear();
}

void MolData::parse_mol_data(utility::vector1<std::string> const & file_lines)
{
	bool inside_block = false;

	std::string data_name = "";
	std::string data = "";

	BOOST_FOREACH(std::string line, file_lines){
		if(line[0] == '>') //we've found a data header
		{
			data_name = line.substr(3,line.size());
			utility::trim(data_name,">");
			data.clear();
			inside_block = true;
		}else if(inside_block && (line[0] == '\n' || line.size() == 0 ))
		{
			//data_name = "";
			mol_data_map_.insert(std::pair<std::string,std::string>(data_name,data));
			inside_block = false;
		}else if(inside_block)
		{
			data.append(line+"\n");
		}
	}
}

std::string MolData::get_mol_data(std::string const & key) const
{
	std::map<std::string,std::string>::const_iterator data_it = mol_data_map_.find(key);
	if(data_it == mol_data_map_.end())
	{
		return "";
	}else
	{
		return data_it->second;
	}
}

utility::vector1<std::string> MolData::get_mol_data_string_vector(std::string const & key, char const & splitter) const
{
	std::map<std::string,std::string>::const_iterator data_it = mol_data_map_.find(key);
	if(data_it == mol_data_map_.end())
	{
		utility::vector1<std::string> data_vector;
		return data_vector;
	}else
	{
		std::string data = data_it->second;
		utility::vector1<std::string> data_vector(utility::string_split(data,splitter));
		return data_vector;
	}
}

void MolData::print() const
{
	std::map<std::string,std::string>::const_iterator data_it;
	for(data_it = mol_data_map_.begin(); data_it != mol_data_map_.end();++data_it)
	{
		MolDataTracer << "\""<< data_it->first <<"\"" << " : " <<data_it->second <<std::endl;
	}
}

}
}
}
