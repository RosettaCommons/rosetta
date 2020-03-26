// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmtf/ExtraDataEnumManager.cc
/// @brief Enum string/enum functions for pose extra data we will be storing/retrieving.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/io/mmtf/ExtraDataEnumManager.hh>
#include <core/io/mmtf/ExtraDataEnum.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.io.mmtf.ExtraDataEnumManager" );


namespace core {
namespace io {
namespace mmtf {

ExtraDataEnumManager::ExtraDataEnumManager():
	utility::VirtualBase()
{
	setup_data_names();
}

void
ExtraDataEnumManager::setup_data_names(){
	string_to_enum_["pose_cache_string_data"] = pose_cache_string_data;
	string_to_enum_["pose_cache_real_data"] = pose_cache_real_data;
	string_to_enum_["pdb_comments"] = pdb_comments;
	string_to_enum_["simple_metric_string_data"] = simple_metric_string_data;
	string_to_enum_["simple_metric_real_data"] = simple_metric_real_data;
	string_to_enum_["simple_metric_composite_string_data"] = simple_metric_composite_string_data;
	string_to_enum_["simple_metric_composite_real_data"] = simple_metric_composite_real_data;
	string_to_enum_["simple_metric_per_residue_string_data"] = simple_metric_per_residue_string_data;
	string_to_enum_["simple_metric_per_residue_real_data"] = simple_metric_per_residue_real_data;
	string_to_enum_["simple_metric_per_residue_string_output"] = simple_metric_per_residue_string_output;
	string_to_enum_["simple_metric_per_residue_real_output"] = simple_metric_per_residue_real_output;


	if ( string_to_enum_.size() != ExtraDataEnum_total ) {
		utility_exit_with_message("ExtraDataEnumManager is missing at lease one enum for string_to_enum conversion");
	}

	enum_to_string_.resize( ExtraDataEnum_total );
	for ( auto const & enum_pair : string_to_enum_ ) {
		enum_to_string_[ enum_pair.second ] = enum_pair.first;
	}
}

ExtraDataEnum
ExtraDataEnumManager::string_to_enum(std::string const & data_name) const {
	if ( ! is_data_type(data_name) ) {
		utility_exit_with_message("Extra Data type not not present " + data_name);
	}

	return string_to_enum_.at(data_name);
}

std::string
ExtraDataEnumManager::enum_to_string( core::io::mmtf::ExtraDataEnum const data_name) const{
	return enum_to_string_[data_name];
}


bool
ExtraDataEnumManager::is_data_type(const std::string & data_name) const{
	return string_to_enum_.count(data_name);
}


ExtraDataEnumManager::~ExtraDataEnumManager(){}

ExtraDataEnumManager::ExtraDataEnumManager( ExtraDataEnumManager const & )  = default;



} //core
} //io
} //mmtf






