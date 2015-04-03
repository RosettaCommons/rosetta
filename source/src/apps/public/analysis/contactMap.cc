// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file contactMap.cc
/// @brief simple application for creating a contact map from multiple models of the same protein.
/// @author Joerg Schaarschmidt

//include block
// C++ headers
#include <string>
#include <sstream>

// libRosetta headers

#include <devel/init.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <protocols/filters/Filter.hh>
//#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/contact_map/ContactMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/contactMap.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

typedef protocols::contact_map::ContactMap ContactMap;


// declare variables that need to be accessed outside the main routine
static thread_local basic::Tracer tr( "contactMap" );

// Method to process region definition string into a vector of initialized ContactMap movers
utility::vector1<ContactMap> processRegions(std::string region_def, core::pose::Pose const & pose) {

	utility::vector1<ContactMap> maps;
	std::string current_region = "";

	// parameters for parse_my_tag routine
	protocols::moves::Movers_map movers;
	protocols::filters::Filters_map filters;
	basic::datacache::DataMap data;

	if(region_def == ""){
		std::ostringstream oss;
		oss << "1-" << pose.n_residue();
		region_def = oss.str();
	}
	// split region_def into single regions and process  each separately
	while(region_def != ""){
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace utility::tag;
		TagOP tag( new Tag() );
		core::Size pos = region_def.find(',');
		if (pos == std::string::npos){
			current_region = region_def;
			region_def = "";
		}else{
			current_region = region_def.substr(0, pos);
			region_def = region_def.substr(pos+1);
		}
		pos=current_region.find(':');
		// set regionX and ligand tags corresponding to region definiton
		if (pos == std::string::npos){
			tag->setOption("region1", current_region);
		}else if (current_region.find("ligand=")==std::string::npos){
			tag->setOption("region1", current_region.substr(0, pos));
			tag->setOption("region2", current_region.substr(pos+1));
		}else{
			if(current_region.find("ligand=")>pos){
				tag->setOption("region1", current_region.substr(0, pos));
				tag->setOption("ligand", current_region.substr(pos+8));
			}else{
				tag->setOption("region1", current_region.substr(pos+1));
				tag->setOption("ligand", current_region.substr(7,pos-7));
			}
		}
		// instantiate ContactMap and call parse_my_tag routine
		protocols::contact_map::ContactMap cm= protocols::contact_map::ContactMap();
		tag->setOption("distance_cutoff", option[contactMap::distance_cutoff]());
		tag->setOption("models_per_file", 0);
		if(option[contactMap::row_format]()) tag->setOption("row_format", "true");
		if(option[contactMap::distance_matrix]()) tag->setOption("distance_matrix", "true");
		cm.parse_my_tag( tag, data,	filters , movers, pose);
		// set output prefix for region
		std::string output_prefix(option[contactMap::prefix]());
		output_prefix +=current_region;
		cm.set_output_prefix(output_prefix);
		// add ContactMap to maps vector
		maps.push_back(cm);
	}
	return maps;
}

// main function
int main(int argc, char* argv[]) {

	try {

	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	// define protocol specific options
	OPT(contactMap::prefix);
	OPT(contactMap::distance_cutoff);
	OPT(contactMap::region_def);
	OPT(contactMap::row_format);
	OPT(contactMap::distance_matrix);

	// define relevant standard options
	OPT(in::path::database);
	OPT(in::file::s);
	OPT(in::file::l);
	OPT(in::file::silent);
	OPT(in::file::tags);
	OPT(in::file::silent_struct_type);
	OPT(in::file::residue_type_set);
	OPT(in::file::extra_res);
	OPT(in::file::extra_res_fa);
	OPT(in::file::extra_res_cen);

	// options, random initialization
	devel::init(argc, argv);

	// initialize variables
	core::chemical::ResidueTypeSetCOP rsd_set =
			ChemicalManager::get_instance()->residue_type_set(
					option[in::file::residue_type_set]());

	core::pose::Pose pose;

	// instantiate PoseInputStream based on input option
	MetaPoseInputStream input( streams_from_cmd_line() );

	// exit if no poses are supplied
	if (! input.has_another_pose()) {
		tr.Info << "No pose supplied - exiting without further action! "<< std::endl;
		return 0;
	}

	// use first pose to initialize ContactMaps
	input.fill_pose(pose, *rsd_set);

	// process region_def option and initialize ContactMaps
	utility::vector1<ContactMap> maps;
	maps = processRegions(option[contactMap::region_def](), pose);

	// reset PoseInputStream
	input.reset();

	// process poses
	while (input.has_another_pose()) {
		input.fill_pose(pose, *rsd_set);
		// loop over ContactMaps and call apply function
		for (core::Size i=1 ; i <= maps.size(); i++) {
			maps[i].apply(pose);
		}
	}
	// loop over ContactMaps and write each one to the corresponding file
	for (core::Size i=1 ; i <= maps.size(); i++) {
		maps[i].write_to_file();
	}

	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

} // main
