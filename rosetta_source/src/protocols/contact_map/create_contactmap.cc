// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file create_contactmap.cc
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
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/ContactMap.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

typedef protocols::moves::ContactMap ContactMap;




// declare variables that need to be accessed outside the main routine
basic::Tracer tr("create_contactmap");
core::Real DISTANCE_CUTOFF=10.0;

// declare protocol specific options
OPT_1GRP_KEY(File, contactmap, prefix)
OPT_1GRP_KEY(Real, contactmap, distance_cutoff)
OPT_1GRP_KEY(Real, contactmap, energy_cutoff)
OPT_1GRP_KEY(String, contactmap, region_def)
OPT_1GRP_KEY(Boolean, contactmap, row_format)

// Method to process region definition string into a vector of initialized ContactMap movers
utility::vector1<ContactMap> processRegions(std::string region_def, core::pose::Pose const & pose) {

	utility::vector1<ContactMap> maps;
	std::string current_region = "";

	// parameters for parse_my_tag routine
	protocols::moves::Movers_map movers;
	protocols::filters::Filters_map filters;
	protocols::moves::DataMap data;

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
		TagPtr tag = TagPtr( new Tag() );;
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
		protocols::moves::ContactMap cm= protocols::moves::ContactMap();
		tag->setOption("distance_cutoff", DISTANCE_CUTOFF);
		tag->setOption("models_per_file", 0);
		if(option[contactmap::row_format]())tag->setOption("row_format", "true");
		cm.parse_my_tag( tag, data,	filters , movers, pose);
		// set output prefix for region
		std::string output_prefix(option[contactmap::prefix]());
		output_prefix +=current_region;
		cm.set_output_prefix(output_prefix);
		// add ContactMap to maps vector
		maps.push_back(cm);
	}
	return maps;
}

// main function
int main(int argc, char* argv[]) {
	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	// add and process protocol specific options


	NEW_OPT(contactmap::prefix, "Prefix of contactmap filename","contact_map_");
	NEW_OPT(contactmap::distance_cutoff,
							"Cutoff Backbone distance for two atoms to be considered interacting", DISTANCE_CUTOFF);
	NEW_OPT(contactmap::energy_cutoff,
							"Energy_Cutoff (percentage value - only affecting silent file input)", 1.0);
	option[contactmap::energy_cutoff].lower(0.0).upper(1.0);
	NEW_OPT(contactmap::region_def,
			"Region definition for comparison eg: 1-10:20-30,40-50,A:ligand=X", "");
	NEW_OPT(contactmap::row_format,
				"Flag whether to output in row instead of matrix format", false);

	// define relevant standard options
	OPT(in::path::database);
	OPT(in::file::s);
	OPT(in::file::silent);
	OPT(in::file::tags);
	OPT(in::file::silent_struct_type);
	OPT(in::file::residue_type_set);

	// options, random initialization
	devel::init(argc, argv);

	// usage declaration
	std::string usage("");
	usage
			+= "\n\nusage:  create_heatmap [options] -in::file::silent <silent_files>\n";
	usage += "\tTo see a list of other valid options, use the option -help.\n";

	// check if files to process have been specified
	if (!option[in::file::silent].user() && !option[in::file::s].user()) {
		std::cerr << usage << std::endl;
		std::exit(1);
	}

	// initialize variables
	utility::vector1<ContactMap> maps;
	DISTANCE_CUTOFF = option[contactmap::distance_cutoff]();

	core::chemical::ResidueTypeSetCAP rsd_set =
			ChemicalManager::get_instance()->residue_type_set(
					option[in::file::residue_type_set]());

	PoseInputStreamOP input;
	core::pose::Pose pose;

	// instantiate PoseInputStream based on input option
	if (option[in::file::silent].user()) {
		if (option[in::file::tags].user()) {
			input = new SilentFilePoseInputStream(option[in::file::silent](),
					option[in::file::tags](), option[contactmap::energy_cutoff]());
		} else {
			input = new SilentFilePoseInputStream(option[in::file::silent](),
					option[contactmap::energy_cutoff]());
		}
	} else if (option[in::file::s].user()) {
		input = new PDBPoseInputStream(option[in::file::s]());
	}

	// use first pose initialize ContactMaps
	input->fill_pose(pose, *rsd_set);

	// process region_def option and initialize ContactMaps
	maps = processRegions(option[contactmap::region_def](), pose);

	// reset PoseInputStream
	input->reset();

	// process poses
	while (input->has_another_pose()) {
		input->fill_pose(pose, *rsd_set);
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
} // main
