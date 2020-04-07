// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file contactMap.cc
/// @brief simple application for creating a contact map from multiple models of the same protein.
/// @author Joerg Schaarschmidt

//include block
// C++ headers
#include <string>

// libRosetta headers

#include <devel/init.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>

#include <protocols/filters/Filter.fwd.hh>
//#include <protocols/moves/Mover.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/contact_map/ContactMap.hh>

#include <core/pose/Pose.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/contactMap.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

using ContactMap = protocols::contact_map::ContactMap;


// declare variables that need to be accessed outside the main routine
static basic::Tracer tr( "contactMap" );

// Method to process region definition string into a vector of initialized ContactMap movers
utility::vector1<ContactMap> processRegions(std::string region_def, core::pose::Pose const & pose) {

	utility::vector1<ContactMap> maps;
	std::string current_region = "";

	if ( region_def == "" ) {
		std::ostringstream oss;
		oss << "1-" << pose.size();
		region_def = oss.str();
	}
	// split region_def into single regions and process  each separately
	while ( region_def != "" ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		core::select::residue_selector::ResidueSelectorCOP region1, region2;
		core::Size pos = region_def.find(',');
		if ( pos == std::string::npos ) {
			current_region = region_def;
			region_def = "";
		} else {
			current_region = region_def.substr(0, pos);
			region_def = region_def.substr(pos+1);
		}

		// instantiate ContactMap and initialize
		protocols::contact_map::ContactMap cm;

		// set regionX and ligand tags corresponding to region definiton
		pos=current_region.find(':');
		if ( pos == std::string::npos ) {
			region1 = cm.parse_region_string( current_region );
			cm.set_regions( region1 );
		} else if ( current_region.find("ligand=")==std::string::npos ) {
			region1 = cm.parse_region_string( current_region.substr(0, pos) );
			region2 = cm.parse_region_string( current_region.substr(pos+1) );
			cm.set_regions( region1, region2 );
		} else {
			if ( current_region.find("ligand=")>pos ) {
				region1 = cm.parse_region_string( current_region.substr(0, pos) );
				region2 = cm.parse_region_string( current_region.substr(pos+8) );
				cm.set_regions( region1, region2 );
				cm.set_region2_all_atom( true );
			} else {
				region1 = cm.parse_region_string( current_region.substr(pos+1) );
				region2 = cm.parse_region_string( current_region.substr(7,pos-7) );
				cm.set_regions( region1, region2 );
				cm.set_region2_all_atom( true );
			}
		}
		cm.set_reference_pose( pose.clone() );
		cm.set_distance_cutoff( option[contactMap::distance_cutoff]() );
		cm.set_models_per_file(0);
		cm.set_row_format( option[contactMap::row_format].value() );
		cm.set_distance_matrix( option[contactMap::distance_matrix].value() );

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
		if ( ! input.has_another_pose() ) {
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
		while ( input.has_another_pose() ) {
			input.fill_pose(pose, *rsd_set);
			// loop over ContactMaps and call apply function
			for ( core::Size i=1 ; i <= maps.size(); i++ ) {
				maps[i].apply(pose);
			}
		}
		// loop over ContactMaps and write each one to the corresponding file
		for ( core::Size i=1 ; i <= maps.size(); i++ ) {
			maps[i].write_to_file();
		}

		return 0;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

} // main
