// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// // vi: set ts=2 noet:
// //
// // (c) Copyright Rosetta Commons Member Institutions.
// // (c) This file is part of the Rosetta software suite and is made available under license.
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// // (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    apps/public/analysis/parcs_ccs_calc
/// @brief   Application to calculate average collision cross section (in squared Angstroms) given a structure
/// @author  SM Bargeen Alam Turzo <turzo.1@osu.edu>
///

// Devel header
#include <devel/init.hh>
// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/energy_methods/CCS_IMMSEnergy.hh>
// Output file header
#include <utility/io/ozstream.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
// Datache headers
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <utility/exit.hh>

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::pose::datacache;

// Local Options
basic::options::IntegerOptionKey const n_rots( "ccs_nrots" ); // Options to give the number of rotations for PARCS
basic::options::RealOptionKey const p_rad( "ccs_prad" ); // Options to give the size of Probe radius

static basic::Tracer TR( "apps.public.analysis.parcs_ccs_calc" );

int main(int argc, char * argv[])
{
	try{
		option.add( n_rots, "Number of Rotations" ).def(300); // Def rotations fixed at 25
		option.add( p_rad, "Probe Radius in Angstroms" ).def(1.0); // Def probe radius size 1.0 Angstroms for He buffer gas
		devel::init (argc,argv);
		core::Size nrot( option[ n_rots ].value() ); // Gets value from user input for number of rotations
		core::Real prad( option[ p_rad ].value() ); // Gets value from user input for the size of probe radius in Angstroms
		if ( nrot <=0 ) utility_exit_with_message("Number of rotations cannot be less than or equal to 0");
		if ( prad <=0 ) utility_exit_with_message("Probe cannot be less than or equal to 0");
		core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
		core::Size pose_counter(0);
		std::ostringstream out;
		out <<"File_Name\tCCS_PARCS" << std::endl;
		while ( input.has_another_pose() ) {
			core::pose::PoseOP mypose( utility::pointer::make_shared < core::pose::Pose >() );
			input.fill_pose ( *mypose );
			pose_counter +=1;
			//TR << "Pose Number from list: " << pose_counter << std::endl;
			core::Real parcs_ccs( core::energy_methods::parcs_ccs( *mypose, nrot, prad) );
			std::string decoy_name("Empty_Tag");
			if ( (*mypose).data().has( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
				decoy_name = static_cast< basic::datacache::CacheableString const & >
					( (*mypose).data().get( core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
			}
			TR << pose_counter << "\t" << decoy_name << "\t" << parcs_ccs << std::endl;
			out << decoy_name << "\t" << parcs_ccs << std::endl;
		}
		if ( pose_counter == 0 ) utility_exit_with_message("No structures were given to the PARCS applications. This must be given with either -in:file:s your_pdb.pdb, -in:file:l your_pdb_list, or -in:file:silent your_rosetta_generated_silent_file.");
		std::string outfile = "CCS_default.out";
		if ( option[ out::file::o ].user() ) {
			outfile = option[ out::file::o ]();
		}
		{
			utility::io::ozstream outz( outfile.c_str() );
			outz << out.str();
		}
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}

