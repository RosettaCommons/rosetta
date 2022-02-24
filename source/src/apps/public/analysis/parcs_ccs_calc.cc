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
// Core headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/energy_methods/CCS_IMMSEnergy.hh>
// Output file header
#include <utility/io/ozstream.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>


using namespace basic::options;
using namespace basic::options::OptionKeys;
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
		core::Real prad( option[ p_rad ].value() ); // Gets value from user input for the size of probe radius in A
		utility::vector1< std::string > pdbnames( basic::options::start_files() ); // Option to read in pdb with in:file:s or in:file:l
		std::ostringstream out;
		out <<"File_Name\tCCS_PARCS" << std::endl;
		TR << "PDB_File\tPARCS_CCS" << std::endl;
		for ( core::Size pdb = 1; pdb <= pdbnames.size(); ++pdb ) {
			if ( pdbnames[ pdb ].size() <= 0 ) {
				TR << " PDB file not found, did you use -in::file::s or -in::file::l option to enter files?" << std::endl;
				return 1;
			}
			core::pose::Pose mypose = *core::import_pose::pose_from_file(pdbnames[ pdb ]);
			core::Real parcs_ccs = core::energy_methods::parcs_ccs( mypose, nrot, prad );
			TR << pdbnames[ pdb ]<<"\t"<< parcs_ccs << std::endl;
			out << pdbnames[ pdb ] << "\t " << parcs_ccs << std::endl;
		}
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

