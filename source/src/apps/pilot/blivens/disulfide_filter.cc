// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file disulfide_filter.cc
/// @brief Outputs a list of residue pairs which are likely candidates for disulfide bonds.
/// @author Spencer Bliven
/// @date Created November 2008
/// @details
/// @section cli Command Line
/// @code disulfide_filter -s input.pdb -o output.txt -database db @endcode


//Utility
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <fstream>

//Options
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>

//Poses
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <apps/pilot/blivens/disulfides.hh>


using namespace core;
using namespace std;
using utility::vector1;
using namespace basic::options;
using namespace basic::options::OptionKeys;

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <basic/Tracer.hh>
static basic::Tracer TR( "pilot_apps.blivens.disulfide_filter" );

int
usage(char* msg)
{
	TR << "usage: disulfide_filter -s in.pdb -o out.txt -database db" << endl
		<< msg << endl;
	exit(1);
}

int main( int argc, char * argv [] )
{
	try {

		//init options system
		devel::init(argc, argv);

		vector1<string> pdbs;
		if ( option[ in::file::s ].user() || option[in::file::l].user() ) {
			pdbs = basic::options::start_files();
		} else return usage("No in file given: Use -s or -l to designate pdb files to search for disulfides");


		ofstream out;
		if ( option[ out::file::o ].user() ) {
			string outfile = option[ out::file::o ]();
			out.open(outfile.c_str());
		} else return usage("No out file given: Use -o to designate an output file");

		pose::Pose pose;
		core::import_pose::pose_from_file( pose, pdbs[1] , core::import_pose::PDB_file);

		//Assign secodary structure
		core::scoring::dssp::Dssp dssp(pose);
		dssp.insert_ss_into_pose(pose);

		TR << pose.secstruct() << endl;


		//for each residue pair
		for ( Size i(1); i<= pose.size()-1; ++i ) {
			for ( Size j(i+1); j<= pose.size(); ++j ) {}
		//end residue pair
		}

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // end main
