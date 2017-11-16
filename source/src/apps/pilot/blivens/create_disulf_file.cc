// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file create_disulf_file.cc
/// @brief Given a pdb file, detect disulfides and output to a disulfide file
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date Created 3/22/2009
/// @details
/// Creates a file pairs of numbers giving rosetta residue indices
/// for each disulfide bond. This file can then be used with -fix_disulf
///
/// Only detects disulfides which have the right rotamer for each Cys.
/// @note requires a full atom input pdb
/// @section cli Command Line
/// @code create_disulf_file -database db -s in.pdb -o out.disulf @endcode

#include <apps/pilot/blivens/disulfides.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
//#include <core/chemical/VariantType.hh>

#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

#include <fstream>
#include <utility>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace std;
using utility::vector1;

static basic::Tracer TR( "apps.pilot.blivens.create_disulf_file" );

int
usage(string msg)
{
	TR << "usage: create_disulf_file -database db -s in.pdb -o out.disulf" << endl
		<< msg << endl;
	exit(1);
}

int main( int argc, char * argv [] )
{
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		//init options system
		devel::init(argc, argv);

		//get i/o files
		string pdbfile, outfile;
		if ( option[ in::file::s ].user() ) {
			pdbfile = start_file();
		} else return usage("No in file given: Use -s to designate pdb files to search for disulfides");

		if ( option[ out::file::o ].user() ) {
			outfile = option[ out::file::o ]();
		} else return usage("No out file given: Use -o to designate a file or directory to output to");

		//read in pose
		pose::Pose pose;
		core::import_pose::pose_from_file( pose, pdbfile , core::import_pose::PDB_file);

		ofstream out(outfile.c_str());
		if ( !out.good() ) {
			return usage( "Error opening "+outfile+" for writing." );
		}

		//find disulfides
		vector1< std::pair<Size,Size> > disulf_bonds(0); //zero needed to disambiguate constructor?
		find_disulfides(pose,disulf_bonds);

		//output
		for ( vector1<pair<Size,Size> >::const_iterator disulf_it = disulf_bonds.begin();
				disulf_it != disulf_bonds.end(); ++disulf_it ) {
			out << disulf_it->first << "\t" << disulf_it->second << endl;
		}

		out.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // end main

