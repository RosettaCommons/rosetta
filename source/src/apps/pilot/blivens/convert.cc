// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file convert.cc
/// @brief Converts PDBs between residue type sets. Use -in:file:residue_type_set to specify type sets
/// @author Spencer Bliven <blivens@u.washington.edu>
/// @date Created October 2008
/// @details
/// @section cli Command Line
/// @code convert2centroid -s fa.pdb -o centroid.pdb -database db @endcode


#include <core/chemical/ChemicalManager.hh>

//Utility
#include <devel/init.hh>
#include <utility/vector1.hh>

//Options
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>


using namespace core;
using namespace std;
using utility::vector1;
using basic::options::option;
using namespace basic::options::OptionKeys;

#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>


static basic::Tracer TR( "pilot_apps.blivens.convert" );

int
usage(char* msg)
{
	TR << "usage: convert -s input.pdb -o output.pdb -database db" << endl
		<< "Example: convert -s in.pdb -in:file:fullatom -o output.pdb -out:file:residue_type_set centroid" <<endl
		<< msg << endl;
	exit(1);
}

int main( int argc, char * argv [] )
{
	try {
		//init options system
		option.add_relevant( in::file::s );
		option.add_relevant( out::file::o );
		option.add_relevant( in::file::residue_type_set );
		option.add_relevant( out::file::residue_type_set );
		option.add_relevant( in::file::fullatom );
		option.add_relevant( out::file::fullatom );
		option.add_relevant( in::file::centroid_input );
		devel::init(argc, argv);


		if ( ! option[ in::file::s ].user() ) {
			return usage("No in file given: Use -s or -l to designate pdb files to search for disulfides");
		}
		if ( ! option[ out::file::o ].user() ) {
			return usage("No out file given: Use -o to designate an output file");
		}

		string pdb = basic::options::start_file();
		string out = option[ out::file::o ]();

		string in_rsd_set;
		if ( option[ in::file::centroid_input ].user() && option[ in::file::centroid_input]() ) {
			if ( option[ in::file::fullatom ].user() && option[ in::file::fullatom ]() ) {
				utility_exit_with_message("Conflicting values for -in:file:centroid_input and -in:file:fullatom");
			}
			in_rsd_set = chemical::CENTROID;
		} else {
			in_rsd_set = option[ in::file::residue_type_set ](); //default fa_standard
		}

		string out_rsd_set;
		if ( option[ out::file::fullatom ].user() && option[ out::file::fullatom ]() ) {
			out_rsd_set = chemical::FA_STANDARD;
		} else {
			out_rsd_set = option[ out::file::residue_type_set ](); //default fa_standard
		}

		TR.Info << "Converting from " << in_rsd_set <<" to " << out_rsd_set << endl;
		chemical::ResidueTypeSetCAP rsd_set =
			chemical::ChemicalManager::get_instance()->residue_type_set(in_rsd_set);
		pose::Pose pose;

		core::import_pose::pose_from_file( pose, *rsd_set, pdb, false, core::import_pose::PDB_file);
		core::util::switch_to_residue_type_set( pose, out_rsd_set);
		pose.dump_pdb(out);

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // end main

