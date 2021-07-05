// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk
#include <string>
#include <sstream>

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

//Protocol Headers
#include <protocols/pockets/PocketGrid.hh>
//#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace basic::options::OptionKeys;
using namespace conformation;
using namespace core::pose::datacache;
using namespace core::id;


static basic::Tracer TR( "apps.pilot.david_pocket_compare.main" );

/// General testing code
int main( int argc, char * argv [] ) {

	try {

		//NEW_OPT( pocket1_fname, "pocket", "fname" );

		//initializes Rosetta functions
		devel::init(argc, argv);

		std::string const resid (option[ OptionKeys::pocket_grid::central_relax_pdb_num ]  );

		for ( core::Size f=1; f <= basic::options::start_files().size(); f++ ) {
			std::string const fname1 = basic::options::start_files().at(f);
			pose::Pose input_pose;

			//read in pdb file from command line
			core::import_pose::pose_from_file( input_pose, fname1 , core::import_pose::PDB_file);
			std::vector< conformation::ResidueCOP > residues = protocols::pockets::PocketGrid::getRelaxResidues(input_pose, resid);


			protocols::pockets::PocketGrid pg( residues );
			pg.autoexpanding_pocket_eval( residues, input_pose );


			std::stringstream out_fname;

			out_fname << fname1 << "."<< option[ OptionKeys::out::output_tag ]()<<".min.pdb";
			pg.dumpTargetPocketsToPDB ( out_fname.str(), true);
		}

		TR << "Done!" << std::endl;
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}

