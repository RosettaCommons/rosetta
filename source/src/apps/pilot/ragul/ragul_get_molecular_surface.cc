// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ragul Gowthaman

//GPU enabling is not default
//To test how many threads are fastest for your computer,
//use -gpu:threads 1024 (or other number) on the command line

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

// Protocol Headers
#include <devel/init.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <basic/options/option_macros.hh>

// Utility Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <utility/options/StringOption.hh>
#include <core/scoring/sc/MolecularSurfaceCalculator.hh>
#include <core/pack/task/TaskFactory.hh>
#include <vector>


using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;
using namespace core::scoring::sc;

OPT_KEY( String, protein )

int main( int argc, char * argv [] ) {

	try {

		NEW_OPT( protein, "protein file name", "protein.pdb" );

		devel::init(argc, argv);

		std::string const input_protein = option[ protein ];
		pose::Pose protein_pose;
		core::import_pose::pose_from_file( protein_pose, input_protein , core::import_pose::PDB_file);

		//create 'tag' for eggshell output filename
		int pfounddir = input_protein.find_last_of("/\\");
		int pfounddot = input_protein.find_last_of(".");
		std::string protein_name = input_protein.substr((pfounddir+1),(pfounddot-(pfounddir+1)));
		std::string pro_output_pdbname = protein_name+"_surface.pdb";

		utility::io::ozstream outPDB_stream;
		outPDB_stream.open(pro_output_pdbname, std::ios::out);

		core::scoring::sc::MolecularSurfaceCalculator calculator;
		// calculator.settings.density = surface_density_;
		calculator.Init();
		calculator.Calc(protein_pose);

		std::vector<DOT> surface_dots = calculator.GetDots(0);

		std::cout << "Generated surface dots: " << surface_dots.size() << std::endl;

		//core::Size zeronormal_dots = 0;
		//core::Size skipped_dots = 0;
		//core::Size selected_dots = 0;
		core::pack::task::PackerTaskOP task;
		task = core::pack::task::TaskFactory::create_packer_task( protein_pose );
		for ( core::Size i = 0; i < surface_dots.size(); i++ ) {
			if ( task->pack_residue(surface_dots[i].atom->nresidue) ) {
				/*
				std::cout<<i << " " <<
				surface_dots[i].atom->nresidue << " " <<
				surface_dots[i].atom->natom << " " <<
				surface_dots[i].atom->residue << " " <<
				surface_dots[i].atom->atom << " " <<
				surface_dots[i].area << " " <<
				surface_dots[i].coor.x() << " " <<
				surface_dots[i].coor.y() << " " <<
				surface_dots[i].coor.z() << " " <<
				surface_dots[i].outnml.x() << " " <<
				surface_dots[i].outnml.y() << " " <<
				surface_dots[i].outnml.z() << " " <<
				"\n";
				*/
				outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   COM X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<surface_dots[i].coor.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<surface_dots[i].coor.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<surface_dots[i].coor.z()<<std::endl;
			}
		}

		outPDB_stream.close();
		outPDB_stream.clear();

	}//try
catch ( utility::excn::EXCN_Base const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}

	return 0;

}
