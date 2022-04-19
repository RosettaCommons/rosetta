// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

//with this application I want to load in a membrane protein using the MPFramework
//and output the value of the transition function to the bfactor column in the pdb file.

#include <iostream>
#include <cmath>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/variant_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>

#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>

//membrane specific
#include <protocols/membrane/AddMembraneMover.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/MembraneGeometry.hh>


int main( int argc, char ** argv ) {
	try {
		devel::init( argc, argv );

		//read in pdb provided in command line
		utility::vector1< std::string > filenames = basic::options::option[basic::options::OptionKeys::in::file::s ].value();
		if ( filenames.size() > 0 ) {
			std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
		} else {
			std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
			return 1;
		}

		//construct a pose object from pdb file
		core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );

		//this is done before AddMembraneMover so if there is not a MEM residue this should not include it
		//TODO fix this for if there is a MEM residue in pose when read in that you don't try to calculate f(z) for MEM residue
		core::Size num_res = mypose->size();
		//add MEM virtual residue to pose
		//Does this depend on the span file or does it automatically place the center at (0,0,0)?
		//^look into this
		protocols::membrane::AddMembraneMoverOP addmem( new protocols::membrane::AddMembraneMover() );
		addmem->apply( *mypose );

		//create conformation object from pose
		const core::conformation::Conformation& myconf = mypose->conformation();
		core::conformation::membrane::MembraneGeometryCOP mp_geometry( myconf.membrane_info()->membrane_geometry() );

		//loop through residues in pose
		std::cout << "looping through resdiues" << std::endl;
		for ( core::Size aa = 1; aa <= num_res; aa++ ) {
			//set residue and number of atoms in that residue
			core::conformation::Residue const & myres = myconf.residue(aa);
			core::Size num_atoms = myres.natoms();

			//loop through atoms in a residue
			for ( core::Size atom = 1; atom <= num_atoms; atom++ ) {
				//get value of transition function
				core::Real f = mp_geometry->f_transition( myconf, aa, atom );
				//save transition value as bfactor for atom
				mypose->pdb_info()->bfactor( aa, atom, f );

			}
		}
		//save pdb with updated bfactor values
		mypose->dump_pdb( "mp_transition_"+filenames[1] );
		std::cout << "Finished writing to new pdb file with transition function value as b factor. File name mp_transition_" << filenames[1] << std::endl;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
