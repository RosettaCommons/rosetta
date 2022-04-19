// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

//Load in a membrane protein and calculate the default inner radius for a micelle/bicelle.
//To set default radius, we first measure the max distance between any CA atoms within 3A
//of membrane center plane (in membrane coordinates this is the xy plane), called here "protein_core".
//The default radius is then 3*protein_core/2.
//This application was written as a step in development and is not meant to be used at this point.
//The code is now in core/conformation/membrane/membrane_geometry/Bicelle.cc

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <core/pose/Pose.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/variant_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/scoring/etable/EtableOptions.hh>
#include <core/conformation/Conformation.hh>


#include <numeric/xyzVector.hh>
#include <core/id/AtomID.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/moves/MoverContainer.hh>
//membrane specific
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/visualize/VisualizeEmbeddingMover.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <core/energy_methods/FaMPEnvEnergy.hh>
#include <core/scoring/memb_etable/MembEtable.hh>


//This function calculates the distance between every CA atom within 3A of the membrane center plane.
//In membrane coordinates the membrane center plane is the xy plane with the z axis as the membrane normal.
//Instead of returning the absolute max distance found, it returns the 95th percentile of the CA distances.
core::Real
membrane_center_protein_core( core::pose::PoseOP mypose ) {
	//create conformation object from pose
	const core::conformation::Conformation& myconf = mypose->conformation();

	//Add membrane
	protocols::membrane::AddMembraneMoverOP addmem( new protocols::membrane::AddMembraneMover() );
	addmem->apply( *mypose );

	//std::vector<core::Real> dist_array;
	core::Size num_res = mypose->size();
	//array to store residue positions that are within 3A of membrane center
	utility::vector1<core::Vector> res_xyz;
	//array to store distances between CA atoms
	utility::vector1<core::Real> dist_array;

	//loop through residues
	for ( core::Size i = 1; i <= num_res; i++ ) {
		//set residue and number of atoms in that residue
		core::conformation::Residue const & res_i = myconf.residue(i);
		//grab xyz of CA (this is atom 2)
		core::Vector i_ca_xyz = res_i.xyz(2);

		//if CA within 3 A of membrane center continue
		core::Real dist_from_center_i =  myconf.membrane_info()->atom_z_position( myconf, i, 2 );
		if ( std::abs(dist_from_center_i) <= 3 ) {

			res_xyz.push_back( i_ca_xyz );
		}

		//loop through residues within 3 A of membrane center
		for ( core::Size ai = 1; ai <= res_xyz.size(); ai++ ) {
			//loop through list again
			for ( core::Size aj = 1; aj <= res_xyz.size(); aj++ ) {
				//skip if same resiue
				if ( aj != ai ) {
					core::Vector i_xyz = res_xyz[ai];
					core::Vector j_xyz = res_xyz[aj];
					//calcuate distance between CA atoms in xy dimensions
					core::Real dist = pow( pow(i_xyz.x() - j_xyz.x(), 2) + pow(i_xyz.y() - j_xyz.y(), 2), 0.5);
					//store in array
					dist_array.push_back( dist );
				}
			}
		}
	}
	//find 95 percentile of distances
	std::sort(dist_array.begin(), dist_array.end());
	//multiply k percent by the total number of values,n
	core::Size index = 0.95 * dist_array.size();
	core::Real protein_core = dist_array[index];
	return protein_core;
}


//Based on the protein "diameter" at the membrane center, this returns a value
//for the micelle/bicelle inner radius to be set as.
core::Real
micelle_raidus_from_protein_core( core::Real protein_core ) {
	core::Real micelle_radius = 3*protein_core/2;
	return micelle_radius;
}


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
		//calculate protein core
		core::Real protein_core = membrane_center_protein_core( mypose );
		std::cout << "Protein core: " << protein_core << std::endl;
		core::Real micelle_radius = micelle_raidus_from_protein_core( protein_core);
		std::cout << "Micelle radius: " << micelle_radius << std::endl;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}


	return 0;
}
