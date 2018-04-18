// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/kkappel/score_rnp.cc
/// @brief Setup coordinate constraints based on the coordinates of an input structure
/// @author Kalli Kappel kappel@stanford.edu


//Unit Headers
//Package Headers
//Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
//Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
//Numeric Headers
#include <numeric/random/random.hh>
//C++ Headers
#include <iostream>
#include <fstream>

#include <devel/init.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// New options for this application
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "apps.pilot.kkappel.setup_coor_csts" );

/////////////////////////////////////////////////////////////////////////////////

void get_csts_from_pose( core::pose::Pose const & pose, std::string name ) {

	// File to write constraints to
	std::ofstream cst_outfile;
	std::string fname = name + "_coord_cst.cst";
	cst_outfile.open( fname.c_str() );

	utility::vector1< std::string > atom_names;
	// List of atoms that we're going to constrain
	atom_names.push_back( " P  ");
	atom_names.push_back( " C1'");

	// Loop through all the residues in the pose
	for ( core::Size i = 1; i<= pose.total_residue(); ++i ) {
		// Loop through all the atoms to constrain in this residue
		for ( core::Size atom_i = 1; atom_i <= atom_names.size(); ++atom_i ) {
			//// Coordinate constraints, not working well
			//core::Vector const pos( pose.residue( i ).xyz( atom_names[atom] ) );
			//cst_outfile << "CoordinateConstraint " << atom_names[atom] << " ";
			//cst_outfile << i << "  " << " P  1  " ;
			//cst_outfile << pos[1] << "  " <<  pos[2] << "  " << pos[3];
			//cst_outfile << "  HARMONIC 0.0 5.0 " << std::endl;

			// Try atom pair constraints instead
			core::Vector const pos_1( pose.residue( i ).xyz( atom_names[atom_i] ) );
			// Loop through
			for ( core::Size j = 1; j < i; ++j ) {
				// Again loop through all the atoms to constrain
				for ( core::Size atom_j = 1; atom_j <= atom_names.size(); ++atom_j ) {
					core::Vector const pos_2( pose.residue( j ).xyz( atom_names[atom_j] ) );
					core::Real const dist = (pos_1 - pos_2).length();
					cst_outfile << "AtomPair " << atom_names[atom_j] << " " << i << " " << atom_names[atom_i] << " " << j << " "
						<< "HARMONIC " << dist << " " << 2.0 << std::endl;
				}
			}

		}
	}
	cst_outfile.close();
}

///////////////////////////////////////////////////////////////////////////////
void setup_coord_csts() {

	using namespace core::scoring;
	using namespace basic::options;

	core::pose::Pose pose;

	//Get the input structure(s)
	utility::vector1< std::string > input_structs = option[ OptionKeys::in::file::s ]();
	for ( core::Size i = 1; i <= input_structs.size(); ++i ) {
		core::import_pose::pose_from_file( pose, input_structs[ i ] );
		std::string out_tag( utility::string_split( utility::string_split( input_structs[i], '/').back(), '.' ).front() );
		get_csts_from_pose( pose, out_tag );
		// if the number of input structures is greater than 1, then we can index the
		// output files, otherwise don't bother
	}

}

/////////////////////////////////////////////////////////////////////////////////
int main( int argc, char ** argv ) {

	try {
		using namespace basic::options;

		devel::init( argc, argv );
		setup_coord_csts();
	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
