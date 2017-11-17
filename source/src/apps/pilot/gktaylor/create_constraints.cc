// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author James Thompson

#include <protocols/viewer/viewers.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <core/types.hh>
#include <core/scoring/sasa.hh>
#include <basic/prof.hh> // profiling

#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include "homolog_cst.hh"

#include <utility/excn/Exceptions.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace ObjexxFCL::format;

utility::vector1< int > calculate_burial(
	core::pose::Pose mypose
);

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char* argv [] )
{
    try {
	// options, random initialization
	devel::init( argc, argv );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//int min_sequence_sep = option[ constraints::min_sequence_sep];
	//double max_distance = option[ constraints::max_distance];
	int min_sequence_sep = 3; // use basic::options later
	double max_distance = 14.0;
	utility::vector1< std::string > pdbfiles = option[ james::pdbfile ]();


	utility::vector1< std::string >::const_iterator iter;
	for ( iter = pdbfiles.begin(); iter != pdbfiles.end(); ++iter ) {
		std::string pdbfile = *iter;
		if ( pdbfile == "" ) {
			utility_exit_with_message( "Unable to open file: " + pdbfile + '\n' );
		}
		core::pose::PoseOP mypose ( new core::pose::Pose );
		std::cerr << "READING " << pdbfile << '\n';
		core::import_pose::pose_from_file( *mypose, pdbfile , core::import_pose::PDB_file); // default is standard fullatom residue_set

		std::string outfile  = pdbfile + ".distances";
		std::ofstream output( outfile.c_str() );
		if ( ! output.is_open() ) {
			utility_exit_with_message( "Unable to open file: " + outfile + '\n' );
		}

		// iterate over all atom pairs
		output 	<< A( 13, "[ atompairs ]" )
						<< '\n';

		// calculate burial
		utility::vector1< int > burial = calculate_burial( *mypose );
		for ( unsigned int i = 1; i <= mypose->size(); ++i ) {
			for ( unsigned int j = i + 1; j <= mypose->size(); ++j ) {
				core::conformation::Residue resi = mypose->residue(i);
				core::conformation::Residue resj = mypose->residue(j);

				for ( unsigned int m = 1; m <= resi.natoms(); ++m ) {
					for ( unsigned int n = 1; n <= resj.natoms(); ++n ) {

						// skip hydrogen atoms
						if ( resi.atom_type(m).is_hydrogen() || resj.atom_type(n).is_hydrogen() ) {
							continue;
						}

						// skip atoms that don't share a type
						if ( resi.atom_type(m).name() != resj.atom_type(n).name() ) {
							continue;
						}

						core::Real distance  = mypose->residue(i).xyz(m).distance( resj.xyz(n) );
						// core::Real bfactor_i = mypose->residue(i).atom(m).temperature();
						// core::Real bfactor_j = mypose->residue(j).atom(n).temperature();

						int sequence_sep = j - i;
						// sanity check
						if (distance == 0) {
							std::cout << "0.0 distance detected" << '\n';
							std::cout
										<< A(  8, resi.atom_name(m) )
										<< I( 10, i )
										<< A(  8, resj.atom_name(n) )
										<< I( 10, j )
										<< A( 10, "HARMONIC")
										<< F( 10, 4, distance ) // distance
										<< F( 10, 4, 1.5 )      // std_dev
										<< A( 10, "TAG")
										<< '\n';
							continue;
						}

						if ( sequence_sep > min_sequence_sep && distance < max_distance &&
								( ! resi.atom_name(m).compare(" CA ") && ! resj.atom_name(n).compare(" CA ") ) ) {
								//! resi.atom_name(m).compare(" CB ") && ! resj.atom_name(n).compare(" CB ") )    ) { // use basic::options later
						output
										<< A(  8, resi.atom_name(m) )
										<< I( 10, i )
										<< A(  8, resj.atom_name(n) )
										<< I( 10, j )
										<< A( 10, "HARMONIC")
										<< F( 10, 4, distance ) // distance
										<< F( 10, 4, 1.5 )      // std_dev
										<< A( 10, "TAG")
										<< '\n';
						}
					} 	// for ( unsigned int n = 1; n <= resj.natoms(); ++n )
				} // 	for ( unsigned int m = 1; m <= resi.natoms(); ++m )
			} // 	for ( unsigned int j = i + 1; j <= mypose->size(); ++j )
		}		// for ( unsigned int i = 1; i <= mypose->size(); ++i )
		output.close();
	} // 	for ( iter = pdbfiles.begin(); iter != pdbfiles.end(); ++iter )
    } catch (utility::excn::Exception const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
            return 0;
        } // int main( int argc, char * argv [] )
