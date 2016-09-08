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


#include <core/types.hh>


#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/options/option.hh>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <ObjexxFCL/format.hh>

#include <utility/excn/Exceptions.hh>


using namespace ObjexxFCL::format;

///////////////////////////////////////////////////////////////////////////////

// burial is represented as number of neighbor atoms within 15 angstroms
utility::vector1< int > calculate_burial(
	core::pose::Pose mypose
) {

	utility::vector1< int > burial;
	burial.resize( mypose.size() );

	for ( unsigned int i = 1; i <= mypose.size(); ++i ) {
		for ( unsigned int j = i + 1; j <= mypose.size(); ++j ) {
			core::conformation::Residue resi = mypose.residue(i);
			core::conformation::Residue resj = mypose.residue(j);

			if ( resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) ) < 15 ) {
				burial[ i ]++;
				burial[ j ]++;
			}
		}
	}

	return burial;
}

int
main( int argc, char* argv [] )
{
    try {
	// options, random initialization
	devel::init( argc, argv );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< std::string > pdbfiles = option[ in::file::s ]();

	utility::vector1< std::string >::const_iterator iter;
	for ( iter = pdbfiles.begin(); iter != pdbfiles.end(); ++iter ) {
		std::string pdbfile = *iter;
		if ( pdbfile == "" ) {
			utility_exit_with_message( "Unable to open file: " + pdbfile + '\n' );
		}
		core::pose::PoseOP mypose ( new core::pose::Pose );
		std::cerr << "READING " << pdbfile << std::endl;
		core::import_pose::pose_from_file( *mypose, pdbfile , core::import_pose::PDB_file); // default is standard fullatom residue_set

		std::string outfile  = pdbfile + ".burial";
		std::ofstream output( outfile.c_str() );
		if ( ! output.is_open() ) {
			utility_exit_with_message( "Unable to open file: " + outfile + '\n' );
		}

		// iterate over all atom pairs
		output 	<< A( 10, "resi_idx" )
						<< A( 10, "resj_idx" )
						<< A(  6, "resi"     )
						<< A(  6, "resj" 		 )
						<< A(  8, "atomi"    )
						<< A(  8, "atomj"    )
						<< A( 10, "burial_i" )
						<< A( 10, "burial_j" )
						// << A( 10, "temp_i"   )
						// << A( 10, "temp_j"   )
						<< A( 10, "dist"     )
						<< std::endl;

		// calculate burial
		utility::vector1< int > burial = calculate_burial( *mypose );
		for ( unsigned int i = 1; i <= mypose->size(); ++i ) {
//			for ( unsigned int j = i + 1; j <= mypose->size(); ++j ) {
				core::conformation::Residue resi = mypose->residue(i);
//				core::conformation::Residue resj = mypose->residue(j);

				for ( unsigned int m = 1; m <= resi.natoms(); ++m ) {
//					for ( unsigned int n = 1; n <= resj.natoms(); ++n ) {

						// skip hydrogen atoms
//						if ( resi.atom_type(m).is_hydrogen() || resj.atom_type(n).is_hydrogen() ) {
						if ( resi.atom_type(m).is_hydrogen() ) {
							continue;
						}

						// skip atoms that don't share a type
//						if ( resi.atom_type(m).name() != resj.atom_type(n).name() ) {
//							continue;
//						}

						//if ( resi.atom_type(m).name() != "CA" ) {
						//	continue;
						//}

//						core::Real distance  = mypose->residue(i).xyz(m).distance( resj.xyz(n) );
						// core::Real bfactor_i = mypose->residue(i).atom(m).temperature();
						// core::Real bfactor_j = mypose->residue(j).atom(n).temperature();

						output 	<< I( 10, i )
//										<< I( 10, j )
										<< A(  6, resi.name1() )
//										<< A(  6, resj.name1() )
										<< A(  8, resi.atom_name(m) )
//										<< A(  8, resj.atom_name(n) )
										<< I( 10, burial[i] ) // burial_i
//										<< I( 10, burial[j] ) // burial_j
										// << F( 10, 4, bfactor_i )
										// << F( 10, 4, bfactor_j )
//										<< F( 10, 4, distance ) // distance
										<< std::endl;
//					} 	// for ( unsigned int n = 1; n <= resj.natoms(); ++n )
				} // 	for ( unsigned int m = 1; m <= resi.natoms(); ++m )
//			} // 	for ( unsigned int j = i + 1; j <= mypose->size(); ++j )
		}		// for ( unsigned int i = 1; i <= mypose->size(); ++i )
		output.close();
	} // 	for ( iter = pdbfiles.begin(); iter != pdbfiles.end(); ++iter )
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
            return 0;
        } // int main( int argc, char * argv [] )

