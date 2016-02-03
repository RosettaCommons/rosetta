// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson

#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <basic/options/option.hh>
#include <basic/database/open.hh>
#include <devel/init.hh>

#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>

// Utility headers
#include <numeric/xyzVector.hh>

// C++ headers
#include <string>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <apps/pilot/james/james_util.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>

using namespace core;
using ObjexxFCL::string_of;

///////////////////////////////////////////////////////////////////////////////

char torsion2big_bin(float const phi,
	float const psi,
	float const omega) {
	if ( std::abs( omega ) < 90 ) {
		return 'O'; // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100 < psi && psi <= 100 ) {
			return 'G'; // alpha-L
		} else {
			return 'E'; // E
		}
	} else {
		if ( -125 < psi && psi <= 50 ) {
			return 'A'; // helical
		} else {
			return 'B'; // extended
		}
	}
	return 'X';
}


int main( int argc, char* argv [] ) {
	try {

		// options, random initialization
		devel::init( argc, argv );

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring::constraints;
		using namespace core::pack::dunbrack;
		using namespace ObjexxFCL::format;
		using core::Size;
		using core::Real;
		using std::string;
		using utility::vector1;

		vector1< string > pdbfiles = option[ in::file::s ]();
		vector1< string >::const_iterator iter;
		for ( iter = pdbfiles.begin(); iter != pdbfiles.end(); ++iter ) {
			string pdbfile = *iter;
			if ( pdbfile == "" ) {
				utility_exit_with_message( "Unable to open file: " + pdbfile + '\n' );
			}
			core::pose::PoseOP mypose( new core::pose::Pose );
			core::import_pose::pose_from_file( *mypose, pdbfile , core::import_pose::PDB_file); // default is standard fullatom residue_set

			std::ostream & output( std::cout );

			output
				<< A( 6, "idx"    )
				<< A( 4, "res"     )
				<< A( 3, "ss"      )
				<< A( 9, "phi"     )
				<< A( 9, "psi"     )
				<< A( 9, "omega"   )
				<< A( 9, "CA_x"    )
				<< A( 9, "CA_y"    )
				<< A( 9, "CA_z"    )
				//<< A( 9, "chi1"    )
				//<< A( 9, "chi2"    )
				//<< A( 9, "chi3"    )
				//<< A( 9, "chi4"    )
				<< A( 9, "bb_bin"  )
				//<< A( 9, "chi_bin" )
				<< A( 9, "burial"  )
				<< std::endl;

			vector1< char > pose_ss( get_ss( *mypose ) );
			vector1< int >  burial( calculate_burial( *mypose ) );

			for ( Size i = 1; i <= mypose->total_residue(); ++i ) {
				core::conformation::Residue resi = mypose->residue(i);
				core::Real phi   = resi.mainchain_torsion( 1 );
				core::Real psi   = resi.mainchain_torsion( 2 );
				core::Real omega = resi.mainchain_torsion( 3 );

				// CA cartesian coordinates
				numeric::xyzVector<core::Real> coords = resi.xyz("CA");

				vector1< core::Real > chis( 4, 0.0 );
				for ( Size jj = 1; jj <= Size(std::min( 4, (int) resi.nchi() )); ++jj ) {
					chis[jj] = resi.chi(jj);
				}

				RotVector rot_vector;
				rotamer_from_chi( resi, rot_vector );
				char bb_bin = torsion2big_bin( phi, psi, omega );
				string chi_bin("");
				for ( Size jj = 1; jj <= rot_vector.size(); ++jj ) {
					chi_bin += string_of( rot_vector[jj] );
				}

				output
					<< I(  6, i            )
					<< A(  4, resi.name1() )
					<< A(  3, pose_ss[i]   )
					<< F(  9, 3, phi       )
					<< F(  9, 3, psi      )
					<< F(  9, 3, omega    )
					<< F(  9, 3, coords.x())
					<< F(  9, 3, coords.y())
					<< F(  9, 3, coords.z())
					//<< F(  9, 3, chis[1]   )
					//<< F(  9, 3, chis[2]   )
					//<< F(  9, 3, chis[3]   )
					//<< F(  9, 3, chis[4]   )
					<< A(  9, bb_bin       )
					//<< A(  9, chi_bin      )
					<< I( 10, burial[i]    )
					<< std::endl;
			} // for ( unsigned int i = 1; i <= mypose->total_residue(); ++i )
		} // for ( iter = pdbfiles.begin(); iter != pdbfiles.end(); ++iter )

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
