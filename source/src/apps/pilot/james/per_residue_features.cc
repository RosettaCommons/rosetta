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

#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <protocols/jumping/util.hh>

// Utility headers
#include <numeric/xyzVector.hh>
#include <utility/string_util.hh>

// C++ headers
#include <string>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>

#include <platform/types.hh>
#include <core/id/types.hh>

#include <utility/vector1.hh>
#include <basic/options/keys/OptionKeys.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>

#include <utility/excn/Exceptions.hh>

///////////////////////////////////////////////////////////////////////////////

char torsion2big_bin(
	float const phi,
	float const psi,
	float const omega
) {
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

utility::vector1< int > calculate_burial(
	core::pose::Pose const & mypose,
	core::Real const dist_cutoff
) {
	utility::vector1< int > burial;
	burial.resize( mypose.total_residue(), 0 );

	using core::Size;
	for ( Size i = 1; i <= mypose.total_residue(); ++i ) {
		for ( Size j = i + 1; j <= mypose.total_residue(); ++j ) {
			core::conformation::Residue const & resi = mypose.residue(i);
			core::conformation::Residue const & resj = mypose.residue(j);

			core::Real const dist(
				resi.xyz( resi.nbr_atom() ).distance( resj.xyz( resj.nbr_atom() ) )
			);
			if ( dist < dist_cutoff ) {
				burial[ i ]++;
				burial[ j ]++;
			}
		}
	}
	return burial;
}

utility::vector1< char > get_ss( core::pose::Pose & pose ) {
	using core::Size;
	protocols::jumping::assign_ss_dssp( pose );

	utility::vector1< char > ss( pose.total_residue(), 'L' );
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		ss[ii] = pose.secstruct(ii);
	}

	return ss;
}

int main( int argc, char* argv [] ) {
	try {

	// options, random initialization
	devel::init( argc, argv );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;
	using namespace ObjexxFCL::format;
	using namespace core::import_pose::pose_stream;
	using namespace core::chemical;
	using core::Size;
	using core::Real;
	using std::string;
	using utility::vector1;
	using utility::file_basename;

	MetaPoseInputStream input = streams_from_cmd_line();
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();

	std::ostream & output( std::cout );
	while ( input.has_another_pose() ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );
		string const pose_tag( file_basename( core::pose::tag_from_pose( pose ) ) );

		output
			<< A( 6,  "resi"	   )
			<< A( 4,  "res"      )
			<< A( 3,  "ss"       )
			<< A( 9,  "phi"      )
			<< A( 9,  "psi"      )
			<< A( 9,  "omega"    )
			<< A( 9,  "bb_bin"   )
			<< A( 10,  "burial6"  )
			<< A( 10,  "burial8"  )
			<< A( 10,  "burial10" )
			<< A( 14, "pose_tag" )
			<< std::endl;

		vector1< int >  burial6( calculate_burial(pose,6) );
		vector1< int >  burial8( calculate_burial(pose,8) );
		vector1< int >  burial10( calculate_burial(pose,10) );
		vector1< char > pose_ss( get_ss( pose ) );

		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			core::conformation::Residue resi = pose.residue(i);
			core::Real phi   = resi.mainchain_torsion(1);
			core::Real psi   = resi.mainchain_torsion(2);
			core::Real omega = resi.mainchain_torsion(3);

			char bb_bin = torsion2big_bin( phi, psi, omega );
			Size const pose_tag_width( std::max( pose_tag.length()+1, (Size)14 ) );
			output
				<< I(  6, i            )
				<< A(  4, resi.name1() )
				<< A(  3, pose_ss[i]   )
				<< F(  9, 3, phi       )
				<< F(  9, 3, psi	     )
				<< F(  9, 3, omega	   )
				<< A(  9, bb_bin       )
				<< I( 10, burial6[i]   )
				<< I( 10, burial8[i]   )
				<< I( 10, burial10[i]  )
				<< A( pose_tag_width, pose_tag )
				<< std::endl;
		}	// for ( unsigned int i = 1; i <= mypose->total_residue(); ++i )
	} // while input.has_another_pose()

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
