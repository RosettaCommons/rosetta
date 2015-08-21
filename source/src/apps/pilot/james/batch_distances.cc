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
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <apps/pilot/james/james_util.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <devel/init.hh>

// Objexx
#include <ObjexxFCL/format.hh>

// Utility headers
#include <basic/Tracer.fwd.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>

// C++ headers
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <string>

#include <utility/excn/Exceptions.hh>

using namespace ObjexxFCL::format;

int main( int argc, char* argv [] ) {
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using core::Real;
		using core::Size;
		using std::string;
		using utility::file_basename;

		// options, random initialization
		devel::init( argc, argv );

		// maximum Euclidean distance between a pair of atoms to be in contact
		Real distance_cutoff = std::numeric_limits<Real>::max();
		if ( option[ james::dist_thresholds ].user() ) {
			distance_cutoff = option[ james::dist_thresholds ]().front();
		}

		// minimum sequence separation between a pair of atoms to be in contact
		Size min_seqsep = option[ james::min_seqsep ]();
		bool const show_all_atoms( option[ james::debug ]() );
		utility::vector1< std::string > atoms_wanted( option[ james::atom_names ]() );

		std::ostream& output( std::cout );

		core::import_pose::pose_stream::MetaPoseInputStream input
			= core::import_pose::pose_stream::streams_from_cmd_line();
		ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );

		output  << A( 10, "resi_idx" )
			<< A( 10, "resj_idx" )
			<< A(  6, "resi"     )
			<< A(  6, "resj"     )
			<< A(  8, "atomi"    )
			<< A(  8, "atomj"    )
			<< A( 10, "burial_i" )
			<< A( 10, "burial_j" )
			<< A( 10, "dist"     )
			<< A( 10, "seqsep"   )
			<< A( 14, "pose_tag" )
			<< std::endl;

		while ( input.has_another_pose() ) {
			core::pose::Pose pose;
			input.fill_pose( pose, *rsd_set );
			string const pose_tag( file_basename( core::pose::tag_from_pose( pose ) ) );
			utility::vector1< int > burial = calculate_burial( pose );
			for ( Size i = 1; i <= pose.total_residue(); ++i ) {
				for ( Size j = i+1; j <= pose.total_residue(); ++j ) {
					const core::conformation::Residue& resi = pose.residue(i);
					const core::conformation::Residue& resj = pose.residue(j);

					for ( Size m = 1; m <= resi.natoms(); ++m ) {
						for ( Size n = 1; n <= resj.natoms(); ++n ) {
							if ( atoms_wanted.size() > 0 ) {
								std::string const & atom_m( resi.atom_type(m).name() );
								std::string const & atom_n( resj.atom_type(n).name() );
								utility::vector1< std::string >::const_iterator find_m(
									std::find(atoms_wanted.begin(),atoms_wanted.end(),atom_m)
								);
								utility::vector1< std::string >::const_iterator find_n(
									std::find(atoms_wanted.begin(),atoms_wanted.end(),atom_n)
								);
								if ( find_m == atoms_wanted.end() || find_n == atoms_wanted.end() ) {
									continue;
								}
							}

							// skip hydrogen atoms
							if ( resi.atom_type(m).is_hydrogen() ||
									resj.atom_type(n).is_hydrogen() ) {
								continue;
							}

							// skip atoms that don't share a type
							if ( resi.atom_type(m).name() != resj.atom_type(n).name() &&
									!show_all_atoms ) {
								continue;
							}

							// determine whether Euclidean distance is below threshold
							core::Real const distance(
								pose.residue(i).xyz(m).distance( resj.xyz(n) ));
							if ( distance > distance_cutoff ) continue;

							// determine whether sequence separation exceeds threshold
							Size pos_i = pose.pdb_info()->number(i);
							Size pos_j = pose.pdb_info()->number(j);
							Size seqsep = std::abs( int(pos_j) - int(pos_i));
							if ( seqsep < min_seqsep ) continue;

							// sometimes pose_tag is too big.
							Size const pose_tag_width( std::max( pose_tag.length(), (Size)14 ) );

							output << I( 10, pos_i)
								<< I( 10, pos_j)
								<< A(  6, resi.name1() )
								<< A(  6, resj.name1() )
								<< A(  8, resi.atom_name(m) )
								<< A(  8, resj.atom_name(n) )
								<< I( 10, burial[i] )
								<< I( 10, burial[j] )
								<< F( 10, 4, distance )
								<< I( 10, seqsep)
								<< A( pose_tag_width, pose_tag )
								<< std::endl;
						}
					}
				}
			}
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
