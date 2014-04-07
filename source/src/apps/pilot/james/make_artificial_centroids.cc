// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file make_artificial_centroids.cc
/// @brief takes an input set of PDBs, creates artificial centroid by averaging
/// coordinates
/// @author James Thompson

// libRosetta headers

#include <devel/init.hh>
#include <core/types.hh>

#include <core/chemical/util.hh>

#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/scoring/rms_util.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <basic/prof.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>
#include <boost/dynamic_bitset.hpp>
#include <core/id/AtomID.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

core::Size count_backbone_atoms(
	core::pose::Pose const & pose
) {
	core::Size count(0);
	for ( core::Size ii = 1, end = pose.total_residue(); ii <= end; ++ii ) {
		core::conformation::Residue const & res(pose.residue(ii));
		for ( int jj = 1; jj <= res.natoms(); ++jj ) {
			if ( res.atom_is_backbone(jj) ) ++count;
		}
	}
	return count;
}

core::Real calc_mean(
	utility::vector1< core::Real > const & vals
) {
	using core::Size;
	using core::Real;
	Real mean(0.0);
	for ( Size ii = 1; ii <= vals.size(); ++ii ) {
		mean += vals[ii] / vals.size();
	}
	return mean;
}

core::Real calc_sdev(
	utility::vector1< core::Real > const & vals
) {
	using core::Size;
	using core::Real;
	Real const mean( calc_mean(vals) );
	Real sum(0.0);
	for ( Size ii = 1; ii <= vals.size(); ++ii ) {
		sum += std::pow( mean - vals[ii], 2 );
	}
	return std::sqrt( sum );
}

int
main( int argc, char * argv [] ) {
	try {

	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::import_pose::pose_stream;
	using namespace core::io::silent;

	using core::Real;
	using core::Size;
	using std::string;
	using utility::vector1;
	using core::pose::Pose;
	using core::PointPosition;

	devel::init( argc, argv );

	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	MetaPoseInputStream input = streams_from_cmd_line();

	vector1< vector1< vector1< PointPosition > > > all_coords; // all_coords[struct_i][resi][atom_i] = coord
	core::pose::Pose pose; // will use the last pose later
	while( input.has_another_pose() ) {
		input.fill_pose( pose, *rsd_set );

		//Size const n_backbone( count_backbone_atoms(pose) );
		//if ( all_coords.size() != 0 ) {
		//	runtime_assert( n_backbone == all_coords.front().size() );
		//}

		vector1< vector1< PointPosition > > pose_coords( pose.total_residue() );
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			core::conformation::Residue const & res(pose.residue(ii));
			vector1< PointPosition > res_coords;
			for ( int jj = 1; jj <= res.natoms(); ++jj ) {
				//if ( res.atom_is_backbone(jj) ) res_coords.push_back(res.xyz(jj));
				res_coords.push_back(res.xyz(jj));
			}
			pose_coords.push_back( res_coords );
		}

		all_coords.push_back( pose_coords );
	}

	Size const nstruct( all_coords.size() );
	//vector1< vector1< PointPosition > > avg_coords(
	//	all_coords.front().size(), vector1< PointPosition >(
	//		all_coords.front().front().size(), PointPosition(0,0,0)
	//	)
	//); // avg_coords[res_jj][atom_kk] = coord
	vector1< vector1< PointPosition > > avg_coords;
	for ( Size jj = 1; jj <= all_coords.front().size(); ++jj ) { // res jj
		avg_coords.push_back( vector1< PointPosition >( all_coords.front().size() ) );
	}

	for ( Size ii = 1; ii <= all_coords.size(); ++ii ) { // struct ii
		for ( Size jj = 1; jj <= all_coords[ii].size(); ++jj ) { // res jj
			for ( Size kk = 1; kk <= all_coords[ii][jj].size(); ++kk ) { // atom kk
				avg_coords[jj][kk] += all_coords[ii][jj][kk].x() / nstruct;
				avg_coords[jj][kk] += all_coords[ii][jj][kk].y() / nstruct;
				avg_coords[jj][kk] += all_coords[ii][jj][kk].z() / nstruct;
			}
		}
	}

	pose.dump_pdb( "debug.pdb" );

	for ( Size jj = 1; jj <= pose.total_residue(); ++jj ) {
		core::conformation::Residue const & res(pose.residue(jj));
		for ( Size kk = 1; kk <= res.natoms(); ++kk ) {
			using core::id::AtomID;
			pose.set_xyz( AtomID(kk,jj), avg_coords[jj][kk] );
			//if ( res.atom_is_backbone(jj) ) pose.set_xyz( id, avg_coords[ii+jj] );
		}
	}
	pose.dump_pdb( "artificial_centroid.pdb" );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
