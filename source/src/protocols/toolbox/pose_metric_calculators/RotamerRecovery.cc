// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin RotamerRecovery
///
/// @brief
/// Compare the rotamer recovery between a native protein and a list of other proteins
///
/// @detailed
/// This is an implementation taken from James Thompson. I am not even sure he knows I stole it
/// from him. The main function that is called is the get_rotamer_recovery() function. You can
/// pass this function a native pdb and a list of altered pdbs, or just 1 native and 1
/// alterd pdb. The rotamer recovery will be output to the screen. Output looks like:
/// # total = 1
///    resi_idx  nat_bb_bin      pct_bb    nat_rot1    pct_rot1    nat_rot2    pct_rot2    nat_rot3    pct_rot3    nat_rot4    pct_rot4
///           1           E      1.0000           1      1.0000           2      1.0000           1      1.0000         999      0.0000
///           2           B      1.0000           2      1.0000           1      1.0000         999      0.0000         999      0.0000
/// Where the # total is how many proteins compared.
/// resi_idx = residue index
/// nat_bb_bin = dssp naming for bb
/// pct_bb = how many match the bb bins?
/// nat_rot1 = chi 1
/// pct_rot1 = how many are correct
/// If 999 appears, that means that the amino acid does not have that chi angle
///
///
///
/// @authors
/// @author James Thompson (original author)
/// @author Steven Combs (moved it to protocols for general use)
///
/// @last_modified October 20 2010
/////////////////////////////////////////////////////////////////////////

#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
//#include <devel/init.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>

// AUTO-REMOVED #include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/util.hh>

// XRW-REMOVE #include <protocols/jumping/util.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/toolbox/pose_metric_calculators/RotamerRecovery.hh>

// option key includes




namespace protocols{
namespace toolbox{
namespace pose_metric_calculators{



char RotamerRecovery::torsion2big_bin(
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

utility::vector1< char > RotamerRecovery::get_ss( core::pose::Pose & pose ) {
	core::scoring::dssp::Dssp dssp(pose);
	dssp.insert_ss_into_pose(pose);

	utility::vector1< char > ss( pose.total_residue(), 'L' );
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		ss[ii] = pose.secstruct(ii);
	}

	return ss;
}
void RotamerRecovery::print_rot_vec(
	core::pack::dunbrack::RotVector rot_vec,
	std::ostream & out
) {
	using core::Size;
	for ( Size ii = 1; ii <= rot_vec.size(); ++ii ) {
		out << rot_vec[ii];
		if ( ii != rot_vec.size() ) {
			out << ",";
		}
	}
}

utility::vector1< char > RotamerRecovery::bb_bins_from_pose(
	core::pose::Pose const & pose
) {
	utility::vector1< char > bb_bins;
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		bb_bins.push_back(
			torsion2big_bin(
				pose.residue(ii).mainchain_torsion(1),
				pose.residue(ii).mainchain_torsion(2),
				pose.residue(ii).mainchain_torsion(3)
			)
		);
	}
	return bb_bins;
}

utility::vector1< core::pack::dunbrack::RotVector >
RotamerRecovery::rots_from_pose(
	core::pose::Pose const & pose
) {


	utility::vector1< core::pack::dunbrack::RotVector > rots;
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		core::pack::dunbrack::RotVector rot;
		core::pack::dunbrack::rotamer_from_chi( pose.residue(ii), rot );
		rots.push_back( rot );
	}
	return rots;
}

utility::vector1< utility::vector1< core::Real > >
RotamerRecovery::chis_from_pose(
	core::pose::Pose const & pose
) {


	utility::vector1<utility:: vector1< core::Real > > chis;
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		utility::vector1< core::Real > chi_vec;
		for ( core::Size jj = 1; jj <= pose.residue(ii).nchi(); ++jj ) {
			chi_vec.push_back( pose.chi(jj,ii) );
		}
		chis.push_back( chi_vec );
	}
	return chis;
}


//overflowed function so that you can give this only one pose
void RotamerRecovery::get_rotamer_recovery(core::pose::Pose & native, core::pose::Pose & compared_pose){
	utility::vector1<core::pose::Pose> poses;
	poses.push_back(compared_pose);
	get_rotamer_recovery(native, poses);
}

//main body. this is where the magic happens
void RotamerRecovery::get_rotamer_recovery(core::pose::Pose & native, utility::vector1<core::pose::Pose> & compared_poses){
	utility::vector1< char > nat_bb_bins( bb_bins_from_pose(native) );
	utility::vector1< core::pack::dunbrack::RotVector > nat_rots( rots_from_pose(native) );

	utility::vector1< core::Size > bb_bins_correct( native.total_residue(), 0 );
	utility::vector1< utility::vector1< core::Size > > rots_correct(
			native.total_residue(), utility::vector1< core::Size >( native.total_residue(), 0 )
	);



	utility::vector1< utility::vector1< core::Size > > chis_correct(
			native.total_residue(), utility::vector1< core::Size >( native.total_residue(), 0 )
	);


	core::Size total(0);
	utility::vector1< core::pose::Pose >::iterator other_poses_itr( compared_poses.begin() ), other_poses_last( compared_poses.end() );
	while(  (other_poses_itr != other_poses_last ) ) {
		core::pose::Pose pose( *other_poses_itr );

		utility::vector1< char > pose_bb_bins( bb_bins_from_pose(pose) );
		utility::vector1< core::pack::dunbrack::RotVector > pose_rots( rots_from_pose(pose) );
		utility::vector1< utility::vector1< core::Real > > pose_chis( chis_from_pose(pose) );

		runtime_assert( pose_bb_bins.size() == nat_bb_bins.size() );
		runtime_assert( pose_rots.size() == pose_rots.size() );
		for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose_bb_bins[ii] == nat_bb_bins[ii] ) {
				bb_bins_correct[ii]++;
			}
			for ( core::Size jj = 1; jj <= pose_rots[ii].size(); ++jj ) {
				//std::cout << "rot(" << ii << "," << jj << ")" << std::endl;
				if ( pose_rots[ii][jj] == nat_rots[ii][jj] ) {
					rots_correct[ii][jj]++;
				}
			}
		}
		++total;
		other_poses_itr++;
	}


	// print out stats on bb_bins_correct and rots_correct
		std::ostream & output( std::cout );
		core::Size const width(12);
		core::Size const prec(4);

		output << "# total = " << ObjexxFCL::string_of(total) << std::endl;
		output
			<< ObjexxFCL::format::A( width, "resi_idx" )
			<< ObjexxFCL::format::A( width, "nat_bb_bin" )
			<< ObjexxFCL::format::A( width, "pct_bb" )
			<< ObjexxFCL::format::A( width, "nat_rot1" )
			<< ObjexxFCL::format::A( width, "pct_rot1" )
			<< ObjexxFCL::format::A( width, "nat_rot2" )
			<< ObjexxFCL::format::A( width, "pct_rot2" )
			<< ObjexxFCL::format::A( width, "nat_rot3" )
			<< ObjexxFCL::format::A( width, "pct_rot3" )
			<< ObjexxFCL::format::A( width, "nat_rot4" )
			<< ObjexxFCL::format::A( width, "pct_rot4" )
			<< std::endl;

		for ( core::Size ii = 1; ii <= native.total_residue(); ++ii ) {
			core::Real pct_bb(0.0);
			if ( bb_bins_correct[ii] > 0 ) {
				pct_bb = (
					static_cast< core::Real > ( bb_bins_correct[ii] ) /
					static_cast< core::Real > ( total )
				);
			}
			utility::vector1< core::Real > pct_natrots;
			utility::vector1< core::Size > out_nat_rots;
			for ( core::Size jj = 1; jj <= 4; ++jj ) {
				if ( jj <= nat_rots[ii].size() ) {
					core::Real pct_natrots_correct(0.0);
					if ( rots_correct[ii][jj] > 0 ) {
						pct_natrots_correct = (
							static_cast< core::Real > ( rots_correct[ii][jj] ) /
							static_cast< core::Real > ( total )
						);
					}
					pct_natrots.push_back( pct_natrots_correct );
					out_nat_rots.push_back( nat_rots[ii][jj] );
				} else {
					//std::cout << "no rot " << jj << std::endl;
					pct_natrots.push_back( 0.0 );
					out_nat_rots.push_back( 999 );
				}
			}

			//std::cout << "nat_rots[" << ii << "][1] = "
			//	<< nat_rots[ii][1] << std::endl;
			//std::cout << "nat_rots[" << ii << "].size() = "
			//	<< nat_rots[ii].size() << std::endl;

			output
				<< ObjexxFCL::format::A( width, ObjexxFCL::string_of(ii) )
				<< ObjexxFCL::format::A( width, nat_bb_bins[ii] )
				<< ObjexxFCL::format::F( width, prec, pct_bb )
				<< ObjexxFCL::format::A( width, ObjexxFCL::string_of(out_nat_rots[1]) )
				<< ObjexxFCL::format::F( width, prec, pct_natrots[1] )
				<< ObjexxFCL::format::A( width, ObjexxFCL::string_of(out_nat_rots[2]) )
				<< ObjexxFCL::format::F( width, prec, pct_natrots[2] )
				<< ObjexxFCL::format::A( width, ObjexxFCL::string_of(out_nat_rots[3]) )
				<< ObjexxFCL::format::F( width, prec, pct_natrots[3] )
				<< ObjexxFCL::format::A( width, ObjexxFCL::string_of(out_nat_rots[4]) )
				<< ObjexxFCL::format::F( width, prec, pct_natrots[4] )
				<< std::endl;
		}



}




/*
int
main( int argc, char* argv [] ) {
	// options, random initialization
	devel::init( argc, argv );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring::constraints;
	using namespace ObjexxFCL::format;
	using namespace core::import_pose::pose_stream;
	using namespace core::chemical;
	using namespace core::pack::dunbrack;

	using core::Size;
	using core::Real;
	using std::string;
	using utility::vector1;
	using ObjexxFCL::string_of;
	using core::pose::Pose;

	core::import_pose::pose_stream::MetaPoseInputStream input
		= core::import_pose::pose_stream::streams_from_cmd_line();
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();

	Pose native;
	core::io::pdb::pose_from_pdb(
		native, option[ in::file::native ]()
	);

	vector1< char > nat_bb_bins( bb_bins_from_pose(native) );
	vector1< RotVector > nat_rots( rots_from_pose(native) );

	vector1< Size > bb_bins_correct( native.total_residue(), 0 );
	vector1< vector1< Size > > rots_correct(
		native.total_residue(), vector1< Size >( native.total_residue(), 0 )
	);


	vector1< vector1< Size > > chis_correct(
		native.total_residue(), vector1< Size >( native.total_residue(), 0 )
	);




	Size total(0);
	while( input.has_another_pose() ) {
		Pose pose;
		input.fill_pose( pose, *rsd_set );

		vector1< char > pose_bb_bins( bb_bins_from_pose(pose) );
		vector1< RotVector > pose_rots( rots_from_pose(pose) );
		vector1< vector1< Real > > pose_chis( chis_from_pose(pose) );

		runtime_assert( pose_bb_bins.size() == nat_bb_bins.size() );
		runtime_assert( pose_rots.size() == pose_rots.size() );
		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
			if ( pose_bb_bins[ii] == nat_bb_bins[ii] ) {
				bb_bins_correct[ii]++;
			}
			for ( Size jj = 1; jj <= pose_rots[ii].size(); ++jj ) {
				//std::cout << "rot(" << ii << "," << jj << ")" << std::endl;
				if ( pose_rots[ii][jj] == nat_rots[ii][jj] ) {
					rots_correct[ii][jj]++;
				}
			}
		}
		++total;
	}



	//////
	/////
	/////

	// print out stats on bb_bins_correct and rots_correct
	std::ostream & output( std::cout );
	Size const width(12);
	Size const prec(4);

	output << "# total = " << string_of(total) << std::endl;
	output
		<< A( width, "resi_idx" )
		<< A( width, "nat_bb_bin" )
		<< A( width, "pct_bb" )
		<< A( width, "nat_rot1" )
		<< A( width, "pct_rot1" )
		<< A( width, "nat_rot2" )
		<< A( width, "pct_rot2" )
		<< A( width, "nat_rot3" )
		<< A( width, "pct_rot3" )
		<< A( width, "nat_rot4" )
		<< A( width, "pct_rot4" )
		<< std::endl;

	for ( Size ii = 1; ii <= native.total_residue(); ++ii ) {
		Real pct_bb(0.0);
		if ( bb_bins_correct[ii] > 0 ) {
			pct_bb = (
				static_cast< Real > ( bb_bins_correct[ii] ) /
				static_cast< Real > ( total )
			);
		}
		vector1< Real > pct_natrots;
		vector1< Size > out_nat_rots;
		for ( Size jj = 1; jj <= 4; ++jj ) {
			if ( jj <= nat_rots[ii].size() ) {
				Real pct_natrots_correct(0.0);
				if ( rots_correct[ii][jj] > 0 ) {
					pct_natrots_correct = (
						static_cast< Real > ( rots_correct[ii][jj] ) /
						static_cast< Real > ( total )
					);
				}
				pct_natrots.push_back( pct_natrots_correct );
				out_nat_rots.push_back( nat_rots[ii][jj] );
			} else {
				//std::cout << "no rot " << jj << std::endl;
				pct_natrots.push_back( 0.0 );
				out_nat_rots.push_back( 999 );
			}
		}

		//std::cout << "nat_rots[" << ii << "][1] = "
		//	<< nat_rots[ii][1] << std::endl;
		//std::cout << "nat_rots[" << ii << "].size() = "
		//	<< nat_rots[ii].size() << std::endl;

		output
			<< A( width, string_of(ii) )
			<< A( width, nat_bb_bins[ii] )
			<< F( width, prec, pct_bb )
			<< A( width, string_of(out_nat_rots[1]) )
			<< F( width, prec, pct_natrots[1] )
			<< A( width, string_of(out_nat_rots[2]) )
			<< F( width, prec, pct_natrots[2] )
			<< A( width, string_of(out_nat_rots[3]) )
			<< F( width, prec, pct_natrots[3] )
			<< A( width, string_of(out_nat_rots[4]) )
			<< F( width, prec, pct_natrots[4] )
			<< std::endl;
	}
	return 0;
} // int main( int argc, char * argv [] )
*/



}
}
}
