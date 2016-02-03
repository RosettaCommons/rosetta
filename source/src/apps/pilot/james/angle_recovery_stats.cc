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

#include <basic/options/option.hh>
#include <devel/init.hh>

#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/id/SequenceMapping.hh>

#include <protocols/jumping/util.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/sasa.hh>

#include <string>

#include <utility/excn/Exceptions.hh>

char
torsion2big_bin(
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

utility::vector1< char > get_ss( core::pose::Pose & pose ) {
	using core::Size;
	protocols::jumping::assign_ss_dssp( pose );

	utility::vector1< char > ss( pose.total_residue(), 'L' );
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		ss[ii] = pose.secstruct(ii);
	}

	return ss;
}
void print_rot_vec(
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

utility::vector1< char > bb_bins_from_pose(
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
rots_from_pose(
	core::pose::Pose const & pose
) {
	using namespace core::pack::dunbrack;
	using utility::vector1;

	vector1< RotVector > rots;
	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		RotVector rot;
		rotamer_from_chi( pose.residue(ii), rot );
		rots.push_back( rot );
	}
	return rots;
}

utility::vector1< utility::vector1< core::Real > >
chis_from_pose(
	core::pose::Pose const & pose
) {
	using core::Size;
	using core::Real;
	using utility::vector1;

	vector1< vector1< Real > > chis;
	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		vector1< Real > chi_vec;
		for ( Size jj = 1; jj <= pose.residue(ii).nchi(); ++jj ) {
			chi_vec.push_back( pose.chi(jj,ii) );
		}
		chis.push_back( chi_vec );
	}
	return chis;
}

void print_rot_vector(
	std::ostream & output,
	utility::vector1< core::pack::dunbrack::RotVector > rots
) {
	using core::Size;
	for ( Size ii = 1; ii <= rots.size(); ++ii ) {
		output << ii;
		for ( Size jj = 1; jj <= rots[ii].size(); ++jj ) {
			output << " " << rots[ii][jj];
		}
		output << std::endl;
	}
}

int
main( int argc, char* argv [] ) {
	try {

		// options, random initialization
		devel::init( argc, argv );

		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring::constraints;
		using namespace ObjexxFCL::format;
		using namespace core::import_pose::pose_stream;
		using namespace core::chemical;
		using namespace core::pack::dunbrack;
		using namespace core::sequence;

		using core::Size;
		using core::Real;
		using std::string;
		using utility::vector1;
		using ObjexxFCL::string_of;
		using core::pose::Pose;
		using core::id::SequenceMapping;

		core::import_pose::pose_stream::MetaPoseInputStream input
			= core::import_pose::pose_stream::streams_from_cmd_line();
		ResidueTypeSetCOP rsd_set( rsd_set_from_cmd_line() );

		Pose native;
		core::import_pose::pose_from_file(
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
		while ( input.has_another_pose() ) {
			Pose pose;
			input.fill_pose( pose, *rsd_set );

			vector1< char > pose_bb_bins( bb_bins_from_pose(pose) );
			vector1< RotVector > pose_rots( rots_from_pose(pose) );
			vector1< vector1< Real > > pose_chis( chis_from_pose(pose) );

			//runtime_assert( pose_bb_bins.size() == nat_bb_bins.size() );
			//runtime_assert( pose_rots.size() == pose_rots.size() );

			SequenceAlignment aln   = align_poses_naive(pose,native);
			SequenceMapping mapping = aln.sequence_mapping(2,1);

			std::cerr << "alignment: " << aln.to_string() << std::endl;
			for ( Size ii = 1; ii <= native.total_residue(); ++ii ) {
				Size const query_res(mapping[ii]);
				if ( query_res == 0 ) continue;
				std::cerr << "examining residue " << ii << " => " << query_res << std::endl;

				if ( pose_bb_bins[query_res] == nat_bb_bins[ii] ) {
					bb_bins_correct[ii]++;
				}
				bool last_rot_correct(true);

				runtime_assert( query_res <= pose_rots.size() );
				std::cout << "have " << pose_rots[query_res].size() << " rotamers" << std::endl;
				for ( Size jj = 1; jj <= pose_rots[query_res].size(); ++jj ) {
					std::cout << "looking at rot " << jj << ", residue (" << ii << "," << query_res << ")" << std::endl;
					utility::vector1< RotVector > rots( pose_rots[query_res] );
					if ( pose_rots[query_res][jj] == nat_rots[ii][jj] && last_rot_correct ) {
						rots_correct[ii][jj]++;
					} else {
						last_rot_correct = false;
					}
				} // jj in pose_rots[ii].size
				std::cerr << "finished looking at residue " << ii << " => " << query_res << std::endl;
			} // ii in total_residue
			++total;
			std::cout << "finished pose " << total << std::endl;
		}

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
			<< A( width, "sasa" )
			<< A( width, "hsasa" )
			<< std::endl;

		// calculate native per-residue sasa
		utility::vector1< core::Real > rsd_sasa( native.total_residue(), 0.0 );
		utility::vector1< core::Real > rsd_hydrophobic_sasa( native.total_residue(), 0.0 );

		core::scoring::calc_per_res_hydrophobic_sasa(
			native, rsd_sasa, rsd_hydrophobic_sasa, 1.4, false
		);

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
			// << nat_rots[ii][1] << std::endl;
			//std::cout << "nat_rots[" << ii << "].size() = "
			// << nat_rots[ii].size() << std::endl;

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
				<< F( width, prec, rsd_sasa[ii] )
				<< F( width, prec, rsd_hydrophobic_sasa[ii] )
				<< std::endl;
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
