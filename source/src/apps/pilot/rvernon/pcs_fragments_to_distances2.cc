// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file calc_pair_stats.cc
/// @brief
/// @author Robert Vernon

#include <core/types.hh>
#include <devel/init.hh>
// AUTO-REMOVED #include <core/util/prof.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/scoring/rms_util.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/pose/Pose.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>


#include <core/kinematics/RT.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>

#include <numeric/constants.hh>

#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/robert.OptionKeys.gen.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>

#include <numeric/random/random.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzMatrix.fwd.hh>

#include <core/sequence/util.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <protocols/jumping/StrandPairing.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

///////////////////////////////////////////////////////////////////////////////

using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using ObjexxFCL::FArray2A_float;
using ObjexxFCL::FArray2D_float;
using namespace core;
using namespace pose;
using namespace conformation;
using namespace basic::options;
using namespace basic::options::OptionKeys;

utility::vector1< pose::Pose >
largest_cluster( utility::vector1< pose::Pose > lowscore_poses, utility::vector1< pose::Pose > pose_list )//, pose::Pose native )
{

	Real maxsub_filter( basic::options::option[ robert::pcs_maxsub_filter ]());//Default 0.9
	Real maxsub_rms( basic::options::option[ robert::pcs_maxsub_rmsd ]() );//Default 4.0

	Size most_hits(0), most_id(0);
	for ( Size ii = 1; ii <= lowscore_poses.size(); ++ii) {
		Size hits(0);
		for ( Size jj = 1; jj <= pose_list.size(); ++jj) {

			Real maxsub = core::scoring::CA_maxsub( lowscore_poses[ii], pose_list[jj], maxsub_rms );

			Real maxsub_coverage = maxsub / ( static_cast< Real >(lowscore_poses[ii].total_residue() ) );

			if ( maxsub_coverage >= maxsub_filter ) {
				++hits;
			}
		}

		if (hits > most_hits) {
			most_hits = hits;
			most_id = ii;
		}

		//		std::cout << "COMPARISON " << ii << " " << core::scoring::CA_maxsub( lowscore_poses[ii], native, maxsub_rms ) << " " << hits << std::endl;
	}

	utility::vector1< pose::Pose > topcluster;
	for ( Size jj = 1; jj <= pose_list.size(); ++jj) {
		Real maxsub = core::scoring::CA_maxsub( lowscore_poses[most_id], pose_list[jj], maxsub_rms );

		Real maxsub_coverage = maxsub / ( static_cast< Real >(pose_list[most_id].total_residue() ) );

		if ( maxsub_coverage >= maxsub_filter ) {
			topcluster.push_back(pose_list[jj]);
		}
	}

	return topcluster;
}


int
main( int argc, char* argv [] )
{

	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::sequence;
	using namespace utility;
	// options, random initialization
	devel::init( argc, argv );

	ResidueTypeSetCAP rsd_set =
		ChemicalManager::get_instance()->residue_type_set( chemical::CENTROID );

	std::string const pdb_file_list( option[ robert::pairdata_input_pdb_list ]() );
	//std::string const pdb_file_list("list");

	utility::io::ozstream outstream(pdb_file_list+".distances");

	utility::io::izstream infile(pdb_file_list);
	//utility::io::ozstream outfile(pdb_file_list+".FAangle.data");

	std::string fasta;
	bool has_native(false);
	pose::Pose native;
	if (basic::options::option[ basic::options::OptionKeys::in::file::native ].user()) {
		core::import_pose::pose_from_pdb( native, *rsd_set, basic::options::option[ in::file::native ]()
		);
		has_native = true;
		fasta = native.sequence();

		if ( basic::options::option[ basic::options::OptionKeys::in::file::fasta ].user()) {
			std::string fasta_file( read_fasta_file_str( option[ OptionKeys::in::file::fasta ]()[1] )[1] );
			runtime_assert( fasta == fasta_file );
		}
	} else {
		fasta = read_fasta_file_str( option[ OptionKeys::in::file::fasta ]()[1] )[1];
	}
	runtime_assert( fasta != "");


	std::cout << "SEQUENCE " << fasta << std::endl;

	Size res1, res2;
	std::string outfile_location, native_location;
	infile >> outfile_location;// >> res1 >> res2 >> native_location;

	//vector1< core::Real > temp(0, 0.0);
	//vector1< vector1< vector1< core::Real > > > vec2d;

	//bool iterate = true;
	while (infile.good()) {

		core::io::silent::SilentFileData sfd;

		sfd.read_file(outfile_location);

		//utility::vector1< Real > pcs_scores;
		//for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end();
		//			iter != end; ++iter ) {
		//	pcs_scores.push_back( iter->get_energy("pcs") );
		//}
		//std::sort( pcs_scores.begin(), pcs_scores.end() );

		//Size index = static_cast< Size >( static_cast< Real >( pcs_scores.size() ) * 0.05 );

		pose::Pose pose;

		utility::vector1< pose::Pose > allpose_list;
		utility::vector1< pose::Pose > pose_list;
		utility::vector1< pose::Pose > lowscore_list;
		utility::vector1< Real > score_list;
		//utility::vector1< Real > rms_list;

		for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end();
					iter != end; ++iter ) {
			iter->fill_pose( pose, *rsd_set );
			allpose_list.push_back( pose );
			score_list.push_back( iter->get_energy("pcs") );
			//rms_list.push_back( iter->get_energy("rms") );
		}

		utility::vector1< Real > pcs_scores( score_list );
		std::sort( pcs_scores.begin(), pcs_scores.end() );

		if (	basic::options::option[ robert::pcs_cluster_lowscoring ]() ) {//Default true
			Size index_20 = static_cast< Size >( static_cast< Real >( pcs_scores.size() ) * 0.2 );
			Size index_median = static_cast< Size >( static_cast< Real >( pcs_scores.size() ) * 0.5 );

			for ( Size pp=1; pp <= allpose_list.size(); ++pp ) {
				if (score_list[pp] <= pcs_scores[index_20]) {
					if ( lowscore_list.size() <= index_20 ) {
						lowscore_list.push_back( allpose_list[pp] );
					}
				}

				if (score_list[pp] <= pcs_scores[index_median]) {
					pose_list.push_back( allpose_list[pp] );
				}
			}
		} else {

			for ( Size pp=1; pp <= allpose_list.size(); ++pp ) {
				lowscore_list.push_back( allpose_list[pp] );
				pose_list.push_back( allpose_list[pp] );
			}
		}


		std::string outfile_fasta( pose_list[1].sequence() );
		res1 = 0;
		res2 = 0;
		Size matches(0);
		for ( int ff=0; ff < static_cast<int>(fasta.length()); ++ff ) {
			if ( (ff+outfile_fasta.length()) <= fasta.length() ) {
				bool all_match(true);
				for (int oo=0; oo < static_cast<int>(outfile_fasta.length()); ++oo ) {
					if (fasta[ff+oo] != outfile_fasta[oo]) {
						all_match = false;
						std::cout << "FASTA_MISMATCH " << ff << " " << oo << std::endl;
					}
				}

				if ( all_match == true ) {
					res1 = static_cast< Size >( ff+1 );
					res2 = static_cast< Size >( ff+outfile_fasta.length() );
					matches += 1;
				}
			}
		}

		if (matches == 0 ) {
			std::cout << "FASTA:   " << fasta << std::endl;
			std::cout << "OUTFILE: " << outfile_fasta << std::endl;
		}

		runtime_assert( matches == 1 );




		utility::vector1< pose::Pose > cluster_list( largest_cluster( lowscore_list, pose_list ) );//, fragnative ) );

		if ( has_native ) {
			pose::Pose fragnative(cluster_list[1]);

			for ( Size fn=1; fn <= fragnative.total_residue(); ++fn ) {
				fragnative.set_phi(fn, native.phi(res1+fn-1));
				fragnative.set_psi(fn, native.psi(res1+fn-1));
				fragnative.set_omega(fn, native.omega(res1+fn-1));
			}

			Real maxsub_rms( basic::options::option[ robert::pcs_maxsub_rmsd ]() );//Default 4.0
			Real natmaxsub = core::scoring::CA_maxsub( cluster_list[1], fragnative, maxsub_rms );

			std::cout << "CLUSTERING: " << outfile_location << " " << cluster_list.size() << " / " << pose_list.size() << "   " << natmaxsub << " / " << fragnative.total_residue() << std::endl;
		} else {
			std::cout << "CLUSTERING: " << outfile_location << " " << cluster_list.size() << " / " << pose_list.size() << std::endl;
		}

		if (basic::options::option[ basic::options::OptionKeys::robert::pcs_dump_cluster ] ) {
			for ( Size pp=1; pp <= cluster_list.size(); ++pp ) {
				cluster_list[pp].dump_pdb("cluster_" + ObjexxFCL::right_string_of(pp,2,'0') + ".pdb");
			}
		}

		Real cluster_coverage( static_cast< Real >( cluster_list.size() ) / static_cast< Real >( pose_list.size() ) );

		std::cout << "CLUSTER_COVERAGE = " << cluster_coverage << std::endl;

		if (cluster_coverage >= basic::options::option[ robert::pcs_cluster_coverage ]() ) {
			for ( Size pp=1; pp <= cluster_list.size(); ++pp ) {
				pose = cluster_list[pp];

				for ( Size ii=1; ii <= pose.total_residue(); ++ii ) {
					for ( Size jj=1; jj <= pose.total_residue(); ++jj ) {
						Real pose_distance( distance( pose.residue(ii).atom("CA").xyz(), pose.residue(jj).atom("CA").xyz()) );

						//Real pose_dotproduct( dot_product( (pose.residue(ii).xyz("N")-pose.residue(ii).xyz("CA")).normalize(),
						//																	 (pose.residue(jj).xyz("N")-pose.residue(jj).xyz("CA")).normalize() ) );

				    //Real pose_dotproduct2( dot_product( (pose.residue(ii).xyz("N")-pose.residue(ii).xyz("C")).normalize(),
						//																	 (pose.residue(jj).xyz("N")-pose.residue(jj).xyz("C")).normalize() ) );
				    //Real pose_dotproduct3( dot_product( (pose.residue(ii).xyz("CEN")-pose.residue(ii).xyz("CA")).normalize(),
						//																	 (pose.residue(jj).xyz("CEN")-pose.residue(jj).xyz("CA")).normalize() ) );



						outstream << "DISTANCE " << res1 << " " << res2 << " " << ii + res1 - 1 << " " << jj + res1 - 1<< " " << pose_distance;// << " " << pose_dotproduct << " " << pose_dotproduct2 << " " << pose_dotproduct3;

						if ( has_native ) {
							Real native_distance( distance( native.residue(ii+res1-1).atom("CA").xyz(), native.residue(jj+res1-1).atom("CA").xyz()) );
							//Real native_dotproduct( dot_product( (native.residue(ii).xyz("N")-native.residue(ii).xyz("CA")).normalize(),
							//																		 (native.residue(jj).xyz("N")-native.residue(jj).xyz("CA")).normalize() ) );

							//Real native_dotproduct2( dot_product( (native.residue(ii).xyz("N")-native.residue(ii).xyz("C")).normalize(),
							//																		 (native.residue(jj).xyz("N")-native.residue(jj).xyz("C")).normalize() ) );
							//Real native_dotproduct3( dot_product( (native.residue(ii).xyz("CEN")-native.residue(ii).xyz("CA")).normalize(),
							//																		 (native.residue(jj).xyz("CEN")-native.residue(jj).xyz("CA")).normalize() ) );
							outstream << " " << native_distance;// << " " << native_dotproduct << " " << native_dotproduct2 << " " << native_dotproduct3;
						}
						outstream << std::endl;

						//vec2d[ii][jj].push_back(pose_distance);
					}
				}
			}
		}

		/*
		for ( Size ii=1; ii <= pose.total_residue(); ++ii ) {
			for ( Size jj=ii+1; jj <= pose.total_residue(); ++jj ) {
				Real average(0.0);
				for ( Size nn=1; nn <= vec2d[ii][jj].size(); ++nn ) {
					average += vec2d[ii][jj][nn];
				}
				average /= static_cast< Real >( vec2d[ii][jj].size() );
				Real sdev(0.0);
				for ( Size nn=1; nn <= vec2d[ii][jj].size(); ++nn ) {
					sdev += std::pow( vec2d[ii][jj][nn] - average, 2);
				}
				sdev = std::sqrt( sdev / static_cast< Real >( vec2d[ii][jj].size() ) );

				std::cout << "AVERAGE " << ii << " " << jj << " " << average << " " << sdev;

				if ( has_native ) {
					Real native_distance( distance( native.residue(ii).atom("CA").xyz(), native.residue(jj).atom("CA").xyz()) );
					std::cout << " " << native_distance;
				}

				std::cout << std::endl;
			}
			}*/

		//std::cout << "HEYO " << outfile_location << " " << res1 << " " << res2 << " " << index << " " << pcs_scores[index] << std::endl;


		infile >> outfile_location;// >> res1 >> res2 >> native_location;
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

} // int main( int argc, char * argv [] )

//}
//}
