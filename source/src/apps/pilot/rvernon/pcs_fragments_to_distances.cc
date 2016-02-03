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
#include <core/io/pdb/pdb_writer.hh>


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

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
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
largest_cluster( utility::vector1< pose::Pose > pose_list )
{

	//mjo commenting out 'maxsub_filter' because it is unused and causes a warning
	//Real maxsub_filter( basic::options::option[ robert::pcs_maxsub_filter ]());//Default 0.9
	Real maxsub_rms( basic::options::option[ robert::pcs_maxsub_rmsd ]() );//Default 4.0

	Size most_hits(0), most_id(0);
	for ( Size ii = 1; ii <= pose_list.size(); ++ii) {
		Size hits(0);
		for ( Size jj = 1; jj <= pose_list.size(); ++jj) {

			Real maxsub = core::scoring::CA_maxsub( pose_list[ii], pose_list[jj], maxsub_rms );

			Real maxsub_coverage = maxsub / ( static_cast< Real >(pose_list[ii].total_residue() ) );

			if ( maxsub_coverage >= 0.9 ) {
				++hits;
			}
		}

		if (hits > most_hits) {
			most_hits = hits;
			most_id = ii;
		}
	}

	utility::vector1< pose::Pose > topcluster;
	for ( Size jj = 1; jj <= pose_list.size(); ++jj) {
		Real maxsub = core::scoring::CA_maxsub( pose_list[most_id], pose_list[jj], maxsub_rms );

		Real maxsub_coverage = maxsub / ( static_cast< Real >(pose_list[most_id].total_residue() ) );

		if ( maxsub_coverage >= 0.9 ) {
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

	bool has_native(false);
	pose::Pose native;
	if (basic::options::option[ basic::options::OptionKeys::in::file::native ].user()) {
		core::import_pose::pose_from_file( native, *rsd_set, basic::options::option[ in::file::native ]()
		);
		has_native = true;
	}

	Size res1, res2;
	std::string outfile_location, native_location;
	infile >> outfile_location >> res1 >> res2 >> native_location;

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

		utility::vector1< pose::Pose > pose_list;
		utility::vector1< pose::Pose > forcluster_list;
		utility::vector1< Real > score_list;

		for ( core::io::silent::SilentFileData::iterator iter = sfd.begin(), end = sfd.end();
					iter != end; ++iter ) {
			iter->fill_pose( pose, *rsd_set );
			pose_list.push_back( pose );
			score_list.push_back( iter->get_energy("pcs") );
		}

		utility::vector1< Real > pcs_scores( score_list );
		std::sort( pcs_scores.begin(), pcs_scores.end() );
		Size index = static_cast< Size >( static_cast< Real >( pcs_scores.size() ) * 0.5 );

		for ( Size pp=1; pp <= pose_list.size(); ++pp ) {
			if (score_list[pp] <= pcs_scores[index]) {
				forcluster_list.push_back( pose_list[pp] );
			}
		}

		utility::vector1< pose::Pose > cluster_list( largest_cluster( forcluster_list ) );

		pose::Pose fragnative;
		core::import_pose::pose_from_file(fragnative, native_location, core::import_pose::PDB_file);
		Real natmaxsub = core::scoring::CA_maxsub( cluster_list[1], fragnative );

		std::cout << "CLUSTERING: " << outfile_location << " " << cluster_list.size() << " / " << forcluster_list.size() << "   " << natmaxsub << " / " << fragnative.total_residue() << std::endl;


		if (basic::options::option[ basic::options::OptionKeys::robert::pcs_dump_cluster ] ) {
			for ( Size pp=1; pp <= cluster_list.size(); ++pp ) {
				cluster_list[pp].dump_pdb("cluster_" + ObjexxFCL::right_string_of(pp,2,'0') + ".pdb");
			}
		}

		pose_list = cluster_list;


		for ( Size pp=1; pp <= pose_list.size(); ++pp ) {
			pose = pose_list[pp];

			//if (score_list[pp] <= pcs_scores[index]) {
			if (pose_list.size() >= 10 ) {

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


		infile >> outfile_location >> res1 >> res2 >> native_location;
	}
} // int main( int argc, char * argv [] )


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

//}
//}
