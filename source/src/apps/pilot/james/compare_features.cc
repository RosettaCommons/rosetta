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

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <algorithm>

#include <utility/excn/Exceptions.hh>


boost::dynamic_bitset<>
get_dm_features(
	core::pose::Pose & pose,
	core::Real const dist_threshold,
	std::string const & atom_name,
	core::Size const min_seqsep = 12,
	core::Size const skip_res   = 0
) {
	using core::Size;
	using core::Real;
	using utility::vector1;
	using boost::dynamic_bitset;

	Real const dist_threshold_sq( dist_threshold * dist_threshold );

	Size const N( pose.total_residue() );
	dynamic_bitset<> features( (N*N - N) / 2 );
	Size feat_idx(1);
	Size const step_size( 1 + skip_res );
	for ( Size ii = 1; ii <= pose.total_residue(); ii += step_size ) {
	for ( Size jj = ii + min_seqsep; jj <= pose.total_residue(); jj += step_size ) {
		Real const dist_sq(
			pose.residue(ii).xyz(atom_name).distance_squared(
				pose.residue(jj).xyz(atom_name)
			)
		);

		if ( dist_sq < dist_threshold_sq ) {
			features[feat_idx] = 1;
		}
		feat_idx++;
	} // jj
	} // ii

	return features;
}

std::string dynamic_bitset_to_hex_string(
	boost::dynamic_bitset<> const & bits
) {
	using core::Size;
	using std::string;
	string str("");
	static string const hex("0123456789ABCDEF");
	for ( Size ii = 0; ii < bits.size(); ii += 4 ) {
		Size idx(0);
		if ( ii < bits.size() && bits[ii] )       idx += 8;
		if ( ii + 1 < bits.size() && bits[ii+1] ) idx += 4;
		if ( ii + 2 < bits.size() && bits[ii+2] ) idx += 2;
		if ( ii + 3 < bits.size() && bits[ii+3] ) idx += 1;
		str += hex[idx];
	}
	return str;
}

core::Real pct_features_in_common(
	boost::dynamic_bitset<> const & set1,
	boost::dynamic_bitset<> const & set2
) {
	using core::Real;

	boost::dynamic_bitset<> result = (set1 & set2);
		//( set1 &  set2) |
		//(~set1 & ~set2);
	//std::string s1, s2;
	//boost::to_string(set1,s1);
	//boost::to_string(set2,s2);
	//std::cout << s1 << " cmp " << s2 << std::endl;
	Real const dist(
		static_cast< Real >( result.count() ) /
		//static_cast< Real >( std::max( set1.count(),set2.count() ) )
		static_cast< Real >( set2.count() )
	);
	//std::cout << result.count() << " " << result.size() << std::endl;
	return dist;
}

int
main( int argc, char * argv [] ) {
	try {

	using namespace core::chemical;
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	using namespace core::import_pose::pose_stream;
	using namespace core::io::silent;

	using std::string;
	using utility::vector1;
	using core::pose::Pose;

	devel::init( argc, argv );

	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();
	MetaPoseInputStream input = streams_from_cmd_line();

  Pose native;
  core::import_pose::pose_from_file(
    native, option[ in::file::native ]()
  );

	boost::dynamic_bitset<> native_features = get_dm_features( native, 12, "CA", 5 );
	vector1< boost::dynamic_bitset<> > feature_strings;

	vector1< SilentStructOP > structs;
	std::map< boost::dynamic_bitset<>, unsigned long > counts;

	std::cout << "ca_rmsd ca_gdtmm feat_score" << std::endl;
	while( input.has_another_pose() ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );

		PROF_START( basic::STAGE1 );
		core::Real const ca_rmsd( core::scoring::CA_rmsd(pose,native) );
		PROF_STOP ( basic::STAGE1 );

		PROF_START( basic::STAGE2 );
		core::Real const ca_gdtmm( core::scoring::CA_gdtmm(pose,native) );
		PROF_STOP ( basic::STAGE2 );

		PROF_START( basic::STAGE3 );
		boost::dynamic_bitset<> features = get_dm_features( pose, 12, "CA", 5 );
		core::Real feat_score = pct_features_in_common(features,native_features);
		PROF_STOP ( basic::STAGE3 );
		//std::cout << gdtmm << " " << dist << std::endl;
		std::cout << ca_rmsd << " " << ca_gdtmm << " " << feat_score << std::endl;
	}

	basic::prof_show();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
