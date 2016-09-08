// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/ContactMapEvaluator.cc
/// @brief
/// @author James Thompson

// libRosetta headers

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/io/silent/SilentStruct.hh>

#include <protocols/simple_filters/ContactMapEvaluator.hh>

#include <ObjexxFCL/string.functions.hh>
#include <boost/dynamic_bitset.hpp>
#include <algorithm>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_filters {

boost::dynamic_bitset<>
get_contact_features(
	core::pose::Pose const & pose,
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

	Size const N( pose.size() );
	dynamic_bitset<> features( (N*N - N) / 2 );
	Size feat_idx(1);
	Size const step_size( 1 + skip_res );
	for ( Size ii = 1; ii <= pose.size(); ii += step_size ) {
		for ( Size jj = ii + min_seqsep; jj <= pose.size(); jj += step_size ) {
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

ContactMapEvaluator::ContactMapEvaluator(
	core::pose::Pose const & native_pose,
	core::Real const max_dist,
	core::Size const min_seqsep
) :
	evaluation::SingleValuePoseEvaluator< core::Real >( "contact_map" ),
	max_dist_(max_dist),
	min_seqsep_(min_seqsep),
	native_(native_pose)
{}

void ContactMapEvaluator::apply(
	core::pose::Pose & pose,
	std::string,
	core::io::silent::SilentStruct & ss
) const {
	using ObjexxFCL::string_of;
	using core::pose::Pose;

	boost::dynamic_bitset<> native_features = get_contact_features(
		native_, max_dist_, "CA", min_seqsep_
	);
	boost::dynamic_bitset<> features = get_contact_features( pose, max_dist_, "CA", min_seqsep_ );
	core::Real const cm_score(
		pct_features_in_common(features,native_features)
	);
	std::string const output_name(
		"contact_score_" + string_of(max_dist_) + "_" + string_of(min_seqsep_)
	);
	ss.add_energy( output_name, cm_score );
}

} // simple_filter
} // protocols
