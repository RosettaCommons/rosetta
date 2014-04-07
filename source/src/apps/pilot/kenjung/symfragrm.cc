// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ken Jung
/// @brief
//	Takes as input an index and offset and returns the corresponding frag in a pdb
//	OR a bin key and returns the pdbs of all the frags in the bin
//	should take only one loopsize

// libRosetta headers
//#include <basic/options/keys/in.OptionKeys.gen.hh>
//j#include <basic/options/keys/wum.OptionKeys.gen.hh>
#include <basic/options/keys/lh.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>

#include <core/io/pdb/file_data.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>


using basic::T;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;
using namespace protocols::loophash;
using namespace core::io::pdb;
static basic::Tracer TR("main");



int
main( int argc, char * argv [] )
{
    try {

	// initialize core
	devel::init(argc, argv);
	utility::vector1 < std::string > pdblist;
	pdblist = option[lh::symfragrm::pdblist]();
	if ( pdblist.size() == 0 ) {
				 TR << "pdblist empty" << std::endl;
					return 1;
	}
	// Set up hash bins that mimic loophash

	core::Size loop_size = 10;

	TR.Info << "Setting up hash_: Size:  " << loop_size << std::endl;

		const int HASH_POSITION_GRID_SIZE= 75;

	BoundingBox bounding_box( core::Vector( -HASH_POSITION_GRID_SIZE,
				-HASH_POSITION_GRID_SIZE,
				-HASH_POSITION_GRID_SIZE),
			core::Vector( HASH_POSITION_GRID_SIZE,
				HASH_POSITION_GRID_SIZE,
				HASH_POSITION_GRID_SIZE ) );
	protocols::match::Size3 euler_offsets;
	euler_offsets[1] = 0;
	euler_offsets[2] = 0;
	euler_offsets[3] = 0;
	protocols::match::Real6 bin_widths;

	core::Real space_multiplier = 1.0;
  core::Real angle_multiplier =  15.0/6.0;

	bin_widths[1] = space_multiplier*loop_size;
	bin_widths[2] = space_multiplier*loop_size;
	bin_widths[3] = space_multiplier*loop_size;
	bin_widths[4] = angle_multiplier*loop_size;
	bin_widths[5] = angle_multiplier*loop_size;
	bin_widths[6] = angle_multiplier*loop_size;

	protocols::match::SixDCoordinateBinnerOP hash_ = new protocols::match::SixDCoordinateBinner( bounding_box, euler_offsets, bin_widths );


	// now go through and bin the input frags
	core::pose::Pose decoy;
	protocols::match::Real6 rt;
	for ( core::Size i = 1; i <= pdblist.size(); i++ ) {
					build_pose_from_pdb_as_is( decoy, pdblist[i] );
					get_rt_over_leap_fast(decoy, 2,12, rt);
					TR << rt[1] << " " << rt[2] << " " << rt[3] << " " << rt[4] << " " << rt[5] << " " << rt[6] << std::endl;
					boost::uint64_t bin_index = hash_->bin_index( rt );
					TR << bin_index << std::endl;
	}

    } catch ( utility::excn::EXCN_Base const & e ) {
			std::cerr << "caught exception " << e.msg() << std::endl;
			return -1;
    }
    return 0;
}

