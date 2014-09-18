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
/// testing sixdtree functions, probably will convert this to a unit test later

#include <basic/Tracer.hh>
#include <numeric/geometry/hashing/SixDHasher.hh>
#include <protocols/match/Hit.hh>
// C++ headers
#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>


using basic::T;
using namespace core;
using namespace protocols::match;
using namespace numeric::geometry;
static thread_local basic::Tracer TR( "main" );


int
main( int argc, char * argv [] )
{
    try {

	// initialize core
	const int HASH_POSITION_GRID_SIZE= 75;
	devel::init(argc, argv);
	SixDCoordinateBinnerOP  hash_;

	Size loop_size = 10;
	BoundingBox< Vector > bounding_box( core::Vector( -HASH_POSITION_GRID_SIZE,
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

	hash_ = new SixDCoordinateBinner( bounding_box, euler_offsets, bin_widths );
	hash_->tree_init(15);
	Real6 center(0);
	center[4]=60;
	center[5]=60;
	center[6]=60;
/*std::vector < boost::uint64_t > offset_list = hash_->radial_bin_index(0,2,center);
	for( std::vector< boost::uint64_t >::iterator itr = offset_list.begin(); itr != offset_list.end(); itr++ ) {
		TR << *itr << std::endl;
	}
	*/
    } catch ( utility::excn::EXCN_Base const & e ) {
			std::cerr << "caught exception " << e.msg() << std::endl;
			return -1;
    }
    return 0;
}

