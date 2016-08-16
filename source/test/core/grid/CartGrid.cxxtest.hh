// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>
#include <test/UTracer.hh>

#include <core/grid/CartGrid.hh>
#include <utility/json_spirit/json_spirit_reader.h>
#include <utility/json_spirit/json_spirit_writer.h>
#include <numeric/xyzVector.hh>


//Auto Headers
#include <utility/vector1.hh>

#include <sstream>
#include <iostream>

class CartGridTest : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
		grid_.read("core/grid/1dwc_3.gridbb.gz");
	}

	void tearDown() {}

	void test_grid_basics()
	{
		using namespace core::grid;


		TS_ASSERT_EQUALS(40, grid_.longestSide());

		TS_ASSERT(grid_.is_in_grid(32.838, 19.878, 25.359));
		TS_ASSERT(grid_.is_in_grid(30.493, 14.834, 23.936));
		TS_ASSERT(grid_.is_in_grid(30.826, 12.247, 27.300));

		TS_ASSERT_EQUALS(0, grid_.getValue(32.838, 19.878, 25.359));
		TS_ASSERT_EQUALS(0, grid_.getValue(30.493, 14.834, 23.936));
		TS_ASSERT_EQUALS(0, grid_.getValue(30.826, 12.247, 27.300));

		TS_ASSERT(grid_.is_in_grid(34.469, 16.389, 27.181));
		TS_ASSERT(grid_.is_in_grid(30.272, 18.897, 19.196));

		TS_ASSERT_EQUALS(1, grid_.getValue(34.469, 16.389, 27.181));
		TS_ASSERT_EQUALS(1, grid_.getValue(30.272, 18.897, 19.196));

		core::Vector const expected_base(21.829, 3.668, 14.403);
		core::Vector const expected_top(21.829+20, 3.668+20, 14.403+20);
		TS_ASSERT( 1e-2 > (expected_base - grid_.getBase()).length() );
		TS_ASSERT( 1e-2 > (expected_top - grid_.getTop()).length() );
	}

	void test_grid_json_serialization()
	{
		std::ostringstream ss;
		ss << utility::json_spirit::write(grid_.serialize());

		utility::json_spirit::mValue grid_data;
		utility::json_spirit::read(ss.str(),grid_data);
		core::grid::CartGrid<int> new_grid;
		new_grid.deserialize(grid_data.get_obj());

		//see if the data we get out is the same as the data we put in.
		TS_ASSERT(new_grid.is_in_grid(32.838, 19.878, 25.359));
		TS_ASSERT(new_grid.is_in_grid(30.493, 14.834, 23.936));
		TS_ASSERT(new_grid.is_in_grid(30.826, 12.247, 27.300));

		TS_ASSERT_EQUALS(0, new_grid.getValue(32.838, 19.878, 25.359));
		TS_ASSERT_EQUALS(0, new_grid.getValue(30.493, 14.834, 23.936));
		TS_ASSERT_EQUALS(0, new_grid.getValue(30.826, 12.247, 27.300));

		TS_ASSERT(new_grid.is_in_grid(34.469, 16.389, 27.181));
		TS_ASSERT(new_grid.is_in_grid(30.272, 18.897, 19.196));

		TS_ASSERT_EQUALS(1, new_grid.getValue(34.469, 16.389, 27.181));
		TS_ASSERT_EQUALS(1, new_grid.getValue(30.272, 18.897, 19.196));

	}

private:
	core::grid::CartGrid<int> grid_;

};

