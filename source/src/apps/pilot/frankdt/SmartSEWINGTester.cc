// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/frankdt/SmartSEWINGTester.cc
/// @brief testing app for SEWING file IO
/// @author frankdt (frankdt@email.unc.edu)

#include <apps/pilot/frankdt/SmartSEWINGTester.hh>
#include <basic/Tracer.hh>

#include <devel/init.hh>
#include <core/pose/Pose.hh>

#include <iostream>
#include <string>
#include <stdlib.h>

#include <protocols/sewing/hashing/ModelFileReader.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>

static basic::Tracer TR( "apps.pilot.frankdt.SmartSEWINGTester" );


namespace apps {
namespace pilot {
namespace frankdt {

SmartSEWINGTester::SmartSEWINGTester():
	utility::pointer::ReferenceCount()
{

}

SmartSEWINGTester::~SmartSEWINGTester(){}

SmartSEWINGTester::SmartSEWINGTester( SmartSEWINGTester const & ) {

}



SmartSEWINGTesterOP
SmartSEWINGTester::clone() const {
	return SmartSEWINGTesterOP( new SmartSEWINGTester( *this ) );
}


} //apps
} //pilot
} //frankdt


int
main( int argc, char * argv [] )
{

	devel::init(argc,argv);

	protocols::sewing::hashing::ModelFileReaderOP modelfilereader = protocols::sewing::hashing::ModelFileReaderOP(new protocols::sewing::hashing::ModelFileReader);
	utility::vector1<protocols::sewing::data_storage::SmartSegmentOP> segment_list = modelfilereader->read_model_file("test_model_file").second;

	for ( protocols::sewing::data_storage::SmartSegmentOP current_segment : segment_list ) {
		TR<< "Testing segment " << current_segment->get_segment_id() << std::endl;
		if ( current_segment->is_n_terminus_fixed() ) {
			TR << "N terminus of segment " << current_segment->get_segment_id() << " is fixed." << std::endl;
			if ( current_segment->get_n_terminal_neighbor()->get_c_terminal_neighbor() != current_segment ) {
				TR << "Bad segment list!" << std::endl;
			}
		}
		if ( current_segment->is_c_terminus_fixed() ) {
			TR << "C terminus of segment " << current_segment->get_segment_id() << " is fixed." << std::endl;
			if ( current_segment->get_c_terminal_neighbor()->get_n_terminal_neighbor() != current_segment ) {
				TR << "Bad segment list!" << std::endl;
			}
		}

	}
	TR << "Worked!" << std::endl;
	return 0;

}




