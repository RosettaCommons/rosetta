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


#include <core/types.hh>
#include <devel/init.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <devel/rbsegment_Moves/RBSegmentMover.hh>
#include <devel/rbsegment_Moves/RBSegmentRelax.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

//silly using/typedef
#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>



using basic::T;
using basic::Error;
using basic::Warning;

basic::Tracer TR("rbsegmove_test");

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::scoring;

////////////////////////////////////////////////
////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// Parses command line options and inits RNG.
	devel::init(argc, argv);

	//	std::string inpdb = option[ in::file::s ]()[1];  // only care about first input file
	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, option[ OptionKeys::RBSegmentRelax::input_pdb ]().name() );

	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );

	//HelixRegisterShiftMoverOP reg_shift( new protocols::moves::HelixRegisterShiftMover( 3,28 ) );
	//	std::string function_tag("cen_std"),patch_tag("score4L");
	//	ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( function_tag, patch_tag ) );

	//	core::scoring::ScoreFunctionOP scorefxn( core::scoring::ScoreFunctionFactory::create_score_function( "cen_std" ) );

	core::scoring::ScoreFunctionOP scorefxn( new ScoreFunction() );
	scorefxn->set_weight( vdw, 1.0 );
	scorefxn->set_weight( env, 1.0 );
	scorefxn->set_weight( cbeta, 1.0 );
	scorefxn->set_weight( pair, 1.0 );

	std::cout << "Scoring the pose" << std::endl;
	(*scorefxn)(pose);
	scorefxn->show( std::cout );

	utility::vector1< protocols::moves::RBSegment > rbsegs;
	std::string filename( option[ OptionKeys::RBSegmentRelax::rb_file ]().name() );
	std::ifstream data( filename.c_str() );

	std::vector< std::vector< int > > fix_segments;
	fix_segments.clear();

	std::string line;
	while( getline( data, line ) ) {
		std::istringstream line_stream( line );
		int index, seg_begin, seg_end;
		char ss;
		std::string identifier;
		//		int seg1, seg2;
		//		line_stream >> identifier >> index >> seg_begin >> seg_end >> ss;
		line_stream >> identifier;
		if( !line_stream.fail() ) {
			if ( identifier == "SEG" ) {
				line_stream >> index >> seg_begin >> seg_end >> ss;
				TR << identifier << " " << index << " " << seg_begin << " " << seg_end << " " << ss << "\n";
				rbsegs.push_back( protocols::moves::RBSegment( index, seg_begin, seg_end, ss ) );
			}
			/*			if ( identifier == "FIX" ) {
				line_stream >> seg1 >> seg2;
				TR << identifier << " " << seg1 << " " << seg2 << "\n";
				fix_segments.push_back( std::make_pair( seg1, seg2 ) );
				} */
 		}
	}

	/*
	int i( 3 ), j( 5 ), k( 11 );
	int a( 1 ), b( 7 );
	std::vector< int > tmp_segment1;
	tmp_segment1.push_back( i );
	tmp_segment1.push_back( j );
	tmp_segment1.push_back( k );
	fix_segments.push_back( tmp_segment1 );
	std::vector< int > tmp_segment2;
	tmp_segment2.push_back( a );
	tmp_segment2.push_back( b );
	fix_segments.push_back( tmp_segment2 );
	*/
	protocols::moves::RBSegmentRelax rb_mover(scorefxn , rbsegs, fix_segments );
	rb_mover.apply( pose );

	std::string outfile = option[ out::file::o ]();
	core::io::pdb::dump_pdb( pose , outfile );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

