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
/// @author Nobuyasu Koga

// Project headers
#include <core/types.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/forge/build/Interval.hh>

#include <core/pose/Pose.hh>
#include <protocols/moves/DsspMover.hh>
#include <core/import_pose/import_pose.hh>

#include <core/fragment/Frame.hh>
#include <core/fragment/IndependentBBTorsionSRFD.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/picking_old/vall/util.hh>

// Utility Headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <devel/init.hh>
#include <utility/excn/Exceptions.hh>

static thread_local basic::Tracer TR( "pick_fragments" );

typedef core::Size Size;
typedef std::string String;

using protocols::forge::build::Interval;
using core::fragment::FragmentIO;
using core::fragment::FrameList;
using core::fragment::FragSetOP;
using core::fragment::ConstantLengthFragSet;


class ThisApplication  {
public:
	ThisApplication(){};
	static void register_options();
};

OPT_KEY( File, blueprint )
OPT_KEY( File, output )
OPT_KEY( Integer, size )
OPT_KEY( Integer, nfrags )
OPT_KEY( Boolean, use_abego )

void ThisApplication::register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	NEW_OPT( blueprint, "blueprint file", "" );
	NEW_OPT( output, "output filename", "fragset" );
	NEW_OPT( size, "size of fragment to pick ", 9 );
	NEW_OPT( nfrags, "number of fragments", 200 );
	NEW_OPT( use_abego, "use abego description if blueprint has abego", false );
}


FrameList pick_fragments(
			   String const & complete_ss,
			   String const & complete_aa,
			   utility::vector1< String > const & complete_abego,
			   Interval const & interval,
			   Size const frag_length,
			   Size const n_frags )
{
	using core::fragment::Frame;
	using core::fragment::FrameOP;
	using core::fragment::IndependentBBTorsionSRFD;

	using core::fragment::picking_old::vall::pick_fragments;

	FrameList frames;

	for ( Size j = 0, je = interval.length(); j < je; ++j ) {
		TR << "picking " << n_frags << " " << frag_length << "-mers for position " << ( interval.left + j ) << std::endl;

		String ss_sub = complete_ss.substr( interval.left + j - 1, frag_length );
		if ( ss_sub.length() < frag_length ) {
			ss_sub.append( frag_length - ss_sub.length(), 'D' );
		}

		String aa_sub;
		if ( !complete_aa.empty() ) {
			aa_sub = complete_aa.substr( interval.left + j - 1, frag_length );
			if ( aa_sub.length() < frag_length ) {
				aa_sub.append( frag_length - aa_sub.length(), '.' );
			}
		} else {
			aa_sub = "";
		}

		utility::vector1< String > abego_sub;
		if ( complete_abego.size() > 0 ) {
			runtime_assert( complete_ss.length() == complete_abego.size() );
			Size pos( 1 );
			abego_sub.resize( frag_length );
			for( Size ii = interval.left + j; ii <= interval.left + j + frag_length - 1; ++ii, ++pos ) {
				if ( ii > complete_abego.size() ) {
					abego_sub[ pos ] = "X";
				} else {
					abego_sub[ pos ] = complete_abego[ ii ];
				}
			}
		} else {
			abego_sub.clear(); // make sure it is empty
		}

		FrameOP frame = new Frame( interval.left + j, frag_length );

		frame->add_fragment( pick_fragments( ss_sub, aa_sub, abego_sub, n_frags, true, IndependentBBTorsionSRFD() ) );

		frames.push_back( frame );
	}

	return frames;
}


int
main( int argc, char * argv [] )
{
	try{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	ThisApplication::register_options();
  devel::init(argc, argv);

	Size frag_size = option[ size ];
	Size n_frags = option[ nfrags ];


	if( option[ blueprint ].user() && option[ in::file::s ].user() ){
		TR << "You need to choose either options of -blueprint or -s " << std::endl;
		runtime_assert( false );
	}


	Size naa( 0 );
	String ss( "" );
    utility::vector1< String >abego;

	if( option[ blueprint ].user() ){
		// read secondary structure from blueprint
		protocols::jd2::parser::BluePrintOP blue = new protocols::jd2::parser::BluePrint( option[ blueprint ] );
		naa = blue->total_residue();
		ss = blue->secstruct();
		if( option[ use_abego ].user() ) {
			abego = blue->abego();
		}
	}else if( option[ in::file::s ].user() ){
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, option[ in::file::s ].value().at( 1 ) );
		protocols::moves::DsspMover dsm;
		dsm.apply( pose );
		naa = pose.total_residue();
		ss = pose.secstruct();
	}else{
		TR << "You need to choose either options of -blueprint or -s " << std::endl;
		runtime_assert( false );
	}

	// define the region to pick frag
	Interval ival( 1, naa );

	// pick fragments
	FragSetOP fragset = new ConstantLengthFragSet( frag_size );
	fragset->add( pick_fragments( ss, "", abego, ival, frag_size, n_frags ) );

	// output fragments
	FragmentIO().write_data( option[ output ], *fragset );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

