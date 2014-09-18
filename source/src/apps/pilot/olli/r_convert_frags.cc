// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/r_frag_quality.cc
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange


#include <core/fragment/FragSet.hh>
#include <core/fragment/OrderedFragSet.hh>
#include <core/fragment/FragmentIO.hh>
#include <core/fragment/FragID_Iterator.hh>
#include <core/fragment/FrameIterator.hh>
#include <basic/options/option_macros.hh>
#include <devel/init.hh>
#include <numeric/random/random.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/util.hh>
#include <utility/excn/Exceptions.hh>
// option key includes





OPT_KEY( File, i )
OPT_KEY( File, o )
OPT_KEY( Integer, ntop )
OPT_KEY( Real, skip_rate )
OPT_KEY( File, loop )

using namespace core;
using namespace fragment;
using namespace pose;
using namespace kinematics;



using namespace basic::options;
using namespace basic::options::OptionKeys;

//using namespace ObjexxFCL::format;

void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  NEW_OPT( i, "fragment file input", "");
  NEW_OPT( o, "fragment file in new fornat", "frags.dat");
	NEW_OPT( ntop, "read top n frags " , 0 );
  NEW_OPT( skip_rate,"if >0 fragments will be stochastically chucked out", 0.0 );
	NEW_OPT( loop, "only accept frags within loop definition", "loop.dat" );
}


static thread_local basic::Tracer tr( "main" );
int main( int argc, char** argv ) {
	try{
  register_options();
  devel::init( argc, argv );

  FragSetOP orig_frags,new_frags;
  orig_frags = FragmentIO( option[ ntop ] ).read_data( option[ OptionKeys::i ]() );

	protocols::loops::Loops filter_loops;
	if ( option[ loop ].user() ) {
		filter_loops.read_loop_file( option[ loop ]() );
		FragSetOP loop_frags = new OrderedFragSet;
		protocols::loops::select_loop_frags( filter_loops, *orig_frags, *loop_frags, false );
		orig_frags = loop_frags;
	}

	new_frags = orig_frags;
  if ( option[ skip_rate ]() > 0.0  ) {
    new_frags = new OrderedFragSet;
    FragID_Iterator it = orig_frags->begin();
    FragID_Iterator eit= orig_frags->end();
    for ( ; it!=eit; ++it ) {
      if ( option[ skip_rate ]() <= 1.0 ) { //throw out fragments
				Real r = numeric::random::rg().uniform();
				if ( r > option[ skip_rate ]() ) { //keep fragment
					new_frags->add( *it );
				}
      } else {
				utility_exit_with_message( "skip_rate >= 1 ... nothing written..");
      }
    }
  }

  FragmentIO().write_data( option[ OptionKeys::o ](), *new_frags );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
