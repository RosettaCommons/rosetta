// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief

// libRosetta headers
#include <protocols/domain_assembly/AssembleLinkerMover.hh>
#include <protocols/domain_assembly/AddAssemblyConstraints.hh>
#include <protocols/domain_assembly/PostDockAssemblyScorer.hh>
//#include <protocols/domain_assembly/CombineChainsMover.hh>

#include <protocols/docking/stateless/SaneDockingProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/loops/loops_main.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/relax/FastRelax.hh>


#include <protocols/moves/Mover.hh>
#include <protocols/moves/CompositionMover.hh>

#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>

#include <core/types.hh>
#include <utility/vector1.hh>

// option key includes

#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <core/kinematics/Jump.hh>
#include <basic/options/option.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] ) {
	try {

		using namespace protocols;
		using namespace protocols::domain_assembly;
		using namespace protocols::jd2;
		using namespace protocols::moves;
		using namespace protocols::relax;
		using namespace protocols::docking::stateless;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring;

		devel::init(argc, argv);
		CompositionMoverOP container( new CompositionMover );

		get_score_function();

		// docking
		container->add_mover( MoverOP( new AddAssemblyConstraints ) );
		container->add_mover( MoverOP( new SaneDockingProtocol ) );

		// container->add_mover( new CombineChainsMover() );

		// scoring
		container->add_mover( MoverOP( new PostDockAssemblyScorer("pre_rebuild_dist") ) );

		// loop remodeling
		if ( option[ OptionKeys::loops::frag_files ].user() ) {
			using core::Size;
			Size const min_loop_size( option[ cm::min_loop_size ]() );
			utility::vector1< core::fragment::FragSetOP > frag_libs;
			protocols::loops::read_loop_fragments( frag_libs );
			MoverOP builder( new AssembleLinkerMover( "quick_ccd", min_loop_size, frag_libs ) );
			container->add_mover(builder);
			container->add_mover( MoverOP( new FastRelax( get_score_function() ) ) );
			container->add_mover( MoverOP( new PostDockAssemblyScorer( "post_rebuild_dist" ) ) );
		}

		// execution
		JobDistributor::get_instance()->go(container);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
