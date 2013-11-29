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

// AUTO-REMOVED #include <protocols/comparative_modeling/util.hh>
#include <protocols/docking/stateless/SaneDockingProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMoverFactory.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover.hh>
// AUTO-REMOVED #include <protocols/loops/Loop.hh>
// AUTO-REMOVED #include <protocols/loops/Loops.hh>

#include <protocols/loops/loops_main.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/relax/FastRelax.hh>


// AUTO-REMOVED #include <core/chemical/ChemicalManager.fwd.hh>

// AUTO-REMOVED #include <core/scoring/func/Func.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Constraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/func/LinearPenaltyFunction.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>


#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/moves/CompositionMover.hh>
// AUTO-REMOVED #include <protocols/simple_moves/ConstraintSetMover.hh>

#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>

// AUTO-REMOVED #include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

// option key includes

#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/constraints.OptionKeys.gen.hh>

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

	getScoreFunction();

	// docking
	container->add_mover( new AddAssemblyConstraints );
	container->add_mover( new SaneDockingProtocol );

//	container->add_mover( new CombineChainsMover() );

	// scoring
	container->add_mover( new PostDockAssemblyScorer("pre_rebuild_dist") );

	// loop remodeling
	if ( option[ OptionKeys::loops::frag_files ].user() ) {
		using core::Size;
		Size const min_loop_size( option[ cm::min_loop_size ]() );
		utility::vector1< core::fragment::FragSetOP > frag_libs;
		protocols::loops::read_loop_fragments( frag_libs );
		MoverOP builder(
			new AssembleLinkerMover( "quick_ccd", min_loop_size, frag_libs )
		);
		container->add_mover(builder);
		container->add_mover( new FastRelax( getScoreFunction() ) );
		container->add_mover( new PostDockAssemblyScorer( "post_rebuild_dist" ) );
	}

	// execution
	JobDistributor::get_instance()->go(container);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	
	return 0;
}
