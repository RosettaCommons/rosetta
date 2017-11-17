// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief

// libRosetta headers
#include <protocols/domain_assembly/AssembleLinkerMover.hh>
#include <protocols/domain_assembly/AddAssemblyConstraints.hh>

#include <protocols/comparative_modeling/util.hh>
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/DockingLowRes.hh>
#include <protocols/docking/stateless/SaneDockingProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/loops/LoopMoverFactory.hh>
#include <protocols/loops/LoopMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/loops/loops_main.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/relax/FastRelax.hh>


#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/LinearPenaltyFunction.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>


#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/moves/CompositionMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>

// C++ headers
#include <iostream>
#include <string>

#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

// option key includes

#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

//Auto Headers
#include <basic/options/option.hh>

#include <core/scoring/saxs/SAXSEnergy.hh>

int
main( int argc, char * argv [] ) {
    try {
	using namespace protocols;
	using namespace protocols::domain_assembly;
	using namespace protocols::jd2;
	using namespace protocols::moves;
	using namespace protocols::relax;
	using namespace protocols::docking;
	using namespace protocols::docking::stateless;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;

	devel::init(argc, argv);
	CompositionMoverOP container( new CompositionMover );

	core::scoring::ScoreFunctionOP score3 = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
	score3->set_weight( core::scoring::saxs_score, 1.0 );

	// docking
	container->add_mover( new AddAssemblyConstraints );
	DockingProtocolOP docking = new SaneDockingProtocol;
//	DockingLowResOP docking = new DockingLowRes( score3 );
	docking->set_lowres_scorefxn( score3 );
	docking->set_highres_scorefxn( scorefxn );
//	docking->set_fullatom( false );
	container->add_mover( docking );

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
		container->add_mover( new FastRelax( get_score_function() ) );
	}

	// execution
	JobDistributor::get_instance()->go(container);
    } catch (utility::excn::Exception const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;

}
