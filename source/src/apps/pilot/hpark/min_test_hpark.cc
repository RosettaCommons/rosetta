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


// libRosetta headers
#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/scoring/electron_density/util.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <devel/init.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>


#include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>

#include <core/scoring/constraints/util.hh>
#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/cryst.OptionKeys.gen.hh>

#include <protocols/constraints_additional/AdditionalConstraintCreators.hh>
#include <protocols/constraints_additional/COMCoordinateConstraint.hh>
#include <core/scoring/constraints/Constraints.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef


//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;
using namespace core;
using namespace protocols;
using utility::vector1;
using io::pdb::dump_pdb;


OPT_1GRP_KEY(Boolean, min, debug)
OPT_1GRP_KEY(Boolean, min, debug_verbose)
OPT_1GRP_KEY(Boolean, min, cartesian)
OPT_1GRP_KEY(String, min, minimizer)


class MinTestMover : public protocols::moves::Mover {
public:
	MinTestMover(){}
	void apply( core::pose::Pose & pose) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		using namespace protocols::moves;
		using namespace scoring;

		// steal relax flags
		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		kinematics::MoveMap mm;
		mm.set_bb  ( true );
		mm.set_chi ( true );
		mm.set_jump( true );
		mm.set( core::id::THETA, option[ OptionKeys::relax::minimize_mainchain_bond_angles ]() );
		mm.set( core::id::D, option[ OptionKeys::relax::minimize_mainchain_bond_lengths ]() );

		if ( option[ OptionKeys::relax::jump_move ].user() )
			mm.set_jump( option[ OptionKeys::relax::jump_move ]() );
		if ( option[ OptionKeys::relax::bb_move ].user() )
			mm.set_bb( option[ OptionKeys::relax::bb_move ]() );
		if ( option[ OptionKeys::relax::chi_move ].user() )
			mm.set_chi( option[ OptionKeys::relax::chi_move ]() );


		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::simple_moves::symmetry::SetupForSymmetryMoverOP symm( new protocols::simple_moves::symmetry::SetupForSymmetryMover );
			symm->apply( pose );
			core::pose::symmetry::make_symmetric_movemap( pose, mm );
		}

		/*
		// csts
		if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
				protocols::simple_moves::ConstraintSetMoverOP loadCsts( new protocols::simple_moves::ConstraintSetMover );
				loadCsts->constraint_file( core::scoring::constraints::get_cst_fa_file_option() );
				loadCsts->apply(pose);
		}
		*/

		// now add density scores from cmd line
		if ( option[ edensity::mapfile ].user() ) {
			core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn );
		}

		// set pose for density scoring if a map was input
		//   + (potentially) dock map into density
		if ( option[ edensity::mapfile ].user() || option[ cryst::mtzfile ].user() ) {
			protocols::electron_density::SetupForDensityScoringMoverOP edens
											 ( new protocols::electron_density::SetupForDensityScoringMover );
			edens->apply( pose );
		}

		Vector COM( 0.0, 0.0, 0.0 );
		utility::vector1< id::AtomID > IDs;
		for( Size ires = 1; ires <= pose.total_residue(); ++ires){
			for( Size iatm = 1; iatm <= pose.residue(ires).natoms(); ++iatm){
				//Size atmno = pose.residue(ires).atom_index(" CA ");
				//std::cout << ires << " " << atmno << std::endl;
				IDs.push_back( id::AtomID(iatm, ires) );
			}
		}

		pose::PoseOP start_pose ( new pose::Pose(pose) );

		core::scoring::constraints::ConstraintCOP cst = new protocols::constraints_additional::COMCoordinateConstraint( IDs, COM );
		pose.add_constraint( cst );

		/*
		utility::vector1< core::scoring::constraints::ConstraintCOP > csts;
		core::scoring::constraints::ConstraintCOP cst = new protocols::constraints_additional::COMCoordinateConstraint( IDs, COM );
		csts.push_back( cst );

		core::scoring::constraints::ConstraintSetOP cstset;
		cstset->add_constraints( csts );

		protocols::simple_moves::ConstraintSetMoverOP loadCsts( new protocols::simple_moves::ConstraintSetMover );
		loadCsts->constraint_set( cstset );
		loadCsts->apply(pose);
		*/

		(*scorefxn)(pose);
		scorefxn->show(std::cout, pose);

		bool debug_verbose = option[ OptionKeys::min::debug_verbose ]();
		bool debug_derivs = option[ OptionKeys::min::debug ]() | debug_verbose;
		std::string minimizer_name = option[ OptionKeys::min::minimizer ]();

		// setup the options
		if ( !option[ OptionKeys::min::cartesian ]() )  {
			if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
				core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_derivs, 200 );
				core::optimization::symmetry::SymAtomTreeMinimizer minimizer;
				std::cout << "SYMTORSION MINTEST: " << "\n";
				std::cout << "start score: " << (*scorefxn)(pose) << "\n";
				long t1=clock();
				minimizer.run( pose, mm, *scorefxn, options );
				long t2=clock();
				double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
				std::cout << "end score: " << (*scorefxn)(pose) << "\n";
				std::cout << "MIN TIME: " << time << " sec \n";
			} else {
				core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_verbose, 200 );
				core::optimization::AtomTreeMinimizer minimizer;
				std::cout << "TORSION MINTEST: " << "\n";
				std::cout << "start score: " << (*scorefxn)(pose) << "\n";
				long t1=clock();
				minimizer.run( pose, mm, *scorefxn, options );
				long t2=clock();
				double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
				std::cout << "end score: " << (*scorefxn)(pose) << "\n";
				std::cout << "MIN TIME: " << time << " sec \n";
				(*scorefxn)(pose);
				scorefxn->show(std::cout, pose);
			}
		} else {
			core::optimization::MinimizerOptions options( minimizer_name, 0.00001, true, debug_derivs, debug_verbose, 200 );
			core::optimization::CartesianMinimizer minimizer;
			std::cout << "CART MINTEST: " << "\n";
			std::cout << "start score: " << (*scorefxn)(pose) << "\n";
			long t1=clock();
			minimizer.run( pose, mm, *scorefxn, options );
			long t2=clock();
			double time = ((double)t2 - t1) / CLOCKS_PER_SEC;
			std::cout << "end score: " << (*scorefxn)(pose) << "\n";
			std::cout << "MIN TIME: " << time << " sec \n";
			(*scorefxn)(pose);
			scorefxn->show(std::cout, pose);
		}
	}
	virtual std::string get_name() const {
		return "MinTestMover";
	}
};

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols::simple_moves::symmetry;

	try{
		protocols::jd2::JobDistributor::get_instance()->go( new MinTestMover() );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}



///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	NEW_OPT(min::debug, "debug derivs?", false);
	NEW_OPT(min::debug_verbose, "debug derivs verbose?", false);
	NEW_OPT(min::cartesian, "cartesian minimization?", false);
	NEW_OPT(min::minimizer, "minimizer?", "lbfgs_armijo_nonmonotone");

	devel::init(argc, argv);

	protocols::viewer::viewer_main( my_main );
}
