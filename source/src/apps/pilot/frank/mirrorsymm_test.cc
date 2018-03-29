// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymAtomTreeMultifunc.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/atom_tree_minimize.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/util.hh>

#include <core/scoring/electron_density/util.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <devel/init.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/jd2/JobDistributor.hh>


#include <core/io/pdb/pdb_writer.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <protocols/constraint_movers/ConstraintSetMover.hh>

#include <core/scoring/constraints/util.hh>

#include <utility/excn/Exceptions.hh>


#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef


//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>
using basic::Error;
using basic::Warning;
using namespace core;
using namespace protocols;
using utility::vector1;



static basic::Tracer TR( "min_test" );

bool rama_list_pred( const std::pair < core::Size, core::Real > &left, const std::pair < core::Size, core::Real > &right ) {
	return left.second > right.second;
}


class MirrorSymmTest : public protocols::moves::Mover {
public:
	MirrorSymmTest(){}

	void apply( core::pose::Pose & pose) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core;
		using namespace core::pack;
		using namespace core::conformation::symmetry;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;

		using namespace protocols::moves;
		using namespace scoring;

		// steal relax flags
		scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
		kinematics::MoveMap mm;
		mm.set_bb  ( true );
		mm.set_chi ( true );
		mm.set_jump( true );
		mm.set( core::id::THETA, option[ OptionKeys::relax::minimize_bond_angles ]() );
		mm.set( core::id::D, option[ OptionKeys::relax::minimize_bond_lengths ]() );

		if ( option[ OptionKeys::relax::jump_move ].user() ) {
			mm.set_jump( option[ OptionKeys::relax::jump_move ]() );
		}
		if ( option[ OptionKeys::relax::bb_move ].user() ) {
			mm.set_bb( option[ OptionKeys::relax::bb_move ]() );
		}
		if ( option[ OptionKeys::relax::chi_move ].user() ) {
			mm.set_chi( option[ OptionKeys::relax::chi_move ]() );
		}

		// test 1: setup symmetric system
		TR << "*** setup_for_symmetry ***" <<std::endl;
		protocols::symmetry::SetupForSymmetryMoverOP symm( new protocols::symmetry::SetupForSymmetryMover );
		symm->apply( pose );
		core::pose::symmetry::make_symmetric_movemap( pose, mm );
		TR << "*** setup_for_symmetry done ***" <<std::endl;

		pose.dump_pdb("test1.pdb");

		// set lattice
		SymmetricConformation & SymmConf ( dynamic_cast<SymmetricConformation &> ( pose.conformation()) );
		Size Ajump_=0, Bjump_=0, Cjump_=0, SUBjump_=0;
		for ( Size i=1; i<=pose.fold_tree().num_jump(); ++i ) {
			std::string jumpname = SymmConf.Symmetry_Info()->get_jump_name( i );
			if ( jumpname == "A" ) Ajump_=i;
			else if ( jumpname == "B" ) Bjump_=i;
			else if ( jumpname == "C" ) Cjump_=i;
			else if ( jumpname == "SUB" ) SUBjump_=i;
		}
		core::id::TorsionID Ax(Ajump_, id::JUMP, 1);
		TR << "*** set_dof ***" <<std::endl;
		core::id::DOF_ID Adofx = pose.conformation().dof_id_from_torsion_id( Ax );
		pose.set_dof( Adofx, 0 );
		TR << "*** set_dof_done ***" <<std::endl;
	}

	virtual std::string get_name() const {
		return "MirrorSymmTest";
	}
};

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols::symmetry;

	try{
		protocols::jd2::JobDistributor::get_instance()->go( protocols::moves::MoverOP( new MirrorSymmTest() ) );
	} catch (utility::excn::Exception& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
