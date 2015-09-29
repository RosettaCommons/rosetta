// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// awatkins: based heavily on kdrew/oop_creator.cc

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/mm/MMLJEnergyTable.hh>
#include <core/scoring/mm/MMLJScore.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>
#include <protocols/simple_moves/hbs/HbsPatcher.hh>
#include <protocols/simple_moves/a3b_hbs/A3BHbsPatcher.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>


//Basic headers
#include <basic/resource_manager/ResourceManager.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::a3b_hbs;
using namespace protocols::simple_moves::chiral;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

void
examine_interval( Size i, Size j, Size path_distance, Real min, Real max, Real step, core::scoring::mm::MMLJEnergyTable const & mmljet ) {

	for ( Real d = min; d < max; d += step ) {

		// arbitrary factor, just so we don't spend SO much of our time zooming in
		Size mult = 2;
		Real rep1, atr1, drep1, datr1;
		Real rep2, atr2, drep2, datr2;

		mmljet.score( i, j, path_distance, d, rep1, atr1 );
		mmljet.deriv_score( i, j, path_distance, d, drep1, datr1 );

		Real d2 = d+step;
		mmljet.score( i, j, path_distance, d2, rep2, atr2 );
		mmljet.deriv_score( i, j, path_distance, d2, drep2, datr2 );

		Real ndrep = ( rep2-rep1 )/step;
		Real ndatr = ( atr2-atr1 )/step;
		Real drep = ( drep1 + drep2 ) / 2;
		Real datr = ( datr1 + datr2 ) / 2;
		if ( ndrep - drep > step*mult || ndrep - drep < -1*step*mult ) {
			std::cout << "Confirming possible issue ( " << ndrep << " vs " << drep << " ) between " << d << " " << d2 << " at " << step << " level." << std::endl;
			examine_interval( i, j, path_distance, d-step/mult, d+step, step/mult, mmljet );
		}

		if ( ndatr - datr > step*mult || ndatr - datr < -1*step*mult ) {
			std::cout << "Confirming possible issue ( " << ndrep << " vs " << drep << " ) between " << d << " " << d2 << " at " << step << " level." << std::endl;
			examine_interval( i, j, path_distance, d-step/mult, d+step, step/mult, mmljet );
		}
	}
}

int
main( int argc, char* argv[] )
{
	try {

		devel::init(argc, argv);
		std::cout.precision(10);
		std::cout.setf( std::ios::fixed, std:: ios::floatfield );

		core::scoring::mm::MMLJEnergyTable mmljet;
		core::scoring::mm::MMLJScore mmljs;

		for ( Size i = 1; i <= 104; ++i ) {
			for ( Size j = 1; j <= 104; ++j ) {
				for ( Size path_distance = 3; path_distance <= 4; ++path_distance ) {
					Real step = 1e-8;
					//Real orig_step = 1e-8;
					//Real orig_d;
					//std::cout << "Testing " << i << " " << j << " " << path_distance << " switch at " << ( 0.6*mmljs.min_dist( i, j, path_distance ) ) << " min at " << mmljs.min_dist( i, j, path_distance ) << std::endl;
					//std::cout << "(in dist_squared that's " << (0.6*mmljs.min_dist( i, j, path_distance )*0.6*mmljs.min_dist( i, j, path_distance )) << " and  " << mmljs.min_dist( i, j, path_distance )*mmljs.min_dist( i, j, path_distance )<< std::endl;
					for ( Real d = mmljs.min_dist( i, j, path_distance )*mmljs.min_dist( i, j, path_distance )*.6; d <= 64; d += step*1e5 ) {
						Real rep1, atr1, drep1, datr1;
						Real rep2, atr2, drep2, datr2;
						mmljet.score( i, j, path_distance, d, rep1, atr1 );
						mmljet.score( i, j, path_distance, d+step, rep2, atr2 );
						mmljet.deriv_score( i, j, path_distance, d, drep1, datr1 );
						mmljet.deriv_score( i, j, path_distance, d+step, drep2, datr2 );
						Real ndrep = ( rep2-rep1 )/step;
						Real ndatr = ( atr2-atr1 )/step;
						Real drep = ( drep1 + drep2 ) / 2;
						Real datr = ( datr1 + datr2 ) / 2;
						if ( ndrep - drep > 1e-6 || ndrep - drep < -1e-6 || ndatr - datr > 1e-6 || ndatr - datr < -1e-6 ) {
							if ( ndrep - drep > 1e-6 || ndrep - drep < -1e-6 ) {
								std::cout << i << " " << j << " rep" << "\t" << mmljs.min_dist( i, j, path_distance ) << "\t" << d << "\t" << ndrep << "\t" << drep << std::endl;
							}
							if ( ndatr - datr > 1e-6 || ndatr - datr < -1e-6 ) {
								std::cout << i << " " << j << " atr" << "\t" << mmljs.min_dist( i, j, path_distance ) << "\t" << d << "\t" << ndatr << "\t" << datr << std::endl;
							}
						} else if ( ndrep - drep > 1e-7 || ndrep - drep < -1e-7 || ndatr - datr > 1e-7 || ndatr - datr < -1e-7 ) {
							mmljet.score( i, j, path_distance, d+step/10, rep2, atr2 );
							mmljet.deriv_score( i, j, path_distance, d+step/10, drep2, datr2 );
							ndrep = ( rep2-rep1 )/step*10;
							ndatr = ( atr2-atr1 )/step*10;
							drep = ( drep1 + drep2 ) / 2;
							datr = ( datr1 + datr2 ) / 2;
							if ( ndrep - drep > 1e-7 || ndrep - drep < -1e-7 ) {
								std::cout << i << " " << j << " rep" << "\t" << mmljs.min_dist( i, j, path_distance ) << "\t" << d << "\t" << ndrep << "\t" << drep << std::endl;
							}
							if ( ndatr - datr > 1e-7 || ndatr - datr < -1e-7 ) {
								std::cout << i << " " << j << " atr" << "\t" << mmljs.min_dist( i, j, path_distance ) << "\t" << d << "\t" << ndatr << "\t" << datr << std::endl;
							}
						}

						//examine_interval( i, j, path_distance, d, d+step, step, mmljet );
					}
				}
			}
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main
