// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/VariantType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/ncbb/util.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>

// Filter headers
#include <numeric/random/random.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

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
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static basic::Tracer TR("test_terpenes" );

int
main( int argc, char* argv[] )
{
	try {

		//utility::vector1< core::Size > empty_vector(0);

		// init command line options
		devel::init(argc, argv);

		ResidueTypeSetCOP rts( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );

		ScoreFunctionOP score_fxn = get_score_function();

		Pose pose;
		//make_pose_from_sequence( pose, "X[isoprene:terpene_lower_full]X[isoprene]X[isoprene]X[isoprene]X[isoprene]X[isoprene]X[isoprene]X[isoprene:terpene_upper_full]", *rts );
		ResidueOP lower = conformation::ResidueFactory::create_residue( rts->name_map( "isoprene:terpene_lower_full" ) );
		ResidueOP inner = conformation::ResidueFactory::create_residue( rts->name_map( "isoprene" ) );
		ResidueOP upper = conformation::ResidueFactory::create_residue( rts->name_map( "isoprene:terpene_upper_full" ) );

		pose.append_residue_by_jump( *lower, 1, "", "", true );
		pose.append_residue_by_bond( *inner, true );
		pose.append_residue_by_bond( *inner, true );
		pose.append_residue_by_bond( *inner, true );
		pose.append_residue_by_bond( *inner, true );
		pose.append_residue_by_bond( *inner, true );
		pose.append_residue_by_bond( *upper, true );

		kinematics::MoveMapOP mm( new kinematics::MoveMap );
		mm->set_bb( true );
		mm->set_chi( true );

		pose.dump_pdb( "done.pdb" );
		for ( Size i = 1; i <= pose.total_residue(); ++i ) {
			pose.set_torsion( id::TorsionID( i, id::BB, pose.residue( i ).mainchain_torsions().size()-1 ), 60 );
			pose.set_torsion( id::TorsionID( i, id::BB, pose.residue( i ).mainchain_torsions().size() ), 180 );
		}
		pose.dump_pdb( "set.pdb" );

		Real before_score = ( *score_fxn )( pose );
		//Real best_score = before_score;
		Pose mcpose = pose;
		Size naccept = 0, nthermal = 0, nreject = 0;
		for ( Size i = 1; i <= 100000; ++i ) {
			Size resi = numeric::random::rg().uniform() * pose.total_residue() + 1;
			// skip torsion 1, which is immobile.
			Size tori = numeric::random::rg().uniform() * ( pose.residue( resi ).mainchain_torsions().size()-1 ) + 2;

			id::TorsionID torid( resi, id::BB, tori );
			Real torval = mcpose.torsion( torid );
			mcpose.set_torsion( torid, torval + numeric::random::rg().gaussian() * 10 );
			Real mcscore = ( *score_fxn )( mcpose );
			if ( mcscore < before_score ) {
				before_score = mcscore;
				//best_score = mcscore;
				pose = mcpose;
				//TR << "Accepted." << std::endl;
				++naccept;
			} else if ( Real(exp( -( mcscore - before_score ) )) > numeric::random::rg().uniform() ) {
				TR << "Accepted thermally (scores " << before_score << " then " << mcscore << "; prob " << Real(exp( -( mcscore - before_score ) ) ) << "). " << std::endl;
				before_score = mcscore;
				pose = mcpose;
				++nthermal;
			} else {
				//TR << "Rejected." << std::endl;
				++nreject;
			}
			if ( i % 1000 == 0 ) {
				TR << "Report." << std::endl;
				TR << i << " moves. " << naccept << " accepted ( " << ( 100 * Real( naccept/i ) ) << " )." << std::endl;
				TR << nthermal << " accepted thermally ( " << ( 100 * Real( nthermal/i ) ) << " ).  ";
				TR << nreject << " rejected ( " << ( 100 * Real( nreject/i ) ) << " )." << std::endl;
			}
		}

		pose.dump_pdb( "mc.pdb" );

		MinMoverOP min( new MinMover( mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, true ) );
		min->apply( pose );

		pose.dump_pdb( "min.pdb" );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}//main

