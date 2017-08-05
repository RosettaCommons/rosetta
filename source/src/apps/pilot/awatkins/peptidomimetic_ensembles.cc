// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   oop_conformation.cc
/// @brief  Creates an OOP dimer and toys around with its dihedrals.
/// @author Watkins


// Package headers
#include <devel/init.hh>

// Project headers
#include <core/types.hh>
#include <core/io/carbohydrates/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/annotated_sequence.hh>
//#include <core/pose/PDBInfo.hh>
//#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
//#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/TorsionID.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/task/TaskFactory.hh>

#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/docking/util.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>

// C++ headers
#include <iostream>
//#include <algorithm>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <protocols/ncbb/oop/OopCreatorMover.hh>
#include <protocols/simple_moves/oop/OopRandomSmallMover.hh>

using namespace std;
using namespace utility;
using namespace core;
using namespace core::chemical;
using namespace kinematics;
using namespace scoring;
using namespace import_pose;
using namespace pose;
using namespace protocols;
using namespace simple_moves;
using namespace oop;
using namespace basic;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR("PeptidomimeticEnsembles");

// application specific options
namespace peptidomimetic_ensembles {
// pert options
StringOptionKey const mimetic_sequence ( "peptidomimetic_ensembles::mimetic_sequence" );
BooleanOptionKey const cartesian_min ( "peptidomimetic_ensembles::cartesian_min" );
}

bool margin (
	core::Real val,
	core::Real comp,
	core::Real range
) {
	return ( ( val > comp - range && val < comp + range ) );
}

void gradual_minimization(
	Pose & pose,
	MinMoverOP min,
	ScoreFunctionOP sfxn
) {
	for ( Real wt = 0.02; wt <= 1; wt += 0.02 ) {
		sfxn->set_weight( core::scoring::atom_pair_constraint, wt );
		sfxn->set_weight( core::scoring::angle_constraint,     wt);
		sfxn->set_weight( core::scoring::dihedral_constraint,  wt );

		min->apply( pose );
	}
}

int
main( int argc, char *argv[] )
{
	try {
		option.add( peptidomimetic_ensembles::mimetic_sequence, "Sequence. Default ''." ).def("");
		option.add( peptidomimetic_ensembles::cartesian_min, "Cart? Default false." ).def(false);

		// initialize core
		devel::init( argc, argv );

		std::string seq = option[ peptidomimetic_ensembles::mimetic_sequence ].value();
		scoring::ScoreFunctionOP score_fxn = scoring::get_score_function();
		core::chemical::ResidueTypeSetCOP residue_set_cap = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

		Pose pose;
		core::pose::make_pose_from_sequence( pose, seq, *residue_set_cap );

		pose.conformation().detect_bonds();
		pose::ncbb::initialize_ncbbs( pose );
		kinematics::MoveMapOP pert_mm( new kinematics::MoveMap() );

		for ( Size i = 1; i <= pose.size(); ++i ) {
			pert_mm->set_bb( i, true );

			for ( Size j = 1; j <= pose.residue(i).mainchain_torsions().size(); ++j ) {
				pose.conformation().set_torsion( id::TorsionID( i, id::BB, j ), 180 );
			}
		}
		// Okay, now sample.
		pose.dump_pdb( "init.pdb" );

		protocols::simple_moves::MinMoverOP desn_min( new simple_moves::MinMover( pert_mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.0001, true ) );
		desn_min->cartesian( option[ peptidomimetic_ensembles::cartesian_min ].value() );

		gradual_minimization( pose, desn_min, score_fxn );
		pose.dump_pdb( "min.pdb" );


		score_fxn->set_weight( core::scoring::atom_pair_constraint, 1 );
		score_fxn->set_weight( core::scoring::angle_constraint, 1.0 );
		score_fxn->set_weight( core::scoring::dihedral_constraint, 1 );//10.0 );


		utility::vector1< Pose > poses;
		utility::vector1< Real > scores;
		Size ns = option[ out::nstruct ].value();
		for ( Size i = 1; i <= ns; ++i ) {
			Pose copy_pose = pose;
			for ( Size i = 1; i <= pose.size(); ++i ) {
				for ( Size j = 1; j <= pose.residue(i).mainchain_torsions().size()-1; ++j ) {
					copy_pose.conformation().set_torsion( id::TorsionID( i, id::BB, j ), numeric::random::rg().uniform() * 360 );
				}
				if ( numeric::random::rg().uniform() < 0.1 ) {
					copy_pose.conformation().set_torsion( id::TorsionID( i, id::BB, pose.residue(i).mainchain_torsions().size() ), 0 );
				} else {
					copy_pose.conformation().set_torsion( id::TorsionID( i, id::BB, pose.residue(i).mainchain_torsions().size() ), 180 );
				}

			}

			desn_min->apply( copy_pose );
			poses.push_back( copy_pose );
			scores.push_back( ( *score_fxn )( copy_pose ) );
		}

		// Find best score
		Pose ref_pose = poses[ arg_min( scores ) ];
		utility::vector1< Real > rmsds;
		for ( Size i = 1; i <= poses.size(); ++i ) {
			rmsds.push_back( CA_rmsd( poses[i], ref_pose, 1, ref_pose.size() ) );
		}
		for ( Size i = 1; i <= poses.size(); ++i ) {
			std::cout << rmsds[i] << " " << scores[i] << std::endl;
		}

		ref_pose.dump_pdb( "done.pdb" );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
