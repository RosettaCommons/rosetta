// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <protocols/relax/FastRelax.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMOLMover.hh>
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
#include <protocols/rigid/RigidBodyMover.hh>

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
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.functions.hh>

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

// tracer - used to replace cout
static basic::Tracer TR("EvaluateBetaMutants");

// application specific options
namespace evaluate_beta_mutants {
// pert options
}

Real evaluate_interface(
	Pose pose,
	ScoreFunctionOP const & score_fxn
) {
	Real v1 = ( *score_fxn )( pose );

	protocols::rigid::RigidBodyTransMover trans_mover( pose, 1 );
	trans_mover.step_size( 1000 );
	trans_mover.apply( pose );

	// repack min in unbound state.
	TaskFactoryOP tf( new TaskFactory() );
	tf->push_back( operation::TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
	operation::RestrictToRepackingOP rtrp( new operation::RestrictToRepacking() );
	tf->push_back( rtrp );
	simple_moves::PackRotamersMoverOP packer( new protocols::simple_moves::PackRotamersMover() );
	packer->task_factory( tf );
	packer->score_function( score_fxn );
	packer->apply( pose );

	kinematics::MoveMapOP separate_min_mm( new kinematics::MoveMap() );
	separate_min_mm->set_bb( true );
	separate_min_mm->set_chi( true );
	separate_min_mm->set_jump( 1, true );
	simple_moves::MinMoverOP separate_min( new simple_moves::MinMover( separate_min_mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, true ) );
	separate_min->apply( pose );

	Real v2 = ( *score_fxn )( pose );

	return ( v1 - v2 );
}

Real
score_ensemble(
	Pose & pose,
	ScoreFunctionOP const & score_fxn,
	std::string const & name
) {
	protocols::relax::FastRelaxOP frlx( new protocols::relax::FastRelax( score_fxn ) );
	kinematics::MoveMapOP mm( new kinematics::MoveMap );
	mm->set_bb( true );
	mm->set_chi( true );
	mm->set_jump( true );
	frlx->set_movemap( mm );

	// We actually want the multiple trials to occur on an unrelaxed pose
	// Otherwise we don't actually sample much.
	//frlx->apply( pose );

	Real score = evaluate_interface( pose, score_fxn );
	Real sum = 0;
	Size ntrials = 3;
	TR << "Starting score of " << name << " pose is: " << score << std::endl;

	utility::vector1< Real > scores;
	for ( Size ii = 1; ii <= ntrials; ++ii ) {
		Pose copypose = pose;

		frlx->apply( copypose );

		bool acceptable = true;

		// filters...

		// reset if failed filters
		if ( !acceptable ) {
			--ii;
			continue;
		}

		Real newscore = evaluate_interface( copypose, score_fxn );

		scores.push_back( newscore );

		sum += newscore;
		if ( newscore < score ) {
			score = newscore;
			copypose.dump_pdb( name );
		}
		TR << "Trial " << ii << " score " << newscore << std::endl;
	}
	TR << "Other stats: average was " << ( sum / Real(ntrials) ) << std::endl;

	// Av of best ten
	if ( ntrials > 10 ) {
		utility::vector1< Size > indices_best_n( 10 );
		utility::arg_least_several( scores, indices_best_n );
		Real sum_10 = 0;
		for ( Size i = 1; i <= 10; ++i ) {
			sum_10 += scores[ indices_best_n[ i ] ];
		}

		TR << "Average of best 10: " << sum_10/10.0 << std::endl;
	}

	return score;
}

int
main( int argc, char* argv[] )
{
	try {

		// init command line options
		//you MUST HAVE THIS CALL near the top of your main function, or your code will crash when you first access the command line options
		devel::init(argc, argv);

		using namespace core::scoring::constraints;

		// create score function
		scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );

		score_fxn->set_weight( fa_rep, 0.3 );

		Pose pose;
		std::string file = option[ in::file::s ].value()[1];
		core::import_pose::pose_from_file( pose, file, core::import_pose::PDB_file );

		bool is_1ycr = ( file.find("1ycr") != std::string::npos );

		bool is_mdm4 = ( file.find("2z5t") != std::string::npos );

		chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

		// Hard coded because who cares.
		// For mdm2 based on 1ycr
		Size const leu_pos = is_mdm4 ? 96  : ( is_1ycr ? 88 : 91 );
		Size const trp_pos = is_mdm4 ? 99  : ( is_1ycr ? 91 : 94 );
		Size const phe_pos = is_mdm4 ? 102 : ( is_1ycr ? 94 : 97 );

		Residue const leu = pose.residue( leu_pos );
		Residue const trp = pose.residue( trp_pos );
		Residue const phe = pose.residue( phe_pos );

		kinematics::MoveMapOP trp_mm( new kinematics::MoveMap );
		trp_mm->set_bb( true, trp_pos );
		trp_mm->set_chi( true, trp_pos );
		protocols::simple_moves::MinMoverOP trp_min( new protocols::simple_moves::MinMover( trp_mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.001, true ) );

		kinematics::MoveMapOP leu_mm( new kinematics::MoveMap );
		leu_mm->set_bb( true, leu_pos );
		leu_mm->set_chi( true, leu_pos );
		protocols::simple_moves::MinMoverOP leu_min( new protocols::simple_moves::MinMover( leu_mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.001, true ) );

		kinematics::MoveMapOP phe_mm( new kinematics::MoveMap );
		phe_mm->set_bb( true, phe_pos );
		phe_mm->set_chi( true, phe_pos );
		protocols::simple_moves::MinMoverOP phe_min( new protocols::simple_moves::MinMover( phe_mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.001, true ) );


		// This starting pose is the first
		TR << "Native." << std::endl;
		Real b53_8_score = score_ensemble( pose, score_fxn, "native_"+file );

		// Next, meta-trifluoromethyl.
		TR << "Mutant 1, meta-trifluoromethyl." << std::endl;
		Pose b53_12_pose = pose;
		ResidueType const & tfm_rt = restype_set->name_map( "B3F:B3F-CE1-trifluoromethylated" );
		ResidueOP tfm = ResidueFactory::create_residue( tfm_rt, pose.residue( trp_pos ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( trp_pos ), *tfm, b53_12_pose.conformation() );
		b53_12_pose.conformation().replace_residue( trp_pos, *tfm, false );
		// Minimize the substituted residue.
		trp_min->apply( b53_12_pose );
		Real b53_12_score = score_ensemble( b53_12_pose, score_fxn, "meta_trifluoromethyl_phe_"+file );

		// Next, meta-trifluoromethyl.
		TR << "Mutant 1b, other meta-trifluoromethyl." << std::endl;
		Pose b53_12b_pose = pose;
		ResidueType const & tfmb_rt = restype_set->name_map( "B3F:B3F-CE2-trifluoromethylated" );
		ResidueOP tfmb = ResidueFactory::create_residue( tfmb_rt, pose.residue( trp_pos ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( trp_pos ), *tfmb, b53_12b_pose.conformation() );
		b53_12b_pose.conformation().replace_residue( trp_pos, *tfmb, false );
		trp_min->apply( b53_12b_pose );
		Real b53_12b_score = score_ensemble( b53_12b_pose, score_fxn, "meta_trifluoromethyl_pheb_"+file );

		// Next, ch2-chloro trp
		TR << "Mutant 2, ch2-chloro trp." << std::endl;
		Pose b53_13_pose = pose;
		ResidueType const & ch2_rt = restype_set->name_map( "B3W:B3W-CH2-chlorinated" );
		ResidueOP ch2 = ResidueFactory::create_residue( ch2_rt, pose.residue( trp_pos ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( trp_pos ), *ch2, b53_13_pose.conformation() );
		b53_13_pose.conformation().replace_residue( trp_pos, *ch2, false );
		trp_min->apply( b53_13_pose );
		Real b53_13_score = score_ensemble( b53_13_pose, score_fxn, "chloro_trp_"+file );

		// Next, phe
		TR << "Mutant 3, phe." << std::endl;
		Pose b53_14_pose = pose;
		ResidueType const & phe2_rt = restype_set->name_map( "B3F" );
		ResidueOP phe2 = ResidueFactory::create_residue( phe2_rt, pose.residue( trp_pos ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( trp_pos ), *phe2, b53_14_pose.conformation() );
		b53_14_pose.conformation().replace_residue( trp_pos, *phe2, false );
		trp_min->apply( b53_14_pose );
		Real b53_14_score = score_ensemble( b53_14_pose, score_fxn, "phe_"+file );

		// Next, meta-chloro
		TR << "Mutant 4, m-chlorophe." << std::endl;
		Pose b53_15_pose = pose;
		ResidueType const & mcl_rt = restype_set->name_map( "B3F:B3F-CE1-chlorinated" );
		ResidueOP mcl = ResidueFactory::create_residue( mcl_rt, pose.residue( trp_pos ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( trp_pos ), *mcl, b53_15_pose.conformation() );
		b53_15_pose.conformation().replace_residue( trp_pos, *mcl, false );
		trp_min->apply( b53_15_pose );
		Real b53_15_score = score_ensemble( b53_15_pose, score_fxn, "meta_chloro_phe_"+file );

		// Next, meta para dichloro
		TR << "Mutant 5, m,p dichloro phe." << std::endl;
		Pose b53_16_pose = pose;
		ResidueType const & mpcl_rt = restype_set->name_map( "B3F:B3F-CE1-chlorinated:B3F-CZ-chlorinated" );
		ResidueOP mpcl = ResidueFactory::create_residue( mpcl_rt, pose.residue( trp_pos ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( trp_pos ), *mpcl, b53_16_pose.conformation() );
		b53_16_pose.conformation().replace_residue( trp_pos, *mpcl, false );
		trp_min->apply( b53_16_pose );
		Real b53_16_score = score_ensemble( b53_16_pose, score_fxn, "meta_para_dichloro_phe_"+file );

		// Next, leu -> ile from 12
		TR << "Mutant 6, ile and also m-trifluoromethyl." << std::endl;
		Pose b53_17_pose = b53_12_pose;
		ResidueType const & ile_rt = restype_set->name_map( "B3I" );
		ResidueOP ile = ResidueFactory::create_residue( ile_rt, pose.residue( leu_pos ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( leu_pos ), *ile, b53_17_pose.conformation() );
		b53_17_pose.conformation().replace_residue( leu_pos, *ile, false );
		leu_min->apply( b53_17_pose );
		Real b53_17_score = score_ensemble( b53_17_pose, score_fxn, "ile_"+file );

		// Next, phe -> 4-fl phe from 16
		TR << "Mutant 7, m,p dichloro phe and also 4-fl phe." << std::endl;
		Pose b53_18_pose = b53_16_pose;
		ResidueType const & fp_rt = restype_set->name_map( "B3F:B3F-CZ-fluorinated" );
		ResidueOP fp = ResidueFactory::create_residue( fp_rt, pose.residue( phe_pos ), pose.conformation() );
		core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( pose.residue( phe_pos ), *fp, b53_18_pose.conformation() );
		b53_18_pose.conformation().replace_residue( phe_pos, *fp, false );
		phe_min->apply( b53_18_pose );
		Real b53_18_score = score_ensemble( b53_18_pose, score_fxn, "para_fluoro_phe_"+file );

		TR << "Starting score beta-53-8 for mdm2: " << b53_8_score << std::endl;
		TR << "meta trifluoromethyl phe beta-53-12 for mdm2: " << b53_12_score << std::endl;
		TR << "meta trifluoromethyl phe beta-53-12 for mdm2: " << b53_12b_score << std::endl;
		TR << "CH2-chloro trp beta-53-13 for mdm2: " << b53_13_score << std::endl;
		TR << "Phe beta-53-14 for mdm2: " << b53_14_score << std::endl;
		TR << "meta chloro phe beta-53-15 for mdm2: " << b53_15_score << std::endl;
		TR << "meta para dichloro phe beta-53-16 for mdm2: " << b53_16_score << std::endl;
		TR << "meta trifluoromethyl phe ile beta-53-17 for mdm2: " << b53_17_score << std::endl;
		TR << "para fluoro phe meta para dichloro phe beta-53-18 for mdm2: " << b53_18_score << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main
