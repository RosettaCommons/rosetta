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
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Patch.hh>
#include <core/chemical/PatchOperation.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/MMAtomType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Mover headers
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/ncbb/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <protocols/residue_optimization/MetapatchEnumeration.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <utility/vector1.functions.hh>

// Filter headers
#include <basic/MetricValue.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
//#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
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
using namespace protocols::rigid;
using namespace protocols::toolbox;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static basic::Tracer TR("test_metapatch" );

namespace metapatch {
// pert options
IntegerOptionKey const first_ligand_res ( "metapatch::first_ligand_res" );
IntegerOptionKey const last_ligand_res ( "metapatch::last_ligand_res" );
IntegerOptionKey const scoring_mode ( "metapatch::scoring_mode" );
RealOptionKey const score_threshold ( "metapatch::score_threshold" );
IntegerOptionKey const top_n ( "metapatch::top_n" );
BooleanOptionKey const also_pack ( "metapatch::also_pack" );
RealOptionKey const taboo_threshold ( "metapatch::taboo_threshold" );
}

int
main( int argc, char* argv[] )
{
	try {

		using namespace protocols::residue_optimization;

		//utility::vector1< core::Size > empty_vector(0);
		option.add( metapatch::first_ligand_res, "first res" ).def( 1 );
		option.add( metapatch::last_ligand_res, "last res" ).def( 2 );
		option.add( metapatch::scoring_mode, "Scoring mode. 1 = bound state, 2 = unbound minus bound. Default 1" ).def( 1 );
		option.add( metapatch::score_threshold, "Scoring threshold. Default -2" ).def( -2 );
		option.add( metapatch::top_n, "Number to report at the end. Default 20" ).def( 20 );
		option.add( metapatch::also_pack, "Also pack. Default false" ).def( false );
		option.add( metapatch::taboo_threshold, "Score above which to discard a modification from further consideration." ).def( 5 );

		// init command line options
		devel::init(argc, argv);

		ResidueTypeSetCAP rts( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );
		ScoreFunctionOP score_fxn = get_score_function();
		// Soft repulsion.
		score_fxn->set_weight( fa_rep, 0.1 );

		Pose pose;
		import_pose::pose_from_file( pose, option[ in::file::s ].value()[1] , core::import_pose::PDB_file);

		Size ligand_begin = option[ metapatch::first_ligand_res ].value();
		Size ligand_end = option[ metapatch::last_ligand_res ].value();
		Real threshold = option[ metapatch::score_threshold ].value();
		Size mode = option[ metapatch::scoring_mode ].value();
		bool pack = option[ metapatch::also_pack ].value();
		Real taboo = option[ metapatch::taboo_threshold ].value();

		Size n_to_save = option[ metapatch::top_n ].value();

		utility::vector1< std::string > metapatch_names;

		//metapatch_names.push_back( "brominated" );
		metapatch_names.push_back( "chlorinated" );
		//metapatch_names.push_back( "fluorinated" );
		//metapatch_names.push_back( "iodinated" );
		metapatch_names.push_back( "methylated" );

		//metapatch_names.push_back( "methoxylated" );
		//metapatch_names.push_back( "ethylated" );
		//metapatch_names.push_back( "hydroxylated" );

		kinematics::MoveMapOP mm( new kinematics::MoveMap() );
		for ( Size ii = 1; ii <= pose.size(); ++ii ) {
			mm->set_bb(  ii, true );
			mm->set_chi( ii, true );
		}

		//utility::vector1< Real > final_scores;
		//utility::vector1< std::string > summary_lines;

		MetapatchEnumeration metapatch_tester( pose, mm, score_fxn, mode, pack, threshold, metapatch_names, taboo );
		for ( Size resi = ligand_begin; resi <= ligand_end; ++resi ) {
			TR << "Starting score is " << metapatch_tester.score() << std::endl;
			metapatch_tester.generate_metapatched_variants( resi );
		}

		utility::vector1< Real > final_scores = metapatch_tester.final_scores();;
		utility::vector1< std::string > summary_lines =  metapatch_tester.summary_lines();;

		if ( final_scores.size() >= n_to_save ) {
			TR << "Best: " << std::endl;
			// Okay, now report the best 10 mutants.
			utility::vector1< Size > indices_best_n( n_to_save );
			utility::arg_least_several( final_scores, indices_best_n );
			for ( Size i = 1; i <= n_to_save; ++i ) {
				TR << std::setw(8) << final_scores[ indices_best_n[ i ] ] << "\t" << summary_lines[ indices_best_n[ i ] ] << std::endl;
			}

			TR << "Worst: " << std::endl;
			// Okay, now report the best n mutants.
			utility::vector1< Size > indices_worst_n( n_to_save );
			utility::arg_greatest_several( final_scores, indices_worst_n );
			for ( Size i = 1; i <= n_to_save; ++i ) {
				TR << std::setw(8) << final_scores[ indices_worst_n[ i ] ] << "\t" << summary_lines[ indices_worst_n[ i ] ] << std::endl;
			}

		} else {
			TR << "All: " << std::endl;
			for ( Size i = 1; i <= final_scores.size(); ++i ) {
				TR << std::setw(8) << final_scores[ i ] << "\t" << summary_lines[ i ] << std::endl;
			}
		}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}//main

