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
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/ncbb/util.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>

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
using namespace protocols::simple_moves::hbs;
using namespace protocols::rigid;
using namespace protocols::toolbox;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static basic::Tracer TR("TrivialAlascanMover");

class TrivialAlascanMover : public Mover {

public:

	//default ctor
	TrivialAlascanMover(): Mover("TrivialAlascanMover"){}

	//default dtor
	virtual ~TrivialAlascanMover(){}

	//methods
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "TrivialAlascanMover"; }

};

typedef utility::pointer::shared_ptr< TrivialAlascanMover > TrivialAlascanMoverOP;
typedef utility::pointer::shared_ptr< TrivialAlascanMover const > TrivialAlascanMoverCOP;


int
main( int argc, char* argv[] )
{
	try {

		//utility::vector1< core::Size > empty_vector(0);

		// init command line options
		devel::init(argc, argv);

		//create mover instance
		TrivialAlascanMoverOP TAM_mover( new TrivialAlascanMover() );

		protocols::ncbb::setup_filter_stats();

		//call job distributor
		protocols::jd2::JobDistributor::get_instance()->go( TAM_mover );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}//main

void
mutate_residue_to_ala( Pose & pose, Size i ) {

	core::chemical::ResidueTypeSetCOP residue_type_set( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );

	core::conformation::Residue res = pose.residue( i );
	core::chemical::ResidueTypeCOPs possible_types = residue_type_set->name3_map_DO_NOT_USE( "ALA" );
	utility::vector1< std::string > variant_types = res.type().properties().get_list_of_variants();

	// Run through all possible new residue types.
	for ( chemical::ResidueTypeCOPs::const_iterator
			type_iter = possible_types.begin(), type_end = possible_types.end();
			type_iter != type_end; ++type_iter ) {
		bool perfect_match( true ); // indicates this type has all the same variant types as the old residue

		//TR << "contemplating " << (*type_iter)->name() << std::endl;
		for ( Size kk = 1; kk <= variant_types.size(); ++kk ) {
			//TR << "checking for variant type " << variant_types[ kk ]<< std::endl;
			if ( ! (*type_iter)->has_variant_type( variant_types[ kk ] ) ) {
				perfect_match = false;
				break;
			}
		}

		if ( perfect_match ) { // Do replacement.
			//TR << (*type_iter)->name() << " success!" << std::endl;
			ResidueOP new_res = ResidueFactory::create_residue( **type_iter, res, pose.conformation() );
			core::conformation::copy_residue_coordinates_and_rebuild_missing_atoms( res, *new_res, pose.conformation() );
			pose.conformation().replace_residue( i, *new_res, false );

			return;
		}
	}
}

void
TrivialAlascanMover::apply(
	core::pose::Pose & pose
)
{

	scoring::ScoreFunctionOP score_fxn = get_score_function();
	scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);
	score_fxn->set_weight( ref, 0 );
	score_fxn->set_weight( unfolded, 0 );

	Real wt_score = (*score_fxn)( pose );
	TR << "Wildtype scores " << wt_score << std::endl;

	// Don't even bother to calculate the interface yet
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {

		Pose pose_copy( pose );
		mutate_residue_to_ala( pose_copy, i );

		Real mut_score = (*score_fxn)( pose_copy );
		if ( mut_score - wt_score > 2 ) {
			TR << pose.pdb_info()->pose2pdb( i ) << (mut_score - wt_score) << std::endl;
		}


	}
}

