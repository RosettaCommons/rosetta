/// @file
/// @brief


#include <devel/init.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <protocols/moves/PackRotamersMover.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>


#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <utility/excn/Exceptions.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>


static basic::Tracer TR("cen_relax");


///////////////////////////////////////////////////////////////////////////////

class CenRelaxMover : public protocols::moves::Mover {
public:
	CenRelaxMover(){}
	void apply( core::pose::Pose & pose) {
		using namespace core;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		core::kinematics::MoveMap mm;
		mm.set_bb  ( true ); mm.set_chi ( true ); mm.set_jump( true );

		core::optimization::CartesianMinimizer minimizer;
		core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 1e-3, true, false, false );
		minoptions.nblist_auto_update( true );
		minoptions.max_iter( 200 );

		// scorefunction0 -- weakly restrain bondlengths
		core::scoring::ScoreFunctionOP scorefxn0 = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth_cart");
		if ( option[ score::weights ].user() ) {
			scorefxn0 = core::scoring::getScoreFunction();
		}

		core::scoring::methods::EnergyMethodOptions options0(scorefxn0->energy_method_options());
		options0.set_cartesian_bonded_linear(true);
		options0.set_cartesian_bonded_parameters(2,0,0,0);
		scorefxn0->set_energy_method_options(options0);

		minimizer.run( pose, mm, *scorefxn0, minoptions );

		core::scoring::ScoreFunctionOP scorefxn1 = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth_cart");
		if ( option[ score::weights ].user() ) {
			scorefxn1 = core::scoring::getScoreFunction();
		}
		core::scoring::methods::EnergyMethodOptions options1(scorefxn1->energy_method_options());
		options1.set_cartesian_bonded_linear(true);
		scorefxn1->set_energy_method_options(options1);

		// ramp chainbreak term
		core::Real max_cart = scorefxn1->get_weight( core::scoring::cart_bonded );
		int ncyc = option[ relax::default_repeats ]();

		for (int i=1; i<=ncyc; ++i) {
			scorefxn1->set_weight( core::scoring::cart_bonded, ((Real)i)/(ncyc) * max_cart );
			minimizer.run( pose, mm, *scorefxn1, minoptions );
			scorefxn1->show( TR.Debug, pose  );
		}
	}

	virtual std::string get_name() const {
		return "CenRelaxMover";
	}
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;

	SequenceMoverOP seq( new SequenceMover() );
	seq->add_mover( new protocols::simple_moves::symmetry::SetupForSymmetryMover() );
	seq->add_mover( new protocols::simple_moves::SwitchResidueTypeSetMover("centroid") );
	seq->add_mover( new CenRelaxMover() );

	try{
		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch ( utility::excn::EXCN_Base& excn ) {
		TR.Error << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {

	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
}
