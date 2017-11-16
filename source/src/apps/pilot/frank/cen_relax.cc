/// @file
/// @brief


#include <devel/init.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/simple_moves/MissingDensityToJumpMover.hh>

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

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>


static basic::Tracer TR( "cen_relax" );

OPT_1GRP_KEY(Boolean, min, debug)


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

		bool debug_derivs = option[ OptionKeys::min::debug ]();
		core::optimization::CartesianMinimizer minimizer;
		core::optimization::MinimizerOptions minoptions( "lbfgs_armijo_nonmonotone", 1e-3, true, debug_derivs, debug_derivs );
		minoptions.nblist_auto_update( true );
		minoptions.max_iter( 200 );

		// scorefunction0 -- fix bad chainbreaks
		core::scoring::ScoreFunctionOP scorefxn0 = core::scoring::get_score_function();
		scorefxn0->reset();
		scorefxn0->set_weight( core::scoring::vdw, 0.1 );
		scorefxn0->set_weight( core::scoring::cart_bonded, 0.1 );
		core::scoring::methods::EnergyMethodOptions options0(scorefxn0->energy_method_options());
		options0.set_cartesian_bonded_linear(true);
		options0.set_cartesian_bonded_parameters(10,1,0,0,0);
		scorefxn0->set_energy_method_options(options0);
		minimizer.run( pose, mm, *scorefxn0, minoptions );

		core::scoring::ScoreFunctionOP scorefxn1 = core::scoring::ScoreFunctionFactory::create_score_function("score4_smooth_cart");
		if ( option[ score::weights ].user() ) {
			scorefxn1 = core::scoring::get_score_function();
		}
		core::scoring::methods::EnergyMethodOptions options1(scorefxn1->energy_method_options());
		options1.set_cartesian_bonded_linear(true);
		scorefxn1->set_energy_method_options(options1);

		// ramp vdw term
		core::Real max_vdw = scorefxn1->get_weight( core::scoring::vdw );
		int ncyc = option[ relax::default_repeats ]();

		for ( int i=1; i<=ncyc; ++i ) {
			scorefxn1->set_weight( core::scoring::vdw, ((Real)i)/(ncyc) * max_vdw );
			(*scorefxn1)(pose);
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
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	SequenceMoverOP seq( new SequenceMover() );
	if ( option[ symmetry::symmetry_definition ].user() ) seq->add_mover( new protocols::simple_moves::symmetry::SetupForSymmetryMover() );
	if ( option[ constraints::cst_file ].user() ) seq->add_mover( new protocols::simple_moves::ConstraintSetMover(  ) );
	seq->add_mover( new protocols::simple_moves::MissingDensityToJumpMover );
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
	try {
		NEW_OPT(min::debug, "debug derivs?", false);

		// initialize option and random number system
		devel::init( argc, argv );

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
