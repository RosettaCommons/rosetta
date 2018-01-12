/// @file
/// @brief


#include <devel/init.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/minimization_packing/PackRotamersMover.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/minimization_packing/symmetry/SymMinMover.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/Remarks.hh>
#include <core/id/AtomID.hh>


#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/rbsegment_relax/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/AtomTreeMultifunc.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh>
#include <core/optimization/symmetry/SymAtomTreeMultifunc.hh>

#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/optimization/CartesianMinimizer.hh>

#include <utility/excn/Exceptions.hh>


#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/kinematics/MoveMap.hh>


#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>


///////////////////////////////////////////////////////////////////////////////

class CustomMover : public protocols::moves::Mover {
public:
	CustomMover(){}
	void apply( core::pose::Pose & pose) {
		using namespace core;
		using namespace core::optimization;
		using namespace core::optimization::symmetry;

		core::scoring::ScoreFunctionOP score_function_ref = core::scoring::get_score_function();
		core::scoring::ScoreFunctionOP rosetta_scorefxn = new core::scoring::ScoreFunction();

		core::kinematics::MoveMap move_map;
		move_map.set_bb  ( true );
		move_map.set_chi ( true );
		move_map.set_jump( true );
		if (core::pose::symmetry::is_symmetric(pose)) {
			core::conformation::symmetry::SymmetricConformation const & symm_conf (
					dynamic_cast<core::conformation::symmetry::SymmetricConformation const & > ( pose.conformation() ) );
			core::conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

			// symmetrize scorefunct & movemap
			rosetta_scorefxn = core::scoring::symmetry::symmetrize_scorefunction( *rosetta_scorefxn );
			core::pose::symmetry::make_symmetric_movemap( pose, move_map );
		}

		// compute gradients using both scorefunctions
		CartesianMinimizerMap min_map;
		min_map.setup( pose, move_map );
		Multivec vars( min_map.ndofs() ), dExtal_dvars;
		min_map.copy_dofs_from_pose( pose, vars );

		utility::vector1< Multivec > dEros_dvars(16);
		for (int ii=0; ii<16; ++ii) {
			rosetta_scorefxn->set_weight( core::scoring::fa_atr      , (ii==0 || ii==1)? score_function_ref->get_weight(core::scoring::fa_atr) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::fa_rep      , (ii==0 || ii==2)? score_function_ref->get_weight(core::scoring::fa_rep) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::fa_sol      , (ii==0 || ii==3)? score_function_ref->get_weight(core::scoring::fa_sol) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::fa_intra_rep      , (ii==0 || ii==2)? score_function_ref->get_weight(core::scoring::fa_intra_rep) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::fa_sol      , (ii==0 || ii==3)? score_function_ref->get_weight(core::scoring::fa_sol) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::fa_elec      , (ii==0 || ii==4)? score_function_ref->get_weight(core::scoring::fa_elec) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::hbond_lr_bb , (ii==0 || ii==5)? score_function_ref->get_weight(core::scoring::hbond_lr_bb) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::hbond_sr_bb , (ii==0 || ii==6)? score_function_ref->get_weight(core::scoring::hbond_sr_bb) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::hbond_bb_sc , (ii==0 || ii==7)? score_function_ref->get_weight(core::scoring::hbond_bb_sc) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::hbond_sc    , (ii==0 || ii==8)? score_function_ref->get_weight(core::scoring::hbond_sc) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::p_aa_pp     , (ii==0 || ii==9)? score_function_ref->get_weight(core::scoring::p_aa_pp) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::pro_close   , (ii==0 || ii==10)? score_function_ref->get_weight(core::scoring::pro_close) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::rama     , (ii==0 || ii==11)? score_function_ref->get_weight(core::scoring::rama) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::omega   , (ii==0 || ii==12)? score_function_ref->get_weight(core::scoring::omega) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::cart_bonded_angle, (ii==0 || ii==13)? score_function_ref->get_weight(core::scoring::cart_bonded) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::cart_bonded_length , (ii==0 || ii==14)? score_function_ref->get_weight(core::scoring::cart_bonded) : 0.0 );
			rosetta_scorefxn->set_weight( core::scoring::cart_bonded_torsion , (ii==0 || ii==15)? score_function_ref->get_weight(core::scoring::cart_bonded) : 0.0 );

			(*rosetta_scorefxn)(pose);  // score pose first
			rosetta_scorefxn->setup_for_minimizing( pose, min_map );
			CartesianMultifunc f_ros( pose, min_map, *rosetta_scorefxn, false, false );
			f_ros.dfunc( vars, dEros_dvars[ii+1] );
		}

		// report
		for (core::Size counter=0; counter<dEros_dvars[1].size()/3; ++counter) {
			id::AtomID id = min_map.get_atom( counter+1 );
			core::conformation::Residue const & rsd_i = pose.residue( id.rsd() );
			if (id.atomno() <= rsd_i.nheavyatoms()) {
				std::cerr << "GRAD   " << rsd_i.name3() << " " << id.rsd() << " " << rsd_i.atom_name( id.atomno() );
				for (int ii=0; ii<16; ++ii) {
					core::Real grad_k = std::sqrt (
						dEros_dvars[ii+1][3*counter+1]*dEros_dvars[ii+1][3*counter+1] +
						dEros_dvars[ii+1][3*counter+2]*dEros_dvars[ii+1][3*counter+2] +
						dEros_dvars[ii+1][3*counter+3]*dEros_dvars[ii+1][3*counter+3] );
					std::cerr << " " << grad_k;
				}
				std::cerr << std::endl;
			}
		}
	}

	virtual std::string get_name() const {
		return "CustomMover";
	}
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;
	using namespace protocols::simple_moves::symmetry;

	SequenceMoverOP seq( new SequenceMover() );
	//seq->add_mover( new SetupForSymmetryMover() );
	seq->add_mover( new CustomMover() );

	try{
		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch (utility::excn::Exception& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
    try {
	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
    } catch (utility::excn::Exception const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
    }
